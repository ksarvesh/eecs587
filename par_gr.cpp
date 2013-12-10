/****************************************************************************************
 * 
 * Program: par_qprf_gr_v2.cpp
 * Authors: Adrian Montero and Sarversh Kirthivasan
 *
 * DESCRIPTION:
 * This program performs the paralel (shared memory) version of the Goldberg's maximum
 * flow algorithm. 
 *
 * The general approach is to maintain a global work queue which keeps track of the
 * active vertices, applies push operations to them, possibly relables them, and puts
 * any newly activated vertices back in the queue. This version also dynamically changes
 * the size of the inQueue and the outQueue there by making it more efficient.
 * 
 * This is different from version 1 in that this does not do a busy wait while trying 
 * to do the global relabel operation. It acquires all vertexLocks instead of the 
 * queueLock. This approach is more parallel. 
 *
 * ACTIVE QUEUE:
 * The queue for active vertices is divided in two: one shared, and the other local to 
 * each processor. The local queue is further subdivided into inqueue and outqueue. A
 * processor takes a vertex from inqueue and discharges it. When inqueue is empty, a 
 * processor gets a new batch of b vertices from the shared queue and stores them in its
 * local inqueue. When outqueue gets full, the processor places the entire content in
 * back to the shared queue. The size of b of inqueue and outqueue is dynamically 
 * adjusted based on the number of idle processors. 
 *
 * LOCKS:
 * There are two main omp locks used in the algorithm:
 *
 * 1.) active queue lock which guarantees mutually exclusive access to the queue. Only
 *     one processor can push and/or pop from the queue at a time.
 *     Lock name: queueLock
 *
 * 2.) vertex locks is an array of locks that prevent access to a specific vertex during
 *     push and relabel operations.
 *     Lock name: vertexLock
 *
 * GLOBAL RELABELING:
 * Done by performing a backwards breadth-first search (BFS) on the residual graph
 * from the sink and from the source. The labels are updated to be the exact 
 * distance in number of edges from each vertex to the sink or, when the vertex 
 * cannot reach the sink, the exact distance to the source plus n. 
 * 
 *
 ***************************************************************************************/
#include <cstdlib>
#include <iostream>
#include <vector>
#include <queue>
#include <fstream>
#include <assert.h>
#include <omp.h>

using namespace std;

#define MIN( X,Y ) X < Y ? X : Y
#define TRUE  1
#define FALSE 0
#define INFINITE 10000000
#define NUM_THREADS 12
#define DEBUG 0
#define DEBUG_temp 0
#define TIMING 0
#define BUSY_WAIT 0
#define GR_Q_SIZE 4

// {s,t} read from input file
int sourceId, sinkId, numIdleProcessors=0, isCompleted=0;
int numVertices, numEdges;

//This is the lock used for debugging
omp_lock_t printLock;

//This basically says what the size of the inQueue should be.
//This will be changed by holding the queueLock since we do
//not want anything to access the queue while the number of
//vertices to be grabbed from the queue is being modified.
int inputQueueSize = 10;

//This is to keep a track of how many discharges have happened.
//Once the num of discharges (over all processors) has reached 
//a threshold amount, it will reset the count to 0 and call 
//changeBufferSize and doGlobalRelabeling.
int totalNumDischarges= 0;


// data structure to represent graph vertices 
struct vertex{
    int height;					
    int excessFlow;
    int waveNumber;
    //int isActive;			// FALSE - vertex is not in active queue
											// TRUE  - vertex already in active queue
};

// data structure to represent graph edges
struct edge{
    int capacity;
    int flow;
    int from;
    int to;
};

// global relabeling data structures
omp_lock_t bfsLock;
queue<int> 	bfsQueue;
vector<int> marked;
bool isGlobalRelabelingInProgress= false;
bool bfsQueueInUse= false;
bool isSecondBfs= false;
int waveNumber= 0;

/*****************************************************************
 * The globalRelabel function performs a backwards breadth-first
 * search (BFS) on the residual graph from the sink and from the 
 * source. The labels are updated to be the exact distance in number
 * of edges from each vertex to the sink or, when a vertex cannot 
 * reach the sink, the exact distance to the source plus n. 
 *
 * Global relabeling is applied periodically to improve the 
 * performance of the push/relabel algorithm.
 *
 ****************************************************************/
void globalRelabel( vector<vector<int> >& adjList, vector<edge>& edgeList, vector<vertex>& vertexList, omp_lock_t* queueLock, vector<omp_lock_t>& vertexLock)
{
	int currentEdgeId, residue, v, w, e;	
    int numExplored= 0;

    omp_set_lock(&bfsLock);

    while(numExplored < GR_Q_SIZE)
	{
		// reset breadth-first search after done with sink 
		if( bfsQueue.empty() ){
			omp_set_lock(queueLock);
            if(isSecondBfs){
                isGlobalRelabelingInProgress= false;
                isSecondBfs= false;
                omp_unset_lock(queueLock);
                break;
            }else{
                isSecondBfs= true;
                bfsQueue.push(sourceId);
                omp_unset_lock(queueLock);
           }
        }

		currentEdgeId = 0;
		v = bfsQueue.front();
		bfsQueue.pop();
        
        omp_set_lock(&vertexLock[v]);
        int vHeight= vertexList[v].height;
        omp_unset_lock(&vertexLock[v]);

		while( currentEdgeId < adjList[v].size() )
		{ // iterate through each v node edge
			e = adjList[v][currentEdgeId];
			
			// need to know if forward or backward edge
			if(edgeList[e].from == v)
			{ // backward edge
				w= edgeList[e].to;
				omp_set_lock(&vertexLock[w]);
                residue = edgeList[e].flow;
			}
			else
			{ // forward edge
				w= edgeList[e].from;
				omp_set_lock(&vertexLock[w]);
				residue = edgeList[e].capacity - edgeList[e].flow;
			}

			// if vertex w has not been visited, push into queue and mark
			if( !marked[w] && residue>0 )
			{
				// relabel w/ depth
				vertexList[w].height = vHeight+1;
				bfsQueue.push( w );
				marked[w] = TRUE;
			}

            //release the vertex locks:
			omp_unset_lock(&vertexLock[w]);

			// next edge
			currentEdgeId++;
		}

        numExplored++;
	}

    omp_unset_lock(&bfsLock);
}

/*****************************************************************
 * This funtion changes the inputQueue size depending on the 
 * number of idle processors and the number of vertices in the
 * shared queue (activeVertexQueue). The queueLock should be
 * unest while calling this function. This will grab all vertex
 * locks before relabeling since no discharging should occur while
 * global relabeling is occuring.
 ****************************************************************/
void startGlobalRelabel(vector<omp_lock_t>& vertexLock, vector<vector<int> >& adjList, vector<edge>& edgeList, vector<vertex>& vertexList, omp_lock_t* queueLock){
    
    omp_set_lock(&bfsLock);
    
    for(int i=0; i<numVertices; i++){
        marked[i]=false;
    }

    marked[sinkId]=true;
    marked[sourceId]=true;

    waveNumber+=1;
    
    vertex_lock(vertexLock, sourceId, sinkId, true);
    vertexList[sourceId].waveNumber= waveNumber;
    vertexList[sinkId].waveNumber= waveNumber;
    vertex_lock(vertexLock, sourceId, sinkId, false);
    
    bfsQueue.push(sinkId);

    omp_unset_lock(&bfsLock);

    globalRelabel(adjList, edgeList, vertexList, queueLock, vertexLock);

    return;
}

/*****************************************************************
 * This funtion changes the inputQueue size depending on the 
 * number of idle processors and the number of vertices in the
 * shared queue (activeVertexQueue). The queueLock should be
 * held while calling this function.
 ****************************************************************/
void changeBufferSize(queue<int>& activeVertexQueue){
    
    //Time the changeBufferSize function:
    double startTime= omp_get_wtime();

    //This function should be called with the queueLock being held
    int numActiveProcessors = NUM_THREADS - numIdleProcessors;
    int numActiveVertices = activeVertexQueue.size();
    
    if(DEBUG){
        omp_set_lock(&printLock);
        cout<<omp_get_thread_num()<<" is changing the buffer size"<<endl;
        omp_unset_lock(&printLock);
    }

    
    if(numIdleProcessors >= 2 || numIdleProcessors >= .15*NUM_THREADS){
        inputQueueSize/=2;
    }else if( numActiveProcessors + ((double)numActiveVertices/inputQueueSize) > 1.5*NUM_THREADS){
        inputQueueSize*=2;
    }

    //If the inputQueueSize is made zero, change it to size 1.
    //This is to prevent pre-mature termination of the algorithm.
    inputQueueSize= inputQueueSize == 0 ? 1 : inputQueueSize;
    inputQueueSize= inputQueueSize > 16 ? 16 : inputQueueSize;

    //Time the global relabel operation and print:
    double endTime= omp_get_wtime();

    if(TIMING){
        omp_set_lock(&printLock);
        cout<<omp_get_thread_num()<<" took "<<endTime-startTime<<" time to change the bufferSize"<<endl;
        omp_unset_lock(&printLock);
    }


}

/*****************************************************************
 * This function takes to vertex indexes and locks the vertex 
 * with the minimum index first or unlocks the vertex with the 
 * minimum index last. Avoids dealocks among processor
 * trying to push the same edge.
 *
 ****************************************************************/
void vertex_lock( vector<omp_lock_t>& vertexLock, int v, int w, int lock )
{
	// get locks v and w, to avoid deadlock grab first vertex with the
	// minimum index
	if( MIN( v,w ) == v )
	{ // v is min
		// set or unset lock
		if( lock==TRUE )
		{ // set lock of min vertex first
			omp_set_lock( &vertexLock[v] );
			omp_set_lock( &vertexLock[w] );
		}
		else
		{ // unset lock of min vertex last
			omp_unset_lock( &vertexLock[w] );
			omp_unset_lock( &vertexLock[v] );
		}
	}
	else
	{ // w is min
		if( lock==TRUE )
		{ // set lock 
			omp_set_lock( &vertexLock[w] );
			omp_set_lock( &vertexLock[v] );
		}
		else
		{ // unset lock
			omp_unset_lock( &vertexLock[v] );
			omp_unset_lock( &vertexLock[w] );
		}
	}
}


/*****************************************************************
 * Function createNewEdge is used by the initialization subroutine
 * to create edges dynamically as they are read from the input file
 * NOTE: after initialization edge data structure becomes static
 *
 ****************************************************************/
edge createNewEdge(int from, int to, int capacity, int flow){
    edge newEdge;
    newEdge.from= from;
    newEdge.to= to;
    newEdge.capacity= capacity;
    newEdge.flow= flow;

    return(newEdge);
}

/*****************************************************************
 * Populate data structures and initialize maximum flow in graph
 * from input text file.
 *
 ****************************************************************/
int initialize(string fileName, vector<vertex>& vertexList, vector<edge>& edgeList, vector<vector<int> >& adjList, queue<int>& activeVertexQueue, vector<omp_lock_t>& vertexLock, omp_lock_t *queueLock_ptr ){

    //Open file with name fileName:
    ifstream inFile(fileName.c_str());

    if(!inFile.is_open()){
        cerr<<"error opening file"<<endl;
        return(-1);
    }

    string line;
    
    //Get num vertices and num edges from the first two lines of the file:
    //Also, get the indices of the source and the sink (global variables)
    inFile >> numVertices >> numEdges >> sourceId >> sinkId;

    //Initilaize the marked vector to be of size numVertices
    marked.resize(numVertices,0);
   
	// Initilize queue and vertex locks
	omp_init_lock( queueLock_ptr );
    omp_init_lock(&bfsLock);
	vertexLock.resize(numVertices);
	for( int i=0;i>numVertices;i++ )
	{
		omp_init_lock( &vertexLock[i] );
	}
 
    //Initialize the vertexList such that all vertices except for the source
    //has a height of zero. Set all excess flow to zero initially.
    //Set height of the source node to be numVertices
    vertex zeroVertex;
    //zeroVertex.height = zeroVertex.excessFlow = zeroVertex.isActive = 0;
    zeroVertex.height = zeroVertex.excessFlow = zeroVertex.waveNumber= 0;
    vertexList.resize(numVertices, zeroVertex);
    vertexList[sourceId].height = numVertices;
    
    //To remove the ending \n char in the previous line (check test1.txt)
    getline(inFile, line);
    
    //Initialize the edgeList and the adjList:
    adjList.resize(numVertices);
    edgeList.resize(numEdges);

    int edgeId=0;

    while(getline(inFile, line)){
        int from, to, cap;
        sscanf(line.c_str(), "%d %d %d", &from, &to, &cap);
        
        //Assert that the vertices are valid indices:
        assert(from < numVertices && to < numVertices);
        
        //If the edge is an outgoing edge from the source, set the flow on that edge to be equal to the capacity
        //Set excess flow on the "to" vertex to be equal to the capacity, and push it onto the activeVertexQueue
        //if the to vertex is not the sink. If it is, then do not push it. Else, just populate the edgeList and adjList:
        if(from != sourceId){
            edgeList[edgeId]= createNewEdge(from, to, cap, 0);
        }else{
            edgeList[edgeId]= createNewEdge(from, to, cap, cap);
            vertexList[to].excessFlow = cap;
            //Do not push sink onto the activeVerticesQueue
            if(to!=sinkId){
                //vertexList[to].isActive = 1;
                activeVertexQueue.push(to);
            }
        }

        adjList[from].push_back(edgeId);
        adjList[to].push_back(edgeId);
    
        //Increment the edgeId count
        edgeId++;
    }

    //Assert that the number of edges is equal to the one specified in the file initially:
    assert(edgeId == numEdges);

    return(0);
}

/*****************************************************************
 * The push operation consists of pushing excess from vertex v
 * to w.
 ****************************************************************/
void push( queue<int>& outQueue, vector<vertex>& vertexList, vector<edge>& edgeList, int v, int w, int e )
{
	int send;
	int capacity = edgeList[e].capacity;
	int flow		 = edgeList[e].flow;
	int excess   = vertexList[v].excessFlow;

	if( edgeList[e].from==v )
	{ // forward push
		send = MIN( excess, capacity-flow );
		edgeList[e].flow += send;	
	}
	else
	{ // backward push
		send = MIN( excess, flow );
		edgeList[e].flow -= send;
	}

    if(DEBUG){
        omp_set_lock(&printLock);
        cout<<"thread "<<omp_get_thread_num()<<" is pushing "<<send<< " from "<<v<<" to "<<w<<endl;
	    omp_unset_lock(&printLock);
    }

	// update excessFlow 
	vertexList[v].excessFlow -= send;
	vertexList[w].excessFlow += send;

	// add vertex w if it is not already in active queue 
	// and vertex - {s,t}
	//if( w!=sinkId && w!=sourceId &&  vertexList[w].isActive == 0 ) 
	if( w!=sinkId && w!=sourceId &&  vertexList[w].excessFlow == send ) 
	{
        //cout<< "pushing v onto the queue " << w<<endl;
		outQueue.push(w);
		//vertexList[w].isActive= 1;
	}
		
}

/*****************************************************************
 * The relabel operation consists of setting node v's height 
 * to 1 higher than its lowest neighbor that has residual flow.
 ****************************************************************/
void relabel( vector<vector<int> >& adjList, vector<vertex>& vertexList, vector<edge>& edgeList, int v )
{
	int e, w, residue;

	// number of edges adjancent to vertex v
	int size = adjList[v].size();
	int min_height = INFINITE;

	// iterate through each adjacent edge and check:
	// 		- edge is in residual graph
	for( int i=0; i<size;i++ )
	{   
        e= adjList[v][i];

		// calculate edge residue from v->e
		if( edgeList[e].from==v )
		{ // forward edge
			residue = edgeList[e].capacity - edgeList[e].flow;
			w = edgeList[e].to;
		}
		else
		{ // backward edge
			residue = edgeList[e].flow;
			w = edgeList[e].from;
		}
		
		// if edge has residue, set vertex v height to 1 higher than w's 	
        if( residue>0 )
		{ // edge is in residual graph
           // cout<<"relabeling heights: "<<min_height<<" "<<vertexList[w].height<<endl;
			min_height = MIN( min_height,vertexList[w].height );
		}
	}
	vertexList[v].height = min_height + 1;
	
    if(DEBUG){
        omp_set_lock(&printLock);
        cout<<"thread "<<omp_get_thread_num()<<" has relabeled height for "<<v<<" : "<<vertexList[v].height<<endl;
	    omp_unset_lock(&printLock);
    }
}

/*****************************************************************
 * The discharge operation consists applying push/relabel to vertex
 * at least until the excess becomes zero or the label of the 
 * vertext increases. 
 ****************************************************************/
void discharge( queue<int>& outQueue, vector<vertex>& vertexList, vector<edge>& edgeList, vector<vector<int> >& adjList, int v, vector<omp_lock_t>& vertexLock )
{	
	int e, residue, w, currentEdgeId=0;

	// try to push to all edges until excess becomes 0 or there are no more
	// edges to push to. 
	while(currentEdgeId < adjList[v].size() && vertexList[v].excessFlow>0){
		 e= adjList[v][currentEdgeId];
	
		// calculate residue, need to know if forward or backward edge
		if(edgeList[e].from == v){
			w= edgeList[e].to;

			// lock vertices before calculating residue
			vertex_lock( vertexLock,v,w,TRUE );	
			residue= edgeList[e].capacity - edgeList[e].flow;

		}else{
			w= edgeList[e].from;
		
			// lock vertices before calculating residue
			vertex_lock( vertexLock,v,w,TRUE );
			residue= edgeList[e].flow;
		}
		
		// push if edge has residue and height of v > w, else go to next edge
		// in vertex's v list 
		if(residue > 0  &&  vertexList[v].height > vertexList[w].height && vertexList[v].waveNumber==vertexList[w].waveNumber){
			push( outQueue, vertexList, edgeList, v, w, e);
		}else{
			currentEdgeId++;
		}

		// release locks associated with previous edge
		vertex_lock( vertexLock,v,w,FALSE );
	}

	// if there is any excess left in v, then relabel and put back at the rear
	// of the queue and mark vertex as active.
	if(vertexList[v].excessFlow > 0){

		// get lock of vertex v before relabel
		omp_set_lock( &vertexLock[v] );
		// DEBUG
        if(DEBUG){
            omp_set_lock(&printLock);
            cout<< "thread "<<omp_get_thread_num()<<" is relabeling vertex "<<v<<endl;
            omp_unset_lock(&printLock);
        }
		relabel(adjList, vertexList, edgeList, v);
		// release lock of vertex v after relabel
		//vertexList[v].isActive= 1;
		omp_unset_lock( &vertexLock[v] );

  	    outQueue.push(v);
	}
	
	return;
}

/*****************************************************************
 * This function takes in two queues as inputs and pushes from the
 * second queue to the frist queue untill the frist queue reaches
 * its max possible size set as a #define, or until the second queue
 * becomes empty
 * ****************************************************************/


void getNewVertex(queue<int>& inQueue, queue<int>&activeVertexQueue, vector<omp_lock_t>& vertexLock, vector<vertex>& vertexList){
    int minSize= inputQueueSize - inQueue.size();
    minSize= minSize > 0 ? minSize : 0;
    
   
   int numNewVertices= MIN(minSize, activeVertexQueue.size());

    for(int i=0; i<numNewVertices; i++){
        int v= activeVertexQueue.front();
		activeVertexQueue.pop();

        if(DEBUG){
         omp_set_lock(&printLock);
         cout<<omp_get_thread_num()<<" is popping vertex "<<v<<" from the shared queue"<<endl;
         omp_unset_lock(&printLock);
        }
		
        inQueue.push(v);
   }

    return;
}



/*****************************************************************
 * This function pushes from queue 1 to queue 2 all its elements.
 * This is implemented to push all the excess active vertics from
 * the outQueue of a processor to the shared activeVertexQueue to
 * enable other possibly free processors to take up the work and 
 * start processing those vertices.
 * ****************************************************************/


void pushNewVertex(queue<int>& outQueue, queue<int>& activeVertexQueue){
   while(!outQueue.empty()){
       int v= outQueue.front();
       if(DEBUG){
       	omp_set_lock(&printLock);
        cout<<omp_get_thread_num()<<" is pushing vertex "<<v<<" into the shared queue"<<endl;
        omp_unset_lock(&printLock);
       }
			 outQueue.pop();
       activeVertexQueue.push(v);
   }
   return;
}


/*****************************************************************
 * The startParallelAlgo takes as input the initialized
 * vertex, edge, adjacency list and queue data structures.
 * Each processor tries to grab a set of vertices from the active
 * vertex queue and works on that set. Once all its out queue is 
 * full, it will push the out queue onto the active vertex queue.
 * If its inqueue is empty, it will push its non empty out queue
 * and put itself on the idleProcessor waiting queue. The elements
 * of the idleProcessor wait queue are woken up when ever there
 * are new elements pushed onto the activeVertexQueue.
 ****************************************************************/
void startParallelAlgo(queue<int>& activeVertexQueue, vector<vertex>& vertexList, vector<edge>& edgeList, vector<vector<int> >& adjList, vector<omp_lock_t>& vertexLock, omp_lock_t* queueLock){
    queue<int> inQueue, outQueue;
    
    omp_set_lock(queueLock);
    
    while(1){

        /*if(DEBUG){
            omp_set_lock(&printLock);
            cout<<omp_get_thread_num()<<" is the thread in progess"<<endl;
            omp_unset_lock(&printLock);
        }*/
        
        //Get new vertices from the shared queue if the size of the current input buffer
        //is smaller than the size we set in inputQueueSize
        getNewVertex(inQueue, activeVertexQueue, vertexLock, vertexList);
        
        if(DEBUG){
            omp_set_lock(&printLock);
            cout<<"thread "<<omp_get_thread_num()<<" has inQueueSize: "<<inQueue.size()<<" "<<endl; 
            omp_unset_lock(&printLock);
        }
        
	//Spin loop (busy wait) implementation of cpu sleep
        //that is to be woken up by an interprocessor interrupt
        while(inQueue.empty()){
            numIdleProcessors++;
            if(DEBUG){
                omp_set_lock(&printLock);
                cout<<omp_get_thread_num()<<" is stuck here with size " << inQueue.size()<<endl;
                omp_unset_lock(&printLock);
            }
            
            /*if(DEBUG){
                omp_set_lock(&printLock);
                cout<<isCompleted<<" "<<endl; 
                omp_unset_lock(&printLock);
            }*/
            
            if(numIdleProcessors == NUM_THREADS || isCompleted){
                isCompleted= 1;
        	    if(DEBUG){
                	omp_set_lock(&printLock);
            		cout<<"size of shared: "<<activeVertexQueue.size()<<" size of inQueue: "<<inQueue.size()<<" outQueue.size() "<<outQueue.size()<<endl;
            		omp_unset_lock(&printLock);
                }
				omp_unset_lock(queueLock);
                return;
            }
            
            //Busy wait loop
            double startTime= omp_get_wtime();
            omp_unset_lock(queueLock);
            //TODO:sleep(10ms);
            if(BUSY_WAIT){
		omp_set_lock(&printLock);
		cout<<"thread: "<<omp_get_thread_num()<<" busyWaiting"<<endl;
		omp_unset_lock(&printLock);
 	    }
            // while(omp_get_wtime()-startTime > .020);
            omp_set_lock(queueLock);
            
            //TODO: What if we put the -- outside the while loop,
            //will that still work. Which will be more efficient.
   	    //Decrement the number of sleeping cpus:
            numIdleProcessors--;
            getNewVertex(inQueue, activeVertexQueue, vertexLock, vertexList);
        }
	

        //Check if there have been more than 200 discharges. If so,
        //do a global relabel and check if the buffer size has to be 
        //modified depending on the number of idle processors.
        bool isGlobalRelabelingReq= false;

        if(totalNumDischarges > 200){
            totalNumDischarges = 0;
            changeBufferSize(activeVertexQueue);
            isGlobalRelabelingReq= !(isGlobalRelabelingInProgress);
	    
   	        if(DEBUG_temp){
   	         omp_set_lock(&printLock);
   	         cout<<omp_get_thread_num()<<" : inputQueueSize "<<inputQueueSize<<" activeVertices: "<<activeVertexQueue.size();
   	         cout<<" inQueue.size: "<<inQueue.size()<<" outQueue.size: "<<outQueue.size()<<" numActive: "<<NUM_THREADS - numIdleProcessors<<endl;
	         omp_unset_lock(&printLock);
   	        }
        }

    	//Keep a local copy of the inputQueueSize so that you do not have to
    	//grab the lock unneccesarily.
    	int tempInputQueueSize= inputQueueSize;

        //At this point, this processor still holds the shared activeVertexQueue 
        //lock. Also, at this point, the inQueue is non zero in size.
        //We can give up the shared queue lock at this point since we are not
        //acessing the shared variable any more.
        omp_unset_lock(queueLock);
        
        //store the current size of the inQueue so that once all discharges 
        //are complete, we can update the total number of discharges that 
        //happened in this cycle.
        int numDischarges= inQueue.size();

        //Start processing the vertices on the inQueue
        while(!inQueue.empty()){
            
            //Concurrent GR:
            omp_set_lock(queueLock);
            if(bfsQueueInUse){
                omp_unset_lock(queueLock);
            }else{
                if(isGlobalRelabelingReq){
                    isGlobalRelabelingInProgress= true;
                    bfsQueueInUse= true;
                    omp_unset_lock(queueLock);
                    //put sink node into the queue:
                    startGlobalRelabel(vertexLock, adjList, edgeList, vertexList, queueLock);
                    //set the inUse to false once 
                    //done with the global relabel 
                    omp_set_lock(queueLock);
                    bfsQueueInUse= false;
                    omp_unset_lock(queueLock);
                }else if(isGlobalRelabelingInProgress){
                    bfsQueueInUse= true;
                    omp_unset_lock(queueLock);
                    globalRelabel(adjList, edgeList, vertexList, queueLock, vertexLock);
                    //set the inUse to false once 
                    //done with the global relabel 
                    omp_set_lock(queueLock);
                    bfsQueueInUse= false;
                    omp_unset_lock(queueLock);
                }else{
                    omp_unset_lock(queueLock);
                }    
            }


            int v= inQueue.front();
            inQueue.pop();

        	if(DEBUG){
            	omp_set_lock(&printLock);
            	cout<<"thread " << omp_get_thread_num() << "is discharging vertex " <<v<<endl;
            	omp_unset_lock(&printLock);
            }

            discharge(outQueue, vertexList, edgeList, adjList, v, vertexLock);
            
            //Don't wait for the entire inQueue to be discharged 
            //before pushing out vertices to the outQueue
	        if(outQueue.size() >= tempInputQueueSize){
               	omp_set_lock(queueLock);
           	 	pushNewVertex(outQueue, activeVertexQueue);
	        	omp_unset_lock(queueLock);
	        }
        }
        
        //Re-grab the queueLock before making changes to any of the 
        //shared variables
        //
        if(BUSY_WAIT){
	    omp_set_lock(&printLock);
	    cout<<"thread: "<<omp_get_thread_num()<<" trying to grab lock"<<endl;
	    omp_unset_lock(&printLock);
 	}
        
        omp_set_lock(queueLock);
        
        if(BUSY_WAIT){
	    omp_set_lock(&printLock);
	    cout<<"thread: "<<omp_get_thread_num()<<" grabs lock"<<endl;
	    omp_unset_lock(&printLock);
 	}

	getNewVertex(inQueue, outQueue, vertexLock, vertexList);
        pushNewVertex(outQueue, activeVertexQueue);
            
        //Increment the total number of discharge operations.
        totalNumDischarges+= numDischarges;

   }

}

/*****************************************************************
 * The preflow_push function takes as input a text file, dynamically
 * initializes vertex, edge, adjacency list and queue data structures.
 * All possible flow is pushed out of source node to applicable neighboring
 * nodes. Those source adjacents are added to active queue.
 *
 ****************************************************************/
void preflow_push(string fileName)
{

	// vertices, edges, adjacency list, active queue 
	// data structures declarion
	vector<vertex> 				vertexList;
	vector<edge> 					edgeList;
	vector<vector<int> > 	adjList;
	queue<int> 						activeVertexQueue;

	// vertices and queue locks
	omp_lock_t          queueLock;
	vector<omp_lock_t> 	vertexLock;

	// populate data structures from input text file, push all possible flow out of 
	// source vertex
	int ret= initialize(fileName, vertexList, edgeList, adjList, activeVertexQueue, vertexLock, &queueLock);
	if(ret == -1){
		return;
	}

	// apply discharge opeartion until active queue is empty
	/*while(activeVertexQueue.size() > 0){

		// pop next vertex from queue, mark as not active
		int v= activeVertexQueue.front();
		vertexList[v].isActive= 0;		
		activeVertexQueue.pop();

		//DEBUG
        cout<<"discharging: "<<v<<endl;

		// discharge vertex until excess becomes 0 or vertex is relabelled
		discharge(activeVertexQueue, vertexList, edgeList, adjList, v, vertexLock);
	}*/
    
    double begin= omp_get_wtime();

    #pragma omp parallel num_threads(NUM_THREADS) 
    {
        startParallelAlgo(activeVertexQueue, vertexList, edgeList, adjList, vertexLock, &queueLock);
        
        //No thread would reach this spot unless all other threads are already sleeping
        //i.e. only when the last thread tries to go to sleep, it will realize that all 
        //other threads are already asleep, and so it will reach this spot. Hence we need 
        //a NO_WAIT, since the other threads are never going to reach this spot ever.    
    	#pragma omp single
        {
           // if(DEBUG){
                double end= omp_get_wtime();
                omp_set_lock(&printLock);
                cout<<vertexList[sinkId].excessFlow<<" is the maximum flow value"<<endl;
                printf("Elapsed time: %f secs\n", end - begin);
                omp_unset_lock(&printLock);
           // }
        }
   }
	
   return;
}
/*****************************************************************
 * Program starts here. Data structures are dynamically initialized
 * given an input text file. 
 *
 ****************************************************************/
int main(int argc, char** argv)
{
	//The second argument is the name of the file from which 
    //the function should read the data regarding the graph
    if(argc < 2){
        cout << "Usage: ./a.out <fileName>" << endl;
        return(-1);
    }

    string fileName(argv[1]);
    

    // read input text file containing graph with vertices, edges, and 
	// corresponding edges capacities
    omp_init_lock(&printLock);
	preflow_push(fileName.c_str());
    
    return(0);
}

