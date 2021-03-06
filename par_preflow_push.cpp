/****************************************************************************************
 * 
 * Program: par_preflow_push.cpp
 * Authors: Adrian Montero and Sarversh Kirthivasan
 *
 * DESCRIPTION:
 * This program performs the paralel (shared memory) version of the Goldberg's maximum
 * flow algorithm. 
 *
 * The general approach is to maintain a global work queue which keeps track of the
 * active vertices, applies push operations to them, possibly relables them, and puts
 * any newly activated vertices back in the queue. 
 *
 * ACTIVE QUEUE:
 *
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
 * 
 *
 * 
 *
 ***************************************************************************************/
#include <omp.h>
#include <iostream>
#include <vector>
#include <queue>
#include <fstream>
#include <assert.h>
#include <cstdio>

using namespace std;

#define MIN( X,Y ) X < Y ? X : Y
#define TRUE  1
#define FALSE 0
#define INFINITE 10000000
#define NUM_THREADS 2
#define IN_OUT_QUEUE_SIZE 10
#define DEBUG 0

// {s,t} read from input file
int sourceId, sinkId, numIdleProcessors=0, isCompleted=0;
omp_lock_t printLock;

// data structure to represent graph vertices 
struct vertex{
    int height;					
    int excessFlow;		 
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

/*****************************************************************
 * This function takes to vertex indexes and locks the vertex 
 * with the minimum index first or unlocks the vertex with the 
 * minimum index last. Avoids dealocks among processor
 * trying to push the same edge.
 *
 ****************************************************************/
void vertex_lock( vector<omp_lock_t>& vertexLock, int v, int w, int lock )
{  
			if(DEBUG && lock){
       	omp_set_lock(&printLock);
        cout<<omp_get_thread_num()<<" trying to acquire locks " << v << " and " << w<<endl;
        omp_unset_lock(&printLock);
       }

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

	if(DEBUG && lock){
   	omp_set_lock(&printLock);
    cout<<omp_get_thread_num()<<" has acquired locks " << v << " and " << w<<endl;
    omp_unset_lock(&printLock);
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
    int numVertices, numEdges;
    inFile >> numVertices >> numEdges >> sourceId >> sinkId;
   
		// Initilize queue and vertex locks
		omp_init_lock( queueLock_ptr );
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
    zeroVertex.height = zeroVertex.excessFlow = 0;
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
		if(residue > 0  &&  vertexList[v].height > vertexList[w].height){
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
    
    int numNewVertices= MIN(IN_OUT_QUEUE_SIZE - inQueue.size(), activeVertexQueue.size());

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
    
    if(DEBUG){
        omp_set_lock(&printLock);
        cout<<"thread "<<omp_get_thread_num()<<" is trying to acquire queue lock" <<endl;
	    omp_unset_lock(&printLock);
    }
    omp_set_lock(queueLock);
    if(DEBUG){
        omp_set_lock(&printLock);
        cout<<"thread "<<omp_get_thread_num()<<" has acquired queue lock" <<endl;
	    omp_unset_lock(&printLock);
    }
    
    while(1){

        /*if(DEBUG){
            omp_set_lock(&printLock);
            cout<<omp_get_thread_num()<<" is the thread in progess"<<endl;
            omp_unset_lock(&printLock);
        }*/
        
        //Get new vertices from the shared queue if the size of the current input buffer
        //is smaller than the size we set in IN_OUT_QUEUE_SIZE
        getNewVertex(inQueue, activeVertexQueue, vertexLock, vertexList);
        
        //Spin loop (busy wait) implementation of cpu sleep
        //that is to be woken up by an interprocessor interrupt
        while(inQueue.empty()){
            numIdleProcessors++;
            
            if(DEBUG){
                omp_set_lock(&printLock);
                cout<<"Thread: "<<omp_get_thread_num()<<" is Completed" <<endl; 
                omp_unset_lock(&printLock);
            }
            
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
            omp_unset_lock(queueLock);
            //TODO:sleep(10ms);
   				   if(DEBUG){
   				       omp_set_lock(&printLock);
   				       cout<<"thread "<<omp_get_thread_num()<<" is trying to acquire queue lock" <<endl;
	 				     omp_unset_lock(&printLock);
   				   }
            omp_set_lock(queueLock);
   				  if(DEBUG){
   				      omp_set_lock(&printLock);
   				      cout<<"thread "<<omp_get_thread_num()<<" has acquired queue lock" <<endl;
	 				    omp_unset_lock(&printLock);
   				  }
            
            //TODO: What if we put the -- outside the while loop,
            //will that still work. Which will be more efficient.
            numIdleProcessors--;
            getNewVertex(inQueue, activeVertexQueue, vertexLock, vertexList);
        }

        //At this point, this processor still holds the shared activeVertexQueue 
        //lock. Also, at this point, the inQueue is non zero in size.
        //We can give up the shared queue lock at this point since we are not
        //acessing the shared variable any more.
        omp_unset_lock(queueLock);

        //Start processing the vertices on the inQueue
        while(!inQueue.empty()){
            int v= inQueue.front();
            inQueue.pop();
            discharge(outQueue, vertexList, edgeList, adjList, v, vertexLock);
						//TODO: don't wait for the entire inQueue to be discharged before pushing out 
						//vertices to the outQueue
        }

        getNewVertex(inQueue, outQueue, vertexLock, vertexList);

   		   if(DEBUG){
   			     omp_set_lock(&printLock);
   				   cout<<"thread "<<omp_get_thread_num()<<" is trying to acquire queue lock" <<endl;
	 				   omp_unset_lock(&printLock);
   			 }
        omp_set_lock(queueLock);
   		  if(DEBUG){
   		      omp_set_lock(&printLock);
   		      cout<<"thread "<<omp_get_thread_num()<<" has acquired queue lock" <<endl;
	 		    omp_unset_lock(&printLock);
   		  }
        pushNewVertex(outQueue, activeVertexQueue);
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

	// wall clock time variable
	double begin;
	double end;

	begin = omp_get_wtime();
	
    #pragma omp parallel num_threads(NUM_THREADS) 
    {
        startParallelAlgo(activeVertexQueue, vertexList, edgeList, adjList, vertexLock, &queueLock);
      
        //No thread would reach this spot unless all other threads are already sleeping
        //i.e. only when the last thread tries to go to sleep, it will realize that all 
        //other threads are already asleep, and so it will reach this spot. Hence we need 
        //a NO_WAIT, since the other threads are never going to reach this spot ever.    
    	#pragma omp single
        {
								end = omp_get_wtime(); 
                omp_set_lock(&printLock);
                cout<<vertexList[sinkId].excessFlow<<" is the maximum flow value"<<endl;
								printf("Elapsed: %f seconds\n", end-begin);
                omp_unset_lock(&printLock);
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
	// read input text file containing graph with vertices, edges, and 
	// corresponding edges capacities
	
	

   omp_init_lock(&printLock);
	preflow_push(argv[1]);
  return(0);
}

