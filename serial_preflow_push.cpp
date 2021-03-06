#include <iostream>
#include <vector>
#include <queue>
#include <fstream>
#include <assert.h>
#include <time.h>

using namespace std;

#define MIN( X,Y ) X < Y ? X : Y
#define TRUE  1
#define FALSE 0
#define INFINITE 10000000

// {s,t} read from input file
int sourceId, sinkId;

// data structure to represent graph vertices 
struct vertex{
    int height;					
    int excessFlow;		 
    int isActive;			// FALSE - vertex is not in active queue
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
int initialize(string fileName, vector<vertex>& vertexList, vector<edge>& edgeList, vector<vector<int> >& adjList, queue<int>& activeVertexQueue){

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
    
    //Initialize the vertexList such that all vertices except for the source
    //has a height of zero. Set all excess flow to zero initially.
    //Set height of the source node to be numVertices
    vertex zeroVertex;
    zeroVertex.height = zeroVertex.excessFlow = zeroVertex.isActive = 0;
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
                vertexList[to].isActive = 1;
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
void push( queue<int>& activeVertexQueue, vector<vertex>& vertexList, vector<edge>& edgeList, int v, int w, int e )
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


   // cout<<"pushing "<<send<< " from "<<v<<" to "<<w<<endl;
	
	// update excessFlow 
	vertexList[v].excessFlow -= send;
	vertexList[w].excessFlow += send;

	// add vertex w if it is not already in active queue 
	// and vertex - {s,t}
	//if( w!=sinkId && w!=sourceId &&  vertexList[w].isActive == 0 ) 
	if( w!=sinkId && w!=sourceId &&  vertexList[w].excessFlow == send ) 
	{
        //cout<< "pushing v onto the queue " << w<<endl;
		activeVertexQueue.push(w);
		vertexList[w].isActive= 1;
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
	// DEBUG
  //	 cout<<"relabeled height for "<<v<<" : "<<vertexList[v].height<<endl;
}


/*****************************************************************
 * The discharge operation consists applying push/relabel to vertex
 * at least until the excess becomes zero or the label of the 
 * vertext increases. 
 ****************************************************************/
void discharge( queue<int>& activeVertexQueue, vector<vertex>& vertexList, vector<edge>& edgeList, vector<vector<int> >& adjList, int v )
{	
	int e, residue, w, currentEdgeId=0;

	// try to push to all edges until excess becomes 0 or there are no more
	// edges to push to. 
	while(currentEdgeId < adjList[v].size() && vertexList[v].excessFlow>0){
		 e= adjList[v][currentEdgeId];
	
		// calculate residue, need to know if forward or backward edge
		if(edgeList[e].from == v){
			residue= edgeList[e].capacity - edgeList[e].flow;
			w= edgeList[e].to;
		}else{
			residue= edgeList[e].flow;
			w= edgeList[e].from;
		}
		
		// push if edge has residue and height of v > w, else go to next edge
		// in vertex's v list 
		if(residue > 0  &&  vertexList[v].height > vertexList[w].height){
			push( activeVertexQueue, vertexList, edgeList, v, w, e);
		}else{
			currentEdgeId++;
		}
	}

	// if there is any excess left in v, then relabel and put back at the rear
	// of the queue and mark vertex as active.
	if(vertexList[v].excessFlow > 0){

		// DEBUG
    //    cout<< "relabeling vertex "<<v<<endl;
		relabel(adjList, vertexList, edgeList, v);
  	activeVertexQueue.push(v);
		vertexList[v].isActive= 1;
	}
	
	return;
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
	vector<vertex> vertexList;
	vector<edge> edgeList;
	vector<vector<int> > adjList;
	queue<int> activeVertexQueue;

	clock_t begin, end;
	double time_spent;
	begin = clock();

	// populate data structures from input text file, push all possible flow out of 
	// source vertex
	int ret= initialize(fileName, vertexList, edgeList, adjList, activeVertexQueue);
	if(ret == -1){
		return;
	}

	// apply discharge opeartion until active queue is empty
	while(activeVertexQueue.size() > 0){

		// pop next vertex from queue, mark as not active
		int v= activeVertexQueue.front();
		vertexList[v].isActive= 0;		
		activeVertexQueue.pop();

		//DEBUG
    //    cout<<"discharging: "<<v<<endl;

		// discharge vertex until excess becomes 0 or vertex is relabelled
		discharge(activeVertexQueue, vertexList, edgeList, adjList, v);
	}
	end = clock();
	cout<<vertexList[sinkId].excessFlow<<" is the maximum flow value"<<endl;
	printf("Elapsed %f seconds\n", (double) (end-begin)/CLOCKS_PER_SEC);
}

/*****************************************************************
 * Program starts here. Data structures are dynamically initialized
 * given an input text file. 
 *
 ****************************************************************/
int main(int argc, char** argv)
{

	if(argc < 2){
		cout << "usage: ./a.out <filename> " << endl;
	}
	
	string fileName(argv[1]);
	
	// read input text file containing graph with vertices, edges, and 
	// corresponding edges capacities
	//
	preflow_push(fileName.c_str());
  return(0);
}

