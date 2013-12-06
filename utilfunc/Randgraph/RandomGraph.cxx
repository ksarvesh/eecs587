
//
// RandomGraph:
//
//    Compile with: g++ RandomGraph.cxx UF.cxx -o RandomGraph
//

#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>

#include "UF.h"

using namespace std;

int gDistType;

struct Point {
  double x;
  double y;
  int id;
  Point *next;
};

double random_weight_uniform(double scale) {
  return(drand48() * scale);
}

// returns weight e^R where R is uniform in [0,E)
//
double random_weight_loguniform(double E) {
  double R = drand48() * E;
  return(exp(R));
}

// returns weight e^e^R where R is uniform in [0,E)
//
double random_weight_logloguniform(double E) {
  double R = drand48() * E;
  return(exp(exp(R)));
}

// mean should be in [0,1]
//
double random_weight_normal(double mean, int precision) {
  // to be implemented later ...
  return(random_weight_uniform(2*mean));
}

double random_weight() {  
  switch(gDistType) {
  case 0:
    return(1.0);
    break;
  case 1:
    return(random_weight_uniform(1.0));
    break;
  case 2:
    return(random_weight_loguniform(10.0));
    break;
  case 3:
    return(random_weight_loguniform(25.0));
    break;
  case 4:
    return(random_weight_loguniform(50.0));
    break;
  case 5:
    return(random_weight_loguniform(100.0));
    break;
  case 6:
    return(random_weight_logloguniform(4.0));
    break;
  case 7:
    return(random_weight_logloguniform(5.0));
    break;
  case 8:
    return(random_weight_normal(0.5,200));
    break;
  case 9:
    return(random_weight_normal(0.2,200));
    break;
  case 10:
    return(random_weight_normal(0.8,200));
  case 11:
    return(floor(random_weight_uniform(100000.0)));
    break;
  case 12:
    return(floor(random_weight_loguniform(11.0)));
    break;
  default:
    cout << "gDistType out of range\n";
    exit(0);
  }
}

void random_edge(int n, int* v1, int* v2) {
  *v1 = lrand48() % n;
  do {
    *v2 = lrand48() % n;
  } while(*v1 == *v2);
}






void swap(int* A, int i, int j) {
  int tmp = A[i];
  A[i]=A[j];
  A[j]=tmp;
}

void random_quicksort(int* A, int lo, int hi) {
  int i1,i2;

  if(lo >= hi) return;
  i1 = lo; // first unprocessed element.
  i2 = hi+1; // first element right of pivot.
  while(i1<i2) {
    if(lrand48() % 2) {
      i1++;
    } else {
      i2--;
      swap(A,i1,i2);
    }
  }
  random_quicksort(A, lo, i1-1);
  random_quicksort(A, i1, hi);
}

void random_d_regular(int n, int d, ofstream &out) {
  int *A = new int[n*d];
  int i;
  double weight;

  int source =0, sink = n-1;

  if(n*d % 2) { cout << "d or n must be even\n"; return; }
  for(i=0;i<n*d;i++)
    A[i] = i/d;
  random_quicksort(A,0,n*d-1);
  out << n << "\n";
  for(i=0;i<n*d;i+=2) {
    weight = random_weight();
    
    //Make it compatible to the network flow graph:
    if(A[i]==sink || A[i+1]==source || A[i]==A[i+1]){
        continue;
    }

    out << A[i] << " " << A[i+1] << " " << weight << "\n";
  }
  //out << "-1";
  delete[] A;
}

void random_bullseye(int n, int d, ofstream &out) {
  int orbits,i,j;
  int per_orbit, endpt1, endpt2;
  double base, scale, weight;

  cout << "How many orbits? ";
  cin >> orbits;
  if(orbits <= 0 || orbits > 100)
    {cout << "too many or too few\n"; exit(0);}
  
  out << n << endl;

  per_orbit = n/orbits;
  base = 1.0;
  scale = 1.0;
  for(i=0; i<orbits; i++) {

    // GENERATE per_orbit * d / 2 INTRA-ORBIT RANDOM EDGES

    for(j=0; j<per_orbit * d / 2; j++) {
      endpt1 = (lrand48() % per_orbit) + i * per_orbit;
      endpt2 = (lrand48() % per_orbit) + i * per_orbit;
      weight = base + scale * drand48();
      out << endpt1 << " " << endpt2 << " " << weight << endl;
    }

    // GENERATE SOME INTER-ORBIT EDGES
    if(i > 0)
      for(j=0; j<(per_orbit / 10)+2; j++) {
	endpt1 = (lrand48() % per_orbit) + i * per_orbit;    
	endpt2 = (lrand48() % per_orbit) + (i-1) * per_orbit;
	weight = base + scale * drand48();
	out << endpt1 << " " << endpt2 << " " << weight << endl;
      }

    base = base * 2;
    scale = scale * 2;
  }
  out << "-1\n";
}

// GENERATES edges EDGES BETWEEN RANDOM VERTICES IN [lo,hi]
// WITH RANDOM WEIGHTS IN [scale,2* scale]
void random_hier(int lo, int hi, int edges, double scale, ofstream &out) {
  int i, endpt1, endpt2;
  int range = (hi - lo + 1);
  double weight;

  for(i=0; i<edges; i++) {
    endpt1 = lrand48() % range + lo;
    endpt2 = lrand48() % range + lo;
    weight = scale + drand48() * scale;
    out << endpt1 << " " << endpt2 << " " << weight << endl;
  }
}

//
void random_hierarchical(int n, int d, ofstream &out) {
  int i,j;
  int initial_step_size = 20;
  int step;
  int hi,edges;
  double scale;
  double weight;
  
  out << n << endl;

  scale = 1.0;
  step = initial_step_size;

  for(j=0; j < n-1; j+= step) {
    for(i=0; i<step -1 && (i+j)<n-1; i++) {
      out << j + i << " "
	  << j + i + 1 << " "
	  << scale * drand48() << endl;
    }
  }

  scale *= 2;
  step *= d;
  i=2;
  while(1) {
    for(j=0; j < n; j += step) {
      hi = (j+step-1 < n) ? (j+step-1) : (n-1);
      edges = (hi-j+1)/i;
      random_hier(j,hi,edges,scale, out);            // add the random edges
    }
    i++;
    scale *= 2;          // up the level
    step *= d;           // and the step
    if(step > n) break;
  }

  random_hier(0,n-1,n/5,scale,out);
  out << "-1\n"; 
}

// RETURNS THE DISTANCE BETWEEN TO POINTS
double distance(Point *pt1, Point *pt2) {
  double xdiff = pt1->x - pt2->x;
  double ydiff = pt1->y - pt2->y;
  return(sqrt(xdiff * xdiff + ydiff * ydiff));
}

// throws points into the unit square, connects two points
// if they are within distance 1/d of each other.
void random_geometric(int n, int d, ofstream &out) {
  Point ***A = new Point**[d];
  Point *randPt, *tmp, *tmp2;
  int i,j,i1,j1,ret;
  double x,y,dist;

  // init a grid of (empty) lists of points
  for(i=0; i<d; i++) {
    A[i] = new Point*[d];
    for(j=0; j<d; j++)
      A[i][j] = 0;
  }

  out << n << "\n";

  // generate n points.   for each point put it in some A[i][j],
  // then find all points in A[i-1][j-1] .. A[i+1][j+1] that are
  // within distance d of the pt.

  for(n=n-1; n>=0; n--) {
    x = drand48();
    y = drand48();
    i = (int) (x * d);
    j = (int) (y * d);
    randPt = new Point;
    randPt->x = x;
    randPt->y = y;
    randPt->id = n;

    if(i==0) i1 = 0; else i1 = i-1;
    for( ; (i1 <= i+1) && (i1<d); i1++) {
      if (j==0) j1 = 0; else j1 = j-1;
      for( ; (j1<= j+1) && (j1<d); j1++) {
	tmp = A[i1][j1];
	while(tmp != 0) {
	  dist = distance(randPt, tmp);
	  if(dist < 1.0/d) {
	    // ADD A NEW EDGE FROM randPt TO tmp
	    out << randPt->id << " " << tmp->id << " " << dist << "\n";
	  }
	  tmp = tmp->next;
	}
      }
    }
    randPt->next = A[i][j];
    A[i][j] = randPt;	
  }

  out << "-1";

  // NOW FREE ALL THE MEMORY
  for(i=0; i<d; i++) 
    for(j=0; j<d; j++) {
      tmp = A[i][j];
      while(tmp != 0) { 
	tmp2 = tmp->next;
	delete tmp;
	tmp = tmp2;
      }
    }
}

// ONLY THE EDGE WEIGHTS ARE RANDOM IN THIS GRAPH.
// IT PRODUCES A "STRING", A GRAPH WITH n-1 EDGES AND DIAMETER n-1.
void random_string(int n, int excess, ofstream &out) {
  double weight;
  int i,v1,v2;
  
  out << n << "\n";
  for(i=0; i<n-1; i++) {
    weight = random_weight();
    out << i << " " << i+1 << " " << weight << "\n";
  }
  for(i=0; i<excess; i++) {
    random_edge(n,&v1,&v2);
    weight = random_weight();
    out << v1 << " " << v2 << " " << weight << endl;
  }

  out << "-1";
}

// ONLY EDGE WEIGHTS ARE RANDOM.
// TREE IS JUST A BINARY TREE
void random_binary_tree(int n, int excess, ofstream &out) {
  double weight;
  int i, v1, v2;
  
  out << n << "\n";
  for(i=2; i<n; i++) {
    weight = random_weight();
    out << i << " " << (i/2) << " " << weight << "\n";
  }
  weight = random_weight();
  out << 0 << " " << (n/2) << " " << weight << "\n";

  for(i=0; i<excess; i++) {
    random_edge(n, &v1, &v2);
    weight = random_weight();
    out << v1 << " " << v2 << " " << weight << endl;
  }
  out << "-1";
}

//
void random_tree(int n, int excess, ofstream &out) {
  int i=0, v1,v2;
  UnionFind UF;
  double weight;

  UF.init(n);
  
  out << n << "\n";

  // i counts the number of edges successfully added to the graph
  while(i<n-1) {
    random_edge(n, &v1, &v2);
    if(!UF.equiv(v1,v2)) {
      i++;
      weight = random_weight();
      out << v1 << " " << v2 << " " << weight << "\n";
    }
  }
  for(i=0; i<excess; i++) {
    random_edge(n, &v1, &v2);
    weight = random_weight();
    out << v1 << " " << v2 << " " << weight << endl;
  }
  out << "-1";
}

void random_connected_Gnm(int n, int density, ofstream &out) {
  int i, v1,v2;
  int excess = density*n/2 - (n-1);
  UnionFind UF;
  double weight;

  UF.init(n);
  
  out << n << "\n";

  // i counts the number of tree edges
  i=0;
  while(i<n-1) {
    random_edge(n, &v1, &v2);
    if(!UF.equiv(v1,v2)) {
      i++;
      weight = random_weight();
      out << v1 << " " << v2 << " " << weight << "\n";
    }
  }

  for(i=0; i< excess; i++) {
    random_edge(n, &v1, &v2);
    weight = random_weight();
    out << v1 << " " << v2 << " " << weight << "\n";
  }

  out << "-1";
}


void random_Gnp(int n, int density, ofstream &out) {

  int i,j;
  double weight;
  double p = ((double) density) / ((double) n); 

  out << n << "\n";    // write first line

  for(i=0;i<n;i++)
    for(j=i+1;j<n;j++)
      if(drand48() <= p) {     // generate (i,j) with prob. p
	weight = random_weight();
	out << i << " " << j << " " << weight << "\n";
      }
  out << "-1";
}

void random_Gnm(int n, int density, ofstream &out) {

  int i,v1,v2;
  int m = n * density / 2;
  double weight;

  out << n << "\n";    // write first line

  for(i=0;i<m;i++) {   // generate a bunch of random edges
    random_edge(n,&v1,&v2);
    weight = random_weight();
    out << v1 << " " << v2 << " " << weight << "\n";
  }
  out << "-1";
}

void random_d_ary(int n, int d, ofstream &out) {
  double weight;
  int i,j,k;
  
  out << n << endl;

  for(i=0; i<n; i++)                  // for each vertex
    for(j=0; j<d; j++) {              // pick d neighbors
      weight = random_weight();       // pick rand edge weight, add the edge
      k = lrand48() % n;
      out << i << " " << k << " " << weight << endl;
    }

  out << "-1" << endl;

}

void random_grid(int n, int type, ofstream &out) {

  int x,y;
  int i,j;

  if(type == 1)
    x = y = (int) floor(sqrt(n));
  else
    x = 16, y = n/16;

  out << (x*y) << endl;

  for(i=0; i<x; i++)
    for(j=0; j<y; j++) {
      if(i != x-1)
	out << i*x + j << " " << (i+1)*x + j << " " << random_weight() << endl;
      if(j != y-1)
	out << i*x + j << " " << i*x + j + 1 << " " << random_weight() << endl;
    }
  
  out << "-1\n";

}

void random_2dtorus(int n, int pct, ofstream &out) {
  double p = pct/100.0;  // probability of an edge being chosen
  int i,j;
  double weight;

  out << n*n << endl;

  for(i=0; i<n; i++) 
    for(j=0; j<n; j++) {
      if(drand48() <= p) {    // add edge going from (i,j) to (i+1,j)
	weight = random_weight();
	out << i*n + j << " " << ((i+1) % n)*n + j << " " << weight << endl;
      }
      if(drand48() <= p) {    // add edge going from (i,j) to (i,j+1)
	weight = random_weight();
	out << i*n + j << " " << i*n + ((j+1) % n) << " " << weight << endl;
      }
    }

  out << "-1" << endl;
}


void random_3dtorus(int n, int pct, ofstream &out) {
  double p = pct/100.0;  // probability of an edge being chosen
  int i,j,k;
  double weight;
  
  out << n*n*n << endl;

  for(i=0; i<n; i++) 
    for(j=0; j<n; j++)
      for(k=0; k<n; k++) {
	if(drand48() <= p) {    // add edge going from (i,j,k) to (i+1,j,k)
	  weight = random_weight();
	  out << i*n*n           + j*n + k << " " 
	      << ((i+1) % n)*n*n + j*n + k << " " 
	      << weight << endl;
	}
	if(drand48() <= p) {    // add edge going from (i,j,k) to (i,j+1,k)
	  weight = random_weight();
	  out << i*n*n + j*n           + k << " " 
	      << i*n*n + ((j+1) % n)*n + k << " " 
	      << weight << endl;
	}
	if(drand48() <= p) {    // add edge going from (i,j,k) to (i,j,k+1)
	  weight = random_weight();
	  out << i*n*n + j*n + k       << " " 
	      << i*n*n + j*n + (k+1)%n << " " 
	      << weight << endl;
	}
      }
  out << "-1" << endl;
}

int main(int argv, char** args) {
  
  char* graphfile;
  int graphsize;       // n - number of vertices
  int graphtype;       
  int disttype;        // weight distribution
  int graphdensity;    
  int seed;            // one may optionally supply the random seed

  // TYPE             INTERPRETATION OF "DENSITY"
  // ----             ---------------------------
  // 0 (d-regular)    all vertices have degree d (defaults to 3)
  // 1 (Gnm)          pick m=dn/2 random edges (defaults to 3)
  // 2 (Gnp)          choose each edge with prob p=d/n (defaults to 3)
  // 3 (geometric)    connect pts within dist 1/d (defaults to ~sqrt(n))
  // 4 (random tree)  n/a
  // 5 (string)       n/a
  // 6 (binary tree)  n/a
  // 7 (connected Gnm) generates random tree, adds dn/2 - (n-1) edges
  // 8 (d-ary)        each vertex chooses d random neighbors
  // 9 (2dtorus)      all edges in the n x n 2dtorus are chosen with probability d/100
  // 10 (3dtorus)     same.  n x n x n 3dtorus.

  
  

  ofstream out;
  
  if(argv < 2) {
    cout << "Arguments:\n"
	 << "name - name of output graph\n"
	 << "size - size of graph (default = 20)\n"
	 << "type - 0 - d-regular (default d = 3)\n"
	 << "     - 1 - Gnm (default m = dn/2)\n"
	 << "     - 2 - Gnp (default p = d/n)\n"
	 << "     - 3 - geometric (default dist = sqrt(n))\n"
	 << "     - 4 - random tree + d edges (default d=0)\n"
	 << "     - 5 - string + d edges\n"
	 << "     - 6 - binary tree + d edges\n"
	 << "     - 7 - connected Gnm (tree + dn/2 - (n-1) edges)\n"
	 << "     - 8 - undirected d-ary (default d=3)\n"
	 << "     - 9 - 2D torus (p = d/100, default d=40)\n"
	 << "     - 10 - 3D torus (p = d/100, default d=60)\n"
	 << "     - 11 - Bullseye (d = density)\n"
	 << "     - 12 - Hierarchical (log_d n levels, ~ (log_d n)/2 dense)\n"
	 << "     - 13 - Grid (d=1 is square, d=2 is 16 x n/16)\n"
	 << "dist - distribution of edge weights (default = 1)\n"
	 << "     - 0 - unit weight edges\n"
	 << "     - 1 - uniform in [0,1)\n"
	 << "     - 2 - = e^R, R uniform in [0,10)\n"
	 << "     - 3 - = e^R, R uniform in [0,25)\n"
	 << "     - 4 - = e^R, R uniform in [0,50)\n"
	 << "     - 5 - = e^R, R unifrom in [0,100)\n"
	 << "     - 6 - = e^{e^R}, R uniform in [0,4)\n"
         << "     - 7 - = e^{e^R}, R uniform in [0,5)\n"
	 << "     - 8 - normal dist, mean = M * .5\n"
	 << "     - 9 - normal dist, mean = M * .2\n"
	 << "     - 10 - normal dist, mean = M * .8\n"
	 << "     - 11 - INTEGER uniform in [0,100000]\n"
	 << "     - 12 - INTEGER log-uniform in [0,17]\n"
	 << "dens - graph density (interpretation depends on type)\n"
	 << "seed - random seed (default = current time)\n\n"
	 << "  e.g., RandomGraph grph 20 0 1 5\n\n"
	 << "will produce a 20 vertex, 5-regular graph named grph\n"
	 << "whose edge weights are uniformly distributed in [0,1)\n\n";
    return 0;
  }

  // COPY PARAMETERS
  graphfile = args[1];
  if(argv >= 3) graphsize = atoi(args[2]); else graphsize = 20;
  if(argv >= 4) graphtype = atoi(args[3]); else graphtype = 0;
  if(argv >= 5) disttype = atoi(args[4]); else disttype = 1;
  if(argv >= 6) graphdensity = atoi(args[5]); else graphdensity = 0;
  if(argv >= 7) seed = atoi(args[6]); else seed=(int) time(0);

  if(graphtype > 13 || disttype > 12) {
    cout << "Try zero arguments.  (if you need help).\n";
    return 0;
  }

  // RANDOM SEEDING 
  srand48(seed);

  // OPEN THE FILE
  out.open(graphfile);
    
  // SET THE GLOBAL VARIABLE INDICATING DISTRIBUTION TYPE
  gDistType = disttype;

  switch (graphtype) {
  case 0: 
    if(graphdensity==0) graphdensity=3;
    random_d_regular(graphsize,graphdensity,out); 
    break;
  case 1:
    if(graphdensity==0) graphdensity=3;
    random_Gnm(graphsize, graphdensity, out);
    break;
  case 2:
    if(graphdensity==0) graphdensity=3;
    random_Gnp(graphsize, graphdensity, out);
    break;
  case 3:
    if(graphdensity==0) graphdensity= (int) sqrt((double) graphsize);
    random_geometric(graphsize, graphdensity, out);
    break;
  case 4:
    random_tree(graphsize,graphdensity,out);
    break;
  case 5:
    random_string(graphsize,graphdensity,out);
    break;
  case 6:
    random_binary_tree(graphsize,graphdensity,out);
    break;
  case 7:
    if(graphdensity ==0) graphdensity = 3;
    random_connected_Gnm(graphsize, graphdensity, out);
    break;
  case 8:
    if(graphdensity ==0) graphdensity = 3;
    random_d_ary(graphsize, graphdensity, out);
    break;
  case 9:
    if(graphdensity == 0) graphdensity = 60;  // 60% 2dtorus edges chosen
    random_2dtorus(graphsize, graphdensity, out);
    break;    
  case 10:
    if(graphdensity == 0) graphdensity = 40;  // 40% 2dtorus edges chosen
    random_3dtorus(graphsize, graphdensity, out);
    break;
  case 11:
    if(graphdensity == 0) graphdensity = 3;
    random_bullseye(graphsize, graphdensity, out);
    break;
  case 12:
    if(graphdensity == 0) graphdensity = 8;  // actual density ~ 1 + ((log_d n) / 2)
    random_hierarchical(graphsize, graphdensity, out);
    break;
  case 13:
    if(graphdensity == 0) graphdensity = 1;  // 1 == square grid, 2 == 16 x n/16 grid
    random_grid(graphsize, graphdensity, out);
    break;
  default:
    cout << "graphtype unexpected.\n";
  }
 
  out.close();
}



