#include <stdlib.h>
#include <stdio.h>
#include <time.h>
 
#define NODES 1000
#define MIN(X,Y) X < Y ? X : Y
#define INFINITE 10000000
 
void push(const int * const * C, int ** F, int *excess, int u, int v) {
  int send = MIN(excess[u], C[u][v] - F[u][v]);
  F[u][v] += send;
  F[v][u] -= send;
  excess[u] -= send;
  excess[v] += send;
}
 
void relabel(const int * const * C, const int * const * F, int *height, int u) {
  int v;
  int min_height = INFINITE;
  for (v = 0; v < NODES; v++) {
    if (C[u][v] - F[u][v] > 0) {
      min_height = MIN(min_height, height[v]);
      height[u] = min_height + 1;
    }
  }
}
 
void discharge(const int * const * C, int ** F, int *excess, int *height, int *seen, int u) {
  while (excess[u] > 0) {
    if (seen[u] < NODES) {
      int v = seen[u];
      if ((C[u][v] - F[u][v] > 0) && (height[u] > height[v])){
        push(C, F, excess, u, v);
      }
      else
        seen[u] += 1;
    } else {
      relabel(C, F, height, u);
      seen[u] = 0;
    }
  }
}
 
void moveToFront(int i, int *A) {
  int temp = A[i];
  int n;
  for (n = i; n > 0; n--){
    A[n] = A[n-1];
  }
  A[0] = temp;
}
 
int pushRelabel(const int * const * C, int ** F, int source, int sink) {
  int *excess, *height, *list, *seen, i, p;
 
  excess = (int *) calloc(NODES, sizeof(int));
  height = (int *) calloc(NODES, sizeof(int));
  seen = (int *) calloc(NODES, sizeof(int));
 
  list = (int *) calloc((NODES-2), sizeof(int));
 
  for (i = 0, p = 0; i < NODES; i++){
    if((i != source) && (i != sink)) {
      list[p] = i;
      p++;
    }
  }
 
  height[source] = NODES;
  excess[source] = INFINITE;
  for (i = 0; i < NODES; i++)
    push(C, F, excess, source, i);
 
  p = 0;
  while (p < NODES - 2) {
    int u = list[p];
    int old_height = height[u];
    discharge(C, F, excess, height, seen, u);
    if (height[u] > old_height) {
      moveToFront(p,list);
      p=0;
    }
    else
      p += 1;
  }
  int maxflow = 0;
  for (i = 0; i < NODES; i++)
    maxflow += F[source][i];
 
  free(list);
 
  free(seen);
  free(height);
  free(excess);
 
  return maxflow;
}
 
void printMatrix(const int * const * M) {
  int i,j;
  for (i = 0; i < NODES; i++) {
    for (j = 0; j < NODES; j++)
      printf("%d\t",M[i][j]);
    printf("\n");
  }
}
 
int main(void) {
  int **flow, **capacities, i;
  flow = (int **) calloc(NODES, sizeof(int*));
  capacities = (int **) calloc(NODES, sizeof(int*));
  for (i = 0; i < NODES; i++) {
    flow[i] = (int *) calloc(NODES, sizeof(int));
    capacities[i] = (int *) calloc(NODES, sizeof(int));
  }



//  printf("Capacity:\n");
//  printMatrix(capacities);
	clock_t begin, end;
	double time_spent;
	begin = clock(); 
  int maxFlow = pushRelabel(capacities, flow, 0, 999);
 	end = clock();
	printf("Elapsed: %f seconds\n", (double) (end-begin) / CLOCKS_PER_SEC);
  printf("Max Flow:\n%d\n", maxFlow);

//  printf("Flows:\n");
//  printMatrix(flow);
 
  return 0;
}
