#include <stdlib.h>
#include <stdio.h>
 
#define NODES 20
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

capacities[9][16]=13609;
capacities[17][7]=53914;
capacities[16][13]=13938;
capacities[6][3]=75971;
capacities[9][18]=49675;
capacities[9][4]=32780;
capacities[11][13]=3453;
capacities[4][13]=63402;
capacities[14][13]=18707;
capacities[10][1]=2575;
capacities[18][12]=37443;
capacities[13][15]=85653;
capacities[8][17]=26807;
capacities[1][8]=34653;
capacities[1][10]=77117;
capacities[8][1]=70827;
capacities[2][7]=86388;
capacities[6][5]=41992;
capacities[4][2]=36246;
capacities[0][1]=99995;
capacities[0][5]=12822;
capacities[11][5]=72223;
capacities[14][11]=64740;
capacities[11][12]=16476;
capacities[16][19]=93378;
capacities[10][15]=45138;
capacities[18][18]=53917;
capacities[11][15]=26325;
capacities[16][2]=93713;
capacities[14][12]=3320;
capacities[17][19]=13454;
capacities[3][17]=15321;
capacities[16][14]=49249;
capacities[9][10]=75876;
capacities[6][8]=61920;
capacities[17][4]=25897;
capacities[9][10]=46919;
capacities[5][3]=76088;
capacities[2][6]=66781;
capacities[3][7]=76351;
capacities[0][15]=81337;
capacities[14][12]=94557;
capacities[15][6]=5284;
capacities[7][5]=16616;
capacities[3][8]=24716;
 
 
  printf("Capacity:\n");
  printMatrix(capacities);
 
  printf("Max Flow:\n%d\n", pushRelabel(capacities, flow, 0, 19));
 
  printf("Flows:\n");
  printMatrix(flow);
 
  return 0;
}
