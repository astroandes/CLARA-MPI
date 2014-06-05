#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "struct.h"
#include "io.h"
#include "grid.h"
/*
  xlc -q64 mtwist.c endian.c io.c grid.c SFR.c
*/
int main(int argc, char **argv){
  char FileName[MAX_FILENAME_SIZE];
  FILE *out;
  float total;
  int i, j;

  grid *G;
  i = 10;
  if(!(out=fopen("SFR_list.dat", "w"))){
    fprintf(stdout, "problem opening file\n");
    exit(0);
  }

  for(i=10;i<=4387;i++){
    sprintf(FileName, "/homea/hpo08/hpo089/LYMAN_MN_50/RAW_CUBES/grid_object_%d_GasSFR.dat", i);    
    G = GridCreate();
    G->cargo = LoadCubeField(FileName, G, "GasSFR");
    total = 0.0;
    for(j=0;j<G->N_cells;j++){
      total+=G->cargo[j];    
    }
    total *= G->cell_size*G->cell_size*G->cell_size;
    fprintf(stdout, "%d %e\n", i, total);
    fprintf(out, "%d %e\n", i, total);
    GridFree(G);
  }

  fclose(out);

}
