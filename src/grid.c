#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "struct.h"
#include "grid.h"


float GridGetValue(grid *G, float x, float y, float z){
    int pos_cell[3];
    int index;
    float value;
    value = 0.0;
    if(G->cargo!=NULL){
	pos_cell[0] = (int)((x - G->x_level)/G->cell_size);
	pos_cell[1] = (int)((y - G->y_level)/G->cell_size);
	pos_cell[2] = (int)((z - G->z_level)/G->cell_size);
	
	/*sanity check*/
	if((pos_cell[0]>=0 && pos_cell[0]<G->N_grid_x) &&
	   (pos_cell[1]>=0 && pos_cell[1]<G->N_grid_y) &&
	   (pos_cell[2]>=0 && pos_cell[2]<G->N_grid_z)){
	    index =   pos_cell[2] + G->N_grid_z*(pos_cell[1] + G->N_grid_y*pos_cell[0]);	
	    value = G->cargo[index];
	}else{
	    fprintf(stderr, "GridGetValue: Problem with the boundaries\n");fflush(stdout);
	    fprintf(stderr, "x y z %f %f %f\n", x, y, z);
	    fprintf(stderr, "i j k %d %d %d \n", pos_cell[0], pos_cell[1], pos_cell[2]);
	    fprintf(stderr, "gridx gridy gridz %d %d %d \n", G->N_grid_x, G->N_grid_y, G->N_grid_z);
	    value = 0.0;
	}
    }
    
    return value;
}

float GridGetValueIndex(grid *G, int index){
    float value;
    value = 0.0;
    if(G->cargo!=NULL){
      value = G->cargo[index];
    }    
    return value;
}


void GridFree(grid *C){
  free(C->cargo);
  C->cargo = NULL;
  free(C);
}

void GridNullify(grid *C){
  int i;
  for(i=0;i<C->N_cells;i++){
    C->cargo[i] = 0.0;
  }
}

grid * GridCreate(void)
{
  grid  *C;

  if(!(C = malloc(sizeof(grid)))){
    fprintf(stderr, "CubeCreate: Problem with data allocation\n");fflush(stdout);
    exit(0);
  } 

  C->N_cells = 0;
  C->cell_size = 0.0;
  C->x_level  = 0.0;
  C->y_level  = 0.0;
  C->z_level  = 0.0;
  C->N_grid_x = 0;
  C->N_grid_y = 0;
  C->N_grid_z = 0;
  C->redshift = 0.0;
  C->cargo = NULL;
  return C;
}
