#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include "struct.h"
#include "parse.h"
#include "mtwist.h"
#include "lyman.h"
#include "propagation.h"

int main(int argc, char **argv)
{
  
  /* start MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisProc);
  MPI_Comm_size(MPI_COMM_WORLD, &SizeProc); 
  
  

  /*initialize Mersene Twister Random Generator*/
  mts_seed32(&(RND_MT_State), All.RandomSeed + ThisProc);


  /* read the simulation setup */
  ReadParameters(argv[1]);


  /*Make basic tests*/
  if(All.TestParallelVel || All.TestParallelVelFast || All.TestFirstScatter ||
      All.TestRND || All.TestPerpVel){
      TestAll();
  }

  /*Make science tests*/
  if(All.NeufeldSlab || All.NeufeldCube || All.ExpandingSphere || All.SimulationCube || All.RotatingSphere){
      PropagateAll();
  }


  /* finish MPI */
  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;    
}
