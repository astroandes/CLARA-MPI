#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "struct.h"
#include "lyman.h"
#include "RND_lyman.h"
#include "RT_ungrid.h"
#include "vector.h"
#include "propagation.h"


void RT_LymanUngrid(double *DirIn, double *DeltaPos, double *x_in, double a, double tau_HI, double n_HI){
    double PosIn[3];
    int n_scatter;

    PosIn[0] = 0.0;
    PosIn[1] = 0.0;
    PosIn[2] = 0.0;

    n_scatter = PropagatePackage(PosIn, DirIn, x_in, a, tau_HI, n_HI);
    
    DeltaPos[0] = PosIn[0]/All.UnitLength_in_cm;
    DeltaPos[1] = PosIn[1]/All.UnitLength_in_cm;
    DeltaPos[2] = PosIn[2]/All.UnitLength_in_cm;
}
