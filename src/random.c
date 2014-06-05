#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "myrand.h"
#include "mtwist.h"
#include "struct.h"



void RandomDirection(float *direction)
     /*produces a random direction in space, homegeneously 
       sampling over the unit sphere*/
{
    double theta;
    double phi;
    theta = acos(2.0*RandFloatUnit()-0.5));      
    phi   = 2.0*PI*(RandFloatUnit());

    direction[0] = sin(theta)*cos(phi);
    direction[1] = sin(theta)*sin(phi);
    direction[2] = cos(theta);
}
