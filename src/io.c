#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "struct.h"
#include "io.h"

void SavePhotonListAscii(char *fname, lyman_RT_photons *Ph)
{
    FILE *f;
    int i;

    if(!(f = fopen(fname,"w")))
    {
	fprintf(stderr,"DumpPhotonList: Problem opening file %s\n",fname);
	MPI_Abort(MPI_COMM_WORLD, 0);

    }
    fprintf(f, "# %d\n", Ph->N_photons);
    for(i=0;i<Ph->N_photons;i++)
    {
	fprintf(f,"%e %e %e %e %e %e %e %d\n", Ph->Pos[3*i], Ph->Pos[3*i + 1], Ph->Pos[3*i + 2], 
		Ph->Dir[3*i], Ph->Dir[3*i + 1], Ph->Dir[3*i + 2], 
		Ph->x_out[i], Ph->Active[i]);      
    }
    fclose(f);
}


void OpenPhotonListAscii(char *fname, lyman_RT_photons *Ph)
{
    FILE *f;
    if(!(f = fopen(fname,"w")))
    {
	fprintf(stderr,"DumpPhotonList: Problem opening file %s\n",fname);
	MPI_Abort(MPI_COMM_WORLD, 0);

    }
    fprintf(f, "# %d %e %e %e %e %e\n", Ph->N_photons, All.Tau, All.Temperature, All. InputFrequency, All.TauDust, All.DustAbsorptionProb);
    fclose(f);
}

void AppendPhotonListAscii(char *fname, lyman_RT_photons *Ph, int index)
{
    FILE *f;
    int i;
    i = index;
    if(!(f = fopen(fname,"a")))
    {
	fprintf(stderr,"DumpPhotonList: Problem opening file %s\n",fname);
	MPI_Abort(MPI_COMM_WORLD, 0);

    }

    fprintf(f,"%e %e %e %e %e %e %e %d %d\n", Ph->Pos[3*i], Ph->Pos[3*i + 1], Ph->Pos[3*i + 2], 
	    Ph->Dir[3*i], Ph->Dir[3*i + 1], Ph->Dir[3*i + 2], 
	    Ph->x_out[i], Ph->Active[i], 
	    Ph->ScatterHI[i]);      
    fclose(f);
}
