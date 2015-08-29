#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include "vector.h"
#include "struct.h"
#include "lyman.h"
#include "propagation.h"
#include "ray.h"
#include "RND_lyman.h"
#include "io.h"
#include "myrand.h"

int PropagateIsInside(double *Pos)
/*depending on the geometrical consideration of the problem at hand,
  decides if the photon is still inside*/
{
    int is_in;
    double radius;
    is_in = 1;

    if(All.SimulationCube){
	if(fabs(Pos[0])<All.SlabLength && 
	   fabs(Pos[1])<All.SlabLength &&
	   Pos[2] < (All.SlabLength*2.0) && Pos[2]>0.0){
	    is_in = 1;
	}else{
	    is_in = 0 ;
	}
    }

    if(All.NeufeldSlab){
	if(fabs(Pos[0])<All.SlabLength){
	    is_in = 1;
	}else{
	    is_in = 0 ;
	}
    }

    if(All.NeufeldCube){
	if(fabs(Pos[0])<All.SlabLength && 
	   fabs(Pos[1])<All.SlabLength &&
	   fabs(Pos[2])<All.SlabLength){
	    is_in = 1;
	}else{
	    is_in = 0 ;
	}
    }

    if(All.ExpandingSphere||All.RotatingSphere||All.OutflowRotationSphere){
	radius = Pos[0]*Pos[0] + Pos[1]*Pos[1] + Pos[2]*Pos[2];
	radius = sqrt(radius);
	if(radius< All.SlabLength){
	    is_in = 1;
	}else{
	    is_in = 0 ;
	}
    }


    
    return is_in;
}

void PropagateAllSetup(void)
/*makes the general problem setup from the input values*/
{
    double a, nu_doppler, ly_sigma_0, column_HI;
    
    nu_doppler = CONSTANT_NU_DOPPLER*sqrt(All.Temperature/10000.0); /* in s^-1 */
    a = Lya_nu_line_width_CGS/(2.0*nu_doppler);
    ly_sigma_0 =  LyaCrossSection(0,a);
    
    All.v_thermal = (nu_doppler/Lya_nu_center_CGS)*C_LIGHT;/*In cm/s*/

    column_HI = All.Tau / ly_sigma_0;
    All.SlabLength = All.Tau/(All.NumberDensityHI*ly_sigma_0);
    
#ifdef DEBUG
    fprintf(stdout, "ColumnDensinty %e\n", column_HI);
    fprintf(stdout, "Temp %g\n", All.Temperature);
    fprintf(stdout, "a %g\n", a);
    fprintf(stdout, "ly_sigma_0 %g\n", ly_sigma_0);
    fprintf(stdout, "Slab Length [cm] %g\n", All.SlabLength);
    fprintf(stdout, "maximum at (a*tau)^1/3 %g\n", pow(a*All.Tau,0.333));
    fprintf(stdout, "scatterings %g\n", 1.612*All.Tau);
#endif

}

int PropagateStep(double *x, double *k_in_photon, double *r_travel, double a, double n_HI)
/* Calculates: 
   1. new frequency 
   2. new propagation direction, 
   3. distance to travel to next interaction

   Note:
   At this point the frequency must be already the comoving one (not lab frame!)
   If the photon is absorbed by dust, then the status is set to -1;
*/

{
    double k_out_photon[3];
    double u_atom[3];
    double lya_sigma;
    double tau;
    double x_in, x_out;
    double r;
    double N;
    double a_in;
    int i;
    int status;
    double HIInteractionProb; /*probability of interacting with HI*/
    double dust_sigma;
    double rand_interaction;
    double rand_absorption;
    double v_parallel;;
    double nu_doppler;
    double v_thermal;
    double g_recoil;
    double temperature;
    

    nu_doppler = Lya_nu_line_width_CGS/(2.0*a);
    v_thermal = (nu_doppler/Lya_nu_center_CGS)*C_LIGHT;/*In cm/s*/
    temperature = (nu_doppler/CONSTANT_NU_DOPPLER)*(nu_doppler/CONSTANT_NU_DOPPLER)*10000.0;
    HIInteractionProb = 1.0;
    rand_absorption = 0.0;
    /*get a random number to decide if the thing interact with HI*/
    rand_interaction  = RandFloatUnit();

    g_recoil = PLANCK*nu_doppler/(2.0*BOLTZMANN*temperature);

    status = ACTIVE;
    /* some initializations */
    x_in = *x;
    N    = n_HI;
    a_in = a;


    /* basic status at the interacion point*/
    lya_sigma  = LyaCrossSection(x_in, a_in);
    tau        = LyaTau();
    r          = tau/(lya_sigma*N); /*the distance is in cm*/

    if(All.UseDust){
	dust_sigma = PI*All.GrainSize*All.GrainSize;
	All.NumberDensityDust = All.TauDust/(All.SlabLength*dust_sigma);
	/*probability of interacting with HI*/
	HIInteractionProb = 
	  All.NumberDensityHI*lya_sigma/(All.NumberDensityHI*lya_sigma + All.NumberDensityDust*dust_sigma);
/*	fprintf(stdout, "Interaction prob %e \n", HIInteractionProb);*/
	/*get a random number to decide if the thing will be absorbed*/
	rand_absorption = RandFloatUnit();
	
	/* The traveled distance is modifyed by the presence of dust*/
	r          = tau/(lya_sigma*N + dust_sigma*All.NumberDensityDust); /*the distance is in cm*/
    }

    if(rand_interaction <= HIInteractionProb){

	/* generate the atom velocity the comoving frame*/
	RND_lyman_atom(&(u_atom[0]), &(k_in_photon[0]), 
		       &(k_out_photon[0]), x_in, a_in);


	/*find the new frequency (in the observer system)*/
	x_out = x_in - point_product(u_atom, k_in_photon) +
	  point_product(u_atom, k_out_photon) +	  
	  g_recoil * (point_product(k_in_photon, k_out_photon) - 1.0);

    }else{
	x_out = x_in;
    }

    for(i=0;i<3;i++){
	k_in_photon[i] = k_out_photon[i];
    }

    *x = x_out;
    *r_travel = r ;

    if((rand_interaction > HIInteractionProb) && All.UseDust && (rand_absorption < All.DustAbsorptionProb)){
	*r_travel = 0.0;
	status = ABSORBED;
    }
    
    return status;
}

void PropagateGetBulkVel(double *BulkVel, double *Pos){
    int i;

    for(i=0;i<3;i++){
	BulkVel[i] = 0.0;
    }

    if(All.ExpandingSphere){
	for(i=0;i<3;i++){
	    BulkVel[i] = (Pos[i]/All.SlabLength)*(All.VmaxSphere*1.0e5); /*in cm/s*/
	}
	
    }
    if(All.RotatingSphere){
      BulkVel[0]=(-Pos[1]/All.SlabLength)*All.VmaxSphere*1.0e5;
      BulkVel[1]=(Pos[0]/All.SlabLength)*All.VmaxSphere*1.0e5;
      BulkVel[2]=0.0;
    }	
    if(All.OutflowRotationSphere){
      /*first the outflow component*/
      	for(i=0;i<3;i++){
	    BulkVel[i] = (Pos[i]/All.SlabLength)*(All.VmaxOutflow*1.0e5); /*in cm/s*/	      
	}
	/*then the rotating component*/
	BulkVel[0] = BulkVel[0] + (-Pos[1]/All.SlabLength)*All.VmaxRotation*1.0e5;
	BulkVel[1] = BulkVel[1] + (Pos[0]/All.SlabLength)*All.VmaxRotation*1.0e5;
	BulkVel[2] = BulkVel[2] + 0.0;
    }
    
}


void PropagateGetTemperature(double *Temp, double *Pos){
    *Temp = All.Temperature;
}

void PropagateGetNumberDensity(double *n_HI, double *Pos){
    *n_HI = All.NumberDensityHI;
}


void PropagateLorentzDirChange(double *Dir, double *BulkVel, int sign){
    int i;
    for(i=0;i<3;i++){
	Dir[i] = Dir[i]*(1.0 + sign*(BulkVel[i]/C_LIGHT)); 
    }
    
    /*renormalize, I just want to know the direction change, bust still be normalized*/
    normalize(Dir);
}

void PropagateLorentzFreqChange(double *x, double *Dir, 
			    double *BulkVel, double ThermalVel, int sign){
    double lorentz_factor;
    lorentz_factor = (BulkVel[0]*Dir[0] + BulkVel[1]*Dir[1] + BulkVel[2]*Dir[2])/ThermalVel;
    
    *x  = *x + sign*lorentz_factor;
}

int PropagatePackage(double *PosIn, double *DirIn, double *x_in, int *status){
    int i;
    double Pos[3];
    double Dir[3];
    double r_travel, x, x_comoving;
    double a, nu_doppler, n_HI, v_thermal, BulkVel[3], temperature;
    int n_iter;
    int stat;
    float last_x;
    n_iter=0;

    /*Make the initialization*/
    for(i=0;i<3;i++){
	Pos[i] = PosIn[i];
	Dir[i] = DirIn[i];
    }
    x = *x_in;
    

#ifdef DEEP_DEBUG
    fprintf(stdout, "starting propagation (x=%f)\n", *x_in);
      fflush(stdout);
#endif
    stat = *status;
    /*difuse the photon in space and frequency until it gets out*/
    while(PropagateIsInside(&(Pos[0]))&&(stat==ACTIVE)&&n_iter<MAX_ITER){
      /* get the temperature at this point*/
      PropagateGetTemperature(&temperature, Pos);
      
      /* get the number density at this point*/
      PropagateGetNumberDensity(&n_HI, Pos);
      
      /*get the bulk velocity of the fluid at this point*/
      PropagateGetBulkVel(BulkVel, Pos);
      
      /*Get the thermal velocity and doppler broadening*/
      nu_doppler = CONSTANT_NU_DOPPLER*sqrt(temperature/10000.0); /* in cm/s */
      a = Lya_nu_line_width_CGS/(2.0*nu_doppler);
      v_thermal = (nu_doppler/Lya_nu_center_CGS)*C_LIGHT;/*In cm/s*/
      
      /*change the value of the frequency to one comoving with the fluid*/
      PropagateLorentzFreqChange(&x, Dir, BulkVel, v_thermal, -1); 
            
      /*change the direction of the photon to the fluid frame*/
      PropagateLorentzDirChange(&(Dir[0]), BulkVel, -1);
      


      /*--------------------------------------------------------------------------*/
      /*Change the frequency and the Propagation direction, find the displacement*/	
      stat = PropagateStep(&x, Dir, &r_travel, a, n_HI);	    	
      /*--------------------------------------------------------------------------*/
      


      /*Change the new direction to the lab frame value*/
      PropagateLorentzDirChange(Dir, BulkVel, 1);
      
      /*Change the frequency comoving to the lab frame value*/
      PropagateLorentzFreqChange(&x, Dir, BulkVel, v_thermal, 1); 
      
      
      /*Update the position*/
      for(i=0;i<3;i++){
	Pos[i] += r_travel*Dir[i];	    
      }
      
      n_iter++;
      last_x = x;
    }
#ifdef DEEP_DEBUG
      fprintf(stdout, "finishing propagation (x=%f)\n", x);
      fflush(stdout);
#endif
    
    /*update the final status of the photon, just to know if it was absorbed, or what*/
    if(stat==ACTIVE){
      stat = OUT_OF_BOX;
    }
    if(n_iter>= MAX_ITER){
      stat = SATURATED_ITERATIONS;
    }
    
    /*Update the final values*/
    for(i=0;i<3;i++){
	PosIn[i] = Pos[i];
	DirIn[i] = Dir[i];
    }

    *x_in = x;
    *status = stat;
    return n_iter;
}
    
void PropagateAll(void)
{
    int n_packages;
    int i;
    lyman_RT_photons *Ph;
    int n_iter;
    char FileName[MAX_FILENAME_SIZE];


    /*update some geometrical values*/
    PropagateAllSetup();
    
    /*get the number of photons to propagate*/
    n_packages = (int)((All.TotalLuminosity/SizeProc)/All.LuminosityPerPackage);
#ifdef DEBUG
    fprintf(stdout, "%d packages to propagate\n", n_packages);
#endif 

    /*create the packages*/
    Ph = PhotonListCreate(n_packages);

    
    /*initialize the packages*/
    PhotonListInitialize(Ph);
    MPI_Barrier (MPI_COMM_WORLD);
    
    if(All.OutputInitList){
	sprintf(FileName, "%s/%s_in.proc.%d.ascii", All.OutputDir, All.OutputFile, ThisProc);
	SavePhotonListAscii(FileName, Ph);
    }
    
    if(All.OutputFinalList){
	sprintf(FileName, "%s/%s_out.proc.%d.ascii", All.OutputDir, All.OutputFile, ThisProc);
	OpenPhotonListAscii(FileName, Ph);
    }
    
    /*propagate each package*/
    for(i=0;i<n_packages;i++){
#ifdef DEEP_DEBUG

      fprintf(stdout, "Propagating packages %d/%d\n pro(%d/%d)\n", i, n_packages, ThisProc, SizeProc);
      fflush(stdout);
#endif

	Ph->ScatterHI[i] = 
	    PropagatePackage(&(Ph->Pos[3*i]), &(Ph->Dir[3*i]), &(Ph->x_out[i]), &(Ph->Active[i]));	
	if(All.OutputFinalList){
	    sprintf(FileName, "%s/%s_out.proc.%d.ascii", All.OutputDir, All.OutputFile, ThisProc);
	    AppendPhotonListAscii(FileName, Ph, i);    
	}    
    }
    fprintf(stdout, "finished writing (prc %d)\n", ThisProc);
    /*free the memory*/
    PhotonFree(Ph);
}

void TestFirstScatter(double x, double a)
{
    double x_in;
    double r_travel;
    double k_in_photon[3];
    double n_HI;
    int status;
    int i;
    char FileName[MAX_FILENAME_SIZE];
    FILE *out;

    status = 0;
    PropagateAllSetup();


    sprintf(FileName, "%s/%s_FirstScatter.proc.%d", All.OutputDir, All.OutputTestFile, ThisProc);
    if(!(out=fopen(FileName, "w"))){
	fprintf(stderr, "TestFirstScatter problem opening file %s\n", FileName);
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD, 0);
    }
    fprintf(out, "%d %e %e\n", N_POINTS_IN_TEST, x, a);
    for(i=0;i<N_POINTS_IN_TEST;i++){
	x_in = x;
	n_HI = All.NumberDensityHI;
	RND_spherical(&(k_in_photon[0]));    
	status = PropagateStep(&x_in, &(k_in_photon[0]), &r_travel, a, n_HI);
	fprintf(out, "%e\n", x_in);
    }
    fclose(out);
} 

void TestRND(void)
{
    double x_in;
    double r_travel;
    double k_in_photon[3];
    double n_HI;
    int status;
    int i;
    char FileName[MAX_FILENAME_SIZE];
    FILE *out;

    sprintf(FileName, "%s/%s_RND.proc.%d", All.OutputDir, All.OutputTestFile, ThisProc);
    if(!(out=fopen(FileName, "w"))){
	fprintf(stderr, "TestRND problem opening file %s\n", FileName);
	MPI_Abort(MPI_COMM_WORLD, 0);
    }

    fprintf(out, "%d %e %e\n", N_POINTS_IN_TEST, 0.0, 0.0);
    for(i=0;i<N_POINTS_IN_TEST;i++){
	r_travel = RandFloatUnit();
	fprintf(out, "%e\n", r_travel);
    }

    fclose(out);
} 

void TestAll(void){
    if(All.TestParallelVel){
	RND_lyman_parallel_vel_test(All.Test_x,All.Test_a);
    }

    if(All.TestParallelVelFast){
	RND_lyman_parallel_vel_fast_test(All.Test_x,All.Test_a);
    }    

    if(All.TestFirstScatter){
	TestFirstScatter(All.Test_x, All.Test_a);
    }

    if(All.TestRND){
	TestRND();
    }

    if(All.TestPerpVel){
       RND_lyman_perp_vel_test();
    }
}
