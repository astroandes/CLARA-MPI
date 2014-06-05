#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include "struct.h"
#include "grid.h"
#include "io.h"
#include "RT_grid.h"
#include "RND_lyman.h"
#include "lyman.h"
#include "ray.h"


void RT_LoadStars(lyman_RT_stars *S){
    char FileName[MAX_FILENAME_SIZE];
    char FieldName[30];
    FILE *in;
    int i;
    float dumm;

    strcpy(FieldName, "Stars");
    sprintf(FileName , "%s/%s_%s.dat", All.InputDir, All.CubeName, FieldName);

    if(!(in=fopen(FileName,"r"))){
	fprintf(stderr, "Problem opening file %s\n", FileName);
	exit(1);
    }


    fscanf(in, "%d\n", &(S->N_stars));
    fprintf(stdout, "Stars to read %d\n", S->N_stars);

    if(!(S->pos_x = malloc(S->N_stars*sizeof(float)))){
	fprintf(stderr, "Problem in the allocation of star pos_x\n");
	exit(1);
    }

    if(!(S->pos_y = malloc(S->N_stars*sizeof(float)))){
	fprintf(stderr, "Problem in the allocation of star pos_y\n");
	exit(1);
    }

    if(!(S->pos_z = malloc(S->N_stars*sizeof(float)))){
	fprintf(stderr, "Problem in the allocation of star pos_z\n");
	exit(1);
    }

    if(!(S->LyaLum = malloc(S->N_stars*sizeof(double)))){
	fprintf(stderr, "Problem in the allocation of star LyaLum\n");
	exit(1);
    }

    
    for(i=0;i<S->N_stars;i++){
	fscanf(in, "%e %e %e %e %e %e %lf\n", &(S->pos_x[i]), &(S->pos_y[i]), &(S->pos_z[i]),
	       &dumm, &dumm, &dumm, &(S->LyaLum[i]));
	fprintf(stdout, "Leido %e %e %e %lf \n", (S->pos_x[i]), (S->pos_y[i]), (S->pos_z[i]),
	       (S->LyaLum[i]));
    }

    fclose(in);    
}

void RT_LoadGrids(grid * GasTauDust, grid *GasTauHI, grid *GasTemp, grid *GasVelX, grid *GasVelY, grid *GasVelZ){
    char FileName[MAX_FILENAME_SIZE];
    char FieldName[30];    
    float min, max;
    /*read the input information*/
    strcpy(FieldName, "GasTauHI");
    sprintf(FileName , "%s/%s_%s.dat", All.InputDir, All.CubeName, FieldName);
    GasTauHI->cargo = LoadCubeField(FileName, GasTauHI, FieldName);
    print_min_max(GasTauHI->cargo, GasTauHI->N_cells, 1, "TauHI", &min, &max);    

    strcpy(FieldName, "GasTemp");
    sprintf(FileName , "%s/%s_%s.dat", All.InputDir, All.CubeName, FieldName);
    GasTemp->cargo = LoadCubeField(FileName, GasTemp, FieldName);
    
    print_min_max(GasTemp->cargo, GasTemp->N_cells, 1, "Temp", &min, &max);

    
    if(All.UseDust){
	strcpy(FieldName, "GasTauDust");
	sprintf(FileName , "%s/%s_%s.dat", All.InputDir, All.CubeName, FieldName);
	GasTauDust->cargo = LoadCubeField(FileName, GasTauDust, FieldName);
    }
    
    if(All.UseVelocities){
	strcpy(FieldName, "GasVelX");
	sprintf(FileName , "%s/%s_%s.dat", All.InputDir, All.CubeName, FieldName);
	GasVelX->cargo = LoadCubeField(FileName, GasVelX, FieldName);
	
	strcpy(FieldName, "GasVelY");
	sprintf(FileName , "%s/%s_%s.dat", All.InputDir, All.CubeName, FieldName);
	GasVelY->cargo = LoadCubeField(FileName, GasVelY, FieldName);
	
	strcpy(FieldName, "GasVelZ");
	sprintf(FileName , "%s/%s_%s.dat", All.InputDir, All.CubeName, FieldName);
	GasVelZ->cargo = LoadCubeField(FileName, GasVelZ, FieldName);
    }
    
    /*Make some consistency checks*/
    if(GasTemp->N_cells!=GasTauHI->N_cells){
	fprintf(stderr, "inconsistent array dimension\n");
	MPI_Abort(MPI_COMM_WORLD, 0);
    }
    
    if(All.UseDust){
	if(GasTemp->N_cells!=GasTauDust->N_cells){
	    fprintf(stderr, "inconsistent array dimension\n");
	    MPI_Abort(MPI_COMM_WORLD, 0);
	}
    }
    
    if(All.UseVelocities){
	if(GasTemp->N_cells!=GasVelX->N_cells || GasVelY->N_cells!=GasVelZ->N_cells || 
	   GasVelZ->N_cells!=GasVelX->N_cells){
	    fprintf(stderr, "inconsistent array dimension\n");
	    MPI_Abort(MPI_COMM_WORLD, 0);
	}
    }
    
    
#ifdef DEBUG
    fprintf(stdout, "Finished loading basic inputs\n");fflush(stdout);
#endif
    
}

int PhotonNumber(lyman_RT_stars *S)
{
    int n_total;
    int i;
    int n_per_star;
    n_total = 0;
    for(i=0;i<S->N_stars;i++){
	n_per_star = (int)(S->LyaLum[i]/(All.LuminosityPerPackage*SizeProc));
	if(n_per_star<1){
	    n_per_star = MIN_N_PHOT/SizeProc;
	    //	    fprintf(stderr, "WARNING the Luminosity PerPackage is small\n");
	}
	n_total +=  n_per_star;
    }
  
#ifdef DEBUG
  fprintf(stdout, "The total number of packages in each proc %d\n",n_total);fflush(stdout);
#endif
  return n_total; 
}



 void  InitializePhotonList(lyman_RT_stars *S, grid *GasVelX, grid *GasVelY, grid *GasVelZ, grid *GasTemp, 
			   lyman_RT_photons *PhList, int n_procs){
  int i, j, k, l;
  int n_here;
  int i_photon;
  int index;
  float direction[3];
  float nu_doppler, v_thermal;
  int n_photons_per_star;
  
  i_photon=0;
  index = 0;
  fprintf(stdout, "npro %d %e\n", n_procs, All.LuminosityPerPackage);

  for(i=0;i<S->N_stars;i++){
      n_photons_per_star = (int)(S->LyaLum[i]/(All.LuminosityPerPackage*SizeProc));
      fprintf(stdout, "nphotons per star %d %e %e\n", n_photons_per_star, S->LyaLum[i], All.LuminosityPerPackage*SizeProc);
      if(n_photons_per_star<1){
	  n_photons_per_star = MIN_N_PHOT/SizeProc;
      }

      for(j=0;j<n_photons_per_star;j++){
	PhList->Lum[i_photon] = S->LyaLum[i]/(1.0*SizeProc*n_photons_per_star);


	PhList->Pos[3*i_photon + 0] = S->pos_x[i];
	PhList->Pos[3*i_photon + 1] = S->pos_y[i];
	PhList->Pos[3*i_photon + 2] = S->pos_z[i];
	
	RND_sphericalf(direction);
	
	PhList->Dir[3*i_photon + 0] = direction[0];
	PhList->Dir[3*i_photon + 1] = direction[1];
	PhList->Dir[3*i_photon + 2] = direction[2];
	
	PhList->Cell[3*i_photon + 0] = (int)((S->pos_x[i] - GasTemp->x_level)/GasTemp->cell_size);
	PhList->Cell[3*i_photon + 1] = (int)((S->pos_y[i] - GasTemp->y_level)/GasTemp->cell_size);
	PhList->Cell[3*i_photon + 2] = (int)((S->pos_z[i] - GasTemp->z_level)/GasTemp->cell_size);
	
	PhList->x_out[i_photon] = 0.0;	      
	
	i_photon++;
      }
  }


  
  /*Minor sanity check*/
  if(i_photon>PhList->N_photons){
      fprintf(stderr, "The photon list is full\n");fflush(stdout);
      MPI_Abort(MPI_COMM_WORLD, 0);    
  }
  PhList->N_photons = i_photon;
  
#ifdef DEBUG
  fprintf(stdout, "Finished creation and initialization of photon list (xProc %d/%d)(%d photons)\n", 
	  ThisProc, SizeProc, PhList->N_photons);fflush(stdout);
#endif
 }

 
	    
	    

	    
	    
	    /*Initialize the frequencies*/
/*
  if(All.UseVelocities){
  nu_doppler = CONSTANT_NU_DOPPLER*sqrt(GasTemp->cargo[index]/10000.0); //in cm/s 
  v_thermal = (nu_doppler/Lya_nu_center_CGS)*C_LIGHT; //In cm/s
  
  
  PhList->x_out[i_photon] = 
  direction[0]*GasVelX->cargo[index] +
  direction[1]*GasVelY->cargo[index] +
  direction[2]*GasVelZ->cargo[index] ;  
  PhList->x_out[i_photon] *= 1.0/(v_thermal);
  }
*/

int InsideBox(float x, float y, float z, grid *G)
{
 
    float min_x, min_y, min_z;
    float max_x, max_y, max_z;
    int is_in;
    int i, j, k;
    min_x = G->x_level;
    max_x = min_x  + (G->N_grid_x-1)  * G->cell_size;

    min_y = G->y_level;
    max_y = min_y  + (G->N_grid_y-1) * G->cell_size;

    min_z = G->z_level;
    max_z = min_z  + (G->N_grid_z-1) * G->cell_size;

    i = (int)((x - G->x_level)/G->cell_size);
    j = (int)((y - G->y_level)/G->cell_size);
    k = (int)((z - G->z_level)/G->cell_size);
    is_in = k + G->N_grid_z * (j + G->N_grid_y * i);


    if((x>min_x)&&(y>min_y)&&(z>min_z)&&(x<max_x)&&(y<max_y)&&(z<max_z)){	
	is_in = k + G->N_grid_z * (j + G->N_grid_y * i);
    }else{
	is_in = -1;
    }    
    return is_in;
}

float EscapeFraction(float a, float tau_HI, float effective_theta, float tau_dust){
  float escape;
  float xsi;
  xsi = 0.525;/*fitting factor from verhamme, schaerer and maselli*/
  escape = pow((a*tau_HI), (1.0/3.0))*tau_dust;
  escape = (sqrt(3)/(pow(PI,5.0/12.0)*xsi))*sqrt(escape);
  escape = 1.0/cosh(escape);
  return escape;
}




void RT_LymanSimple(double *Dir, double *DeltaPos, double *x, double a, double tau, double effective_eta){
    double x_in;
    double x_out;
    double tmp1, tmp2, tmp3;
    double sign;
    x_in = *x;



    tmp1 = mts_drand(&RND_MT_State);
    while(tmp1<1e-18)
    {
	tmp1 = mts_drand(&RND_MT_State);
    }

    tmp2 = (double)sqrt(54.0/pow(PI,3.0))* effective_eta * a * tau * log(tan(PI*tmp1/2.0));    
    tmp3   = x_in*x_in*x_in;
    if((tmp2+tmp3)>0){
	sign = 1;
    }else{
	sign = -1;
    }

    x_out = sign*pow(fabs(tmp2 + tmp3), (1.0/3.0));
    *x = x_out;
}


float RT_LocalTau(float x_value, float tau_HI, float temp, 
		  float dir_x, float dir_y, float dir_z,
		  float vel_x, float vel_y, float vel_z, 
		  float delta_step, float cell_size){
    double nu_doppler, a, thermal_vel;    
    double lya_sigma, lya_sigma_0;
    float x_in;
    float tau;

    thermal_vel = LyaThermalVel((double)(temp));
    nu_doppler = LyaNuDoppler((double)(temp));
    a = Lya_nu_line_width_CGS/(2.0*nu_doppler);
   
    x_in = x_value - (vel_x*dir_x + vel_y*dir_y + vel_z*dir_z)/thermal_vel;



    lya_sigma_0  = LyaCrossSection(0.0, a);
    lya_sigma  = LyaCrossSection(x_in, a);

    tau   =  (delta_step/cell_size)*(lya_sigma/lya_sigma_0)*tau_HI;

    return tau;
}

void RT_DustAbsorption(float *intensity, float tau_HI, float temp, float tau_dust){
    float I;
    float escape;
    float xsi;
    double nu_doppler;
    double a;
    I = *intensity;

    nu_doppler = LyaNuDoppler((double)(temp));
    a = Lya_nu_line_width_CGS/(2.0*nu_doppler);

    xsi = 0.525;/*fitting factor from verhamme, schaerer and maselli*/
    escape = pow((a*tau_HI), (1.0/3.0))*tau_dust;
    escape = (sqrt(3)/(pow(PI,5.0/12.0)*xsi))*sqrt(escape);
    escape = 1.0/cosh(escape);
    
    *intensity = I*escape;
}

void RT_LymanCube(float *x_value, float *Intensity, float tau_HI, float temp, float tau_dust, 
	float *dir_x, float *dir_y, float *dir_z, float vel_x, float vel_y, float vel_z, float cube_size){

    double nu_doppler, a, thermal_vel;    
    double lya_sigma, lya_sigma_0;
    double n_HI;
    double x_in;
    float tau;
    double Dir[3];
    double DeltaPos[3];
    Dir[0] = *dir_x;
    Dir[1] = *dir_y;
    Dir[2] = *dir_z;

//    fprintf(stdout, "Lyman cube vel %e %e %e\n", vel_x, vel_y, vel_z);
    /*get the thermal properties*/

    thermal_vel = LyaThermalVel((double)(temp));
    nu_doppler = LyaNuDoppler((double)(temp));
    a = Lya_nu_line_width_CGS/(2.0*nu_doppler);
    
    lya_sigma_0 =  LyaCrossSection(0,a);
    n_HI = tau_HI/(0.5*cube_size*All.UnitLength_in_cm*lya_sigma_0);


    /*move the frequency to frame comoving with the fluid*/

    x_in = (*x_value) - (vel_x*Dir[0] + vel_y*Dir[1] + vel_z*Dir[2])/thermal_vel ;
//    fprintf(stdout, "atau %e tau %e\n", a*tau_HI/0.5, tau_HI/0.5);

    RND_spherical(Dir);
    RT_LymanSimple(Dir, DeltaPos, &x_in, a, tau_HI, All.EffectiveEta);
    /*
    if(a*tau_HI > 1.0){
	if(a*tau_HI > All.ThresholdCube){

	}else{
	    RT_LymanUngrid(Dir, DeltaPos, &x_in, a, tau_HI, n_HI);
	}
    }
    */

    /*move the frequency again to the observers value*/
    *x_value = x_in + (vel_x*Dir[0] + vel_y*Dir[1] + vel_z*Dir[2])/thermal_vel;

    *dir_x = Dir[0];
    *dir_y = Dir[1];
    *dir_z = Dir[2];
}

void RecordHistoryPhoton(float *pos_x, float *pos_y, float *pos_z, float *d_x, float *d_y, float *d_z, int *n_scatter, 
			 grid *GasTauHI,     lyman_RT_record *H){
    float old_dir[3],old_pos[3];
    float delta_step;
    int n_cells_scatter;
    int is_in;
    /*Initialize photon variables*/

    delta_step = GasTauHI->cell_size;
    old_pos[0] = *pos_x;
    old_pos[1] = *pos_y;
    old_pos[2] = *pos_z;
    old_dir[0] = *d_x;
    old_dir[1] = *d_y;
    old_dir[2] = *d_z;
    n_cells_scatter = 0;
    
    is_in = InsideBox(old_pos[0], old_pos[1], old_pos[2], GasTauHI);  
    //    fprintf(stdout, "is in %d\n", is_in);fflush(stdout);

    /*First check that the initial position is inside the domain*/
    if(is_in>0){
	/*Make the loop*/
	while(is_in>0 && n_cells_scatter<MAX_N_SCATTER){    
	    is_in = InsideBox(old_pos[0], old_pos[1], old_pos[2], GasTauHI);  	    
	    H->index_grid[n_cells_scatter] = is_in;
	    RND_sphericalf(old_dir);	    
	    old_pos[0] = old_pos[0] + old_dir[0]*delta_step;
	    old_pos[1] = old_pos[1] + old_dir[1]*delta_step;
	    old_pos[2] = old_pos[2] + old_dir[2]*delta_step;
	    n_cells_scatter++;
	}	
    }
    fprintf(stdout, "n_cells_scatter %d\n", n_cells_scatter);fflush(stdout);
    H->N_steps = n_cells_scatter;
    *n_scatter = n_cells_scatter;
}

void PropagatePhoton(float *pos_x, float *pos_y, float *pos_z, float *d_x, float *d_y, float *d_z, float *x, float *I,
		     int *nscatterHI, int *ncatterDust, int *nwrong, 
		     grid *GasTauHI, grid *GasTauDust, grid *GasTemp,
		     grid *GasVelX, grid *GasVelY, grid *GasVelZ){

  int n_scatter_HI, n_wrong, n_scatter_dust;
  int new_pos_i, new_pos_j, new_pos_k;
  float old_dir[3],old_pos[3];
  int old_pos_cell[3];
  int is_in;
  int index;
  float tau_HI, temp, tau_dust, vel_x, vel_y, vel_z;
  float effective_tau_HI;
  float a;
  float x_lab, x_old, x_new;
  
  float intensity;
  float delta_step;
  float move_delta_step;
  double CellSize;
  double v_thermal, nu_doppler;
  double new_nu_doppler;
  double old_nu_doppler, old_v_thermal, old_vel_x, old_vel_y, old_vel_z;
  double tau, local_tau;
  float x_value;
  int n_cells_scatter;
  

  /*Initialize photon variables*/
  delta_step = GasTauHI->cell_size;
  
  old_pos[0] = *pos_x;
  old_pos[1] = *pos_y;
  old_pos[2] = *pos_z;
  old_dir[0] = *d_x;
  old_dir[1] = *d_y;
  old_dir[2] = *d_z;
  x_value = *x;
  intensity = *I;
  
  /* Initialize bookkeep variables*/
  n_scatter_HI = 0;
  n_scatter_dust = 0;
  n_wrong = 0;
  n_cells_scatter = 0;
  
  is_in = InsideBox(old_pos[0], old_pos[1], old_pos[2], GasTauHI);  

  /*First check that the initial position is inside the domain*/
  if(is_in){
      
      /* Initialize the physical state of the cell*/
      tau_HI     = GridGetValue(GasTauHI,    old_pos[0], old_pos[1], old_pos[2]);
      temp       = GridGetValue(GasTemp,     old_pos[0], old_pos[1], old_pos[2]);
      vel_x      = GridGetValue(GasVelX,     old_pos[0], old_pos[1], old_pos[2]);
      vel_y      = GridGetValue(GasVelY,     old_pos[0], old_pos[1], old_pos[2]);
      vel_z      = GridGetValue(GasVelZ,     old_pos[0], old_pos[1], old_pos[2]);
      tau_dust   = GridGetValue(GasTauDust,  old_pos[0], old_pos[1], old_pos[2]);
      
      v_thermal  = LyaThermalVel((double)(temp));
      nu_doppler = LyaNuDoppler((double)(temp));
      
      
      
      /*values at entry*/
      is_in = InsideBox(old_pos[0], old_pos[1], old_pos[2], GasTauHI);  
      /*Make the loop*/
      while(is_in){    
	  /*keep some values of the old cell*/
	  old_v_thermal = v_thermal;
	  old_nu_doppler = nu_doppler;
	  old_vel_x = vel_x;
	  old_vel_y = vel_y;
	  old_vel_z = vel_z;
	  
	  
	  RT_LymanCube(&x_value, &intensity, tau_HI, temp, tau_dust, 
		       &old_dir[0], &old_dir[1], &old_dir[2], vel_x, vel_y, vel_z, delta_step);
	  RT_DustAbsorption(&intensity, tau_HI, temp, tau_dust);      
	  
	  old_pos[0] = old_pos[0] + old_dir[0]*delta_step;
	  old_pos[1] = old_pos[1] + old_dir[1]*delta_step;
	  old_pos[2] = old_pos[2] + old_dir[2]*delta_step;
	  
	  is_in = InsideBox(old_pos[0], old_pos[1], old_pos[2], GasTauHI);      
	  if(is_in){	  	  
	      /* find the updated physical values at the new position*/
	      tau_HI     = GridGetValue(GasTauHI,    old_pos[0], old_pos[1], old_pos[2]);
	      temp       = GridGetValue(GasTemp,     old_pos[0], old_pos[1], old_pos[2]);
	      tau_dust   = GridGetValue(GasTauDust,  old_pos[0], old_pos[1], old_pos[2]);
	      vel_x      = GridGetValue(GasVelX,     old_pos[0], old_pos[1], old_pos[2]);
	      vel_y      = GridGetValue(GasVelY,     old_pos[0], old_pos[1], old_pos[2]);
	      vel_z      = GridGetValue(GasVelZ,     old_pos[0], old_pos[1], old_pos[2]);	    
	      v_thermal  = LyaThermalVel((double)(temp));
	      nu_doppler = LyaNuDoppler((double)(temp));
	      
	      /*modify the value of x, accounting for the change in temperature, the frequency here 
		is always the seen in the lab */
	      x_value = x_value*(old_nu_doppler/nu_doppler);	  
	      
	  }
	  n_cells_scatter++;
      }
  }
  fprintf(stdout, "N scatter %d\n", n_cells_scatter);
  
  *nscatterHI = n_scatter_HI;
  *nwrong = n_wrong;
  *pos_x = old_pos[0];
  *pos_y = old_pos[1];
  *pos_z = old_pos[2];
  *d_x   = old_dir[0];
  *d_y   = old_dir[1];
  *d_z   = old_dir[2];
  *x     = (x_value*nu_doppler)*Lya_lambda_line_center_CGS/C_LIGHT;/*to get things out in units of wavelength */
  *I     = intensity;
}


void TracePhotonDirty(lyman_RT_record *R, float *pos_x, float *pos_y, float *pos_z, float *d_x, float *d_y, float *d_z, 
		      float *total_dust_escape, float *total_H_tau, int *nscatterHI,
		      grid *GasTauHI, grid *GasTauDust, grid *GasTemp,
		      grid *GasVelX, grid *GasVelY, grid *GasVelZ){
  
  int n_scatter_HI, n_wrong, n_scatter_dust;
  int new_pos_i, new_pos_j, new_pos_k;
  float old_dir[3],old_pos[3];
  int old_pos_cell[3];
  int is_in;
  int index;
  float tau_HI, temp, tau_dust, vel_x, vel_y, vel_z;
  float effective_tau_HI;
  float a;
  float x_lab, x_old, x_new;
  
  float intensity;
  float delta_step;
  float move_delta_step;
  double CellSize;
  double v_thermal, nu_doppler;
  double new_nu_doppler;
  double old_nu_doppler, old_v_thermal, old_vel_x, old_vel_y, old_vel_z;
  double tau, local_tau;
  double total_tau, total_tau_a_0, total_tau_dust, a_0, max_ext;
  float x_value;
  int n_cells_scatter;
  int i;
  float escape_fraction, total_escape_fraction;

  /*Initialize photon variables*/
  delta_step = GasTauHI->cell_size;
  
  old_pos[0] = *pos_x;
  old_pos[1] = *pos_y;
  old_pos[2] = *pos_z;
  old_dir[0] = *d_x;
  old_dir[1] = *d_y;
  old_dir[2] = *d_z;


  
  /* Initialize bookkeep variables*/
  n_scatter_HI = 0;
  n_cells_scatter = 0;

  total_tau = 0.0;
  total_tau_a_0 = 0.0;
  total_escape_fraction = 1.0; /* The one that goes in the extinction formulae */
  max_ext = 0.0;

  //  for(i=0;i<R->N_steps;i++){      
  for(i=0;i<R->N_steps - 1;i++){      
    /* Initialize the physical state of the cell*/
    tau_HI     = GridGetValueIndex(GasTauHI,   R->index_grid[i]);
    temp       = GridGetValueIndex(GasTemp,    R->index_grid[i]);
    vel_x      = GridGetValueIndex(GasVelX,    R->index_grid[i]);
    vel_y      = GridGetValueIndex(GasVelY,    R->index_grid[i]);
    vel_z      = GridGetValueIndex(GasVelZ,    R->index_grid[i]);
    tau_dust   = GridGetValueIndex(GasTauDust, R->index_grid[i]);

    if(temp<0.0)
      temp =1.0e4;

    nu_doppler = LyaNuDoppler((double)(temp));
    a_0 = Lya_nu_line_width_CGS/(2.0*nu_doppler);

    escape_fraction = 1.0;
    RT_DustAbsorption(&escape_fraction, tau_HI, temp, tau_dust);
    total_escape_fraction *= escape_fraction;      

    if(total_escape_fraction < 1.0e-10){
      total_escape_fraction = 1.0e-10;
    }
      
    total_tau += tau_HI * a_0;


    /*keep some values of the old cell*/
    old_v_thermal = v_thermal;
    old_nu_doppler = nu_doppler;
    old_vel_x = vel_x;
    old_vel_y = vel_y;
    old_vel_z = vel_z;
    n_cells_scatter++;
  }
  fprintf(stdout, "OUT %d %g %g %g %g\n", i, tau_HI, a_0, total_tau, escape_fraction);
  

  *pos_x = old_pos[0];
  *pos_y = old_pos[1];
  *pos_z = old_pos[2];
  *d_x   = old_dir[0];
  *d_y   = old_dir[1];
  *d_z   = old_dir[2];
  *total_dust_escape =  total_escape_fraction;
  *nscatterHI = n_cells_scatter;
  
}


void TracePhotonQuick(lyman_RT_record *R, float *pos_x, float *pos_y, float *pos_z, float *d_x, float *d_y, float *d_z,
		      float *total_dust_escape, float *x, int *nscatterHI,
		      grid *GasTauHI, grid *GasTauDust, grid *GasTemp,
		      grid *GasVelX, grid *GasVelY, grid *GasVelZ){
  
  int n_scatter_HI, n_wrong, n_scatter_dust;
  int new_pos_i, new_pos_j, new_pos_k;
  float old_dir[3],old_pos[3];
  int old_pos_cell[3];
  int is_in;
  int index;
  float tau_HI, temp, tau_dust, vel_x, vel_y, vel_z;
  float effective_tau_HI;
  float a;
  float x_lab, x_old, x_new;
  
  float intensity;
  float delta_step;
  float move_delta_step;
  double CellSize;
  double v_thermal, nu_doppler;
  double new_nu_doppler;
  double old_nu_doppler, old_v_thermal, old_vel_x, old_vel_y, old_vel_z;
  double tau, local_tau;
  double total_tau, total_tau_a_0, total_tau_dust, a_0, max_ext;
  float x_value;
  int n_cells_scatter;
  int i;
  float escape_fraction, total_escape_fraction;

  /*Initialize photon variables*/
  delta_step = GasTauHI->cell_size;
  
  old_pos[0] = *pos_x;
  old_pos[1] = *pos_y;
  old_pos[2] = *pos_z;
  old_dir[0] = *d_x;
  old_dir[1] = *d_y;
  old_dir[2] = *d_z;
  x_value = *x;

  
  /* Initialize bookkeep variables*/
  n_scatter_HI = 0;
  n_cells_scatter = 0;

  total_tau = 0.0;
  total_tau_a_0 = 0.0;
  total_escape_fraction = 1.0; /* The one that goes in the extinction formulae */
  max_ext = 0.0;

  //  for(i=0;i<R->N_steps;i++){      
  for(i=0;i<R->N_steps - 1;i++){      
    /* Initialize the physical state of the cell*/
    tau_HI     = GridGetValueIndex(GasTauHI,   R->index_grid[i]);
    temp       = GridGetValueIndex(GasTemp,    R->index_grid[i]);
    vel_x      = GridGetValueIndex(GasVelX,    R->index_grid[i]);
    vel_y      = GridGetValueIndex(GasVelY,    R->index_grid[i]);
    vel_z      = GridGetValueIndex(GasVelZ,    R->index_grid[i]);
    tau_dust   = GridGetValueIndex(GasTauDust, R->index_grid[i]);

    nu_doppler = LyaNuDoppler((double)(temp));
    a_0 = Lya_nu_line_width_CGS/(2.0*nu_doppler);





    /*keep some values of the old cell*/
    old_v_thermal = v_thermal;
    old_nu_doppler = nu_doppler;
    old_vel_x = vel_x;
    old_vel_y = vel_y;
    old_vel_z = vel_z;
    n_cells_scatter++;

    RT_LymanCube(&x_value, &intensity, tau_HI, temp, tau_dust, 
		 &old_dir[0], &old_dir[1], &old_dir[2], vel_x, vel_y, vel_z, delta_step);

  }
  fprintf(stdout, "OUT %e %d %g %g %g %g\n", i, tau_HI, a_0, total_tau, escape_fraction);
  

  *pos_x = old_pos[0];
  *pos_y = old_pos[1];
  *pos_z = old_pos[2];
  *d_x   = old_dir[0];
  *d_y   = old_dir[1];
  *d_z   = old_dir[2];
  *total_dust_escape =  total_escape_fraction;
  *x = x_value;
  *nscatterHI = n_cells_scatter;
  
}

void TracePhoton(lyman_RT_record *R, float *pos_x, float *pos_y, float *pos_z, float *d_x, float *d_y, float *d_z, float *x, float *I,
		 int *nscatterHI, int *ncatterDust, int *nwrong, 
		 grid *GasTauHI, grid *GasTauDust, grid *GasTemp,
		 grid *GasVelX, grid *GasVelY, grid *GasVelZ){
  
  int n_scatter_HI, n_wrong, n_scatter_dust;
  int new_pos_i, new_pos_j, new_pos_k;
  float old_dir[3],old_pos[3];
  int old_pos_cell[3];
  int is_in;
  int index;
  float tau_HI, temp, tau_dust, vel_x, vel_y, vel_z;
  float effective_tau_HI;
  float a;
  float x_lab, x_old, x_new;
  
  float intensity;
  float delta_step;
  float move_delta_step;
  double CellSize;
  double v_thermal, nu_doppler;
  double new_nu_doppler;
  double old_nu_doppler, old_v_thermal, old_vel_x, old_vel_y, old_vel_z;
  double tau, local_tau;
  double total_tau, total_tau_a_0, total_tau_dust, a_0, max_ext;
  float x_value;
  int n_cells_scatter;
  int i;

  /*Initialize photon variables*/
  delta_step = GasTauHI->cell_size;
  
  old_pos[0] = *pos_x;
  old_pos[1] = *pos_y;
  old_pos[2] = *pos_z;
  old_dir[0] = *d_x;
  old_dir[1] = *d_y;
  old_dir[2] = *d_z;
  x_value = *x;
  intensity = *I;
  
  /* Initialize bookkeep variables*/
  n_scatter_HI = 0;
  n_scatter_dust = 0;
  n_wrong = 0;
  n_cells_scatter = 0;

  total_tau = 0.0;
  total_tau_a_0 = 1.0e8;
  total_tau_dust = 0.0;
  max_ext = 0.0;

  for(i=0;i<R->N_steps;i++){      
    /* Initialize the physical state of the cell*/
    tau_HI     = GridGetValueIndex(GasTauHI,   R->index_grid[i]);
    temp       = GridGetValueIndex(GasTemp,    R->index_grid[i]);
    vel_x      = GridGetValueIndex(GasVelX,    R->index_grid[i]);
    vel_y      = GridGetValueIndex(GasVelY,    R->index_grid[i]);
    vel_z      = GridGetValueIndex(GasVelZ,    R->index_grid[i]);
    tau_dust   = GridGetValueIndex(GasTauDust, R->index_grid[i]);

    if(temp<0.0)
      temp =1.0e4;


    nu_doppler = LyaNuDoppler((double)(temp));
    a_0 = Lya_nu_line_width_CGS/(2.0*nu_doppler);
    //    fprintf(stdout, "%g\n", a_0);
    
    if(pow(tau_HI * a_0,0.333)*tau_dust>max_ext)
      max_ext =       pow(tau_HI * a_0,0.333)*tau_dust;
    
    total_tau_dust += tau_dust;
    total_tau += tau_HI;

    if(nu_doppler < total_tau_a_0)
      total_tau_a_0 = nu_doppler;


    /*keep some values of the old cell*/
    old_v_thermal = v_thermal;
    old_nu_doppler = nu_doppler;
    old_vel_x = vel_x;
    old_vel_y = vel_y;
    old_vel_z = vel_z;
    
    /*
    RT_LymanCube(&x_value, &intensity, tau_HI, temp, tau_dust, 
		 &old_dir[0], &old_dir[1], &old_dir[2], vel_x, vel_y, vel_z, delta_step);
    RT_DustAbsorption(&intensity, tau_HI, temp, tau_dust);      
    */
    n_cells_scatter++;
  }
  fprintf(stdout, "N scatter %d %g %g %g %g\n", n_cells_scatter, total_tau, total_tau_a_0, max_ext, total_tau_dust);
  
  *nscatterHI = n_scatter_HI;
  *nwrong = n_wrong;
  *pos_x = old_pos[0];
  *pos_y = old_pos[1];
  *pos_z = old_pos[2];
  *d_x   = old_dir[0];
  *d_y   = old_dir[1];
  *d_z   = old_dir[2];
  *x     = (x_value*nu_doppler)*Lya_lambda_line_center_CGS/C_LIGHT;/*to get things out in units of wavelength */
  *I     = intensity;
}


lyman_RT_record * AllocateHistory(grid *G){
    int n_max_points;
    lyman_RT_record  *H;
    H =  malloc(sizeof(lyman_RT_record));
    
    H->N_max_steps = MAX_N_SCATTER;
    H->N_steps = 0;
    if(!(H->index_grid=malloc(sizeof(int)*MAX_N_SCATTER))){
	fprintf(stderr, "problem allocating index_grid for the history\n");
	exit(1);
    }

    return H;
}

void NullHistory(lyman_RT_record *H){
    int i;
    H->N_steps = 0;
    for(i=0;i<H->N_max_steps;i++){
	H->index_grid[i] = -1;
    }
}



void RecordPhotons(void){
    char FileName[MAX_FILENAME_SIZE];
    int i;
    int N_photons_per_proc;
    int n_scatt;
    FILE *out;
    lyman_RT_photons *PhList;

    grid *GasTauDust;
    grid *GasTauHI;
    grid *GasTemp;
    grid *GasVelX;
    grid *GasVelY;
    grid *GasVelZ;
    lyman_RT_stars *StarList;    
    lyman_RT_record  *PhHistory;
   
    /*create the grids*/

    GasTauDust = GridCreate(); 
    GasTauHI   = GridCreate();
    GasTemp    = GridCreate();
    GasVelX    = GridCreate();
    GasVelY    = GridCreate();
    GasVelZ    = GridCreate();
    
    StarList = malloc(sizeof(lyman_RT_stars));
    
    /*load the data*/
    RT_LoadGrids(GasTauDust, GasTauHI, GasTemp, GasVelX, GasVelY, GasVelZ);

    /*load the star luminosity data*/
    RT_LoadStars(StarList);

    /*Initialize the photon list*/
    N_photons_per_proc = PhotonNumber(StarList);
    PhList  = PhotonListCreate(N_photons_per_proc);

    InitializePhotonList(StarList, GasVelX, GasVelY, GasVelZ, GasTemp, PhList, SizeProc);

    /* Allocate the array for the individual photon histories */
    PhHistory = AllocateHistory(GasTauHI);

    
    /*write the input list*/
    if(All.OutputInitList){
	sprintf(FileName , "%s/%s_%s.%d.dat", All.OutputDir, All.OutputFile, "ph_in", ThisProc);
	if(All.OutputBinary){
	    SavePhotonListBinary(FileName, PhList);
	}else{
	    sprintf(FileName , "%s/%s_%s.%d.dat.ascii", All.OutputDir, All.OutputFile, "ph_in", ThisProc);
	    SavePhotonListAscii(FileName, PhList);
	}
    }

    
    /*open the file with the history*/
    sprintf(FileName , "%s/%s_%s.proc_%d.dat", All.OutputDir, All.HistoryName, "hist", ThisProc);
    if(!(out=fopen(FileName, "w"))){
	fprintf(stdout, "Problem opening file %s\n", FileName);
	exit(1);
    }
    
    /*Trace the random path of the photons*/
    fwrite(&N_photons_per_proc, sizeof(int), 1, out);
    
    fprintf(stdout, "N_photons_per_proc %d\n", N_photons_per_proc);fflush(stdout);
    MPI_Barrier (MPI_COMM_WORLD);

    for(i=0;i<N_photons_per_proc;i++){	
	NullHistory(PhHistory);
	RecordHistoryPhoton(&(PhList->Pos[3*i + 0]), &(PhList->Pos[3*i + 1]), &(PhList->Pos[3*i + 2]), 
			    &(PhList->Dir[3*i + 0]), &(PhList->Dir[3*i + 1]), &(PhList->Dir[3*i + 2]), &(PhList->ScatterHI[i]),
			    GasTauHI, PhHistory);	    
	
	fwrite(&i, sizeof(int), 1, out);
	fwrite(&(PhList->Pos[3*i]), sizeof(float), 3, out);
	fwrite(&(PhList->Lum[i]), sizeof(float), 1, out);
	fwrite(&(PhHistory->N_steps), sizeof(int), 1, out);
	if(PhHistory->N_steps > 0){
	    fwrite(&(PhHistory->index_grid[0]), sizeof(int), PhHistory->N_steps, out);
	}
	//	if(!(i%(N_photons_per_proc/100)))
	fprintf(stdout, "Done %d : Proc %d\n", i, ThisProc);fflush(stdout);
    }
    
    /* Close the History file */
    fclose(out);
    fprintf(stdout, "Done!\n");fflush(stdout);
    
    
    /*write the output list*/
    if(All.OutputFinalList){
	sprintf(FileName , "%s/%s_%s.%d.dat", All.OutputDir, All.OutputFile, "ph_out", ThisProc);
	if(All.OutputBinary){
	    SavePhotonListBinary(FileName, PhList);
	}else{
	    sprintf(FileName , "%s/%s_%s.%d.dat.ascii", All.OutputDir, All.OutputFile, "ph_out", ThisProc);
	    SavePhotonListAscii(FileName, PhList);
	}
    }  
    
    
    PhotonFree(PhList);
    GridFree(GasTauHI);
    GridFree(GasTemp);
    GridFree(GasVelX);
    GridFree(GasVelY);
    GridFree(GasVelZ);
}

void TracePhotonsFromHistory(void){  
    char FileName[MAX_FILENAME_SIZE];
    int i;
    int N_photons_per_proc;
    int n_scatt;
    int index;
    FILE *in;
    lyman_RT_photons *PhList;

    grid *GasTauDust;
    grid *GasTauHI;
    grid *GasTemp;
    grid *GasVelX;
    grid *GasVelY;
    grid *GasVelZ;
    lyman_RT_stars *StarList;    
    lyman_RT_record  *PhHistory;

    

    /*create the grids*/
    GasTauDust = GridCreate(); 
    GasTauHI   = GridCreate();
    GasTemp    = GridCreate();
    GasVelX    = GridCreate();
    GasVelY    = GridCreate();
    GasVelZ    = GridCreate();
    
    StarList = malloc(sizeof(lyman_RT_stars));
    
    /*load the data*/
    RT_LoadGrids(GasTauDust, GasTauHI, GasTemp, GasVelX, GasVelY, GasVelZ);

    /*open the file with the history*/
    sprintf(FileName , "%s/%s_%s.proc_%d.dat", All.OutputDir, All.HistoryName, "hist", ThisProc);
    if(!(in=fopen(FileName, "r"))){
      fprintf(stdout, "Problem opening file %s\n", FileName);
      exit(1);
    }
    
    
    fflush(stdout);
    fread(&N_photons_per_proc, sizeof(int), 1, in);	
    fprintf(stdout, "N_photons_per_proc to read (proc %d):%d\n", ThisProc, N_photons_per_proc); fflush(stdout);
    /*Initialize the photon list*/
    PhList  = PhotonListCreate(N_photons_per_proc);

    /* Allocate the array for the individual photon histories */
    PhHistory = AllocateHistory(GasTauHI);


    for(i=0;i<N_photons_per_proc;i++){	
      NullHistory(PhHistory);
      fread(&index, sizeof(int), 1, in);
      fread(&(PhList->Pos[3*i]), sizeof(float), 3, in);
      fread(&(PhList->Lum[i]), sizeof(float), 1, in);
      fread(&(PhHistory->N_steps), sizeof(int), 1, in);
#ifdef DEBUG
      fprintf(stdout, "posx  %g \n", PhList->Pos[3*i]);
      fprintf(stdout, "lum  %g \n", PhList->Lum[i]);
      fprintf(stdout, "index %d\n", index);
      fprintf(stdout, "N_steps %d\n", PhHistory->N_steps);
#endif
      if(PhHistory->N_steps>0){
	fread(&(PhHistory->index_grid[0]), sizeof(int), PhHistory->N_steps, in);


	if(All.DirtyRT){
	  TracePhotonDirty(PhHistory, 
			   &(PhList->Pos[3*i + 0]), &(PhList->Pos[3*i + 1]), &(PhList->Pos[3*i + 2]), 
			   &(PhList->Dir[3*i + 0]), &(PhList->Dir[3*i + 1]), &(PhList->Dir[3*i + 2]), 
			   &(PhList->Intensity[i]), &(PhList->x_out[i]), &(PhList->ScatterHI[i]), 
			   GasTauHI, GasTauDust, GasTemp, GasVelX, GasVelY, GasVelZ);    	
#ifdef DEBUG
	  fprintf(stdout, "%e %e %d\n", PhList->Intensity[i], PhList->x_out[i], PhList->ScatterHI[i]);
#endif
	}

	if(All.QuickRT){
	  TracePhotonQuick(PhHistory, 
			   &(PhList->Pos[3*i + 0]), &(PhList->Pos[3*i + 1]), &(PhList->Pos[3*i + 2]), 
			   &(PhList->Dir[3*i + 0]), &(PhList->Dir[3*i + 1]), &(PhList->Dir[3*i + 2]), 
			   &(PhList->Intensity[i]), &(PhList->x_out[i]), &(PhList->ScatterHI[i]), 
			   GasTauHI, GasTauDust, GasTemp, GasVelX, GasVelY, GasVelZ);    	
#ifdef DEBUG
	  fprintf(stdout, "%e %e %d\n", PhList->Intensity[i], PhList->x_out[i], PhList->ScatterHI[i]);
#endif
	}
	
	if(All.ArchiveRT){
	/*
	TracePhoton(PhHistory, 
		    &(PhList->Pos[3*i + 0]), &(PhList->Pos[3*i + 1]), &(PhList->Pos[3*i + 2]), 
		    &(PhList->Dir[3*i + 0]), &(PhList->Dir[3*i + 1]), &(PhList->Dir[3*i + 2]), 
		    &(PhList->x_out[i]), &(PhList->Intensity[i]),
		    &(PhList->ScatterHI[i]), &(PhList->ScatterDust[i]), &(PhList->Wrong[i]), 
		    GasTauHI, GasTauDust, GasTemp, GasVelX, GasVelY, GasVelZ);    	
	*/	  
	}



      }
      fprintf(stdout, "Done read and trace (proc %d):%d / %d \n", ThisProc, i, N_photons_per_proc);fflush(stdout);
      if(index>N_photons_per_proc){
	fprintf(stdout, " Hared Limit index Done read and trace (proc %d):%d / %d \n", ThisProc, index, N_photons_per_proc);fflush(stdout);
	exit(1);
      }
    }
    fclose(in);
    fprintf(stdout, "Done!\n");fflush(stdout);

    /*write the output list*/
    if(All.OutputFinalList){
	sprintf(FileName , "%s/%s_%s.proc.%d.dat", All.OutputDir, All.OutputFile, "ph_out", ThisProc);
	if(All.OutputBinary){
	    SavePhotonListBinary(FileName, PhList);
	}else{
	  sprintf(FileName , "%s/%s_%s.proc.%d.ascii", All.OutputDir, All.OutputFile, "ph_out", ThisProc);
	    SavePhotonListAscii(FileName, PhList);
	}
    }  

    
    PhotonFree(PhList);
    GridFree(GasTauHI);
    GridFree(GasTemp);
    GridFree(GasVelX);
    GridFree(GasVelY);
    GridFree(GasVelZ);
}
