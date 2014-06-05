#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "struct.h"
#include "myrand.h"
#include "lyman.h"
#include "RND_lyman.h"
#include "vector.h"



void RND_spherical(double *vec)
/*vector randomly distributed over the sphere*/
{
    double theta, phi;
    theta = acos(2.0*(RandFloatUnit()-0.5));
    phi = 2.0*PI*RandFloatUnit();
    vec[0] = sin(theta)*cos(phi);
    vec[1] = sin(theta)*sin(phi);
    vec[2] = cos(theta);
}

void RND_pair(double *r_1, double *r_2){
  int finished;
  double rand_1, rand_2;
  finished = 0;
  while (finished == 0) {
    rand_1 = 2.0*(RandFloatUnit() - 0.5);
    rand_2 = 2.0*(RandFloatUnit() - 0.5);

    if(rand_1*rand_1 + rand_2*rand_2 < 1.0){
      finished = 1;
    }
  }    
  *r_1 = rand_1;
  *r_2 = rand_2;
}

void RND_lyman_perp_vel_test(void)
/*generates a series  of velocities for a set of values of x and a */
{
    char FileName[MAX_FILENAME_SIZE];
    FILE *out;
    int i;
    double vel_1, vel_2;

    sprintf(FileName, "%s/%s_perp_vel.proc_%d.dat",All.OutputDir, All.OutputTestFile, ThisProc);
    if(!(out=fopen(FileName, "w"))){
	fprintf(stderr, "RND_lyman_perp_vel_testL problem opening file %s\n", FileName);
	MPI_Abort(MPI_COMM_WORLD, 0);
    }
    fprintf(out, "%d %e %e\n", N_POINTS_IN_TEST, 0.0, 0.0);
    for(i=0;i<N_POINTS_IN_TEST;i++){
	RND_lyman_perp_vel(&vel_1,&vel_2);
	fprintf(out, "%e %e\n", vel_1, vel_2);
    }
    
    fclose(out);
}


void RND_lyman_parallel_vel_test(double x, double a)
/*generates a series  of velocities for a set of values of x and a */
{
    char FileName[MAX_FILENAME_SIZE];
    FILE *out;
    int i;
    double vel;
    sprintf(FileName, "%s/%s_par_vel.proc_%d.dat",All.OutputDir, All.OutputTestFile, ThisProc);
    if(!(out=fopen(FileName, "w"))){
	fprintf(stderr, "RND_lyman_parallel_vel_testL problem opening file %s\n", FileName);
	MPI_Abort(MPI_COMM_WORLD, 0);
    }
    fprintf(out, "%d %e %e\n", N_POINTS_IN_TEST, All.Test_x, All.Test_a);
    for(i=0;i<N_POINTS_IN_TEST;i++){
	vel = RND_lyman_parallel_vel(x,a);
	fprintf(out, "%e\n", vel);
    }
    
    fclose(out);
}




void RND_lyman_parallel_vel_fast_test(double x, double a)
/*generates a series  of velocities for a set of values of x and a */
{
    char FileName[MAX_FILENAME_SIZE];
    FILE *out;
    int i;
    double vel;

    sprintf(FileName, "%s/%s_par_vel_fast.proc_%d.dat",All.OutputDir, All.OutputTestFile, ThisProc);
    if(!(out=fopen(FileName, "w"))){
	fprintf(stderr, "RND_lyman_parallel_vel_testL problem opening file %s\n", FileName);
	MPI_Abort(MPI_COMM_WORLD, 0);
    }
    fprintf(out, "%d %e %e\n", N_POINTS_IN_TEST, All.Test_x, All.Test_a);
    for(i=0;i<N_POINTS_IN_TEST;i++){
	vel = RND_lyman_parallel_vel(x,a);
	fprintf(out, "%e\n", vel);
    }
    
    fclose(out);
}


double RND_lyman_parallel_vel(double x, double a) 
/* 
   it generates a random number \theta between -pi/2 and pi/2, then 
   it generates u_parallel through u_parallel = a \tan\theta + x
   then this value of u_parallel is kept if another random number 
   between [0,1] is smaller than \exp^{-u_parallel^2}. 
   At the end the value of u_parallel is multiplied by the sign of x.       
*/
{
    int finished = 0;
    double tmp0, tmp1, tmp2;        
    int counter = 0;

    while (finished == 0) {
	tmp0 = ((RandFloatUnit() - 0.5))*PI ;
	tmp1 = (a*tan(tmp0)) + fabs(x); 
	tmp2  = RandFloatUnit();
	if(tmp2 <= (exp(-(tmp1*tmp1)))) finished = 1;		
	counter++;		
	if(counter > MAX_VEL_ITER) {
	    finished = 1;		    
	    fprintf(stderr, "Warning: RND_MT_parallel_vel MAX_VEL_ITER reached... continuing, but PDF is now biased.\n");
	}	
    }
    /* fprintf(stdout, "counts %d %e %e \n", counter, x, a); */
    /* now change the sign accordingly*/

//   fprintf(stdout, "scale %e\n", scale);
    if(x > 0.0)
    {
	return tmp1;
    }
    else
    {
	return -tmp1;
    }
}



double RND_lyman_parallel_vel_fast(double my_x, double a) 
    /* it generates a random number \theta between -pi/2 and pi/2, then 
       it generates u_parallel through u_parallel = a \tan\theta + x
       then this value of u_parallel is kept if another random number 
       between [0,1] is smaller than \exp^{-u_parallel^2}. 
       At the end the value of u_parallel is multiplied by the sign of x.       
       It uses the speeding mechanism of Zheng&Miralda-Scude, with
       a critical value taken from the paper of Semelin, Combes & Baek.
    */
{
    int finished = 0;
    double tmp0, tmp1, tmp2, tmp3;

    double u_critical,  theta_0, p_ratio;
    int counter = 0;
    double final, x;
    x = fabs(my_x);
    /* first find the value of u_critical */
    if(x > 3)
    {
	u_critical = 1.85 - log(a)/6.73 + log(x)*log(x);    
	theta_0 = atan((u_critical + x)/a);        
	p_ratio = (theta_0 + PI*0.5)/((1.0 - exp(-u_critical*u_critical))*theta_0 + 
				      (1.0 + exp(-u_critical*u_critical))*PI*0.5);        
	/* printf("%g\n", u_critical); */
    }
    else
    {
	u_critical = 0.0;
	theta_0  = PI*0.5;        
	p_ratio = 1.0; 
    }

    tmp0  = RandFloatUnit();       
    if(tmp0 <= p_ratio)/*use one side of the comparation function*/
    {        
	while (finished == 0) {
	    tmp1 = (RandFloatUnit()*(theta_0 + PI*0.5)) - PI*0.5;
	    tmp2 = a*tan(tmp1) + x;        
	    tmp3  = RandFloatUnit();
	    if(tmp3 <= (exp(-(tmp2*tmp2)))) finished = 1;
	    counter++;
	    if(counter > MAX_VEL_ITER) {
		finished = 1;    
		printf("Warning: RND_MT_parallel_vel_fast (< p_ratio)-- MAX_VEL_ITER reached... continuing, but PDF is now biased.\n");
	    }
	}
    }
    else/*use the other side of the comparation function*/
    {       
	while (finished == 0) {     
	    tmp1 = (RandFloatUnit()*(PI*0.5 - theta_0)) + theta_0;
	    tmp2 = a*tan(tmp1) + x;
	    tmp3  = RandFloatUnit();
	    if(tmp3 <= (exp(-(tmp2*tmp2))/exp(-(u_critical*u_critical)))) finished = 1;
	    counter++;
	    if(counter > MAX_VEL_ITER) {
		finished = 1;    
		printf("Warning: RND_MT_parallel_vel_fast (> p_ratio)-- MAX_VEL_ITER reached... continuing, but PDF is now biased.\n");
	    }
	}
    }
    
    
    /*now change the sign accordingly*/    

    if(my_x > 0.0)
    {
	final = tmp2;
    }
    else
    {
	final  = -tmp2;
    }    
  

    
    return final;
}

void RND_lyman_perp_vel(double *u_1, double *u_2)
/* here I generate the magnitudes of the perpendicular velocities using Box&Muller method*/
{
    double tmp1;
    double tmp2;
    double vel_1, vel_2;


    
    tmp1 = RandFloatUnit();
    tmp2 = RandFloatUnit();
    vel_1 = sqrt(-log(tmp1))*cos(2.0*PI*tmp2);
    vel_2 = sqrt(-log(tmp1))*sin(2.0*PI*tmp2);
    *u_1 = vel_1;
    *u_2 = vel_2;
    
    return;
}




void RND_lyman_atom(double *Vel, double *DirPhoton, 
		    double *DirOutPhoton, double x, double a)
/* obtains a random velocity for the atom.
 the velocity is in units of the termal velocity.*/
{
    int i;
    double LocalVel[3];
    double x_axis[3];
    double y_axis[3];
    double z_axis[3];
    double rand_axis[3];
    double R_1, R_2, R_3, T, mu, iso;
    double x_corewing;


    /*get first the parallel velocity*/
    LocalVel[2] = RND_lyman_parallel_vel(x,a);

    /*get the perpendicular velocity*/
    RND_lyman_perp_vel(&(LocalVel[0]), &(LocalVel[1]));

    /*get the axis in the coordinate system of the atom, where 
      the z direction is the propagation direction of the photon*/
    for(i=0;i<3;i++){
	z_axis[i] = DirPhoton[i];
    }

    /*get another random vector*/
    RND_spherical(rand_axis);

    /*make the cross product and get y_axis*/
    cross_product(z_axis, rand_axis, y_axis);

    /*make the cross product and get x_axis*/
    cross_product(y_axis, z_axis, x_axis);

    /*normalize the vectors*/
    normalize(x_axis);
    normalize(y_axis);
    normalize(z_axis);

    /*see if they are perpendicular*/
    rand_axis[0] = point_product(x_axis, z_axis);
    rand_axis[1] = point_product(x_axis, y_axis);
    rand_axis[2] = point_product(y_axis, z_axis);
    if(fabs(rand_axis[0]) + fabs(rand_axis[1])+ fabs(rand_axis[2])>1.0e-10){
	fprintf(stderr, "not orthogonal\n");
	exit(1);
    }

    /*Now make the transformation into the coordinate frame of the lab*/
    for(i=0;i<3;i++){
	Vel[i] = LocalVel[0]*x_axis[i] + LocalVel[1]*y_axis[i] + LocalVel[2]*z_axis[i];
    }

    /*now get the outgoing direction of the photon, 
      taking advantage of the vector basis I have just generated
      Here I just take the value of dijkstra for the wing.
    */

    /*first define if it's in the core or not*/    

    x_corewing = 1.59 - 0.60*log10(a) - 0.03*log10(a)*log10(a);

    R_1 = RandFloatUnit();

    iso = RandFloatUnit();
    if(iso<(1.0/3.0)){/*isotropic*/
      mu = (2.0*R_1 - 1.0);
    }else{
      if(fabs(x)<x_corewing){			/*In the core*/
	/*now we make the decision if it's isotropic or not*/
	T = (1.0/7.0)*(14.0  - 24.0*R_1  + sqrt(245.0 - 672.0*R_1 + 576*R_1*R_1));	  
	mu = 1.0/(pow(T,1.0/3.0)) - pow(T, 1.0/3.0);	
      }else{
	T = 2.0 - 4.0*R_1  + sqrt(5.0 -16.0*R_1 + 16*R_1*R_1);
	mu = 1.0/(pow(T,1.0/3.0)) - pow(T, 1.0/3.0);	
      }
    }



    
    RND_pair(&R_1, &R_2);    
    R_3 = R_1*R_1 + R_2*R_2;
    for(i=0;i<3;i++){
	DirOutPhoton[i] = 
	    sqrt((1.0-(mu*mu))/R_3)*R_1*x_axis[i] + 
	    sqrt((1.0-(mu*mu))/R_3)*R_2*y_axis[i] + 
	    mu*z_axis[i];
    }

}
