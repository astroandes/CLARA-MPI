#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define USE "./dijkstra.x tau temperature x_in"
#define N_POINTS  10000
#define PI  3.14159
int main(int argc, char **argv){
    float tau;
    float temperature;
    float x_arr[N_POINTS];
    float neufeld_arr[N_POINTS];
    float min_x, max_x, delta_x;
    int i;
    float Lya_nu_line_width_CGS = 9.936E7;
    float nu_doppler;
    float a;
    float x_in;
    float tmp;
    float total;
    if(argc!=4){
	fprintf(stderr, "%s\n", USE);
	exit(1);
    }

    /*read the inputs*/
    tau = atof(argv[1]);
    temperature = atof(argv[2]);
    x_in = atof(argv[3]);

    
    nu_doppler = 1.057E11*sqrt(temperature/10000.0); 
    a =  Lya_nu_line_width_CGS/(2.0*nu_doppler);                                                          
    
    min_x = -150.0;
    max_x = 150.0;
    delta_x = (max_x - min_x)/N_POINTS;
    for(i=0;i<N_POINTS;i++){
        x_arr[i] = min_x + i*delta_x;
    }

    fprintf(stdout, "# %d\n", N_POINTS);
    total = 0.0;
    for(i=0;i<N_POINTS;i++){
        neufeld_arr[i] = x_arr[i]*x_arr[i];
        neufeld_arr[i] *= sqrt(PI)/(4.0*sqrt(6)*a*tau);
	tmp  = sqrt(2.0*PI*PI*PI/27.0);
	tmp *= (fabs(pow(x_arr[i],3.0) - pow(x_in,3.0))/(a*tau));
        neufeld_arr[i] *= 1.0/(1.0+cosh(tmp));
	//	fprintf(stdout, "%e %e\n", x_arr[i], neufeld_arr[i]);
	total += neufeld_arr[i]*delta_x;
    }    
    for(i=0;i<N_POINTS;i++){
      neufeld_arr[i] *= 1.0/(total*4.0*PI);
    }

    fprintf(stderr, "# total %e [1/(4pi) =%e]\n", total, 1.0/(4.0*PI));
    total= 0.0;
    for(i=0;i<N_POINTS;i++){
      total += neufeld_arr[i]*delta_x;
	fprintf(stdout, "%e %e\n", x_arr[i], neufeld_arr[i]);
    }    
    fprintf(stderr, "# new total %e [1/(4pi) =%e]\n", total, 1.0/(4.0*PI));
    return  0;    
}
