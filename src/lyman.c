#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "struct.h"
#include "lyman.h"
#include "myrand.h"


double LyaTau(void)
/*generates a value of optical depth*/
{
    double tau;
    double tmp;
    tmp = RandFloatUnit();

    tau = -log(tmp);
    return tau;
}

double LyaH(double x, double a)
{
    double z, P, q, H;
    z = 0.0; P=0.0; q=0.0; H=0.0;

    z = (x*x - 0.855)/(x*x + 3.42);
    
    P = 5.674*(z*z*z*z) - 9.207*(z*z*z)  + 4.421*(z*z) + 0.1117*z;

    if(z<=0)
    {
	q = 0.0;
    }
    else
    {
	q = (1.0 + (21.0/(x*x)))*(a/(PI*(x*x + 1.0)))*P;
    }

    H = q*sqrt(PI);
    H = H + exp(-(x*x));
    return H;
}

double LyaCrossSection(double x, double a){
    double sigma;
    double nu_doppler;


    nu_doppler = Lya_nu_line_width_CGS/(2.0*a);
    sigma = 0.4162*sqrt(PI)*(CHARGEELECTRON*CHARGEELECTRON)/(ELECTRONMASS*C_LIGHT*nu_doppler);
    sigma = sigma*LyaH(x,a);
        
    return sigma;
}

void LyaHTest(void)
{
    int i;    
    char FileName[MAX_FILENAME_SIZE];
    FILE * output;
    double  a, x, H, delta_x = 0.01, min_x = -10;
    a  = 0.000479;
    x = 0.1;
    sprintf(FileName, "%s/lya_H_test.dat", All.OutputDir);

    output = fopen(FileName, "w");   
    for(i=0; i<(2.*(fabs(min_x)/delta_x)); i++) {
	x = min_x + delta_x*i;
	H = LyaH(x, a);	
	fprintf(output, "%f %f \n", x,H);	
    }
    fclose(output);
}

