#include <stdio.h>
#include <math.h>
#include "vector.h"
double point_product(double *vec_1, double *vec_2){
    double point;
    int i;
    point = 0.0;
    for(i=0;i<3;i++){
        point+=vec_1[i]*vec_2[i];
    }
    return point;
}


void cross_product(double *vec_1, double *vec_2, double *result){
    result[0] = vec_1[1]*vec_2[2] - vec_1[2]*vec_2[1];
    result[1] = vec_1[2]*vec_2[0] - vec_1[0]*vec_2[2];
    result[2] = vec_1[0]*vec_2[1] - vec_1[1]*vec_2[0];
}


void normalize(double *vec){
    double norm;
    int i;
    norm = 0.0;
    norm = point_product(vec, vec);
    norm = sqrt(norm);
    for(i=0;i<3;i++){
        vec[i] = vec[i]/norm;
    }
}
