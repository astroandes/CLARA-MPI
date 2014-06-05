#ifndef VECTOR_H
#define VECTOR_H
void cross_product(double *vec_1, double *vec_2, double *result);
double point_product(double *vec_1, double *vec_2);
void normalize(double *vec);
void PropagateLorentzDirChange(double *Dir, double *BulkVel, int sign);
#endif
