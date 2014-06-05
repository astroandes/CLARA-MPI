#ifndef RND_LYMAN_H
#define RND_LYMAN_H
void RND_pair(double *r_1, double *r_2);
void RND_lyman_perp_vel_test(void);
void RND_lyman_parallel_vel_test(double x, double a);
void RND_lyman_parallel_vel_fast_test(double x, double a);
void RND_spherical(double *vec);
double RND_lyman_parallel_vel(double x, double a);
double RND_lyman_parallel_vel_fast(double my_x, double a);
void RND_lyman_perp_vel(double *u_1, double *u_2);
void RND_lyman_atom(double *Vel, double *DirPhoton, 
		    double *DirOutPhoton, double x, double a);
#endif 

