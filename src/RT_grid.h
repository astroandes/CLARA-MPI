#ifndef RT_GRID_H
#define RT_GRID_H
lyman_RT_record * AllocateHistory(grid *G);
void TracePhotonsFromHistory(void);
void RecordPhotons(void);
void NullHistory(lyman_RT_record *H);
void RecordHistoryPhoton(float *pos_x, float *pos_y, float *pos_z, float *d_x, float *d_y, float *d_z, int *n_scatter, 
			 grid *GasTauHI,     lyman_RT_record *H);
 void  InitializePhotonList(lyman_RT_stars *S, grid *GasVelX, grid *GasVelY, grid *GasVelZ, grid *GasTemp, 
			    lyman_RT_photons *PhList, int n_procs);
void RT_LoadStars(lyman_RT_stars *S);
void RT_LoadGrids(grid * GasTauDust, grid *GasTauHI, 
		  grid *GasTemp, grid *GasVelX, grid *GasVelY, grid *GasVelZ);
int PhotonNumber(lyman_RT_stars *S);
void TracePhotons(void);
int InsideBox(float x, float y, float z, grid *G);
float EscapeFraction(float a, float tau_HI, float effective_theta, float tau_dust);
void PropagatePhoton(float *pos_x, float *pos_y, float *pos_z, float *d_x, float *d_y, float *d_z, float *x, float *I,
		     int *nscatterHI, int *ncatterDust, int *nwrong, 
		     grid *GasTauHI, grid *GasTauDust, grid *GasTemp,
		     grid *GasVelX, grid *GasVelY, grid *GasVelZ);
void RT_LymanCube(float *x_value, float *Intensity, float tau_HI, float temp, float tau_dust, 
		  float *dir_x, float *dir_y, float *dir_z, float vel_x, float vel_y, float vel_z, float cube_size);
void RT_DustAbsorption(float *intensity, float tau_HI, float temp, float tau_dust);
float RT_LocalTau(float x_value, float tau_HI, float temp, 
		  float dir_x, float dir_y, float dir_z,
		  float vel_x, float vel_y, float vel_z, 
		  float delta_step, float cell_size);
void RT_LymanSimple(double *Dir, double *DeltaPos, double *x, double a, double tau, double effective_eta);
void TracePhotonDirty(lyman_RT_record *R, float *pos_x, float *pos_y, float *pos_z, float *d_x, float *d_y, float *d_z, 
		      float *total_dust_escape, float *total_H_tau, int *nscatterHI,
		      grid *GasTauHI, grid *GasTauDust, grid *GasTemp,
		      grid *GasVelX, grid *GasVelY, grid *GasVelZ);
void TracePhotonQuick(lyman_RT_record *R, float *pos_x, float *pos_y, float *pos_z, float *d_x, float *d_y, float *d_z, 
		      float *total_dust_escape, float *x, int *nscatterHI,
		      grid *GasTauHI, grid *GasTauDust, grid *GasTemp,
		      grid *GasVelX, grid *GasVelY, grid *GasVelZ);

#endif
