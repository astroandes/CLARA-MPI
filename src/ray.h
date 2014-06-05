#ifndef RAY_H
#define RAY_H
lyman_RT_photons * PhotonListCreate(int N_packages);
void PhotonFree(lyman_RT_photons *Ph);
void PhotonListInitialize(lyman_RT_photons *Ph);
#endif
