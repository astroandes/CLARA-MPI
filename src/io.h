#ifndef IO_H
#define IO_H
void SavePhotonListAscii(char *fname, lyman_RT_photons *Ph);
void OpenPhotonListAscii(char *fname, lyman_RT_photons *Ph);
void AppendPhotonListAscii(char *fname, lyman_RT_photons *Ph, int index);
#endif
