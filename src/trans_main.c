#include <stdio.h>
#include <string.h>
#include "struct.h"
#include "grid.h"
#include "io.h"
#define USE "./translate.c fieldname fileout filein [filein2]"

int main(int argc, char **argv){
    char FileNameIn[MAX_FILENAME_SIZE];
    char FileNameInAux[MAX_FILENAME_SIZE];
    char FileNameOut[MAX_FILENAME_SIZE];
    char FieldName[30];
    grid *G;
    grid *G_aux;
    grid *G_new;
    float min_val, max_val;
    float Temp;
    double unit_l, unit_d, unit_t;
    double tconv, dconv;
    int i;

    /**/
    unit_l      =  0.439345258038073E+26;
    unit_d      =  0.346996787056857E-27;
    unit_t      =  0.175476585435706E+17;
    tconv=PROTONMASS * pow((1e-2 * unit_l/unit_t), 2.0)/ BOLTZMANN;
    dconv=(unit_d*1e-3 / PROTONMASS);
    fprintf(stdout, "tconv %e\n", tconv);

    if(argc>5){
	fprintf(stderr, "%s\n", USE);
	exit(1);
    }

    /*initialize all*/
    strcpy(FileNameIn, argv[3]);
    strcpy(FileNameOut, argv[2]);
    sprintf(FieldName, "%s", argv[1]);
    
    fprintf(stdout, "FileIn %s\n", FileNameIn);
    fprintf(stdout, "FileOut %s\n", FileNameOut);
    fprintf(stdout, "FieldName %s\n", FieldName);


    /*create the grid*/
    G = GridCreate();
    G_new = GridCreate();
    G_aux = GridCreate();
    
    /*load the original data in Pierre's format*/
    G->cargo = LoadField(FileNameIn,G);
    G_new->cargo = LoadField(FileNameIn,G_new);


    if(strcmp(FieldName, "GasTemp")==0 || 
       strcmp(FieldName, "GasTauHI")==0){
	G_aux->cargo = LoadField(argv[4],G_aux);
    }
    




    /*Make the new data cube*/
    for(i=0;i<G->N_cells;i++){
	if(strcmp(FieldName, "GasTemp")==0){
	    G_new->cargo[i] = G->cargo[i]/G_aux->cargo[i];
	    G_new->cargo[i] *= tconv/dconv;
	}

	if(strcmp(FieldName, "GasTauHI")==0){
	    Temp = G->cargo[i]/G_aux->cargo[i];
	    G_new->cargo[i] *= tconv/dconv;
	}

	if(strcmp(FieldName, "GasTauDust")==0){

	}

	if(strcmp(FieldName, "GasVelX")==0){

	}

	if(strcmp(FieldName, "GasVelY")==0){

	}

	if(strcmp(FieldName, "GasVelZ")==0){

	}
    }
    print_min_max(G_new->cargo, G->N_cells, 1, FieldName, &min_val, &max_val);

    return 0;
}


