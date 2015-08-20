#ifndef STRUCT_H
#define STRUCT_H


#include "mtwist.h"
#define MAX_FILENAME_SIZE  1024
#define MAX_VEL_ITER       1000000000
#define MAX_ITER           1000000000
#define ABSORBED            -1
#define OUT_OF_BOX           0
#define ACTIVE               1
#define SATURATED_ITERATIONS 100000
#define N_POINTS_IN_TEST     1000000

/*some units in (cgs) sistem*/
#define G_GRAVITY         6.672e-8
#define HUBBLE            3.2407789e-18 /* in h/sec */
#define HUBBLE_TIME       3.09e+17   /*in sec/h*/
#define C_LIGHT           2.9979e+10   /*cm/s*/
#define CM_PER_MPC        3.085678e24
#define HUBBLE            3.2407789e-18	/* in h/sec */
#define LOG10_L_SUN       22.583     /*log10(3.83E22) = L_bol_sun [in ergs/s] */
#define PI                3.14159265 
#define BOLTZMANN         1.3806e-16 /*cgs units*/
#define PROTONMASS        1.6726e-24 /*cgs units*/
#define GAMMA             1.66666 
#define CHARGEELECTRON	  4.8032E-10  /*E.S.E			H.Scheffler/Elsässer Bau und Physik der Galaxies*/
#define ELECTRONMASS	  9.109382616E-28			/*g				wikipedia*/
#define PLANCK            6.626075E-27 

/*Lyman alpha related constants*/
#define Lya_nu_center_CGS 2.466E15
#define Lya_nu_line_width_CGS 9.936E7 
#define Lya_lambda_line_center_CGS 121.6E-7 /* cm wikipedia */
#define CONSTANT_NU_DOPPLER 1.057E11


typedef struct global_setup
{
  /*input files and formats*/
    char InputDir[MAX_FILENAME_SIZE];  /* directory where all base cubes are*/
    char CubeName[MAX_FILENAME_SIZE];  /* base name of all the cubes*/
    
    /*output files*/
    char OutputDir[MAX_FILENAME_SIZE];            /* directory where all the outputs are written*/
    char OutputFile[MAX_FILENAME_SIZE];            /* file were the outputs are written*/
    char OutputTestFile[MAX_FILENAME_SIZE];            /* file were the tests are written*/
    
    /*parameters controling the algorithm*/
    double LuminosityPerPackage;
    double EffectiveEta;
    double ThresholdCube;
    double EffectiveDust;

    double VmaxSphere;
    double VmaxRotation;
    double VmaxOutflow;

    /*parameters for the dust model*/
    double GrainSize;         /*in cm*/
    double NumberDensityDust;
    double TauDust;
    double DustAbsorptionProb;   /*when interacting with dust the probability to be absorbed*/
    int UseDust;
    int UseVelocities;
    int UseAtomSpeedUp;
    
    /*define the problem*/
    int NeufeldSlab;
    int NeufeldCube;
    int SimulationCube;
    int ExpandingSphere;
    int RotatingSphere;
    int TestParallelVel;
    int TestParallelVelFast;
    int TestFirstScatter;
    int TestPerpVel;
    int TestRND;
    int OutflowRotationSphere;
    int HomogeneousInit; /*specifies if the photons are to be homogeneously distributed over the volume*/

    /*Define some physical characteristics of the problem*/
    double Temperature;
    double Tau;
    double InputFrequency;
    double TotalLuminosity;
    double NumberDensityHI;
  double v_thermal;

    /*units*/
    double UnitMass_in_g;  
    double UnitVelocity_in_cm_per_s;
    double UnitLength_in_cm;
    double UnitLymanLuminosity;
    double UnitTime_in_s;
    double UnitDensity_in_cgs;
    double UnitEnergy_in_cgs;
    
    /*output options*/
    int OutputInitList;
    int OutputFinalList;
    int OutputBinary;

  int RandomSeed;
    /*parameters to be used in case of test*/
    double Test_a;
    double Test_x;

    /*quantities to be calculated afterwards*/
    double SlabLength;

} setup;


typedef struct grid_str
{    
  int N_cells;
  float *cargo;
  float cell_size;
  float x_level;
  float y_level;
  float z_level;
  int N_grid_x; /*Size of the grid*/
  int N_grid_y;
  int N_grid_z;
  float redshift;
  
} grid;


typedef struct lyman_RT_photons_str
{
  /*All the arrays have size N_photons*/
  int  N_photons;    

  /*geometrical and physical properties*/
  int     *Cell;    /* marks the [i,j,k] where the cube is (it's a 3D array)*/  
  double   *x_out;    /* output frequency */
  double   *Pos;     /* marks the position*/
  double   *Dir;     /* last recorded direction of movement along x (it's a 3D array)*/
  double   *Intensity; /* Gives a weight to the intensity by dust absorption*/
  
  /*Boookkeeping data*/
  int *Active;        /* signals if the photon is still active (i.e. inside the box and not absorbed)*/
  int *ScatterHI;    /* Number of scatters with HI*/
  int *ScatterDust;  /* Number of scatters with dust particles*/
  int *Wrong;        /* Number of scatters where the approximation is wrongly used*/
    
  /*parent grid data*/
  int N_grid_x;
  int N_grid_y;
  int N_grid_z;
  
} lyman_RT_photons;



extern double G_INTERNAL_UNITS, LENGTH_INTERNAL_UNITS, MASS_INTERNAL_UNITS, TIME_INTERNAL_UNITS, C_INTERNAL_UNITS;
extern mt_state RND_MT_State;
extern setup All;
extern int ThisProc, SizeProc;


#endif
