#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "struct.h"
#include "parse.h"

/*! This function parses the parameterfile in a simple way.  Each paramater
 *  is defined by a keyword (`tag'), and can be either of type double, int,
 *  or character string.  The routine makes sure that each parameter
 *  appears exactly once in the parameterfile, otherwise error messages are
 *  produced that complain about the missing parameters.
 */
void ReadParameters(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[400];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int  errorFlag = 0;

  
  nt = 0;
  
  strcpy(tag[nt], "InputDir");
  addr[nt] = All.InputDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "CubeName");
  addr[nt] = All.CubeName;
  id[nt++] = STRING;

  strcpy(tag[nt], "OutputDir");
  addr[nt] = All.OutputDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "OutputFile");
  addr[nt] = All.OutputFile;
  id[nt++] = STRING;

  strcpy(tag[nt], "OutputTestFile");
  addr[nt] = All.OutputTestFile;
  id[nt++] = STRING;

  strcpy(tag[nt], "LuminosityPerPackage");
  addr[nt] = &All.LuminosityPerPackage;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "EffectiveEta");
  addr[nt] = &All.EffectiveEta;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "EffectiveDust");
  addr[nt] = &All.EffectiveDust;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "GrainSize");
  addr[nt] = &All.GrainSize;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "TauDust");
  addr[nt] = &All.TauDust;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "DustAbsorptionProb");
  addr[nt] = &All.DustAbsorptionProb;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ThresholdCube");
  addr[nt] = &All.ThresholdCube;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "VmaxSphere");
  addr[nt] = &All.VmaxSphere;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "HomogeneousInit");
  addr[nt] = &All.HomogeneousInit;
  id[nt++] = INT;

  strcpy(tag[nt], "TestParallelVel");
  addr[nt] = &All.TestParallelVel;
  id[nt++] = INT;

  strcpy(tag[nt], "TestParallelVelFast");
  addr[nt] = &All.TestParallelVelFast;
  id[nt++] = INT;

  strcpy(tag[nt], "TestFirstScatter");
  addr[nt] = &All.TestFirstScatter;
  id[nt++] = INT;

  strcpy(tag[nt], "TestRND");
  addr[nt] = &All.TestRND;
  id[nt++] = INT;

  strcpy(tag[nt], "TestPerpVel");
  addr[nt] = &All.TestPerpVel;
  id[nt++] = INT;

  strcpy(tag[nt], "ExpandingSphere");
  addr[nt] = &All.ExpandingSphere;
  id[nt++] = INT;

  strcpy(tag[nt], "RotatingSphere");
  addr[nt] = &All.RotatingSphere;
  id[nt++] = INT;

  strcpy(tag[nt], "UseDust");
  addr[nt] = &All.UseDust;
  id[nt++] = INT;

  strcpy(tag[nt], "UseAtomSpeedUp");
  addr[nt] = &All.UseAtomSpeedUp;
  id[nt++] = INT;

  strcpy(tag[nt], "UseVelocities");
  addr[nt] = &All.UseVelocities;
  id[nt++] = INT;

  strcpy(tag[nt], "NeufeldSlab");
  addr[nt] = &All.NeufeldSlab;
  id[nt++] = INT;

  strcpy(tag[nt], "NeufeldCube");
  addr[nt] = &All.NeufeldCube;
  id[nt++] = INT;

  strcpy(tag[nt], "SimulationCube");
  addr[nt] = &All.SimulationCube;
  id[nt++] = INT;

  strcpy(tag[nt], "Temperature");
  addr[nt] = &All.Temperature;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Tau");
  addr[nt] = &All.Tau;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "NumberDensityHI");
  addr[nt] = &All.NumberDensityHI;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "InputFrequency");
  addr[nt] = &All.InputFrequency;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "TotalLuminosity");
  addr[nt] = &All.TotalLuminosity;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "RandomSeed");
  addr[nt] = &All.RandomSeed;
  id[nt++] = INT;

  strcpy(tag[nt], "OutputInitList");
  addr[nt] = &All.OutputInitList;
  id[nt++] = INT;

  strcpy(tag[nt], "OutputFinalList");
  addr[nt] = &All.OutputFinalList;
  id[nt++] = INT;

  strcpy(tag[nt], "OutputBinary");
  addr[nt] = &All.OutputBinary;
  id[nt++] = INT;

  strcpy(tag[nt], "UnitMass_in_g");
  addr[nt] = &All.UnitMass_in_g;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
  addr[nt] = &All.UnitVelocity_in_cm_per_s;
  id[nt++] = DOUBLE;
  
  strcpy(tag[nt], "UnitLength_in_cm");
  addr[nt] = &All.UnitLength_in_cm;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "UnitLymanLuminosity");
  addr[nt] = &All.UnitLymanLuminosity;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Test_a");
  addr[nt] = &All.Test_a;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Test_x");
  addr[nt] = &All.Test_x;
  id[nt++] = DOUBLE;


  if((fd = fopen(fname, "r")))
  {
      sprintf(buf, "%s%s", fname, "-usedvalues");
      if(!(fdout = fopen(buf, "w")))
      {
	  printf("error opening file '%s' \n", buf);
	  errorFlag = 1;
      }
      else
      {
	  while(!feof(fd))
	  {
	      *buf = 0;
	      fgets(buf, 200, fd);
	      if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
		  continue;
	      
	      if(buf1[0] == '%')
		  continue;
	      
	      for(i = 0, j = -1; i < nt; i++)
		  if(strcmp(buf1, tag[i]) == 0)
		  {
		      j = i;
		      tag[i][0] = 0;
		      break;
		  }
	      
	      if(j >= 0)
	      {
		  switch (id[j])
		  {
		      case DOUBLE:
			  *((double *) addr[j]) = atof(buf2);
			  fprintf(fdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  break;
		      case STRING:
			  strcpy(addr[j], buf2);
			  fprintf(fdout, "%-35s%s\n", buf1, buf2);
			  break;
		      case INT:
			  *((int *) addr[j]) = atoi(buf2);
			  fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  break;
		  }
	      }
	      else
	      {
		  fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
			  fname, buf1);
		  errorFlag = 1;
	      }
	  }
	  fclose(fd);
	  fclose(fdout);
	  
	  i = strlen(All.OutputDir);
	  if(i > 0)
	      if(All.OutputDir[i - 1] != '/')
		  strcat(All.OutputDir, "/");
	  
	  sprintf(buf1, "%s%s", fname, "-usedvalues");
	  sprintf(buf2, "%s%s", All.OutputDir, "parameters-usedvalues");
	  sprintf(buf3, "cp %s %s", buf1, buf2);
	  system(buf3);
      }
  }
  else
  {
      printf("\nParameter file %s not found.\n\n", fname);
      errorFlag = 2;
  }
  
  if(errorFlag != 2)
      for(i = 0; i < nt; i++)
      {
	  if(*tag[i])
	  {
	      printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	      errorFlag = 1;
	  }
      }


  /*set all the units*/  
  All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitDensity_in_cgs =  All.UnitMass_in_g/ pow(All.UnitLength_in_cm,3);
  All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm,2);
  All.UnitEnergy_in_cgs /= pow(All.UnitTime_in_s,2);

}
