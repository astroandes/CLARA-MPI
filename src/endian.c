#include <stdio.h>
#include "endian.h"

int fread_sw(void *data,int size,int nb,FILE *f,int swap)
{
    
    fread(data,size,nb,f);

    if (swap)
    {
      switch (size)
	{
	case 8: Dswap8BArr(data,nb);break;
	case 4: Dswap4BArr(data,nb);break;
	case 2: Dswap2BArr(data,nb);break;
	}
    }

    return 0;
}

int swapI(int val)
{
    int out;
    const char *i=(const char *)&val;
    char *o=(char *)&out;

    o[3]=i[0];
    o[2]=i[1];
    o[1]=i[2];
    o[0]=i[3];

    return out; 
}

float swapF(float val)
{
    float out;
    const char *i=(const char *)&val;
    char *o=(char *)&out;

    o[3]=i[0];
    o[2]=i[1];
    o[1]=i[2];
    o[0]=i[3];

    return out; 
}

double swapD(double val)
{
    double out;
    const char *i=(const char *)&val;
    char *o=(char *)&out;

    o[7]=i[0];
    o[6]=i[1];
    o[5]=i[2];
    o[4]=i[3];
    o[3]=i[4];
    o[2]=i[5];
    o[1]=i[6];
    o[0]=i[7];


    return out; 
}

void Dswap2B(void *val)
{
    char *c=(char *)val;
    char a;
    
   a=c[0];c[0]=c[1];c[1]=a; 
}

void Dswap4B(void *val)
{
    char *c=(char *)val;
    char a;
    
   a=c[0];c[0]=c[3];c[3]=a;
   a=c[1];c[1]=c[2];c[2]=a; 
 
}
void Dswap8B(void *val)
{
   char *c=(char *)val;
   char a;
    
   a=c[0];c[0]=c[7];c[7]=a;
   a=c[1];c[1]=c[6];c[6]=a;
   a=c[2];c[2]=c[5];c[5]=a;
   a=c[3];c[3]=c[4];c[4]=a;
}

void Dswap2BArr(void *val,int n)
{
    int i;
    char a;

    char *c=(char *)val;

    for (i=0;i<2*n;i+=2)
    {
	a=c[i];
	c[i]=c[i+1];
	c[i+1]=a;
    }

}


void Dswap4BArr(void *val,int n)
{
    int i;
    char a,b;

    char *c=(char *)val;

    for (i=0;i<4*n;i+=4)
    {
	a=c[i];
	b=c[i+1];
	c[i]=c[i+3];
	c[i+1]=c[i+2];
	c[i+2]=b;
	c[i+3]=a;
    }

}

void Dswap8BArr(void *val,int n)
{
    int i;
    char a,b,u,v;

    char *c=(char *)val;

    for (i=0;i<8*n;i+=8)
    {
	a=c[i];
	b=c[i+1];
	u=c[i+2];
	v=c[i+3];
	c[i]=c[i+7];
	c[i+1]=c[i+6];
	c[i+2]=c[i+5];
	c[i+3]=c[i+4];
	c[i+4]=v;
	c[i+5]=u;
	c[i+6]=b;
	c[i+7]=a;
    }

}
