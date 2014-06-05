#include <stdio.h>
#include <stdlib.h>
#include "struct.h"
#include "mtwist.h"
#include "myrand.h"

double RandFloatUnit(void){
    double ll;
    int dat;
    ll = mts_ldrand(&RND_MT_State);   
    return ll;
}
