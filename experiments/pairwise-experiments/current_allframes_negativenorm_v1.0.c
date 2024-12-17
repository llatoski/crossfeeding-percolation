/*************************************************************************
*             Hilbert Transform Luis C. F. Latoski              *
*                             V1.0 26/06/2021                            *
*************************************************************************/


/***************************************************************
 *                            INCLUDES                      
 **************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>
#include <complex.h>
#include "mc.h"

#define NCELLS 158
#define RADIUS 5 
#define ULIM 475
#define LLIM -25
#define PRECISION 10000
#define PI 3.14159265359
/***************************************************************
 *                            FUNCTIONS                       
 **************************************************************/

/***************************************************************
 *                         GLOBAL VARIABLES                   
 **************************************************************/

double *x,*y,*r,*flux,*fluxerr,normalization=0,**multiplier;
double **current,*norm,**dist;
int *absorb,*emit,**cont;
unsigned long seed;

int main(void)  {

    seed = time(0);
    if (seed%2==0) ++seed;
    start_randomic(seed);
    x = (double*)malloc(NCELLS*sizeof(double));
    y = (double*)malloc(NCELLS*sizeof(double));
    r = (double*)malloc(NCELLS*sizeof(double));
    flux = (double*)malloc(NCELLS*sizeof(double));
    fluxerr = (double*)malloc(NCELLS*sizeof(double));
    absorb = (int*)malloc(NCELLS*sizeof(int));
    emit = (int*)malloc(NCELLS*sizeof(int));
    cont = (int**)malloc(NCELLS*sizeof(int*));
    current = (double**)malloc(NCELLS*sizeof(double*));
    dist = (double**)malloc(NCELLS*sizeof(double*));
    norm = (double*)malloc(NCELLS*sizeof(double));

    int nemit=0,nabsorb=0;
    for(int i=0; i<NCELLS; i++){
        scanf("%lf %lf %lf",&x[i],&y[i],&flux[i]);
        flux[i]=flux[i]/7;
        if(flux[i]<=0){
            absorb[nabsorb]=i;
            nabsorb++;
        }
        else{
            emit[nemit]=i;
            nemit++;
        }
        dist[i] = (double*)malloc(NCELLS*sizeof(double));
        current[i] = (double*)malloc(NCELLS*sizeof(double));
    }   

    for(int j=0; j<nabsorb; j++){
        int cellj = absorb[j];
        for(int i=0; i<nemit; i++){
            int celli = emit[i];
            dist[celli][cellj]=sqrt((x[cellj]-x[celli])*(x[cellj]-x[celli]) + (y[cellj]-y[celli])*(y[cellj]-y[celli]));
            norm[cellj]+=flux[celli]/dist[celli][cellj];
        }
    }

    //Current calculations
    for(int i=0; i<nemit; i++){
        int celli = emit[i];
        for(int j=0; j<nabsorb; j++){
            int cellj = absorb[j];
            current[celli][cellj]=-(flux[celli]*flux[cellj])/(norm[cellj]*dist[celli][cellj]);
            printf("%f ", current[celli][cellj]);
        }
        printf("\n");
    }
}