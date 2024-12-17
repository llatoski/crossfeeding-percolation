/*************************************************************************
*             Flux Network Reconstruction using RW method                *
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

/***************************************************************
 *                            FUNCTIONS                       
 **************************************************************/

void initialize(void);
void errorpropagation(void);
void generate_stochastic_network(void);
/***************************************************************
 *                         GLOBAL VARIABLES                   
 **************************************************************/

double *x,*y,*emit,*absorb,*flux,*fluxerr,**dist,*norm,**properror;
double **current;
unsigned long seed;

void main(void)  {
    seed=SEED;
    start_randomic(seed);

    initialize();
    errorpropagation();
    generate_stochastic_network();
        
}

/*######################################## 
               Functions
#########################################*/

void initialize(void){
    int NCELLS=NEMIT+NABSORB;
    
    x = (double*)malloc(NCELLS*sizeof(double));
    y = (double*)malloc(NCELLS*sizeof(double));
    flux = (double*)malloc(NCELLS*sizeof(double));
    fluxerr = (double*)malloc(NCELLS*sizeof(double));
    absorb = (double*)malloc(NABSORB*sizeof(double));
    emit = (double*)malloc(NEMIT*sizeof(double));
    norm = (double*)malloc(NCELLS*sizeof(double));
    current = (double**)malloc(NEMIT*sizeof(double*));
    dist = (double**)malloc(NEMIT*sizeof(double*));
    properror = (double**)malloc(NEMIT*sizeof(double*));

    for(int i=0; i<NEMIT; i++){
        current[i] = (double*)malloc(NABSORB*sizeof(double));
        dist[i] = (double*)malloc(NABSORB*sizeof(double));
        properror[i] = (double*)malloc(NABSORB*sizeof(double));
        for(int j=0; j<NABSORB; j++)scanf("%lf", &current[i][j]);
    }
    int nemit=0,nabsorb=0;
    for(int i=0; i<NCELLS; i++){
        scanf("%lf %lf",&x[i],&y[i]);
        scanf("%lf",&flux[i]);
        scanf("%lf",&fluxerr[i]);
        if(flux[i]>0){
            emit[nemit]=i;
            nemit++;
        }
        else{
            absorb[nabsorb]=i;
            nabsorb++;
        }
    }
}

void errorpropagation(void){
    for(int j=0; j<NABSORB; j++){
        int cellj = absorb[j];
        for(int i=0; i<NEMIT; i++){
            int celli = emit[i];
            dist[i][j]=sqrt((x[cellj]-x[celli])*(x[cellj]-x[celli]) + (y[cellj]-y[celli])*(y[cellj]-y[celli]));
            norm[cellj]+=flux[celli]/dist[i][j];
        }
    }
    for(int i=0; i<NEMIT; i++){
        int celli = emit[i];
        for(int j=0; j<NABSORB; j++){
            double inc = 0;
            int cellj = absorb[j];
            for(int k=0; k<NEMIT; k++){
                int cellk = emit[k];
                if(cellk!=celli)inc += pow(((current[i][j]*fluxerr[cellk])/dist[k][j]),2);
            }
            inc += pow(((current[i][j] + flux[cellj])*fluxerr[celli]/dist[i][j]),2);
            inc += pow((flux[celli]*fluxerr[cellj]/dist[i][j]),2);
            properror[i][j]=(double)sqrt(inc)/norm[cellj];
            // if(current[i][j]>=THRESHOLD)if(celli==CELL1 && cellj==CELL2)printf("%d %d %d %.8f %.8f\n",FRAME,celli,cellj,current[i][j],properror[i][j]);
            // if(celli==CELL2 && cellj==CELL1)printf("%d %d %d %.8f %.8f\n",FRAME,celli,cellj,current[i][j],properror[i][j]);
            // if(current[i][j]>=THRESHOLD)printf("%d %d %d %.8f %.8f\n",FRAME,celli,cellj,current[i][j],properror[i][j]);
        }
    }
}

void printview(void){
    for(int i=0; i<NEMIT; i++){ 
        for(int j=0; j<NABSORB; j++){
            if( current[i][j] > THRESHOLD ){
                int celli = emit[i];
                int cellj = absorb[j];
                //if(sqrt(pow((x[celli]-x[cellj]),2) + pow((y[celli]-y[cellj]),2)) < 200.);
                printf("set arrow nohead lw %f from %f,%f to %f,%f front\n",current[celli][cellj],x[celli],y[celli],x[cellj],y[cellj]);
            }
        }
    }
}

void generate_stochastic_network(void){
    for(int i=0; i<NEMIT; i++){ 
        for(int j=0; j<NABSORB; j++){
            double signalnoise = current[i][j]/properror[i][j];
            //printf("%f\n",signalnoise);
            double PROBABILITY = 1 - exp(-signalnoise*DT);
            double RAND = FRANDOM;
            //printf("PROB = %lf, RAND = %lf\n",PROBABILITY,RAND);
            if( RAND <= PROBABILITY && current[i][j] >= THRESHOLD){
                int celli = emit[i];
                int cellj = absorb[j];
                //printf("%d %d %f %f\n",celli,cellj,current[i][j],signalnoise);
                printf("%d %d\n",celli,cellj);
            }
        }
    }
}
