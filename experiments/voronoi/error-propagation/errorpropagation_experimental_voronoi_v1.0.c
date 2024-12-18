/*
    Supporting material software to the article 
    Cross-feeding percolation phase transitions of inter-cellular metabolic networks
    
    Copyright (C) 2064 L.C.F. Latoski, D.De Martino, A.De Martino

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    a with this program.  If not, see <http://www.gnu.org/licenses/>.
                                                                        */


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

/***************************************************************
 *                           MAIN CODE                   
 **************************************************************/

void main(void)  {
    seed=SEED;
    start_randomic(seed);

    initialize();
    errorpropagation();
    generate_stochastic_network();
        
}

/***************************************************************
 *              Initialization and File reading                   
 **************************************************************/
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
    for(int i=0; i<NCELLS; i++){//Read auxiliar file
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

/***************************************************************************
    Propagate error using analytical expression (pairwise approximation)
****************************************************************************/
void errorpropagation(void){
    for(int j=0; j<NABSORB; j++){//Calculate distances and normalization factor
        int cellj = absorb[j];
        for(int i=0; i<NEMIT; i++){
            int celli = emit[i];
            dist[i][j]=sqrt((x[cellj]-x[celli])*(x[cellj]-x[celli]) + (y[cellj]-y[celli])*(y[cellj]-y[celli]));
            norm[cellj]+=flux[celli]/dist[i][j];
        }
    }
    for(int i=0; i<NEMIT; i++){//Error propagation (Eq.28 of the manuscript)
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
        }
    }
}

/************************************************************************
          Generate stochastic network based on signal/noise ratio
*************************************************************************/
void generate_stochastic_network(void){
    for(int i=0; i<NEMIT; i++){ 
        for(int j=0; j<NABSORB; j++){
            double signalnoise = current[i][j]/properror[i][j];
            double PROBABILITY = 1 - exp(-signalnoise*DT); //Estimate probability based on signal noise ratio
            double RAND = FRANDOM; //Generate random number
            if( RAND <= PROBABILITY && current[i][j] >= THRESHOLD){ //If random is smaller than probability and current above threshold we say that there is a link between the two cells
                int celli = emit[i];
                int cellj = absorb[j];
                printf("%d %d\n",celli,cellj);
            }
        }
    }
}
