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

#define NCELLS 158
#define RADIUS 5 
#define ULIM 475
#define LLIM -25
#define PRECISION 10000
#define PI 3.14159265359
#define STEP 10
/***************************************************************
 *                            FUNCTIONS                       
 **************************************************************/

/***************************************************************
 *                         GLOBAL VARIABLES                   
 **************************************************************/

double *x,*y,*r,*flux,*fluxerr,normalization=0,**multiplier;
double **current,**dist;
int *absorb,*emit,*norm,**cont;
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
    norm = (int*)malloc(NCELLS*sizeof(int));

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
            normalization+=flux[i];
        }
        cont[i] = (int*)malloc(NCELLS*sizeof(int));
    }   
    
    double dt=10;
    for(int i=0; i<nemit; i++){
        int celli = emit[i];
        printf("Running for emitting cell %d\n",i);
        for(int k=0; k<PRECISION; k++){
            double xwalk=x[celli];
            double ywalk=y[celli];
            // printf("Walker %d leaving from %.2f %.2f\n",k,(double)xwalk/100,(double)ywalk/100);
            int step=0;
            int found=0;
            int alive=1;
            while(alive==1){
                double theta = FRANDOM*(2*PI);
                xwalk+=STEP*sin(theta);
                ywalk+=STEP*cos(theta);
                if(xwalk>ULIM | ywalk>ULIM | xwalk<LLIM | ywalk<LLIM){
                    alive=0;
                    // printf("Saiu do frame na posicao x = %.2f,y=%.2f\n",(double)xwalk,(double)ywalk);
                }
                else{
                    for(int j=0; j<nabsorb; j++){
                        // printf("Testando\n");
                        int cellj = absorb[j];
                        double relativedist = sqrt( (pow((xwalk-x[cellj]),2)) + (pow((ywalk-y[cellj]),2)) );
                        double RAND = (double)FRANDOM;
                        double probability = (double)flux[celli]*dt/normalization;
                        // if( relativedist < RADIUS)printf("Probabilidade %f compara com o n aleatorio %f\n",probability,RAND);
                        if( ( relativedist<RADIUS && RAND<probability) ){
                            alive=0;
                            found=1;
                            cont[celli][cellj]++;
                            break;
                        }                
                    }
                }
                step++;
            }
        }
    }
    for(int j=0; j<nabsorb; j++){
        int celli=absorb[j];
        for(int k=0; k<nemit; k++){
            int cellj=emit[k];
            norm[celli]+=cont[cellj][celli];
        }
    }
    for(int j=0; j<nemit; j++){
        int celli=emit[j];
        for(int k=0; k<nabsorb; k++){
            int cellj=absorb[k];
            printf("%f ",-(double)cont[celli][cellj]*flux[cellj]/norm[cellj]);
        }
        printf("\n");
    }
}
