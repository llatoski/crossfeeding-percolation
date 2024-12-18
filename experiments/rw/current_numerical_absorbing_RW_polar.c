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
 *                            DEFINITIONS
 **************************************************************/
#define NCELLS 158 //Number of cells in the culture
#define RADIUS 5 //Cell radius
#define ULIM 475 //Upper limit of the lattice
#define LLIM -25 //Lower limit of the lattice
#define PRECISION 10000 //Number of walkers emitted for each emitting cell
#define PI 3.14159265359 
#define STEP 10 //Radius of movement of each random walk

/***************************************************************
 *                            FUNCTIONS                       
 **************************************************************/
void readfile(void);
void evolution(void);
void normalization_calc(void);
void estimates(void);

/***************************************************************
 *                         GLOBAL VARIABLES                   
 **************************************************************/
double *x,*y,*r,*flux,*fluxerr,normalization=0,**multiplier;
double **current,**dist;
int *absorb,*emit,*norm,**cont, nemit, nabsorb;
unsigned long seed;

/***************************************************************
 *                           MAIN CODE                    
 **************************************************************/
int main(void)  {

    seed = time(0);
    if (seed%2==0) ++seed;
    start_randomic(seed);

    readfile();
    evolution();
    normalization_calc();
    estimates();

}

/***************************************************************
 *        Read file containing cells position and flux                    
 **************************************************************/
void readfile(){
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
}

/***************************************************************
 *                  Random walk evolution
 *    ( PRECISION walkers are released for each emitter )
 **************************************************************/    
void evolution(){
    double dt=10;
    for(int i=0; i<nemit; i++){
        int celli = emit[i];
        for(int k=0; k<PRECISION; k++){ //Initialize the walker at the center of a cell and starts the evolution
            double xwalk=x[celli];
            double ywalk=y[celli];
            int step=0;
            int found=0;
            int alive=1;
            while(alive==1){ //While the walker is inside the frame and unarbsobed
                double theta = FRANDOM*(2*PI); //Move the walker
                xwalk+=STEP*sin(theta);
                ywalk+=STEP*cos(theta);
                if(xwalk>ULIM | ywalk>ULIM | xwalk<LLIM | ywalk<LLIM){//If the walker leaves the frame a new emission begins
                    alive=0;
                }
                else{ //Test for absorbtion
                    for(int j=0; j<nabsorb; j++){
                        int cellj = absorb[j];
                        double relativedist = sqrt( (pow((xwalk-x[cellj]),2)) + (pow((ywalk-y[cellj]),2)) );
                        double RAND = (double)FRANDOM;
                        double probability = (double)flux[celli]*dt/normalization; //Absorbption probability is proportional to the emiting flux and normalized
                        if( ( relativedist<RADIUS && RAND<probability) ){ //If inside the range of an absorbing cell, test for for absorbtion
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
}

/***************************************************************
 *              Calculate normalization factor       
 *      (number of walkers absorbed at a given position)
 **************************************************************/
void normalization_calc(){
    for(int j=0; j<nabsorb; j++){
        int celli=absorb[j];
        for(int k=0; k<nemit; k++){
            int cellj=emit[k];
            norm[celli]+=cont[cellj][celli];
        }
    }
}

/***************************************************************
 *         Current estimates using random walk method
 **************************************************************/
void estimates(){
    for(int j=0; j<nemit; j++){
        int celli=emit[j];
        for(int k=0; k<nabsorb; k++){
            int cellj=absorb[k];
            printf("%f ",-(double)cont[celli][cellj]*flux[cellj]/norm[cellj]);
        }
        printf("\n");
    }
}
