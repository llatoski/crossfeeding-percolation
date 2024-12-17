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

#define PRECISION 1e6
//#define RADIUS 5. 


/***************************************************************
 *                            FUNCTIONS                       
 **************************************************************/

void openfiles(void);
void initialize(void);
int test_absorption(int,int);
void matrixview(void);
void evolution(int);
void change_box(int*);

/***************************************************************
 *                         GLOBAL VARIABLES                   
 **************************************************************/

FILE *fp2;
double *EmiterBox, *EmiterPa;
double *AbsorbPe,*AbsorbFlux,*AbsorbNNeigh,**AbsorbNeighborhood,**AbsorbNeighborhoodProbs;
double **current,**AbsorbingProbability;
int **cont;
unsigned long seed;

int main(void)  {

    initialize();   
    
    int where[NEMIT];    
    int k=0;
    for(int i=0; i<NEMIT; i++)if(EmiterPa[i]>0.025){
        where[k]=i;
        k++;
    }

    for(int j=0; j<k; j++){
        int i = where[j];
        printf("Running for emitting cell %d with Pa = %f\n",i,EmiterPa[i]);
        for(int k=0; k<PRECISION; k++){
            evolution(i);
        }
    }

    // /*Generating hit matrix*/
    // for(int i=0; i<NEMIT; i++){
    //         printf("Running for emitting cell %d with Pa = %f\n",i,EmiterPa[i]);
    //         for(int k=0; k<PRECISION*EmiterPa[i]; k++){
    //             evolution(i);
    //         }
    // }

 

    /*Estimating currents*/
    for(int i=0; i<NABSORB; i++){
        int hits=0;
        for(int j=0; j<NEMIT; j++){
            hits+=cont[j][i];
        }
        for(int j=0; j<NEMIT; j++){
            if(hits!=0)current[j][i] = -(double)cont[j][i]*AbsorbFlux[i]/hits;
        }       
    }

    matrixview();

}

/*######################################## 
               Functions
#########################################*/

void initialize(void){
    
    #ifdef DEBUG
        seed = 1;
    #else
        seed = time(0);
        if (seed%2==0) ++seed;
    #endif
    start_randomic(seed);

    EmiterBox = (double*)malloc(NEMIT*sizeof(double));
    EmiterPa = (double*)malloc(NEMIT*sizeof(double));
    AbsorbPe = (double*)malloc(NABSORB*sizeof(double));
    AbsorbNNeigh = (double*)malloc(NABSORB*sizeof(double));
    AbsorbFlux = (double*)malloc(NABSORB*sizeof(double));
    for(int i=0; i<NEMIT; i++)scanf("%*f %lf %lf",&EmiterBox[i],&EmiterPa[i]);
    for(int i=0; i<NABSORB; i++)scanf("%*f %lf %lf %lf",&AbsorbNNeigh[i],&AbsorbPe[i],&AbsorbFlux[i]);
    AbsorbNeighborhood = (double**)malloc(NABSORB*sizeof(double*));
    AbsorbNeighborhoodProbs = (double**)malloc(NABSORB*sizeof(double*));
    AbsorbingProbability = (double**)malloc(NEMIT*sizeof(double*));
    cont = (int**)malloc(NEMIT*sizeof(int*));
    current = (double**)malloc(NEMIT*sizeof(double*));     
    for(int i=0; i<NABSORB; i++){
        AbsorbNeighborhood[i]=(double*)malloc(AbsorbNNeigh[i]*sizeof(double));
        AbsorbNeighborhoodProbs[i]=(double*)malloc(AbsorbNNeigh[i]*sizeof(double));
        for(int j=0; j<AbsorbNNeigh[i]; j++)scanf("%*f %lf %lf",&AbsorbNeighborhood[i][j],&AbsorbNeighborhoodProbs[i][j]);
    }
    for(int i=0; i<NEMIT; i++){
        cont[i]=(int*)malloc(NABSORB*sizeof(int));
        current[i]=(double*)malloc(NABSORB*sizeof(double));
        AbsorbingProbability[i] = (double*)malloc(NABSORB*sizeof(double));
        double Pa = EmiterPa[i];
        for(int j=0; j<NABSORB; j++){
            double Pe = AbsorbPe[j];
            AbsorbingProbability[i][j]=Pa*(1-Pe)/(1-(1-Pa)*(1-Pe));
	    //            printf("%d %d %f %f\n",i,j,Pa*(1-Pe),AbsorbingProbability[i][j]);
        }
    }
    openfiles();
}

int test_absorption(int emitter, int absorber){
    double rand = FRANDOM;
    double probability = AbsorbingProbability[emitter][absorber];
    if(rand<probability){
        return 0;
    }
    else {
        return 1;
    }
}

void evolution(int celli){
    int alive=1;
    int currentbox = EmiterBox[celli];
    do{
        alive=test_absorption(celli,currentbox);
        cont[celli][currentbox]+=(1-alive);
        if(alive==1){
            int teste=currentbox;
            change_box(&teste);
            while(teste==-1 & alive==1){
                alive=test_absorption(celli,currentbox);
                cont[celli][currentbox]+=(1-alive);
                teste=currentbox;
                change_box(&teste);
            }
            currentbox=teste;
        }
    }while(alive==1);
}

void change_box(int *r){
    double rand = FRANDOM;
    int k=0;
    double probability=0;
    int newbox = *r;
    //printf("New exchange\n");
    do{
        probability+=AbsorbNeighborhoodProbs[*r][k];
        //printf("Accumulated probability: %f\n",probability);
        if(rand>probability)newbox=AbsorbNeighborhood[*r][k];
        k++;
    }while(k<AbsorbNNeigh[*r]);
    *r=newbox;
}


void matrixview(void){
    for(int i=0; i<NEMIT; i++){
        for(int j=0; j<NABSORB; j++){
            fprintf(fp2,"%f ",current[i][j]);
        }
        fprintf(fp2,"\n");
    }
}

void openfiles(void) {
    char output_file[300];
    char prefix[250];

    sprintf(prefix,"networks/frame%d",FRAME);
    sprintf(output_file,"%s_Matrix.dat",prefix);
    fp2 = fopen(output_file,"w");
    fflush(fp2);
    
    return;
}

