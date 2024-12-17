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
void initialize2(void);
void printview(void);
void printview2(void);
void printview3(void);
void properties(void);
void probabilitymatrix(void);
/***************************************************************
 *                         GLOBAL VARIABLES                   
 **************************************************************/

double *x,*y,*emit,*absorb,*flux;
double **current;

void main(void)  {
   
        initialize2();
        printview2();      
}

/*######################################## 
               Functions
#########################################*/

void initialize2(void){
    int NCELLS=NEMIT+NABSORB;
    x = (double*)malloc(NCELLS*sizeof(double));
    y = (double*)malloc(NCELLS*sizeof(double));
    flux = (double*)malloc(NCELLS*sizeof(double));
    absorb = (double*)malloc(NABSORB*sizeof(double));
    emit = (double*)malloc(NEMIT*sizeof(double));
    current = (double**)malloc(NEMIT*sizeof(double*));
    for(int i=0; i<NEMIT; i++){
        current[i] = (double*)malloc(NABSORB*sizeof(double));
        for(int j=0; j<NABSORB; j++)scanf("%lf", &current[i][j]);
    }
    int nemit=0,nabsorb=0;
    for(int i=0; i<NCELLS; i++){
        #ifdef TC
            scanf("%lf %lf",&x[i],&y[i]);
            for(int j=0; j<7; j++){
                if(j==FRAME)scanf("%lf",&flux[i]);
                else scanf("%*f");
            }
        #else
            scanf("%lf %lf",&x[i],&y[i]);
            scanf("%lf",&flux[i]);
        #endif
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

void printview2(void){
    for(int i=0; i<NEMIT; i++){ 
        for(int j=0; j<NABSORB; j++){
            if( current[i][j] > THRESHOLD ){
                int celli = emit[i];
                int cellj = absorb[j];
                // printf("Aqui tem corrente na linha %d coluna %d\n",i,j);
                printf("%d %d\n",celli,cellj);//,current[i][j]);
            }
        }
    }
}