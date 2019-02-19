// Numerator for z_ig in logs  
// (Likelihood as a function of theta, pi, and data)
// Author : Roy Costilla
// Version : Feb19

#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

extern "C" void Zrc_tt(double *theta, int *y, double *pi, int *n, int *p, int *q, int *G, double *Z)
{
  int i, g, j, y_prev=0, y_now=0, counter=0; 
  double sum_jk, log_mg;

  for (i=0 ; i<*n ; ++i){
    for (g=0 ; g<*G ; ++g){
      sum_jk=0;
      log_mg=0;
      // Transitional models start in 2nd occasion
      for (j=1; j<*p; ++j){
            y_prev = y[i + (j-1) * (*n)]   ;  
            y_now = y[i + j * (*n)]  ;
            counter = g + (y_prev-1)*(*G)  + (y_now-1) * (*G) * (*q) + j*(*q)*(*G)*(*q) ;
            sum_jk+= log(theta[counter]); 
       }
    log_mg = log(pi[g]) + sum_jk;
    Z[g + i*(*G)] =log_mg;
    }
  }
}
