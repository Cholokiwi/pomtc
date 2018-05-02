// Numerator for z_ig in logs  
// (Likelihood as a function of theta, pi, and data)
// Author : Roy Costilla
// Version : Nov2016
// z_ig's formulation is from Eleni and Ivy's (Trickwithlogs.pdf, changed mr for mg)
// z_ig=mg/sum(mg for all row groups)

#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

extern "C" void Zrc_tt(double *theta, int *y, double *pi, int *n, int *p, int *q, int *G, double *Z)
{
  int i, g, j, k,kdash, max, y_prev=0, y_now=0,counter1=0, counter=0;  // Dont forget to initialise all vars in C++!!
  double sum_jk, log_mg;
  // theta coming from R is Rqqp

  for (i=0 ; i<*n ; ++i){
    //cout << " i=" << i << endl;
    for (g=0 ; g<*G ; ++g){
      //cout << " r=" << g << endl;
      sum_jk=0;
      log_mg=0;
      // Transitional models start in 2nd occasion
      for (j=1; j<*p; ++j){
            //cout << " i=" << i << " G=" << *G << " g=" << g << " j=" << j  << " q=" << *q << endl;
            //counter1 = i + j * (*n) ;
            y_prev = y[i + (j-1) * (*n)]   ;  // In C++ counters start at 0!!
            y_now = y[i + j * (*n)]  ;  // In C++ counters start at 0!!
            //counter = (y_prev-1)  + (y_now-1) *(*q) + g * (*q) * (*q);
            counter = g + (y_prev-1)*(*G)  + (y_now-1) * (*G) * (*q) + j*(*q)*(*G)*(*q) ;
            //cout << "j=" << j << " y_prev=" << y_prev << " y_now="<< y_now << " counter=" << counter << " theta=" << theta[counter] << endl;
            //cout << "theta and log=" << theta[counter] << " " <<  log(theta[counter]) << endl;
            if (counter >= (*p)*(*G)*(*q)*(*q) )  cout << " counter=" << counter << " which is higher than R*q*q=" << max << " THETA doesn't make sense!!" << endl;
            sum_jk+= log(theta[counter]); //USE += with no spaces!!
       }
    log_mg = log(pi[g]) + sum_jk;
    Z[g + i*(*G)] =log_mg;
    //cout << " i=" << i  << " g=" << g << " Z_ig=" << Z[g+i*(*G)] << endl;
    }
  }
  
}
