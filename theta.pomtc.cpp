// POMtt: theta Rpqp (pdf) for the pomtt  (Transitional POM+column effects)
// Author : Roy Costilla
// Version : Nov16

#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

extern "C" void theta_pomtc(double *mu, double *alpha, double *beta, double *gamma, int *R, int *p, int *q, double *thetacdf, double *theta)
{
  int r,k,kdash,j, counter=0;  
// theta.temp <- array(1, dim=c(p,q,q,R))  matrix coming from R

for (r=0 ; r<*R ; ++r){
  for (kdash=0 ; kdash<*q ; ++kdash){
    for (k=0 ; k< *q ; ++k){
      for( j=0 ; j < (*p) ; ++j ) {
        //cout << "j=" << j << endl;
            counter = r + kdash*(*R) + k*(*q)*(*R) + j*(*q)*(*q)*(*R)  ;
             // Calculate thetacdf for k<q (cdf=1 if k=q)
            if ( k < (*q-1) ) {
              //thetacdf[counter] = 1/(1+exp(-(mu[k] - alpha[r] - beta[kdash] - gamma[j] )));
              thetacdf[counter] = 1/(1+exp(-(mu[k] - alpha[r] - beta[r+kdash*(*R)] - gamma[j] )));
              //thetacdf[counter] = 1/(1+exp(-(mu[k] - alpha[r] - beta[r+j*(*R)] )));
              
              
            }
            // Calculate theta pdf
            if (k==0) {
              theta[counter] = thetacdf[counter] ;          
            }
            else {
              // this doesnt changed as it within each j 
              theta[counter] = thetacdf[counter] - thetacdf[counter - (*q)*(*R)] ;
            }
           //cout << " r=" << r+1 << " kdash=" << kdash+1 <<" k=" << k+1 << " j=" << j+1 <<  
          //    " nu=" << mu[k] - alpha[r] - beta[kdash] - gamma[j] <<" theta.cdf=" << thetacdf[counter] <<" theta.pdf="<< theta[counter] << endl;;
      }
    }
  }
 }
}