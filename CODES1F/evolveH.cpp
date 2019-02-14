/* evolve states in real-time with the Hamiltonian */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

// Notation: eigenstate |n> = \sum_k \alpha_k |k>, |k> is a cartoon state
void evolve_cartoons(std::vector<double> &evals, std::vector<std::vector<double>> &evecs){
   extern int scan(std::vector<bool>&);
   int ix,iy,parity,p,q1,q2;
   double t;
   double temp1,temp2;
   double ampl1_RE, ampl1_IM, ampl2_RE, ampl2_IM;
   double lam1,lam2;
   std::vector<bool> cart1(2*VOL), cart2(2*VOL);
   std::vector<double> alpha1,alpha2;
   FILE *outf;
   outf = fopen("overlap.dat","w");
   /* construct cartoon state */
   for(iy=0;iy<LY;iy++){
   for(ix=0;ix<LX;ix++){
    parity=(ix+iy)%2;
    p = 2*(iy*LX+ix);
    if(parity){
       cart1[p]=false; cart1[p+1]=true; cart2[p]=true;  cart2[p+1]=false;
    }
    else{
       cart1[p]=true; cart1[p+1]=false; cart2[p]=false;  cart2[p+1]=true;
    }
   }}
   for(p=0;p<NH;p++){
     q1=scan(cart1);
     q2=scan(cart2);
     alpha1.push_back(evecs[p][q1]);
     alpha2.push_back(evecs[p][q2]);
   }
   // check that the coefficients are correctly set
   //for(p=0;p<NH;p++) std::cout<<alpha1[p]<<" "<<alpha2[p]<<std::endl;
   // compute the real-time evolution
   temp1 = 0.0; temp2 = 0.0;
   for(p=0;p<NH;p++){
    if(fabs(evals[p]) < 1e-10){
      temp1 = temp1 + alpha1[p]*alpha1[p];
      temp2 = temp2 + alpha2[p]*alpha1[p];
    }
   }
   //std::cout<<ampl1<<" "<<ampl2<<std::endl;
   for(t=Ti;t<Tf;t=t+dT){
     ampl1_RE=0.0; ampl2_RE=0.0;
     ampl1_IM=0.0; ampl2_IM=0.0;
     for(p=0;p<NH;p++){
        if(fabs(evals[p]) > 1e-10){
          ampl1_RE = ampl1_RE + alpha1[p]*alpha1[p]*cos(evals[p]*t);
          ampl2_RE = ampl2_RE + alpha2[p]*alpha1[p]*cos(evals[p]*t);
          ampl1_IM = ampl1_IM + alpha1[p]*alpha1[p]*sin(evals[p]*t);
          ampl2_IM = ampl2_IM + alpha2[p]*alpha1[p]*sin(evals[p]*t);
        } // close if-loop
     } // close for-loop
     ampl1_RE = ampl1_RE + temp1; ampl2_RE = ampl2_RE + temp2;
     //std::cout<<t<<" "<<ampl1<<std::endl;
     lam1 = -log(ampl1_RE*ampl1_RE + ampl1_IM*ampl1_IM)/VOL;
     lam2 = -log(ampl2_RE*ampl2_RE + ampl2_IM*ampl2_IM)/VOL;
     fprintf(outf,"%.2f % lf % lf % lf % lf % lf % lf\n",t,ampl1_RE,ampl1_IM,ampl2_RE,ampl2_IM,lam1,lam2);
   }
   fclose(outf);
}
// Note: I was confused here about a factor of 2 before cos(evals[p]*t)
// This is simply because I go over the whole set of eigenvalues in
// calculating the overlap.

// Note2: For lambda=0, the spectrum is symmetric about E=0, and hence
// the imaginary part is exactly zero. The imaginary part no longer
// vanishes for lambda!=0, and hence needs to be computed explicitly.
// As in the paper https://arxiv.org/pdf/1709.07461.pdf, eq (10), 
// L(t) = |G(t)|^2, and hence has contributions from both the real
// and the imaginary part. 
