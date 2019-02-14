/* calculate the local operator Oflip */
/* Oflip calculates the expectation value of the number of flippable
   plaquettes in a given eigenstate, normalized with the total number
   of plaquettes   */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

// Notation: eigenstate |n> = \sum_k \alpha_k |k>, |k> is a cartoon state
void calc_Oflip(std::vector<double> &evals, std::vector<std::vector<double>> &evecs){
  int p,q;
  double Oflip_avg;
  FILE *outf;
  outf = fopen("Oflip.dat","w");
  // scan through all the eigenvalues
  for(p=0;p<NH;p++){
    // calculate the expectation value in each eigenstate
    Oflip_avg = 0.0;
    for(q=0;q<NH;q++){
      Oflip_avg += evecs[p][q]*evecs[p][q]*Oflip[q];
    }
    Oflip_avg = Oflip_avg/((double)VOL);
    fprintf(outf,"%lf %lf\n",evals[p],Oflip_avg);
  }
 fclose(outf);
}
