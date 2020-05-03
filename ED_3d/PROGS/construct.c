  /* construct gauge inequivalent states of the U(1) model */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#include "define.h"
#include "routines.h"

void conststates(){
  FILE *fptr,*fptr1;
  FILE *tfptr,*tfptr1;
  int d,p,dum;
  int count,fluxcount;
  int flagGI, flagFLUX;
  long long int NST,n,k;

  fptr =fopen("FLUXSTATES","w");
  fptr1=fopen("SPINSTATES","w");
  //tfptr =fopen("TFLUXSTATES","w");
  //tfptr1=fopen("TSPINSTATES","w");

  /* count states */
  NST = pow(2,DIM*VOL);
  count = 0;
  fluxcount = 0;
  printf("No of states: %lld\n",NST);
  for(n=0;n<NST;n++){
   k=n;
   for(p=0;p<VOL;p++){
    for(d=0;d<DIM;d++){
     dum = k%2;
     if(dum==0) conf[d][p]= 1;
     else conf[d][p] = -1;
      k = k/2;
     }
   }
   flagGI = checkconf();
   count += flagGI;
   if(flagGI){
   //storeconf(fptr1);
   flagFLUX = checkflux(fptr,tfptr);
   /* This option will only store states with non-trivial flux */
   if(flagFLUX) storeconf(fptr1);
   //else storeconf(tfptr1);
   fluxcount += flagFLUX;
   }
  }

  fclose(fptr);
  fclose(fptr1);
  //fclose(tfptr);
  //fclose(tfptr1);
  printf("No of gauge invariant states: %d\n",count);
  printf("No of gauge invariant states with non-trivial flux: %d\n",fluxcount);

  /* Initialize */
  N=count;
  N1=fluxcount;

}
