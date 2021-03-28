  /* construct gauge inequivalent states of the U(1) model */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#include "define.h"
#include "routines.h"

void conststates(){
  FILE *fptr,*fptr1;
  int d,p,dum;
  int count,fluxcount;
  int flagGI, flagFLUX;
  long long int NST,n,k;
  int wx, wy, wz;

  fptr =fopen("FLUXSTATES","w");
  fptr1=fopen("SPINSTATES","w");

  /* count states */
  NST = pow(2,DIM*VOL);
  count = 0;
  fluxcount = 0;
  printf("counting: %lld\n",NST);
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

   flagFLUX = checkflux(fptr);
   /* Store states with matching windings */
   if(flagFLUX) storeconf(fptr1);
   //else storeconf(tfptr1);
   fluxcount += flagFLUX;
   }
  }

  fclose(fptr);
  fclose(fptr1);
  //fclose(tfptr);
  //fclose(tfptr1);
  printf("Total GI states: %d\n",count);
  printf("States in sector (%d,%d,%d) is: %d\n",(int)Wx,(int)Wy,(int)Wz,fluxcount);

  /* Initialize */
  N=count;
  N1=fluxcount;

}
