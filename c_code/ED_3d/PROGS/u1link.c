#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#include "define.h"
#include "routines.h"

int main(){
  FILE *fptr;
  char string[50];
  int i,d,p,F;
  fptr = fopen("QUEUE","r");
  if(fptr == NULL)
    {
      printf("could not open QUEUE FILE to open\n");
      exit(1);
    }
  fscanf(fptr,"%s %d\n",string,&LX);
  fscanf(fptr,"%s %d\n",string,&LY);
  fscanf(fptr,"%s %d\n",string,&LZ);
  fscanf(fptr,"%s %d\n",string,&F);
  fscanf(fptr,"%s %lf\n",string,&lam);
  fclose(fptr);
  VOL = LX*LY*LZ;
  Wx=0; Wy=0; Wz=0;

  if(F==1) printf("Going to diagonalize the fermionic QLM. \n");
  else     printf("Going to diagonalize the bosonic QLM. \n");
  /* Allocate memory for configurations */
  for(d=0;d<DIM;d++)
    conf[d] = (int *)malloc(VOL*sizeof(int));

  /* Initialize nearest neighbours */
  for(i=0;i<=2*DIM;i++)
    next[i] = (int *)malloc(VOL*sizeof(int));

  initneighbor();
  /* construct gauge inequivalent states with given (Wx,Wy,Wz) */
  conststates();
  /* construct the hamiltonian for the non-triv states */
  H=allocate2d(N1,N1);
  if(F==1)  constFermH(H,N1);
  else      constH(H,N1);
  deallocate2d(H,N1,N1);

  /* Clear memory */
  for(i=0;i<=2*DIM;i++)
    free(next[i]);

  /* Free memory for configs */
  for(p=0;p<DIM;p++)
   free(conf[p]);



  return 0;
}
