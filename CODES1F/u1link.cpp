#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "define.h"
#include<iostream>
#include<algorithm>
#include<vector>
#include<iterator>

/* according to the rules of cpp, the variables are declared here 
 * and also in the header file as extern such that they are avl to
 * the other functions.
 */
int *next[2*DIM+1];
int *nextCHK[2*DIM+1];
int *chk2lin,*lin2chk;
int LX,LY,VOL,VOL2;
double lam,Ti,Tf,dT;
int NTOT,NH;
std::vector<std::vector<bool>> basis;
std::vector<std::vector<bool>> basis_nonflip;
std::vector<std::vector<bool>> basis_flip;
std::vector<int> Oflip;

int main(){
  FILE *fptr;
  char string[50];
  int i,d,p,q;
  int x,y;
  extern void initneighbor(void);
  extern void conststates(void);
  extern void printbasis(void);
  // eigenvalues and eigenvectors
  std::vector<double> evals;
  std::vector<std::vector<double>> evecs;

  fptr = fopen("QUEUE","r");
  if(fptr == NULL)
    {
      printf("could not open QUEUE FILE to open\n");
      exit(1);
    }
  fscanf(fptr,"%s %d\n",string,&LX);
  fscanf(fptr,"%s %d\n",string,&LY);
  fscanf(fptr,"%s %lf\n",string,&lam);
  fscanf(fptr,"%s %lf\n",string,&Ti);
  fscanf(fptr,"%s %lf\n",string,&Tf);
  fscanf(fptr,"%s %lf\n",string,&dT);
  fclose(fptr);
  VOL = LX*LY;
  VOL2 = VOL/2;

  /* Initialize nearest neighbours */
  for(i=0;i<=2*DIM;i++){
    next[i] = (int *)malloc(VOL*sizeof(int)); 
    nextCHK[i] = (int *)malloc(VOL*sizeof(int));
  }

  /* Initialize checkerboard co-ordinates */
  lin2chk = (int *)malloc(VOL*sizeof(int));
  chk2lin = (int *)malloc(VOL*sizeof(int));
  initneighbor();

  /* construct states */
  conststates();
  // for(i=0;i<NTOT;i++) printf("%d\n",(int)basis[i].size());
  // for(i=0;i<NTOT;i++) std::cout<<basis[i].size()<<std::endl;
  //printbasis();
  /* construct and diagonalize the matrix */
  constH(evals, evecs);
  //if(basis[3][4]==true) printf("spin = 1\n");
  //else printf("spin = -1\n");
  /* real-time evolution of cartoon states and Locshmidt Echo */
  evolve_cartoons(evals, evecs);
  //std::cout<<basis[2][3]<<std::endl;
  //std::cout<<basis[3][4]<<std::endl;
  /* calculate the expectation value of Oflip for every eigenstate */
  calc_Oflip(evals,evecs);

  /* Clear memory */
  for(i=0;i<=2*DIM;i++){  free(next[i]); free(nextCHK[i]); }
  free(chk2lin); free(lin2chk);
  return 0;
}
