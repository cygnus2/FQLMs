#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_eigen.h>


#include "define.h"
#include "routines.h"

void constFermH(int **H,int size){
   extern int calc_sign(int*, int, int);
   FILE *spin,*flux,*ham;
   int newstate[DIM][VOL];
   int oldstate[DIM*VOL];
   // array oldstate is used to compute the signs. The storage: l0x, l0y, l0z, l1x, l1y, l1z, etc
   // first index is the site and then the three links in three directions
   int q;
   int i,j,k,p,sx,sy,sz;
   int fx,fy,fz;
   int sign;
   int l1,l2,l3,l4;
   int flag;

   int ***stspin;
   int stflux[size][DIM][VOL];
   
   spin = fopen("SPINSTATES","r");
   flux = fopen("FLUXSTATES","r");
 
   // the flux variables are written in the order: xy, xz, yz 
   stspin = allocate3d(size,DIM,VOL);
   for(i=0;i<size;i++){
   for(j=0;j<VOL;j++){
     fscanf(spin,"%d %d %d",&sx,&sy,&sz);
     stspin[i][0][j]=sx; stspin[i][1][j]=sy; stspin[i][2][j]=sz;
     fscanf(flux,"%d %d %d",&fx,&fy,&fz);
     stflux[i][0][j]=fx; stflux[i][1][j]=fy; stflux[i][2][j]=fz;
   }
   }

   fclose(spin);
   fclose(flux);

   /* Act the hamiltonian state by state to construct the full matrix */
   for(i=0;i<size;i++){
    for(j=0;j<VOL;j++){
      /* xy plaquette */
      if(stflux[i][0][j]!=0){ 
        /* copy the state */
        q=0;
        for(k=0;k<VOL;k++){
          newstate[0][k]=stspin[i][0][k]; newstate[1][k]=stspin[i][1][k]; newstate[2][k]=stspin[i][2][k];
          oldstate[q]   =stspin[i][0][k]; oldstate[q+1] =stspin[i][1][k]; oldstate[q+2] =stspin[i][2][k];
          q=q+3;
        }
        if(q != 3*VOL){ printf("Error"); exit(0); } 
        /* act the hamiltonian on the xy plaquette at site j */
        newstate[0][j]=-newstate[0][j];
        newstate[1][j]=-newstate[1][j];
        newstate[1][next[DIM+1][j]]=-newstate[1][next[DIM+1][j]];
        newstate[0][next[DIM+2][j]]=-newstate[0][next[DIM+2][j]];
        /* compare which state it is */
        p=scan(newstate,stspin,size);
        /* check if the flip really exists at the specified location */
        flag=0;
        l1=newstate[0][j]; l2=newstate[1][next[DIM+1][j]]; l3=newstate[0][next[DIM+2][j]]; l4=newstate[1][j];
        if((l1==l2)&&(l3==l4)&&(l1!=l4)) flag=1;
        if(flag==0){ printf("Why no xy flip exists? Check! State=%d, site=%d\n",i,j); exit(0); }
        /* find the sign. j is the lattice site, 0=xy plaquette */
        sign=calc_sign(oldstate, j, 0);  
        if(p>=0){
          if(H[i][p]!=0){ printf("xy plaq! i=%d, p=%d, H(i,p)=%d \n",i,p,H[i][p]); exit(0);  }
          else            H[i][p]=-sign;  // J = -1
        }
      }
      /* xz plaquette */
      if(stflux[i][1][j]!=0){  
        /* copy the state */
        q=0;
        for(k=0;k<VOL;k++){
          newstate[0][k]=stspin[i][0][k]; newstate[1][k]=stspin[i][1][k]; newstate[2][k]=stspin[i][2][k];
          oldstate[q]   =stspin[i][0][k]; oldstate[q+1] =stspin[i][1][k]; oldstate[q+2] =stspin[i][2][k];
          q=q+3;
        }
        if(q!= 3*VOL){ printf("Error"); exit(0); } 
        /* act the hamiltonian on the xz plaquette at site j */
        newstate[0][j]=-newstate[0][j];
        newstate[2][j]=-newstate[2][j];
        newstate[2][next[DIM+1][j]]=-newstate[2][next[DIM+1][j]];
        newstate[0][next[DIM+3][j]]=-newstate[0][next[DIM+3][j]];
        /* compare which state it is */
        p=scan(newstate,stspin,size);
        /* check if the flip really exists at the specified location */
        flag=0;
        l1=newstate[0][j]; l2=newstate[2][next[DIM+1][j]]; l3=newstate[0][next[DIM+3][j]]; l4=newstate[2][j];
        if((l1==l2)&&(l3==l4)&&(l1!=l4)) flag=1;
        if(flag==0){ printf("Why no xz flip exists? Check! State=%d, site=%d\n",i,j); exit(0); }
        /* find the sign. j is the lattice site, 1=xz plaquette */
        sign=calc_sign(oldstate, j, 1);  
        if(p>=0){
          if(H[i][p]!=0){ printf("xz plaq! i=%d, p=%d, H(i,p)=%d \n",i,p,H[i][p]); exit(0); }
          else            H[i][p]=-sign;  // J = -1
        }
      }
      /* yz plaquette */
      if(stflux[i][2][j]!=0){  
        /* copy the state */
        q=0;
        for(k=0;k<VOL;k++){
          newstate[0][k]=stspin[i][0][k]; newstate[1][k]=stspin[i][1][k]; newstate[2][k]=stspin[i][2][k];
          oldstate[q]   =stspin[i][0][k]; oldstate[q+1] =stspin[i][1][k]; oldstate[q+2] =stspin[i][2][k];
          q=q+3;
        }
        if(q!= 3*VOL){ printf("Error"); exit(0); } 
        /* act the hamiltonian on the yz plaquette at site j */
        newstate[1][j]=-newstate[1][j];
        newstate[2][j]=-newstate[2][j];
        newstate[2][next[DIM+2][j]]=-newstate[2][next[DIM+2][j]];
        newstate[1][next[DIM+3][j]]=-newstate[1][next[DIM+3][j]];
        /* compare which state it is */
        p=scan(newstate,stspin,size);
        /* find the sign. j is the lattice site, 2=yz plaquette */
        sign=calc_sign(oldstate, j, 2);   
        if(p>=0){
          if(H[i][p]!=0){ printf("yz plaq! i=%d, p=%d, H(i,p)=%d \n",i,p,H[i][p]); exit(0); }
                      H[i][p]=-sign;   // J = -1  
        }
      }
   }}
   
   deallocate3d(stspin,size,DIM,VOL);
   printf("Hamiltonian construction done. \n");

   /* Print the Hamiltonian */
   //ham=fopen("HAMILTONIAN.dat","w");
   //for(i=0;i<size;i++){
   //  for(j=0;j<size;j++)
   //   fprintf(ham,"% d ",H[i][j]);
   //  fprintf(ham,"\n");
   //}
   //fclose(ham);

   diagH(H,size);

}


/* The arrangment of the plaquette is
          p3  
       o-------o
       |       |
    p4 |       | p2
       |       |
       o-------o
          p1       
 */
int calc_sign(int oldstate[DIM*VOL], int j, int or){
  int p, sign, count, count1;
  int l1, l2, l3, l4;
  int p1, p2, p3, p4;
  sign=1;
  /* collect the links for the relevant plaquette */
  if(or==0){   // xy plaquette
    p1=3*j; p2=3*next[DIM+1][j]+1; p3=3*next[DIM+2][j]; p4=p1+1; 
    l1=oldstate[p1]; l2=oldstate[p2]; l3=oldstate[p3]; l4=oldstate[p4];
  }
  else if(or==1){ // xz plaquette
    p1=3*j; p2=3*next[DIM+1][j]+2; p3=3*next[DIM+3][j]; p4=p1+2;
    l1=oldstate[p1]; l2=oldstate[p2]; l3=oldstate[p3]; l4=oldstate[p4];
  }
  else if(or==2){ // yz plaquette
    p1=3*j+1; p2=3*next[DIM+2][j]+2; p3=3*next[DIM+3][j]+1; p4=p1+1;
    l1=oldstate[p1]; l2=oldstate[p2]; l3=oldstate[p3]; l4=oldstate[p4];
  }
  else{
    printf("wroing orientation! exit(0) \n");
  }
  if((l1==1)&&(l2==1)&&(l3==-1)&&(l4==-1)){ // interaction: cdag(p4) cdag(p3) c(p2) c(p1)
      // action of c(p1)
      count=0;
      for(p=p1+1;p<(DIM*VOL);p++){ if(oldstate[p]==1) count++; }
      if(count%2 ==1) sign = -sign;    // flip sign
      oldstate[p1] = -oldstate[p1];    // flip link; otherwise next counts can go wrong
      // action of c(p2)
      count=0;
      for(p=p2+1;p<(DIM*VOL);p++){ if(oldstate[p]==1) count++; }
      if(count%2 ==1) sign = -sign;    // flip sign
      oldstate[p2] = -oldstate[p2];    // flip link; otherwise next counts can go wrong
      // action of cdag(p3)
      count=0;
      for(p=p3+1;p<(DIM*VOL);p++){ if(oldstate[p]==1) count++; }
      if(count%2 ==1) sign = -sign;    // flip sign
      count1=0;
      for(p=p3;p<(DIM*VOL);p++){  if(oldstate[p]==1) count1++; }
      if(count != count1){ printf("Error in sign routine! \n");  exit(0); }
      oldstate[p3] = -oldstate[p3];    // flip link; otherwise next counts can go wrong
      // action of cdag(p4)
      count=0;
      for(p=p4+1;p<(DIM*VOL);p++){ if(oldstate[p]==1) count++; }
      if(count%2 ==1) sign = -sign;    // flip sign
      count1=0;
      for(p=p4;p<(DIM*VOL);p++){  if(oldstate[p]==1) count1++; }
      if(count != count1){ printf("Error in sign routine p1! \n");  exit(0); }
      oldstate[p4] = -oldstate[p4];    // flip link; otherwise next counts can go wrong
  }
  else if((l1==-1)&&(l2==-1)&&(l3==1)&&(l4==1)){ // interaction: cdg(p1) cdg(p2) c(p3) c(p4)
      // action of c(p4)
      count=0;
      for(p=p4+1;p<(DIM*VOL);p++){ if(oldstate[p]==1) count++; }
      if(count%2 ==1) sign = -sign;     // flip sign
      oldstate[p4] = -oldstate[p4];     // flip link; otherwise next counts will go wrong
      // action of c(p3)
      count=0;
      for(p=p3+1;p<(DIM*VOL);p++){ if(oldstate[p]==1) count++; }
      if(count%2 ==1) sign = -sign;
      oldstate[p3] = -oldstate[p3];
      // action of cdg(p2)
      count=0;
      for(p=p2+1;p<(DIM*VOL);p++){ if(oldstate[p]==1) count++; }
      if(count%2 ==1) sign = -sign;
      count1=0;
      for(p=p2;p<(DIM*VOL);p++){  if(oldstate[p]==1) count1++; }
      if(count != count1){ printf("Error in sign routine p2! \n");  exit(0); }
      oldstate[p2] = -oldstate[p2];
      // action of cdg(p1)
      count=0;
      for(p=p1+1;p<(DIM*VOL);p++){ if(oldstate[p]==1) count++; }
      if(count%2 ==1) sign = -sign;
      count1=0;
      for(p=p1;p<(DIM*VOL);p++){  if(oldstate[p]==1) count1++; }
      if(count != count1){ printf("Error in sign routine! \n");  exit(0); }
      oldstate[p1] = -oldstate[p1];
    }
  else{
       printf("Error in some interaction! \n"); exit(0);
  }
  return sign;
}
