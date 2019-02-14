#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "define.h"

void initneighbor(void)
{
  int p,x,y,q;
  int lin;

  /* lattice linear co-ordinates */
  for(p=0;p<VOL;p++){
      x = (p%LX);
      y = (p/LX)%LY;

      next[DIM][p] = p;
      next[DIM+1][p] = y*LX + ((x+1)%LX);
      next[DIM+2][p] = ((y+1)%LY)*LX + x;
      next[DIM-1][p] = y*LX + ((x-1+LX)%LX);
      next[DIM-2][p] = ((y-1+LY)%LY)*LX + x;
  }
 
  /* checkerboard to linear and vice-versa */
  p=0;
  for(y=0;y<LY;y++){
  for(x=0;x<LX;x++){
    if((x+y)%2==0){
       q = y*LX + x;
       chk2lin[p] = q; lin2chk[q] = p;
       printf("checkerboard site = %d(%d,%d); Lin site = %d\n",p,x,y,q);
       p++;
    }
  }}
  for(y=0;y<LY;y++){
  for(x=0;x<LX;x++){
    if((x+y)%2==1){
       q = y*LX + x;
       chk2lin[p]= q; lin2chk[q] = p;
       printf("checkerboard site = %d(%d,%d); Lin site = %d\n",p,x,y,q);
       p++;
    }
  }} 

  // nextCHK targets the neighboring sites in the string 
  // representing the GL implementation on the even sites
  for(p=0;p<VOL;p++){
      lin = chk2lin[p];
      nextCHK[DIM][p]   = p;
      nextCHK[DIM+1][p] = lin2chk[next[DIM+1][lin]];
      nextCHK[DIM-1][p] = lin2chk[next[DIM-1][lin]]; 
      nextCHK[DIM+2][p] = lin2chk[next[DIM+2][lin]];
      nextCHK[DIM-2][p] = lin2chk[next[DIM-2][lin]];
      //printf("checkerboard NN. site = %d (in chkrboard = %d)\n",lin,p);
      //printf("+x=%d,   -x=%d,   +y=%d,   -y=%d\n",nextCHK[DIM+1][p],
      //nextCHK[DIM-1][p],nextCHK[DIM+2][p],nextCHK[DIM-2][p]); 
   }
}

/* configurations are stored such that i=2*p and i=2*p+1
   denotes the x and the y links of the site p=y*LX + x.
 */
void printconf(bool *conf){
  int p,d,x,y;
  for(p=0;p<VOL;p++){
   x=p%LX; y=p/LX;
   if(conf[2*p]==true)         printf("site (x,y)=(%d,%d); link-x=% d \n",x,y,1);
   else if(conf[2*p]==false)   printf("site (x,y)=(%d,%d); link-x=% d \n",x,y,-1);
   if(conf[2*p+1]==true)       printf("site (x,y)=(%d,%d); link-y=% d \n",x,y,1);
   else if(conf[2*p+1]==false) printf("site (x,y)=(%d,%d); link-y=% d \n",x,y,-1);
  }
}


/* Prints all the basis; use only for checking */
void printbasis(){
 FILE *fptr;
 int i,p;
 fptr=fopen("BASIS_FLIP.txt","w");
 fprintf(fptr,"Printing the (flippable) basis states for (Lx,Ly)=(%d,%d) lattice.\n",LX,LY);
 for(i=0;i<NH;i++){
   fprintf(fptr,"basis %d ",i);
   for(p=0;p<2*VOL;p++){
    fprintf(fptr,"%d ",(int)basis_flip[i][p]);
   }
   fprintf(fptr,"\n");
  }
 fclose(fptr);
}

double **allocatedouble2d(int row, int col){
  int i,j;
  double **mat;
  mat = (double **)malloc(row*sizeof(double*));
  if(mat==NULL) {printf("Out of memory\n"); exit(0);}

  for(i=0;i<row;i++){
   mat[i]=(double *)malloc(col*sizeof(double));
   if(mat[i]==NULL)  {printf("Out of memory\n"); exit(0);}
  }

 for(i=0;i<row;i++)
 for(j=0;j<col;j++)
  mat[i][j]=0.0;

 return mat;
}

void deallocate2d(int **mat, int row, int col){
  int i;
  for(i=0;i<row;i++)
   free(mat[i]);

 free(mat);
}


