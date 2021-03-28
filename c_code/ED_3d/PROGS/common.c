#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#include "define.h"
#include "routines.h"

int **allocate2d(int row, int col)
{
  int i,j;
  int **mat;
  mat = (int **)malloc(row*sizeof(int*));
  if(mat==NULL) {printf("Out of memory\n"); exit(0);}

  for(i=0;i<row;i++){
   mat[i]=(int *)malloc(col*sizeof(int));
   if(mat[i]==NULL)  {printf("Out of memory\n"); exit(0);}
  }

 for(i=0;i<row;i++)
 for(j=0;j<col;j++)
  mat[i][j]=0;

 return mat;

}

double **allocatedouble2d(int row, int col)
{
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

int ***allocate3d(int size1, int size2, int size3)
{
   int i,j,k;
   int ***mat;
   mat = (int ***)malloc(size1*sizeof(int**));

   for(i=0;i<size1;i++)
     mat[i]= (int **)malloc(size2*sizeof(int*));

   for(i=0;i<size1;i++)
   for(j=0;j<size2;j++)
     mat[i][j] = (int *)malloc(size3*sizeof(int));

   for(i=0;i<size1;i++)
   for(j=0;j<size2;j++)
   for(k=0;k<size3;k++)
     mat[i][j][k]=0;

   return mat;
}

void deallocate2d(int **mat, int row, int col)
{
  int i;
  for(i=0;i<row;i++)
   free(mat[i]);

 free(mat);
}

void deallocatedouble2d(double **mat, int row, int col)
{
  int i;
  for(i=0;i<row;i++)
   free(mat[i]);

 free(mat);
}

void deallocate3d(int ***mat, int size1, int size2, int size3)
{
  int i,j;
  for(i=0;i<size1;i++)
  for(j=0;j<size2;j++)
    free(mat[i][j]);

  for(i=0;i<size1;i++)
    free(mat[i]);

  free(mat);
}


void initneighbor(void)
{
  int p,x,y,z;
  for(p=0;p<VOL;p++)
    {
      x = (p%LX);
      y = (p/LX)%LY;
      z = (p/(LX*LY))%LZ;

      next[DIM][p] = p;
      next[DIM+1][p] = z*LX*LY + y*LX + ((x+1)%LX);
      next[DIM+2][p] = z*LX*LY + ((y+1)%LY)*LX + x;
      next[DIM+3][p] = ((z+1)%LZ)*LX*LY + y*LX + x;
      next[DIM-1][p] = z*LX*LY + y*LX + ((x-1+LX)%LX);
      next[DIM-2][p] = z*LX*LY + ((y-1+LY)%LY)*LX + x;
      next[DIM-3][p] = ((z-1+LZ)%LZ)*LX*LY + y*LX + x;
    }
}

void initconf(void)
{
  int p,d;

  for(p=0;p<VOL;p++)
      for(d=0;d<DIM;d++) conf[d][p] = 0;
}

void printconf(){
  int p,d;

  for(p=0;p<VOL;p++){
   printf("site %d(x,y):",p+1);
   for(d=0;d<DIM;d++)
    printf("% d ",conf[d][p]);
  }
  printf("\n");
}

void storeconf(FILE *fp){
  int p,d;
  for(p=0;p<VOL;p++){
   for(d=0;d<DIM;d++)
    fprintf(fp,"% d ",conf[d][p]);
  }
  fprintf(fp,"\n");
}

int checkflux(FILE *fptr)
{
  int p,d,flux[DIM][VOL],flag;
  int link1,link2,link3,link4;
  double wx,wy,wz,flagW;

  flag=0; flagW=0;
  /* initialize flux variable */
  for(p=0;p<VOL;p++) for(d=0;d<DIM;d++) flux[d][p]=0;
  /* initialize the winding */
  wx=0.0; wy=0.0; wz=0.0;
  /* look for non-trivial flux */
  // The link variables are labelled as
  //           link3
  //        O--------O
  //        |        |
  //  link4 |        | link2
  //        |        |
  //        O--------O
  //          link1
  //
  for(p=0;p<VOL;p++){
  /* xy plaquette */
  link1 = conf[0][p];              link2 = conf[1][next[DIM+1][p]];
  link3 = conf[0][next[DIM+2][p]]; link4 = conf[1][p];
  if((link1==link2)&&(link3==link4)&&(link1!=link4)){
      flux[0][p]=1; flag=1;
  }

  /* xz plaquette */
  link1 = conf[0][p];              link2 = conf[2][next[DIM+1][p]];
  link3 = conf[0][next[DIM+3][p]]; link4 = conf[2][p];
  if((link1==link2)&&(link3==link4)&&(link1!=link4)){
      flux[1][p]=1; flag=1;
  }

  /* yz plaquette */
  link1 = conf[1][p];              link2 = conf[2][next[DIM+2][p]];
  link3 = conf[1][next[DIM+3][p]]; link4 = conf[2][p];
  if((link1==link2)&&(link3==link4)&&(link1!=link4)){
      flux[2][p]=1; flag=1;
  }

  // compute the windings
  wx = wx + conf[0][p];   wy = wy + conf[1][p];  wz = wz + conf[2][p];
 }

 wx /= LX; wy /= LY; wz/= LZ;

 /* if-option, if activated will store only states matching specified winding */
 // flux states indicate the flippabilites of xy, xz and yz plaquettes
 // at the site p
 if(wx==Wx && wy==Wy && wz==Wz){
  for(p=0;p<VOL;p++){
      for(d=0;d<DIM;d++)
        fprintf(fptr,"% d ",flux[d][p]);
  }
  fprintf(fptr,"\n");
  flagW=1;
 }

 return flagW;
}

int checkconf(void)
{
  int p,d,k;
  int flag=1;

  /* check Gauss Law at each site */

  for(p=0;p<VOL;p++)
    {
      k=0;
      for(d=0;d<DIM;d++)
        {
          k += conf[d][p];
          k -= conf[d][next[DIM-(d+1)][p]];
        }
      if(k!=0) flag=0;
    }
   return flag;
}


int scan(int newstate[][VOL],int ***stspin,int num)
{
   int p,flag,i,j;
   for(i=0;i<num;i++){
     p=-1;
     flag=0;
     for(j=0;j<VOL;j++){
      if(newstate[0][j]==stspin[i][0][j]) flag++;
      if(newstate[1][j]==stspin[i][1][j]) flag++;
      if(newstate[2][j]==stspin[i][2][j]) flag++;
     }
     if(flag==3*VOL) { p=i; break; }
   }

   return p;
}
