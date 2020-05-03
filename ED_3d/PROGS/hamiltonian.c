#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_eigen.h>


#include "define.h"
#include "routines.h"

void constH(int **H,int size){
   FILE *spin,*flux,*ham;
   int newstate[DIM][VOL];
   int i,j,k,p,sx,sy,sz;
   int fx,fy,fz;
   int size1;
   size1=N-size;

   int ***stspin;
   int stflux[size][DIM][VOL];

   printf("Hamiltonian for bosonic links. \n");
   
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
        for(k=0;k<VOL;k++){
          newstate[0][k]=stspin[i][0][k]; newstate[1][k]=stspin[i][1][k]; newstate[2][k]=stspin[i][2][k];
        }
        /* act the hamiltonian on the xy plaquette at site j */
        newstate[0][j]=-newstate[0][j];
        newstate[1][j]=-newstate[1][j];
        newstate[1][next[DIM+1][j]]=-newstate[1][next[DIM+1][j]];
        newstate[0][next[DIM+2][j]]=-newstate[0][next[DIM+2][j]];
        /* compare which state it is */
        p=scan(newstate,stspin,size);
        if(p>=0){
          if(H[i][p]!=0){ printf("xy plaq! i=%d, p=%d, H(i,p)=%d \n",i,p,H[i][p]); exit(0);  }
          else            H[i][p]=-1;
        }
      }
      /* xz plaquette */
      if(stflux[i][1][j]!=0){  
        /* copy the state */
        for(k=0;k<VOL;k++){
          newstate[0][k]=stspin[i][0][k]; newstate[1][k]=stspin[i][1][k]; newstate[2][k]=stspin[i][2][k];
        }
        /* act the hamiltonian on the xz plaquette at site j */
        newstate[0][j]=-newstate[0][j];
        newstate[2][j]=-newstate[2][j];
        newstate[2][next[DIM+1][j]]=-newstate[2][next[DIM+1][j]];
        newstate[0][next[DIM+3][j]]=-newstate[0][next[DIM+3][j]];
        /* compare which state it is */
        p=scan(newstate,stspin,size);
        if(p>=0){
          if(H[i][p]!=0){ printf("xz plaq! i=%d, p=%d, H(i,p)=%d \n",i,p,H[i][p]); exit(0); }
          else            H[i][p]=-1;
        }
      }
      /* yz plaquette */
      if(stflux[i][2][j]!=0){  
        /* copy the state */
        for(k=0;k<VOL;k++){
          newstate[0][k]=stspin[i][0][k]; newstate[1][k]=stspin[i][1][k]; newstate[2][k]=stspin[i][2][k];
        }
        /* act the hamiltonian on the yz plaquette at site j */
        newstate[1][j]=-newstate[1][j];
        newstate[2][j]=-newstate[2][j];
        newstate[2][next[DIM+2][j]]=-newstate[2][next[DIM+2][j]];
        newstate[1][next[DIM+3][j]]=-newstate[1][next[DIM+3][j]];
        /* compare which state it is */
        p=scan(newstate,stspin,size);
        if(p>=0){
          if(H[i][p]!=0){ printf("yz plaq! i=%d, p=%d, H(i,p)=%d \n",i,p,H[i][p]); exit(0); }
                      H[i][p]=-1;
        }
      }
   }}
   
   deallocate3d(stspin,size,DIM,VOL);

   printf("Hamiltonian for bosonic links done. \n");

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

void diagH(int **H,int size){

  FILE *vals,*vecs; 
  double *mat;
  //double norm;
    
  int i,j,k;
  int Nsq=size*size;
  mat = (double *) malloc(Nsq*sizeof(double));
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
 	mat[i+size*j]=H[i][j];
  
  gsl_matrix_view m=gsl_matrix_view_array(mat,size,size);

  gsl_vector *eval=gsl_vector_alloc(size);
  gsl_matrix *evec=gsl_matrix_alloc(size,size);
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(size);
  gsl_eigen_symmv(&m.matrix,eval,evec,w);
  gsl_eigen_symmv_free(w);
  gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_ASC);

  vals=fopen("evals.dat","w");
  vecs=fopen("evecs.dat","w");

  for(i=0;i<size;i++){
   double eval_i=gsl_vector_get(eval,i);
   //gsl_vector_view evec_i=gsl_matrix_column(evec,i);

   //printf("eigenvalue=% g\n",eval_i);
   //printf("eigenvector=\n");
   //gsl_vector_fprintf(stdout, &evec_i.vector, "%g");

        fprintf(vals,"%d % .18le\n",i,eval_i);
        for(k=0;k<size;k++)
         fprintf(vecs,"% .6e ",gsl_matrix_get(evec,k,i));
        fprintf(vecs,"\n");
      }
    
   fclose(vals);
   fclose(vecs);  

   /* Check for normalization */
   /*
   for(i=0;i<size;i++){
    norm=0;
    for(j=0;j<size;j++) 
      norm += pow(gsl_matrix_get(evec,j,i),2);
    printf("evec %d, norm=%f\n",i,norm);
   }
   */
   /* Check for orthogonality */
   /*
   for(i=0;i<size;i++){
    for(j=i+1;j<size;j++){
      norm=0;
      for(k=0;k<size;k++)
        norm += gsl_matrix_get(evec,k,i)*gsl_matrix_get(evec,k,j);
      printf("evecs %d %d, norm=%f\n",i,j,norm);
     }
    }
   */

  gsl_vector_free(eval);
  gsl_matrix_free(evec);
}


