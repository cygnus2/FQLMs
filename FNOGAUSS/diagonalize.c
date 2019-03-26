#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "diagonalize.h"

void exponentiate( double * m, double * evec, double * eval, size_t length )
{
  void diagonalize( double * m, double * evec, double * eval, size_t length );
  
  int i,j;
  double * h;
  double alpha = 1.0;
  double beta = 0.0;
  char transA;
  char transB;

  h = (double *)malloc(  length * length * sizeof(double) ); //placeholder

  diagonalize(m, evec, eval, length);

  printf("\n Energies are: \n");

  for(i = 0; i < length; i++)
    {
      printf("%.8f ", eval[i]);
    }
  printf("\n\n");
  
  /*for( i = 0; i < length; i++ )
    {
      eval[i] = exp( - BETA * eval[ i ] );
    }

  for( i = 0; i < length; i++ )
    {
      for( j = 0; j < length; j++ )
	{
	  if( i != j )
	    m[ i + length * j ] = 0.0;
	  else
	    m[ i + length * j ] = eval[ i ];
	  //printf("%10.3f", m[i + length * j]);
	}
      //printf("\n");
    }
    printf("\n");

  transA = 'N';
  transB = 'N';

  dgemm_( &transA, &transB, &length, &length, &length, &alpha, evec, &length, m, &length, &beta, h, &length);

  transB = 'T';
  dgemm_( &transA, &transB, &length, &length, &length, &alpha, h, &length, evec, &length, &beta, m, &length);*/
  
  free( h );
}

void diagonalize( double * m, double * evec, double * eval, size_t length )
{
  char comp;
  int i;
  int lWork, info, lg, liWork;
  double *e, *tau, *work;
  int * iWork;

  for( i = 0; i < length * length; i++ )
    {
      evec[ i ] = m[ i ];
    }

  comp = 'U';
  lWork = length * length;
  info = 0;
  lg = (int)ceil( log(length) / log(2.0) );
  liWork = 6 + 6 * length + 5 * length * lg;

  e = (double *)malloc( ( length - 1 ) * sizeof(double) );
  tau = (double *)malloc( ( length - 1 ) * sizeof(double) );
  work = (double *)malloc( lWork * sizeof(double) );
  iWork = (int *)malloc( liWork * sizeof(int) );

  /*Finding eigenvectors and eignevalues*/
  printf("length is %zu \n", length);
  dsytrd_( &comp, &length, evec, &length, eval, e, tau, work, &lWork, &info );

  dorgtr_( &comp, &length, evec, &length, tau, work, &lWork, &info );

  comp = 'V';
  lWork = 1 + 3 * length + 2 * length * lg + 4 * length * length;
  work = (double *)realloc( work, lWork * sizeof(double) );
  dstedc_( &comp, &length, eval, e, evec, &length, work, &lWork, iWork, &liWork, &info );

  free(e);
  free(tau);
  free(work);
  free(iWork);
}

void test( double * m, double * s, double * e, int n) //testing to make sure diagonalization works
{
  int i, j;
  double * h;
  double alpha = 1.0;
  double beta = 0.0;
  char transA;
  char transB;

  h = (double *)malloc(  n * n * sizeof(double) ); //placeholder
  
  for( i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      h[ i + n * j] = 0.0;

  transA = 'T';
  transB = 'N';
  
  dgemm_( &transA, &transB, &n, &n, &n, &alpha, s, &n, m, &n, &beta, h, &n); 
  
  transA = 'N';
  dgemm_( &transA, &transB, &n, &n, &n, &alpha, h, &n, s, &n, &beta, m, &n);
  /*for( i = 0; i < n; i++ )
    {
      for( j = 0; j < n; j++ )
	{
	  printf("%.3f ", m[i + n * j]);
	}
      printf("\n");
      }*/

  free( h );
}

void expBlocks(double * m, size_t *pointers, size_t *dimensions, int n)
{
  size_t i;
  double *e;
  double *s;
  for( i = 0; i < n + 1; i++ )
    printf("dimensions[ i ] = %zu \n", dimensions[ i ]);
  for( i = 0; i < n + 1; i++ )
    {
      printf("i is %d \n", i);
      e = (double *)malloc( dimensions[i] * sizeof(double) );
      s = (double *)malloc( dimensions[i] * dimensions[i] * sizeof(double) );
      
      exponentiate( &(m[pointers[i]]), s, e, dimensions[i] );
      
      free(e);
      free(s);
    }
}
