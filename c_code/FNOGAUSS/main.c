#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "constants.h"
#include "diagonalize.h"
#include "categorize.h"
#include "observable.h"

size_t n = 8;


void matrixInit( double *h, size_t *pointers, size_t *dimensions, size_t *states, int dim );
int isEven( int n);


int main(void)
{
  size_t dim;
  dim = (size_t) 1 << n;
  size_t length = 0;
  size_t k, l;
  size_t i, j;
  double *h;
  size_t *states;
  size_t *pointers;
  size_t *dimensions;
  double ans;

  pointers = (size_t *)malloc( (n + 1 ) * sizeof( size_t ) );
  dimensions = (size_t *)malloc( (n + 1) * sizeof( size_t) );
  states = (size_t *)malloc( dim * sizeof(size_t));

  sort( states, n );
  dims( dimensions, n );
  locations( pointers, dimensions, n );

  for( i = 0; i < n + 1; i++ )
    {
      printf("dimensions[ %d ] = %zu \n", i, dimensions[ i ]);
    }
  
  for( i = 0; i < n + 1; i++ )
    {
      length += dimensions[i] * dimensions[i];
    }

  //printf("the dimension is %zu", length);
  
  h = malloc( length * sizeof(double) );

  /*for( i = 0; i < dim; i++)
    {
      printf("states[%d] = %d\n", i, states[i]);
    }


  for( i = 0; i < dim; i++ )
    {
      printf("state %d has %d particles \n", i, bin(i , n));
    }
  printf("\n");
  for( i = 0; i < dim; i++ )
    {
      printf("the %d element of states is %zu\n", i, states[i]);
    }*/

  matrixInit(h, pointers, dimensions, states, dim);
  printf("matrix init completed.\n");
  /*for( l = 0; l < n + 1; l++)
    {
      for(i = 0; i < dimensions[ l ]; i++)
	{
	  for(j = 0; j < dimensions[ l ]; j++ )
	    {
	      printf("%10.3f ", h[pointers[l] + i + dimensions[l] * j]);
	    }
	  printf("\n");
	}
      printf("\n");
      }*/

  expBlocks(h, pointers, dimensions, n);

  //printf("\n");

  /*for( l = 0; l < n + 1; l++)
    {
      for(i = 0; i < dimensions[ l ]; i++)
	{
	  for(j = 0; j < dimensions[ l ]; j++ )
	    {
	      printf("%10.3f ", h[pointers[l] + i + dimensions[l] * j]);
	    }
	  printf("\n");
	}
      printf("\n");
      }*/

  //ans = observable(h, pointers, dimensions, states, n, dim);

  //printf("\nthe average observable is %.12f\n", ans);

  free( h );
  free( states );
  free( dimensions );
  free( pointers );
  
  return 0;
}


void matrixInit( double *h, size_t *pointers, size_t *dimensions, size_t *states, int dim)
{
  int parity( int k, int i );
  int annihilate( int k, int i);
  int create( int i, int j );
  double couplingV(int i);
  int nbr( int i, int q);
  int i, j, k, l;
  int count = 0;
  int counti = 0;
  int countj = 0;

  /*This part is unfinished!!*/
  for( l = 0; l < n + 1; l++ )
    {
      printf("l = %d \n", l);
      for( i = 0; i < dimensions[ l ]; i++ )
	{
	  for( j = 0; j < dimensions[ l ]; j++ )
	    {
	      for( k = 0; k < n; k++ )
		{
		  if(create( 0, create(3, annihilate( 4, annihilate(1, states[i + count] ) ) ) ) == states[j + count] )
		    {
		      //h[pointers[ l ] + i + dimensions[ l ] * j] += (double) T * parity(1, states[i + count]) * parity( 4, annihilate(1, states[i + count]) ) * parity( 3, annihilate( 4, annihilate(1, states[i + count]) ) ) * parity( 0, create( 3, annihilate( 4, annihilate(1, states[i + count]) ) ) );
		      h[pointers[ l ] + i + dimensions[ l ] * j] += 1.0;
		    }
		  if(create( 2, create(1, annihilate( 6, annihilate(3, states[i + count] ) ) ) ) == states[j + count] )
		    {
		      //h[pointers[ l ] + i + dimensions[ l ] * j] += (double) T * parity(3, states[i + count]) * parity( 6, annihilate(3, states[i + count]) ) * parity( 1, annihilate( 6, annihilate(3, states[i + count]) ) ) * parity( 2, create( 1, annihilate( 6, annihilate(3, states[i + count]) ) ) );
		      h[pointers[ l ] + i + dimensions[ l ] * j] += 1.0;
		    }
		  if(create( 4, create(7, annihilate( 0, annihilate(5, states[i + count] ) ) ) ) == states[j + count] )
		    {
		      //h[pointers[ l ] + i + dimensions[ l ] * j] += (double) T * parity(5, states[i + count]) * parity( 0, annihilate(5, states[i + count]) ) * parity( 7, annihilate( 0, annihilate(5, states[i + count]) ) ) * parity( 4, create( 7, annihilate( 0, annihilate(5, states[i + count]) ) ) );
		      h[pointers[ l ] + i + dimensions[ l ] * j] += 1.0;
		    }
		  if(create( 6, create(5, annihilate( 2, annihilate(7, states[i + count] ) ) ) ) == states[j + count] )
		    {
		      //h[pointers[ l ] + i + dimensions[ l ] * j] += (double) T * parity(7, states[i + count]) * parity( 2, annihilate(7, states[i + count]) ) * parity( 5, annihilate( 2, annihilate(7, states[i + count]) ) ) * parity( 6, create( 5, annihilate( 2, annihilate(7, states[i + count]) ) ) );
		      h[pointers[ l ] + i + dimensions[ l ] * j] += 1.0;
		    }
		  if(create( 1, create(4, annihilate( 3, annihilate(0, states[i + count] ) ) ) ) == states[j + count] )
		    {
		      //h[pointers[ l ] + i + dimensions[ l ] * j] += (double) T * parity(0, states[i + count]) * parity( 3, annihilate(0, states[i + count]) ) * parity( 4, annihilate( 3, annihilate(0, states[i + count]) ) ) * parity( 1, create( 4, annihilate( 3, annihilate(0, states[i + count]) ) ) );
		      h[pointers[ l ] + i + dimensions[ l ] * j] += 1.0;
		    }
		  if(create( 3, create(6, annihilate( 1, annihilate(2, states[i + count] ) ) ) ) == states[j + count] )
		    {
		      //h[pointers[ l ] + i + dimensions[ l ] * j] += (double) T * parity(2, states[i + count]) * parity( 1, annihilate(2, states[i + count]) ) * parity( 6, annihilate( 1, annihilate(2, states[i + count]) ) ) * parity( 3, create( 6, annihilate( 1, annihilate(2, states[i + count]) ) ) );
		      h[pointers[ l ] + i + dimensions[ l ] * j] += 1.0;
		    }
		  if(create( 5, create(0, annihilate( 7, annihilate(4, states[i + count] ) ) ) ) == states[j + count] )
		    {
		      //h[pointers[ l ] + i + dimensions[ l ] * j] += (double) T * parity(4, states[i + count]) * parity( 7, annihilate(4, states[i + count]) ) * parity( 0, annihilate( 7, annihilate(4, states[i + count]) ) ) * parity( 5, create( 0, annihilate( 7, annihilate(4, states[i + count]) ) ) );
		      h[pointers[ l ] + i + dimensions[ l ] * j] += 1.0;
		    }
		  if(create( 7, create(2, annihilate( 5, annihilate(6, states[i + count] ) ) ) ) == states[j + count] )
		    {
		      //h[pointers[ l ] + i + dimensions[ l ] * j] += (double) T * parity(6, states[i + count]) * parity( 5, annihilate(6, states[i + count]) ) * parity( 2, annihilate( 5, annihilate(6, states[i + count]) ) ) * parity( 7, create( 2, annihilate( 5, annihilate(6, states[i + count]) ) ) );
		      h[pointers[ l ] + i + dimensions[ l ] * j] += 1.0;
		    }
		}
	    }
	}
      count += dimensions[ l ];
    }

  count = 0;
  
  /*for( l = 0; l < n + 1; l++)
    {
      for(i = 0; i < dimensions[ l ]; i++)
	{
	  h[pointers[l] + i + dimensions[ l ] * i] = couplingV( states[ i + count ] );
	}
      count += dimensions[ l ];
      }*/
}

/*k is a site, i is a state*/
int parity( int k, int i)
{
  int step( int n, int k );
  int j;
  int par = 0;
  for( j = 0; j < n-1; j++ )
    {
      par += (i / ( (int) pow(2.0, (double)j)) ) % 2 * step(j, k);
    }
  return 1 - 2 * (par % 2);
}

/*"q" specifies which neighbor we want for coordinate i*/
int nbr( int i, int q)
{
  int i_x = i % LX;
  int i_y = (i / LX) % LY;
  if( q == 0 )
    {
      return (i_x + 1) % LX + LX * i_y;
    }
  else if( q == 1)
    {
      return (i_x - 1 + LX) % LX + LX * i_y;
    }
  else if( q == 2)
    {
      return i_x + LX * ((i_y + 1) % LY );
    }
  else
    {
      return i_x + LX * ((i_y - 1 + LY) % LY);
    }
}

/*Fermion sign functions*/
int step( int n, int k )
{
  if(n < k)
    {
      return 1;
    }
  else
    {
      return 0;
    }
}

int annihilate( int k, int i )
{
  if( i == -1 )
    {
      return -1;
    }
  else if( ( ( i / (int) pow( 2.0, (double) k ) ) % 2 ) == 0 )
    {
      return -1;
    }
  else
    {
      return i - (int) pow( 2.0, (double) k );
    }
}

int create( int k, int i )
{
  if( i == -1 )
    {
      return -1;
    }
  else if( ( ( i / (int) pow(2.0, (double) k) ) % 2 ) == 1 )
    {
      return -1;
    }
  else
    return i + (int) pow( 2.0, (double) k );
}

double couplingV(int i)
{
  int j;
  double sum = 0.0;
  for( j = 0; j < n; j++ )
    {
      sum += V * (((i / (int) pow(2.0, (double) nbr(j,0))) % 2 - 1.0/2.0) *
		  ((i / (int) pow(2.0, (double) j)) % 2 - 1.0/2.0)  );
      sum += V * (((i / (int) pow(2.0, (double) nbr(j,2))) % 2 - 1.0/2.0) *
		  ((i / (int) pow(2.0, (double) j)) % 2 - 1.0/2.0)  );
    }
  if( LX == 2 && LY == 2 )
    sum /= 2.0;
  return sum;
}

int isEven(int n)
{
  int n_x;
  int n_y;
  n_x = n % LX;
  n_y = (n / LX) % LY;
  if( (n_x % 2) == (n_y % 2) )
    return 0;
  else
    return 1;
}
