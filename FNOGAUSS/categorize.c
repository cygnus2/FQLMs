#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "constants.h"
#include "categorize.h"

/*A state is input, and its number of particles is output*/
int bin( size_t k, size_t n )
{
  int sum = 0;
  size_t i;
  size_t dim = (size_t) 1 << n;
  for( i = 0; i < n; i++ )
    {
      sum += (k / (1 << i)) % 2;
    }

  return sum;
}

/*Sorts the states by particle number*/
void sort( size_t *states, size_t n )
{
  size_t dim;
  dim = (size_t) 1 << n;
  int i;
  size_t *sector;
  /*Sector is an array of n-1 index values. Each index gives the
    location in the states array for the next sector of states
    there's a zero particle sector, one particle sector, etc.*/
  sector = (size_t *)malloc( (n + 1) * sizeof(size_t) );

  initSector( sector, n );

  for( i = 0; i < dim; i++)
    {
      int num;
      int loc;
      num = bin( i, n );
      loc = sector[ num ];
      states[ loc ] = i;
      sector[ num ] += 1;
    }

  free(sector);
}

/*the numbers of the first states in each different particle number category*/
void initSector( size_t *sector, size_t n )
{
  size_t fact( int i );
  int i;
  size_t index = 0;

  sector[ 0 ] = 0;
  
  /*Initialize the boundary-marking indices*/
  for( i = 1; i <= n; i++)
    {
      sector[ i ] = index + 1;
      
      index += fact( n ) / fact( n - i ) / fact( i );

    }
}

/*factorial function*/
size_t fact( int i )
{
  int j;
  size_t result = 1;
  for( j = 1; j <= i; j++ )
    {
      result *= j;
    }
  return result;
}

/*the dimensions of each particle number sector*/
void dims( size_t *dimensions, size_t n )
{
  size_t i;
  size_t *sector;
  sector = (size_t *)malloc( (n + 1) * sizeof( size_t ) );

  initSector( sector, n );

  for( i = 0; i < n; i++ )
    {
      dimensions[ i ] = sector[ i + 1 ] - sector[ i ];
      //printf("dimensions[%d] = %zu\n", i, dimensions[ i ]);
    }
  dimensions[ n ] = 1;
  
  free( sector );
}

void locations( size_t *pointers, size_t *dimensions, size_t n )
{
  size_t loc = 0;
  size_t i;
  for( i = 0; i <= n; i++ )
    {
      pointers[ i ] = loc;
      loc += dimensions[ i ] * dimensions[ i ];
      //printf("pointers[%zu] = %zu\n", i, pointers[ i ]);
    }
  
}
