#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "constants.h"

double observable( double *h , size_t *pointers, size_t *dimensions, size_t *states, size_t n, size_t dim)
{
  int sublat( int i );

  size_t i, j, k, l;
  double denom = 0.0;
  double num = 0.0;
  int count = 0;

  for( l = 0; l < n + 1; l++)
    {
      printf("l is %d \n", l);
      for( i = 0; i < dimensions[ l ]; i++)
	{
	  denom += h[pointers[l] + i + dimensions[l] * i ];
	  //printf("\n added %.3f\n", h[pointers[l] + i + dimensions[l] * i ]);
	}
    }

  printf("\n\nthe partition function is %10.3f\n\n", denom);
  for( l = 0; l < n + 1; l++ )
    {
      printf("l is %d \n", l);
      for( k = 0; k < dimensions[ l ]; k++ )
	{
	  for( i = 0; i < 1; i++ )
	    {
	      for( j = LX/2; j < LX/2+1; j++ )
		{
		  num += (1 - 2 * sublat((int)i) ) * (1 - 2 * sublat((int)j) )
		    * ( (states[count + k] / (1 << (i % n))) % 2 - 1.0 / 2)
		    * ( (states[count + k] / (1 << (j % n))) % 2 - 1.0 / 2 )
		    * h[ pointers[ l ] + k + dimensions[l] * k ];
		 
		}
	    }
	}
      count += dimensions[ l ];
    }
  printf("\n denominator is %.3f \n", denom);
  printf("\n numerator is %.3f \n", num);

  return num / denom;
}

int sublat( int i )
{
  int i_x, i_y;

  i_x = i % LX;
  i_y = ( i / LX ) % LY;

  //printf("sublat(%d) = %d \n", i, ((i_x + i_y) % 2));

  return ((i_x + i_y) % 2);
}
