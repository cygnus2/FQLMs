#define DIM 3
#define PI 2*asin(1.0)

int *conf[DIM];
int *next[2*DIM+1];
int LX,LY,LZ,VOL;
double lam;
/* N is the total no of states while N1 is the *
 * no of states which have non-trivial flux    */
int N,N1,N2;
/* winding number sectors */
double Wx, Wy, Wz;
/* allocate memory for states */
int *states;
/* The Hamiltonian */
int **H;
/* Charge conjugation matrix */
//int **cnj,**cnj1;
/* Translation matrix */
//int **tx, **ty;
/* states and spins */
//int ***stspin;
