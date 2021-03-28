/* routines in common.c */

int **allocate2d(int,int);
int ***allocate3d(int,int,int);
double **allocatedouble2d(int,int);
void deallocate2d(int**,int,int);
void deallocate3d(int***,int,int,int);
void deallocatedouble2d(double**,int,int);
void initneighbor(void);
void initconf(void);
void printconf(void);
int checkconf(void);
void storeconf(FILE *);
int checkflux(FILE *);
int scan(int A[][VOL],int ***B,int);
int calc_sign(int A[], int, int);


void conststates(void);
void constH(int **,int);
void constFermH(int **,int);
void teval(int );
void diagH(int **,int);
void ch(int **,int);
void ch1(int **,int);
void trans(int **,int **,int);
void diagTrans(double *,double *,double **,int **,int **,int);
void measurePOL(int, int);
void measureWLOOP(int, int);
double getloopvalue(int, int);
