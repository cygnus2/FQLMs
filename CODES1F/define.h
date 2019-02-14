#ifndef __DEFINE_H_INCLUDED__
#define __DEFINE_H_INCLUDED__

#include <vector>
#define DIM 2
#define PI 2*asin(1.0)

extern int *next[2*DIM+1];
extern int *nextCHK[2*DIM+1];
extern int *chk2lin,*lin2chk;
extern int LX,LY,VOL,VOL2;
extern double lam,Ti,Tf,dT;
/* NTOT = total no of basis states 
 * NH   = states not killed by H  */
extern int NTOT,NH;

extern std::vector<std::vector<bool>> basis;
extern std::vector<std::vector<bool>> basis_nonflip;
extern std::vector<std::vector<bool>> basis_flip;
extern std::vector<int> Oflip;

/* routines */
void initneighbor(void);
void conststates(void);
void constH(std::vector<double>&, std::vector<std::vector<double>>&);
void constH_zeroL(std::vector<double>&, std::vector<std::vector<double>>&);
void evolve_cartoons(std::vector<double>&, std::vector<std::vector<double>>&);
void calc_Oflip(std::vector<double>&, std::vector<std::vector<double>>&);
int scan(std::vector<bool>&);
#endif 
