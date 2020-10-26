/*******************************************************************************
* Copyright 2005-2016 Intel Corporation All Rights Reserved.
*
* The source code,  information  and material  ("Material") contained  herein is
* owned by Intel Corporation or its  suppliers or licensors,  and  title to such
* Material remains with Intel  Corporation or its  suppliers or  licensors.  The
* Material  contains  proprietary  information  of  Intel or  its suppliers  and
* licensors.  The Material is protected by  worldwide copyright  laws and treaty
* provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
* modified, published,  uploaded, posted, transmitted,  distributed or disclosed
* in any way without Intel's prior express written permission.  No license under
* any patent,  copyright or other  intellectual property rights  in the Material
* is granted to  or  conferred  upon  you,  either   expressly,  by implication,
* inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
* property rights must be express and approved by Intel in writing.
*
* Unless otherwise agreed by Intel in writing,  you may not remove or alter this
* notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
* suppliers or licensors in any way.
*******************************************************************************/

/*
!   Content: Example for Intel MKL Extended Eigensolvers (sparse format,
!            double precision)
!
!*******************************************************************************
!
! The following routines are used in the example:
!          DGEMM  DFEAST_SCSREV  DFEAST_SCSRGV  FEASTINIT
!
!  The test performs the following operations:
!
!  Step 1. Calls FEASTINIT to define the default values for the input
!          FEAST parameters.
!
!  Step  2. The code solves the standard eigenvalue problem Ax=ex using
!          DFEAST_SCSREV.
!
!*******************************************************************************/
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl_solvers_ee.h"
#include <fstream>

#define max(a, b) (a) < (b) ? (b): (a)

int main()
{
    char          UPLO = 'F'; /* Type of matrix: (F means full matrix, L/U - lower/upper triangular part of matrix) */
    /* Matrix A of size N in CSR format. Size N and 3 arrays are used to store matrix in CSR format */
    
    //MKL_INT rows[7] = {1,3,7,9,11,15,17};
    // MKL_INT cols[16] = {2,5,1,3,4,6,2,5,2,5,1,3,4,6,2,5};
    // double  val[16] = {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};

    unsigned int k;
    /* read matrix A from a binary file */
    std::ifstream inFile ("sparse_mat.bin", std::ios::in | std::ios::binary); 
    int nrows, ncols;
    inFile.read((char*)&nrows,sizeof(int));
    inFile.read((char*)&ncols,sizeof(int));
    std::cout<<"No of rows = "<<nrows<<std::endl;
    std::cout<<"No of cols = "<<ncols<<std::endl;

    const MKL_INT N = nrows-1;
    MKL_INT *rows,*cols;
    double *val;
    rows = (MKL_INT*)(calloc(nrows,sizeof(MKL_INT)));
    cols = (MKL_INT*)(calloc(ncols,sizeof(MKL_INT)));
    val =  (double*)(calloc(ncols,sizeof(double)));
    inFile.read((char*)rows, nrows*sizeof(MKL_INT));
    inFile.read((char*)cols, ncols*sizeof(MKL_INT));

    for(k=0;k<ncols;k++) val[k]=-1.0;

    printf("row[10]=%d\n",rows[10]);
    printf("column[10]=%d\n",cols[10]);

    /* Declaration of FEAST variables */
    MKL_INT      fpm[128];      /* Array to pass parameters to Intel MKL Extended Eigensolvers */
    double       Emin, Emax;    /* Lower/upper bound of search interval [Emin,Emax] */

    double       epsout;        /* Relative error on the trace */
    MKL_INT      loop;          /* Number of refinement loop */
    //MKL_INT      L = N;
    MKL_INT      L = 100;
    MKL_INT      M0;            /* Initial guess for subspace dimension to be used */
    MKL_INT      M;             /* Total number of eigenvalues found in the interval */
    double   *E;         /* Eigenvalues */
    double   *X;       /* Eigenvectors */
    double   *res;       /* Residual */
    E = (double*)(calloc(N,sizeof(double)));
    X = (double*)(calloc(N*N,sizeof(double)));
    res = (double*)(calloc(N,sizeof(double)));


    /* Declaration of local variables */
    MKL_INT      info;          /* Errors */
    double       one = 1.0;     /* alpha parameter for GEMM */
    double       zero = 0.0;    /* beta  parameter for GEMM */
    MKL_INT      i, j;
    double       trace, smax, eigabs;

    printf("\n FEAST DFEAST_SCSREV AND DFEAST_SCSRGV EXAMPLE\n");
    /* Initialize matrix X */
    for (i=0; i<N*N; i++){
        X[i] = zero;
    }

    printf("Sparse matrix size %i\n", (int)N);

    /* Search interval [Emin,Emax] */
    Emin = -16.0;
    Emax = -8.0;
    printf("Search interval [ %.15e, %.15e  ]  \n", Emin, Emax);

    M0   = L;
    M    = L;
    loop = 2;
    info = 0;
    epsout = 0.0;

    /* Step 1. Call  FEASTINIT to define the default values for the input FEAST parameters */
    feastinit(
        fpm /* OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
        );

    fpm[0] =  1; /* Extended Eigensolver routines print runtime status to the screen. */

    /* Step 2. Solve the standard Ax = ex eigenvalue problem. */
    printf("Testing dfeast_scsrev routine:\n");
    dfeast_scsrev(
        &UPLO,   /* IN: UPLO = 'F', stores the full matrix */
        &N,      /* IN: Size of the problem */
        val,     /* IN: CSR matrix A, values of non-zero elements */
        rows,    /* IN: CSR matrix A, index of the first non-zero element in row */
        cols,    /* IN: CSR matrix A, columns indices for each non-zero element */
        fpm,     /* IN/OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
        &epsout, /* OUT: Relative error of on the trace */
        &loop,   /* OUT: Contains the number of refinement loop executed */
        &Emin,   /* IN: Lower bound of search interval */
        &Emax,   /* IN: Upper bound of search interval */
        &M0,     /* IN: The initial guess for subspace dimension to be used. */
        E,       /* OUT: The first M entries of Eigenvalues */
        X,       /* IN/OUT: The first M entries of Eigenvectors */
        &M,      /* OUT: The total number of eigenvalues found in the interval */
        res,     /* OUT: The first M components contain the relative residual vector */
        &info    /* OUT: Error code */
        );
    printf("FEAST OUTPUT INFO %d \n",info);
    if ( info != 0 )
    {
        printf("Routine dfeast_scsrev returns code of ERROR: %i", (int)info);
        return 1;
    }
    printf("Number of eigenvalues found %d \n", M);
    //for (i=0; i<M; i++){
    //    printf("%.15e \n", E[i]);
    //}
    
    free(rows); free(cols);
    return 0;
}
