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
! Consider the matrix A
!                 |  5   2   1   1   0   0   0   0   0   0   0  |
!                 |  2   6   3   1   1   0   0   0   0   0   0  |
!                 |  1   3   6   3   1   1   0   0   0   0   0  |
!                 |  1   1   3   6   3   1   1   0   0   0   0  |
!                 |  0   1   1   3   6   3   1   1   0   0   0  |
!    A    =       |  0   0   1   1   3   6   3   1   1   0   0  |,
!                 |  0   0   0   1   1   3   6   3   1   1   0  |
!                 |  0   0   0   0   1   1   3   6   3   1   1  |
!                 |  0   0   0   0   0   1   1   3   6   3   1  |
!                 |  0   0   0   0   0   0   1   1   3   6   2  |
!                 |  0   0   0   0   0   0   0   1   1   2   5  |
!
! stored as sparse matrix.
! B is a unit matrix:
!                 |  1   0   0   0   0   0   0   0   0   0   0  |
!                 |  0   1   0   0   0   0   0   0   0   0   0  |
!                 |  0   0   1   0   0   0   0   0   0   0   0  |
!                 |  0   0   0   1   0   0   0   0   0   0   0  |
!                 |  0   0   0   0   1   0   0   0   0   0   0  |
!    B    =       |  0   0   0   0   0   1   0   0   0   0   0  |.
!                 |  0   0   0   0   0   0   1   0   0   0   0  |
!                 |  0   0   0   0   0   0   0   1   0   0   0  |
!                 |  0   0   0   0   0   0   0   0   1   0   0  |
!                 |  0   0   0   0   0   0   0   0   0   1   0  |
!                 |  0   0   0   0   0   0   0   0   0   0   1  |
!
!  In what follows the symbol ' represents a transpose operation.
!
!  The test performs the following operations:
!
!  Step 1. Calls FEASTINIT to define the default values for the input
!          FEAST parameters.
!
!  Step  2. The code solves the standard eigenvalue problem Ax=ex using
!          DFEAST_SCSREV.
!
!  Step  3. The code computes the residual R(i) = | E(i) - Eig(i) | where Eig(i)
!          are the expected eigenvalues and E(i) are eigenvalues computed
!          by DFEAST_SCSREV().
!
!  Step 4. The code computes the maximum absolute value of elements
!          of the matrix Y = X' * X - I, where X is the matrix of eigenvectors
!          computed by DFEAST_SCSREV.
!          DGEMM (BLAS Level 3 Routine) is called to compute X' * X.
!
!  Step  5. The code solves the generalized eigenvalue problem Ax=eBx using
!          DFEAST_SCSRGV.
!
!  Step  6. The code computes the residual R(i) = | E(i) - Eig(i) | where Eig(i)
!          are the expected eigenvalues  and E(i) are eigenvalues computed
!          by DFEAST_SCSRGV().
!
!  Step 7. The code computes the maximum absolute value of the elements of
!          the matrix Y = X' * X - I, where X is the matrix of eigenvectors
!          computed by DFEAST_SCSRGV.
!          DGEMM (BLAS Level 3 Routine) is called to compute X' * X.
!
!*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl_solvers_ee.h"

#define max(a, b) (a) < (b) ? (b): (a)

int main()
{
    char          UPLO = 'F'; /* Type of matrix: (F means full matrix, L/U - lower/upper triangular part of matrix) */
    /* Matrix A of size N in CSR format. Size N and 3 arrays are used to store matrix in CSR format */
/*    const MKL_INT N = 11;
    MKL_INT       rows[12] = { 1, 5, 10, 16, 23, 30, 37, 44, 51, 57, 62, 66 };
    MKL_INT       cols[65] = {    1,   2,   3,   4,
                                  1,   2,   3,   4,   5,
                                  1,   2,   3,   4,   5,   6,
                                  1,   2,   3,   4,   5,   6,   7,
                                       2,   3,   4,   5,   6,   7,   8,
                                            3,   4,   5,   6,   7,   8,   9,
                                                 4,   5,   6,   7,   8,   9,  10,
                                                      5,   6,   7,   8,   9,  10,  11,
                                                           6,   7,   8,   9,  10,  11,
                                                                7,   8,   9,  10,  11,
                                                                     8,   9,  10,  11
                            };

    double        val[65] = {   5.0, 2.0, 1.0, 1.0,
                                2.0, 6.0, 3.0, 1.0, 1.0,
                                1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                     1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                          1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                               1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                                    1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                                         1.0, 1.0, 3.0, 6.0, 3.0, 1.0,
                                                              1.0, 1.0, 3.0, 6.0, 2.0,
                                                                   1.0, 1.0, 2.0, 5.0
                            };
 */
     const MKL_INT N = 6;
     MKL_INT rows[7] = {1,3,7,9,11,15,17};
     MKL_INT cols[16] = {2,5,1,3,4,6,2,5,2,5,1,3,4,6,2,5};
     double  val[16] = {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};


    /* Matrix B of size N in CSR format. Size N and 3 arrays are used to store matrix in CSR format */
/*    MKL_INT       rowsb[12] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
    MKL_INT       colsb[11] = { 1,
                                    2,
                                        3,
                                            4,
                                                5,
                                                    6,
                                                        7,
                                                            8,
                                                                9,
                                                                    10,
                                                                        11
                              };

    double        valb[11] = {  1.0,
                                    1.0,
                                        1.0,
                                            1.0,
                                                1.0,
                                                    1.0,
                                                        1.0,
                                                            1.0,
                                                                1.0,
                                                                    1.0,
                                                                        1.0
                            };
*/
    /* Declaration of FEAST variables */
    MKL_INT      fpm[128];      /* Array to pass parameters to Intel MKL Extended Eigensolvers */
    double       Emin, Emax;    /* Lower/upper bound of search interval [Emin,Emax] */

    double       epsout;        /* Relative error on the trace */
    MKL_INT      loop;          /* Number of refinement loop */
    MKL_INT      L = 6;
    MKL_INT      M0;            /* Initial guess for subspace dimension to be used */
    MKL_INT      M;             /* Total number of eigenvalues found in the interval */

    double       E[6];         /* Eigenvalues */
    double       X[36];        /* Eigenvectors */
    double       res[6];       /* Residual */

    /* Declaration of local variables */
    MKL_INT      info;          /* Errors */
    //double       Eig[11];       /* Eig - array for storing exact eigenvalues */
    //double       R[11];         /* R = |E-Eig| */
    double        Y[6][6];     /* Y=(X')*X-I */

    char         DGEMMC = 'T';  /* Character for GEMM routine, transposed case */
    char         DGEMMN = 'N';  /* Character for GEMM routine, non-transposed case */
    double       one = 1.0;     /* alpha parameter for GEMM */
    double       zero = 0.0;    /* beta  parameter for GEMM */
    MKL_INT      ldx = 6;      /* Leading dimension for source arrays in GEMM */
    MKL_INT      ldy = 6;      /* Leading dimension for destination array in GEMM */

    MKL_INT      i, j;
    double       trace, smax, eigabs;

    /* Exact eigenvalues in range (3.0, 7.0) */
    //Eig[0] = 3.1715728752538100;
    //Eig[1] = 4.0000000000000000;
    //Eig[2] = 4.0000000000000000;
    //Eig[3] = 4.1292484841890931;
    //Eig[4] = 4.4066499006731521;
    //Eig[5] = 6.0000000000000000;
    //for (i=6; i<N; i++)
    // {
    //    Eig[i] = 0.0;
    //}

    printf("\n FEAST DFEAST_SCSREV AND DFEAST_SCSRGV EXAMPLE\n");
    /* Initialize matrix X */
    for (i=0; i<N*N; i++)
    {
        X[i] = zero;
    }

    printf("Sparse matrix size %i\n", (int)N);

    /* Search interval [Emin,Emax] */
    Emin = -5.0;
    Emax = 5.0;
    printf("Search interval [ %.15e, %.15e  ]  \n", Emin, Emax);

    M0   = L;
    M    = L;
    loop = 0;
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

    /* Step 3. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
    * are the expected eigenvalues and E(i) are eigenvalues computed by DFEAST_SCSREV(). */
    printf("Number of eigenvalues found %d \n", M);
    //printf("Computed      |    Expected    \n");
    //printf("Eigenvalues   |    Eigenvalues \n");
    eigabs = 0.0;
    for (i=0; i<M; i++)
    {
        //R[i] = fabs(E[i] - Eig[i]);
        //eigabs = max(eigabs, R[i]);
        //printf("%.15e %.15e \n", E[i], Eig[i]);
        printf("%.15e \n", E[i]);
    }
    //printf("Max value of | computed eigenvalue - expected eigenvalues | %.15e \n", eigabs);

    /* Step 4. The code computes the maximum absolute value of elements
     * of the matrix Y = X' *X - I, where X is the matrix of eigenvectors
     *  computed by DFEAST_SCSREV.
     *
     * Call BLAS to compute Y = X' * X  */

    dgemm(
        &DGEMMC, /* IN: 'T', transposed case*/
        &DGEMMN, /* IN: 'N', non-transposed case*/
        &M,      /* IN: Number of rows in matrix Y */
        &M,      /* IN: Number of columns in matrix X */
        &N,      /* IN: Number of columns in matrix Y */
        &one,    /* IN: alpha = 1.0 */
        X,       /* IN: Source #1 for GEMM, will be transposed */
        &ldx,    /* IN: Leading dimension of Source 1 */
        X,       /* IN: Source #2 for GEMM */
        &ldx,    /* IN: Leading dimension of Source 2 */
        &zero,   /* IN: beta = 0.0 */
        Y,       /* OUT: Destination */
        &ldy     /* IN: Leading dimension of Destination */
        );

    /* Compute Y = Y - I */
    for (i=0; i<M; i++)
    {
        Y[i][i] -= 1.0;
    }

    printf("*************************************************\n");
    printf("************** REPORT ***************************\n");
    printf("*************************************************\n");
    printf("#Search interval [Emin,Emax] %.15e %.15e\n", Emin, Emax);
    printf("#mode found/subspace %d %d \n", M, M0);
    printf("#iterations %d \n", loop);
    trace = 0.0;
    for (i=0; i<M; i++)
    {
        trace += E[i];
    }
    printf("TRACE %.15e \n", trace);
    printf("Relative error on the Trace %.15e \n", epsout);
    printf("Index/Eigenvalues/Residuals\n");
    for (i=0; i<M; i++)
    {
        printf("   %d  %.15e %.15e \n",i, E[i], res[i]);
    }
    smax = 0.0;
    for (i=0; i<M; i++)
    {
        for (j=0; j<M; j++)
        {
            smax = max(smax, fabs(Y[i][j]));
        }
    }
    printf( "Max(X' * X - I) = %.15e \n", smax);
    
    /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
    /*!!!!!!!!!!!!!!! GENERALIZED EIGENVALUE PROBLEM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
    /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

    /* Reset initial parameters for new problem */
 /*   M0 = L;
    for (i=0; i<N; i++)
    {
        E[i] = 0.0;
    }
    for (i=0; i<N*N; i++)
    {
        X[i] = zero;
    }
*/
    /* Step 5. Solve the generalized eigenvalue problem Ax=eBx by DFEAST_SCSRGV */
//    printf("Testing dfeast_scsrgv  \n");

//    dfeast_scsrgv(
//        &UPLO,   /* IN: UPLO = 'F', stores the full matrix */
//        &N,      /* IN: Size of the problem */
//        val,     /* IN: CSR matrix A, values of non-zero elements */
//        rows,    /* IN: CSR matrix A, index of the rist non-zero in row */
//        cols,    /* IN: CSR matrix A, columns indeces for each non-zero element */
//        valb,    /* IN: CSR matrix B, values of non-zero elements */
//        rowsb,   /* IN: CSR matrix B, index of the rist non-zero in row */
//        colsb,   /* IN: CSR matrix A, columns indeces for each non-zero element */
//        fpm,     /* IN: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
//        &epsout, /* OUT: Relative error of on the trace */
//        &loop,   /* OUT: Contains the number of refinement loop executed */
//        &Emin,   /* IN: Lower bound of search interval */
//        &Emax,   /* IN: Upper bound of search interval */
//        &M0,     /* IN/OUT: The initial guess for subspace dimension to be used. */
//        E,       /* OUT: The first M entries of Eigenvalues */
//        X,       /* OUT: The first M entries of Eigenvectors */
//        &M,      /* OUT: The total number of eigenvalues found in the interval */
//        res,     /* OUT: The first M components contain the relative residual vector */
//        &info    /* OUT: Code of error */
//        );

//    printf("FEAST OUTPUT INFO %d \n" ,info);
//    if ( info != 0 )
//    {
//        printf("Routine dfeast_scsrgv return error: %i", (int)info);
//        return 1;
//    }
    /* Task 6. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
     *  are the expected eigenvalues  and E(i) are eigenvalues computed
     *  by  DFEAST_SCSRGV */
/*    printf("Number of eigenvalues found %d \n", M);
    printf("Computed      |    Expected    \n");
    printf("Eigenvalues   |    Eigenvalues \n");
    eigabs = 0.0;
    for (i=0; i<M; i++)
    {
        R[i] = fabs(E[i] - Eig[i]);
        eigabs = max(eigabs, R[i]);
        printf("%.15e %.15e \n", E[i], Eig[i]);
    }
    printf("Max value of | computed eigenvalue - expected eigenvalues | %.15e \n", eigabs);
*/
    /* Step 7. The code computes the maximum absolute value of the elements of
     * the matrix  Y = X' * X - I, where X is the matrix of eigenvectors
     * computed by DFEAST_SCSRGV.
     *
     * Call BLAS to compute X' * X */
//    dgemm(
//        &DGEMMC, /* IN: 'T', transposed case*/
//        &DGEMMN, /* IN: 'N', non-transposed case*/
//        &M,      /* IN: Number of rows in matrix Y */
//        &M,      /* IN: Number of columns in matrix X */
//        &N,      /* IN: Number of columns in matrix Y */
//        &one,    /* IN: alpha = 1.0 */
//        X,       /* IN: Source #1 for GEMM, will be transposed */
//        &ldx,    /* IN: Leading dimension of Source 1 */
//        X,       /* IN: Source #2 for GEMM */
//        &ldx,    /* IN: Leading dimension of Source 2 */
//        &zero,   /* IN: beta = 0.0 */
//        Y,       /* OUT: Destination */
//        &ldy     /* IN: Leading dimension of Destination */
//        );

    /* Compute Y = Y - I */
//    for (i=0; i<M; i++)
//    {
//        Y[i][i] -= 1.0;
//    }

    /* Check the orthogonality of X' * X */
/*    smax = 0.0;
    for (i=0; i<M; i++)
    {
        for (j=0; j<M; j++)
        {
            smax = max(smax, fabs(Y[i][j]));
        }
    }

    printf("*************************************************\n");
    printf("************** REPORT ***************************\n");
    printf("*************************************************\n");
    printf("# Search interval [Emin,Emax] %.15e %.15e\n", Emin, Emax);
    printf("# mode found/subspace %d %d \n",M,M0);
    printf("# iterations %d \n",loop);
    printf("Max(X' * X - I) = %.15e \n", smax);
    trace = 0.0;
    for (i=0; i<M; i++)
    {
        trace += E[i];
    }
    printf("TRACE %.15e \n", trace);
    printf("Relative error on the Trace %.15e \n", epsout);
    printf("Index/Eigenvalues/Residuals\n");
    for (i=0; i<M; i++)
    {
        printf("   %d  %.15e %.15e \n", i, E[i], res[i]);
    }
*/
    return 0;
}
