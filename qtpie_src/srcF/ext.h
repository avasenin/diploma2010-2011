/***************************************************************************
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/*         
  *         -------------------------------------------------
  *                LINEAR LEAST SQUARE SOLVER
  *         -------------------------------------------------
  *
  * 
  * This file has implementation of the LAPACK routine 'dgels', 'dgelsy' and 'dgelss'
  * for C/C++. This program solves the least square problem defined as follows:
  *
  * minimize  || Ax - B ||
  *     x 
  *     
  *  ---------------------------------------------   
  *  Date        : 30 December 2006, Saturday
  *  ---------------------------------------------
  * please refer to man pages of 'dgels' for details.   
  *
  * void dgels(double **A, double *b, int m, int n, double *x);
  * void dgelsy(double **A, double *b, int m, int n, double *x);
  * void dgelss(double **A, double *b, int m, int n, double *x);
  * void dgelss_pi(double **A, double *b, int m, int n, double *x);
  *
  * 
  * A    : m by n matrix
  * b    : m by 1 matrix 
  * x    : n by 1 matrix for returned values
  *
  *
  *  dgels     : computes minimum norm solution when rank(A) = min(m,n)
  *  dgelsy    : computes minimum norm solution when rank(A) <= min(m,n) 
  *              uses complete orthogonal factorization of A
  *  dgelss    : computes minimum norm solution when rank(A) <= min(m,n) uses SVD
  *  dgelss_pi : computes the solution x = pinv(A)*b which is the NOT the minimum norm solution
  *
  *
  * ----------------------------------
  * Date: 04 September 2007 Tuesday
  * ------------------------------------
  * Least Square Solution for Multiple columns
  *
  * minimize || A X - B ||
  *    X
  * 
  *
  *  Author      : Swagat Kumar
  *  Email       : swagatk@iitk.ac.in / swagat.kumar@gmail.com
  *  Affiliation : Indian Institute of Technology Kanpur, India 
  *  License     : GPL (GNU Public License)
  *  --------------------------------------------------------------*/
#include <math.h>
#include<cblas.h>

#define MAXIMUM(m, n) (m) < (n) ? (n) : (m)
#define MINIMUM(m, n) (m) < (n) ? (m) : (n)

void dgels(double **A, double *b, int m, int n, double *x);
void dgelsy(double **A, double *b, int m, int n, double *x);
void dgelss(double **A, double *b, int m, int n, double *x);
void dgelss_pi(double **A, double *b, int m, int n, double *x);
void dgelss(double **A, double **B, int m, int n, int p, double **X);

double* dgels_ctof(double **in, int rows, int cols);
 
extern "C" void dgels_(char *trans, int *m, int *n, int *nrhs, double *a, int *lda,
		       double *b, int *ldb, double *work, int *lwork, int *info);

extern "C" void dgelsy_(int *m, int *n, int *nrhs, double *a, int *lda,
		       double *b, int *ldb, int *JPVT, double *RCOND, int *RANK,
             double *work, int *lwork, int *info);
 
extern "C" void dgelss_(int *m, int *n, int *nrhs, double *a, int *lda,
		       double *b, int *ldb, double *s, double *RCOND, int *RANK,
             double *work, int *lwork, int *info);
 
//------------------------------------------------------------------------
// Maximum of 3 numbers
int max3(int a, int b, int c)
{
  if( a >= b)
  {
    if(a >= c)
      return a;
    else
      return c;
  }
  else
  {
    if(b >= c)
      return b;
    else
      return c;
  }
}
  
//--------------------------------------------------------------
void dgels(double **A, double *b, int m, int n, double *x)
{
  int nrhs, lda, ldb, info, lwork;
  double *a;                      
  double *work;

  
  nrhs = 1;

  lwork = m + n;


  lda = m; // The leading dimension of A
  
  a = dgels_ctof(A, lda, n); // Convert A to a Fortran style matrix

  ldb = MAXIMUM(m, n); // leading dimension of B 


  work = new double [lwork];
  
  /* The Fortran routine replaces the input information in b with the
     results.  Since we are interested in the results, we put the
     initial conditions in the output array and use that in the
     call. */

  for (int i = 0; i < ldb; i++)
      x[i] = b[i];

  // Now the function call

  dgels_("N", &m, &n, &nrhs, a, &lda, x, &ldb, work, &lwork, &info);

  if(info < 0)
    printf("Lapack routine dgels returns error ...\n");
  
  // Clean up the memory before returning to the calling program
delete [] work;
}
//---------------------------------------------------------------------
void dgelsy(double **A, double *b, int m, int n, double *x)
{
  int nrhs, lda, ldb, info, lwork, rank, *jpvt;
  double *a, *work, rcond;


  jpvt = new int[n];
  rcond = 0.00001;
  
  nrhs = 1;

  lwork = 3*(m + n);


  lda = m; // The leading dimension of A
  
  a = dgels_ctof(A, lda, n); // Convert A to a Fortran style matrix

  ldb = MAXIMUM(m, n); // leading dimension of B 


  work = new double [lwork];
  
  /* The Fortran routine replaces the input information in b with the
     results.  Since we are interested in the results, we put the
     initial conditions in the output array and use that in the
     call. */

  for (int i = 0; i < ldb; i++)
      x[i] = b[i];

  // Now the function call

  dgelsy_(&m, &n, &nrhs, a, &lda, x, &ldb, jpvt, &rcond, &rank, work, &lwork, &info);

  if(info < 0)
    printf("Lapack routine dgels returns error ...\n");
  
  // Clean up the memory before returning to the calling program

  delete a;              
  delete b;
  delete jpvt;
}
// ===================================================================
void dgelss(double **A, double *b, int m, int n, double *x)
{
  int nrhs, lda, ldb, info, lwork, rank, sdim;
  double *a, *work, *s, rcond;
  double **X, *B;

  int i, j;

  sdim = MINIMUM(m,n);
  s = new double [sdim];         // Singular values of A returned by dgelss_

  rcond = 1e-6;               // Condition number ... choose a small value 

  ldb = MAXIMUM(m, n);               // leading dimension of B 

  nrhs = 1;                      // no. of columns of B

  lwork = 30 * m * n;            // Allocate memory of the workspace ..
  // In case you get memory related error, refer to manpage
  // or increase this number.


  lda = m;                       // The leading dimension of A

  a = dgels_ctof(A, lda, n);     // Convert A to a Fortran style matrix



  work = new double [lwork];     // workspace memory to be used by LAPACK routine

  /* The Fortran routine replaces the input information in b with the
     results.  Since we are interested in the results, we put the
     initial conditions in the output array and use that in the
     call. */

  for(i = 0; i < ldb; i++)
    x[i] = b[i];

  // Now the function call

  dgelss_(&m, &n, &nrhs, a, &lda, x, &ldb, s, &rcond, &rank, work, &lwork, &info);

  if(info < 0)
  {
    printf("Lapack routine dgelss returns error ...\n");
    printf("Invalid argument: %d\n", info);

  }
  else if(info > 0)
  {
    printf("SVD fails to converge and returns: %d\n", info);
    printf("see man page for details ...\n");
  }
  else 
  {
    printf("Rank of matrix A: %d\n", rank);
    printf("Singular values of A: ");
    for(i = 0; i < sdim; i++)
      printf("%f\t", s[i]);
  }  
  printf("\n\n");    


  // Clean up the memory before returning to the calling program

  delete [] s;
  delete [] work;
}
// -------------------------------------------------------------------
// Solves Ax = B such that x = pinv(A) * B
// -----------------------------------------------------------------
void dgelss_pi(double **A, double *b, int m, int n, double *x)
{
  int nrhs, lda, ldb, info, lwork, rank, sdim, k;
  double *a, *work, *s, rcond;
  double **B, *br;
  
  int i, j;
  
  sdim = MINIMUM(m,n);

  rcond = 1e-6;                  // Condition number ... choose a small value 
  
  ldb = MAXIMUM(m, n);               // leading dimension of B 
  lda = m;                       // The leading dimension of A
  
  nrhs = m;                      // no. of columns of B

  s = new double [sdim];         // Singular values of A returned by dgelss_

  B = new double *[ldb];
  for(i = 0; i < ldb; i++)
    B[i] = new double [nrhs];
  
  for(i = 0; i < ldb; i++)
    for(j = 0; j < nrhs; j++)
    {
      if(i == j)
         B[i][j] = 1.0;
      else 
        B[i][j] = 0.0;
    }
      
  

  lwork = 30 * m * n;            // Allocate memory of the workspace ..
                                 // refer manpage for details


  
  a = dgels_ctof(A, lda, n);     // Convert A to a Fortran style matrix
  br = dgels_ctof(B, ldb, nrhs); // Convert B to Fortran style matrix



  work = new double [lwork];     // workspace memory to be used by LAPACK routine
  
  
  // compute pinv(A) using SVD

  dgelss_(&m, &n, &nrhs, a, &lda, br, &ldb, s, &rcond, &rank, work, &lwork, &info);
  /* 'br' contains the pseudo-inverse of A */
  
  if(info < 0)
  {
    printf("Lapack routine dgelss returns error ...\n");
    printf("Invalid argument: %d\n", info);
    
  }
  else if(info > 0)
  {
    printf("SVD fails to converge and returns: %d\n", info);
    printf("see man page for details ...\n");
    exit(-1);
  } 

  // COmpute the solution x = pinv(A) * b
  for(i = 0; i < ldb; i++)
    for(j = 0; j < nrhs; j++)
      B[i][j] = br[i + j * ldb];

  for(i = 0; i < n; i++)
  {
    x[i] = 0.0;
    for(j = 0; j < m; j++)
      x[i] += B[i][j] * b[j];
  } 


  
  // Clean up the memory before returning to the calling program

  delete [] s;
  delete [] work;
  delete [] B;
}
// ===================================================================
// Date: 04 September 2007 Tuesday
// A : M x N
// X : N x P
// B : M x P
// solve : A * X = B
//
// NOTE: This routine Overwrites B Matrix.
// ================================================================
void dgelss(double **A, double **B, int m, int n, int p, double **X)
{
  int nrhs, lda, ldb, info, lwork, rank, sdim;
  double *a, *work, *s, rcond;
  double *x;
  
  sdim = MINIMUM(m,n);
  s = new double [sdim];         // Singular values of A returned by dgelss_

  rcond = 1e-16;                  // Condition number ... choose a small value 
  
  ldb = MAXIMUM(1, MAXIMUM(m, n));  // leading dimension of B 
  
  nrhs = p;                      // no. of columns of B and X

  lwork = 3 * MINIMUM(m, n) + max3(2*MINIMUM(m,n), MAXIMUM(m,n), nrhs) + 10;      

                                 // Allocate memory of the workspace ..
                                 // In case you get memory related error, refer to manpage
                                 // or increase this number.


  lda = m;                       // The leading dimension of A
  
  a = dgels_ctof(A, lda, n);     // Convert A to a Fortran style matrix



  work = new double [lwork];     // workspace memory to be used by LAPACK routine
  

  /* Note: B is M x P matrix. dgelss_ overwrites B with the resulting
   * N x P solution matrix X. We initialize X with B and hence DGELSS
   * replaces B with the solution X. */

  // Input Matrix B
  
  x = dgels_ctof(B, ldb, nrhs);     // Convert X to a Fortran style matrix

  // Now the function call

  dgelss_(&m, &n, &nrhs, a, &lda, x, &ldb, s, &rcond, &rank, work, &lwork, &info);

  if(info < 0)
  {
    printf("Lapack routine dgelss returns error ...\n");
    printf("Invalid argument: %d\n", info);
    
  }
  else if(info > 0)
  {
    printf("SVD fails to converge and returns: %d\n", info);
    printf("see man page for details ...\n");
  }
 
  /*for(int i = 0; i < 10; i++)
    cout << x[i] << "\t";
  cout << endl;
  getchar();*/

  // Return the Solution Matrix X (N x NRHS)
  for(int i = 0; i < ldb; i++)
    for(int j = 0; j < nrhs; j++)
      B[i][j] = x[i + j * ldb];  // Column Major mode

  for(int i = 0; i < n; i++)
    for(int j = 0; j < nrhs; j++)
      X[i][j] = B[i][j];
  
  // Clean up the memory before returning to the calling program

  delete [] s;
  delete [] work;
  delete [] a;
  delete [] x;
}
// -------------------------------------------------------------------
 
// ------------------------------------------------------
// Lapack routines uses Matrices in Column-Major mode and hence need to be
// transformed accordingly. 
//
double* dgels_ctof(double **in, int rows, int cols)
{
  double *out;
  int i, j;

  out = new double[rows*cols];
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++) 
      out[i+j*rows] = in[i][j]; // column-major mode
  return(out);
} 

