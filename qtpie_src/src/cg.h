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
#ifndef CG__H
#define CG__H

#include "parameters.h"

inline void lapack_svdsolver(int N, double *A, double *b, double *x)
{
	printf("[ERROR] unexpecting call");

/*	double *S = new double [N]; // Matrix of singular values
	int Rank, stat, WorkSize;
	double *WORK;
	double RCond = 1.0e-8;
	
	for (int i=0; i<N; i++) x[i] = b[i];
	
	// First find optimal workspace size
	WORK = new double [1];
	dgelss(N, N, 1, A, N, x, N, S, RCond, Rank, WORK, -1, stat);
	WorkSize = (int)WORK[0];
	delete [] WORK;
	
	WORK = new double [WorkSize];
	dgelss(N, N, 1, A, N, x, N, S, RCond, Rank, WORK, WorkSize, stat);
	delete [] WORK;
	delete [] S;*/
}

/**
* @brief Use LAPACK routine to calculate condition number of the NxN matrix A
*/
inline double ConditionNumber(int N, double *A)
{
	printf("[ERROR] unexpecting call");
/*	double *c = new double [N];
	double *S = new double [N];
	int Rank, stat, WorkSize;
	double *WORK;
	double RCond = 1.0e-8;
	double ConditionNumber_ = 0.;
	
	for (int i=0; i<N; i++) c[i] = 0.;
	
	// First find optimal workspace size
	WORK = new double [1];
	dgelss(N, N, 1, A, N, c, N, S, RCond, Rank, WORK, -1, stat);
	WorkSize = (int)WORK[0];
	delete [] WORK;
	
	WORK = new double [WorkSize];
	dgelss(N, N, 1, A, N, c, N, S, RCond, Rank, WORK, WorkSize, stat);
	delete [] WORK;
	
	if (Rank == N) ConditionNumber_ = S[0] / S[N-1];
	else
	{
		printf("Matrix found to be singular\n");
		ConditionNumber_ = 0.;
	}
	
	delete [] S;
	delete [] c;
	return ConditionNumber_;*/
}

/**
* @brief Calculates the solution x of the approximate preconditioned problem 
* \f[
*  \mathbf{P}(\mathbf{A}) \vec{x} = \vec{b}
* \f]
*/
inline void Precondition(int N, double *A, double *b, double *x)
{
	double ReciprocalSumOfDiagonals, MatrixElement;
	
	for (int i=0; i<N-1; i++) { x[i] = b[i] / _E(N, A, i, i);
	//printf("i=%d x=%10.5e b=%10.5e A=%10.5e\n", i+1, x[i], b[i], _E(N, A, i, i));
	}
	
	// Add in exact solution for last column and last row
	//   ( 0 v ) ( x ) = (   y v   )
	//    ( v w ) ( y ) = ( x.v + wy) 
	ReciprocalSumOfDiagonals = 0.;
	
	for (int i=0; i<N-1; i++) ReciprocalSumOfDiagonals += 1./ _E(N, A, i, i);
	//printf("ReciprocalSumOfDiagonals=%10.5e\n", ReciprocalSumOfDiagonals);
	x[N-1] = - b[N-1] / ReciprocalSumOfDiagonals;
	//printf("N=%d x=%10.5e\n", N-1, x[N-1]);
	
	for (int i=0; i<N-1; i++)
	{
		MatrixElement = 1./(ReciprocalSumOfDiagonals * _E(N, A, i, i));
		x[i] += b[N-1] * MatrixElement;
		x[N-1] += b[i] * MatrixElement;
	}
	//exit(0);
}

/**
* @brief Double precision conjugate gradient solver with Jacobi preconditioner
*      Solves the matrix problem Ax = b for x
*     Implemented from Golub and van Loan's stuff
* @author Jiahao Chen
* @date   2008-01-28
* @param N : an integer specifying the size of the problem
* @param A : a real, positive definite NxN matrix
* @param b : a real vector with N elements
* @param x : (Output) solution to matrix equation
*     On input, contains initial guess
*/
inline void dcg(int N, double *A, double *b, double *x)
{
	int max_k = 100000; // maximum number of iterations
	double tol = 1.0e-10; // Convergence tolerance
	
	double *r = new double [N];
	double *p = new double [N];
	double *q = new double [N];
	double *z = new double [N];
	for (int i=0; i<N; i++) r[i] = p[i] = q[i] = z[i] = 0.;
	
	double alpha, norm, critical_norm, gamma, gamma0;
	bool Verbose = false;
	
	// Termination criterion norm
	//printf("N=%d tol=%lf %lf %lf %lf\n",  N, tol, b[0], b[1], b[N-1]);
	critical_norm = tol * dnrm2(N, b, 1);
	//printf("critical_norm=%10.5e\n", critical_norm);
	// Calculate initial guess x from diagonal part P(A) x = b 
	// The secret code to want an initial guess calculated is to pass an
	// initial guess with the first entry equal to floating-point zero.
	// If not, we'll just use the pre-specified initial guess that's already in x
//printf("precond N=%d A[0][0]=%10.5e x[0]=%10.5e x[1]=%10.5e\n",  N, _E(N, A, 0, 0), x[0], x[1]);
	if (x[0] == 0.)
	{
//printf("precond N=%d A[0][0]=%10.5e x[0]=%10.5e x[1]=%10.5e\n",  N, _E(N, A, 0, 0), x[0], x[1]);
	Precondition(N, A, b, x);
//printf("precond N=%d A[0][0]=%10.5e x[0]=%10.5e x[1]=%10.5e\n",  N, _E(N, A, 0, 0), x[0], x[1]);
	}
	
	// Calculate residual r = b - Ax
	dcopy(N, b, 1, r, 1);
	// r = r - Ax
	dgemv(CblasNoTrans, N, N, -1., A, N, x, 1, 1., r, 1);
	
	// Calculate norm
	norm = dnrm2(N, r, 1);
	if (Verbose) printf("Iteration 0 : %lf %lf\n", norm, norm/critical_norm);
	
	for (int k=1; k<max_k; k++)
	{
		// Generate preconditioned z from P(A) z = r
		Precondition(N, A, r, z);
		// Propagate old vectors
		gamma0 = gamma;
		gamma  = ddot(N, r, 1, z, 1);
		if (k != 1)
		{
			// p = z + gamma/gamma0 * p
			// z = z + gamma/gamma0 * p
			daxpy(N, gamma/gamma0, p, 1, z, 1);
		}
		dcopy(N, z, 1, p, 1);
		// Form matrix-vector product
		// q = A p
		dgemv(CblasNoTrans, N, N, 1., A, N, p, 1, 0., q, 1);
		// Calculate step size alpha = gamma / p.q
		alpha = gamma / ddot(N, p, 1, q, 1);
		// Propagate by step size
		//c$omp    sections
		//c$omp    section
		// x = x + alpha * p
		daxpy(N, alpha, p, 1, x, 1);
		//c$omp    section
		// r = r - alpha * q
		daxpy(N,-alpha,q,1,r,1);
		//c$omp    end sections
	
		// Calculate new norm of residual
		norm = dnrm2(N, r, 1);
		// If requested, print convergence information
		if (Verbose) printf("Iteration %d : %lf %lf\n", k, norm, norm/critical_norm);
		
		// Check termination criterion. Done if || r || < tol || b ||
		if (norm < critical_norm) goto label1;
	}
	// Oops, reached maximum iterations without convergence
	printf("dcg: maximum iterations reached.\n");
	printf("WARNING: Solution may not be converged.\n");
	
label1:
	if (Verbose) printf("dcg: solution found with residual %lf\n", norm);
	
	delete [] r;
	delete [] p;
	delete [] q;
	delete [] z;
}

/**
* @brief Calculates an approximate inverse to a matrix of the form
*  \f[
*    \mathbf{M}=\left(\begin{array}{cc}\mathbf{J} & 1\\
*        1 & 0\end{array}\right)
*  \f]
*  The inverse is calculated by approximating J by its diagonal, in which
*  case an exact inverse can be constructed.
* @param N Size of matrix
* @param M Matrix to invert
* @param W Approximate inverse matrix
*/
inline void ApproximateInverse(double *M, double *W, int N)
{
	double MatrixElement;
 
	for (int i=0; i<N; i++) 
	for (int j=0; j<N; j++) 
		_E(N, W, i, j) = 0.;
	
	double ReciprocalSum = 0.;
	for (int i=0; i<N; i++)
		ReciprocalSum += 1. / _E(N, M, i, i);
	_E(N, W, N-1, N-1) = -1. / ReciprocalSum;
	
	for (int i=0; i<N; i++)
	{
		MatrixElement = 1. / (ReciprocalSum * _E(N, M, i, i));
		_E(N, W, i, N-1) = MatrixElement;
		_E(N, W, N-1, i) = MatrixElement;
		// This is an approximation to the approximate problem, replacing that which follows
		_E(N, W, i, i) = 1. / _E(N, M, i, i);
	}
}

/**
* @param N      : an integer specifying the size of the problem
* @param A(N,N) : a real, positive definite NxN matrix
* @param b(N)   : a real vector with N elements
* @param x(N)   : (Output) solution to matrix equation
*/
inline void solver(int N, double *A, double *b, double *x)
{
	// Use conjugate gradients routine
	dcg(N, A, b, x);
}

#endif
