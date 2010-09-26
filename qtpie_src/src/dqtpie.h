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
#ifndef DGTPIE__H
#define DGTPIE__H

#include "parameters.h"

/**
* @brief Calculates the matrix to be solved in QTPIE
* \f[
*  \left(\begin{array}{cc}
*  \mathbf{J} & \vec{1}\\
*  \vec{1}^{T} & 0\end{array}\right)\left(\begin{array}{c}
*  \vec{x}\\
*  y\end{array}\right)=\left(\begin{array}{c}
*  \mathbf{J}\vec{x}+y\\
*  \vec{1}\cdot\vec{x}\end{array}\right)
*  \f]
*  where \f$J\f$ is the classical Coulomb matrix.
*  To speed this up, employ prescreening
*/
void MultiplyByA(const Molecule &Mol, double *X, int N, double *V, double CoulIntMaxR)
{
	double R, MatrixElement;
	for (unsigned i=0; i<N; i++) V[i] = 0.;
	
	//   (A 1) (x) = (Ax + y)
	//  (1 0) (y) = ( 1.x  )
	for (int i=0; i<N-1; i++)
	{
		for (int j=0; j<N-1; j++)
		{
			if (i == j) MatrixElement = Mol.Atoms[i].Element.Hardness;
			else
			{
				R = Distance(Mol.Atoms[i].Basis.Position, Mol.Atoms[j].Basis.Position);
				if (R > CoulIntMaxR) MatrixElement = 1. / R;
				MatrixElement = sGTOCoulInt(Mol.Atoms[i].Basis.zeta, Mol.Atoms[j].Basis.zeta, R);
			}
			V[i] += MatrixElement * X[j];
		}
		V[i] += X[N-1];
	}
	for (int i=0; i<N-1; i++) V[N-1] += X[i];
}

inline void dsolver(const Molecule &Mol, double *b, int N, double *x, double CoulIntMaxR)
{
	int max_k = 10; // maximum number of iterations
	max_k = max(max_k, N+1);
	
	double tol = 1.0e-8; // Convergence tolerance
	int k,i;
	double *r = new double[N];
	double *p = new double[N];
	double *q = new double[N];
	double *z = new double[N];
	double *tmp = new double[N];
	double alpha, norm, gamma, gamma0;
	bool Verbose = true;
	
	// Termination criterion norm
	double critical_norm = tol * dnrm2(N, b, 1);
	
	// Calculate initial guess from diagonal part
	for (int i=0; i<N-1; i++)
		x[i] = b[i] / Mol.Atoms[i].Element.Hardness;
	x[N-1] = ONE;
	
	// Calculate residual r = b - Ax
	MultiplyByA(Mol, x, N, tmp, CoulIntMaxR);
	for (int i=0; i<N; i++) r[i] = b[i] - tmp[i];
	
	// Calculate norm
	norm = dnrm2(N, r, 1);
	if (Verbose)
		printf("Iteration 0 : norm = %lf %lf\n", norm, norm/critical_norm);
	
	for (int k=0; k<max_k; k++)
	{
		// Generate preconditioned P z = r
		for (int i=0; i<N-1; i++)
			z[i] = r[i] / Mol.Atoms[i].Element.Hardness;
		z[N-1] = ONE;
		
		// Propagate old vectors
		gamma0 = gamma;
		gamma  = ddot(N, r, 1, z, 1);
		
		if (k != 1) daxpy(N, gamma/gamma0, p, 1, z, 1);
		dcopy(N, z, 1, p, 1);
		
		// Form matrix-vector product q = A p
		MultiplyByA(Mol, p, N, q, CoulIntMaxR);

		// Calculate step size
			alpha = gamma / ddot(N, p, 1, q, 1);
		
		// Propagate by step size x = x + alpha * p
		daxpy(N, alpha, p, 1, x, 1);
		// r = r - alpha * q
		daxpy(N, -alpha, q, 1, r, 1);
		
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
	if (Verbose)
		printf("dcg: solution found with residual %lf\n", norm);
		
	delete [] r;
	delete [] p;
	delete [] q;
	delete [] z;
	delete [] tmp;
}

/**
* @brief Populates atomic charges according to the QTPIE charge model
* @param Mol : of the Molecule data type Mol.Atoms(i).Charge are computed
* @note The model is described in the paper below:
*   J. Chen and T. J. Martinez, Chem. Phys. Lett., 438 (4-6), 2007, 315-320
*/
inline void dQTPIE(Molecule &Mol)
{
	double ThisCharge; // Temporary atomic charge variable
	double Overlap; // Temporary overlap integral
	double R; // Temporary distance vairiable
	
	int i1, i2; // i1-i2 loop over atoms
	
	// Define size of matrix problem
	int N = Mol.NumAtoms;
	double *Voltage = new double [N+1];
	
	// Calculate integral pre-screening threshold
	double SmallestGaussianExponentInSystem = 1.0e40;
	for (int i1=0; i1<Mol.NumAtoms; i1++)
		SmallestGaussianExponentInSystem = min(
			SmallestGaussianExponentInSystem,
			Mol.Atoms[i1].Basis.zeta);
	
	/// Store pre-calculated thresholds for prescreening
	double OvIntMaxR = sqrt(
		log( (pi/pow(2*SmallestGaussianExponentInSystem,3))
					/ pow(OvIntThreshold,2))
			/SmallestGaussianExponentInSystem);
	
	/// Store pre-calculated thresholds for prescreening
	double CoulIntMaxR = 4.32*2/sqrt(SmallestGaussianExponentInSystem);
		// Hard coded to threshold of 1d-9
	
	// Construct voltages
//c$omp   parallel do private(R, Overlap, ThisCharge) schedule(static, 32)
	for (int i1=0; i1<N; i1++)
	{
		ThisCharge = ZERO;
		for (int i2=0; i2<N; i2++)
		{
			if (i1 != i2)
			{
				R = Distance(Mol.Atoms[i1].Basis.Position, Mol.Atoms[i2].Basis.Position);
				if (R < OvIntMaxR)
				{
					Overlap = sGTOOvInt(Mol.Atoms[i1].Basis.zeta, Mol.Atoms[i2].Basis.zeta, R);
					ThisCharge -= Overlap
						* (Mol.Atoms[i1].Element.Electronegativity - Mol.Atoms[i2].Element.Electronegativity);
				}
			}
		}
		Voltage[i1] = ThisCharge / N;
	}
//c$omp   end parallel do
	
	// Put in charge constraint
	Voltage[N+1] = Mol.TotalCharge;
	
	// Use internal conjugate gradients routine
	dsolver(Mol, Voltage, N+1, 0, CoulIntMaxR);
	
	// Calculate energy
	Mol.Energy = 0.;
	for (int i1=0; i1<N; i1++)
	{
		ThisCharge = Mol.Atoms[i1].Charge;
		Mol.Energy = Mol.Energy - ThisCharge * Voltage[i1];
		for (int i2=0; i2<i1; i2++)
			Mol.Energy += ThisCharge * Mol.Atoms[i2].Charge * _E(N+1, Mol.Coulomb, i1, i2);
		Mol.Energy += 0.5 * ThisCharge * ThisCharge * _E(N+1, Mol.Coulomb, i1, i1);
	}
	
	delete [] Voltage;
}

#endif
