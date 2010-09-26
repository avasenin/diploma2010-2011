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
#ifndef QTPIE__H
#define QTPIE__H

#include "parameters.h"

/**
* @brief Populates integral matrices in Mol data type
*   Mol.Coulomb and Mol.Overlap are initialized
* @param Mol : of the Molecule data type
*/
inline void DosGTOIntegrals(Molecule &Mol)
{
	double R; // Temporary distance
	double Integral, Norm; // Temporary integrals
	Point Pos;
	
	int N = Mol.NumAtoms;
	// Calculate integral pre-screening thresholds
	double SmallestGaussianExponentInSystem = 1e40;
	for (int i=0; i<N; i++)
		SmallestGaussianExponentInSystem = min(SmallestGaussianExponentInSystem, Mol.Atoms[i].Basis.zeta);
	double OvIntMaxR = sqrt(
		log( (pi/cube(2*SmallestGaussianExponentInSystem))
				/ sqr(OvIntThreshold))
		/SmallestGaussianExponentInSystem);


	// An asymptotic expansion of erfc-1(x) gives this formula
	double CoulIntMaxR = 2 * sqrt(-log(CoulIntThreshold) / SmallestGaussianExponentInSystem);
	
	// Populate integral matrices
	//  Note: Only the lower (i2<i1) triangle of Mol.Overlap is populated
	//   because the rest of the code is written to take advantage of
	//   its symmetry.
	
//c$omp   parallel do private(R, Integral, Pos) schedule(dynamic, 32)
	for (int i1=0; i1<N; i1++)
	{
		Pos  = Mol.Atoms[i1].Basis.Position;
		for (int i2=0; i2<i1; i2++)
		{
			R = Distance(Pos, Mol.Atoms[i2].Basis.Position);
			if (R > CoulIntMaxR) { Integral = 1. / R; }
			else Integral = sGTOCoulInt(Mol.Atoms[i1].Basis.zeta, Mol.Atoms[i2].Basis.zeta, R);
			_E(N+1, Mol.Coulomb, i1, i2) = Integral;
			_E(N+1, Mol.Coulomb, i2, i1) = Integral;
			if (R > OvIntMaxR) Integral = 0.;
			else Integral = sGTOOvInt(Mol.Atoms[i1].Basis.zeta, Mol.Atoms[i2].Basis.zeta, R);
			_E(N, Mol.Overlap, i1, i2) = Integral;
		}
		// For the diagonal elements, use hardness
		_E(N+1, Mol.Coulomb, i1, i1) = Mol.Atoms[i1].Element.Hardness;
	}
//c$omp   end parallel do

	// Calculate due normalization
//c$omp   parallel do private(Norm) schedule(static, 32)
	for (int i1=0; i1<N; i1++)
	{
		Norm = ONE;
		for (int i2=0; i2<i1; i2++) Norm += _E(N, Mol.Overlap, i1, i2);
		for (int i2=i1+1; i2<N; i2++) Norm += _E(N, Mol.Overlap, i2, i1);
		Mol.OvNorm[i1] = N / Norm;
	}
//c$omp   end parallel do
}

/**
* @brief Populates integral matrices in Mol data type
*   Mol.Coulomb and Mol.Overlap are initialized
* @param Mol : of the Molecule data type
*/
inline void DosSTOIntegrals(Molecule &Mol)
{
	double *RefOverlap = new double [Mol.NumAtoms * Mol.NumAtoms];
	
	// Now compute Coulomb matrix
	int N = Mol.NumAtoms;
	for (int i1=0; i1<N; i1++)
	{
		for (int i2=1; i2<i1; i2++)
		{
			_E(N+1, Mol.Coulomb, i1, i2) = sSTOCoulInt(
				Mol.Atoms[i1].Basis.zeta, Mol.Atoms[i2].Basis.zeta,
				Mol.Atoms[i1].Basis.n , Mol.Atoms[i2].Basis.n,
				Distance(Mol.Atoms[i1].Basis.Position, Mol.Atoms[i2].Basis.Position));
			_E(N+1, Mol.Coulomb, i2, i1) = _E(N+1, Mol.Coulomb, i1, i2);
		}
		// For the diagonal elements, use hardness
		_E(N+1, Mol.Coulomb, i1, i1) = Mol.Atoms[i1].Element.Hardness;
	}
	
	// Now compute Overlap and RefOverlap matrices
	for (int i1=0; i1<Mol.NumAtoms; i1++)
	{
		for (int i2=0; i2<i1; i2++)
		{
			_E(N, Mol.Overlap, i1, i2) = sSTOOvInt(
					Mol.Atoms[i1].Basis.zeta, Mol.Atoms[i2].Basis.zeta,
					Mol.Atoms[i1].Basis.n , Mol.Atoms[i2].Basis.n,
					Distance(Mol.Atoms[i1].Basis.Position, Mol.Atoms[i2].Basis.Position));
			
			// Calculate the same quantity but referenced to an intrinsic length scale
			_E(N, RefOverlap, i1, i2) = sSTOOvInt(
					Mol.Atoms[i1].Basis.zeta, Mol.Atoms[i2].Basis.zeta,
					Mol.Atoms[i1].Basis.n , Mol.Atoms[i2].Basis.n,
					ExpectR(Mol.Atoms[i1].Basis) + ExpectR(Mol.Atoms[i2].Basis));
			
			// Fill in the other triangle
			_E(N, Mol.Overlap, i2, i1) = _E(N, Mol.Overlap, i1, i2);
			_E(N, RefOverlap, i2, i1) = _E(N, RefOverlap, i1, i2);
		}
		// For the diagonal elements, the overlap is just the orbital normalization
		_E(N, Mol.Overlap, i1, i1) = 1.;
		_E(N, RefOverlap, i1, i1) = 1.;
	}
	
	// Now compute normalization of Attenuation (overlap) matrix
	for (int i1=0; i1<Mol.NumAtoms; i1++)
	{
		Mol.OvNorm[i1] = 0.;
		for (int i2=1; i2<Mol.NumAtoms; i2++)
			Mol.OvNorm[i1] += _E(N, RefOverlap, i1, i2);
		Mol.OvNorm[i1] /= Mol.NumAtoms;
	}
	
	// Multiply in Norm
	for (int i1=1; i1<Mol.NumAtoms; i1++)
	{
		for (int i2=1; i2<Mol.NumAtoms; i2++)
			_E(N, Mol.Overlap, i1, i2) /= Mol.OvNorm[i1];
	}
	
	delete [] RefOverlap;
}

/**
* @brief Populates atomic charges according to the QEq(-H) charge model
* @param Mol : of the Molecule data type Mol.Atoms(i)%Charge are computed
* @note The model is described in the seminal paper below:
*  "Charge equilibration for Molecular dynamics simulations"
*  A. K. Rappe and W. A. Goddard, J. Phys. Chem., 1991, 95(8), 3358-3363
* @note This implementation does not do the additional procedure for H atoms
*  nor does it check for overly large charges that exceed the principal
*  quantum number of the given atom.
*/
inline void QEq(Molecule &Mol)
{
	int N = Mol.NumAtoms;
	double *Capacitance = new double[(N + 1) * (N + 1)];
	double *Voltage = new double[N + 1];
	double *Charge  = new double[N + 1];
	for (int i=0; i<N+1; i++) Voltage[i] = Charge[i] = 0.;
	
	// Calculate problem
	// Construct the capacitance matrix
	for (int i=0; i<N; i++)
	for (int j=0; j<N; j++)
		_E(N+1, Capacitance, i, j) = _E(N+1, Mol.Coulomb, i, j);
	// Add charge conservation constraints
	for (int i=0; i<N; i++) 
	{
		_E(N+1, Capacitance, i, N) = ONE;
		_E(N+1, Capacitance, N, i) = ONE;
	}
	_E(N+1, Capacitance, N, N) = ZERO;
	
	// Construct voltage matrix
	for (int i=0; i<N; i++) Voltage[i] = -Mol.Atoms[i].Element.Electronegativity;
	
	// The last row expresses charge conservation
	Voltage[N] = Mol.TotalCharge;
	
	// Solve for Charge, where Capacitance * Charge = Voltage
	// i.e. an A*x = b problem where A is symmetric and positive nonnegative
	solver(N + 1, Capacitance, Voltage, Charge);
	for (int i=0; i<N; i++)
	{
	Mol.Atoms[i].Charge = Charge[i];
	printf("i=%d CHARGE EQ=%10.5e\n", i, Charge[i]);
	}
	
	// Calculate energy
	Mol.Energy = 0.;
	for (int i1=0; i1<N; i1++)
	{
		Mol.Energy += Mol.Atoms[i1].Charge * Mol.Atoms[i1].Element.Electronegativity;
		for (int i2=0; i2<N; i2++)
			Mol.Energy += 0.5 * Mol.Atoms[i1].Charge * Mol.Atoms[i2].Charge * _E(N+1, Mol.Coulomb, i1, i2);
	}
	printf("EQ Mol.Energy=%10.5e\n", Mol.Energy);
	
	delete [] Capacitance;
	delete [] Voltage;
	delete [] Charge;
}

/**
* @brief  Populates atomic charges according to the QTPIE charge model
* @param Mol : of the Molecule data type
*   Mol.Atoms(i)%Charge are computed
* @note The model is described in the paper below:
*    J. Chen and T. J. Martinez, Chem. Phys. Lett., 438 (4-6), 2007, 315-320
*/
inline void QTPIE(Molecule &Mol)
{
	int N = Mol.NumAtoms;
	double ThisCharge; // Temporary atomic charge variable
	double Overlap;    // Temporary integral
	double VoltageDifference, Norm;
	
	double *Voltage = new double [N+1];
	double *Charge  = new double [N+1];
	
	// Construct problem
	// Add charge conservation constraints
	for (int i=0; i<N; i++)
	{
		_E(N+1, Mol.Coulomb, i, N) = ONE;
		_E(N+1, Mol.Coulomb, N, i) = ONE;
		Charge[i] = 0.;
	}
	_E(N+1, Mol.Coulomb, N, N) = ZERO;
	Charge[N] = 0.;
	
// Construct voltages
// The code here is a little convoluted but knows that Overlap is sparse and symmetric

//c$omp   parallel do private(ThisCharge, Overlap, Norm,
//c$omp&     VoltageDifference][i2] schedule(dynamic, 32)
	for (int i1=0; i1<N; i1++)
	{
		ThisCharge = ZERO;
		Norm = Mol.OvNorm[i1];
		for (int i2=0; i2<i1; i2++)
		{
			Overlap = _E(N, Mol.Overlap, i1, i2);
			if (Overlap != ZERO)
			{
				VoltageDifference = (Mol.Atoms[i1].Element.Electronegativity - Mol.Atoms[i2].Element.Electronegativity);
				if (VoltageDifference != ZERO)
					ThisCharge -= Overlap * VoltageDifference * Norm;
			}
		}
		for (int i2=i1+1; i2<N; i2++)
		{
			Overlap = _E(N, Mol.Overlap, i2, i1);
			if (Overlap != ZERO)
			{
				VoltageDifference = Mol.Atoms[i1].Element.Electronegativity - Mol.Atoms[i2].Element.Electronegativity;
				if (VoltageDifference != ZERO)
					ThisCharge -= Overlap * VoltageDifference * Norm;
			}
		}
		Voltage[i1] = ThisCharge / N;
		//printf("i1=%d volt=%lf\n", i1+1, Voltage[i1]);
	}
//c$omp   end parallel do
	
	// Put in charge constraint
	Voltage[N] = Mol.TotalCharge;
	//printf("%d volt=%lf\n", N+1, Voltage[N]);
	
	// Use internal conjugate gradients routine
	solver(N + 1, Mol.Coulomb, Voltage, Charge);
	
	// Copy out solution from work into charges
	for (int i=0; i<N; i++)
	{
	Mol.Atoms[i].Charge = Charge[i];
printf("i=%d Charge=%lf\n", i, Charge[i]);
	}

	
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
	printf("Mol.Energy=%lf\n", Mol.Energy);
	
	delete [] Voltage;
	delete [] Charge;
}

/**
*  Computes energy gradients numerically
*  Calculates energy gradients using the method of finite differences
*  using forward gradients
*  As you can imagine, this is pretty slow
*  You should not use this routine!
* @param Mol : of the Molecule data type
*  Mol.EGradient is calculated
*/
inline void DoGradientsByFiniteDifference(Molecule &Mol)
{
	double OriginalEnergy = Mol.Energy; // Save current energy
	double Eps = 1.0e-3;
	
	// Calculate energy gradients
	int N = Mol.NumAtoms;
	for (int i1=0; i1<Mol.NumAtoms; i1++)
	{
		for (int i2=0; i2<3; i2++)
		{
			// Perturb Geometry
			Mol.Atoms[i1].Position[i2] += Eps;
			Mol.Atoms[i1].Basis.Position[i2] += Eps;
			// Redo QTPIE
			DosGTOIntegrals(Mol);
			QTPIE(Mol);
			// Calculate gradient
			_E(3, (double*)Mol.EGradient, i1, i2) = (Mol.Energy - OriginalEnergy) / Eps;
			// Perturb Geometry
			Mol.Atoms[i1].Position[i2] -= Eps;
			Mol.Atoms[i1].Basis.Position[i2] -= Eps;
		}
	}
	// Redo integrals
	DosGTOIntegrals(Mol);
}

/**
*  Computes energy gradients analytically
*  Calculates energy gradients using analytic derivatives
* @param Mol : of the Molecule data type
*  Mol.EGradient is calculated
*/
inline void DoGradientsAnalytically(Molecule &Mol)
{
	double a, b, R, Force;
	int m, n;
	
	// Initialize gradients
	for (int i=0; i<Mol.NumAtoms*3; i++)
		((double*)Mol.EGradient)[i] = 0.;
	
	for (int i1=0; i1<Mol.NumAtoms; i1++)
	{
		// Obtain basis set parameters for atom i1
		a = Mol.Atoms[i1].Basis.zeta;
		m = Mol.Atoms[i1].Basis.n;
		
		// Calculate contribution to gradient from voltage term
		int N = Mol.NumAtoms;
		for (int i2=0; i2<Mol.NumAtoms; i2++)
		{
			// Diagonal part has no contribution to gradient
			if (i1 != i2)
			{
				// Obtain basis set parameters for atom i2
				b = Mol.Atoms[i2].Basis.zeta;
				n = Mol.Atoms[i2].Basis.n;
						// Calculate pairwise distance
				R = Distance(Mol.Atoms[i1].Basis.Position, Mol.Atoms[i2].Basis.Position);
				Force = 2 * Mol.Atoms[i1].Charge / Mol.NumAtoms
					* (Mol.Atoms[i1].Element.Electronegativity - Mol.Atoms[i2].Element.Electronegativity)
					* sGTOOvIntGrad(a,b,R) / Mol.OvNorm[i1];
				
				Force -= Mol.Atoms[i1].Charge / Mol.NumAtoms
					* (Mol.Atoms[i1].Element.Electronegativity - Mol.Atoms[i2].Element.Electronegativity)
					* _E(N, Mol.Overlap, i1, i2) / sqr(Mol.OvNorm[i1])
					* sGTOOvIntGrad(a,b,R);
					
				Force += Mol.Atoms[i1].Charge * Mol.Atoms[i2].Charge * sGTOCoulIntGrad(a,b,R);
				// Calculates projection onto direction vector
				//  $Temp*\frac{\partial R_{i1,i2}}{\partial R_{k,i3}}
				//  * (\delta_{i1,k} - \delta_{i2,k})$
				Force /= R;
				for (int i3=0; i3<3; i3++)
					_E(3, (double*)Mol.EGradient, i1, i3) += (Mol.Atoms[i1].Basis.Position[i3] 
						- Mol.Atoms[i2].Basis.Position[i3]) * Force;
			}
		}
	}
}

#endif
