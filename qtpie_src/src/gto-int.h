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
#ifndef GTO_INT__H
#define GTO_INT__H

#include "parameters.h"
#include "sto-int.h"

/**
* @brief Assigns a Gaussian-type orbital to the atom
* @note See research notes dated 2008-03-14
*/
inline void AssignsGTOBasis(Atom &theAtom)
{
	// Assign position
	theAtom.Basis.Position = theAtom.Position;
	
	// Assign Gaussian orbital exponent
	for (int i=0; i<numParameterizedAtoms; i++)
	{
		if (theAtom.Element.Z == ParameterizedAtoms[i].Z)
			theAtom.Basis.zeta = GaussianExponent[i];
	}
	// Scaling - does not work well with QEq parameters
}

/**
*  Calculates best-fit GTO exponent given best-fit STO exponent
* @param n: principal quantum number
* @param zeta: exponent for s-type Slater orbital
* @return the best-fit exponent for the s-type Gaussian orbital
* @note See research notes dated 2007-08-31
* @note deprecated
*/
double sSTO2sGTO(int n, double zeta)
{
	static double conversion[7] =
	{ 
		0.2709498089, 0.2527430925, 0.2097635701,
		0.1760307725, 0.1507985107, 0.1315902101,
		0.1165917484 
	};
	return conversion[n] * zeta * zeta;
}

/**
*  Calculates a best-fit Gaussian-type orbital (STO-1G) to
*  the Slater-type orbital defined from the hardness parameters
* @param Hardness: chemical hardness in atomic units
* @param n: principal quantum number
* @note See research notes dated 2007-08-30
*  deprecated
*/
inline void AssignsSTO1GBasis(Atom &theAtom)
{
	// Approximate the exact value of the constant of proportionality
	//  by its value at a very small distance epsilon
	//  since the exact R = 0 case has not be programmed into STOIntegrals
	double epsilon = 1.0e-6;
	
	// Assign position
	theAtom.Basis.Position = theAtom.Position;
	
	// Assign principal quantum number
	int n = pqn(theAtom.Element);
	theAtom.Basis.n = n;
	
	// Assign orbital exponent
	double zeta = pow(sSTOCoulInt(1., 1., n, n, epsilon)
		/ theAtom.Element.Hardness, -1./(3. + 2.*n));
	
	// Rewrite it with best-fit Gaussian
	theAtom.Basis.zeta = sSTO2sGTO(n, zeta);
}

/**
*  Computes Coulomb integral analytically over s-type GTOs
* 
*  Computes the two-center Coulomb integral over Gaussian-type orbitals
*  of s symmetry.
* 
* @param a: Gaussian exponent of first atom in atomic units (inverse squared Bohr)
* @param b: Gaussian exponent of second atom in atomic units (inverse squared Bohr)
* @param R: internuclear distance in atomic units (bohr)
* @return the value of the Coulomb potential energy integral
* @note Reference: T. Helgaker, P. Jorgensen, J. Olsen, Molecular Electronic Structure Theory
*    Wiley, NY, 2000, Equations (9.7.21) and (9.8.23)
*/
inline double sGTOCoulInt(double a, double b, double R)
{
	double p = sqrt(a * b / (a + b));
	return erf(p * R) / R;
}

/** Computes overlap integral analytically over s-type GTOs
*  Computes the overlap integral over two Gaussian-type orbitals of s symmetry.
* @param a: Gaussian exponent of first atom in atomic units (inverse squared Bohr)
* @param b: Gaussian exponent of second atom in atomic units (inverse squared Bohr)
* @param R: internuclear distance in atomic units (bohr)
* @note Reference: T. Helgaker, P. Jorgensen, J. Olsen, Molecular Electronic Structure Theory
*    Wiley, NY, 2000, Equation (9.2.41)
* @note With normalization constants added, calculates
*  \f[
*  S = \left(\frac{4\alpha\beta}{(\alpha + \beta)^2}\right)^\frac{3}{4}
*      \exp\left(-\frac{\alpha\beta}{\alpha+\beta} R^2 \right)
*  \f]
*/
inline double sGTOOvInt(double a, double b, double R)
{
	double p, q;
	p = a + b;
	q = a * b / p;
	return pow(4*q/p, 0.75) * exp(-q*R*R);
}

/**
* @brief Computes derivative of Coulomb integral wrt R
* @param a: Gaussian exponent of first atom in atomic units (inverse squared Bohr)
* @param b: Gaussian exponent of second atom in atomic units (inverse squared Bohr)
* @param R: internuclear distance in atomic units (bohr)
* @return the derivative of the Coulomb potential energy integral
*/
inline double sGTOCoulIntGrad(double a, double b, double R)
{
	double pi =  3.141592653589793;
	
	if (fabs(R) == 0.)
	{
		printf("FATAL ERROR: R = 0 in sGTOCoulIntGrad\n");
		exit(0);
	}
	
	double p = sqrt(a * b / (a + b));
	return 2. * p / (R * sqrt(pi)) * exp(-sqr(p*R)) - sGTOCoulInt(a,b,R) / R;
}

/**
*  Computes gradient of overlap integral wrt R
*  Computes the derivative of the overlap integral over two Gaussian-type orbitals of s symmetry.
* @param a: Gaussian exponent of first atom in atomic units (inverse squared Bohr)
* @param b: Gaussian exponent of second atom in atomic units (inverse squared Bohr)
* @param R: internuclear distance in atomic units (bohr)
* @return the derivative of the sGTOOvInt integral
*/
inline double sGTOOvIntGrad(double a, double b, double R)
{
	return -2 * (a*b)/(a+b)* R * sGTOOvInt(a, b, R);
}

#endif
