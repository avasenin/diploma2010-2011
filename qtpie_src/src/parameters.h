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
#ifndef PARAMETERS__H
#define PARAMETERS__H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomicunits.h"

#undef max
#undef min
extern "C" {
//#include "clapack.h"
#include "cblas.h"
}

double dnrm2(const int N, const double *X, const int incX)
{ return cblas_dnrm2(N, X, incX); }

void dcopy(unsigned n, const double *x, unsigned incx, double *y, unsigned incy)
{ cblas_dcopy(n, x, incx, y, incy); }

void dgemv(const enum CBLAS_TRANSPOSE trans, unsigned m, unsigned n,
			double alpha, double *a, unsigned lda, const double *x, unsigned incx,
			double beta, double *y, unsigned incy)
{ cblas_dgemv(CblasRowMajor, trans, m, n, alpha, a, lda, x, incx, beta, y, incy); }

double ddot(unsigned n, const double *x, unsigned incx, const double *y, unsigned incy)
{ return cblas_ddot(n, x, incx, y, incy); }

void daxpy(unsigned n, double alpha, const double *x, unsigned incx, double *y, unsigned incy)
{ return cblas_daxpy(n, alpha, x, incx, y, incy); }


#define _E(N, A, i, j)    (*(A + (i)*(N) + (j)))
#define _D(N, A, i)       (*(A + (i)*((N) + 1)))

template <typename T> inline T sqr(T t) { return t*t; }
template <typename T> inline T cube(T t) { return t*t*t; }
template <typename T> inline T max(T t1, T t2) { return t1>t2? t1 : t2; }
template <typename T> inline T min(T t1, T t2) { return t1>t2? t2 : t1; }
template <typename T> inline T mod(T x, T m)
{ return x<0 ? m - 1 - ((-x) - 1)%m : x%m; }

struct Point
{
	double &operator[](unsigned i) { return _[i]; }
	double operator[](unsigned i) const { return _[i]; }
	Point &operator=(const Point &p)
	{
		_[0] = p[0];
		_[1] = p[1];
		_[2] = p[2];
		return *this;
	}
private:
	double _[3];
};

/**
* @brief Computes pairwise distances from Cartesian coordinates
* @param Point1,Point2: 3-vectors of double precisions
* @return Cartesian distance in atomic units
*/
inline double Distance(const Point &Point1, const Point &Point2)
{
	double x, y, z;
	x = Point2[0] - Point1[0];
	y = Point2[1] - Point1[1];
	z = Point2[2] - Point1[2];
	return sqrt(x*x + y*y + z*z);
}


/**
* @brief Stores parameters for our charge models
*/

const double pi =  3.141592653589793;


/**
* @brief Parameters for a s-type Slater type orbital (STO) basis function
* @param n : principal quantum number
* @param zeta : zeta exponent with dimensions of inverse length in atomic units
*/
struct sSTO
{
	Point Position;
	int n;
	double zeta;
};

/**
* @brief Parameters for a s-type Gaussian type orbital (GTO) basis function
* @param Position : an array of three double precisions describing Cartesian coordinates
* @param zeta: exponent with dimensions of inverse square length in atomic units
*/
struct sGTO
{
	Point Position;
	double zeta;
};

/**
* @brief Atomic parameters
* @param Symbol            : elemental symbol
* @param Z                 : atomic number
* @param FormalCharge      : formal charge, integers only
* @param Electronegativity : Mulliken electronegativity in atomic units
* @param Hardness          : Parr-Pearson chemical hardness in atomic units
*/
struct AtomData
{
	char Symbol[3];
	int Z, FormalCharge;
	double Electronegativity, Hardness;
};

//-----------------------------------------------
// Here are a bunch of predefined elements
// As parameterized by Rappe and Goddard for QEq
//-----------------------------------------------
AtomData Hydrogen   = {  "H",  1, 0, 4.528*eV, 13.890*eV };
AtomData Lithium    = { "Li",  3, 0, 3.006*eV,  4.772*eV };
AtomData Carbon     = {  "C",  6, 0, 5.343*eV, 10.126*eV };
AtomData Nitrogen   = {  "N",  7, 0, 7.139*eV, 12.844*eV };
AtomData Oxygen     = {  "O",  8, 0, 8.741*eV, 13.364*eV };
AtomData Fluorine   = {  "F",  9, 0,10.874*eV, 14.948*eV };
AtomData Sodium     = { "Na", 11, 0, 2.843*eV,  4.592*eV };
AtomData Silicon    = { "Si", 14, 0, 4.168*eV,  6.974*eV };
AtomData Phosphorus = {  "P", 15, 0, 5.463*eV,  8.000*eV };
AtomData Sulphur    = {  "S", 16, 0, 6.084*eV, 10.660*eV };
AtomData Chlorine   = { "Cl", 17, 0, 8.564*eV,  9.892*eV };
AtomData Potassium  = {  "K", 19, 0, 2.421*eV,  3.840*eV };
AtomData Bromine    = { "Br", 35, 0, 7.790*eV,  8.850*eV };
AtomData Rubidium   = { "Rb", 37, 0, 2.331*eV,  3.692*eV };
AtomData Iodine     = {  "I", 53, 0, 6.822*eV,  7.524*eV };
AtomData Cesium     = { "Cs", 55, 0, 2.183*eV,  3.422*eV };

const int numParameterizedAtoms = 16; // Number of defined atomic parameters
AtomData ParameterizedAtoms[numParameterizedAtoms] =
{
	Hydrogen, Lithium, Carbon, Nitrogen, Oxygen, Fluorine,
	Sodium, Silicon, Phosphorus, Sulphur, Chlorine,
	Potassium, Bromine, Rubidium, Iodine, Cesium
}; // Array of defined atomic parameters

// Parameters for cations. All experimental values!
/** @brief Sodium cation */
AtomData SodiumCation = { "Na", 11, +1, 4562*kJ_mol, 5.13908*eV };

/**
* @brief Describes an atom in a molecule
* @param  Element : type(AtomData) containing atomic parameters
* @param Basis   : A basis function associated with the atom
* @param Position: double precision(3) vector of Cartesian coordinates describing spatial location
* @param Charge  : double precision, result of charge model calculation
*/
struct Atom
{
	AtomData Element;
	sSTO Basis;
	Point Position;
	double Charge;
};

// Data for newly parameterized Gaussian basis set
double GaussianExponent[numParameterizedAtoms] =
{
	0.534337523756312, 0.166838519142176, 0.206883838259186,
	0.221439796025873, 0.223967308625516, 0.231257590182828,
	0.095892938712585, 0.105219608142377, 0.108476721661715,
	0.115618357843499, 0.113714050615107, 0.060223294377778,
	0.070087547802259, 0.041999054745368, 0.068562697575073,
	0.030719481189777
};

/// Threshold for calculating overlap integrals
const double OvIntThreshold = 1.0e-9;

/// Threshold for calculating Coulomb integrals
const double CoulIntThreshold = 1.0e-9;

/**
* @brief Computes the expectation value of the radial distance over s-type STOs
* @param Basis: s-type STO basis function
* @return the expectation value of the radial distance over s-type STOs
*/
inline double ExpectR(const sSTO &Basis)
{ return (Basis.n + 0.5) / Basis.zeta; }

/**
* @brief Computes the principal quantum number of an atom given its atomic number
* @param theAtom Atom to determine principle quantum number for
* @return the principal quantum number
*/
inline int pqn(const AtomData &theAtom)
{
	static int maxelectrons[7] = { 2, 10, 18, 36, 54, 86, 118 };
		// Lookup table for max number of electrons for that quantum number
	
	int pqn_ = 1;
	// work through each shell
	for (int j=0; j<7; j++) if (theAtom.Z > maxelectrons[j]) pqn_++;
	return pqn_;
}

#endif
