/***************************************************************************
 *   Copyright (C) 2008 by Jiahao Chen                                     *
 *   cjiahao@gmail.com                                                     *
 *   Copyright (C) 2008 by Eduard S. Fomin                                 *
 *   fomin@bionet.nsc.ru                                                   *
 *                                                                         *
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
#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <cstdlib>
#include <map>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "_minimize_lbfgs.h"
using namespace prgkern;

#define _E(N, A, i, j)    (*(A + (i)*(N) + (j)))
class v3dense
{
	double v_[3];
public:
	v3dense() { v_[0] = 0; v_[1] = 0; v_[2] = 0; }
	v3dense(double x, double y, double z) { v_[0] = x; v_[1] = y; v_[2] = z; }
	double &operator[](unsigned i) { return v_[i]; }
	double operator[](unsigned i) const { return v_[i]; }
};

inline double distance1(v3dense &X, v3dense &Y)
{ return sqrt(sqr(X[0] - Y[0]) + sqr(X[1] - Y[1]) + sqr(X[2] - Y[2])); }

const double EV_2_HARTREE          = 3.67493245e-2;      // electron volt to Hartree
const double ANGSTROM_2_BORH       = 1./0.529177249;     // Angstrom to bohr

/// Threshold for calculating Coulomb integrals
const double COULOMB_INTEGRAL_THRESHOLD = 1.0e-9;
const double LOG_COULOMB_INTEGRAL_THRESHOLD = log(COULOMB_INTEGRAL_THRESHOLD);

/// Threshold for calculating overlap integrals
const double OVERLAP_INTEGRAL_THRESHOLD = 1.0e-9;
const double SQR_OVERLAP_INTEGRAL_THRESHOLD = sqr(OVERLAP_INTEGRAL_THRESHOLD);

/**
* @brief Parameters for a s-type Gaussian type orbital (GTO) basis function
* @param zeta exponent with dimensions of inverse square length in atomic units
*/
struct sGTO
{
	sGTO(double exponent=0.) : zeta(exponent) {}
	double zeta;
};

//template <typename> class Basis_;

/*template <>*/ class Basis_/*<sGTO>*/ : public std::map<std::string, sGTO>
{
	typedef std::map<std::string, sGTO> _Base;

	double Coulomn_Integral_Distance_Threshold;
	double Overlap_Integral_Distance_Threshold;
	double smallest_gaussian_exponent;

public:
	typedef std::map<std::string, sGTO>::value_type  value_type;

	Basis_() : smallest_gaussian_exponent(1.e20) {}

	void insert(const value_type &val)
	{
		_Base::insert(val);
		double zeta = val.second.zeta;
		if (smallest_gaussian_exponent > zeta)
		{
			smallest_gaussian_exponent = zeta;
			Coulomn_Integral_Distance_Threshold =
				2 * sqrt(-LOG_COULOMB_INTEGRAL_THRESHOLD / smallest_gaussian_exponent);
			Overlap_Integral_Distance_Threshold =
				sqrt(
					log( (M_PI / cube(2 * smallest_gaussian_exponent))
							/ SQR_OVERLAP_INTEGRAL_THRESHOLD
					) / smallest_gaussian_exponent
				);
		}
	}

	/**
	* @brief Computes the 2-center Coulomb integral over Gaussian-type s-orbitals analytically.
	* @note Reference: T. Helgaker, P. Jorgensen, J. Olsen, Molecular Electronic Structure Theory
	*   Wiley, NY, 2000, Equations (9.7.21) and (9.8.23)
	* @param a,b atom names
	* @param R internuclear distance in atomic units (bohr)
	* @return the value of the Coulomb potential energy integral
	*/
	double coulomb_integral(sGTO a, sGTO b, double R)
	{
		if (R > Coulomn_Integral_Distance_Threshold) return 1. / R;
		double p = sqrt(a.zeta * b.zeta / (a.zeta + b.zeta));
		return erf(p * R) / R;
	}

	/**
	* @brief Computes overlap integral analytically over s-type GTOs
	* @note Reference: T. Helgaker, P. Jorgensen, J. Olsen, Molecular Electronic Structure Theory
	*    Wiley, NY, 2000, Equation (9.2.41)
	* @param a Gaussian exponent of first atom in atomic units (inverse squared Bohr)
	* @param b Gaussian exponent of second atom in atomic units (inverse squared Bohr)
	* @param R internuclear distance in atomic units (bohr)
	*/
	double overlap_integral(sGTO a, sGTO b, double R)
	{
		if (R > Overlap_Integral_Distance_Threshold) return 0.;
		double p = a.zeta + b.zeta;
		double q = a.zeta * b.zeta / p;
		return pow(4 * q / p, 0.75) * exp(-q * R * R);
	}
};

/**
* @brief Atomic parameters
* @param Symbol elemental symbol
* @param nuclear_charge atomic number
* @param formal_charge formal charge
* @param electronegativity Mulliken electronegativity in atomic units
* @param hardness Parr-Pearson chemical hardness in atomic units
*/
struct AtomicParams
{
	std::string symbol;
	int nuclear_charge;
	int formal_charge;
	double electronegativity;
	double hardness;
	double gaussian_exponent;
};

/**
* @brief Predefined Rappe and Goddard elements && gaussians for QEq
*/
const AtomicParams RAPPLE_GODDARD_PARAMS[] =
{
	{ "H",  1, 0,  4.528 * EV_2_HARTREE, 13.890 * EV_2_HARTREE, 0.534337523756312 },
	{ "LI", 3, 0,  3.006 * EV_2_HARTREE,  4.772 * EV_2_HARTREE, 0.166838519142176 },
	{ "C",  6, 0,  5.343 * EV_2_HARTREE, 10.126 * EV_2_HARTREE, 0.206883838259186 },
	{ "N",  7, 0,  7.139 * EV_2_HARTREE, 12.844 * EV_2_HARTREE, 0.221439796025873 },
	{ "O",  8, 0,  8.741 * EV_2_HARTREE, 13.364 * EV_2_HARTREE, 0.223967308625516 },
	{ "F",  9, 0, 10.874 * EV_2_HARTREE, 14.948 * EV_2_HARTREE, 0.231257590182828 },
	{ "NA", 11, 0,  2.843 * EV_2_HARTREE,  4.592 * EV_2_HARTREE, 0.095892938712585 },
	{ "SI", 14, 0,  4.168 * EV_2_HARTREE,  6.974 * EV_2_HARTREE, 0.105219608142377 },
	{ "P",  15, 0,  5.463 * EV_2_HARTREE,  8.000 * EV_2_HARTREE, 0.108476721661715 },
	{ "S",  16, 0,  6.084 * EV_2_HARTREE, 10.660 * EV_2_HARTREE, 0.115618357843499 },
	{ "CL", 17, 0,  8.564 * EV_2_HARTREE,  9.892 * EV_2_HARTREE, 0.113714050615107 },
	{ "K",  19, 0,  2.421 * EV_2_HARTREE,  3.840 * EV_2_HARTREE, 0.060223294377778 },
	{ "BR", 35, 0,  7.790 * EV_2_HARTREE,  8.850 * EV_2_HARTREE, 0.070087547802259 },
	{ "RB", 37, 0,  2.331 * EV_2_HARTREE,  3.692 * EV_2_HARTREE, 0.041999054745368 },
	{ "I",  53, 0,  6.822 * EV_2_HARTREE,  7.524 * EV_2_HARTREE, 0.068562697575073 },
	{ "CS", 55, 0,  2.183 * EV_2_HARTREE,  3.422 * EV_2_HARTREE, 0.030719481189777 }
};
const unsigned RAPPLE_GODDARD_PARAM_COUNT =
	sizeof(RAPPLE_GODDARD_PARAMS) / sizeof(AtomicParams);

/**
* @brief an atom in a molecule
* @param element atomic parameters
* @param basis basis function associated with the atom
* @param X cartesian coordinates of atom (in atom units!)
*/
struct Atom
{
	std::string symbol;
	v3dense X;
	double charge;
	double du__dq;
	double electronegativity;
	double hardness;
};

/**
* @brief a molecular system
*/
class Molecule
{
	typedef Basis_/*<sGTO>*/        _Basis;
	typedef _Basis::value_type  _Orbital;
	typedef Atom                _Atom;

public:
	std::vector<_Atom> atoms_; // array of atoms
	Basis_/*<sGTO>*/ basis_; // basis set
	double total_charge_; // total charge of system
	double energy_; // electrostatic contribution to the potential energy
	int number_of_bonds;
public:
	typedef _Atom       atom_type;

	Molecule(double total_charge=0.) : total_charge_(total_charge), energy_(0.), number_of_bonds(0)
	{
	}

	void insert(std::string symbol, const v3dense &X)
	{
		// look up AtomSymbol to assign parameters
		for (unsigned i=0; i<RAPPLE_GODDARD_PARAM_COUNT; i++)
		{
			if (symbol == RAPPLE_GODDARD_PARAMS[i].symbol)
			{
				Atom atom;
				atom.symbol = symbol;
				atom.electronegativity = RAPPLE_GODDARD_PARAMS[i].electronegativity;
				atom.hardness = RAPPLE_GODDARD_PARAMS[i].hardness;
				atom.X[0] = X[0] * ANGSTROM_2_BORH;
				atom.X[1] = X[1] * ANGSTROM_2_BORH;
				atom.X[2] = X[2] * ANGSTROM_2_BORH;
				atoms_.push_back(atom);
				basis_.insert(_Orbital(symbol, RAPPLE_GODDARD_PARAMS[i].gaussian_exponent));
				return;
			}
		}
		std::string msg = _S("[ERROR] Unknown element type : ") + symbol;
		PRINT_ERR(msg);
	}

	/**
	* @brief loads an external file containing a XYZ geometry
	* @param filename name of the XYZ geometry file
	*/
	void load(const char *filename);

	/**
	* @brief populates integral matrices
	*/
	void make_coulomb_integrals(std::vector<double> &coulomb_matrix);
	std::pair<int, int> get_bond_by_index(int bondIndex);
	void make_topology_matrix();
	void make_geometric_matrix();
	double dU__dQ();

	unsigned get_atom_count() const { return atoms_.size(); }
	const _Atom* get_atoms() const { return &atoms_[0]; }
	_Atom* get_atoms() { return &atoms_[0]; }
	std::vector<double> coulomb_;
	std::vector<int> bonds;
	std::vector<int> topology; //matrix A in article
	std::vector<double> geometric; //matrix K in article
};

inline void Molecule::load(const char *filename)
{
	FILE *file = fopen(filename, "rt");
	if (file == 0)
	{
		printf("Problem loading geometry file %s\n", filename);
		exit(0);
	}

	char buffer[255], symbol[3];
	unsigned count;
	double x, y, z;

	fgets(buffer, 255, file); // atom count
	sscanf(buffer, "%d", &count);
	fgets(buffer, 255, file); // skip comment

	for (unsigned i=0; i<count; i++)
	{
		fgets(buffer, 255, file);
		sscanf(buffer, "%s %lf %lf %lf", &symbol, &x, &y, &z);
		for (unsigned j=0; j<2; j++) symbol[j] = toupper(symbol[j]);
		insert(std::string(symbol), v3dense(x, y, z));
	}
	bonds.resize(count*count);
	for (unsigned i=0; i < count; i++)
	{
		for (unsigned j=0; j < count; j++)
		{
			int bondMark = 0;
			fscanf(file, "%d", &bondMark);
			if (bondMark)
			{
				_E(count, &bonds[0], i, j) = bondMark;
				number_of_bonds++;
			}
		}
	}
	//delete symmetric bounds
	number_of_bonds /= 2;
	fclose(file);
	make_topology_matrix();
	make_geometric_matrix();
	make_coulomb_integrals(coulomb_);
}
inline void Molecule::make_topology_matrix()
{
	int numberOfAtoms = atoms_.size();
	topology.resize(number_of_bonds * numberOfAtoms);
	int bondIndex = 0;
	for (int i=0; i < numberOfAtoms; i++)
	{
		for (int j=0; j < i; j++)
		{
			if (_E(numberOfAtoms, &bonds[0], i, j))
			{
				_E(number_of_bonds, &topology[0], i, bondIndex) = 1;
				_E(number_of_bonds, &topology[0], j, bondIndex) = -1;
				bondIndex++;
			}
		}
	}
}
inline std::pair<int, int> Molecule::get_bond_by_index(int bondIndex)
{
	assert(bondIndex < number_of_bonds);
	int numberOfAtoms = atoms_.size();
	std::pair<int, int> ret = std::pair<int, int>(-1,-1);
	for (int i = 0; i < numberOfAtoms; i++)
	{
		int isIncludeToBond = _E(number_of_bonds, &topology[0], i, bondIndex);
		switch (isIncludeToBond) {
			case 1:
				ret.first = i;
				break;
			case -1:
				ret.second = i;
				break;
		}
	}
	assert(-1 != ret.first && -1 != ret.second);
	return ret;
}
inline void Molecule::make_geometric_matrix()
{
	geometric.resize(number_of_bonds*number_of_bonds);
	for (int i = 0; i < number_of_bonds; i++)
	{
		std::pair<int, int> bond = get_bond_by_index(i);
		sGTO first_sGTO= basis_[atoms_[bond.first].symbol];
		sGTO second_sGTO= basis_[atoms_[bond.second].symbol];
		double r = distance1(atoms_[bond.first].X, atoms_[bond.second].X);
		_E(number_of_bonds, &geometric[0], i, i) = basis_.overlap_integral(first_sGTO, second_sGTO, r);
	}
}
inline void Molecule::make_coulomb_integrals(std::vector<double> &coulomb_matrix)
{
	unsigned n = atoms_.size();
	if (coulomb_matrix.size() != n * n) coulomb_matrix.resize(n * n);

	for (unsigned i=0; i<n; i++)
	{
		const sGTO &sgto= basis_[atoms_[i].symbol];
		for (unsigned j=0; j<i; j++)
		{
			const sGTO &sgto__= basis_[atoms_[j].symbol];
			double r = distance1(atoms_[i].X, atoms_[j].X);
			double s = basis_.coulomb_integral(sgto, sgto__, r);
			_E(n, &coulomb_matrix[0], i, j) = s;
			_E(n, &coulomb_matrix[0], j, i) = s;
		}
		_E(n, &coulomb_matrix[0], i, i) = atoms_[i].hardness;
			// for the diagonal elements, use hardness
	}
}

inline double Molecule::dU__dQ()
{
	unsigned n = atoms_.size();

	double sum_charge = 0.;
	for (unsigned i=0; i<n-1; i++)
	{
		sum_charge -= atoms_[i].charge;
	}
	atoms_[n-1].charge = sum_charge;

	double energy = 0., ch_;
	for (unsigned i=0; i<n; i++)
	{
		ch_ = atoms_[i].charge;
		energy +=  ch_ * atoms_[i].electronegativity;
		energy += 0.5 * ch_ * ch_ * _E(n, &coulomb_[0], i, i);
		for (unsigned j=0; j<i; j++)
			energy += ch_ * atoms_[j].charge * _E(n, &coulomb_[0], i, j);
	}

	for (unsigned i=0; i<n; i++)
	{
		atoms_[i].du__dq = atoms_[i].electronegativity;
		//
		//atoms_[i].du__dq +=  atoms_[i].charge * _E(n, &coulomb_[0], i, i);
		//
		for (unsigned j=0; j<n; j++)
			atoms_[i].du__dq += atoms_[j].charge * _E(n, &coulomb_[0], i, j);
	}

	for (unsigned i=0; i<n-1; i++) // correction
	{
		atoms_[i].du__dq -= atoms_[n-1].du__dq;
	}
	return energy;
}

/**
* @brief geometry optimizer for molecule (or complex)
*/
class QEq_optimizer
{
	std::vector<double> q_, du__dq_;

public:

	double operator()(Molecule *molecule, unsigned max_iter=2000)
	{
		typedef Molecule::atom_type       _Atom;
		_Atom *atoms = molecule->get_atoms();

		unsigned n = molecule->get_atom_count() - 1;

		q_.resize(n);
		du__dq_.resize(n);

		//for (unsigned i=0; i<n; i++) q_[i] = atoms[i].charge;

		double energy = prgkern::Minimizer_<prgkern::LMBFGS_>()(&dU__dQ,
			n, &q_[0], &du__dq_[0], (void*)molecule, max_iter, 0.1, 0.1, 1e-4, 1e-16, 1e-16);

		double full_charge = 0.;
		for (unsigned i=0; i<n+1; i++)
		{
			full_charge += atoms[i].charge;
			printf("i=%d CHARGE OPTIMIZE=%10.5e\n", i, atoms[i].charge);
		}
		printf("FULL CHARGE =%10.5e\n", full_charge);
		printf("FULL energy =%10.5e Hartree\n", energy);

		return energy;
	}

	static double dU__dQ(unsigned n, const double *q, double *du__dq, void *param)
	{
		typedef Molecule::atom_type       _Atom;
		Molecule *molecule = (Molecule *)param;
		_Atom *atoms = molecule->get_atoms();

		n = molecule->get_atom_count() - 1;
		for (unsigned i=0; i<n; i++) atoms[i].charge = q[i];
		double energy = molecule->dU__dQ();
		for (unsigned i=0; i<n; i++) du__dq[i] =atoms[i].du__dq;

		return energy;
	}
};
#undef _E

int main(int argc, char *argv[]) 
{
	if (argc == 1)
	{
		printf("usege : QTPIE directory *.xyz\n");
		exit(0);
	}
	
	Molecule molecule;
	molecule.load(argv[1]);
	QEq_optimizer()(&molecule);
	
	return 0;
}
