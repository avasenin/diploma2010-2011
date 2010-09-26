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
#ifndef MOLECULE__H
#define MOLECULE__H

#include "parameters.h"
#include "gto-int.h"
#include "sto-int.h"

/**
* @brief Describes a molecular system
* @param Description: a text label of 132 characters
* @param    NumAtoms: number of atoms (integer)
* @param TotalCharge: total charge of system (double precision)
* @param       Atoms: array of atoms
* @param     Overlap: overlap matrix
* @param      OvNorm: overlap norm vector (useful temporary variable)
* @param     Coulomb: Coulomb matrix
* @param      Energy: QTPIE contribution to the potential energy
* @param   EGradient: Energy gradients
*/
struct Molecule
{
	char Description[132];
	int NumAtoms;
	double TotalCharge;
	double Energy;
	Atom *Atoms;
	double *Overlap; // dim = 2
	double *OvNorm;
	double *Coulomb; // dim = 2
	Point *EGradient;
	
	Molecule() : NumAtoms(0), TotalCharge(0.), Energy(0.),
		Atoms(0), Overlap(0), OvNorm(0), Coulomb(0), EGradient(0) {}
		
	~Molecule()
	{
		delete [] Atoms;
		delete [] Overlap;
		delete [] OvNorm;
		delete [] Coulomb;
		delete [] EGradient;
	}
	
	void load(const char *fileName);
};

/**
* @brief Reads XYZ file
*  loads an external file containing a XYZ geometry
*  \param fileName: name of the XYZ geometry file
*  \return A Molecule data structure
*/
inline void Molecule::load(const char *fileName)
{
	char buffer[255];
	char AtomSymbol[2];

	bool isParameterized;
	FILE *file = fopen(fileName, "rt");
	if (file == 0)
	{
		printf("Problem loading geometry file %s\n", fileName);
		exit(0);
	}
	
	// First line says how many atoms there are
	fgets(buffer, 255, file);
	sscanf(buffer, "%d", &NumAtoms);
	// Second line may contain a comment, skip it
	fgets(buffer, 255, file);
	
	Atoms = new Atom [NumAtoms];
	Overlap = new double [sqr(NumAtoms)];
	Coulomb = new double [sqr(NumAtoms + 1)];
	OvNorm  = new double [NumAtoms];
	EGradient = new Point [NumAtoms];
	
	for (int j=0; j<NumAtoms; j++)
	{
		fgets(buffer, 255, file);
		sscanf(buffer, "%s %lf %lf %lf", &AtomSymbol, &Atoms[j].Position[0],
			&Atoms[j].Position[1], &Atoms[j].Position[2]);
			
		// Convert units from Angstroms to atomic units (Bohr)
		Atoms[j].Position[0] *= Angstrom;
		Atoms[j].Position[1] *= Angstrom;
		Atoms[j].Position[2] *= Angstrom;
		
		// look up AtomSymbol to assign parameters
		isParameterized = false;
		for (int l=0; l<numParameterizedAtoms; l++)
		{
			if (::strncmp(AtomSymbol, ParameterizedAtoms[l].Symbol, 2) == 0)
			{
				Atoms[j].Element = ParameterizedAtoms[l];
				isParameterized = true;
			}
		}
		
		// assign basis set
		if (isParameterized){
			AssignsGTOBasis(Atoms[j]);
			// Replace with this line to assign STO
			// AssignsSTOBasis(Atoms[j]);
			}
		else
		{
			printf("Error: Unknown element type: %s\n", AtomSymbol);
			exit(0);
		}
	}
	fclose(file);
}

#endif
