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
#ifndef PROPERTIES__H
#define PROPERTIES__H

#include "parameters.h"

/**
* @brief Computes the dipole moment
* @param Centroid Point about which to define dipole moment
* @return the dipole moment vector
*/
Point dipmom(const Molecule &Mol)
{
	double Centroid[3];
	Point dipmom_;
	for (unsigned i=0; i<Mol.NumAtoms; i++)
		for (unsigned j=0; j<3; j++)
			dipmom_[j] += Mol.Atoms[i].Charge 
				* (Mol.Atoms[i].Position[j] - Centroid[j]);
// [?ERROR?] Centroid[j] using without initialization
	return dipmom_;
}

/**
* @brief Computes the dipole polarizability tensor
* @note currently returns 0
*/
double polarizability()
{
	//double precision, dimension(3:3) :: polarizability
	return 0.;
}

#endif
