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
#ifndef ATOMICUNITS__H
#define ATOMICUNITS__H

/**
*  Our QTPIE charge model works exclusively in atomic units
*  The values stored here are conversion factors to convert
*  from that unit into atomic units
*  There is no dimensional checking implemented!
*/
const double ONE         = 1.;
const double ZERO        = 0.;
const double eV          = 3.67493245e-2;      // electron volt to Hartree
const double kJ_mol      = 6.6744644952e-3;    // kilojoule per mole to Hartree
const double kcal_mol    = 1.5952353e-3;       // kilocalorie per mole to Hartree
const double invAngstrom = 455.6335252760;     // inverse Angstrom to Hartree
const double Debye       = 0.3934302014076827; // Debye to atomic unit of dipole moment
const double Angstrom    = 1./0.529177249;     // Angstrom to bohr

#endif
