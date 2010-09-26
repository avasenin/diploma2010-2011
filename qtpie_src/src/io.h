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
#ifndef IO__H
#define IO__H

#include <stdio.h>
#include <stdlib.h>
#include "parameters.h"

/**
* @brief dumps QTPIE calculation results
* @param Mol: molecule data structure
* @param fileName Name of the log file to write or append to
*/
inline void WriteLog(const Molecule &Mol, const char *fileName)
{
	FILE *file = fopen(fileName, "wt");
	if (file == 0)
	{
		printf("Problem open log file %s\n", fileName);
		exit(0);
	}
	fprintf(file, "energy is %lf\n", Mol.Energy);
	for (int i=0; i<Mol.NumAtoms; i++)
		fprintf(file, "%5d %8.3lf\n", i, Mol.Atoms[i].Charge);
	
	fclose(file);
}

/**
* @brief dumps molecular geometry from QTPIE in XYZ formal
* @param Mol: molecule data structure
* @param fileName: Name of geometry file to write or append to 
*/
inline void WriteXYZ(const Molecule &Mol, const char *fileName)
{
	FILE *file = fopen(fileName, "wt");
	if (file == 0)
	{
		printf("Problem open dump file %s\n", fileName);
		exit(0);
	}
	fprintf(file, "%d\n", Mol.NumAtoms);
	fprintf(file, "Written by QTPIE : WriteXYZ()\n");
	for (int j=0; j<Mol.NumAtoms; j++)
		fprintf(file, "%2s %10.3lf %10.3lf %10.3lf\n", Mol.Atoms[j].Element.Symbol, 
			Mol.Atoms[j].Position[0]/Angstrom, Mol.Atoms[j].Position[1]/Angstrom,
			Mol.Atoms[j].Position[2]/Angstrom);
	
	fclose(file);
}

#endif
