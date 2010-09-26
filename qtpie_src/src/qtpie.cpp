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
#include "parameters.h"
#include "molecule.h"
#include "cg.h"
#include "io.h"
#include "qtpie.h"
#include "dqtpie.h"

/// Runs QTPIE for a single XYZ geometry
int main(int argc, char *argv[])
{
	if (argc == 1)
	{
		printf("usege : QTPIE *.xyz\n");
		exit(0);
	}
	
	Molecule Mol;
	double elapsed[2];
	double total = 0;
	
	Mol.load(argv[1]);
	printf("Read file %s\n", argv[1]);
	
	DosGTOIntegrals(Mol);
	QEq(Mol);
	printf("Integrals done...\n");
	
	QTPIE(Mol);
	printf("Calculations done...\n");
	printf("QTPIE Energy is %lf\n", Mol.Energy);
	
	WriteLog(Mol, "qtpie.log");
	printf("Calculated charges written to qtpie.log\n");
	
	return 0;
}
