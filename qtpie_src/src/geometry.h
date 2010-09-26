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
#ifndef GEOMETRY__H
#define GEOMETRY__H

struct Point
{
	double &operator[](unsigned i) { return _[i]; }
	double _[3];
};

/**
* @brief Computes pairwise distances from Cartesian coordinates
* @param Point1,Point2: 3-vectors of double precisions
* @return Cartesian distance in atomic units
*/
inline double Distance(const double Point1[3], const double Point2[3])
{
	double x, y, z;
	x = Point2[0] - Point1[0];
	y = Point2[1] - Point1[1];
	z = Point2[2] - Point1[2];
	return sqrt(x*x + y*y + z*z);
}

#endif
