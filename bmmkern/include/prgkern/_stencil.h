/***************************************************************************
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
#ifndef _STENCIL__0077A726_2C84_57b4_2A66_CE43B6700100__H
#define _STENCIL__0077A726_2C84_57b4_2A66_CE43B6700100__H

#include "prgkern/_prgconfig.h"
#include "prgkern/_assert.h"
#include "prgkern/_math.h"
#include "prgkern/_pproc.h"
#include "prgkern/_index.h"

namespace prgkern
{

	/**
	* @brief tool for smoothing data in multidimensional arrays
	* @param N - dimension of stencil
	* @param M - number of points for smoothing
	*/
	template <unsigned N, unsigned M, typename T> class Stencil_;

	#define __EQ_xyz(T, arr, inc, ndxcast, _Index3, ndx0, i, j, k, op, value) \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i),  (j),  (k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i),  (j), -(k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i), -(j),  (k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i),  (j),  (k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i), -(j), -(k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i),  (j), -(k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i), -(j),  (k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i), -(j), -(k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i),  (k),  (j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i),  (k), -(j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i), -(k),  (j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i),  (k),  (j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i), -(k), -(j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i),  (k), -(j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i), -(k),  (j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i), -(k), -(j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (j),  (i),  (k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (j),  (i), -(k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (j), -(i),  (k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(j),  (i),  (k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (j), -(i), -(k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(j),  (i), -(k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(j), -(i),  (k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(j), -(i), -(k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (j),  (k),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (j),  (k), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (j), -(k),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(j),  (k),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (j), -(k), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(j),  (k), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(j), -(k),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(j), -(k), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (k),  (i),  (j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (k),  (i), -(j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (k), -(i),  (j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(k),  (i),  (j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (k), -(i), -(j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(k),  (i), -(j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(k), -(i),  (j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(k), -(i), -(j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (k),  (j),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (k),  (j), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (k), -(j),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(k),  (j),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (k), -(j), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(k),  (j), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(k), -(j),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(k), -(j), -(i))]) op (value);

	#define __EQ_xxz(T, arr, inc, ndxcast, _Index3, ndx0, i, _, k, op, value) \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i),  (i),  (k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i),  (i), -(k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i), -(i),  (k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i),  (i),  (k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i), -(i), -(k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i),  (i), -(k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i), -(i),  (k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i), -(i), -(k))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i),  (k),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i),  (k), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i), -(k),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i),  (k),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i), -(k), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i),  (k), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i), -(k),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i), -(k), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (k),  (i),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (k),  (i), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (k), -(i),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(k),  (i),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (k), -(i), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(k),  (i), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(k), -(i),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(k), -(i), -(i))]) op (value);

	#define __EQ_xxx(T, arr, inc, ndxcast, _Index3, ndx0, i, _, __, op, value) \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i),  (i),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i),  (i), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i), -(i),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i),  (i),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i), -(i), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i),  (i), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i), -(i),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i), -(i), -(i))]) op (value);

	#define __EQ_xy0(T, arr, inc, ndxcast, _Index3, ndx0, i, j, _, op, value) \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i),  (j),    0)]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i), -(j),    0)]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i),  (j),    0)]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i), -(j),    0)]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i),    0,  (j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i),    0, -(j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i),    0,  (j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i),    0, -(j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (j),  (i),    0)]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (j), -(i),    0)]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(j),  (i),    0)]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(j), -(i),    0)]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (j),    0,  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (j),    0, -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(j),    0,  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(j),    0, -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(   0,  (i),  (j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(   0,  (i), -(j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(   0, -(i),  (j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(   0, -(i), -(j))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(   0,  (j),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(   0,  (j), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(   0, -(j),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(   0, -(j), -(i))]) op (value);

	#define __EQ_xx0(T, arr, inc, ndxcast, _Index3, ndx0, i, _, __, op, value) \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i),  (i),    0)]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i), -(i),    0)]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i),  (i),    0)]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i), -(i),    0)]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i),    0,  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i),    0, -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i),    0,  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i),    0, -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(   0,  (i),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(   0,  (i), -(i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(   0, -(i),  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(   0, -(i), -(i))]) op (value);

	#define __EQ_x00(T, arr, inc, ndxcast, _Index3, ndx0, i, _, __, op, value) \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3( (i),    0,    0)]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(-(i),    0,    0)]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(   0,  (i),    0)]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(   0, -(i),    0)]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(   0,    0,  (i))]) op (value); \
		*(T*)((char*)arr + inc * ndxcast[ndx0 + _Index3(   0,    0, -(i))]) op (value);

	#define __EQ_000(T, arr, inc, ndxcast, _Index3, ndx0, _, __, ___, op, value) \
		*(T*)((char*)arr + inc * ndxcast[ndx0]) op (value);

	/**
	*     Распределяет значение по 27 ближайшим точкам сетки. Сетка является
	*     внешней, потому один стенсил может работать с любым их количеством.
	*
	* (1) Рассматривает одномерный массив как 3-х мерный с границами
	*     [I1..I2), [J1..J2), [K1..K2)
	* (2) распределяет заданное значение у заданной точки (i0, j0, k0) равномерно
	*     по всем направлениям
	*/
	template <typename T> class Stencil_<3, 3, T>
	{
		typedef index_<3, int>  _Index;
		T fn_[4]; // предвычисленные значения f()

	public:

		/**
		* @param fn - функция смазывания, зависящая только от расстояния
		*/
		template <typename F> Stencil_(F fn)
		{
			fn_[0] = fn(T());      // (0, 0, 0) -> 1 element
			fn_[1] = fn(1.);       // (1, 0, 0) -> 6 elements
			fn_[2] = fn(M_SQRT_2); // (1, 1, 0) -> 12 elements
			fn_[3] = fn(M_SQRT_3); // (1, 1, 1) -> 8 elements
			T norm = 1. / (fn_[3] * 8 + fn_[2] * 12 + fn_[1] * 6 + fn_[0]);
			fn_[0] *= norm;
			fn_[1] *= norm;
			fn_[2] *= norm;
			fn_[3] *= norm;
		}

		/**
		*   Возвращает коэффициент смазывания
		* @param ndx индекс смещения от центра
		*/
		T operator[](const _Index &ndx) const
		{
			assert(_GE(ndx[0], -1)); assert(_LE(ndx[0], 1));
			assert(_GE(ndx[1], -1)); assert(_LE(ndx[1], 1));
			assert(_GE(ndx[2], -1)); assert(_LE(ndx[2], 1));
			return fn_[abs(ndx[2]) + abs(ndx[1]) + abs(ndx[0])];
				// this is trick due to the next table
				// (1, 1, 1) -> fn_[3]
				// (1, 1, 0) -> fn_[2]
				// (1, 0, 0) -> fn_[1]
				// (0, 0, 0) -> fn_[0]
		}

		/**
		* @brief вставка значения и его размазка по линейному массиву памяти
		* @note Внешнаяя структура памяти, куда вставляются значения, должна быть обязательно
		*   вектором, который занимает сплошную область памяти. Однако, элементами этого
		*   вектора могут быть любые структуры данных. По этой причине, необходимо указывать
		*   адрес того элемента структуры, которая является первой в массиве, и указывать
		*   смещение между двумя элементами сетки.
		* @param arr адрес поля структуры первого объекта в линейном массиве структур в памяти
		* @param inc смещение в байтах между последовательными элементами &arr[i].f - &arr[i-1].f
		* @param ndxcast объект, который "знает" как преобразовать мультиинджекс в линейный индекс
		* @param ndx точка (i,j,k) сетки, к которой привязывается центр стенсила
		* @param value вставляеме значение
		*/
		template <typename C, typename S>
		void insert(S *arr, unsigned inc, const C &ndxcast, const _Index &ndx, S value) const
		{
			__EQ_xxx(S, arr, inc, ndxcast, _Index, ndx, 1, 1, 1, +=, expr2value(value * fn_[3]));
			__EQ_xx0(S, arr, inc, ndxcast, _Index, ndx, 1, 1, 0, +=, expr2value(value * fn_[2]));
			__EQ_x00(S, arr, inc, ndxcast, _Index, ndx, 1, 0, 0, +=, expr2value(value * fn_[1]));
			__EQ_000(S, arr, inc, ndxcast, _Index, ndx, 0, 0, 0, +=, expr2value(value * fn_[0]));
		}
	};

	template <typename T> class Stencil_<3, 5, T>
	{
		typedef index_<3, int>  _Index;
		T fn_[10]; // precalculated f() values

	public:

		/**
		* @brief construct and precalculate internal normalized arrays
		* @param F - distance dependent functor
		* @param h - step
		*/
		template <typename F> Stencil_(F fn)
		{
			fn_[0] = fn(T(0.));        // 000 ->  1 element
			fn_[1] = fn(T(1.));        // 100 ->  6 elements
			fn_[2] = fn(T(M_SQRT_2));  // 110 -> 12 elements
			fn_[3] = fn(T(M_SQRT_3));  // 111 ->  8 elements
			fn_[4] = fn(T(2.));        // 200 ->  6 elements
			fn_[5] = fn(T(M_SQRT_5));  // 210 -> 24 elements
			fn_[6] = fn(T(M_SQRT_6));  // 211 -> 24 elements
			fn_[7] = fn(T(M_SQRT_8));  // 220 -> 12 elements
			fn_[8] = fn(T(3.));        // 221 -> 24 elements
			fn_[9] = fn(T(M_SQRT_12)); // 222 ->  8 elements

			T norm = (fn_[5] + fn_[6] + fn_[8]) * 24 + (fn_[3] + fn_[9]) * 8
				+ (fn_[2] + fn_[7]) * 12 + (fn_[1] + fn_[4]) * 6 + fn_[0];

			norm = 1./ norm;
			fn_[0] *= norm;
			fn_[1] *= norm;
			fn_[2] *= norm;
			fn_[3] *= norm;
			fn_[4] *= norm;
			fn_[5] *= norm;
			fn_[6] *= norm;
			fn_[7] *= norm;
			fn_[8] *= norm;
			fn_[9] *= norm;
		}

		T operator[](const _Index &ndx) const
		{
			assert(_GE(ndx[0], -1)); assert(_LE(ndx[0], 1));
			assert(_GE(ndx[1], -1)); assert(_LE(ndx[1], 1));
			assert(_GE(ndx[2], -1)); assert(_LE(ndx[2], 1));

			// calculate indexes for searching fn_
			int s = abs(ndx[2]) + abs(ndx[1]) + abs(ndx[0]);
			int p = 1;
			if (ndx[2]) p *= ndx[2];
			if (ndx[1]) p *= ndx[1];
			if (ndx[0]) p *= ndx[0];
			p = abs(p);

			switch (s)
			{
				case 6 : return fn_[9];
				case 5 : return fn_[8];
				case 4 : if (p == 4) return fn_[7]; else return fn_[6];
				case 3 : if (p == 2) return fn_[5]; else return fn_[3];
				case 2 : if (p == 2) return fn_[4]; else return fn_[2];
				case 1 : return fn_[1];
				case 0 : return fn_[0];
			}
			return 0;
		}

		/**
		* @brief inserts the value into the arr according to stensil dispersion
		* @param ndxcast - object for transform multindex [i,j,k] into linear memory offset
		* @param arr - array which are inserted the values
		* @param inc - memory offset (bytes) between elements of array inc = &arr[i] - &arr[i-1]
		* @param ndx - multindex (i,j,k)
		* @param value - value inserted
		*/
		template <typename C, typename S>
		void insert(S *arr, unsigned inc, const C &ndxcast, const _Index &ndx, S value) const
		{
			__EQ_xxx(S, arr, inc, ndxcast, _Index, ndx, 2, 2, 2, +=, expr2value(value * fn_[9]));
			__EQ_xxz(S, arr, inc, ndxcast, _Index, ndx, 2, 2, 1, +=, expr2value(value * fn_[8]));
			__EQ_xx0(S, arr, inc, ndxcast, _Index, ndx, 2, 2, 0, +=, expr2value(value * fn_[7]));
			__EQ_xxz(S, arr, inc, ndxcast, _Index, ndx, 1, 1, 2, +=, expr2value(value * fn_[6]));
			__EQ_xy0(S, arr, inc, ndxcast, _Index, ndx, 2, 1, 0, +=, expr2value(value * fn_[5]));
			__EQ_x00(S, arr, inc, ndxcast, _Index, ndx, 2, 0, 0, +=, expr2value(value * fn_[4]));
			__EQ_xxx(S, arr, inc, ndxcast, _Index, ndx, 1, 1, 1, +=, expr2value(value * fn_[3]));
			__EQ_xx0(S, arr, inc, ndxcast, _Index, ndx, 1, 1, 0, +=, expr2value(value * fn_[2]));
			__EQ_x00(S, arr, inc, ndxcast, _Index, ndx, 1, 0, 0, +=, expr2value(value * fn_[1]));
			__EQ_000(S, arr, inc, ndxcast, _Index, ndx, 0, 0, 0, +=, expr2value(value * fn_[0]));
		}
	};

	#undef __EQ_xxx
	#undef __EQ_xxz
	#undef __EQ_xx0
	#undef __EQ_xxz
	#undef __EQ_xy0
	#undef __EQ_x00
	#undef __EQ_xxx
	#undef __EQ_xx0
	#undef __EQ_x00
	#undef __EQ_000

	/**
	*  Объект, который размазывает введенное значение по ряду точек сетки.
	*  Данный объект отличается от других подобных тем, что размазка идет не по
	*  фиксированному набору точек, а по динамическому набору. Размерность
	*  подсетки для размазывания меняется в зависимости от выбранного радиуса
	*  взаимодействия и длины шага сетки.
	*  Данная реализация ограничена равносторонним характером внутренней подсетки.
	* @param T тип реального числа
	*/
	template <typename T> class Stencil_<3, 0, T>
	{
		typedef index_<3, int>  _Index;
		typedef vdense_<3, T>   _Vector3;

		mutable std::vector<T> fn_; ///< внутренная подсетка
		unsigned ndx_; ///< размерность внутренней сетки (по всем направлениям едина)
		T h_; ///< шаг (под)сетки для каждого направления
		T radius2_; ///< квадрат радиуса взаимодействия
		T tune_param_; ///< параметр тьюнинга внутренней функции

	public:

		/**
		* @brief конструирование объекта
		* @param h - шаг сетки
		* @param radius - радиус сглаживания (взаимодействия)
		*/
		Stencil_(T h, T radius) : h_ (h), radius2_(sqr(radius)), tune_param_(1.)
		{
			ndx_ = (int)::floor(radius / h) * 2 + 1;
			fn_.resize(cube(ndx_));
			smooth_fn_tune_(radius2_); // настройка внутренную функцию на границу
		}

		/**
		* @brief вставляет размазанное значение во внешнюю сетку
		* @param arr - внешняя сетка
		* @param inc - разность (в байтах) между последовательными элементами внешней сетки
		* @param ndxcast - объект, конвертирующий трехмерный индекс в линейный для внешней сетки
		* @param ndx - трехмерный индекс внешней сетки, около которого "складывать" значения
		* @param X - абсолютная координата точки относительно ближайшей точки сетки
		* @param value - вставляемое значение
		*/
		template <typename C>
		void insert(T *arr, unsigned inc, const C &ndxcast,
			const _Index &ndx, const _Vector3 &X, T value) const
		{
			insert_to_internal_net_(X, value);
			copy_internal_net_(arr, inc, ndxcast, ndx);
		}

	private:

		/**
		*  Скалярная функция сглаживания значения f(x) = exp(-x*x/sigma).
		*  В данной реализации она определена как внутренная, а не внешняя поскольку:
		*   (1) она требует настройки под область взаимодействия
		*   (2) не требует установки пользователем
		*  Функция с помощью параметра tune_param_ подстраивается таким образом,
		*  что на границе взаимодействия значение функции становится близким к 0.
		*  При необходимости можно написать любое количество функций и подменять
		*  их на этапе компиляции через #ifdef определения.
		* @param x2 - квадрат расстояния
		* @param param - параметр для подстройки функции под границу взаимодействия
		* @return значение в точке x
		*/
		T smooth_fn_(T x2, T param) const { return exp( x2 * param); }

		/**
		* Настройка внутренней функции сглаживания на границу взаимодействия.
		* Функция согласована с функцией smooth_fn_()
		* @param xlimit2 - радиус взаимодействия (в квадрате)
		* @param accuracy - относительный уровень падения на границе
		*/
		void smooth_fn_tune_(T xlimit2, T accuracy=0.1)
		{ tune_param_ = log(accuracy) / xlimit2; }

		/**
		* @brief вставляет значение и размазывает его по точкам внутренней сетки
		* @note автоматически очищает внутренную сетку от предыдущих значений
		*  Реализация использует способ прохода сетки и построения координаты
		*  текущего узла эффективным способом, исключающим умножения.
		* @param fn - функция размазки (работает с квадратом(!) расстояния)
		* @param X - относительная координата точки относительно ближайшей точки сетки
		* @param value - вставляемое значение
		*/
		void insert_to_internal_net_(const _Vector3 &X, T value) const
		{
			_Vector3 X__ = X * h_; // абсолютный вектор из относительного
			X__ += (ndx_ >> 1) * h_; // расстояние относительно начала внутренней сетки
			unsigned n = 0; // согласованный с порядком обхода линейный индекс массива fn
			T norm = 0.; // коэффициент нормировки значений
			_Vector3 Y__(0., 0., 0.); // абсолютная координата узла сетки
			for (unsigned i=0; i<ndx_; i++)
			{
				Y__[1] = 0.;
				for (unsigned j=0; j<ndx_; j++)
				{
					Y__[2] = 0.;
					for (unsigned k=0; k<ndx_; k++)
					{
						T r2 = distance2(X__, Y__);
						T s = 0.; // очистка от предыдущего значения
						if (r2 < radius2_) s = smooth_fn_(r2, tune_param_); // новое значение
						norm += s;
						fn_[n++] = s;
						Y__[2] += h_;
					}
					Y__[1] += h_;
				}
				Y__[0] += h_;
			}
			norm = value / norm;
			for (unsigned i=0; i<n; i++) fn_[i] *= norm; // нормировка на value
		}

		/**
		* Добавляет значения из внутренней сетки во внешнюю.
		* @param arr - внешняя сетка
		* @param inc - разность (в байтах) между последовательными элементами внешней сетки
		* @param ndxcast - объект, конвертирующий трехмерный индекс в линейный для внешней сетки
		* @param ndx - трехмерный индекс внешней сетки, около которого "складывать" значения
		*/
		template <typename C>
		void copy_internal_net_(T *arr, unsigned inc, const C &ndxcast,
			const _Index &ndx) const
		{
			_Index ndx__(ndx[0] - ndx_/2, ndx[1] - ndx_/2, ndx[2] - ndx_/2);
				// точка привязки внутри внешней сетки
			unsigned n = 0; // линейный индекс внутренней сетки
			for (unsigned i=0; i<ndx_; i++)
			for (unsigned j=0; j<ndx_; j++)
			for (unsigned k=0; k<ndx_; k++)
			{
				*(T*)((char*)arr + inc * ndxcast[ndx__ + _Index(i, j, k)]) += fn_[n++];
			}
		}
	};

	/**
	* @brief smooth the values of array
	* @param stensil - object for dispersing of inserted values
	* @param array - array which are inserted the values
	* @param inc - memory offset (bytes) between elements of array inc = &arr[i] - &arr[i-1]
	* @param ndxcast - object for transform multindex
	*/
	template <typename S, typename C1, typename C2, typename T>
	INLINE void smooth(const S &stencil, T *array_to, unsigned inc_to, const C1 &ndxcast_to,
		const T *array_from, unsigned inc_from, const C2 &ndxcast_from)
	{
		for (unsigned i=0, sz=ndxcast_from.size(); i<sz; i++)
		{
			T value = *(T*)((char*)array_from + inc_from * i);
			typename C2::multi_index_type ndx = ndxcast_from[i];
			stencil.insert(array_to, inc_to, ndxcast_to, ndx, value);
		}
	}

	/**
	* @brief disperse given value at X point along the nearest mesh points
	* @param X - relative coordinates of point (step of net = 1(!))
	*/
	template <typename C, typename T>
	INLINE void insert(T *arr_to, unsigned inc_to, const C &ndxcast_to, const index_<3, int> &ndx_to,
		const vdense_<3, T> &X, T value)
	{
		vdense_<3, T> X__ =  vdense_<3, T>(1.) - X; // inserting part to nearest node

	#define __V(j, i, k)  *(T*)((char*)arr_to + inc_to * ndxcast_to[ndx_to + index_<3, int>((i), (j), (k))])
		T s00 = value * X__[0] * X__[1]; __V(0, 0, 0) += s00 * X__[2]; __V(0, 0, 1) += s00 * X[2];
		T s01 = value * X__[0] * X  [1]; __V(0, 1, 0) += s01 * X__[2]; __V(0, 1, 1) += s01 * X[2];
		T s10 = value * X  [0] * X__[1]; __V(1, 0, 0) += s10 * X__[2]; __V(1, 0, 1) += s10 * X[2];
		T s11 = value * X  [0] * X  [1]; __V(1, 1, 0) += s11 * X__[2]; __V(1, 1, 1) += s11 * X[2];
	#undef __V
	}

	/**
	* @brief insert, disperse and smooth the given value at X point
	* @param X - relative coordinates of point (step of net = 1(!))
	*/
	template <unsigned N, typename C, typename T>
	INLINE void insert(const Stencil_<3, N, T> &stencil, T *arr_to, unsigned inc_to, const C &ndxcast_to,
		const index_<3, int> &ndx_to, const vdense_<3, T> &X, T value)
	{
		typedef index_<3, int>  _Index;
		vdense_<3, T> X__ =  vdense_<3, T>(1.) - X; // inserting part to nearest node

		T s = value * X__[0] * X__[1];
		stencil.insert(arr_to, inc_to, ndxcast_to, ndx_to + _Index(0, 0, 0), s * X__[2]);
		stencil.insert(arr_to, inc_to, ndxcast_to, ndx_to + _Index(0, 0, 1), s * X[2]);

		s = value * X__[0] * X[1];
		stencil.insert(arr_to, inc_to, ndxcast_to, ndx_to + _Index(0, 1, 0), s * X__[2]);
		stencil.insert(arr_to, inc_to, ndxcast_to, ndx_to + _Index(0, 1, 1), s * X[2]);

		s = value * X[0] * X__[1];
		stencil.insert(arr_to, inc_to, ndxcast_to, ndx_to + _Index(1, 0, 0), s * X__[2]);
		stencil.insert(arr_to, inc_to, ndxcast_to, ndx_to + _Index(1, 0, 1), s * X[2]);

		s = value * X[0] * X[1];
		stencil.insert(arr_to, inc_to, ndxcast_to, ndx_to + _Index(1, 1, 0), s * X__[2]);
		stencil.insert(arr_to, inc_to, ndxcast_to, ndx_to + _Index(1, 1, 1), s * X[2]);
	}

}
#endif
