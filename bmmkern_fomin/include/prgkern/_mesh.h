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
#ifndef _MESH__0077A726_56A4_5b9b_E3D3_D1451DA80101__H
#define _MESH__0077A726_56A4_5b9b_E3D3_D1451DA80101__H
#include "prgkern/_prgconfig.h"

#include "prgkern/_string.h"
#include "prgkern/_assert.h"
#include "prgkern/_index.h"
#include "prgkern/_v3dense.h"
#include "prgkern/_stencil.h"
#include "prgkern/_transforms.h"
#include <algorithm>

namespace prgkern
{

	/**
	*   Mesh - это объект, который структурирует набор одинаковых сеток.
	*   Каждая сетка является простым неструктурированным линейным массивом в памяти.
	*   Она может быть задана, например std::vector<_Real>, или через указатель на память.
	*   Эти сетки хранят информацию о пространственно распределенных данных - как
	*   найти данные для точки, соответствующей заданной пространственной точке,
	*   в данном массиве, или как решить обратную задачу и определяет этот объект.
	*   Предполагается цикличность данных за пределами заданной области пространства,
	*   чтобы избежать проблем вставки данных за пределы области.
	* @note Область индексов определена в STL стиле, то есть как [I1, I2[
	* @param N размерность пространства
	* @param _Real тип плаваюего числа
	*/
	template <unsigned N, typename _Real>
	class Mesh_
	{
		typedef vdense_<N, _Real>    _Point; ///< координаты точки
		typedef index_<N, int>       _Index; ///< многомерный индекс точки
		typedef Box_<N, _Real>       _Box;   ///< область пространства

		typedef vdense_<2, _Real>    _Point2;
		typedef vdense_<3, _Real>    _Point3;
		typedef vdense_<4, _Real>    _Point4;

		typedef index_<2, int>       _Index2;
		typedef index_<3, int>       _Index3;
		typedef index_<4, int>       _Index4;

		index_<N, unsigned> sz_; // пользовательская размерность сетки
			// по z она отличается от внутренней за счет паддинга
		unsigned size_; // должное число элементов в каждой внешней сетке
		unsigned pad_; // z-паддинг внешней сетки
		_Real h_, _1h_; // шаг и обратный шаг сетки

	public:

		typedef _Index  index_type;

		/**
		*  Размечает объект, получает размер необходимой памяти.
		* @param interaction_radius длина грани базовой ячейки
		* @param sz число базовых ячеек по всем направлениям
		* @param h приблизительный шаг деления ячейки
		*  Шаг уточняется согласно разрешенному числу делений каждой грани.
		* @param r{ight}pad паддинг по последнему измерению (число дополнительных точек).
		*  Паддинг является величиной определяемой пользователем, зависит от задачи
		*  и используется для эффективной записи элементов. Его наличие позволяет
		*  записывать элементы не по одному, а z-рядами. При этом происходит
		*  запись как бы за границу массива. По окончании такой записи множества
		*  элементов необходимо свернуть массив, то есть перенести элементы,
		*  записанные за границу внутрь массива. Такая техника работает быстрее,
		*  чем запись с постоянным контролем. Паддинга по другим координатам
		*  легко избежать. Внутри массив паддится дважды (впереди и сзади).
		* @return число точек в объекте (необходимый размер памяти под сетку)
		*/
		unsigned init(_Real interaction_radius, index_<N, unsigned> sz, _Real h,
			unsigned rpad=1)
			// паддинг равный 1 по умалчиванию нужен для быстрой работы функции
			// вставки, которая размазывает значение по ближайшим 8 точкам
		{
			pad_ = rpad;
				// размер паддинга слева и справа различаются, слева на 1 меньше,
				// так как справа граница не входит в интервал, а левая входит.

			_S msg("Have been builded -> mesh [");
			unsigned dims_count = STATIC_DIMENSION(dfft_sizes, unsigned);
				// используем только разрешенные размеры сеток, поскольку для
				// других размеров расчеты могут быть неэффективными (dfft ограничения)

			unsigned dim = (unsigned)ceil(interaction_radius / h);
			const unsigned *p = std::upper_bound(&dfft_sizes[0],
				&dfft_sizes[0] + (int)dims_count, dim);
			if ((unsigned)(p - &dfft_sizes[0]) < dims_count)
			{
				dim = *p; h_ = interaction_radius / dim;
				_1h_ = 1. / h_;
			}
			else
			{
				_S msg("[ERROR] mesh size is out of range [1, 10125]");
				PRINT_ERR(msg);
			}

			size_ = 1;
			for (unsigned i=0; i<N-1; i++)
			{
				sz_[i] = sz[i] * dim;
				size_ *= sz_[i];
				msg += make_string("%d, ", sz_[i]);
			}
			sz_[N-1] = sz[N-1] * dim;
			size_ *= sz[N-1] * dim + pad_ + pad_ - 1;
				// на величину паддинга увеличим только размер требуемого внешнего массива,
				// пользователь о нем как бы не знает, но работать с таким массивом он
				// напрямую не должен, а только через mesh, который учитывает наличие
				// паддинга при обращении к элементам массива.

			msg += make_string("%d]", sz_[N-1]);
			PRINT_MESSAGE(msg);

			return size_;
		}

		/**
		* Должное число элементов во внешних массивах для данных сеток.
		* @return размер массивов сеток
		*/
		unsigned size() const { return size_; }

		/**
		* Размерность сеток по всем направлениям.
		* @return многомерный вектор, описывающих размерности
		*/
		_Index dimension() const { return sz_; }

		/**
		*  Вставляет значение в сетку и "распыляет" его по ближайшим соседям
		* @param fn сетка, куда вставляется значение
		* @param X координата точки
		* @param value вставляемое значение
		*/
		void insert(_Real *fn, const _Point3 &X, _Real s)
		{
			assert(_GE(pad_, (unsigned)1));
				// функция коректно работает только при наличии паддинга

			int z = (int)floor(X[2] * _1h_);
			int y = (int)floor(X[1] * _1h_);
			int x = (int)floor(X[0] * _1h_);

			_Point3 X__(X[0] - x * h_, X[1] - y * h_, X[2] - z * h_);
			_Point3 alpha(X__[0] * _1h_, X__[1] * _1h_, X__[2] * _1h_);
				// приводим к единичной длине

			z = mod(z, (int)sz_[2]);
			y = mod(y, (int)sz_[1]);
			x = mod(x, (int)sz_[0]);
				// координаты записи в линейный массив
			int y1 = (y == (int)sz_[1] - 1) ? 0 : y + 1;
			int x1 = (x == (int)sz_[0] - 1) ? 0 : x + 1;
				// координаты записи + 1 в линейный массив
				// z1 отсутствует, так как запись идет с учетом паддинга

			unsigned sz__ = sz_[2] + pad_ + pad_ - 1;
				// реальная длина ряда по z направлению
			unsigned z__ = z + pad_ - 1;
				// реальная смещение по z направлению

			_Real p, q, w; unsigned pos;

			p = s * (1. - alpha[0]);  q = p * (1. - alpha[1]);  w = q * alpha[2];
			pos = (x * sz_[1] + y) * sz__ + z__;
			fn[pos] = q - w; fn[pos + 1] = w;
				// запись в (x, y)

			q = p * alpha[1]; w = q * alpha[2];
			pos = (x * sz_[1] + y1) * sz__ + z__;
			fn[pos] = q - w; fn[pos + 1] = w;
				// запись в (x, y+1)

			p = s * alpha[0]; q = p * (1. - alpha[1]); w = q * alpha[2];
			pos = (x1 * sz_[1] + y) * sz__ + z__;
			fn[pos] = q - w; fn[pos + 1] = w;
				// запись в (x+1, y)

			q = p * alpha[1]; w = q * alpha[2];
			pos = (x1 * sz_[1] + y1) * sz__ + z__;
			fn[pos] = q - w; fn[pos + 1] = w;
				// запись в (x+1, y+1)
		}

		/**
		*   Вызывается после вставки всех точек в сетку и выполняет перенос данных,
		*   попавших в паддинг область, "внутрь" сетки.
		* @param fn сетка, куда вставляется значение
		*/
		void finalize(_Real *fn)
		{
			if (pad_ == 0) return;

			unsigned k, pos = 0, offset = sz_[2] + 2 * pad_ - 1, off = sz_[2];
			for (unsigned x=0; x<sz_[0]; x++)
			for (unsigned y=0; y<sz_[1]; y++)
			{
				for (int z=0; z<(int)pad_-1; z++)
				{
					k = pos + z;
					fn[k + off] += fn[k];
					fn[k] = 0.; // очистка для правильной следующей вставки
				}
				for (int z=0; z<(int)pad_; z++)
				{
					k = pos + pad_ - 1 + z;
					fn[k] += fn[k + off];
					fn[k + off] =0.; // очистка для правильной следующей вставки
				}
				pos += offset;
			}
		}
	};

}
#endif

