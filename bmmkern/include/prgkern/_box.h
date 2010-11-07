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
#ifndef _BOX__0077A726_9B21_5114_65E1_C84410790000__H
#define _BOX__0077A726_9B21_5114_65E1_C84410790000__H

/** @file
 *  Файл содержит описание объекта Box (ящик), который является прямоугольной
 *  областью пространства.
 *
 *  Ящики нужны для ограничения областей простанства для различных алгоритмов.
 *  Наиболее часто ящики используются для алгоритмов поиска ближайших соседей.
 */

#include "prgkern/_prgconfig.h"
#include "prgkern/_math.h"
#include "prgkern/_type.h"
#include "prgkern/_v3dense.h"
#include "prgkern/_for.h"

namespace prgkern
{
	/**
	* @brief Шаблон для описания ящика любой размерности.
	* @note У ящика все ребра параллельны осям координат. Это условие дает
	*  максимально простой и эффективный код.
	* @param N размерность пространства
	* @param T тип вещественного числа
	*/
	template <unsigned N, typename T>
	class Box_
	{
		typedef vdense_<N, T>  _Vector;

	public:
		typedef _Vector  point_type;

		/**
		* @note По умалчиванию создается инвертированный ящик. В нем нижняя граница
		*  больше верхней. Это позволяет подставлять такой ящик как стартовый
		*  объект при поиске максимального/минимального ящика. См. EMPTY ящик.
		*/
		Box_() { T inf = infinity<T>(); d_ = inf; u_ = -inf; }

		/**
		*   Конструктор одномерного ящика.
		* @param x0 левая граница
		* @param x1 правая граница
		*/
		Box_(T x0, T x1) { d_[0] = x0; u_[0] = x1; }

		/**
		*   Конструктор двумерного ящика.
		* @param x0 левая граница x-координаты
		* @param y0 левая граница y-координаты
		* @param x1 правая граница x-координаты
		* @param y1 правая граница y-координаты
		*/
		Box_(T x0, T y0, T x1, T y1)
		{
			d_[0] = x0; d_[1] = y0;
			u_[0] = x1; u_[1] = y1;
		}

		/**
		* @param x0 левая граница x-координаты
		* @param y0 левая граница y-координаты
		* @param z0 левая граница z-координаты
		* @param x1 правая граница x-координаты
		* @param y1 правая граница y-координаты
		* @param z1 правая граница z-координаты
		*/
		Box_(T x0, T y0, T z0, T x1, T y1, T z1)
		{
			d_[0] = x0; d_[1] = y0; d_[2] = z0;
			u_[0] = x1; u_[1] = y1; u_[2] = z1;
		}

		/**
		*   Конструктор ящика любой размерности.
		* @param d левая (нижняя) граница
		* @param u правая (верхняя) граница
		*/
		Box_(const _Vector &d, const _Vector &u) : d_(d), u_(u) {}

		/**
		*  Обменивает размеры ящиков друг с другом
		* @param b другой ящик
		*/
		void swap(Box_ b)
		{
			prgkern::swap(d_, b.d_);
			prgkern::swap(u_, b.u_);
		}

		/**
		* Обеспечивает доступ к границе ящика
		* @return левая нижняя точка ящика
		*/
		const _Vector &bottom() const { return d_; }

		/**
		* Обеспечивает доступ на изменение к границе ящика
		* @return левая нижняя точка ящика
		*/
		_Vector &bottom() { return d_; }

		/**
		* Обеспечивает доступ к границе ящика
		* @return правая верхняя точка ящика
		*/
		const _Vector &top() const { return u_; }

		/**
		* Обеспечивает доступ на изменение к границе ящика
		* @return правая верхняя точка ящика
		*/
		_Vector &top() { return u_; }

		/**
		*   Рассчитывает объем ящика с контролем на его недействительность.
		*/
		T volume() const
		{
			T s = (T) 1; bool sign = true;
			for (unsigned i=0; i<N; i++)
			{
				s *= u_[i] - d_[i];
				sign = sign && (u_[i] > d_[i]);
			}
			return sign ? s : (T) 0.;
		}

		/**
		*  Возвращает длину любого ребра
		* @param i номер ребра
		* @return длина ребра i
		*/
		T length(unsigned i) const { assert(_LT(i, N)); return u_[i] - d_[i]; }

		/**
		*  Возвращает вектор длин ребер
		* @return вектор, содержащий все длины ребер
		*/
		_Vector length() const { return u_ - d_; }

		/**
		*  Расшираяет ящик, чтобы включить в него заданный ящик.
		* @param b - заданный включаемый ящик
		* @return текущий объединенный ящик
		*/
		Box_ &make_union(const Box_ &b)
		{
			bottom () = min(bottom (), b.bottom ());
			top() = max(top(), b.top());
			return *this;
		}

		/**
		*  Расшираяет ящик, чтобы включить в него заданную точку.
		* @param X координаты точки
		* @return текущий объединенный ящик
		*/
		Box_ &make_union(const _Vector &X)
		{
			bottom() = min(bottom (), X);
			top() = max(top(), X);
			return *this;
		}

		/**
		* Увеличивает ящик относительно его центра.
		* @param coef - коэффициент увеличения (1 - не меняется)
		* @return текущий увеличенный ящик
		*/
		Box_ &enlarge(T coef)
		{
			_Vector p = length() * (0.5 * (coef - 1.));
			d_ -= p; u_+= p;
			return *this;
		}

		/**
		* Увеличивает ящик относительно заданной точки X. Если
		* точка вне ящика, то он сдвигается от нее (побочный эффект).
		* @param coef - коэффициент увеличения (1 - не меняется)
		* @param X - центр преобразования
		* @return текущий увеличенный ящик (и сдвинутый)
		*/
		Box_ &enlarge(T coef, const _Vector &X)
		{
			_Vector p = X * (1. - coef);
			d_ *= coef;
			u_ *= coef;
			d_ += p;
			u_ += p;
			return *this;
		}

		/**
		* Увеличивает ящик на заданную длину X.
		* @param X - центр преобразования
		* @return текущий увеличенный ящик (и сдвинутый)
		*/
		Box_ &enlarge(const _Vector &X)
		{
			d_ -= X;
			u_ += X;
			return *this;
		}

		/**
		* Сдвигает ящик на заданную длину X.
		* @param X вектор сдвига
		* @return текущий сдвинутый ящик
		*/
		Box_ &move(const _Vector &X)
		{
			d_ += X;
			u_ += X;
			return *this;
		}

		/**
		* Сдвигает ящик так, чтобы центр совпал с заданной точкой.
		* @param X новое положение центра ящика
		* @return текущий сдвинутый ящик
		*/
		Box_ &moveto(const _Vector &X)
		{
			_Vector X__ = X - (u_ + d_) * (T)0.5;
			d_ += X__;
			u_ += X__;
			return *this;
		}

		/**
		* Проверяет, входит ли заданная точка в ящик.
		* @param X координаты точки
		* @return true/false входит/не всходит
		*/
		bool is_included(const _Vector &X) const
		{
			if (d_ <= X && X < u_) return true;
			return false;
		}

		/**
		* Возвращает центр ящика.
		*/
		_Vector center() const { return (d_ + u_) * (T)0.5; }

	protected:
		_Vector d_; ///< левый нижний край
		_Vector u_; ///< правый верхний край
	};

	/**
	*  Объединение двух ящиков произвольной но единой размерности.
	*  При этом включается не только оба ящика, но и пространство между ними.
	* @param box - результат объединения
	* @param b1 - 1-й ящик
	* @param b2 - 2-й ящик
	*/
	template <unsigned N, typename T>
	INLINE void get_union(Box_<N, T> &box, const Box_<N, T> &b1, const Box_<N, T> &b2)
	{
		box.bottom() = min(b1.bottom(), b2.bottom());
		box.top() = max(b1.top(), b2.top());
	}

	/**
	*  Обмен двух ящиков произвольной но единой размерности.
	* @param b1 - 1-й ящик
	* @param b2 - 2-й ящик
	*/
	template <unsigned N, typename T>
	INLINE void swap(Box_<N, T> &b1, Box_<N, T> &b2)
	{
		swap(b1.bottom(), b2.bottom());
		swap(b1.top(), b2.top());
	}

	/**
	* Строковое представление ящика
	* @param box - ящик
	* @return строка информации о координатах ящика.
	*/
	template <unsigned N, typename T>
	INLINE std::string make_string(const Box_<N, T> &box)
	{
		return _S("[") + make_string(box.bottom()) + _S(",")
			+ make_string(box.top()) + _S("]");
	}

	/**
	*  Вывод строкового представления ящика в стандартный выводной поток
	* @param os - не требуется
	* @param box - ящик
	*/
	template <unsigned N, typename T>
	INLINE std::ostream &operator<<(std::ostream &os, const Box_<N, T> &box)
	{ os << make_string(box); return os; }

}
#endif
