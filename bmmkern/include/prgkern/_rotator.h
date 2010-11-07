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
#ifndef _ROTATOR__BFC42EF6_195A_5a1e_6CF1_6F41AEE10300__H
#define _ROTATOR__BFC42EF6_195A_5a1e_6CF1_6F41AEE10300__H
#include "prgkern/_prgconfig.h"

#include "prgkern/_type.h"
#include "prgkern/_sse.h"
#include "prgkern/_v3dense.h"
#include "prgkern/_m3x3dense.h"

namespace prgkern
{

	/**
	* @brief Prepares 3x3 rotate matrix R(alpha, beta, gamma)
	*        = R(gamma) * R(beta) * R(alpha)
	* @param m - R(alpha, beta, gamma) matrix
	* @param angle - array of angles ([0]alpha, [1]beta, [2]gamma).
	*/
	template <typename T>
	inline void make_rotate0_matrix(mdense_<3, 3, T> *m, const vdense_<3, T> &angle)
	{
		if (m == NULL) return;

		// | cos(g)  sin(g)  0 |   | 1     0       0    |   | cos(a)  sin(a) 0 |
		// | -sin(g) cos(g)  0 | * | 0   cos(b)  sin(b) | * | -sin(a) cos(a) 0 | =
		// |   0       0     1 |   | 0  -sin(b)  cos(b) |   |   0       0    1 |
		//
		// | cg ca - cb sa sg    cg sa + cb ca sg    sg sb |
		// | -sg ca - cb sa cg   -sg sa + cb ca cg   cg sb |
		// |      sb sa               -sb ca         cb    |

		T ca = (T) cos(angle[0]); T sa = (T) sin(angle[0]);
		T cb = (T) cos(angle[1]); T sb = (T) sin(angle[1]);
		T cg = (T) cos(angle[2]); T sg = (T) sin(angle[2]);

		T sasg = sa * sg;
		T cacg = ca * cg;
		T casg = ca * sg;
		T sacg = sa * cg;

		(*m)[0][0] =  cacg - cb*sasg; (*m)[0][1] =  sacg + cb*casg; (*m)[0][2] = sg*sb;
		(*m)[1][0] = -casg - cb*sacg; (*m)[1][1] = -sasg + cb*cacg; (*m)[1][2] = cg*sb;
		(*m)[2][0] =  sb*sa;          (*m)[2][1] = -sb*ca;          (*m)[2][2] = cb;
	}

	const int EULER_ROTATOR_    = 1; // euler R(alpha, beta, gamma) rotator
	const int AXIS_ROTATOR_     = 2; // axis R(phi) rotator
	const int INFINITE_ROTATOR_ = 3; // moving rotator (R = infinity)
	const int XYZ_ROTATOR_      = 4;

	const int NORMAL_AXIS_ROTATOR_ = 5;

	/**
	* @brief Template of rotator class
	* @param ROTATOR_TYPE - ident of rotator type
	*/
	template <int ROTATOR_TYPE, typename T> class Rotator;

	/**
	* @brief Implementation of EULER_ROTATOR_TYPE rotator
	*/
	template <typename T> class Rotator<INFINITE_ROTATOR_, T>
	{
	public:
		/**
		* @brief constructs object and prepares the internal data
		* @note  X - moving the coordinate system
		*       -X - moving the atoms
		* @param X0 - moving vector
		*/
		Rotator(const vdense_<3, T> &X0) : X_(X0) {}

		/**
		* @brief move vector
		* @param R - result of moving
		* @param X - vector to move
		*/
		void operator()(vdense_<3, T> &R, const vdense_<3, T> &X) const { R = X - X_; }

	protected:

		vdense_<3, T> X_;
	};

	/**
	* @brief Implementation of EULER_ROTATOR_ rotator
	*/
	template <typename T> class Rotator<EULER_ROTATOR_, T>
	{
		vdense_<3, T> X0_;
		mdense_<3, 3, T> M_;

	public:

		/**
		 * Конструирование и приготовление внутренних данных
		 * @param X0 центр вращения
		 * @param angles (+++) углы вращения координатной системы, (---) - углы вращения атомов
		 * @return
		 */
		Rotator(const vdense_<3, T> &angles, const vdense_<3, T> &X0) : X0_(X0)
		{ make_rotate0_matrix(&M_, angles); }

		/**
		* @brief rotate vector onto (alpha, beta, gamma)
		* @param X - vector to rotate
		*/
		vdense_<3, T> operator()(vdense_<3, T> X) const { return M_ * (X - X0_) + X0_; }
	};

	/**
	* @brief Implementation of AXIS_ROTATOR_ rotator R(phi)
	*   R = X * cosa + RV * (1 - cosa) * scalar_product(RV, X)
	*     + vector_product(RV, X) * sina
	*/
	template <typename T> class Rotator<AXIS_ROTATOR_, T>
	{
	public:
		/**
		* @param W ось вращения (нормализуется внутри класса)
		* @param w угол вращения
		* @param O точка опоры для вращения
		*/
		Rotator(const vdense_<3, T> &W, T w, const vdense_<3, T> &O)
			: ca_((T)cos(w)), sa_((T)sin(w)), W_(W), O_(O) { W_.normalize(); }

		/**
		* @brief rotate vector about axis by phi angle
		* @note function uses previous calculated data
		* @param R - result of rotation
		* @param X - vector to rotate
		*/
		vdense_<3, T> operator()(const vdense_<3, T> &X) const
		{
			vdense_<3, T> X0 = X - O_;
			vdense_<3, T> XW = W_ * scalar_product(W_, X0); // невращаемая компонента вдоль оси W
			vdense_<3, T> XP = X0 - XW; // вращаемая перпендикулярная компонента
			vdense_<3, T> XT = vector_product(W_, X0);
			return XW + XP * ca_ + XT * sa_ + O_;
		}

	protected:
		T ca_, sa_;
		vdense_<3, T> W_, O_;
	};


	template <typename T> class Rotator<XYZ_ROTATOR_, T>
	{
		Rotator<AXIS_ROTATOR_, T> x_rotator_;
		Rotator<AXIS_ROTATOR_, T> y_rotator_;
		Rotator<AXIS_ROTATOR_, T> z_rotator_;

	public:
		/**
		* @param W ось вращения (нормализуется внутри класса)
		* @param w угол вращения
		* @param X0 точка опоры для вращения
		*/
		Rotator(const vdense_<3, T> &angle, const vdense_<3, T> &X0)
		: x_rotator_(vdense_<3, T>(1., 0., 0.), angle[0], X0),
		  y_rotator_(vdense_<3, T>(0., 1., 0.), angle[1], X0),
		  z_rotator_(vdense_<3, T>(0., 0., 1.), angle[2], X0)
		{}

		/**
		* @brief rotate vector about axis by phi angle
		* @note function uses previous calculated data
		* @param R - result of rotation
		* @param X - vector to rotate
		*/
		vdense_<3, T> operator()(const vdense_<3, T> &X) const
		{ return z_rotator_(y_rotator_(x_rotator_(X))); }

	};

	template <typename T> class Rotator<NORMAL_AXIS_ROTATOR_, T>
	{
		vdense_<3, T> O_, W_;
		T ca_, sa_;
		vdense_<3, T> A_, B_;

	public:

		/**
		 * Инициализация ротатора.
		 * (!) Все проверки на допустимость вращения должны быть выполнены ранее, иначе использование
		 * становится неэффективным из-за проверок, применяемых к множеству вращаемых атомов.
		 * (!) Поведение некорректно при использовании с ненормированными A & B.
		 * @param O - точка опоры для вращения
		 * @param A - вектор, к которому нужно привести
		 * @param B - вектор, который нужно привести
		 */
		Rotator(const vdense_<3, T> &O, const vdense_<3, T> &A, const vdense_<3, T> &B)
			: O_(O), W_(vector_product(B, A)),
			  ca_(scalar_product(B, A)), // проверка на sc == (+|-)1 делается вне
			  sa_(W_.normalize()), // проверка на W != 0 делается вне
			  A_(A), B_(vector_product(A, W_))
		{
		}

		/**
		 * Выполнение вращения любого вектора. (!) Проверка на позиционирование вращаемого вектора на
		 * оси вращения не выполняется, так как и для них формулы дают корректный результат.
		 * @param X - вектор
		 * @return результат вращения
		 */
		vdense_<3, T> operator()(const vdense_<3, T> &X) const
		{
			vdense_<3, T> XX = X - O_; // приведеный к точке O вектор
			vdense_<3, T> XW = W_ * scalar_product(W_, XX); // невращаемая компонента вдоль оси W
			vdense_<3, T> XA = A_ * scalar_product(A_, XX); // вращаемая компонента вдоль оси A
			vdense_<3, T> XB = B_ * scalar_product(B_, XX); // вращаемая компонента вдоль оси B

			return vdense_<3, T>(
				XW[0] + XB[0] * ca_ + XA[0] * sa_ + O_[0],
				XW[1] + XB[1] * ca_ + XA[1] * sa_ + O_[1],
				XW[2] + XB[2] * ca_ + XA[2] * sa_ + O_[2]
			);
		}

	};

}
#endif
