#ifndef _PROTONIZATION__0077A726_81CC_5dce_36BD_A3453EED0200__H
#define _PROTONIZATION__0077A726_81CC_5dce_36BD_A3453EED0200__H

#include "molkern/__moldefs.h"

/*==============================================================================
 *                 ОПИСАНИЕ ОСОБЕННОСТЕЙ РЕАЛИЗАЦИИ
 *
 *    В файле помещаются объекты, которые позволяют восстанавливать геометрию
 *  отсутствующих атомов водорода. В реализации предполагается, что все
 *  координаты атомов имеют тип real_t, независимо от того, являются ли
 *  координаты атомов рассчитываемыми внутри функций, или приходят как
 *  параметры. Это делается специально, чтобы избежать непреднамеренных
 *  преобразований данных и акцентировать внимание программиста на необходимость
 *  для любых координат использования заданного типа.
 =============================================================================*/

namespace molkern
{

	/**
	* @brief builds the position of one atom using coordinates of atom group
	* @param XX[out] array of added atoms (must be 1..3)
	* @param XA array of exist atoms (may be 1..3)
	* @param X0 center atom (new atoms are added to it)
	* @param r0 equalibrium distance for new atoms (from X0)
	*/
	// Это ранний вариант, где водород устанавливается вдоль прямой. Это неверно делать.
	//	INLINE void build_1A_1H_(vector_t *XH, const vector_t *XA,
	//		const vector_t &X0, real_t r0=HYDROGEN_ATOM_DISTANCE)
	//	{
	//		vector_t AX0 = X0 - XA[0];
	//		AX0.normalize(r0);
	//		XH[0] = X0 + AX0;
	//	}

	INLINE void build_1A_1H_(vector_t *XH, const vector_t *XA,
		const vector_t &X0, real_t r0=HYDROGEN_ATOM_DISTANCE)
	{
		vector_t AX0 = X0 - XA[0];
		vector_t AX1 = AX0 - vector_t(ran0(-1., 1.), ran0(-1., 1.), ran0(-1., 1.));

		AX0.normalize();
		AX1 = AX1 - AX0 * scalar_product(AX1, AX0);
		AX1.normalize();

		// В данной функции атом строится таким образом, чтобы образовывать угол 120 градусов
		// между атомами A-A-H. Такой угол практически между любыми атомами (кроме HZ).
		// При установке угла 180 градусов часто возникают nan, так как это вырожденный угол.
		XH[0] = X0 + AX0 * M_COS_60 * r0 + AX1 * M_SIN_60 * r0;
	}

	/**
	* @brief builds the position of two atoms using coordinates of atom group
	* @param XX[out] array of added atoms (must be 1..3)
	* @param XA array of exist atoms (may be 1..3)
	* @param X0 center atom (new atoms are added to it)
	* @param req equalibrium distance for new atoms (from X0)
	*/
	INLINE void build_1A_2H_(vector_t *XH, const vector_t *XA, const vector_t &X0,
		real_t r0=HYDROGEN_ATOM_DISTANCE)
	{
		vector_t AX0 = X0 - XA[0];
		AX0.normalize((real_t) 0.5 * r0);

		vector_t ort = ortogonal(AX0);
		ort.normalize((real_t) M_SIN_60 * r0);

		AX0 = X0 + AX0;
		XH[0] = AX0 + ort;
		XH[1] = AX0 - ort;
	}

	/**
	* @brief builds the position of 3 atoms using coordinates of atom group
	* @param XX[out] array of added atoms (must be 1..3)
	* @param XA array of exist atoms (may be 1..3)
	* @param X0 center atom (new atoms are added to it)
	* @param req equalibrium distance for new atoms (from X0)
	*/
	INLINE void build_1A_3H_(vector_t *XH, const vector_t *XA, const vector_t &X0,
		real_t r0=HYDROGEN_ATOM_DISTANCE)
	{
		// look to mathworld.wolfram.com/tetrahedron.html
		// vertexs are located in (x, 0, 0), (-d, (+|-) a/2, 0), (0, 0, h)
		// x = a / Sqrt(3), d = 0.5 * x, h = a * Sqrt(2)/ Sqrt(3)
		// circumradius R = a * Sqrt(6) / 4, ( circumradius == HYDROGEN_ATOM_DISTANCE)
		// intraradius r = a * Sqrt(6) / 12

		vector_t AX0 = X0 - XA[0];
		real_t a = (real_t) (4 * r0 / M_SQRT_6); // length of tetrahedron edge
		AX0.normalize((real_t) (a * M_SQRT_6 / 12)); // normalize to intraR

		real_t x = (real_t) (a / M_SQRT_3);
		vector_t ort_x = ortogonal(AX0);
		ort_x.normalize(x); // normalize to x
		vector_t ort_a= vector_product(AX0, ort_x);
		ort_a.normalize((real_t) 0.5 * a); // normalize to a

		AX0 = X0 + AX0;
		XH[0] = AX0 + ort_x;
		ort_x *= (real_t) 0.5;
		AX0 -= ort_x;
		XH[1] = AX0 + ort_a;
		XH[2] = AX0 - ort_a;
	}

	/**
	* @brief builds the position of 1 atom using coordinates of atom group
	* @param XX[out] array of added atoms (must be 1..3)
	* @param XA array of exist atoms (may be 1..3)
	* @param X0 center atom (new atoms are added to it)
	* @param req equalibrium distance for new atoms (from X0)
	*/
	INLINE void build_2A_1H_(vector_t *XH, const vector_t *XA, const vector_t &X0,
		real_t r0=HYDROGEN_ATOM_DISTANCE)
	{
		vector_t AX0 = X0 - XA[0];
		AX0.normalize();
		vector_t AX0__ = X0 - XA[1];
		AX0__.normalize();
		AX0 += AX0__;
		AX0.normalize(r0);
		XH[0] = X0 + AX0;
	}

	/**
	* @brief builds the position of 2 atom using coordinates of atom group
	* @param XX[out] array of added atoms (must be 1..3)
	* @param XA array of exist atoms (may be 1..3)
	* @param X0 center atom (new atoms are added to it)
	* @param req equalibrium distance for new atoms (from X0)
	*/
	INLINE void build_2A_2H_(vector_t *XH, const vector_t *XA, const vector_t &X0,
		real_t r0=HYDROGEN_ATOM_DISTANCE)
	{
		vector_t AX0 = X0 - XA[0];
		AX0.normalize();
		vector_t AX0__ = X0 - XA[1];
		AX0__.normalize();
		vector_t ort = vector_product(AX0, AX0__);
		ort.normalize((real_t) M_SIN_60 * r0);
		AX0 += AX0__;
		AX0.normalize((real_t) 0.5 * r0);

		AX0 = X0 + AX0;
		XH[0] = AX0 + ort;
		XH[1] = AX0 - ort;
	}

	/**
	* @brief builds the position of 2 atom using coordinates of atom group
	* @param XX[out] array of added atoms (must be 1..3)
	* @param XA array of exist atoms (may be 1..3)
	* @param X0 center atom (new atoms are added to it)
	* @param req equalibrium distance for new atoms (from X0)
	*/
	INLINE void build_3A_1H_(vector_t *XH, const vector_t *XA, const vector_t &X0,
		real_t r0=HYDROGEN_ATOM_DISTANCE)
	{
		vector_t A = X0 - XA[0];
		A.normalize();
		vector_t B = X0 - XA[1];
		B.normalize();
		vector_t C = X0 - XA[2];
		C.normalize();

		A += B + C;
		A.normalize(r0);
		XH[0] = X0 + A;
	}

	/**
	* @brief построение позиции атома по образцу
	* @param XA[in,out] массив существующих (3) и добавляемых (1) атомов
	* @param XS образец для сравнения (4 атома)
	* @note порядок атомов в массивах XA и XS должен совпадать и первым должен
	*  располагаться цетральный атом, с которым связан восстанавливаемый атом,
	*  далее должен идти атом, определяющий ось ротамера, чтобы улучшить точность.
	*  Например, при восстановлении атома OXT таким атомы должны идти в порядке
	*  C(цетральный), CA(ось ротамера), O(вращающийся вокруш ротамерной оси).
	*/
	INLINE void build_2A_1A_(vector_t *XA, const vector_t *XS)
	{
		//--------------------------------------------------------------------------
		//                        совмещение центров
		//--------------------------------------------------------------------------
		vector_t OO = XA[0] - XS[0];
		vector_t XS__[4] = { XA[0], XS[1] + OO, XS[2] + OO, XS[3] + OO };

		//--------------------------------------------------------------------------
		//                  точное совмещение 01 и 01' осей
		//--------------------------------------------------------------------------
		vector_t A1 = XA[1] - XA[0], B1 = XS__[1] - XS__[0];
		real_t w = get_angle(B1, A1); // угол между осями 01 и 01'
		vector_t W = vector_product(B1, A1); // ось вращения, чтобы совместить 01 и 01'
		{
			Rotator<AXIS_ROTATOR_, real_t> rotator(W, w, XA[0]);
			XS__[1] = rotator(XS__[1]);
			XS__[2] = rotator(XS__[2]);
			XS__[3] = rotator(XS__[3]);
		}

		//--------------------------------------------------------------------------
		//             максимально близкое совмещение 02 и 02' осей
		//--------------------------------------------------------------------------
		vector_t A2 = XA[2] - XA[0], B2 = XS__[2] - XS__[0];
		vector_t Q1 = vector_product(A1, A2);
		vector_t Q2 = vector_product(B1, B2);
		w = get_angle(Q1, Q2); // угол между плоскостями 102 и 102'
		if (scalar_product(A1, vector_product(B2, A2)) < 0.) w = -w;
			// определение направления вращения
		{
			Rotator<AXIS_ROTATOR_, real_t> rotator(A1, w, XA[0]);
			XS__[2] = rotator(XS__[2]);
			XS__[3] = rotator(XS__[3]);
		}

		//--------------------------------------------------------------------------
		//             определение позиции добавляемого атома
		//--------------------------------------------------------------------------
		XA[3] = XS__[3];
	}

}
#endif

