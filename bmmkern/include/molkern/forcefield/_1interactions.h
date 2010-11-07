#ifndef _1INTERACTION__F9ED1116_D4D6_5cf5_C20E_AC432B290700__H
#define _1INTERACTION__F9ED1116_D4D6_5cf5_C20E_AC432B290700__H

#include "molkern/__moldefs.h"

namespace molkern
{
	using namespace prgkern;

	template <int TYPE> struct Interaction1_;


	template <> struct Interaction1_<BOND_>
	{
		/**
		* @brief энергия связи
		* @param ke,q0 параметры подфункции
		* @param R вектор XB-XA
		*/
		static real_t U(real_t ke, real_t q0, const vector_t &R)
		{
			real_t q = (real_t) sqrt(scalar_product(R, R));
			return ke * sqr(q - q0);
		}

		/**
		* @brief расчет силы, действующих на атомы связи
		*   Вектор, определяющие расстояние, заменяются силой.
		* @param ke,q0 параметры подфункции
		* @param R[in] вектор XB-XA
		* @param R[out] сила, действующая на атомы
		* @return энергия
		*/
		static real_t dU__dX(real_t ke, real_t q0, vector_t &R)
		{
			real_t q = (real_t) sqrt(scalar_product(R, R));
			real_t q__ = q - q0;
			R *= (ke + ke) * q__ / q;
			return ke * sqr(q__);
		}
	};

	template <> struct Interaction1_<ANGLE_>
	{
		/**
		* @brief энергия валентных углов
		* @param ke силовая постоянная
		* @param q0 равновесный угол
		* @param A,B координаты векторов XA-XB и XC-XB
		*/
		static real_t U(real_t ke, real_t q0, const vector_t &A, const vector_t &B)
		{
			real_t a2 = scalar_product(A, A);
			real_t b2 = scalar_product(B, B);
			real_t ab = scalar_product(A, B);
			real_t q  = (real_t) safe_acos(ab / sqrt(a2 * b2));
			return ke * sqr(q - q0);
		}

		/**
		* @brief расчет сил, действующих на атомы валентных углов
		*   Вектора, определяющие угол, заменяются силами. Причем возвращаются
		*   силы только для крайних атомов. Для центрального атома F = -(FA + FB).
		* @param ke силовая постоянная
		* @param q0 равновесный угол
		* @param A,B[in] координаты векторов XA-XB и XC-XB
		* @param A,B[out] силы, действующие на краевые атомы A и B
		* @return энергия
		*/
		static real_t dU__dX(real_t ke, real_t q0, vector_t &A, vector_t &B)
		{
			const real_t high_accuracy = (real_t) 10. * std::numeric_limits<real_t>::epsilon();

			real_t a2 = scalar_product(A, A);
			real_t b2 = scalar_product(B, B);
			real_t ab = scalar_product(A, B);
			real_t ab2 = sqr(ab);

			real_t q = (real_t) safe_acos(ab / sqrt(a2 * b2));
			real_t energy = ke * sqr(q - q0);

			if ( fabs(ab2 - a2*b2) < high_accuracy)
			{
				A[0] = 0.; A[1] = 0.; A[2] = 0.;
				B[0] = 0.; B[1] = 0.; B[2] = 0.;
					// При вырождении следует передавать значение энергии != 0,
					// иначе происходит скачок энергии около вырожденного угла.
					// Все силы равны нулю в точке вырождения (неустойчивое равновесие).
				return energy;
			}
			real_t m = (ke + ke) * (q - q0) / (real_t) sqrt(a2 * b2 - ab2);
				// момент сил, поворачивающий атомы, деленный на расстояние

			real_t A_[3] = { A[0], A[1], A[2] };
			real_t B_[3] = { B[0], B[1], B[2] };

			real_t ab_a2 = ab / a2;
			A[0] = (B_[0] - A_[0] * ab_a2) * m;
			A[1] = (B_[1] - A_[1] * ab_a2) * m;
			A[2] = (B_[2] - A_[2] * ab_a2) * m;

			real_t ab_b2 = ab / b2;
			B[0] = (A_[0] - B_[0] * ab_b2) * m;
			B[1] = (A_[1] - B_[1] * ab_b2) * m;
			B[2] = (A_[2] - B_[2] * ab_b2) * m;
				// знак не соответствует определению функции (запись сил, а не производных)

			return energy;
		}
	};

	template <> struct Interaction1_<TORSION_>
	{
		/**
		* @brief энергия дигедральных углов
		* @param curf,v,phi,n параметры подфункций
		* @param XA,XB,XC,XD координаты векторов
		*/
		static real_t U(int curf, const int *n, const real_t *v, const real_t *phi,
			const vector_t &XA, const vector_t &XB, const vector_t &XC, const vector_t &XD)
		{
			vector_t A(XA[0] - XB[0], XA[1] - XB[1], XA[2] - XB[2]);
			vector_t B(XD[0] - XC[0], XD[1] - XC[1], XD[2] - XC[2]);
			vector_t C(XC[0] - XB[0], XC[1] - XB[1], XC[2] - XB[2]);

			vector_t U = vector_product(A, C);
			vector_t T = vector_product(B, C);

			real_t ut = (real_t) sqrt(scalar_product(U, U) * scalar_product(T, T));
			if (ut < ACCURACY) return 0;
				// избегание "раскрытых" углов в системах с H-C=C-M

			real_t q = safe_acos(scalar_product(U, T) / ut);

			_E(real_t) energy = 0;
			// Так как только положительные значения n[i] используются
			// отрицательные значения n[i] говорят о конце цикла.
			for (int i=0; i<curf; i++)
			{
				real_t s = q * n[i] - phi[i];
				energy += v[i] * (1 + (real_t)cos(s));
					// 1/2 коэффициент внутри V
			}
			return (real_t)energy;
		}

		/**
		* @brief расчет сил, действующих на атомы дигедральных углов
		*   Вектора, определяющие угол, заменяются силами. Возвращаются
		*   все силы, в отличие от расчета валентных углов.
		* @param curf число подфункций
		* @param v,phi,n параметры подфункций
		* @param XA,XB,XC,XD[in] координаты векторов
		* @param XA,XB,XC,XD[out] силы, действующие на атомы
		* @return энергия
		*/
		static real_t dU__dX(int curf, const int *n, const real_t *v, const real_t *phi,
			vector_t &XA, vector_t &XB, vector_t &XC, vector_t &XD)
		{
			// Предыдущий вариант этой функции был построен на прямом дифференцировании
			// энергии торсионного угла от смещения координат. Это оказалось неэффективно.
			// Текущий вариант построен на уравновешивании всех моментов вращения,
			// которые образуются при угле, не равном равновесному. Что оказалось
			// заметно проще и эффективнее. Поскольку у меня нет проработанной математики
			// из которой нужные формулы выводятся, то данный алгоритм был построен на
			// основе наблюдения за тем, как распределяются моменты в первоначальной
			// версии алгоритма.

			const real_t high_accuracy = (real_t) 10. * std::numeric_limits<real_t>::epsilon();

			vector_t A(XA[0] - XB[0], XA[1] - XB[1], XA[2] - XB[2]);
			vector_t B(XD[0] - XC[0], XD[1] - XC[1], XD[2] - XC[2]);
			vector_t C(XC[0] - XB[0], XC[1] - XB[1], XC[2] - XB[2]);

			real_t a2 = scalar_product(A, A);
			real_t c2 = scalar_product(C, C);
			real_t b2 = scalar_product(B, B);
			real_t ac = scalar_product(A, C);
			real_t ac2 = sqr(ac);
			real_t bc = -scalar_product(B, C);
			real_t bc2 = sqr(bc);

			vector_t U = vector_product(A, C);
			vector_t T = vector_product(B, C);
			real_t u2 = scalar_product(U, U);
			real_t t2 = scalar_product(T, T);
			real_t ut = scalar_product(U, T);
			if (u2 < ACCURACY || t2 < ACCURACY) return 0;
				// избегание "раскрытых" углов в системах с H-C=C-M

			real_t q = safe_acos(ut / (real_t) sqrt(u2 * t2));
			real_t energy = 0.;
			real_t m = 0.; // момент вращения
			for (int i=0; i<curf; i++)
			{
				real_t s = q * n[i] - phi[i];
				energy += v[i] * ((real_t) 1. + (real_t) cos(s));
				m -= v[i] * n[i] * (real_t) sin(s);
					// коэффициент 1/2 включен в params_.v
			}
			if ( fabs(ac2 - a2*c2) < high_accuracy
				|| fabs(bc2 - b2*c2) < high_accuracy
				|| fabs(sqr(ut) - u2*t2) < high_accuracy
			)
			{
				XA[0] = 0.; XA[1] = 0.; XA[2] = 0.;
				XB[0] = 0.; XB[1] = 0.; XB[2] = 0.;
				XC[0] = 0.; XC[1] = 0.; XC[2] = 0.;
				XD[0] = 0.; XD[1] = 0.; XD[2] = 0.;

				return energy;
				// При вырождении следует передавать значение энергии != 0,
				// иначе происходит скачок энергии около вырожденного угла.
				// Все силы равны нулю в точке вырождения (неустойчивое равновесие).
			}

			vector_t UxT = vector_product(U, T);
			if (scalar_product(C, UxT) > 0) m = -m;
				// отслеживаем направление сил, ои должны идти в том направлении,
				// чтобы сводить угол к равновесному

			real_t __u = -m * (real_t) sqrt(c2 / (u2 * (a2 * c2 - ac2)));
			real_t __t =  m * (real_t) sqrt(c2 / (t2 * (b2 * c2 - bc2)));
			vector_t FA(U[0] * __u, U[1] * __u, U[2] * __u);
			vector_t FB(T[0] * __t, T[1] * __t, T[2] * __t);

			real_t __1_c2 = (real_t) 1. / c2;
			real_t __ac_c2 = ac * __1_c2;
			real_t __bc_c2 = bc * __1_c2;
			real_t s = (real_t) 1. - __ac_c2;
			real_t p = (real_t) 1. - __bc_c2;

			XA[0] = -FA[0]; XA[1] = -FA[1]; XA[2] = -FA[2];

			XB[0] = FA[0] * s + FB[0] * __bc_c2;
			XB[1] = FA[1] * s + FB[1] * __bc_c2;
			XB[2] = FA[2] * s + FB[2] * __bc_c2;

			XC[0] = FB[0] * p + FA[0] * __ac_c2;
			XC[1] = FB[1] * p + FA[1] * __ac_c2;
			XC[2] = FB[2] * p + FA[2] * __ac_c2;

			XD[0] = -FB[0]; XD[1] = -FB[1]; XD[2] = -FB[2];
				// знак не соответствует определению функции (запись сил, а не производных)
			return energy;
		}
	};

}
#endif
