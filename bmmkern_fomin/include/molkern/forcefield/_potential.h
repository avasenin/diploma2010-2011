#ifndef _POTENTIAL__0077A726_C0DA_57a5_2C2E_29442B440B00__H
#define _POTENTIAL__0077A726_C0DA_57a5_2C2E_29442B440B00__H

#include "molkern/__moldefs.h"
#include "molkern/forcefield/_atom.h"

namespace molkern
{
	using namespace prgkern;
	
	/// радиус разделения потенциала на ближнюю и дальнюю части
	const double DEFAULT_COULOMB_SPLIT_RADIUS = 8.;
	
	/**
	* @brief шаблон потенциалов разного типа
	* @param _Type - тип потенциала (COULOMB, VDW, ..)
	* @param _Real - тип вещественного числа
	*/
	template <int _Type, typename _Real=double> class Potential;
	
	/**
	*  Данный потенциал имеет специфичные особенности, не имеющие общего характера.
	*  К ним относится способ сглаживания особенности в 0 (с помощью erf функции)
	*  и связанный только с этим потенциалом способ расчета энергии самодействия.
	*/
	template <typename _Real> class Potential<COUL_, _Real>
	{
		_Real alpha_; ///< параметр, связанный с радиусом разделения потенциала 
		              ///< на ближнюю и дальнюю части
		
	public:
		
		/**
		*   Конструктор объекта.
		* @param split_radius - расстояние разделения потенциала на ближнюю и дальнюю части.
		*  В области (0, split_radius) из дальнего потенциала вырезана кулоновская особенность
		*  в 0 и потенциал близок к константе. В области (split_radius, inf) потенциал
		*  асимптотически приближается к 1/r.
		*/
		Potential(_Real split_radius=DEFAULT_COULOMB_SPLIT_RADIUS)
			: alpha_( M_SQRT_PI / split_radius ) {}
		
		/**
		* Данный параметр необходим для расчета энергии самодействия.
		* @return параметр в функции сглаживания erf(alpha * x)
		*/
		_Real alpha() const { return alpha_; }
		
		/**
		* Рассчитывает ближнюю часть потенциала на заданном расстоянии.
		* Потенциал определен только для r > M_SAFETY
		* @param r расстояние до источника
		* @return значение потенциала на заданном расстоянии
		*/
		_Real near(_Real r) const
		{
			assert(_GT(r, (_Real)M_SAFETY)); // потенциал определен для r >= M_SAFETY
			return (1. - erf(r * alpha_)) / r;
		}
		
		/**
		* Рассчитывает дальнюю часть потенциала на заданном расстоянии.
		* Потенциал определен только для r > M_SAFETY
		* @param r расстояние до источника
		* @return значение потенциала на заданном расстоянии
		*/
		_Real far(_Real r) const
		{
			if (r < M_SAFETY) return _FMULT<2>(alpha_) / M_SQRT_PI;
			return erf(r * alpha_) / r;
		}
		
		/**
		* Рассчитывает полный потенциал на заданном расстоянии.
		* Потенциал определен только для r > M_SAFETY
		* @param r расстояние до источника
		* @return значение потенциала на заданном расстоянии
		*/
		_Real full(_Real r) const
		{
			assert(_GT(r, (_Real)M_SAFETY)); // потенциал определен для r >= M_SAFETY
			return 1. / r;
		}
		
		/**
		* Рассчитывает производную ближней части потенциала на заданном расстоянии.
		* Потенциал определен только для r > M_SAFETY
		* @param r расстояние до источника
		* @return значение производной потенциала на заданном расстоянии
		*/
		_Real near1(_Real r) const
		{
			assert(_GT(r, (_Real)M_SAFETY)); // потенциал определен для r >= M_SAFETY
			
			_Real coef = _FMULT<2>(alpha_) / M_SQRT_PI;
			_Real s = coef * exp(-sqr(alpha_ * r));
			return ((erf(r * alpha_) - 1.) / r - s) / r;
		}
		
		/**
		* Рассчитывает производную дальней части потенциала на заданном расстоянии.
		* Потенциал определен и для r = 0, что необходимо для учета самодействия.
		* @param r расстояние до источника
		* @return значение производной потенциала на заданном расстоянии
		*/
		_Real far1(_Real r) const
		{
			if (r < M_SAFETY) return 0.;
			
			_Real coef = _FMULT<2>(alpha_) / M_SQRT_PI;
			_Real s = coef * exp(-sqr(alpha_ * r));
			return (s - erf(r * alpha_) / r) / r;
		}
		
		/**
		* Рассчитывает производную полного потенциала на заданном расстоянии.
		* Потенциал определен только для r > M_SAFETY
		* @param r расстояние до источника
		* @return значение производной потенциала на заданном расстоянии
		*/
		_Real full1(_Real r) const
		{
			assert(_GT(r, (_Real)M_SAFETY)); // потенциал определен для r >= M_SAFETY
			return -1. / sqr(r);
		}
		
		/**
		* Вычисляет фурье образ дальней части потенциала для 3D-случая.
		* @note Результат неопределен для k = 0, поскольку в дискретных схемах
		*  невозможно интегрировать по дельта-функциям и их требуется избегать.
		*  Для потенциала erf(ar)/r фурье образ достаточно сложная функция и нет
		*  необходимости вычислять его точно, поскольку нужет только предел 
		*  при cutoff->inf 
		* @param k - абсолютное значение вектора |k| > 0
		* @return значение фурье образа потенциала (без нормировочного множителя)
		*/
		_Real fourier_near3(_Real k) const
		{
			assert(_GT(k, (_Real)M_SAFETY)); // фурье образ определен для |k| > 0
			return M_4PI * (1 - exp(-sqr(0.5* k / alpha_))) / sqr(k);
		}
		
		_Real fourier_far3(_Real k) const
		{
			assert(_GT(k, (_Real)M_SAFETY)); // фурье образ определен для |k| > 0
			return M_4PI * exp(-sqr(0.5* k / alpha_)) / sqr(k);
		}
		
		_Real fourier3(_Real k) const
		{
			assert(_GT(k, (_Real)M_SAFETY)); // фурье образ определен для |k| > 0
			return M_4PI / sqr(k);
		}
		/**
		* Считает нормировочный множитель фурье образа для 3D-случая.
		* @return нормировочный множитель
		*/
		_Real fourier_norm3() const { return 1./ (M_2PI * sqrt(M_2PI)); }
	};
	
}
#endif
