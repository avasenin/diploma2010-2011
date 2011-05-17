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
#ifndef _RANDOM__F9ED1116_22BE_5068_B247_C1442B540601__H
#define _RANDOM__F9ED1116_22BE_5068_B247_C1442B540601__H

#include <boost/random/normal_distribution.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/variate_generator.hpp>

#include "prgkern/_prgconfig.h"

#ifdef USE_LINUX_SPECIFIC_CODE
	#include <sys/time.h>
	#include <sys/utsname.h>
#endif

namespace prgkern
{
	/**
	* @brief random number generator
	* @note This is an corrected implementation of the algorithm used in Numerical
	*   Recipe's ran0 generator. It is the same as MINSTD with an XOR mask
	*   of 123459876 on the seed. The period of this generator is 2^31.
	*/
	class randgen
	{
	#define RAN0_MASK 123459876
	#define RAN0_MAX 2147483646

	public:
		/// constructor
		randgen(long seed=1) : idum_(seed_ctrl_(seed)) {}

		/// sets seed of generator
		void seed(long s) { idum_ = seed_ctrl_(s); }

		/// sets random seed of generator
		void seed()
		{
		#ifdef USE_LINUX_SPECIFIC_CODE
			struct timeval tv;
			::gettimeofday(&tv, NULL);
			idum_ = seed_ctrl_(tv.tv_usec);
		#endif
		}

		template <typename T>
		T operator()(T max=1, T min=0) { return (T)(min + make_() * (max-min)); }
			// использовать обратный порядок [min, max] безопасно
			// как плюс, это дает возможность вызывать функцию с одним параметром

	protected:
		/**
		* @note if you choose a seed of 123459876 it would give a degenerate
		*   series 0,0,0,0, ...  I've made the code to avoid that error.
		*/
		long seed_ctrl_(long s)	{	return s == RAN0_MASK ? s : s ^ RAN0_MASK; }

		double make_()
		{
			long k = idum_ / 127773;
			idum_ = 16807 * (idum_ - k * 127773) - 2836 * k;
			if (idum_ < 0) idum_ += (RAN0_MAX + 1);
			return ((double)idum_) / (RAN0_MAX + 1);
		}

	private:
		long idum_;

	#undef RAN0_MAX
	#undef RAN0_MASK
	};

	template <typename _Real>
	INLINE _Real ran0(_Real max=1., _Real min=0.)
	{
		static randgen random;
		return random(max, min);
	}

	INLINE unsigned ran0(unsigned max) { return (unsigned)(max * ran0<double>()); }

	/**
	*  Функция нормального распределения.
	*/
	template <typename _Real>
	class normal_distribution_ : public boost::normal_distribution<_Real>
	{
		typedef boost::normal_distribution<_Real>  _Base;

	public:
		typedef typename boost::normal_distribution<_Real>::result_type result_type;

		/// конструктор
		normal_distribution_(_Real sigma, _Real mean=0.)
		: boost::normal_distribution<_Real>(mean, sigma) {}

		/**
		* Возвращает значение нормального распределения без нормировочного множителя
		* @note НЕЛЬЗЯ в качестве этой функции использовать operator(), поскольку
		*  это закрывает подобный оператор в базовом классе, используемый для
		*  генерации случайных чисел.
		*  Предположительно это дефект boost::math::statistical библиотеки.
		*/
		_Real value(_Real x)
		{
			x -= _Base::mean();
			x /= _Base::sigma();
			return exp(-0.5 * sqr(x));
		}

		/// коэффициент нормировки распределения
		_Real norm() { return 1. / (sqrt(M_PI + M_PI) * _Base::sigma()); }

	};
	typedef boost::rand48 rand48;

	/**
	*  Генератор случайных чисел, вероятность получения которых согласована
	*  с заданной функцией распределения.
	* @param G генератор равномерно распределенного случайного числа
	* @param D функцией распределения, согласно которому равномерное распределение
	*  корректируется
	*/
	template <class Distribution, class Generator=rand48>
	class variate_generator_
	: public boost::variate_generator<Generator, Distribution>
	{
	public:

		/**
		*  Конструктор генератора
		* @param g генератор равномерного распределения
		* @param d целевое распределение
		*/
		variate_generator_(Distribution d, Generator g=Generator())
		: boost::variate_generator<Generator, Distribution>(g, d) {}
	};

}
#endif
