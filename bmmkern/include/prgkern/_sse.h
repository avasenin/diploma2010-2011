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
#ifndef _SSE___0077A726_E6E3_58a1_C16D_CE436AB90502__H
#define _SSE___0077A726_E6E3_58a1_C16D_CE436AB90502__H

#ifdef USE_SSE
	#include <math.h>
	#include <xmmintrin.h>
	#include <emmintrin.h>
#endif
#include "prgkern/_debug.h"
#include "prgkern/_string.h"

/**
 * @note При работе с SSE замечено (для gcc), что преобразования (__m128i)m128 и (__m128)m128i не
 * меняют битовое представление числа, то есть не являются аналогами преобразования float <-> int.
 * Таким образом, незачем использовать преобразование *(__m128i*)&f, поскольку оно более долгое.
 *
 * @note Обертки SSE векторов построены таким образом, чтобы обеспечить сохранение эффективности.
 * По этой причине, обертки не содержат операций с ссылками. Например,
 *     v4_int operator+(v4_int a, v4_int b);  работает эффективно
 *     v4_int operator+(const v4_int &a, const v4_int &b);  работаеи НЕ эффективно, поскольку
 * ссылка реализуется компилятором через указатель и требует дополнительного разыменования
 *
 * @note Некоторые операции отсутствуют, так как их нет в sse intristic. Например, нет операции
 * умножения для целых чисел (включающей перемножение 4-х чисел).
 */

namespace prgkern
{
	// избегаем сокрытия стандартных sqrt(), abs() и т.д., из-за того, что SSE вводит собственные
	using std::sqrt;
	using std::abs;
	using std::floor;
	using std::ceil;

	/**
	 * @brief булевское значение для операций с плавающими числами
	 */
	template <unsigned N, typename T> class vecbool_;
	template <unsigned N, typename T> class vecint_;
	template <unsigned N, typename T> class vecreal_;

	template <typename T> INLINE T summarize(T v) { return v; }

#ifdef USE_SSE

	const unsigned non_empty_element_count[16] =
	{
		0, // 0x0000
		1, // 0x0001
		1, // 0x0010
		2, // 0x0011
		1, // 0x0100
		2, // 0x0101
		2, // 0x0110
		3, // 0x0111
		1, // 0x1000
		2, // 0x1001
		2, // 0x1010
		3, // 0x1011
		2, // 0x1100
		3, // 0x1101
		3, // 0x1110
		4  // 0x1111
	};

	template <> class vecbool_<4, int>
	{
		__m128 b_;

	public:
		enum { size = 4 };
		typedef int value_type;

		vecbool_() : b_(_mm_setzero_ps()) {}
		vecbool_(__m128 b) : b_(b) {}
		vecbool_(__m128i b) : b_((__m128)b) {}
		vecbool_(const vecbool_ &b) : b_(b.b_) {}

		vecbool_(bool b) : b_(_mm_cmpgt_ps(_mm_set1_ps(b), _mm_setzero_ps())) {}
			// данная реализация конструктора из bool корректна, несмотря на использование
			// сравнения ">" с нулем, поскольку bool хранит только два значения {0, 1}

		vecbool_(unsigned b) : b_(_mm_cmpgt_ps(_mm_set1_ps(b), _mm_setzero_ps())) {}
			// данная реализация конструктора из unsigned корректна, несмотря на использование
			// сравнения ">" с нулем, поскольку unsigned хранит положительные значения.

		vecbool_(int b) : b_(_mm_cmpgt_ps(_mm_set1_ps((unsigned)b), _mm_setzero_ps())) {}
			// Заметим, что операция конструирования из int будет работает только через
			// преобразование к unsigned, а без него работает некорректно.

		vecbool_ &operator=(vecbool_ b) { b_ = b.b_; return *this; }
		vecbool_ &operator=(bool b) { *this = vecbool_(b); return *this; }
		vecbool_ &operator=(unsigned b) { *this = vecbool_(b); return *this; }
		vecbool_ &operator=(int b) { *this = vecbool_(b); return *this; }

		int operator[](unsigned i) const { return *((int*)&b_ + i); }
		int &operator[](unsigned i) { return *((int*)&b_ + i); }

		operator bool() { return _mm_movemask_ps(b_) != 0; }
		operator __m128() { return b_; }
		operator __m128i() { return (__m128i)b_; }

		/// число ненулевых компонент
		unsigned count() const { return non_empty_element_count[_mm_movemask_ps(b_)]; }
	};

	INLINE vecbool_<4, int> operator!(vecbool_<4, int> a) { return vecbool_<4, int>(_mm_cmpeq_ps(a, _mm_setzero_ps())); }
	INLINE vecbool_<4, int> operator&(vecbool_<4, int> a, vecbool_<4, int> b) { return vecbool_<4, int>(_mm_and_ps(a, b)); }
	INLINE vecbool_<4, int> operator|(vecbool_<4, int> a, vecbool_<4, int> b) { return vecbool_<4, int>(_mm_or_ps (a, b)); }
	INLINE vecbool_<4, int> operator^(vecbool_<4, int> a, vecbool_<4, int> b) { return vecbool_<4, int>(_mm_xor_ps(a, b)); }
	INLINE bool operator==(vecbool_<4, int> a, bool b) { return (bool)a == b; }
	INLINE bool operator==(bool a, vecbool_<4, int> b) { return a == (bool)b; }

	INLINE vecbool_<4, int> multiple_not(vecbool_<4, int> a, vecbool_<4, int> b) { return vecbool_<4, int>(_mm_andnot_ps(a, b)); }
	INLINE vecbool_<4, int> multiple(vecbool_<4, int> a, vecbool_<4, int> b) { return vecbool_<4, int>(_mm_and_ps(a, b)); }

	/**
	* @brief операция c = !a && b
	* @note важно не забывать, что отрицание относится к первому аргументу
	*/
	INLINE vecbool_<4, int> andnot(vecbool_<4, int> a, vecbool_<4, int> b) { return vecbool_<4, int>(_mm_andnot_ps(a, b)); }
		// операция эмулирует соответствующую SSE операцию и это быстрее чем делать через явные ~a & b

	/**
	 * @brief класс-обертка для операций с целыми числами
	 */
	template <> class vecint_<4, int>
	{
		__m128i i_;

	public:
		enum { size = 4 };
		typedef int value_type;

		vecint_() : i_(_mm_setzero_si128()) {}
		vecint_(__m128i i) : i_(i) {}
		vecint_(vecbool_<4, int> b) : i_((__m128i)b) {}
		vecint_(const vecint_ &i) : i_(i.i_) {}
		vecint_(int i) : i_(_mm_set1_epi32(i)) {}
		vecint_(int i0, int i1, int i2, int i3) { i_ = _mm_set_epi32(i3, i2, i1, i0); }

		vecint_ &operator=(vecint_ i) { i_ = i.i_; return *this; }
		vecint_ &operator=(int i) { i_ = _mm_set1_epi32(i); return *this; }

		int operator[](unsigned i) const { return *((int*)&i_ + i); }
		int &operator[](unsigned i) { return *((int*)&i_ + i); }

		vecint_ &operator++() { i_ = _mm_add_epi32(i_, _mm_set1_epi32(1)); return *this; }
		vecint_ &operator--() { i_ = _mm_sub_epi32(i_, _mm_set1_epi32(1)); return *this; }
		vecint_ &operator+=(vecint_ i) { i_ = _mm_add_epi32(i_, i.i_); return *this; }
		vecint_ &operator-=(vecint_ i) { i_ = _mm_sub_epi32(i_, i.i_); return *this; }

		/**
		* @brief эмулятор булевой операции a = a && b
		* @note Данная операция должна иметь название &&=, поскольку она логическая, а не битовая,
		* но такого оператора в C++ нет. Операция испоьзуется для очистки тех компонент вектора,
		* для которых в соответствующих компонентах булевого аргумента записаны нули.
		*/
		vecint_ &operator&=(vecbool_<4, int> b) { i_ = _mm_and_si128(i_, (__m128i)b); return *this;  }

		/**
		* @brief эмулятор битовой операции a = a & b
		* @note Операция испоьзуется для очистки тех бит вектора, для которых в соответствующих
		* битах аргумента записаны нули.
		*/
		vecint_ &operator&=(vecint_ i) { i_ = _mm_and_si128(i_, i.i_); return *this;  }

		operator __m128i() { return i_; }

		void load(const int *i) { i_ = _mm_loadu_si128((__m128i*)i); }
		void store(int *p) const { _mm_storeu_si128((__m128i*)p, i_); }
	};

	INLINE vecint_<4, int> operator+(vecint_<4, int> a, vecint_<4, int> b) { return _mm_add_epi32(a, b); }
	INLINE vecint_<4, int> operator-(vecint_<4, int> a, vecint_<4, int> b) { return _mm_sub_epi32(a, b); }

	INLINE vecint_<4, int> operator+(vecint_<4, int> a, int b) { return _mm_add_epi32(a, _mm_set1_epi32(b)); }
	INLINE vecint_<4, int> operator-(vecint_<4, int> a, int b) { return _mm_sub_epi32(a, _mm_set1_epi32(b)); }
	INLINE vecint_<4, int> operator+(int a, vecint_<4, int> b) { return _mm_add_epi32(_mm_set1_epi32(a), b); }
	INLINE vecint_<4, int> operator-(int a, vecint_<4, int> b) { return _mm_sub_epi32(_mm_set1_epi32(a), b); }

	INLINE vecbool_<4, int> operator< (vecint_<4, int> a, vecint_<4, int> b) { return (vecbool_<4, int>)_mm_cmplt_epi32(a, b); }
	INLINE vecbool_<4, int> operator> (vecint_<4, int> a, vecint_<4, int> b) { return (vecbool_<4, int>)_mm_cmpgt_epi32(a, b); }
	INLINE vecbool_<4, int> operator==(vecint_<4, int> a, vecint_<4, int> b) { return (vecbool_<4, int>)_mm_cmpeq_epi32(a, b); }

	INLINE vecbool_<4, int> operator< (vecint_<4, int> a, int b) { return (vecbool_<4, int>)_mm_cmplt_epi32(a, _mm_set1_epi32(b)); }
	INLINE vecbool_<4, int> operator> (vecint_<4, int> a, int b) { return (vecbool_<4, int>)_mm_cmpgt_epi32(a, _mm_set1_epi32(b)); }
	INLINE vecbool_<4, int> operator==(vecint_<4, int> a, int b) { return (vecbool_<4, int>)_mm_cmpeq_epi32(a, _mm_set1_epi32(b)); }

	INLINE vecbool_<4, int> operator< (int a, vecint_<4, int> b) { return (vecbool_<4, int>)_mm_cmplt_epi32(_mm_set1_epi32(a), b); }
	INLINE vecbool_<4, int> operator> (int a, vecint_<4, int> b) { return (vecbool_<4, int>)_mm_cmpgt_epi32(_mm_set1_epi32(a), b); }
	INLINE vecbool_<4, int> operator==(int a, vecint_<4, int> b) { return (vecbool_<4, int>)_mm_cmpeq_epi32(_mm_set1_epi32(a), b); }

	INLINE vecint_<4, int> operator&(vecint_<4, int> a, vecint_<4, int> b) { return vecint_<4, int> (_mm_and_si128(a, b)); }
	INLINE vecint_<4, int> operator|(vecint_<4, int> a, vecint_<4, int> b) { return vecint_<4, int> (_mm_or_si128 (a, b)); }
	INLINE vecint_<4, int> operator^(vecint_<4, int> a, vecint_<4, int> b) { return vecint_<4, int> (_mm_xor_si128(a, b)); }

	INLINE vecint_<4, int> operator&(vecint_<4, int> a, vecbool_<4, int> b) { return vecint_<4, int>(_mm_and_si128(a, b)); }
	INLINE vecint_<4, int> operator&(vecbool_<4, int> a, vecint_<4, int> b) { return vecint_<4, int>(_mm_and_si128(a, b)); }

	INLINE vecint_<4, int> multiple_not(vecbool_<4, int> a, vecint_<4, int> b) { return vecint_<4, int>(_mm_andnot_si128(a, b)); }
	INLINE vecint_<4, int> multiple    (vecbool_<4, int> a, vecint_<4, int> b) { return vecint_<4, int>(   _mm_and_si128(a, b)); }

	/**
	 * Сдвиг влево каждого элемента вектора на заданное в параметре число бит.
	 * @note В отличие от SSE операций сдвиг для каждого элемента свой, а не единый для всех.
	 * @param a исходный вектор
	 * @param n вектор смещений
	 * @return вектор со смещенными компонентами
	 */
	INLINE vecint_<4, int> operator<<(vecint_<4, int> a, vecint_<4, int> n)
	{
		return vecint_<4, int>(a[0] << n[0], a[1] << n[1], a[2] << n[2], a[3] << n[3]);
	}

	/**
	* @brief Получение абсолютного значения смещения.
	* @note см. Уоррен. Алгоритмические трюки для программистов
	*/
	INLINE vecint_<4, int> abs(vecint_<4, int> a) { vecbool_<4, int> f = (a < 0); return _mm_sub_epi32(_mm_xor_si128(f, a), f); }

	/**
	* @brief Получение абсолютного значения смещения.
	* @note Данная операция чуть эффективней, чем операция с одним аргументом.
	* Второй аргумент используется для указания на то, какие компоненты вектора меньше нуля.
	*/
	INLINE vecint_<4, int> abs(vecint_<4, int> a, vecbool_<4, int> less_zero)
	{ return _mm_sub_epi32(_mm_xor_si128(less_zero, a), less_zero); }

	/**
	* @brief Циклическое вращение компонент вектора на соответствующее число позиций.
	*/
	INLINE vecint_<4, int> ror0(vecint_<4, int> i) { return i; }
	INLINE vecint_<4, int> ror1(vecint_<4, int> i) { return _mm_shuffle_epi32(i, _MM_SHUFFLE(0, 3, 2, 1)); }
	INLINE vecint_<4, int> ror2(vecint_<4, int> i) { return _mm_shuffle_epi32(i, _MM_SHUFFLE(1, 0, 3, 2)); }
	INLINE vecint_<4, int> ror3(vecint_<4, int> i) { return _mm_shuffle_epi32(i, _MM_SHUFFLE(2, 1, 0, 3)); }

	INLINE vecint_<4, int> rol0(vecint_<4, int> i) { return i; }
	INLINE vecint_<4, int> rol1(vecint_<4, int> i) { return _mm_shuffle_epi32(i, _MM_SHUFFLE(2, 1, 0, 3)); }
	INLINE vecint_<4, int> rol2(vecint_<4, int> i) { return _mm_shuffle_epi32(i, _MM_SHUFFLE(1, 0, 3, 2)); }
	INLINE vecint_<4, int> rol3(vecint_<4, int> i) { return _mm_shuffle_epi32(i, _MM_SHUFFLE(0, 3, 2, 1)); }

	/**
	 * @brief класс-обертка для операций с float числами
	 */
	template<> class vecreal_<4, float>
	{
		__m128 f_;

	public:
		enum { size = 4 };
		typedef float value_type;

		vecreal_() : f_(_mm_setzero_ps()) {}
		vecreal_(__m128 f) : f_(f) {}
		vecreal_(const vecreal_ &f) : f_(f.f_) {}
		vecreal_(float f) : f_(_mm_set1_ps(f)) {}
		vecreal_(float f0, float f1, float f2, float f3) : f_(_mm_set_ps(f3, f2, f1, f0)) {}
		vecreal_(const float *f) : f_(_mm_loadu_ps(f)) {}

		vecreal_ &operator=(const vecreal_ &f) { f_ = f.f_; return *this; }
		vecreal_ &operator=(float f) { f_ = _mm_set1_ps(f); return *this; }

		float operator[](unsigned i) const { return *((float*)&f_ + i); }
		float &operator[](unsigned i) { return *((float*)&f_ + i); }

		vecreal_ &operator+=(vecreal_ f) { f_ = _mm_add_ps(f_, f.f_); return *this; }
		vecreal_ &operator-=(vecreal_ f) { f_ = _mm_sub_ps(f_, f.f_); return *this; }
		vecreal_ &operator*=(vecreal_ f) { f_ = _mm_mul_ps(f_, f.f_); return *this; }
		vecreal_ &operator/=(vecreal_ f) { f_ = _mm_div_ps(f_, f.f_); return *this; }

		vecreal_ &operator+=(float f) { f_ = _mm_add_ps(f_, _mm_set1_ps(f)); return *this; }
		vecreal_ &operator-=(float f) { f_ = _mm_sub_ps(f_, _mm_set1_ps(f)); return *this; }
		vecreal_ &operator*=(float f) { f_ = _mm_mul_ps(f_, _mm_set1_ps(f)); return *this; }
		vecreal_ &operator/=(float f) { f_ = _mm_div_ps(f_, _mm_set1_ps(f)); return *this; }

		/**
		* @brief эмулятор булевой операции a = a && b
		* @note Данная операция должна иметь название &&=, поскольку она логическая, а не битовая,
		* но такого оператора в C++ нет. Операция испоьзуется для очистки тех компонент вектора,
		* для которых в соответствующих компонентах булевого аргумента записаны нули.
		*/
		vecreal_ &operator&=(vecbool_<4, int> b) { f_ = _mm_and_ps(f_, (__m128)b); return *this; }

		operator __m128() { return f_; }

		void load(const float *f) { f_ = _mm_loadu_ps(f); }
		void store(float *p) const { _mm_storeu_ps(p, f_); }
	};

	INLINE vecreal_<4, float> operator+(vecreal_<4, float> a, vecreal_<4, float> b) { return _mm_add_ps(a, b); }
	INLINE vecreal_<4, float> operator-(vecreal_<4, float> a, vecreal_<4, float> b) { return _mm_sub_ps(a, b); }
	INLINE vecreal_<4, float> operator*(vecreal_<4, float> a, vecreal_<4, float> b) { return _mm_mul_ps(a, b); }
	INLINE vecreal_<4, float> operator/(vecreal_<4, float> a, vecreal_<4, float> b) { return _mm_div_ps(a, b); }
	INLINE vecreal_<4, float> operator+(vecreal_<4, float> a, float b) { return _mm_add_ps(a, _mm_set1_ps(b)); }
	INLINE vecreal_<4, float> operator-(vecreal_<4, float> a, float b) { return _mm_sub_ps(a, _mm_set1_ps(b)); }
	INLINE vecreal_<4, float> operator*(vecreal_<4, float> a, float b) { return _mm_mul_ps(a, _mm_set1_ps(b)); }
	INLINE vecreal_<4, float> operator/(vecreal_<4, float> a, float b) { return _mm_div_ps(a, _mm_set1_ps(b)); }
	INLINE vecreal_<4, float> operator+(float a, vecreal_<4, float> b) { return _mm_add_ps(_mm_set1_ps(a), b); }
	INLINE vecreal_<4, float> operator-(float a, vecreal_<4, float> b) { return _mm_sub_ps(_mm_set1_ps(a), b); }
	INLINE vecreal_<4, float> operator*(float a, vecreal_<4, float> b) { return _mm_mul_ps(_mm_set1_ps(a), b); }
	INLINE vecreal_<4, float> operator/(float a, vecreal_<4, float> b) { return _mm_div_ps(_mm_set1_ps(a), b); }

	INLINE vecreal_<4, float> operator+(vecreal_<4, float> a, vecint_<4, int> b)
	{ return _mm_add_ps(a, _mm_cvtepi32_ps(b)); }

	INLINE vecreal_<4, float> operator-(vecreal_<4, float> a, vecint_<4, int> b)
	{ return _mm_sub_ps(a, _mm_cvtepi32_ps(b)); }

	INLINE vecreal_<4, float> operator*(vecreal_<4, float> a, vecint_<4, int> b)
	{ return _mm_mul_ps(a, _mm_cvtepi32_ps(b)); }

	INLINE vecreal_<4, float> operator/(vecreal_<4, float> a, vecint_<4, int> b)
	{ return _mm_div_ps(a, _mm_cvtepi32_ps(b)); }

	INLINE vecreal_<4, float> operator+(vecint_<4, int> a, vecreal_<4, float> b)
	{ return _mm_add_ps(_mm_cvtepi32_ps(a), b); }

	INLINE vecreal_<4, float> operator-(vecint_<4, int> a, vecreal_<4, float> b)
	{ return _mm_sub_ps(_mm_cvtepi32_ps(a), b); }

	INLINE vecreal_<4, float> operator*(vecint_<4, int> a, vecreal_<4, float> b)
	{ return _mm_mul_ps(_mm_cvtepi32_ps(a), b); }

	INLINE vecreal_<4, float> operator/(vecint_<4, int> a, vecreal_<4, float> b)
	{ return _mm_div_ps(_mm_cvtepi32_ps(a), b); }

	INLINE vecreal_<4, float> sqr(vecreal_<4, float> a) { return _mm_mul_ps(a, a); }
	INLINE vecreal_<4, float> cube(vecreal_<4, float> a) { return _mm_mul_ps(_mm_mul_ps(a, a), a); }
	INLINE vecreal_<4, float> sqrt(vecreal_<4, float> a) { return _mm_sqrt_ps(a); }

	INLINE vecreal_<4, float> max(vecreal_<4, float> a, vecreal_<4, float> b) { return _mm_max_ps(a, b); }
	INLINE vecreal_<4, float> min(vecreal_<4, float> a, vecreal_<4, float> b) { return _mm_min_ps(a, b); }

	/// объект, востанавливающих моду округления после использования функций округления
	/// При его отсутствии нарушается работа даже sin(), cos() и т.д.
	struct round_mode_saver_
	{
		round_mode_saver_() : mode_(_MM_GET_ROUNDING_MODE()) {}
		~round_mode_saver_() { _MM_SET_ROUNDING_MODE(mode_); }
		unsigned mode_;
	};

	#define DECLARE_AS_RESTORING_ROUND_MODE \
		round_mode_saver_ mode##0077A726_E6E3_58a1_C16D_CE436AB90502;

	INLINE vecreal_<4, float> floor(vecreal_<4, float> a)
	{
		DECLARE_AS_RESTORING_ROUND_MODE
		_MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
		return vecreal_<4, float>(_mm_cvtepi32_ps(_mm_cvtps_epi32(a)));
	}

	INLINE vecreal_<4, float> ceil(vecreal_<4, float> a)
	{
		DECLARE_AS_RESTORING_ROUND_MODE
		_MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
		return vecreal_<4, float>(_mm_cvtepi32_ps(_mm_cvtps_epi32(a)));
	}

	INLINE vecint_<4, int> truncate(vecreal_<4, float> a)
	{
		DECLARE_AS_RESTORING_ROUND_MODE
		_MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);
		return vecint_<4, int>(_mm_cvtps_epi32(a));
	}

	INLINE vecreal_<4, float> round(vecreal_<4, float> a)
	{
		//DECLARE_AS_RESTORING_ROUND_MODE
		//_MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
		return vecreal_<4, float>(_mm_cvtepi32_ps(_mm_cvtps_epi32(a)));
	}

	INLINE vecint_<4, int> ifloor(vecreal_<4, float> a)
	{
		DECLARE_AS_RESTORING_ROUND_MODE
		_MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);
		return vecint_<4, int>(_mm_cvtps_epi32(a));
	}

	INLINE vecint_<4, int> iceil(vecreal_<4, float> a)
	{
		DECLARE_AS_RESTORING_ROUND_MODE
		_MM_SET_ROUNDING_MODE(_MM_ROUND_UP);
		return vecint_<4, int>(_mm_cvtps_epi32(a));
	}

	INLINE vecreal_<4, float> exp(vecreal_<4, float> a)
	{
		return vecreal_<4, float>(::exp(a[0]), ::exp(a[1]), ::exp(a[2]), ::exp(a[3]));
	}

	INLINE vecreal_<4, float> erf(vecreal_<4, float> a)
	{
		return vecreal_<4, float>(::erf(a[0]), ::erf(a[1]), ::erf(a[2]), ::erf(a[3]));
	}

	// (!) основной источник неточности расчета в SSE это эти две функции
	INLINE vecreal_<4, float> rcp(vecreal_<4, float> a)
	{
		_MSG("Do you really want use this fast but inaccurate function?");
		return _mm_rcp_ps(a);
	}
	INLINE vecreal_<4, float> rsqrt(vecreal_<4, float> a)
	{
		_MSG("Do you really want use this fast but inaccurate function?");
		return _mm_rsqrt_ps(a);
	}

	INLINE vecbool_<4, int> operator< (vecreal_<4, float> a, vecreal_<4, float> b)
	{ return (vecbool_<4, int>)_mm_cmplt_ps(a, b); }

	INLINE vecbool_<4, int> operator> (vecreal_<4, float> a, vecreal_<4, float> b)
	{ return (vecbool_<4, int>)_mm_cmpgt_ps(a, b); }

	INLINE vecbool_<4, int> operator< (vecreal_<4, float> a, float b)
	{ return (vecbool_<4, int>)_mm_cmplt_ps(a, _mm_set1_ps(b)); }

	INLINE vecbool_<4, int> operator> (vecreal_<4, float> a, float b)
	{ return (vecbool_<4, int>)_mm_cmpgt_ps(a, _mm_set1_ps(b)); }

	INLINE vecbool_<4, int> operator< (float a, vecreal_<4, float> b)
	{ return (vecbool_<4, int>)_mm_cmplt_ps(_mm_set1_ps(a), b); }

	INLINE vecbool_<4, int> operator> (float a, vecreal_<4, float> b)
	{ return (vecbool_<4, int>)_mm_cmpgt_ps(_mm_set1_ps(a), b); }

	INLINE vecreal_<4, float> multiple_not(vecbool_<4, int> a, vecreal_<4, float> b)
	{ return vecreal_<4, float>(_mm_andnot_ps(a, b)); }

	INLINE vecreal_<4, float> multiple    (vecbool_<4, int> a, vecreal_<4, float> b)
	{ return vecreal_<4, float>(   _mm_and_ps(a, b)); }

	/**
	* @brief Циклическое вращение компонент вектора на соответствующее число позиций.
	*/
  INLINE vecreal_<4, float> ror0(vecreal_<4, float> f) { return f; }
  INLINE vecreal_<4, float> ror1(vecreal_<4, float> f) { return _mm_shuffle_ps(f, f, _MM_SHUFFLE(0, 3, 2, 1)); }
  INLINE vecreal_<4, float> ror2(vecreal_<4, float> f) { return _mm_shuffle_ps(f, f, _MM_SHUFFLE(1, 0, 3, 2)); }
  INLINE vecreal_<4, float> ror3(vecreal_<4, float> f) { return _mm_shuffle_ps(f, f, _MM_SHUFFLE(2, 1, 0, 3)); }

  INLINE vecreal_<4, float> rol0(vecreal_<4, float> f) { return f; }
  INLINE vecreal_<4, float> rol1(vecreal_<4, float> f) { return _mm_shuffle_ps(f, f, _MM_SHUFFLE(2, 1, 0, 3)); }
  INLINE vecreal_<4, float> rol2(vecreal_<4, float> f) { return _mm_shuffle_ps(f, f, _MM_SHUFFLE(1, 0, 3, 2)); }
  INLINE vecreal_<4, float> rol3(vecreal_<4, float> f) { return _mm_shuffle_ps(f, f, _MM_SHUFFLE(0, 3, 2, 1)); }

	INLINE int summarize(vecint_<4, int> v) { return v[0] + v[1] + v[2] + v[3]; }
	INLINE float summarize(vecreal_<4, float> v) { return v[0] + v[1] + v[2] + v[3]; }

#else // USE_SSE

	#define DECLARE_AS_RESTORING_ROUND_MODE
	#define _MM_SET_ROUNDING_MODE(_)

	#define VECTOR_OPERATOR(v_, op, v) { v_[0] op v[0]; v_[1] op v[1]; v_[2] op v[2]; v_[3] op v[3]; }
	#define SCALAR_OPERATOR(v_, op, v) { v_[0] op v; v_[1] op v; v_[2] op v; v_[3] op v; }
	#define MULTISCALAR_OPERATOR(v_, op, v0, v1, v2, v3) { v_[0] op v0; v_[1] op v1; v_[2] op v2; v_[3] op v3; }

	#define VECTOR_ARG(a, op, b) a[0] op b[0], a[1] op b[1], a[2] op b[2], a[3] op b[3]
	#define SCALAR_ARG1(a, op, b) a op b[0], a op b[1], a op b[2], a op b[3]
	#define SCALAR_ARG2(a, op, b) a[0] op b, a[1] op b, a[2] op b, a[3] op b

	/**
	 * @brief булевское значение для операций с плавающими числами
	 * Хранит все единичные биты (true) и все нулевые биты (false) для совместимости с SSE.
	 */
	template <>	class vecbool_<4, int>
	{
		int v_[4];

	public:
		enum { size = 4 };
		typedef int value_type;

		vecbool_() { SCALAR_OPERATOR(v_, =, 0); }

		vecbool_(bool v0, bool v1, bool v2, bool v3)
		{
			v_[0] = v0 ? -1 : 0;
			v_[1] = v1 ? -1 : 0;
			v_[2] = v2 ? -1 : 0;
			v_[3] = v3 ? -1 : 0;
		}

		vecbool_(const vecbool_ &b) { VECTOR_OPERATOR(v_, =, b); }

		vecbool_ &operator=(vecbool_ b) { VECTOR_OPERATOR(v_, =, b); return *this; }

		vecbool_(bool b)
		{
			v_[0] = b ? -1 : 0;
			v_[1] = v_[0];
			v_[2] = v_[0];
			v_[3] = v_[0];
		}

		value_type operator[](unsigned i) const { return v_[i]; }
		value_type &operator[](unsigned i) { return v_[i]; }

		operator bool() { return (bool)(v_[0] || v_[1] || v_[2] || v_[3]); }

		/// число ненулевых компонент
		unsigned count() const { return std::abs(v_[0] + v_[1] + v_[2] + v_[3]); }
	};

	INLINE vecbool_<4, int> operator!(vecbool_<4, int> a) { return vecbool_<4, int>(SCALAR_ARG2(a, ==, 0)); }
	INLINE vecbool_<4, int> operator&&(vecbool_<4, int> a, vecbool_<4, int> b) { return vecbool_<4, int>(VECTOR_ARG(a, &, b)); }
	INLINE bool operator&&(vecbool_<4, int> a, bool b) { return (bool)a && b; }
	INLINE bool operator&&(bool a, vecbool_<4, int> b) { return a && (bool)b; }

	INLINE vecbool_<4, int> operator||(vecbool_<4, int> a, vecbool_<4, int> b) { return vecbool_<4, int>(VECTOR_ARG(a, |, b)); }
	INLINE bool operator==(vecbool_<4, int> a, bool b) { return (bool)a == b; }
	INLINE bool operator==(bool a, vecbool_<4, int> b) { return a == (bool)b; }

	INLINE vecbool_<4, int> multiple_not(vecbool_<4, int> a, vecbool_<4, int> b) { return (!a) && b; }
	INLINE vecbool_<4, int> multiple(vecbool_<4, int> a, vecbool_<4, int> b) { return a && b; }
	INLINE vecbool_<4, int> andnot(vecbool_<4, int> a, vecbool_<4, int> b) { return (!a) && b; }

	/**
	 * @brief класс-обертка для операций с целыми числами
	 */
	template <> class vecint_<4, int>
	{
		int v_[4];

	public:
		enum { size = 4 };
		typedef int value_type;

		vecint_() { SCALAR_OPERATOR(v_, =, 0); }
		vecint_(const vecint_ &i) { VECTOR_OPERATOR(v_, =, i); }
		vecint_(int i) { SCALAR_OPERATOR(v_, =, i); }
		vecint_(int i0, int i1, int i2, int i3) { MULTISCALAR_OPERATOR(v_, =, i0, i1, i2, i3); }
		vecint_(const vecbool_<4, int> &b) { VECTOR_OPERATOR(v_, = , b); }

		vecint_ &operator=(const vecint_ &i) { VECTOR_OPERATOR(v_, =, i); return *this; }
		vecint_ &operator=(int i) { SCALAR_OPERATOR(v_, =, i); return *this; }

		value_type operator[](unsigned i) const { return v_[i]; }
		value_type &operator[](unsigned i) { return v_[i]; }

		vecint_ &operator+=(vecint_ i) { VECTOR_OPERATOR(v_, +=, i); return *this; }
		vecint_ &operator-=(vecint_ i) { VECTOR_OPERATOR(v_, -=, i); return *this; }
	};

	INLINE vecint_<4, int> operator+(vecint_<4, int> a, vecint_<4, int> b) { return a += b; }
	INLINE vecint_<4, int> operator-(vecint_<4, int> a, vecint_<4, int> b) { return a -= b; }
	INLINE vecint_<4, int> operator+(vecint_<4, int> a, int b) { return a += vecint_<4, int>(b); }
	INLINE vecint_<4, int> operator-(vecint_<4, int> a, int b) { return a -= vecint_<4, int>(b); }
	INLINE vecint_<4, int> operator+(int a, vecint_<4, int> b) { vecint_<4, int> a_(a); return a_ += b; }
	INLINE vecint_<4, int> operator-(int a, vecint_<4, int> b) { vecint_<4, int> a_(a); return a_ -= b; }

	INLINE vecbool_<4, int> operator< (vecint_<4, int> a, vecint_<4, int> b) { return vecbool_<4, int>(VECTOR_ARG(a, <, b)); }
	INLINE vecbool_<4, int> operator> (vecint_<4, int> a, vecint_<4, int> b) { return vecbool_<4, int>(VECTOR_ARG(a, >, b)); }
	INLINE vecbool_<4, int> operator==(vecint_<4, int> a, vecint_<4, int> b) { return vecbool_<4, int>(VECTOR_ARG(a,==, b)); }

	INLINE vecbool_<4, int> operator< (vecint_<4, int> a, int b) { return vecbool_<4, int>(SCALAR_ARG2(a, <, b)); }
	INLINE vecbool_<4, int> operator> (vecint_<4, int> a, int b) { return vecbool_<4, int>(SCALAR_ARG2(a, >, b)); }
	INLINE vecbool_<4, int> operator==(vecint_<4, int> a, int b) { return vecbool_<4, int>(SCALAR_ARG2(a,==, b)); }

	INLINE vecbool_<4, int> operator< (int a, vecint_<4, int> b) { return vecbool_<4, int>(SCALAR_ARG1(a, <, b)); }
	INLINE vecbool_<4, int> operator> (int a, vecint_<4, int> b) { return vecbool_<4, int>(SCALAR_ARG1(a, >, b)); }
	INLINE vecbool_<4, int> operator==(int a, vecint_<4, int> b) { return vecbool_<4, int>(SCALAR_ARG1(a,==, b)); }

	INLINE vecint_<4, int> operator&(vecint_<4, int> a, vecint_<4, int> b) { return vecint_<4, int>(VECTOR_ARG(a, &, b)); }
	INLINE vecint_<4, int> operator|(vecint_<4, int> a, vecint_<4, int> b) { return vecint_<4, int>(VECTOR_ARG(a, |, b)); }
	INLINE vecint_<4, int> operator^(vecint_<4, int> a, vecint_<4, int> b) { return vecint_<4, int>(VECTOR_ARG(a, ^, b)); }

	INLINE vecint_<4, int> multiple_not(vecbool_<4, int> a, vecint_<4, int> b) { return vecint_<4, int>(!a) & b; }
	INLINE vecint_<4, int> multiple    (vecbool_<4, int> a, vecint_<4, int> b) { return vecint_<4, int>( a) & b; }

	/**
	 * Сдвиг влево каждого элемента вектора на заданное в параметре число бит.
	 * @note В отличие от SSE операций сдвиг для каждого элемента свой, а не единый для всех.
	 * @param a исходный вектор
	 * @param n вектор смещений
	 * @return вектор со смещенными компонентами
	 */
	INLINE vecint_<4, int> operator<<(vecint_<4, int> a, vecint_<4, int> n)
	{ return vecint_<4, int>(VECTOR_ARG(a, <<, n)); }

	/**
	* @brief Получение абсолютного значения смещения.
	* @note см. Уоррен. Алгоритмические трюки для программистов
	*/
	INLINE vecint_<4, int> abs(vecint_<4, int> a)
	{
		return vecint_<4, int>
		(
			std::abs(a[0]),
			std::abs(a[1]),
			std::abs(a[2]),
			std::abs(a[3])
		);
	}

	/**
	 * @brief класс-обертка для операций с float числами
	 */
	template <typename T> class vecreal_<4, T>
	{
		T v_[4];

	public:

		enum { size = 4 };
		typedef T value_type;

		vecreal_() { SCALAR_OPERATOR(v_, =, 0.f); }
		template <typename S> vecreal_(S f) { SCALAR_OPERATOR(v_, =, f); }
		template <typename S> vecreal_(const vecreal_<4, S> &f) { VECTOR_OPERATOR(v_, =, f); }
		template <typename S> vecreal_(const vecint_<4, S> &f) { VECTOR_OPERATOR(v_, =, f); }
		template <typename S> vecreal_(S f0, S f1, S f2, S f3) { MULTISCALAR_OPERATOR(v_, =, f0, f1, f2, f3); }
		template <typename S> vecreal_(const S *f) { VECTOR_OPERATOR(v_, =, f); }

		template <typename S> vecreal_ &operator=(S v) { SCALAR_OPERATOR(v_, =, v); return *this; }
		template <typename S> vecreal_ &operator=(vecreal_<4, S> v) { VECTOR_OPERATOR(v_, =, v); return *this; }

		value_type operator[](unsigned i) const { return v_[i]; }
		value_type &operator[](unsigned i) { return v_[i]; }

		template <typename S> vecreal_ &operator+=(vecreal_<4, S> f) { VECTOR_OPERATOR(v_, +=, f); return *this; }
		template <typename S> vecreal_ &operator-=(vecreal_<4, S> f) { VECTOR_OPERATOR(v_, -=, f); return *this; }
		template <typename S> vecreal_ &operator*=(vecreal_<4, S> f) { VECTOR_OPERATOR(v_, *=, f); return *this; }
		template <typename S> vecreal_ &operator/=(vecreal_<4, S> f) { VECTOR_OPERATOR(v_, /=, f); return *this; }

		template <typename S> vecreal_ &operator+=(vecint_<4, S> f) { VECTOR_OPERATOR(v_, +=, f); return *this; }
		template <typename S> vecreal_ &operator-=(vecint_<4, S> f) { VECTOR_OPERATOR(v_, -=, f); return *this; }
		template <typename S> vecreal_ &operator*=(vecint_<4, S> f) { VECTOR_OPERATOR(v_, *=, f); return *this; }
		template <typename S> vecreal_ &operator/=(vecint_<4, S> f) { VECTOR_OPERATOR(v_, /=, f); return *this; }

		template <typename S> vecreal_ &operator+=(S f) { SCALAR_OPERATOR(v_, +=, f); return *this; }
		template <typename S> vecreal_ &operator-=(S f) { SCALAR_OPERATOR(v_, -=, f); return *this; }
		template <typename S> vecreal_ &operator*=(S f) { SCALAR_OPERATOR(v_, *=, f); return *this; }
		template <typename S> vecreal_ &operator/=(S f) { SCALAR_OPERATOR(v_, /=, f); return *this; }

		template <typename S> void store(S *p) const { *p = v_[0]; *(p + 1) = v_[1]; *(p + 2) = v_[2]; *(p + 3) = v_[3]; }
	};

	#define IMPLEMENT_4FUNCTION_2(name, type) \
	INLINE vecreal_<4, type> name(vecreal_<4, type> a, vecreal_<4, type> b) \
	{ \
		return vecreal_<4, type>(std::name(a[0], b[0]), std::name(a[1], b[1]), \
			std::name(a[2], b[2]), std::name(a[3], b[3])); \
	}

	IMPLEMENT_4FUNCTION_2(max, float);
	IMPLEMENT_4FUNCTION_2(max, double);
	IMPLEMENT_4FUNCTION_2(min, float);
	IMPLEMENT_4FUNCTION_2(min, double);

	#undef IMPLEMENT_4FUNCTION_2

	#define IMPLEMENT_4FUNCTION_1(name, type) \
	INLINE vecreal_<4, type> name(vecreal_<4, type> a) \
	{ \
		return vecreal_<4, type> (std::name(a[0]), std::name(a[1]), \
			std::name(a[2]), std::name(a[3])); \
	}

	IMPLEMENT_4FUNCTION_1(sqrt,  float);
	IMPLEMENT_4FUNCTION_1(sqrt,  double);
	IMPLEMENT_4FUNCTION_1(floor, float);
	IMPLEMENT_4FUNCTION_1(floor, double);
	IMPLEMENT_4FUNCTION_1(ceil,  float);
	IMPLEMENT_4FUNCTION_1(ceil,  double);

	#undef IMPLEMENT_4FUNCTION_1

	INLINE vecreal_<4, float> exp(vecreal_<4, float> a)
	{
		return vecreal_<4, float>(::exp(a[0]), ::exp(a[1]), ::exp(a[2]), ::exp(a[3]));
	}

	INLINE vecreal_<4, float> erf(vecreal_<4, float> a)
	{
		return vecreal_<4, float>(::erf(a[0]), ::erf(a[1]), ::erf(a[2]), ::erf(a[3]));
	}


	INLINE vecreal_<4, float> round(vecreal_<4, float> a)
	{
		return vecreal_<4, float> (roundf(a[0]), roundf(a[1]), roundf(a[2]), roundf(a[3]));
	}

	INLINE vecreal_<4, double> round(vecreal_<4, double> a)
	{
		return vecreal_<4, double> (::round(a[0]), ::round(a[1]), ::round(a[2]), ::round(a[3]));
	}

	#define IMPLEMENT_TRUNCATE_FUNCTION(type) \
	INLINE vecreal_<4, int> truncate(vecreal_<4, type> a) \
	{ \
		return vecint_<4, int> ((int)a[0], (int)a[1], (int)a[2], (int)a[3]); \
	}

	IMPLEMENT_TRUNCATE_FUNCTION(float);
	IMPLEMENT_TRUNCATE_FUNCTION(double);

	#undef IMPLEMENT_TRUNCATE_FUNCTION

	#define IMPLEMENT_TRUNCATE_FUNCTION_1(name, type) \
	INLINE vecint_<4, int> i##name(vecreal_<4, type> a) \
	{ \
		return vecint_<4, int> ((int)std::name(a[0]), (int)std::name(a[1]), \
			(int)std::name(a[2]), (int)std::name(a[3])); \
	}

	IMPLEMENT_TRUNCATE_FUNCTION_1(floor, float);
	IMPLEMENT_TRUNCATE_FUNCTION_1(floor, double);
	IMPLEMENT_TRUNCATE_FUNCTION_1(ceil,  float);
	IMPLEMENT_TRUNCATE_FUNCTION_1(ceil,  double);

	#undef IMPLEMENT_TRUNCATE_FUNCTION_1

	#define IMPLEMENT_CMP_FUNCTION(op, type) \
	INLINE vecbool_<4, int> operator op(vecreal_<4, type> a, vecreal_<4, type> b) \
	{ return vecbool_<4, int>(VECTOR_ARG(a, op, b)); }

	IMPLEMENT_CMP_FUNCTION(<, float);
	IMPLEMENT_CMP_FUNCTION(<, double);
	IMPLEMENT_CMP_FUNCTION(>, float);
	IMPLEMENT_CMP_FUNCTION(>, double);

	#undef IMPLEMENT_CMP_FUNCTION

	#define IMPLEMENT_CMP_FUNCTION(op, type) \
	INLINE vecbool_<4, int> operator op(vecreal_<4, type> a, type b) \
	{ return vecbool_<4, int>(SCALAR_ARG2(a, op, b)); }

	IMPLEMENT_CMP_FUNCTION(<, float);
	IMPLEMENT_CMP_FUNCTION(<, double);
	IMPLEMENT_CMP_FUNCTION(>, float);
	IMPLEMENT_CMP_FUNCTION(>, double);

	#undef IMPLEMENT_CMP_FUNCTION

	#define IMPLEMENT_CMP_FUNCTION(op, type) \
	INLINE vecbool_<4, int> operator op(type a, vecreal_<4, type> b) \
	{ return vecbool_<4, int>(SCALAR_ARG1(a, op, b)); }

	IMPLEMENT_CMP_FUNCTION(<, float);
	IMPLEMENT_CMP_FUNCTION(<, double);
	IMPLEMENT_CMP_FUNCTION(>, float);
	IMPLEMENT_CMP_FUNCTION(>, double);

	#undef IMPLEMENT_CMP_FUNCTION

	#define IMPLEMENT_MULTIPLE_NOT_FUNCTION(type) \
	INLINE vecreal_<4, type> multiple_not(vecbool_<4, int> a, vecreal_<4, type> b) \
	{ \
		return vecreal_<4, type>(a[0] ? 0.f : b[0], a[1] ? 0.f : b[1], \
			a[2] ? 0.f : b[2], a[3] ? 0.f : b[3]); \
	}

	IMPLEMENT_MULTIPLE_NOT_FUNCTION(float);
	IMPLEMENT_MULTIPLE_NOT_FUNCTION(double);

	#undef IMPLEMENT_MULTIPLE_NOT_FUNCTION

	#define IMPLEMENT_MULTIPLE_FUNCTION(type) \
	INLINE vecreal_<4, type> multiple(vecbool_<4, int> a, vecreal_<4, type> b) \
	{ \
		return vecreal_<4, type>(a[0] ? b[0] : 0.f, a[1] ? b[1] : 0.f, \
				a[2] ? b[2] : 0.f, a[3] ? b[3] : 0.f); \
	}

	IMPLEMENT_MULTIPLE_FUNCTION(float);
	IMPLEMENT_MULTIPLE_FUNCTION(double);

	#undef IMPLEMENT_MULTIPLE_FUNCTION

	INLINE vecreal_<4, float> operator*(vecint_<4, int> a, float b) { vecreal_<4, float> a_(a); return a_ *= b; }
	INLINE vecreal_<4, float> operator*(float b, vecint_<4, int> a) { vecreal_<4, float> a_(a); return a_ *= b; }
	INLINE vecreal_<4, double> operator*(vecint_<4, int> a, double b) { vecreal_<4, double> a_(a); return a_ *= b; }
	INLINE vecreal_<4, double> operator*(double b, vecint_<4, int> a) { vecreal_<4, double> a_(a); return a_ *= b; }

	/**
	* @brief Циклическое вращение компонент вектора на соответствующее число позиций.
	*/
	INLINE vecint_<4, int> ror0(vecint_<4, int> i) { return i; }
	INLINE vecint_<4, int> ror1(vecint_<4, int> i) { return vecint_<4, int>(i[3], i[0], i[1], i[2]); }
	INLINE vecint_<4, int> ror2(vecint_<4, int> i) { return vecint_<4, int>(i[2], i[3], i[0], i[1]); }
	INLINE vecint_<4, int> ror3(vecint_<4, int> i) { return vecint_<4, int>(i[1], i[2], i[3], i[0]); }

	INLINE vecint_<4, int> rol0(vecint_<4, int> i) { return i; }
	INLINE vecint_<4, int> rol1(vecint_<4, int> i) { return vecint_<4, int>(i[1], i[2], i[3], i[0]); }
	INLINE vecint_<4, int> rol2(vecint_<4, int> i) { return vecint_<4, int>(i[2], i[3], i[0], i[1]); }
	INLINE vecint_<4, int> rol3(vecint_<4, int> i) { return vecint_<4, int>(i[3], i[0], i[1], i[2]); }

	INLINE vecreal_<4, float> ror0(vecreal_<4, float> i) { return i; }
	INLINE vecreal_<4, float> ror1(vecreal_<4, float> i) { return vecreal_<4, float>(i[3], i[0], i[1], i[2]); }
	INLINE vecreal_<4, float> ror2(vecreal_<4, float> i) { return vecreal_<4, float>(i[2], i[3], i[0], i[1]); }
	INLINE vecreal_<4, float> ror3(vecreal_<4, float> i) { return vecreal_<4, float>(i[1], i[2], i[3], i[0]); }

	INLINE vecreal_<4, float> rol0(vecreal_<4, float> i) { return i; }
	INLINE vecreal_<4, float> rol1(vecreal_<4, float> i) { return vecreal_<4, float>(i[1], i[2], i[3], i[0]); }
	INLINE vecreal_<4, float> rol2(vecreal_<4, float> i) { return vecreal_<4, float>(i[2], i[3], i[0], i[1]); }
	INLINE vecreal_<4, float> rol3(vecreal_<4, float> i) { return vecreal_<4, float>(i[3], i[0], i[1], i[2]); }

	#undef VECTOR_OPERATOR
	#undef SCALAR_OPERATOR
	#undef MULTISCALAR_OPERATOR

	#undef VECTOR_ARG
	#undef SCALAR_ARG1
	#undef SCALAR_ARG2

	#define VECTOR_OPERATOR(v_, op, v) { v_[0] op v[0]; }
	#define SCALAR_OPERATOR(v_, op, v) { v_[0] op v; }
	#define MULTISCALAR_OPERATOR(v_, op, v0) { v_[0] op v0; }

	#define VECTOR_ARG(a, op, b) a[0] op b[0]
	#define SCALAR_ARG1(a, op, b) a op b[0]
	#define SCALAR_ARG2(a, op, b) a[0] op b
	/**
	 * @brief булевское значение для операций с плавающими числами
	 * Хранит все единичные биты (true) и все нулевые биты (false) для совместимости с SSE.
	 */
	template <>	class vecbool_<1, int>
	{
		int v_[1];

	public:
		enum { size = 1 };
		typedef int value_type;

		vecbool_() { SCALAR_OPERATOR(v_, =, 0); }

		vecbool_(bool v0)
		{
			v_[0] = v0 ? -1 : 0;
		}

		vecbool_(const vecbool_ &b) { VECTOR_OPERATOR(v_, =, b); }

		vecbool_ &operator=(vecbool_ b) { VECTOR_OPERATOR(v_, =, b); return *this; }

		value_type operator[](unsigned i) const { return v_[i]; }
		value_type &operator[](unsigned i) { return v_[i]; }

		operator bool() { return (bool)(v_[0]); }
		unsigned count() const { return std::abs(v_[0]); }
	};

	INLINE vecbool_<1, int> operator!(vecbool_<1, int> a) { return vecbool_<1, int>(SCALAR_ARG2(a, ==, 0)); }
	INLINE vecbool_<1, int> operator&&(vecbool_<1, int> a, vecbool_<1, int> b) { return vecbool_<1, int>(VECTOR_ARG(a, &, b)); }
	INLINE bool operator&&(vecbool_<1, int> a, bool b) { return (bool)a && b; }
	INLINE bool operator&&(bool a, vecbool_<1, int> b) { return a && (bool)b; }

	INLINE vecbool_<1, int> operator||(vecbool_<1, int> a, vecbool_<1, int> b) { return vecbool_<1, int>(VECTOR_ARG(a, |, b)); }
	INLINE bool operator==(vecbool_<1, int> a, bool b) { return (bool)a == b; }
	INLINE bool operator==(bool a, vecbool_<1, int> b) { return a == (bool)b; }

	INLINE vecbool_<1, int> multiple_not(vecbool_<1, int> a, vecbool_<1, int> b) { return (!a) && b; }
	INLINE vecbool_<1, int> multiple(vecbool_<1, int> a, vecbool_<1, int> b) { return a && b; }
	INLINE vecbool_<1, int> andnot(vecbool_<1, int> a, vecbool_<1, int> b) { return (!a) && b; }

	/**
	 * @brief класс-обертка для операций с целыми числами
	 */
	template <> class vecint_<1, int>
	{
		int v_[1];

	public:
		enum { size = 1 };
		typedef int value_type;

		vecint_() { SCALAR_OPERATOR(v_, =, 0); }
		vecint_(const vecint_ &i) { VECTOR_OPERATOR(v_, =, i); }
		vecint_(int i) { SCALAR_OPERATOR(v_, =, i); }
		vecint_(const vecbool_<1, int> &b) { VECTOR_OPERATOR(v_, = , b); }

		vecint_ &operator=(const vecint_ &i) { VECTOR_OPERATOR(v_, =, i); return *this; }
		vecint_ &operator=(int i) { SCALAR_OPERATOR(v_, =, i); return *this; }

		value_type operator[](unsigned i) const { return v_[i]; }
		value_type &operator[](unsigned i) { return v_[i]; }

		vecint_ &operator+=(vecint_ i) { VECTOR_OPERATOR(v_, +=, i); return *this; }
		vecint_ &operator-=(vecint_ i) { VECTOR_OPERATOR(v_, -=, i); return *this; }
	};

	INLINE vecint_<1, int> operator+(vecint_<1, int> a, vecint_<1, int> b) { return a += b; }
	INLINE vecint_<1, int> operator-(vecint_<1, int> a, vecint_<1, int> b) { return a -= b; }
	INLINE vecint_<1, int> operator+(vecint_<1, int> a, int b) { return a += vecint_<1, int>(b); }
	INLINE vecint_<1, int> operator-(vecint_<1, int> a, int b) { return a -= vecint_<1, int>(b); }
	INLINE vecint_<1, int> operator+(int a, vecint_<1, int> b) { vecint_<1, int> a_(a); return a_ += b; }
	INLINE vecint_<1, int> operator-(int a, vecint_<1, int> b) { vecint_<1, int> a_(a); return a_ -= b; }

	INLINE vecbool_<1, int> operator< (vecint_<1, int> a, vecint_<1, int> b) { return vecbool_<1, int>(VECTOR_ARG(a, <, b)); }
	INLINE vecbool_<1, int> operator> (vecint_<1, int> a, vecint_<1, int> b) { return vecbool_<1, int>(VECTOR_ARG(a, >, b)); }
	INLINE vecbool_<1, int> operator==(vecint_<1, int> a, vecint_<1, int> b) { return vecbool_<1, int>(VECTOR_ARG(a,==, b)); }

	INLINE vecbool_<1, int> operator< (vecint_<1, int> a, int b) { return vecbool_<1, int>(SCALAR_ARG2(a, <, b)); }
	INLINE vecbool_<1, int> operator> (vecint_<1, int> a, int b) { return vecbool_<1, int>(SCALAR_ARG2(a, >, b)); }
	INLINE vecbool_<1, int> operator==(vecint_<1, int> a, int b) { return vecbool_<1, int>(SCALAR_ARG2(a,==, b)); }

	INLINE vecbool_<1, int> operator< (int a, vecint_<1, int> b) { return vecbool_<1, int>(SCALAR_ARG1(a, <, b)); }
	INLINE vecbool_<1, int> operator> (int a, vecint_<1, int> b) { return vecbool_<1, int>(SCALAR_ARG1(a, >, b)); }
	INLINE vecbool_<1, int> operator==(int a, vecint_<1, int> b) { return vecbool_<1, int>(SCALAR_ARG1(a,==, b)); }

	INLINE vecint_<1, int> operator&(vecint_<1, int> a, vecint_<1, int> b) { return vecint_<1, int>(VECTOR_ARG(a, &, b)); }
	INLINE vecint_<1, int> operator|(vecint_<1, int> a, vecint_<1, int> b) { return vecint_<1, int>(VECTOR_ARG(a, |, b)); }
	INLINE vecint_<1, int> operator^(vecint_<1, int> a, vecint_<1, int> b) { return vecint_<1, int>(VECTOR_ARG(a, ^, b)); }

	INLINE vecint_<1, int> multiple_not(vecbool_<1, int> a, vecint_<1, int> b) { return vecint_<1, int>(!a) & b; }
	INLINE vecint_<1, int> multiple    (vecbool_<1, int> a, vecint_<1, int> b) { return vecint_<1, int>( a) & b; }

	/**
	 * Сдвиг влево каждого элемента вектора на заданное в параметре число бит.
	 * @note В отличие от SSE операций сдвиг для каждого элемента свой, а не единый для всех.
	 * @param a исходный вектор
	 * @param n вектор смещений
	 * @return вектор со смещенными компонентами
	 */
	INLINE vecint_<1, int> operator<<(vecint_<1, int> a, vecint_<1, int> n)
	{ return vecint_<1, int>(VECTOR_ARG(a, <<, n)); }

	/**
	* @brief Получение абсолютного значения смещения.
	* @note см. Уоррен. Алгоритмические трюки для программистов
	*/
	INLINE vecint_<1, int> abs(vecint_<1, int> a)
	{ return vecint_<1, int>(std::abs(a[0])); }

	/**
	 * @brief класс-обертка для операций с float числами
	 */
	template <typename T> class vecreal_<1, T>
	{
		T v_[1];

	public:

		enum { size = 1 };
		typedef T value_type;

		vecreal_() { SCALAR_OPERATOR(v_, =, 0.f); }
		template <typename S> vecreal_(S f) { SCALAR_OPERATOR(v_, =, f); }
		template <typename S> vecreal_(const vecreal_<1, S> &f) { VECTOR_OPERATOR(v_, =, f); }
		template <typename S> vecreal_(const vecint_<1, S> &f) { VECTOR_OPERATOR(v_, =, f); }
		template <typename S> vecreal_(const S *f) { VECTOR_OPERATOR(v_, =, f); }

		template <typename S> vecreal_ &operator=(S v) { SCALAR_OPERATOR(v_, =, v); return *this; }
		template <typename S> vecreal_ &operator=(vecreal_<1, S> v) { VECTOR_OPERATOR(v_, =, v); return *this; }

		value_type operator[](unsigned i) const { return v_[i]; }
		value_type &operator[](unsigned i) { return v_[i]; }

		template <typename S> vecreal_ &operator+=(vecreal_<1, S> f) { VECTOR_OPERATOR(v_, +=, f); return *this; }
		template <typename S> vecreal_ &operator-=(vecreal_<1, S> f) { VECTOR_OPERATOR(v_, -=, f); return *this; }
		template <typename S> vecreal_ &operator*=(vecreal_<1, S> f) { VECTOR_OPERATOR(v_, *=, f); return *this; }
		template <typename S> vecreal_ &operator/=(vecreal_<1, S> f) { VECTOR_OPERATOR(v_, /=, f); return *this; }

		template <typename S> vecreal_ &operator+=(vecint_<1, S> f) { VECTOR_OPERATOR(v_, +=, f); return *this; }
		template <typename S> vecreal_ &operator-=(vecint_<1, S> f) { VECTOR_OPERATOR(v_, -=, f); return *this; }
		template <typename S> vecreal_ &operator*=(vecint_<1, S> f) { VECTOR_OPERATOR(v_, *=, f); return *this; }
		template <typename S> vecreal_ &operator/=(vecint_<1, S> f) { VECTOR_OPERATOR(v_, /=, f); return *this; }

		template <typename S> vecreal_ &operator+=(S f) { SCALAR_OPERATOR(v_, +=, f); return *this; }
		template <typename S> vecreal_ &operator-=(S f) { SCALAR_OPERATOR(v_, -=, f); return *this; }
		template <typename S> vecreal_ &operator*=(S f) { SCALAR_OPERATOR(v_, *=, f); return *this; }
		template <typename S> vecreal_ &operator/=(S f) { SCALAR_OPERATOR(v_, /=, f); return *this; }

		template <typename S> void store(S *p) const { *p = v_[0]; }
	};

	template <unsigned N, typename T>
	INLINE vecreal_<N, T> operator+(vecreal_<N, T> a, vecreal_<N, T> b) { return a += b; }

	template <unsigned N, typename T>
	INLINE vecreal_<N, T> operator-(vecreal_<N, T> a, vecreal_<N, T> b) { return a -= b; }

	template <unsigned N, typename T>
	INLINE vecreal_<N, T> operator*(vecreal_<N, T> a, vecreal_<N, T> b) { return a *= b; }

	template <unsigned N, typename T>
	INLINE vecreal_<N, T> operator/(vecreal_<N, T> a, vecreal_<N, T> b) { return a /= b; }

	template <unsigned N, typename T>
	INLINE vecreal_<N, T> operator+(vecreal_<N, T> a, vecint_<N, int> b) { return a += b; }

	template <unsigned N, typename T>
	INLINE vecreal_<N, T> operator-(vecreal_<N, T> a, vecint_<N, int> b) { return a -= b; }

	template <unsigned N, typename T>
	INLINE vecreal_<N, T> operator*(vecreal_<N, T> a, vecint_<N, int> b) { return a *= b; }

	template <unsigned N, typename T>
	INLINE vecreal_<N, T> operator/(vecreal_<N, T> a, vecint_<N, int> b) { return a /= b; }

	template <unsigned N, typename T>
	INLINE vecreal_<N, T> operator+(vecint_<N, int> a, vecreal_<N, T> b) { return b += a; }

	template <unsigned N, typename T>
	INLINE vecreal_<N, T> operator-(vecint_<N, int> a, vecreal_<N, T> b) { vecreal_<N, T> a_(a); return a_ - b; }

	template <unsigned N, typename T>
	INLINE vecreal_<N, T> operator*(vecint_<N, int> a, vecreal_<N, T> b) { return b *= a; }

	template <unsigned N, typename T>
	INLINE vecreal_<N, T> operator/(vecint_<N, int> a, vecreal_<N, T> b) { vecreal_<N, T> a_(a); return a_ / b; }

	template <unsigned N, typename T>
	INLINE vecreal_<N, T> operator+(vecreal_<N, T> a, T b) { return a += b; }

	template <unsigned N, typename T>
	INLINE vecreal_<N, T> operator-(vecreal_<N, T> a, T b) { return a -= b; }

	template <unsigned N, typename T>
	INLINE vecreal_<N, T> operator*(vecreal_<N, T> a, T b) { return a *= b; }

	template <unsigned N, typename T>
	INLINE vecreal_<N, T> operator/(vecreal_<N, T> a, T b) { return a /= b; }

	template <unsigned N, typename T, typename S>
	INLINE vecreal_<N, T> operator+(S a, vecreal_<N, T> b) { return b += a; }

	template <unsigned N, typename T, typename S>
	INLINE vecreal_<N, T> operator-(S a, vecreal_<N, T> b) { vecreal_<N, T> a_(a); return a_ -= b; }

	template <unsigned N, typename T, typename S>
	INLINE vecreal_<N, T> operator*(S a, vecreal_<N, T> b) { return b *= a; }

	template <unsigned N, typename T, typename S>
	INLINE vecreal_<N, T> operator/(S a, vecreal_<N, T> b) { vecreal_<N, T> a_(a); return a_ /= b; }

	template <unsigned N, typename T>
	INLINE vecreal_<N, T> sqr(vecreal_<N, T> a) { return a * a; }

	template <unsigned N, typename T>
	INLINE vecreal_<N, T> cube(vecreal_<N, T> a) { return a * a * a; }

	#define IMPLEMENT_4FUNCTION_2(name, type) \
	INLINE vecreal_<1, type> name(vecreal_<1, type> a, vecreal_<1, type> b) \
	{ \
		return vecreal_<1, type>(std::name(a[0], b[0])); \
	}

	IMPLEMENT_4FUNCTION_2(max, float);
	IMPLEMENT_4FUNCTION_2(max, double);
	IMPLEMENT_4FUNCTION_2(min, float);
	IMPLEMENT_4FUNCTION_2(min, double);

	#undef IMPLEMENT_4FUNCTION_2

	#define IMPLEMENT_4FUNCTION_1(name, type) \
	INLINE vecreal_<1, type> name(vecreal_<1, type> a) \
	{ \
		return vecreal_<1, type> (std::name(a[0])); \
	}

	IMPLEMENT_4FUNCTION_1(sqrt,  float);
	IMPLEMENT_4FUNCTION_1(sqrt,  double);
	IMPLEMENT_4FUNCTION_1(floor, float);
	IMPLEMENT_4FUNCTION_1(floor, double);
	IMPLEMENT_4FUNCTION_1(ceil,  float);
	IMPLEMENT_4FUNCTION_1(ceil,  double);

	#undef IMPLEMENT_4FUNCTION_1

	template <typename T>
	INLINE vecreal_<1, T> exp(vecreal_<1, T> a)
	{
		return vecreal_<1, T>(::exp(a[0]));
	}

	template <typename T>
	INLINE vecreal_<1, T> erf(vecreal_<1, T> a)
	{
		return vecreal_<1, T>(::erf(a[0]));
	}

	template <typename T>
	INLINE vecreal_<1, T> round(vecreal_<1, T> a)
	{
		return vecreal_<1, T> (roundf(a[0]));
	}

	#define IMPLEMENT_TRUNCATE_FUNCTION(type) \
	INLINE vecreal_<1, int> truncate(vecreal_<1, type> a) \
	{ \
		return vecint_<1, int> ((int)a[0]); \
	}

	IMPLEMENT_TRUNCATE_FUNCTION(float);
	IMPLEMENT_TRUNCATE_FUNCTION(double);

	#undef IMPLEMENT_TRUNCATE_FUNCTION

	#define IMPLEMENT_TRUNCATE_FUNCTION_1(name, type) \
	INLINE vecint_<1, int> i##name(vecreal_<1, type> a) \
	{ \
		return vecint_<1, int> ((int)std::name(a[0])); \
	}

	IMPLEMENT_TRUNCATE_FUNCTION_1(floor, float);
	IMPLEMENT_TRUNCATE_FUNCTION_1(floor, double);
	IMPLEMENT_TRUNCATE_FUNCTION_1(ceil,  float);
	IMPLEMENT_TRUNCATE_FUNCTION_1(ceil,  double);

	#undef IMPLEMENT_TRUNCATE_FUNCTION_1

	#define IMPLEMENT_CMP_FUNCTION(op, type) \
	INLINE vecbool_<1, int> operator op(vecreal_<1, type> a, vecreal_<1, type> b) \
	{ return vecbool_<1, int>(VECTOR_ARG(a, op, b)); }

	IMPLEMENT_CMP_FUNCTION(<, float);
	IMPLEMENT_CMP_FUNCTION(<, double);
	IMPLEMENT_CMP_FUNCTION(>, float);
	IMPLEMENT_CMP_FUNCTION(>, double);

	#undef IMPLEMENT_CMP_FUNCTION

	#define IMPLEMENT_CMP_FUNCTION(op, type) \
	INLINE vecbool_<1, int> operator op(vecreal_<1, type> a, type b) \
	{ return vecbool_<1, int>(SCALAR_ARG2(a, op, b)); }

	IMPLEMENT_CMP_FUNCTION(<, float);
	IMPLEMENT_CMP_FUNCTION(<, double);
	IMPLEMENT_CMP_FUNCTION(>, float);
	IMPLEMENT_CMP_FUNCTION(>, double);

	#undef IMPLEMENT_CMP_FUNCTION

	#define IMPLEMENT_CMP_FUNCTION(op, type) \
	INLINE vecbool_<1, int> operator op(type a, vecreal_<1, type> b) \
	{ return vecbool_<1, int>(SCALAR_ARG1(a, op, b)); }

	IMPLEMENT_CMP_FUNCTION(<, float);
	IMPLEMENT_CMP_FUNCTION(<, double);
	IMPLEMENT_CMP_FUNCTION(>, float);
	IMPLEMENT_CMP_FUNCTION(>, double);

	#undef IMPLEMENT_CMP_FUNCTION

	#define IMPLEMENT_MULTIPLE_NOT_FUNCTION(type) \
	INLINE vecreal_<1, type> multiple_not(vecbool_<1, int> a, vecreal_<1, type> b) \
	{ \
		return vecreal_<1, type>(a[0] ? 0.f : b[0]); \
	}

	IMPLEMENT_MULTIPLE_NOT_FUNCTION(float);
	IMPLEMENT_MULTIPLE_NOT_FUNCTION(double);

	#undef IMPLEMENT_MULTIPLE_NOT_FUNCTION

	#define IMPLEMENT_MULTIPLE_FUNCTION(type) \
	INLINE vecreal_<1, type> multiple(vecbool_<1, int> a, vecreal_<1, type> b) \
	{ \
		return vecreal_<1, type>(a[0] ? b[0] : 0.f); \
	}

	IMPLEMENT_MULTIPLE_FUNCTION(float);
	IMPLEMENT_MULTIPLE_FUNCTION(double);

	#undef IMPLEMENT_MULTIPLE_FUNCTION

	INLINE vecreal_<1, float> operator*(vecint_<1, int> a, float b) { vecreal_<1, float> a_(a); return a_ *= b; }
	INLINE vecreal_<1, float> operator*(float b, vecint_<1, int> a) { vecreal_<1, float> a_(a); return a_ *= b; }
	INLINE vecreal_<1, double> operator*(vecint_<1, int> a, double b) { vecreal_<1, double> a_(a); return a_ *= b; }
	INLINE vecreal_<1, double> operator*(double b, vecint_<1, int> a) { vecreal_<1, double> a_(a); return a_ *= b; }

	/**
	* @brief Циклическое вращение компонент вектора на соответствующее число позиций.
	*/
	INLINE vecint_<1, int> ror0(vecint_<1, int> i) { return i; }
	INLINE vecint_<1, int> rol0(vecint_<1, int> i) { return i; }
	INLINE vecreal_<1, float> ror0(vecreal_<1, float> i) { return i; }
	INLINE vecreal_<1, float> rol0(vecreal_<1, float> i) { return i; }

	INLINE int summarize(vecint_<1, int> v) { return v[0]; }
	INLINE int summarize(vecint_<4, int> v) { return v[0] + v[1] + v[2] + v[3]; }
	INLINE float summarize(vecreal_<1, float> v) { return v[0]; }
	INLINE float summarize(vecreal_<4, float> v) { return v[0] + v[1] + v[2] + v[3]; }
	INLINE double summarize(vecreal_<1, double> v) { return v[0]; }
	INLINE double summarize(vecreal_<4, double> v) { return v[0] + v[1] + v[2] + v[3]; }


	template <typename T>
	INLINE T summarize(vecreal_<1, T> v) { return v[0]; }

	#undef VECTOR_OPERATOR
	#undef SCALAR_OPERATOR
	#undef MULTISCALAR_OPERATOR

	#undef VECTOR_ARG
	#undef SCALAR_ARG1
	#undef SCALAR_ARG2

	#endif
	/**
	* @brief Преобразование в строку для отладки.
	* @note Ипользуется порядок компонент, аналогичный порядку битов, то есть справа налево.
	*/
	INLINE std::string make_string(const vecreal_<4, float> &f)
	{
		float *p = (float*)&f;
		return make_string("{%e %e %e %e}", p[0], p[1], p[2], p[3]);
	}
	INLINE std::string make_string(const vecint_<4, int> &f)
	{
		int *p = (int*)&f;
		return make_string("{%d %d %d %d}", p[0], p[1], p[2], p[3]);
	}

	INLINE std::string make_string(const vecbool_<4, int> &f)
	{
		int *p = (int*)&f;
		return make_string("{%d %d %d %d}", p[0], p[1], p[2], p[3]);
	}

	/**
	* @brief Преобразование в строку для отладки.
	* @note Ипользуется порядок компонент, аналогичный порядку битов, то есть справа налево.
	*/
	INLINE std::string make_string(const vecreal_<1, float> &f)
	{
		float *p = (float*)&f;
		return make_string("{%e}", p[0]);
	}
	INLINE std::string make_string(const vecint_<1, int> &f)
	{
		int *p = (int*)&f;
		return make_string("{%d}", p[0]);
	}

	INLINE std::string make_string(const vecbool_<1, int> &f)
	{
		int *p = (int*)&f;
		return make_string("{%d}", p[0]);
	}

	/**
	 * Возращает число векторов (например SSE вектор), для полного размещения заданного числа объектов.
	 * @note при неиспользовании SSE возвращатся исходное число элементов
	 * @param count число объектов
	 * @param n число объектов, которые помещаются в вектор
	 * @return число векторов (SSE)
	 */
	template <typename T> INLINE unsigned velement_count(unsigned count)
	{ return (count - 1) / T::size + 1; }

}
#endif
