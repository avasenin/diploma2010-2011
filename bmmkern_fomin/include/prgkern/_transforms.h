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
#ifndef _TRANSFORMS__0077A726_E6E3_58a1_C16D_CE436AB90500__H
#define _TRANSFORMS__0077A726_E6E3_58a1_C16D_CE436AB90500__H
#include "prgkern/_prgconfig.h"
#include "prgkern/_assert.h"
#include "prgkern/_type.h"
#include "prgkern/_index.h"
#include "prgkern/_v3dense.h"

#ifdef USE_FFTW_LIBRARY
	#include <complex>
	#include <fftw3.h>
#endif

namespace prgkern
{

	// the mesh sizes which are support the standard FFTW implementation (primes - 2, 3, 5, 7)
	const unsigned dfft_sizes[] = {
		   2,     3,     4,     5,     6,     7,     8,     9,    10,    12,    14,    15,
		  16,    18,    20,    21,    24,    25,    27,    28,    30,    32,    35,    36,
		  40,    42,    45,    48,    49,    50,    54,    56,    60,    63,    64,    70,
		  72,    75,    80,    81,    84,    90,    96,    98,   100,   105,   108,   112,
		 120,   125,   126,   128,   135,   140,   144,   147,   150,   160,   162,   168,
		 175,   180,   189,   192,   196,   200,   210,   216,   224,   225,   240,   243,
		 245,   250,   252,   256,   270,   280,   288,   294,   300,   315,   320,   324,
		 336,   343,   350,   360,   375,   378,   384,   392,   400,   405,   420,   432,
		 441,   448,   450,   480,   486,   490,   500,   504,   512,   525,   540,   560,
		 567,   576,   588,   600,   625,   630,   640,   648,   672,   675,   686,   700,
		 720,   729,   735,   750,   756,   768,   784,   800,   810,   840,   864,   875,
		 882,   896,   900,   945,   960,   972,   980,  1000,  1008,  1024,  1029,  1050,
		1080,  1120,  1125,  1134,  1152,  1176,  1200,  1215,  1225,  1250,  1260,  1280,
		1296,  1323,  1344,  1350,  1372,  1400,  1440,  1458,  1470,  1500,  1512,  1536,
		1568,  1575,  1600,  1620,  1680,  1701,  1715,  1728,  1750,  1764,  1792,  1800,
		1875,  1890,  1920,  1944,  1960,  2000,  2016,  2025,  2048,  2058,  2100,  2160,
		2187,  2205,  2240,  2250,  2268,  2304,  2352,  2400,  2401,  2430,  2450,  2500,
		2520,  2560,  2592,  2625,  2646,  2688,  2700,  2744,  2800,  2835,  2880,  2916,
		2940,  3000,  3024,  3072,  3087,  3125,  3136,  3150,  3200,  3240,  3360,  3375,
		3402,  3430,  3456,  3500,  3528,  3584,  3600,  3645,  3675,  3750,  3780,  3840,
		3888,  3920,  3969,  4000,  4032,  4050,  4116,  4200,  4320,  4374,  4375,  4410,
		4480,  4500,  4536,  4608,  4704,  4725,  4800,  4802,  4860,  4900,  5000,  5040,
		5103,  5120,  5145,  5184,  5250,  5292,  5376,  5400,  5488,  5600,  5625,  5670,
		5760,  5832,  5880,  6000,  6048,  6075,  6125,  6144,  6174,  6250,  6272,  6300,
		6400,  6480,  6615,  6720,  6750,  6804,  6860,  6912,  7000,  7056,  7168,  7200,
		7203,  7290,  7350,  7500,  7560,  7680,  7776,  7840,  7875,  7938,  8000,  8064,
		8100,  8232,  8400,  8505,  8575,  8640,  8748,  8750,  8820,  8960,  9000,  9072,
		9216,  9261,  9375,  9408,  9450,  9600,  9604,  9720,  9800, 10000, 10080, 10125
	};
	/**
	* @brief calculate the number of bits needed to represent an integer n
	*/
	inline unsigned bitsizeof(long long x)
	{
		unsigned n = 0;
		if (x >> 32) { n += 32; x >>= 32; }
		if (x >> 16) { n += 16; x >>= 16; }
		if (x >> 8 ) { n +=  8; x >>=  8; }
		if (x >> 4 ) { n +=  4; x >>=  4; }
		if (x >> 2 ) { n +=  2; x >>=  2; }
		if (x >> 1 ) { n +=  1; x >>=  1; }
		return n + x;
	}

	/**
	* @brief calculate the number of bits needed to represent an integer n
	*/
	INLINE unsigned bitsizeof(unsigned x)
	{
		unsigned n = 0;
		if (x >> 16) { n += 16; x >>= 16; }
		if (x >> 8 ) { n +=  8; x >>=  8; }
		if (x >> 4 ) { n +=  4; x >>=  4; }
		if (x >> 2 ) { n +=  2; x >>=  2; }
		if (x >> 1 ) { n +=  1; x >>=  1; }
		return n + x;
	}

	/**
	* @brief calculate integer log2(integer)
	* @note result is cut to the upper(!) bound
	*/
	INLINE unsigned ilog2(unsigned n)
	{
		assert(_GT((unsigned)n, (unsigned)0));
		return bitsizeof(--n);
	}

	/**
	* @brief Fast Discrete Hadamard Transform
	* @param n size must be a power of two
	* @note [out] vector must be normalized by 1./sqrt(n)
	*/
	template <typename _Real> INLINE void dfht(unsigned n, _Real *v)
	{
		unsigned m = ilog2(n);
		assert(_EQ((unsigned)(1 << m), (unsigned)n));

		unsigned l = 1;
		for (unsigned i=0; i<m; ++i)
		{
			n >>= 1;
			unsigned i__ = 0;
			for (unsigned k=0; k<n; ++k)
			{
				for (unsigned j=0; j<l; ++j)
				{
					_Real t = v[i__ + j];
					v[i__ + j] += v[i__ + j + l];
					v[i__ + j + l] = t - v[i__ + j + l];
				}
				i__ += l + l;
			}
			l <<= 1;
		}
	}

	/**
	* @brief Fast Discrete Walsh Hadamard Transform
	* @param n size must be a power of two
	* @note [out] vector must be normalized by 1./sqrt(n)
	*/
	template <typename _Real> INLINE void dfwht(unsigned n, _Real *v)
	{
		fdht(n, v); // make fast Hadamard transform

		// restore ordering of elements
		for (unsigned i=0, j=0, k; i<n-1; ++i)
		{
			if (i < j) std::swap(v[i], v[j]);
			k = n >> 1;
			while (k <= j) { j -= k; k >>= 1; }
			j += k;
		}
	}

	/**
	* @brief Fast 2D Discrete Hadamard Transform
	* @param (n1,n2) sizes must be a power of two
	* @note [out] must be normalized by 1./sqrt(n1*n2) and transpose
	*/
	template <typename _Real>
	INLINE void dfht2d(unsigned n1, unsigned n2, _Real *m)
	{
		for (unsigned i=0, i__=0; i<n1; ++i, i__+=n2)
			fdht(n2, &m[i__]);

		std::vector<_Real> v(n1);
		for (unsigned i=0; i<n2; ++i)
		{
			for (unsigned k=0, k__=i; k<n1; ++k, k__+=n2) v[k] = m[k__];
				// copy column to compress data in memory for efficiency
			fdht(n1, &v[0]);
			for (unsigned k=0, k__=i; k<n1; ++k, k__+=n2) m[k__] = v[k];
				// restore column
		}
	}

	/**
	* @brief Fast Discrete Walsh Hadamard Transform
	* @param (n1,n2) sizes must be a power of two
	* @note [out] must be normalized by 1./sqrt(n1*n2) and transpose
	*/
	template <typename _Real>
	INLINE void dfwht2d(unsigned n1, unsigned n2, _Real *m)
	{
		for (unsigned i=0, i__=0; i<n1; ++i, i__+=n2)
			fdwht(n2, &m[i__]);

		std::vector<_Real> v(n1);
		for (unsigned i=0; i<n2; ++i)
		{
			for (unsigned k=0, k__=i; k<n1; ++k, k__+=n2) v[k] = m[k__];
				// copy column to compress data in memory for efficiency
			fdwht(n1, &v[0]);
			for (unsigned k=0, k__=i; k<n1; ++k, k__+=n2) m[k__] = v[k];
				// restore column
		}
	}

#ifdef USE_FFTW_LIBRARY

	#define FFTW3_CONNECT_BRIDGE_MACRO(_N, _Real, symbol, class_name) \
	protected: \
		typedef fftw##symbol##_complex           fft_complex; \
		typedef _Real                            fft_real; \
		typedef fftw##symbol##_plan              fft_plan; \
		typedef index_<_N, int>                  fft_index; \
		typedef vdense_<_N, fft_real>            fft_v3dense; \
		typedef std::pair<fft_index, fft_index>  fft_pair_index; \
		\
		fft_complex *fft_icdata_; \
		fft_complex *fft_ocdata_; \
		fft_real *fft_irdata_; \
		fft_real *fft_ordata_; \
		fft_plan fft_plan_; \
		fft_index fft_n_; \
		fft_index fft_kdown_; \
		fft_index fft_kup_; \
		fft_index fft_ksz_; \
		fft_v3dense fft_knorm_; \
		unsigned fft_size_; \
		unsigned fft_csize_; \
		\
		class_name(const fft_index &ndx) : fft_icdata_(NULL), fft_ocdata_(NULL), \
			fft_irdata_(NULL), fft_ordata_(NULL), fft_plan_(NULL), \
			fft_n_(ndx), fft_size_(1), fft_csize_(1) \
		{ \
			for (unsigned i=0; i<_N; i++) \
			{ \
				fft_kdown_[i] = ( 1 - fft_n_[i] ) / 2; \
				fft_kup_[i] = fft_n_[i] / 2 + 1; \
				fft_ksz_[i] = fft_n_[i]; \
				fft_size_ *= fft_n_[i]; \
				fft_knorm_[i] = M_2PI / fft_n_[i]; \
			} \
			fft_csize_ = (fft_size_ / fft_n_[_N-1]) * (fft_n_[_N-1] / 2 + 1); \
			fft_ksz_[_N-1] = fft_n_[_N-1] / 2 + 1; \
		} \
		~class_name() { \
			if (fft_icdata_) fft_free(fft_icdata_); \
			if (fft_ocdata_) fft_free(fft_ocdata_); \
			if (fft_irdata_) fft_free(fft_irdata_); \
			if (fft_ordata_) fft_free(fft_ordata_); \
			fft_destroy_plan(); \
		 } \
		fft_pair_index xbox() const { return _Pair(fft_index(), fft_n_); } \
		fft_pair_index kbox() const { return _Pair(fft_kdown_, fft_kup_); } \
		fft_real norm() const { return sqrt(1. / fft_size_); } \
		\
		void *fft_malloc(size_t n) { return fftw##symbol##_malloc(n); } \
		void fft_free(void *p) { fftw##symbol##_free(p); } \
		\
		fft_plan fft_plan_dft(fft_complex *in, fft_complex *out, \
			int sign=FFTW_FORWARD, unsigned flags=FFTW_ESTIMATE) \
		{ return fftw##symbol##_plan_dft(_N, &fft_n_[0], in, out, sign, flags); } \
		\
		fft_plan fft_plan_dft(_Real *in, fft_complex *out, unsigned flags=FFTW_ESTIMATE) \
		{ return fftw##symbol##_plan_dft_r2c(_N, &fft_n_[0], in, out, flags); } \
		\
		fft_plan fft_plan_dft(fft_complex *in, _Real *out, unsigned flags=FFTW_ESTIMATE) \
		{ return fftw##symbol##_plan_dft_c2r(_N, &fft_n_[0], in, out, flags); } \
		\
		fft_plan fft_plan_dft(int howmany, \
			fft_complex *in, const int *inembed, int istride, int idist, \
			fft_complex *out, const int *onembed, int ostride, int odist, \
			int sign=FFTW_FORWARD, unsigned flags=FFTW_ESTIMATE) \
		{ \
			return fftw##symbol##_plan_many_dft(_N, &fft_n_[0], howmany, in, inembed, istride, idist, \
			out, onembed, ostride, odist, sign, flags); \
		} \
		\
		fft_plan fft_plan_dft(int howmany, \
			_Real *in, const int *inembed, int istride, int idist, \
			fft_complex *out, const int *onembed, int ostride, int odist, \
			unsigned flags=FFTW_ESTIMATE) \
		{ \
			return fftw##symbol##_plan_many_dft_r2c(_N, &fft_n_[0], howmany, in, inembed, istride, idist, \
			out, onembed, ostride, odist, flags); \
		} \
		\
		fft_plan fft_plan_dft(int howmany, \
			fft_complex *in, const int *inembed, int istride, int idist, \
			_Real *out, const int *onembed, int ostride, int odist, \
			unsigned flags=FFTW_ESTIMATE) \
		{ \
			return fftw##symbol##_plan_many_dft_c2r(_N, &fft_n_[0], howmany, in, inembed, istride, idist, \
			out, onembed, ostride, odist, flags); \
		} \
		\
		void fft_destroy_plan() { if (fft_plan_) fftw##symbol##_destroy_plan(fft_plan_); } \
		\
		void fft_execute_dft(fft_complex *in, fft_complex *out) \
		{ fftw##symbol##_execute_dft(fft_plan_, in, out); } \
		\
		void fft_execute_dft(_Real *in, fft_complex *out) \
		{ fftw##symbol##_execute_dft_r2c(fft_plan_, in, out); } \
		\
		void fft_execute_dft(fft_complex *in, _Real *out) \
		{ fftw##symbol##_execute_dft_c2r(fft_plan_, in, out); } \

	#define USING_FFTW3_CONNECT_DEFINES(_N, _Real) \
		typedef typename fftw3_connect<_N, _Real>::fft_complex  fft_complex; \
		typedef typename fftw3_connect<_N, _Real>::fft_real     fft_real; \
		typedef typename fftw3_connect<_N, _Real>::fft_plan     fft_plan; \
		typedef typename fftw3_connect<_N, _Real>::fft_index    fft_index; \
		\
		using fftw3_connect<_N, _Real>::fft_icdata_; \
		using fftw3_connect<_N, _Real>::fft_ocdata_; \
		using fftw3_connect<_N, _Real>::fft_irdata_; \
		using fftw3_connect<_N, _Real>::fft_ordata_; \
		using fftw3_connect<_N, _Real>::fft_plan_; \
		using fftw3_connect<_N, _Real>::fft_n_; \
		using fftw3_connect<_N, _Real>::fft_kdown_; \
		using fftw3_connect<_N, _Real>::fft_kup_; \
		using fftw3_connect<_N, _Real>::fft_ksz_; \
		using fftw3_connect<_N, _Real>::fft_knorm_; \
		using fftw3_connect<_N, _Real>::fft_size_; \
		using fftw3_connect<_N, _Real>::fft_csize_; \
		\
		using fftw3_connect<_N, _Real>::fft_malloc; \
		using fftw3_connect<_N, _Real>::fft_free; \
		using fftw3_connect<_N, _Real>::fft_plan_dft; \
		using fftw3_connect<_N, _Real>::fft_destroy_plan; \
		using fftw3_connect<_N, _Real>::fft_execute_dft; \
		using fftw3_connect<_N, _Real>::norm; \
		using fftw3_connect<_N, _Real>::xbox; \
		using fftw3_connect<_N, _Real>::kbox; \

	/**
	* @brief FTTW3 (-lfftw3, -lfftw3f, -lfftw3l ) library connect class
	*/
	template <unsigned _N, typename _Real> struct fftw3_connect;

	template <unsigned _N> struct fftw3_connect<_N, float>
	{ FFTW3_CONNECT_BRIDGE_MACRO(_N, float, f, fftw3_connect) };

	template <unsigned _N> struct fftw3_connect<_N, double>
	{ FFTW3_CONNECT_BRIDGE_MACRO(_N, double, ,fftw3_connect) };

	template <unsigned _N> struct fftw3_connect<_N, long double>
	{ FFTW3_CONNECT_BRIDGE_MACRO(_N, long double, l, fftw3_connect) };

	#define IMPORT_RUN_EXPORT(vi, vo) \
		import_(vi); execute_(); export_(vo);

	/**
	* @brief Fast Discrete Fourier Transform
	* The implementation is built upon FFTW (version 3.0.0 or higher)
	* @note FFTW-based implementation is the fastest for powers of two.
	* Furthermore, the second time you call the routine with the same size,
	* the calculation is much faster due to many things were calculated and
	* stored the first time the routine was called.
	* Achieving maximum runtime efficiency with the FFTW library on some
	* computer architectures requires that data are stored in the memory with
	* a special alignment (to 16-byte boundaries).
	* @note [out] must be normalized by 1./sqrt(2*n)
	*/
	template <unsigned _N, typename _From, typename _To> class Dfft;

	template <unsigned _N, typename _Real>
	class Dfft<_N, std::complex<_Real>, std::complex<_Real> >
	: public fftw3_connect<_N, _Real>
	{
		typedef fftw3_connect<_N, _Real>  _Base;
		USING_FFTW3_CONNECT_DEFINES(_N, _Real);

	public:

		/**
		* @brief конструирует объект и делает предварительные вычисления для эффективности
		* @param ndx мультииндекс, описывающий размерности входного массива
		* @param direction направление преобразования { DFFT_FORWARD, DFFT_BACKWARD }
		*/
		Dfft(const fft_index &ndx, int direction) : _Base(ndx)
		{
			fft_icdata_ = (fft_complex *)fft_malloc(sizeof(fft_complex) * fft_size_);
			fft_ocdata_ = (fft_complex *)fft_malloc(sizeof(fft_complex) * fft_size_);
			fft_plan_ = fft_plan_dft(fft_icdata_, fft_ocdata_, direction, FFTW_ESTIMATE);
		}

		/**
		* @brief вычисляет число элементов в выходном массиве
		* @return число элементов во входном массиве
		*/
		unsigned isize() const { return fft_size_; }
		unsigned osize() const { return fft_size_; }

		/**
		* @brief перекопирует данные и выполняет преобразование
		* @param vo - выходной массив
		* @param vi - входной массив
		*/
		void operator()(std::complex<_Real> *vo, const std::complex<_Real> *vi)
		{ IMPORT_RUN_EXPORT(vi, vo) }

	protected:

		void import_(const std::complex<_Real> *vi)
		{
			for (unsigned i=0; i<fft_size_; i++)
			{
				fft_icdata_[i][0] = vi[i].real();
				fft_icdata_[i][1] = vi[i].imag();
			}
		}

		void execute_() { fft_execute_dft(fft_icdata_, fft_ocdata_); }

		void export_(std::complex<_Real> *vo)
		{
			for (unsigned i=0; i<fft_size_; i++)
			{
				vo[i] = std::complex<_Real>(fft_ocdata_[i][0], fft_ocdata_[i][1]);
			}
		}

	};

	template <unsigned _N, typename _Real>
	class Dfft<_N, _Real, std::complex<_Real> > : public fftw3_connect<_N, _Real>
	{

		typedef fftw3_connect<_N, _Real>  _Base;
		USING_FFTW3_CONNECT_DEFINES(_N, _Real);

	public:

		/**
		* @brief конструирует объект и делает предварительные вычисления для эффективности
		* @param ndx мультииндекс, описывающий размерности входного массива
		*/
		Dfft(const fft_index &ndx) : _Base(ndx)
		{
			fft_irdata_ = (fft_real *)fft_malloc(sizeof(fft_real) * fft_size_);
			fft_ocdata_ = (fft_complex *)fft_malloc(sizeof(fft_complex) * fft_csize_);
			fft_plan_ = fft_plan_dft(fft_irdata_, fft_ocdata_, FFTW_ESTIMATE);
		}

		/**
		* @brief вычисляет число элементов во входном массиве
		* @return число элементов во входном массиве
		*/
		unsigned isize() const { return fft_size_; }
		unsigned osize() const { return fft_csize_; }

		/**
		* @brief перекопирует данные и выполняет преобразование
		* @param vo - выходной массив
		* @param vi - входной массив
		*/
		void operator()(std::complex<_Real> *vo, const _Real *vi)
		{ IMPORT_RUN_EXPORT(vi, vo) }

		/**
		* @brief make dfft (use another order of params)
		* @note this variant used by mesh to avoid creation of 2 functions with different params order
		*/
		void operator()(const _Real *vi, std::complex<_Real> *vo) { operator()(vo, vi); }

		/**
		* @brief находит позицию элемента с индексом k в k-массиве
		* @note в случае, если элемент попадает в отсутствующую половину массива,
		*  возвращается индекс сопряженного элемента и устанавливается флаг=false
		*  При поиске учитывется цикличность вектора k
		* @param flag[out] - true истинный элемент, false - требуется сопряжение
		* @param k - индекс k-вектора
		* @return позиция элемента в k-массиве
		*/
		unsigned find_pos(bool &flag, const fft_index &k) const
		{
			fft_index ndx__ = mod(k, fft_n_); // сдвинем в область определения
			flag = true;
			if (ndx__[_N-1] >= fft_ksz_[_N-1]) // проверим знак
			{
				flag = false;
				ndx__ = mod(-ndx__, fft_n_);
			}
			return memory_offset(ndx__, fft_ksz_);
		}

		/**
		* @brief находит индекс k элемента при известной его позиции
		* @note всегда возвращается k-индекс прямого элемента
		* @param n позиция элемента в k-массиве
		* @return индекс k-вектора
		*/
		fft_index find_pos(unsigned n) const
		{
			fft_index k = memory_offset(n, fft_ksz_);
			for (unsigned i=0; i<_N-1; i++)
				if (k[i] >= fft_kup_[i]) k[i] -= fft_ksz_[i];
			return k;
		}

		/**
		* @brief рассчитывает ненормированный k вектор элемента при известной его позиции
		* Нормировка делается домножением каждой компоненты на множитель 2*pi/h[i].
		* @param n позиция элемента в k-массиве
		* @return ненормированный k-вектор
		*/
		void extract_k(vdense_<3, _Real> &k, unsigned n) const
		{
			fft_index kndx = memory_offset(n, fft_ksz_);
			k[0] = (kndx[0] >= fft_kup_[0] ? kndx[0] - fft_ksz_[0] : kndx[0]) * fft_knorm_[0];
			k[1] = (kndx[1] >= fft_kup_[1] ? kndx[1] - fft_ksz_[1] : kndx[1]) * fft_knorm_[1];
			k[2] = kndx[2] * fft_knorm_[2];
		}

		void extract_k(vdense_<2, _Real> &k, unsigned n) const
		{
			fft_index kndx = memory_offset(n, fft_ksz_);
			k[0] = (kndx[0] >= fft_kup_[0] ? kndx[0] - fft_ksz_[0] : kndx[0]) * fft_knorm_[0];
			k[1] = kndx[1] * fft_knorm_[1];
		}

		void extract_k(vdense_<1, _Real> &k, unsigned n) const
		{
			fft_index kndx = memory_offset(n, fft_ksz_);
			k[0] = kndx[0] * fft_knorm_[0];
		}

	protected:

		void import_(const _Real *vi)
		{ for (unsigned i=0; i<fft_size_; i++) fft_irdata_[i] = vi[i]; }

		void execute_() { fft_execute_dft(fft_irdata_, fft_ocdata_); }

		void export_(std::complex<_Real> *vo)
		{
			for (unsigned i=0; i<fft_csize_; i++)
				vo[i] = std::complex<_Real>(fft_ocdata_[i][0], fft_ocdata_[i][1]);
		}
	};

	template <unsigned _N, typename _Real>
	class Dfft<_N, std::complex<_Real>, _Real> : public fftw3_connect<_N, _Real>
	{
		typedef fftw3_connect<_N, _Real> _Base;
		USING_FFTW3_CONNECT_DEFINES(_N, _Real);

	public:

		/**
		* @brief конструирует объект и делает предварительные вычисления для эффективности
		* @param ndx мультииндекс, описывающий размерности входного массива
		*/
		Dfft(const fft_index &ndx) : _Base(ndx)
		{
			fft_icdata_ = (fft_complex *)fft_malloc(sizeof(fft_complex) * fft_csize_);
			fft_ordata_ = (fft_real *)fft_malloc(sizeof(fft_real) * fft_size_);
			fft_plan_ = fft_plan_dft(fft_icdata_, fft_ordata_, FFTW_ESTIMATE);
		}

		/**
		* @brief перекопирует данные и выполняет преобразование
		* @param vo - выходной массив
		* @param vi - входной массив
		*/
		void operator()(_Real *vo, const std::complex<_Real> *vi)
		{ IMPORT_RUN_EXPORT(vi, vo) }

		/**
		* @brief перекопирует данные и выполняет преобразование
		* @param vo - выходной массив
		* @param vi - входной массив
		*/
		void operator()(const std::complex<_Real> *vi, _Real *vo) { operator()(vo, vi); }

	private:

		void import_(const std::complex<_Real> *vi)
		{
			for (unsigned i=0; i<fft_csize_; i++)
			{
				fft_icdata_[i][0] = vi[i].real();
				fft_icdata_[i][1] = vi[i].imag();
			}
		}

		void execute_() { fft_execute_dft(fft_icdata_, fft_ordata_); }

		void export_(_Real *vo)
		{
			for (unsigned i=0; i<fft_size_; i++)
				vo[i] = fft_ordata_[i];
		}
	};

	template <unsigned _N, unsigned _M, typename _Real>
	class Dfft<_N, vdense_<_M, std::complex<_Real> >, vdense_<_M, _Real> >
	: public fftw3_connect<_N, _Real>
	{
		typedef fftw3_connect<_N, _Real> _Base;
		USING_FFTW3_CONNECT_DEFINES(_N, _Real);

	public:

		/**
		* @brief конструирует объект и делает предварительные вычисления для эффективности
		* @param ndx мультииндекс, описывающий размерности входного массива
		*/
		Dfft(const fft_index &ndx) : _Base(ndx)
		{
			fft_icdata_ = (fft_complex *)fft_malloc(sizeof(fft_complex) * _M * fft_csize_);
			fft_ordata_ = (fft_real *)fft_malloc(sizeof(fft_real) * _M * fft_size_);
			fft_plan_ = fft_plan_dft(_M, fft_icdata_, NULL, _M, 1, fft_ordata_, NULL, _M, 1, FFTW_ESTIMATE);
		}

		/**
		* @brief перекопирует данные и выполняет преобразование
		* @param vo - выходной массив
		* @param vi - входной массив
		*/
		void operator()(vdense_<_M, _Real> *vo, const vdense_<_M, std::complex<_Real> > *vi)
		{ IMPORT_RUN_EXPORT(vi, vo) }

		/**
		* @brief перекопирует данные и выполняет преобразование
		* @param vo - выходной массив
		* @param vi - входной массив
		*/
		void operator()(const vdense_<_M, std::complex<_Real> > *vi, vdense_<_M, _Real> *vo)
		{ operator()(vo, vi); }

	protected:

		/**
		* @brief перекопирует данные внутрь для эффективной обработки
		* @param vi - входной массив
		*/
		void import_(const vdense_<2, std::complex<_Real> > *vi)
		{
			for (unsigned i=0, k=0; i<fft_csize_; i++, k+=2)
			{
				fft_icdata_[k    ][0] = vi[i][0].real(); fft_icdata_[k    ][1] = vi[i][0].imag();
				fft_icdata_[k + 1][0] = vi[i][1].real(); fft_icdata_[k + 1][1] = vi[i][1].imag();
			}
		}

		/**
		* @brief перекопирует данные внутрь для эффективной обработки
		* @param vi - входной массив
		*/
		void import_(const vdense_<3, std::complex<_Real> > *vi)
		{
			for (unsigned i=0, k=0; i<fft_csize_; i++, k+=3)
			{
				fft_icdata_[k    ][0] = vi[i][0].real(); fft_icdata_[k    ][1] = vi[i][0].imag();
				fft_icdata_[k + 1][0] = vi[i][1].real(); fft_icdata_[k + 1][1] = vi[i][1].imag();
				fft_icdata_[k + 2][0] = vi[i][2].real(); fft_icdata_[k + 2][1] = vi[i][2].imag();
			}
		}

		/**
		* @brief перекопирует данные внутрь для эффективной обработки
		* @param vi - входной массив
		*/
		void import_(const vdense_<4, std::complex<_Real> > *vi)
		{
			for (unsigned i=0, k=0; i<fft_csize_; i++, k+=4)
			{
				fft_icdata_[k    ][0] = vi[i][0].real(); fft_icdata_[k    ][1] = vi[i][0].imag();
				fft_icdata_[k + 1][0] = vi[i][1].real(); fft_icdata_[k + 1][1] = vi[i][1].imag();
				fft_icdata_[k + 2][0] = vi[i][2].real(); fft_icdata_[k + 2][1] = vi[i][2].imag();
				fft_icdata_[k + 3][0] = vi[i][3].real(); fft_icdata_[k + 3][1] = vi[i][3].imag();
			}
		}

		/**
		* @brief перекопирует данные внутрь для эффективной обработки
		* @param vi - входной массив
		*/
		template <unsigned M>
		void import_(const vdense_<M, std::complex<_Real> > *vi)
		{
			for (unsigned i=0; i<M * fft_csize_; i+=M)
			for (unsigned k=0; k<M; k++)
			{
				fft_icdata_[i + k][0] = vi[i][k].real();
				fft_icdata_[i + k][1] = vi[i][k].imag();
			}
		}

		/**
		* @brief перекопирует результирующие данные наружу
		* @param vo - выходной массив
		*/
		void export_(vdense_<2, _Real> *vo)
		{
			for (unsigned i=0, k=0; i<fft_size_; i++, k+=2)
			{
				vo[i][0] = fft_ordata_[k    ];
				vo[i][1] = fft_ordata_[k + 1];
			}
		}

		/**
		* @brief перекопирует результирующие данные наружу
		* @param vo - выходной массив
		*/
		void export_(vdense_<3, _Real> *vo)
		{
			for (unsigned i=0, k=0; i<fft_size_; i++, k+=3)
			{
				vo[i][0] = fft_ordata_[k    ];
				vo[i][1] = fft_ordata_[k + 1];
				vo[i][2] = fft_ordata_[k + 2];
			}
		}

		/**
		* @brief перекопирует результирующие данные наружу
		* @param vo - выходной массив
		*/
		void export_(vdense_<4, _Real> *vo)
		{
			for (unsigned i=0, k=0; i<fft_size_; i++, k+=4)
			{
				vo[i][0] = fft_ordata_[k    ];
				vo[i][1] = fft_ordata_[k + 1];
				vo[i][2] = fft_ordata_[k + 2];
				vo[i][3] = fft_ordata_[k + 3];
			}
		}

		/**
		* @brief перекопирует результирующие данные наружу
		* @param vo - выходной массив
		*/
		template <unsigned M>
		void export_(vdense_<M, _Real> *vo)
		{
			for (unsigned i=0, sz=M * fft_size_; i<fft_size_; i+=M)
			for (unsigned k=0; k<M; k++)
				vo[i][k] = fft_ordata_[i + k];
		}

		void execute_() { fft_execute_dft(fft_icdata_, fft_ordata_); }

	};

#endif // USE_FFTW_LIBRARY

}
#endif
