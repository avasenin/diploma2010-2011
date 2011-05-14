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
#ifndef _TIME__0077A726_0EE7_5e96_7E5A_CE432FBA0800__H
#define _TIME__0077A726_0EE7_5e96_7E5A_CE432FBA0800__H

#include <boost/preprocessor/iteration/local.hpp>
#include "prgkern/_prgconfig.h"
#include "prgkern/_string.h"
#include <ctime>

namespace prgkern
{
	struct system_time_t
	{
		double sec_;
		long nsec_;
		system_time_t(double t=0., long ns=0.) : sec_(t), nsec_(ns) {}
		operator double() const { return sec_ + nsec_ * 1.e-9; }
		operator float() const { return float(sec_ + nsec_ * 1.e-9); }
	};

	INLINE system_time_t operator-(system_time_t t1, system_time_t t2)
	{
		return system_time_t(t1.sec_ - t2.sec_, t1.nsec_ - t2.nsec_);
	}

	INLINE std::string make_string(system_time_t tm)
	{ return make_string("%14.6f", (float)tm) + _S(" s"); }

#ifdef USE_LINUX_SPECIFIC_CODE

	//---------------------------------------------------------------------------------------
	// Необходимо подкулючать библиотеку -lrt для того, чтобы использовать часы стандарта
	// POSIX-2001 с точностью до 1 нс.
	//---------------------------------------------------------------------------------------
	#include <sys/time.h>

	// Такие сложные имена идентификаторов выбраны для того, чтобы избежать
	// ситуации, когда эти имена совпадают со внешними. В этом случае, если внешние
	// имена не вызываются внутри TIME_TESTING_START {...} TIME_TESTING_FINISH, то
	// проблемы не возникает, но если вызываются, то результат непредсказуем.

	#define TIME_TESTING_START(msg, iteration) { \
		std::string msg_09267106_A953_52b7_A54A_024279B90300(msg); \
		struct timeval tv_09267106_A953_52b7_A54A_024279B90300; \
		::gettimeofday(&tv_09267106_A953_52b7_A54A_024279B90300, NULL); \
		int iteration_09267106_6E53_5abf_C84A_02425EE50600 = iteration; \
		int time1_09267106_6E53_5abf_C84A_02425EE50600 \
			= tv_09267106_A953_52b7_A54A_024279B90300.tv_sec; \
		long time1ms_09267106_6E53_5abf_C84A_02425EE50600 \
			= tv_09267106_A953_52b7_A54A_024279B90300.tv_usec; \
		for (int i__09267106_6E53_5abf_C84A_02425EE50600=0; \
			i__09267106_6E53_5abf_C84A_02425EE50600<(iteration); \
			i__09267106_6E53_5abf_C84A_02425EE50600++) {

	#define TIME_TESTING_FINISH } \
			::gettimeofday(&tv_09267106_A953_52b7_A54A_024279B90300, NULL); \
			int time2_09267106_6E53_5abf_C84A_02425EE50600 \
				= tv_09267106_A953_52b7_A54A_024279B90300.tv_sec; \
			long time2ms_09267106_6E53_5abf_C84A_02425EE50600 \
				= tv_09267106_A953_52b7_A54A_024279B90300.tv_usec; \
			double diff = iteration_09267106_6E53_5abf_C84A_02425EE50600 > 0 ? \
				double(time2_09267106_6E53_5abf_C84A_02425EE50600 \
				- time1_09267106_6E53_5abf_C84A_02425EE50600 \
				+ (time2ms_09267106_6E53_5abf_C84A_02425EE50600 \
				- time1ms_09267106_6E53_5abf_C84A_02425EE50600) * 0.000001 ) \
				/ iteration_09267106_6E53_5abf_C84A_02425EE50600 : 0; \
			if (!msg_09267106_A953_52b7_A54A_024279B90300.empty()) \
				if (msg_09267106_A953_52b7_A54A_024279B90300.length()) \
					msg_09267106_A953_52b7_A54A_024279B90300 += _S("\n"); \
			msg_09267106_A953_52b7_A54A_024279B90300 += _S("@   time difference is "); \
			if (diff < 60) msg_09267106_A953_52b7_A54A_024279B90300 += prgkern::ftoa(diff) + _S(" s"); \
			else \
			{ \
				int minutes = (int)floor(diff / 60); diff -= minutes * 60; \
				msg_09267106_A953_52b7_A54A_024279B90300 += prgkern::itoa(minutes) \
					+ _S(" m ") + prgkern::ftoa(diff) + _S(" s"); \
			} \
			PRINT_MESSAGE(msg_09267106_A953_52b7_A54A_024279B90300); \
		}
#ifdef USE_POSIX_2001
	INLINE system_time_t current_time()
	{
		timespec tm;
		clock_gettime(CLOCK_REALTIME, &tm);
		return system_time_t(tm.tv_sec, tm.tv_nsec);
	}
#else
	INLINE system_time_t current_time()
	{
		timeval tm;
		gettimeofday(&tm, NULL);
		return system_time_t(tm.tv_sec, tm.tv_usec * 1000);
	}
#endif

#else

	#define TIME_TESTING_START(iteration)
	#define TIME_TESTING_FINISH

	INLINE system_time_t current_time() { time_t tm; time(&tm); return system_time_t(tm, 0); }

#endif

	template <int N, unsigned NT=8> class function_timer_t
	{
		system_time_t start_time_;
		unsigned n_;

	public:

		static double accumulate_time[NT + 1];
		static long unsigned accumulate_count[NT + 1];

		function_timer_t(unsigned n=0) : start_time_(current_time()), n_(n) { accumulate_count[n_]++; }
		~function_timer_t() { accumulate_time[n_] += double(current_time() - start_time_); }

		static void clear()
		{
			for (unsigned i=0; i<NT+1; i++) { accumulate_time[i] = 0; accumulate_count[i] = 0; }
		}

		static unsigned slow_thread()
		{
			double max_time = accumulate_time[0]; unsigned slow_ndx = 0;
			for (unsigned i=1; i<NT; i++)
			{
				if (max_time < accumulate_time[i])
				{
					max_time = accumulate_time[i];
					slow_ndx = i;
				}
			}
			return slow_ndx;
		}

		static void sinc()
		{
			unsigned slow_ndx = slow_thread();
			accumulate_time[NT] += accumulate_time[slow_ndx];
			accumulate_count[NT] += accumulate_count[slow_ndx];
			for (unsigned i=0; i<NT; i++) { accumulate_time[i] = 0; accumulate_count[i] = 0; }
		}

		static unsigned count() { return accumulate_count[NT]; }
		static double time() { return accumulate_time[NT]; }

	};

	template <int N, unsigned NT> double function_timer_t<N, NT>::accumulate_time[NT + 1];
	template <int N, unsigned NT> long unsigned function_timer_t<N, NT>::accumulate_count[NT + 1];
}

#endif
