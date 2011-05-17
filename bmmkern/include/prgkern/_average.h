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
#ifndef _AVERAGE__0077A726_F90F_5e24_13E9_DD43611C0B00__H
#define _AVERAGE__0077A726_F90F_5e24_13E9_DD43611C0B00__H
#include <vector>

#include "prgkern/_prgconfig.h"
#include "prgkern/_assert.h"
#include "prgkern/_type.h"

namespace prgkern
{

	/**
	* Объект для ведения статистики. Он позволяет вводить любое число данных,
	* однако после некоторого порога, данные начинают замещаться таким образом,
	* что выбрасываются наиболее ранние введенные. Объект позволяет усреднять
	* данные по некоторому интервалу из N значений.
	* @param S структура, поля которой накапливаются и вычисляется среднее
	*/
	template <typename S>
	class Average_ : protected std::vector<S>
	{
		typedef std::vector<S> _Base; // массив значений величины
		int pos_; // позиция последней вставки
		unsigned count_; // количество точек вставленных
		typename extended<S>::type average_;

	public:

		/**
		* @param sz число точек по которым происходит усреднение
		*/
		Average_(unsigned sz=0) : pos_(-1), count_(0), average_(0.) { _Base::resize(sz); }

		/**
		*  Изменение размера накопителя
		* @param sz число точек по которым происходит усреднение
		*/
		void resize(unsigned sz)
		{
			pos_ = -1;
			count_ = 0;
			average_ = 0.;
			_Base::resize(sz);
		}

		/**
		* @return последнее введенное значение
		*/
		const S &top() const
		{
			assert(_GT(count_, (unsigned)0));
			assert(_GE(pos_, (unsigned)0));
			return _Base::operator[](pos_);
		}

		/// возврат уже насчитанного среднего значения
		typename extended<S>::type average() const { return average_; }

		/**
		* Запомнить еще одно значение.
		*/
		void push(const S &s)
		{
			pos_++; pos_ %= _Base::size();
			if (count_ < _Base::size())
			{
				average_ = (average_ * count_ + s) / (count_ + 1);
				count_++;
			}
			else
			{
				average_ += (s - _Base::operator[](pos_)) / count_;
			}
			_Base::operator[](pos_) = s;

//			pos_++; pos_ %= _Base::size();
//			_Base::operator[](pos_) = s;
//			if (count_ < _Base::size()) count_++;
		}

		/**
		* @return среднее значение величины
		*/
		S operator()() const
		{
			assert(_GT(count_, (unsigned)0));
			typename extended<S>::type s; s = 0;
			for (unsigned i=0; i<count_; i++) s += _Base::operator[](i);
			s *= ((S)1. / count_);
			return (S) s;
		}

		/**
		* @return среднее значение величины по последним n точкам
		*/
		S operator()(unsigned n) const
		{
			assert(_LT(n, count_));

			int pos = (int )pos_;
			typename extended<S>::type s = 0;
			for (unsigned k=0; k<n; k++)
			{
				s += _Base::operator[](pos--);
				if (pos == -1) pos = count_ - 1;
			}
			s *= ((S)1. / n);
			return (S) s;
		}
	};

}
#endif
