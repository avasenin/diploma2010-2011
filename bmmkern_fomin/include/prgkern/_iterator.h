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
#ifndef _ITERATOR__F9ED1116_5C6E_5fd6_C36C_114687FB0400__H
#define _ITERATOR__F9ED1116_5C6E_5fd6_C36C_114687FB0400__H
#include "prgkern/_prgconfig.h"

namespace prgkern
{

	/**
	* @brief iterator producing sequense such as [-1, 0, 1, 2, 3, ...]
	*/
	template <typename T>
	class range_iterator_
	{
		T pos_;
	public:

    range_iterator_(T pos=0) : pos_(pos) {}
		range_iterator_(const range_iterator_ &t) : pos_(t.pos_) {}

		range_iterator_ &operator=(const range_iterator_ &t) { pos_ = t.pos_; return *this; }

		T operator*() const { return pos_; }
		T &operator*() { return pos_; }
		range_iterator_ &operator++() { ++pos_; return *this; }
		range_iterator_ &operator+=(int n) { pos_ += n; return *this; }

		bool operator!=(const range_iterator_ &t) const { return pos_ != t.pos_; }
	};
	typedef range_iterator_<int>  range_iterator;
	typedef range_iterator_<unsigned>  urange_iterator;

	/**
	*  Итератор, проходящий любое подмножество заданного std::vector<>
	*  в произвольном порядке. Порядок обхода задается внешним итератором
	*  Iterator, который возвращает индексы данного массива.
	*  Если внешний итератор возвращает индексы элементов, отсутствующие в массиве,
	*  результат неопределен.
	* @param S тип элемента заданного std::vector<>
	* @param Iterator внешний итератор, дающий индексы массива
	*/
	template <typename S, typename Iterator>
	class array_iterator
	{
		Iterator it_; // индексный итератор
		S *v_; // опорный массив

	public:

		typedef S value_type; // нужен для внешних пользователей

		/**
		*  Конструктор итератора. Стандартный вариант конструирования следующий:
		*    array_iterator<V, A::iterator> it(v, i.begin());
		* @param v опорный массив, чьи элементы будут возвращаться
		* @param it итератор, генерирующий индексы опорного массива
		*/
		array_iterator() : it_(), v_(0) {}
		array_iterator(S *v, const Iterator &it) : it_(it), v_(v) {}
		array_iterator(const array_iterator &it) : it_(it.it_), v_(it.v_) {}

		array_iterator &operator=(const array_iterator &it)
		{ it_ = it.it_; v_ = it.v_; return *this; }

		array_iterator &operator++() { ++it_; return *this; }

		S &operator*() { return v_[*it_]; }
		const S &operator*() const { return v_[*it_]; }

		S *operator->() { return &v_[*it_]; }
		const S *operator->() const { return &v_[*it_]; }

		bool operator!=(const array_iterator &it) const { return it_ != it.it_; }
	};

	/**
	*  Итератор, проходящий любое подмножество заданного std::vector<>
	*  в произвольном порядке. Порядок обхода задается внешним итератором
	*  Iterator, который возвращает индексы данного массива.
	*  Если внешний итератор возвращает индексы элементов, отсутствующие в массиве,
	*  результат неопределен.
	* @param S тип элемента заданного std::vector<>
	* @param Iterator внешний итератор, дающий индексы массива
	*/
	template <typename S, typename Iterator>
	class const_array_iterator
	{
		Iterator it_; // индексный итератор
		const S *v_; // опорный массив

	public:

		typedef S value_type; // нужен для внешних пользователей

		/**
		*  Конструктор итератора. Стандартный вариант конструирования следующий:
		*    array_iterator<V, A::iterator> it(v, i.begin());
		* @param v опорный массив, чьи элементы будут возвращаться
		* @param it итератор, генерирующий индексы опорного массива
		*/
		const_array_iterator() : it_(), v_(0) {}
		const_array_iterator(const S *v, const Iterator &it) : it_(it), v_(v) {}
		const_array_iterator(const const_array_iterator &it) : it_(it.it_), v_(it.v_) {}

		const_array_iterator &operator=(const const_array_iterator &it)
		{ it_ = it.it_; v_ = it.v_; return *this; }

		const_array_iterator &operator++() { ++it_; return *this; }
		const S &operator*() const { return v_[*it_]; }
		const S *operator->() const { return &v_[*it_]; }

		bool operator!=(const const_array_iterator &it) const { return it_ != it.it_; }
	};

	/**
	*  Итератор, проходящий любое подмножество заданного std::vector<>
	*  в произвольном порядке. Порядок обхода задается внешним итератором
	*  Iterator, который возвращает индексы данного массива.
	*  Если внешний итератор возвращает индексы элементов, отсутствующие в массиве,
	*  результат неопределен.
	* @param S тип элемента заданного std::vector<>
	* @param Iterator внешний итератор, дающий индексы массива
	*/
	template <typename Iterator, typename Member, Member Iterator::value_type :: *MemberPtr>
	class sub_iterator
	{
		Iterator it_; // индексный итератор

	public:

		typedef Member value_type; // нужен для внешних пользователей

		/**
		*  Конструктор итератора. Стандартный вариант конструирования следующий:
		*    array_iterator<V, A::iterator> it(v, i.begin());
		* @param v опорный массив, чьи элементы будут возвращаться
		* @param it итератор, генерирующий индексы опорного массива
		*/
		sub_iterator() : it_() {}
		sub_iterator(const Iterator &it) : it_(it) {}
		sub_iterator(const sub_iterator &it) : it_(it.it_) {}

		sub_iterator &operator=(const sub_iterator &it) { it_ = it.it_; return *this; }
		sub_iterator &operator++() { ++it_; return *this; }

		Member &operator*() { return it_->*MemberPtr; }
		const Member &operator*() const { return it_->*MemberPtr; }

		Member *operator->() { return &it_->*MemberPtr; }
		const Member *operator->() const { return &it_->*MemberPtr; }

		bool operator!=(const sub_iterator &it) const { return it_ != it.it_; }
	};

	/**
	 * Структура данных, хранящая несколько массивов векторов, и позволяющая сделать полный
	 * перебор элементов во всех хранимых векторах. Можно сказать, что это обобщение вектора,
	 * но элементы которого не помещены в одну область памяти, а расбросаны.
	 *
	 * Структура может быть объявлена на константные объекты, например piecewise_vector<const T>,
	 * или на неконстантные, как piecewise_vector<T>. То есть, именно она, а не тип итератора
	 * определяет тип доступа к элементам. Итератор ничего не знает об этой структуре, его задача
	 * сделать обход элементов, попытаться дать доступ, и если доступа нет, то это не его забота.
	 */
	template <typename T>
	class piecewise_vector
	{
		typedef std::pair<T *, unsigned>   _piece;

		std::vector<_piece> all_pieces_;

		class iterator_
		{
			_piece *pieces_;
			unsigned count_;
			_piece *cur_piece_;
			T *elem_;

		public:

			iterator_() : pieces_(0), count_(0), cur_piece_(0), elem_(0) {}

			iterator_(_piece *p, unsigned count) : pieces_(p), count_(count), cur_piece_(p), elem_(0)
			{ if (p && count) elem_ = cur_piece_->first; }

			iterator_(const iterator_ &it) : pieces_(it.pieces_), count_(it.count_),
				cur_piece_(it.cur_piece_), elem_(it.elem_) {}

			T &operator*() const { return *elem_; }

			iterator_ &operator++()
			{
				if (++elem_ - cur_piece_->first == cur_piece_->second)
				{
					if ((unsigned)(++cur_piece_ - pieces_) == count_) elem_ = 0;
					else elem_ = cur_piece_->first;
				}
				return *this;
			}

			iterator_ &operator+=(unsigned n)
			{
				elem_ += n;
				if (elem_ - cur_piece_->first >= cur_piece_->second)
				{
					n = elem_ - (cur_piece_->first + cur_piece_->second);
					++cur_piece_;

					while ( (n >= cur_piece_->second)
						&& ((unsigned)(cur_piece_ - pieces_) < count_)
					)
					{
						n -= cur_piece_->second;
						++cur_piece_;
					}

					if ((unsigned)(cur_piece_ - pieces_) == count_) elem_ = 0;
					else elem_ = cur_piece_->first + n;
				}
				return *this;
			}
			iterator_ operator+(unsigned n) { return iterator_(*this) += n; }

			bool operator!=(const iterator_ &it) const { return elem_ != it.elem_; }
		};

	public:

		typedef iterator_  iterator;

		piecewise_vector() {}
		piecewise_vector(T *arr, unsigned count) { insert(arr, count); }

		void insert(T *arr, unsigned count)
		{ if (arr && count) all_pieces_.push_back(_piece(arr, count)); }

		iterator begin()
		{
			unsigned count = all_pieces_.size();
			if (count) return iterator(&all_pieces_[0], count);
			return iterator();
		}
		iterator end() { return iterator(); }

		unsigned count() const
		{
			unsigned cnt = 0;
			for (unsigned i=0,sz=all_pieces_.size(); i<sz; i++) cnt += all_pieces_[i].second;
			return cnt;
		}

	};

}
#endif
