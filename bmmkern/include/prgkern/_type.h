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
#ifndef _TYPE__0077A726_AC88_55c0_4A44_CE4344760600__H
#define _TYPE__0077A726_AC88_55c0_4A44_CE4344760600__H
#include "prgkern/_prgconfig.h"

namespace prgkern
{

	/// type used to overloading of functions and templates
	template <int n> struct Int2Type_ { enum { result = n }; };
	#define _I2T  Int2Type_

	/// type used to overloading of functions and templates
	template <bool n> struct Bool2Type_ { enum { result = n }; };
	#define _B2T  Bool2Type_

	/**
	*  Расширенные типы, это те типы, которые описываются большим числом байт.
	*  Они необходимы в случае накопления сумм для избежания потери точности.
	*  Так, например, любое скалярное произведение векторов большой размерности
	*  требует использования внутреннего сумматора расширенной точности.
	*  Поскольку расширенный тип всегда имеет всего вдвое больше число байт, чем
	*  его базовый тип, то без потери точности можно суммировать не более чем
	*  то количество объектов, которое может описать базовый тип.
	*
	*  Например, если базовый тип есть int, то в накопителе можно суммировать
	*  не более чем ~65000 чисел. При необходимости большего суммирования нужно
	*  еще раз расширить тип, типа _E(_E(type)). Нельзя полностью полагаться
	*  на расширение типа. Автоматическое расширение типа возможно не во всех
	*  случаях, например, для функций scalar_product(n, T1, T2) нельзя на этапе
	*  компиляции подобрать такое расширение, чтобы оно работало для всех n.
	*  Пользователь сам должен следить достаточно ли стандартное расширение для
	*  его целей.
	*/
	template <typename T> struct extended;

	template <> struct extended<float> { typedef float type; };
	template <> struct extended<double> { typedef double type; };
/*	template <> struct extended<float> { typedef double type; };
	template <> struct extended<double> { typedef double type; };*/
		// это замыкание от дальнейшего расширения

	template <> struct extended<char> { typedef short type; };
	template <> struct extended<short> { typedef int type; };
	template <> struct extended<int> { typedef int type; };
		// это замыкание от дальнейшего расширения

	template <> struct extended<unsigned char> { typedef unsigned short type; };
	template <> struct extended<unsigned short> { typedef unsigned int type; };
	template <> struct extended<unsigned int> { typedef unsigned int type; };
		// это замыкание от дальнейшего расширения

	#define _E(type_) extended<type_>::type

	template <typename T, typename S> struct maxtype;

	template <> struct maxtype<float, float> { typedef float type; };
	template <> struct maxtype<float, double> { typedef double type; };
	template <> struct maxtype<double, float> { typedef double type; };
	template <> struct maxtype<double, double> { typedef double type; };


	/**
	*  Предикат для фильтрации любых контейнеров по значению некоторого поля
	*  элемента контейнера. Предикат выделен в отдельный класс, чтобы его не
	*  писать каждый раз и не помещать где-то в глобальном пространстве имен,
	*  не зная куда помещать и не зная где искать. Предикат позволяет фунциям
	*  поиска настроиться на заданное поле внутри класса и сравнивать только
	*  с ним.
	* @param _T класс элемента любого стандартного контейнера
	* @param _F тип поля внутри элемента
	* @param _mp относительный адрес этого поля
	*/
	template <typename _T, typename _F, _F _T::*_mp>
	class equal_
	{
		_F f_;

	public:

		/**
		*  Инициализация предиката.
		* @param f значение, по которому делается поиск элементов контейнера
		*/
		equal_(_F f) : f_(f) {}

		/**
		*  Функция сравнения предиката с элементом.
		* @param f значение, по которому делается поиск элементов контейнера
		*/
		bool operator()(const _T &t) const { return t.*_mp == f_; }
	};

	/**
	*  Предикат для упорядочения векторов индексов по возрастанию элементов
	*  массивов, на которые ссылаются эти индексы. Вот пример.
	*  Допустим, мы имеем массив вещественных чисел f[N], элементы которого
	*  переупорядочивать запрещено или бессмысленно. Например, необходимы
	*  различные упорядочивания этого массива. Мы создаем массив индексов I[N]
	*  (0, 1, 2 ..). Элементы этого массива индесов мы может переупорядочить
	*  таким образом, чтобы первый элемент этого массива I[0] давал номер самого
	*  минимального элемента массива f[], следующий элемент массива I[1] ссылался
	*  на следующий по величине элемент массива f[] и т.д. Как упорядочить
	*  индексы, чтобы они указывали на упорядоченную последовательность в f[]?
	*  Это и делает указанный предикат. Он используется так.
	*    std::sort(&f[0], &f[N], compare_(&I[0]));
	*
	* @param _T класс элемента любого стандартного контейнера
	* @param _F тип поля внутри элемента
	* @param _mp относительный адрес этого поля
	*/
	template <typename _T>
	class less_than_
	{
		const _T *f_;

	public:
		less_than_(const _T *f) : f_(f) {}

		template <typename _I>
		bool operator() (_I i, _I j) const { return f_[i] < f_[j]; }
	};

	/**
	*  Предикат для фильтрации любых контейнеров по значению некоторого поля
	*  элемента контейнера. Предикат выделен в отдельный класс, чтобы его не
	*  писать каждый раз и не помещать где-то в глобальном пространстве имен,
	*  не зная куда помещать и не зная где искать. Предикат позволяет фунциям
	*  поиска настроиться на заданное поле внутри класса и сравнивать только
	*  с ним.
	* @param _T класс элемента любого стандартного контейнера
	* @param _F тип поля внутри элемента
	* @param _mp относительный адрес этого поля
	*/
	template <typename _T, typename _F, _F _T::*_mp>
	class less_
	{
	public:

		/**
		*  Функция сравнения предиката с элементом.
		* @param f значение, по которому делается поиск элементов контейнера
		*/
		bool operator()(const _T &s1, const _T &s2) const
		{ return s1.*_mp < s2.*_mp; }
	};

}
#endif
