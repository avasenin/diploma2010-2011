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
#ifndef _GRAPH___0077A726_E6E3_58a1_C16D_CE436AB90501__H
#define _GRAPH___0077A726_E6E3_58a1_C16D_CE436AB90501__H

#include <boost/tuple/tuple.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

namespace prgkern
{
	/**
	 * Функция позволяет найти все подграфы в молекуле. Для этого молекула представляется как граф
	 * из заданного числа точек, а валентные связи молекулы как ребра этого графа. Ротамерные связи
	 * представляют собой ребра, которые разделяют группы атомов на подграфы. Каждый такой подграф
	 * и является полным ротамером (совокупностью атомов, вращающихся относительно ротамерных связей
	 * как единое целое).
	 * @param n - число атомов в молекуле
	 * @param dest - начало массива "цветов" атомов, по которым разделяются подграфы
	 * @param first - итератор по связям молекулы (два номера атома, представляющим связь)
	 * @param last - запредельный итератор по связям
	 * @return количество подграфов
	 */
	template <typename _BondIterator>
	unsigned make_subgraphs(unsigned n, int *dest, _BondIterator first, _BondIterator last)
	{
		typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;

		Graph G(first, last, n); // создали граф со всеми ребрами
		return boost::connected_components(G, dest); // разделили граф на подграфы
	}

}
#endif
