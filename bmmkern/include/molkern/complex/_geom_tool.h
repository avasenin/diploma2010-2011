#ifndef _GEOM_TOOL__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H
#define _GEOM_TOOL__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H

#include "molkern/__moldefs.h"

namespace molkern
{
	using namespace prgkern;

	/**
	*   Вращает набор атомов молекулы в 3D пространстве.
	* @param X0 центр вращения
	* @param angle углы Эйлера для вращения
	* @param it,ite итераторы начала и конца последовательности индексов атомов
	*/
	template <typename Atom, typename Iterator>
	INLINE void rotate(_I2T<EULER_ROTATOR_>, Atom *atom, const vector_t &angle, const vector_t &X0,
		Iterator it, Iterator ite)
	{
		Rotator<EULER_ROTATOR_, real_t> rotator(angle, X0);
		for (; it!=ite; ++it) atom[*it].X = rotator(atom[*it].X);
	}

	template <typename Atom, typename Iterator>
	INLINE void rotate(_I2T<XYZ_ROTATOR_>, Atom *atom, const vector_t &angle, const vector_t &X0,
		Iterator it, Iterator ite)
	{
		Rotator<XYZ_ROTATOR_, real_t> rotator(angle, X0);
		for (; it!=ite; ++it) atom[*it].X = rotator(atom[*it].X);
	}

	/**
	*   Вращает набор атомов молекулы в 3D пространстве.
	* @param X0 центр вращения
	* @param angle углы Эйлера для вращения
	* @param it,ite итераторы начала и конца последовательности индексов атомов
	*/
	template <typename Atom, typename __Atom, typename Iterator>
	INLINE void rotate(_I2T<EULER_ROTATOR_>, Atom *atom, const __Atom *atom__,
		const vector_t &angle, const vector_t &X0, Iterator it, Iterator ite)
	{
		Rotator<EULER_ROTATOR_, real_t> rotator(angle, X0);
		for (; it!=ite; ++it) atom[*it].X = rotator(atom__[*it].X);
	}

	/**
	*   Сдвигает набор атомов молекулы в пространстве любой размерности.
	* @param X вектор сдвига
	* @param it,ite итераторы начала и конца последовательности индексов атомов
	*/
	template <typename Atom, typename __Atom, typename Iterator>
	INLINE void move(Atom *atom, const __Atom *atom__, const vector_t &X,
		Iterator it, Iterator ite)
	{
		for (; it!=ite; ++it) atom[*it].X = atom__[*it].X + X;
	}

	/**
	*   Сдвигает набор атомов молекулы в пространстве любой размерности.
	* @param molecule молекула, которая сдвигается
	* @param X вектор сдвига
	* @param it,ite итераторы начала и конца последовательности индексов атомов
	*/
	template <typename Atom, typename Iterator>
	INLINE void move(Atom *atom, const vector_t &X, Iterator it, Iterator ite)
	{
		for (; it!=ite; ++it) atom[*it].X += X;
	}

	/**
	*   Сдвигает набор атомов молекулы таким образом, чтобы их геометрический
	*   центр попал в заданную точку.
	* @param X точка, куда должен быть сдвинут центр группы атомов
	* @param it,ite итераторы начала и конца последовательности индексов атомов
	*/
	template <typename Atom, typename Iterator>
	INLINE void move_to(Atom *atom, const vector_t &X, Iterator it, Iterator ite)
	{
		vector_t cm = calculate(GEOM_CENTER, atom, it, ite);
		move(atom, X - cm, it, ite);
	}

	/**
	*   Рассчиывает минимальный ящик, который охватывает все атомы в заданной
	*   последовательности. Последовательность дает индексы атомов. Она может
	*   быть задана любым типом (вектор, список и т.д.).
	* @param it,ite итераторы начала и конца последовательности индексов атомов
	* @param molecule молекула (комплекс), чьи атомы используются
	* @return ящик
	*/
	template <typename Atom, typename Iterator>
	INLINE Box_<3, real_t> calculate(_I2T<BOX_>, const Atom *atom, Iterator it, Iterator ite)
	{
		Box_<3, real_t> box(vector_t(INFINITY), vector_t(-INFINITY));
		for (; it!=ite; ++it)
		{
			box.bottom() = min(box.bottom(), atom[*it].X);
			box.top()    = max(box.top(),    atom[*it].X);
		}
		return box;
	}

	/**
	*   Рассчитывает геометрический центр группы атомов молекулы
	* @param molecule молекула (комплекс), чьи атомы используются
	* @param it,ite итераторы начала и конца последовательности индексов атомов
	* @return найденный центр группы
	*/
	template <typename Atom, typename Iterator>
	INLINE vector_t	calculate(_I2T<GEOM_CENTER_>, const Atom *atom, Iterator it, Iterator ite)
	{
		vector_t cm = 0; unsigned count = 0;
		for (; it!=ite; ++it, ++count) cm += atom[*it].X;
		cm *= (1./ count);
		return cm;
	}

	template <typename Atom, typename Iterator>
	INLINE vector_t	calculate(_I2T<MASS_CENTER_>, const Atom *atom, Iterator it, Iterator ite)
	{
		vector_t cm = 0.; real_t mass = 0., m;
		for (; it!=ite; ++it) { m = atom[*it]->mass; mass += m; cm += (atom[*it].X) * m; }
		cm *= (1./ mass);
		return cm;
	}

	template <typename Atom, typename Iterator>
	INLINE vector_t	calculate_center_of_mass(const Atom *atom, Iterator it, Iterator ite)
	{
		vector_t cm = 0.; real_t mass = 0., m;
		for (; it!=ite; ++it) { m = atom[*it]->mass; mass += m; cm += (atom[*it].X) * m; }
		cm *= (1./ mass);
		return cm;
	}

	/**
	*   Рассчитывает наибольший Ван-дер-Ваальсовый радиус в группе атомов.
	* @param molecule молекула (комплекс)
	* @param it,ite итераторы начала и конца последовательности индексов атомов
	* @return Ван-дер-Ваальсовый радиус
	*/
	template <typename Atom, typename Iterator>
	INLINE real_t calculate(_I2T<MAX_VDW_RADIUS_>, const Atom *atom,
		Iterator it, Iterator ite)
	{
		real_t max_radius = 0.;
		for (; it!=ite; ++it) if (max_radius < atom[*it]->radius)
			max_radius = atom[*it]->radius;
		return max_radius;
	}

}
#endif
