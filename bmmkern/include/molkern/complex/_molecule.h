#ifndef _MOLECULE__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H
#define _MOLECULE__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H

#include "boost/tuple/tuple.hpp"

#include "molkern/__moldefs.h"
#include "molkern/complex/_geom_tool.h"
#include "molkern/forcefield/_rotamer.h"

namespace molkern
{
	using namespace prgkern;

	/**
	*  Класс Molecule_ является промежутоным типом между _Archetype и Complex_.
	*  Он является служебным и только организует данных по атомам, принадлежащим
	*  одной молекуле вместе. Он выделяет из архетипа переменную часть,
	*  связанную с изменением числа свободных атомов, то есть позволяя иметь
	*  ряд совершенно одинаковых молекул, но с разным внутренним состоянием.
	*
	*  Класс не должен зависить от силового поля и должен обрабатывать любой
	*  тип атомов, который ему предложит комплекс. По этой причине, для этого
	*  класса меняются шаблонные параметры.
	*/
	template <typename _Archetype>
	class Molecule_
	{
		typedef typename _Archetype::atomdata_type  _Atomdata;
		typedef typename _Archetype::rotamer_type   _Rotamer;
		typedef typename _Archetype::edge_type      _Edge;
		typedef typename _Archetype::chain_type     _Chain;
		typedef typename _Rotamer::stick_type       _Stick;

		const _Archetype *archetype_; // прототип молекулы
		unsigned freedom_type_; // тип свободы молекулы
			// (FREEDOM_FIXED_, FREEDOM_CM_, FREEDOM_CHAIN_, FREEDOM_ROTAMER_, FREEDOM_ATOM_)
		unsigned calc_type_; // какие энергетические вклады учитываются

		unsigned_t *free_atoms_;
		unsigned_t *free_bonds_;
		unsigned_t *free_angles_;
		unsigned_t *free_torsions_;
		unsigned_t *free_pair14s_;

		unsigned_t free_atoms_count_;
		unsigned_t free_bonds_count_;
		unsigned_t free_angles_count_;
		unsigned_t free_torsions_count_;
		unsigned_t free_pair14s_count_;

	public:

		enum { dimension = 3 };
		typedef real_t  real_type;

		Molecule_(const _Archetype *archetype, unsigned freedom_type)
		: archetype_(archetype), freedom_type_(freedom_type), calc_type_(CALC_NOTHING_)
		{
			// Работаю с памятью, как с "сырой" памятью, не собираясь ее перемещать, инициализировать
			// при загрузке и т.д. По этой причине использую функции, типа calloc().
			free_atoms_    = (unsigned_t *)calloc(archetype_->count(_I2T<ATOM_>()),    sizeof(unsigned_t));
			free_bonds_    = (unsigned_t *)calloc(archetype_->count(_I2T<BOND_>()),    sizeof(unsigned_t));
			free_angles_   = (unsigned_t *)calloc(archetype_->count(_I2T<ANGLE_>()),   sizeof(unsigned_t));
			free_torsions_ = (unsigned_t *)calloc(archetype_->count(_I2T<TORSION_>()), sizeof(unsigned_t));
			free_pair14s_  = (unsigned_t *)calloc(archetype_->count(_I2T<PAIR14_>()),  sizeof(unsigned_t));

			free_atoms_count_    = archetype_->count(ATOM);
			free_bonds_count_    = archetype_->count(BOND);
			free_angles_count_   = archetype_->count(ANGLE);
			free_torsions_count_ = archetype_->count(TORSION);
			free_pair14s_count_  = archetype_->count(PAIR14);

			for (unsigned i=0; i<free_atoms_count_; i++) free_atoms_[i] = i;
			for (unsigned i=0; i<free_bonds_count_; i++) free_bonds_[i] = i;
			for (unsigned i=0; i<free_angles_count_; i++) free_angles_[i] = i;
			for (unsigned i=0; i<free_torsions_count_; i++) free_torsions_[i] = i;
			for (unsigned i=0; i<free_pair14s_count_; i++) free_pair14s_[i] = i;

			if (freedom_type_ & YES_ATOM_)
				calc_type_ = CALC_BOND_ | CALC_ANGLE_ | CALC_TORSION_ | CALC_PAIR14_;
			else if (freedom_type_ & YES_ROTAMER_)
				calc_type_ = CALC_TORSION_ | CALC_PAIR14_;
			else calc_type_ = CALC_NOTHING_;
		}

		~Molecule_()
		{
			free((void*)free_atoms_);
			free((void*)free_bonds_);
			free((void*)free_angles_);
			free((void*)free_torsions_);
			free((void*)free_pair14s_);
		}

		/**
		* @brief сохранение молекулы в поток с различными форматами
		* @param file поток вывода
		* @param prn_hydrogens печатать (или нет) атомы водорода
		* @param header вставка в поток комментария
		*/
	#define IMPLEMENT_SAVE_FUNCTION(format) \
		template <typename _Atom> \
		unsigned save(_I2T<format>, std::ofstream &file, const _Atom *atoms, bool prn_hydrogens) const \
		{ return archetype_->save(_I2T<format>(), file, atoms, prn_hydrogens); }

		IMPLEMENT_SAVE_FUNCTION(FORMAT_PDB_)
		IMPLEMENT_SAVE_FUNCTION(FORMAT_HIN_)
		IMPLEMENT_SAVE_FUNCTION(FORMAT_MOL2_)
		IMPLEMENT_SAVE_FUNCTION(FORMAT_BMM_)

	#undef IMPLEMENT_SAVE_FUNCTION


		/**   ФУНКЦИЯ ДЛЯ НАХОЖДЕНИЯ ОБЪЕКТОВ, СВЯЗАННЫХ С ЗАДАННЫМИ АТОМАМИ
		*  Данные функции используют, например, для получения подмножетва связей,
		*  углов и т.д., которые связаны со списком свободных атомов. Имея такое
		*  подмножество, при расчетах энергий можно не вычислять вклады, которые
		*  связаны с фиксированными атомами, и которые дают всего лишь константу
		*  в полную энергию.
		*
		* @param start,end итераторы последовательности индексов свободных атомов
		*/
		template <typename _Iterator>	void make_free(_Iterator it, _Iterator ite);

		/// возвращает число свободных атомов
		unsigned count(_I2T<FREE_>) const
		{
		#ifdef USE_FIXED_ATOMS
			return free_atoms_count_;
		#else
			return archetype_->count(_I2T<ATOM_>());
		#endif
		}

		/// текущее число степеней движения молекулы
		unsigned count(_I2T<FREEDOM_>) const { return archetype_->count(FREEDOM, freedom_type_); }

		template <typename _Atom> unsigned read(_I2T<POSITION_>, const real_t *x, _Atom *atoms) const
		{ return archetype_->read(POSITION, x, freedom_type_, atoms); }

		template <typename _Atom> unsigned write(_I2T<GRADIENT_>, real_t *g, const _Atom *atoms) const
		{ return archetype_->write(GRADIENT, g, freedom_type_, atoms); }

		unsigned get(_I2T<FREEDOM_TYPE_>) const { return freedom_type_; }

		/**
		*    ФУНКЦИИ ДЛЯ РАСЧЕТА ПОТЕНЦИАЛЬНОЙ ЭНЕРГИИ И ЕЕ ПРОИЗВОДНЫХ
		* @param region область пространства для расчета
		* @return полная энергия взаимодействия
		*/
		template <typename _Atom> _E(real_t) U(const _Atom *atoms, bool make_print=YES_PRINT) const;
		template <typename _Atom> _E(real_t) dU__dQ(_Atom *atoms) const;
		template <typename _Atom> _E(real_t) dU__dX(_Atom *atoms) const;

		/**
		*      Функция доступа к полному интерфейсу класса Archetype.
		*  Используется, чтобы не писать кучу перенаправляющих функций, типа
		*     count(_I2T<...>){ return archetype->count(_I2T<...>()); }
		*  Например, чтобы вызвать функцию get(_I2T<ATOM_>) класса Archetype
		*  достаточно написать для объекта(!) молекулы вызов:
		*  molecule->get(_I2T<ATOM_>()), хотя такой функции в классе молекулы нет.
		*  Примечание: molecule является не указателем, а объектом.
		*/
		const _Archetype *operator->() const { return archetype_; }

	};

	#define TEMPLATE_HEADER  template <typename _Archetype>
	#define TEMPLATE_ARG     _Archetype

	TEMPLATE_HEADER
	template <typename _Iterator>
	INLINE void Molecule_<TEMPLATE_ARG>
	::make_free(_Iterator it, _Iterator ite)
	{
		free_atoms_count_ = 0;
		for (_Iterator it__=it; it__!=ite; ++it__) free_atoms_[free_atoms_count_++] = *it__;

		free_bonds_count_    = archetype_->free(_I2T<BOND_>(),    &free_bonds_[0],    it, ite);
		free_angles_count_   = archetype_->free(_I2T<ANGLE_>(),   &free_angles_[0],   it, ite);
		free_torsions_count_ = archetype_->free(_I2T<TORSION_>(), &free_torsions_[0], it, ite);
		free_pair14s_count_  = archetype_->free(_I2T<PAIR14_>(),  &free_pair14s_[0],  it, ite);
	}

	TEMPLATE_HEADER
	template <typename _Atom>
	INLINE _E(real_t) Molecule_<TEMPLATE_ARG>
	::U(const _Atom *atoms, bool make_print) const
	{
		if (make_print)
		{
			std::string msg = _S("molecule ") + archetype_->name() + _S(" :");
			PRINT_MESSAGE(msg);
		}

		_E(real_t) energy = 0.;

		if (calc_type_ & CALC_BOND_) energy += archetype_->U(BOND, atoms, &free_bonds_[0],
				&free_bonds_[0] + free_bonds_count_, make_print);

		if (calc_type_ & CALC_ANGLE_) energy += archetype_->U(ANGLE, atoms, &free_angles_[0],
				&free_angles_[0] + free_angles_count_, make_print);

		if (calc_type_ & CALC_TORSION_) energy += archetype_->U(TORSION, atoms, &free_torsions_[0],
				&free_torsions_[0] + free_torsions_count_, make_print);

		if (calc_type_ & CALC_PAIR14_) energy += archetype_->U(PAIR14, atoms, &free_pair14s_[0],
					&free_pair14s_[0] + free_pair14s_count_, make_print);

		return energy;
	}

	TEMPLATE_HEADER
	template <typename _Atom>
	INLINE _E(real_t) Molecule_<TEMPLATE_ARG>
	::dU__dQ(_Atom *atoms) const
	{
            return archetype_->dU__dQ(atoms);
	}

        TEMPLATE_HEADER
	template <typename _Atom>
	INLINE _E(real_t) Molecule_<TEMPLATE_ARG>
	::dU__dX(_Atom *atoms) const
	{
		_E(real_t) energy = 0.;

		if (calc_type_ & CALC_BOND_) energy += archetype_->dU__dX(BOND, atoms, &free_bonds_[0],
				&free_bonds_[0] + free_bonds_count_);

		if (calc_type_ & CALC_ANGLE_) energy += archetype_->dU__dX(ANGLE, atoms, &free_angles_[0],
				&free_angles_[0] + free_angles_count_);

		if (calc_type_ & CALC_TORSION_) energy += archetype_->dU__dX(TORSION, atoms, &free_torsions_[0],
				&free_torsions_[0] + free_torsions_count_);

		if (calc_type_ & CALC_PAIR14_) energy += archetype_->dU__dX(PAIR14, atoms, &free_pair14s_[0],
					&free_pair14s_[0] + free_pair14s_count_);

		return energy;
	}

	#undef TEMPLATE_HEADER
	#undef TEMPLATE_ARG

}
#endif
