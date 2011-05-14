#ifndef LINKCELL__09ED1116_C3FE_55a4_CA63_F84536C80D02__H
#define LINKCELL__09ED1116_C3FE_55a4_CA63_F84536C80D02__H

#include "molkern/__moldefs.h"
#include "molkern/complex/_parallel.h"
#include "molkern/complex/_region.h"

namespace molkern
{
	using namespace prgkern;

	/**
	 * LC объект используется расчетов ближних невалентных взаимодействий. Вместо него можно
	 * использовать Верлет таблицу, которая более эффективна при низкой плотности атомов.
	 */
	template <typename _LPComplex>
	class _Link_cell : public near_range_integrator<_LPComplex>
	{
		typedef typename _LPComplex::region_type  _Region;
		typedef typename _LPComplex::atom_type    _Atom;
		typedef Parall_<_Link_cell>               _Parall;

		typedef __LJAtom<vreal_t, vint_t>         _LJAtom;
		typedef typename _Region::value_type      _GridAtom;
		typedef typename _Region::key_type        _Key;
		typedef typename _Region::index_type      _Index;
		typedef typename _Region::node_type       _Node;
		typedef Interaction                       _Interaction;

		_Parall parall_; // обеспечивает счет в параллельном режиме

		_LPComplex *complex_; // комплекс, который обслуживает Верлет таблица
		_Region *region_; // область взаимодействия
		real_t rcut_; // радиус взаимодействия
		real_t density_; // плотность частиц

	public:

		/**
		 * Создать пустую связанную с заданным комплексом Верлет таблицу
		 * @param complex молекулярный комплекс
		 */
		_Link_cell(_LPComplex *complex) : complex_(complex), region_(complex->get(REGION))
		{ parall_.init(this); }

		/**
		 * Зарезервировать память под Верлет таблицу при заданных параметрах расчета
		 * @param rskin - skin радиус
		 * @param density - плотность атомов (число частиц в 1 A**3)
		 */
		void resize(real_t rcut, real_t density)
		{
			rcut_ = rcut;
			density_ = density;
		}

		/**
		 * Сделать update таблицы, если какие-либо атомы ушли на большое расстояние
		 */
		void update(bool print=NO_PRINT) { update_(); }

		_E(real_t) dU__dX();

		real_t average(_I2T<RSKIN_>) const { return rcut_; }
		unsigned count(_I2T<REBUILD_>) const { return 0; }


	protected:

		void update_();

		_E(real_t) dU__dX_(unsigned pos);
		_E(real_t) dU__dX_(unsigned pos, unsigned pos__);

		/**
		 * Расчет пар для заданной ячейки в параллельном режиме
		 * @param wparam энергия пар в ячейке
		 * @param rparam номер ячейки
		 */
		void parall_dU__dX1_(void *wparam, void *rparam, void *params=0)
		{ *(_E(real_t) *)wparam += dU__dX_(((unsigned *)rparam)[0]); }

		/**
		 * Расчет пар для заданных ячеек в параллельном режиме
		 * @param wparam энергия пар в ячейке
		 * @param rparam номера ячейки
		 */
		void parall_dU__dX2_(void *wparam, void *rparam, void *params=0)
		{ *(_E(real_t) *)wparam += dU__dX_(((unsigned *)rparam)[0], ((unsigned *)rparam)[1]); }

		template <typename T1, typename T2>
		void read_data(__LJAtom<T1, T2> *ljatom1, __LJAtom<T1, T2> *ljatom2, unsigned pos) const;

		template <typename T1, typename T2>
		void save_data(__LJAtom<T1, T2> *ljatom1, __LJAtom<T1, T2> *ljatom2, unsigned pos);

		template <typename T1, typename T2>
		void read_data(__LJAtom<T1, T2> *ljatom1, __LJAtom<T1, T2> *ljatom2,
			unsigned pos1, unsigned pos2) const;

		template <typename T1, typename T2>
		void save_data(__LJAtom<T1, T2> *ljatom1, __LJAtom<T1, T2> *ljatom2,
			unsigned pos1, unsigned pos2);

	#ifdef USE_GONNET
		void parall_gonnet_ordering(void *wparam, void *rparam, void *params=0)
		{
			_Node *node = region_->get(NODE, *(unsigned *)rparam);
			node->make_gonnet_ordering();
		}
	#endif
	};

	#define TEMPLATE_HEADER  template <typename _LPComplex>
	#define TEMPLATE_ARG     _LPComplex

	TEMPLATE_HEADER
	INLINE void _Link_cell<TEMPLATE_ARG>
	::update_()
	{
		region_->clear(); // очистка от предыдущего расчета

		_Atom *atoms_ = complex_->get(ATOM);
		unsigned atom_count = complex_->count(ATOM);

		typename _LPComplex::ljatom_type ljatom;
		for (unsigned i=0,sz=atom_count; i<sz; i++)
		{
			ljatom.make(atoms_[i]);
			region_->insert(YES_PERIODIC, ljatom, atoms_[i].X);
		}

	#ifdef USE_GONNET
		{
			unsigned node_count = region_->count(NODE);
			unsigned nodes[node_count]; // идентификаторы рядов Верлет раблицы
			for (unsigned i=0; i<node_count; i++) nodes[i] = i;

			parall_.start_parallel_running(&_Link_cell::parall_gonnet_ordering, (unsigned *)0,
				node_count, nodes);
			parall_.wait();
		}
	#endif

	}

	TEMPLATE_HEADER
	INLINE _E(real_t) _Link_cell<TEMPLATE_ARG>
	::dU__dX()
	{
		_E(real_t) energy = 0.; // накопители энергии должны иметь расширенный тип

	#if !defined(SKIP_COUL_ENERGY) || !defined(SKIP_VDW_ENERGY)

		update();

		// подготовим данные для параллельного счета
		unsigned thread_count = parall_.thread_count();
		_E(real_t) energies[thread_count];
		for (unsigned i=0; i<thread_count; i++) energies[i] = 0.;

		{
			DECLARE_AS_RESTORING_ROUND_MODE
			_MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);

			unsigned set_count = region_->count(SET);
			for (unsigned iset=0; iset<set_count; iset++)
			{
				unsigned cnt = region_->count(SET, iset);

				if (region_->get(DIRECTION, iset) != region_->get(EMPTY_DIRECTION))
				{
					parall_.start_parallel_running(&_Link_cell::parall_dU__dX2_, &energies[0],
						cnt / 2, (std::pair<unsigned, unsigned>*)region_->get(PAIR, iset), (void*)region_);
				}
				else
				{
					parall_.start_parallel_running(&_Link_cell::parall_dU__dX1_, &energies[0],
						cnt, (unsigned *)region_->get(PAIR, iset), (void*)region_);
				}
				parall_.wait();
			}
		}
		for (unsigned i=0; i<thread_count; i++) energy += energies[i];

		_Atom *atoms_ = complex_->get(ATOM);
		const _Node *nodes = region_->get(NODE);
		for (unsigned pos=0,sz=region_->count(NODE); pos<sz; pos++)
		{
			const _Node &node = nodes[pos];
			typename _Node::const_iterator it = node.begin(), ite = node.end();
			for (; it!=ite; ++it) atoms_[it->insert_data].F += it->F;
		}

	#endif // !defined(SKIP_COUL_ENERGY) || !defined(SKIP_VDW_ENERGY)

		return energy;
	}

	TEMPLATE_HEADER
	INLINE _E(real_t) _Link_cell<TEMPLATE_ARG>
	::dU__dX_(unsigned pos)
	{
		const _Node *node = region_->get(NODE, pos);
		unsigned atom_count = node->count(ATOM);
		if (atom_count < 2) return (_E(real_t))0.;

		real_t rmax2 = sqr(_Interaction::interaction_radius());

		unsigned vsz = velement_count<vreal_t>(atom_count);
			// число элементов в массиве из векторных данных (для SSE оно кратно 4)

		//===============================================================================================
		// Алгоритм использует два массива SSE векторов для расчета. Первый массив хранит
		// значения плотно упакованными - каждый элемент SSE вектора соответствует разным значениям.
		// Второй массив хранит одинаковые сдублированные по 4 раза значения. Такая схема позволяет
		// максимально быстро создать и считать всевозможные пары, поскольку не требует переноса
		// данных по памяти или SHUFFLE операций внутри SSE векторов для получения перекрестных пар,
		// и дает почти 2 кратный выигрыш по сравнению с ними.
		//===============================================================================================
		__LJAtom<vreal_t, vint_t> ljatom1[vsz], ljatom2[atom_count];
		read_data(&ljatom1[0], &ljatom2[0], pos);

		vreal_t energy = 0.f;

		//===============================================================================================
		// Использована сложная схема суммирования, которая суммирует сперва наиболее удаленные пары и
		// в последнюю очередь более близкие пары. Это позволяет несколько увеличить точность расчетов.
		// Однако суммирование выполняется в 3 цикла, вместо обычных 2-х. Эффективность меняется в
		// пределах долей процента, поскольку число ячеек, для которых работает этот алгоритм в 26 раз
		// меньше, чем для алгоритма, рассчитывающего взаимодействия между соседними ячейками.
		// Данная схема написана для Интел процессора и не должна прямо переноситься на графические
		// процессоры типа Теслы. Схема выполняет суммирование по диагоналям матрицы взаимодействий, а
		// не по строкам этой матрицы, как это было бы естественно.
		//===============================================================================================
		for (unsigned i=1; i<vsz; i++) // проход по диагоналям матрицы взаимодействий
		{
			int ii__ = atom_count - vreal_t::size * i;
			for (unsigned j=0,ii=0; j<i; j++,ii++,ii__+=vreal_t::size)
			{
				__LJAtom<vreal_t, vint_t> &atom1 = ljatom1[ii];
				for (unsigned k=0; k<vreal_t::size; k++)
				{
					__LJAtom<vreal_t, vint_t> &atom2 = ljatom2[ii__ + k];
					LJ_CALCULATE(is_interaction, energy__, fx, fy, fz, rmax2, atom1, atom2);
					atom1.fx += fx; atom1.fy += fy; atom1.fz += fz;
					atom2.fx -= fx; atom2.fy -= fy; atom2.fz -= fz;
					energy += energy__;
				}
			}
		}

		for (unsigned ii=1; ii<vsz; ii++)
		{
			// Суммирование остатков, лежащих на диагонали матрицы (1, 2, 3, 4) -> (2, 3, 4, 1) и т.д.
			// делается в последнюю очередь, так как там находятся члены, даюшщие самые большие
			// результаты. Так пытаемся сохранить точность. Использовано дублирование взаимодействий,
			// так что энергия должна уменьшаться вдвое, а силы записываться только для одного атома пары.
			// Второй атом (точнее 4 атома, записанных в SSE вектор) вращаются ror(atom2),
			// чтобы обеспечить расчет всех требуемых пар.

			__LJAtom<vreal_t, vint_t> &atom1 = ljatom1[ii];
			__LJAtom<vreal_t, vint_t> atom2 = ljatom1[ii];
			for (unsigned i=1; i<vreal_t::size; i++)
			{
				ror1(atom2);
				LJ_CALCULATE(is_interaction, energy__, fx, fy, fz, rmax2, atom1, atom2);
				atom1.fx += fx; atom1.fy += fy; atom1.fz += fz;
				energy += 0.5 * energy__;
			}
		}

		unsigned n = atom_count - (vsz - 1) * vreal_t::size;
		for (unsigned l=0; l<n-1; l++)
		{
			__LJAtom<vreal_t, vint_t> &atom1 = ljatom2[l];
			for (unsigned l__=l+1; l__<n; l__++)
			{
				__LJAtom<vreal_t, vint_t> &atom2 = ljatom2[l__];
				LJ_CALCULATE(is_interaction, energy__, fx, fy, fz, rmax2, atom1, atom2);
				atom1.fx[0] += fx[0]; atom1.fy[0] += fy[0]; atom1.fz[0] += fz[0];
				energy[0] += energy__[0];
			}
		}

		// Сохранение насчитанных сил в ячейках региона.
		save_data(&ljatom1[0], &ljatom2[0], pos);

		return summarize(energy);
	}

	TEMPLATE_HEADER
	INLINE _E(real_t) _Link_cell<TEMPLATE_ARG>
	::dU__dX_(unsigned pos, unsigned pos__)
	{
		const _Node *node1 = region_->get(NODE, pos);
		const _Node *node2 = region_->get(NODE, pos__);
		unsigned atom_count1 = node1->count(ATOM);
		unsigned atom_count2 = node2->count(ATOM);

		if (atom_count1 == 0 || atom_count2 == 0) return (_E(real_t))0.;

		unsigned pairs = 0;
		vreal_t energy = 0.f;
		real_t rmax2 = sqr(_Interaction::interaction_radius());

		unsigned vsz = velement_count<vreal_t>(atom_count1);
			// число элементов в массиве из векторных данных (для SSE оно кратно 4)

		// Создаем 2 массива данных, причем в первом массиве начало заполняется пустыми элементами,
		// а во втором все данные дублируются во весь SSE вектор.	Такая техника позволяет избежать
		// проверок контроля последних элементов, при неполном заполнении SSE вектора данных (1 массив).
		__LJAtom<vreal_t, vint_t> ljatom1[vsz], ljatom2[atom_count2];

		// Загрузка данных из ячеек региона. Создание "псевдоатомов" на конце массива, связанных со
		// второй ячейкой, чтобы обеспечить возможность сдвига данных в SSE векторах. Установка
		// псевдоатомов в начале первого массива вместо его конца (более 20%) ускоряет счет при Гоннете.
		read_data(&ljatom1[0], &ljatom2[0], pos, pos__);

		// kk__ является индексом последнего атома, до которого имеет смысл считать взаимодействия.
		// Так как атомы упорядочены по расстояниям, то атомы за kk__ индексом слишком далеки.
		// Параметр настраивается перед счетом в Гоннет варианте через поиск вперед, стартуя от 0.
	#ifdef USE_GONNET
		unsigned kk__ = 0;
	#else
		unsigned kk__ = atom_count2;
	#endif

		for (unsigned ii=0; ii<vsz; ii++)
		{
			__LJAtom<vreal_t, vint_t> &atom1 = ljatom1[ii];
		#ifdef USE_GONNET
			// Делаем поиск вперед, чтобы найти самые дальние взаимодействующие атомы.
			// Хэш проверяем для максимально близких атомов SSE группы, а не для всех ее атомов,
			// что дает небольшой выигрыш по эффективности.
		#ifdef USE_VECTORIZATION
			real_t hash1 = ljatom1[ii].hash[3];
		#else
			real_t hash1 = ljatom1[ii].hash[0];
		#endif
			while (sqr(ljatom2[kk__].hash[0] - hash1) < rmax2 && kk__ < atom_count2) kk__++;

		#endif

			// Счет в обратную сторону от самых дальних атомов, чтобы сохранить точность расчетов.
			for (int ii__=kk__ - 1; ii__>=0; ii__--)
			{
				__LJAtom<vreal_t, vint_t> &atom2 = ljatom2[ii__];
				LJ_CALCULATE(is_interaction, energy__, fx, fy, fz, rmax2, atom1, atom2);
				pairs += is_interaction.count();
				atom1.fx += fx; atom1.fy += fy; atom1.fz += fz;
				atom2.fx -= fx; atom2.fy -= fy; atom2.fz -= fz;
				energy += energy__;
			}
		}

		// Сохранение насчитанных сил в ячейках региона.
		save_data(&ljatom1[0], &ljatom2[0], pos, pos__);

		return summarize(energy);
	}

	TEMPLATE_HEADER
	template <typename T1, typename T2>
	INLINE void _Link_cell<TEMPLATE_ARG>
	::read_data(__LJAtom<T1, T2> *ljatom1, __LJAtom<T1, T2> *ljatom2, unsigned pos) const
	{
		const _Node *node = region_->get(NODE, pos);
		unsigned atom_count = node->count(ATOM);
		const typename _LPComplex::ljatom_type *atom__ = node->get(ATOM);

		unsigned vsz = velement_count<T1>(atom_count);
			// число элементов в массиве из векторных данных (для SSE оно кратно 4)

		unsigned n = T1::size;
		unsigned sz = vsz * n;
		unsigned empty = sz - atom_count;

		// сделаем пустое начало в 1-м массиве
		molkern::read_empty(ljatom1[0], atom__);

		// заполним 1-й массив уникальным атомом в каждой позиции SSE вектора
		for (unsigned i=empty, ndx=0; i<sz; i++, ndx++)
			molkern::read_data(i % n, ljatom1[i / n], atom__ + ndx);

		// заполним 2-й массив дублируя атом по всем позициям SSE вектора
		for (unsigned i=0, ndx=0; i<atom_count; i++, ndx++)
			molkern::read_data(ljatom2[i], atom__ + ndx);
	}

	TEMPLATE_HEADER
	template <typename T1, typename T2>
	INLINE void _Link_cell<TEMPLATE_ARG>
	::save_data(__LJAtom<T1, T2> *ljatom1, __LJAtom<T1, T2> *ljatom2, unsigned pos)
	{
		typedef typename T1::value_type real_t;
		typedef typename T2::value_type int_t;

		_Node *nodes_ = region_->get(NODE);
		unsigned atom_count = nodes_[pos].count(ATOM);

		unsigned n = T1::size; // число компонент в sse векторе
		unsigned vsz = velement_count<T1>(atom_count);
			// число элементов в массиве из векторных данных (для SSE оно кратно 4)

		unsigned sz1 = vsz * n;
		unsigned empty = sz1 - atom_count;

		typename _LPComplex::ljatom_type *atom = nodes_[pos].get(ATOM);
		for (unsigned i=empty, ndx=0; i<sz1; i++, ndx++)
		{
			__LJAtom<T1, T2> &ljatom = ljatom1[i / n];
			unsigned m = i % n;
			real_t *f = &(atom + ndx)->F[0];
			f[0] += (real_t)*((real_t *)&ljatom.fx + m);
			f[1] += (real_t)*((real_t *)&ljatom.fy + m);
			f[2] += (real_t)*((real_t *)&ljatom.fz + m);
		}

		for (unsigned i=0, ndx=0; i<atom_count; i++, ndx++)
		{
			__LJAtom<T1, T2> &ljatom = ljatom2[i];
			real_t *f = &(atom + ndx)->F[0];
			f[0] += (real_t)summarize(ljatom.fx);
			f[1] += (real_t)summarize(ljatom.fy);
			f[2] += (real_t)summarize(ljatom.fz);
		}
	}

	TEMPLATE_HEADER
	template <typename T1, typename T2>
	INLINE void _Link_cell<TEMPLATE_ARG>
	::read_data(__LJAtom<T1, T2> *ljatom1, __LJAtom<T1, T2> *ljatom2,
		unsigned pos1, unsigned pos2) const
	{
		const _Node *node1 = region_->get(NODE, pos1);
		const _Node *node2 = region_->get(NODE, pos2);
		unsigned atom_count1 = node1->count(ATOM);
		unsigned atom_count2 = node2->count(ATOM);

	#ifdef USE_GONNET
		vector_t direction = region_->translation_vector(pos1, pos2);

		unsigned dircode = Gonnet_vector::make_direction_code(direction);
		real_t w = Gonnet_vector::make_direction_shift(direction);

		const std::pair<unsigned, real_t> *gonnet_pair1 = node1->get_gonnet_pair(dircode);
		const std::pair<unsigned, real_t> *gonnet_pair2 = node2->get_gonnet_pair(dircode);
	#endif

		unsigned n = T1::size; // число компонент в sse векторе
		unsigned vsz = velement_count<T1>(atom_count1);
			// число элементов в массиве из векторных данных (для SSE оно кратно 4)

		unsigned sz1 = vsz * n;
		unsigned empty = sz1 - atom_count1;

		const typename _LPComplex::ljatom_type *atom__ = node1->get(ATOM);

		// сделаем пустое начало в 1-м массиве
		molkern::read_empty(ljatom1[0], atom__);

		// заполним 1-й массив уникальным атомом в каждой позиции SSE вектора
		for (unsigned i=empty, ndx=0; i<sz1; i++, ndx++)
		{
			__LJAtom<T1, T2> &ljatom = ljatom1[i / n];
			unsigned m = i % n;
		#ifdef USE_GONNET
			ndx = gonnet_pair1[i - empty].first;
			*((typename T1::value_type *)&ljatom.hash + m) = gonnet_pair1[i - empty].second;
		#endif
			molkern::read_data(m, ljatom, atom__ + ndx);
		}

		for (unsigned i=0; i<empty; i++)
		{
			__LJAtom<T1, T2> &ljatom = ljatom1[i / n];
			unsigned m = i % n;
		#ifdef USE_GONNET
			*((typename T1::value_type *)&ljatom.hash + m) = *((real_t *)&ljatom.hash + empty);
		#endif
		}

		// заполним 2-й массив дублируя атом по всем позициям SSE вектора
		atom__ = node2->get(ATOM);
		for (unsigned i=0, ndx=0; i<atom_count2; i++, ndx++)
		{
			__LJAtom<T1, T2> &ljatom = ljatom2[i];
		#ifdef USE_GONNET
			ndx = gonnet_pair2[i].first;
			ljatom.hash = w + gonnet_pair2[i].second;
		#endif
			molkern::read_data(ljatom, atom__ + ndx);
		}
	}

	TEMPLATE_HEADER
	template <typename T1, typename T2>
	INLINE void _Link_cell<TEMPLATE_ARG>
	::save_data(__LJAtom<T1, T2> *ljatom1, __LJAtom<T1, T2> *ljatom2, unsigned pos1, unsigned pos2)
	{
		_Node *node1 = region_->get(NODE, pos1);
		_Node *node2 = region_->get(NODE, pos2);
		vector_t direction = region_->translation_vector(pos1, pos2);

	#ifdef USE_GONNET
		unsigned dircode = Gonnet_vector::make_direction_code(direction);
		const std::pair<unsigned, real_t> *gonnet_pair1 = node1->get_gonnet_pair(dircode);
		const std::pair<unsigned, real_t> *gonnet_pair2 = node2->get_gonnet_pair(dircode);
	#endif

		unsigned atom_count1 = node1->count(ATOM);
		unsigned atom_count2 = node2->count(ATOM);

		unsigned n = T1::size; // число компонент в sse векторе
		unsigned vsz = velement_count<T1>(atom_count1);
			// число элементов в массиве из векторных данных (для SSE оно кратно 4)

		unsigned sz1 = vsz * n;
		unsigned empty = sz1 - atom_count1;

		typename _LPComplex::ljatom_type *atom = node1->get(ATOM);
		for (unsigned i=empty, ndx=0; i<sz1; i++, ndx++)
		{
			__LJAtom<T1, T2> &ljatom = ljatom1[i / n];
			unsigned m = i % n;
		#ifdef USE_GONNET
			ndx = gonnet_pair1[i - empty].first;
		#endif
			real_t *f = &(atom + ndx)->F[0];
			f[0] += *((real_t *)&ljatom.fx + m);
			f[1] += *((real_t *)&ljatom.fy + m);
			f[2] += *((real_t *)&ljatom.fz + m);
		}

		//atom = nodes_[pos2].get(ATOM);
		atom = node2->get(ATOM);
		for (unsigned i=0, ndx=0; i<atom_count2; i++, ndx++)
		{
			__LJAtom<T1, T2> &ljatom = ljatom2[i];
		#ifdef USE_GONNET
			ndx = gonnet_pair2[i].first;
		#endif
			real_t *f = &(atom + ndx)->F[0];
			f[0] += (real_t)summarize(ljatom.fx);
			f[1] += (real_t)summarize(ljatom.fy);
			f[2] += (real_t)summarize(ljatom.fz);
		}
	}


#undef TEMPLATE_HEADER
#undef TEMPLATE_ARG

}
#endif
