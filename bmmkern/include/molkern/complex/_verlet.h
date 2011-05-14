#ifndef VERLET__09ED1116_C3FE_55a4_CA63_F84536C80D01__H
#define VERLET__09ED1116_C3FE_55a4_CA63_F84536C80D01__H

#include "molkern/__moldefs.h"
#include "molkern/complex/_parallel.h"
#include "molkern/complex/_region.h"

/// определяется при тестировании Верлет таблицы, требует огромного количества времени
// #define VERLET_DEBUG

namespace molkern
{
	using namespace prgkern;

	/**
	 * Среднее число соседей для выбранной частицы (в полусфере).
	 * @param density плотность частиц (частиц в A**3)
	 * @param radius радиус сферы
	 * @return число соседей
	 */
	INLINE unsigned count(_I2T<NEIGHBOR_>, real_t radius, real_t density)
	{
		return (unsigned)ceil((4. / 6.) * M_PI * density * cube(radius));
	}

	/**
	 * Верлет таблица используется для того, чтобы хранить пары атомов, которые находятся на
	 * расстоянии r < rskin друг от друга. Она используется только для ускорения расчетов
	 * ближних невалентных взаимодействий.
	 *
	 * Верлет таблица является редким массивом. Такой массив обычно не имеет богатой функциональности
	 * для того, чтобы быть выделенным отдельно. Но, с другой стороны, эта таблица не является
	 * обязательным атрибутом молекулярного комплекса, то есть и ее реализация должна быть отдельной.
	 *
	 * Таким образом, Верлет таблица может рассматриваться как некий инструмент, оперирующий с
	 * данными комплекса, который может быть использован, а может и нет. Если он используется, то
	 * он использует всю функциональность самого комплекса для эффективности. То есть, он максимально
	 * завязан на комплекс, а комплекс на него не завязан.
	 */
	template <typename _LPComplex>
	class _Verlet_table : public near_range_integrator<_LPComplex>
	{
		typedef typename _LPComplex::region_type  _Region;
		typedef typename _LPComplex::atom_type    _Atom;
		typedef Parall_<_Verlet_table>            _Parall;
		typedef Interaction                       _Interaction;

		typedef __LJAtom<vreal_t, vint_t>         _LJAtom;
		typedef typename _Region::value_type      _GridAtom;
		typedef typename _Region::key_type        _Key;
		typedef typename _Region::index_type      _Index;
		typedef typename _Region::node_type       _Node;

		std::vector<std::vector<unsigned> > verlet_table_; // Верлет таблица
		std::vector<vector_t> verlet_offsets_; // координаты атомов для контроля перестройки таблицы
		std::vector<std::vector<vector_t> > verlet_forces_;
		_Parall parall_; // обеспечивает счет в параллельном режиме

		_LPComplex *complex_; // комплекс, который обслуживает Верлет таблица
		std::auto_ptr<_Region> region_; // область взаимодействия с гридом по rskin

		real_t rskin_; // текущий skin радиус
		real_t density_; // плотность частиц
		unsigned rebuild_verlet_table_count_; // число перестроек Верлет таблицы
		real_t rskin_average_; // среднее <rskin> за прошедшее число шагов

#ifdef USE_VERLET_TABLE_TUNE
		// Верлет таблица позволяет во время счета подстраивать rskin чтобы минимизировать полное
		// время расчета. Для этого производятся измерения времени отработки между циклами перестройки
		// таблицы и минимизируется вреднее время на один шаг динамики.

		system_time_t prev_system_time_; // предыдущее астрономическое время (число секунд)
		model_time_t prev_model_time_; // предыдущее модельное время (число fs)
		real_t prev_average_time_; // предыдущее среднее время на шаг (сек)
		real_t rskin_step_; // текущее направление и величина изменения rskin
#endif

#ifdef USE_VERLET_TABLE_GROUP

		void compress_(bool print=NO_PRINT);
		void parall_compress(void *wparam, void *rparam, void *params);
		void compress_(unsigned *share, const _Region *grid, unsigned igroup);

		struct Group_
		{
			unsigned elem_count_;
			unsigned group_elems_[MAX_GROUP_ELEM]; // объединяемые в группу ряды таблицы
			std::vector<unsigned> share_atoms_; // разделяемые атомы группы
			std::vector<unsigned> uniq_atoms_[MAX_GROUP_ELEM];
		};
		std::vector<Group_> verlet_groups_;
		std::vector<std::vector<unsigned> > share_;
#endif

	public:

		/**
		 * Создать пустую связанную с заданным комплексом Верлет таблицу
		 * @param complex молекулярный комплекс
		 */
		_Verlet_table(_LPComplex *complex) : complex_(complex) { parall_.init(this); }

		/**
		 * Зарезервировать память под Верлет таблицу при заданных параметрах расчета
		 * @param rskin - skin радиус
		 * @param density - плотность атомов (число частиц в 1 A**3)
		 */
		void resize(real_t rskin, real_t density);

		/**
		 * Сделать update таблицы, если какие-либо атомы ушли на большое расстояние
		 */
		void update(bool print=NO_PRINT)
		{
			if (has_to_be_updated_())
			{
				unsigned pairs_count = update_();

			#ifdef USE_VERLET_TABLE_GROUP
				compress_(print);
			#endif

				if (print)
				{
					std::string msg = make_string("Have been rebuilt -> verlet table (%d), rskin = %f",
						pairs_count, (float)rskin_);
					PRINT_MESSAGE(msg);
				}
			}
		}

		_E(real_t) dU__dX();

		real_t average(_I2T<RSKIN_>) const { return rskin_average_; }
		unsigned count(_I2T<REBUILD_>) const { return rebuild_verlet_table_count_; }


	protected:

		unsigned update_();
		unsigned update_(unsigned pos);
		unsigned update_(unsigned pos, unsigned pos__);
		_E(real_t) dU__dX_(vector_t *F, unsigned irow);


		bool has_to_be_updated_() const;

		/**
		 * Построение пар для заданной ячейки в параллельном режиме
		 * @param wparam количество пар в ячейке
		 * @param rparam номер ячейки
		 * @param params грид
		 */
		void parall_build1_(void *wparam, void *rparam, void *params)
		{
			*(unsigned *)wparam += update_(((unsigned *)rparam)[0]);
		}

		/**
		 * Построение пар для заданных ячеек в параллельном режиме
		 * @param wparam количество пар в ячейке
		 * @param rparam номера ячейки
		 * @param params грид
		 */
		void parall_build2_(void *wparam, void *rparam, void *params)
		{
			*(unsigned *)wparam += update_(((unsigned *)rparam)[0], ((unsigned *)rparam)[1]);
		}

		void parall_dU__dX_(void *wparam, void *rparam, void *params)
		{
			unsigned irow = *(unsigned *)rparam;

#ifdef USE_VERLET_TABLE_GROUP
			const Group_ &group = verlet_groups_[irow];
			if (group.elem_count_ == 0) return;
#else
			if (verlet_table_[irow].size() == 0) return;
#endif

			typedef std::pair<_E(real_t), vector_t *>  wparam_t;
			wparam_t &wparams = *(wparam_t *)wparam;
			vector_t *F__ = wparams.second;

			wparams.first += dU__dX_(F__, irow);
		}

		template <typename T1, typename T2>
		void read_data(__LJAtom<T1, T2> *ljatoms, const std::vector<unsigned> &row);

		template <typename T1, typename T2>
		void read_data(__LJAtom<T1, T2> *ljatom1, __LJAtom<T1, T2> *ljatom2,
			unsigned pos1, unsigned pos2) const;

		template <typename T1, typename T2>
		void save_data(vector_t *F, const __LJAtom<T1, T2> *ljatoms, const std::vector<unsigned> &row);

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
	INLINE void _Verlet_table<TEMPLATE_ARG>
	::resize(real_t rskin, real_t density)
	{
		rskin_ = rskin;
		density_ = density;

		unsigned atom_count = complex_->count(ATOM);
		verlet_table_.resize(atom_count);
		unsigned n = molkern::count(NEIGHBOR, rskin, density);
		for (unsigned i=0; i<atom_count; i++) verlet_table_[i].resize(n);

		unsigned thread_count = parall_.thread_count();
		verlet_forces_.resize(thread_count);
		for (unsigned i=0; i<thread_count; i++)
			verlet_forces_[i].resize(atom_count, 0.f);

		verlet_offsets_.resize(atom_count, 0.);
		rebuild_verlet_table_count_ = 0;
		rskin_average_ = rskin;

		//--------------------------------------------------------------------
		//      создадим внутренний грид с ячейками нужного размера
		//--------------------------------------------------------------------
		const _Region *region = complex_->get(REGION);

		vector_t T = region->translation_vector();
		real_t rcut = _Interaction::interaction_radius();
		_Index region_sz = region->box_size();

		real_t compress_factor = rcut / rskin_;
		_Index grid_sz
		(
			(unsigned)ceil(region_sz[0] * compress_factor),
			(unsigned)ceil(region_sz[1] * compress_factor),
			(unsigned)ceil(region_sz[2] * compress_factor)
		);
		region_ = std::auto_ptr<_Region>(new _Region(grid_sz, rskin_, &T));
	}

	TEMPLATE_HEADER
	INLINE bool _Verlet_table<TEMPLATE_ARG>
	::has_to_be_updated_() const
	{
		real_t max2a = 0, max2b = 0; // квадраты максимальных смещений частиц A и B

		unsigned atom_count = complex_->count(ATOM);
		const _Atom *atoms = complex_->get(ATOM);
		for (unsigned i=0; i<atom_count; i++)
		{
			vector_t offs = atoms[i].X - verlet_offsets_[i];
			real_t d2 = scalar_product(offs, offs);
			if (d2 > max2b)
			{
				if (d2 > max2a) { max2b = max2a; max2a = d2; }
				else { max2b = d2; }
			}
		}

		if (sqrt(max2a) + sqrt(max2b) > rskin_ - _Interaction::interaction_radius())
			return true;

		return false;
	}

	TEMPLATE_HEADER
	inline unsigned _Verlet_table<TEMPLATE_ARG>
	::update_()
	{
		typedef typename _Region::value_type          _LJAtom;
		typedef typename _Region::index_type          _Index;
		typedef typename _Region::key_type            _Key;
		typedef typename _Region::region_vector_type  region_vector_t;
		typedef typename _Region::region_key_type     region_key_t;

	#ifdef USE_VERLET_TABLE_TUNE
		//-----------------------------------------------------------------------
		//      минимизация времени выполнения путем оптимизации rskin
		//-----------------------------------------------------------------------
		{
			system_time_t curr_real_time = current_time();
			model_time_t curr_model_time = current_model_time();

			if (curr_model_time == 0)
			{
				prev_model_time_ = curr_model_time;
				prev_system_time_ = curr_real_time;
				prev_average_time_ = infinity<real_t>();
				rskin_step_ = DEFAULT_RSKIN_STEP;
				rskin_average_ = 0;
			}
			else
			{
				unsigned nt = curr_model_time - prev_model_time_;
				double curr_average_time = double(curr_real_time - prev_system_time_) / nt;

				if (curr_average_time > prev_average_time_) rskin_step_ = -rskin_step_;
					// меняем направление в случае замедления

				rskin_ += rskin_step_;
				if (rskin_ < _Interaction::interaction_radius())
					rskin_ = _Interaction::interaction_radius();
					// запрет для избежания потери пар

				rskin_average_ = (rskin_average_ * prev_model_time_ + rskin_ * nt) / curr_model_time;
					// рассчитываем среднее значение rskin для контроля

				prev_model_time_ = curr_model_time;
				prev_system_time_ = curr_real_time;
				prev_average_time_ = curr_average_time;
			}
		}
	#endif

		unsigned pairs_count = 0;
		rebuild_verlet_table_count_++;

		//--------------------------------------------------------------------------
		//                  очистка от предыдущей итерации
		//--------------------------------------------------------------------------
		unsigned atom_count = complex_->count(ATOM);
		for (unsigned i=0; i<atom_count; i++) verlet_table_[i].resize(0);

		const _Atom *atoms = complex_->get(ATOM);
		for (unsigned i=0; i<atom_count; i++) verlet_offsets_[i] = atoms[i].X;

		//--------------------------------------------------------------------------
		//                задание нового грида с расширенными ячейками
		//--------------------------------------------------------------------------
		region_->clear();
		_LJAtom ljatom;
		for (unsigned i=0; i<atom_count; i++)
		{
			ljatom.make(atoms[i]);
			region_->insert(YES_PERIODIC, ljatom, atoms[i].X);
		}

	#ifdef USE_GONNET
		{
			unsigned node_count = region_->count(NODE);
			unsigned nodes[node_count]; // идентификаторы рядов Верлет раблицы
			for (unsigned i=0; i<node_count; i++) nodes[i] = i;

			parall_.start_parallel_running(&_Verlet_table::parall_gonnet_ordering, (unsigned *)0,
				node_count, nodes);
			parall_.wait();
		}
	#endif

		//--------------------------------------------------------------------------
		//                        счет по ячейкам
		//--------------------------------------------------------------------------
		{
			DECLARE_AS_RESTORING_ROUND_MODE
			_MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);

			unsigned thread_count = parall_.thread_count();
			std::vector<unsigned> pairs_counts(thread_count, 0);

			unsigned set_count = region_->count(SET);
			for (unsigned iset=0; iset<set_count; iset++)
			{
				unsigned cnt = region_->count(SET, iset);

				if (region_->get(DIRECTION, iset) != region_->get(EMPTY_DIRECTION))
				{
					parall_.start_parallel_running(&_Verlet_table::parall_build2_, &pairs_counts[0],
						cnt / 2, (std::pair<unsigned, unsigned>*)region_->get(PAIR, iset));
				}
				else
				{
					parall_.start_parallel_running(&_Verlet_table::parall_build1_, &pairs_counts[0],
						cnt, (unsigned *)region_->get(PAIR, iset));
				}
				parall_.wait();
			}
			for (unsigned i=0; i<thread_count; i++) pairs_count += pairs_counts[i];
		}

#ifdef VERLET_DEBUG
		{
			unsigned atom_count = complex_->count(ATOM);

			std::vector<unsigned> sh(atom_count, 0);
			for (unsigned i=0; i<atom_count; i++)
			{
				std::vector<unsigned> &row = verlet_table_[i];
				for (unsigned k=0; k<row.size(); k++) sh[row[k]]++;

				for (unsigned k=0; k<atom_count; k++)
					if (sh[k] > 1)
					{
						std::string msg =
							make_string("the same atom %d meets more then once in verlet table row %d", k, i);
						PRINT_ERR(msg);
					}
				for (unsigned k=0; k<atom_count; k++) sh[k] = 0;
			}
		}
#endif

		return pairs_count;
	}

	TEMPLATE_HEADER
	INLINE unsigned _Verlet_table<TEMPLATE_ARG>
	::update_(unsigned pos)
	{
		typedef typename _Region::node_type  _Node;
		unsigned pairs_count = 0;

		const _Node *node = region_->get(NODE, pos);
		unsigned atom_count = node->count(ATOM);
		if (atom_count < 2) return pairs_count;

		real_t rskin2 = sqr(region_->get(EDGE_LENGTH));
		unsigned vsz = velement_count<vreal_t>(atom_count);
			// число элементов в массиве из векторных данных (для SSE оно кратно 4)

		typedef __LJAtom<vreal_t, vint_t>  _LJAtom;
		typedef typename _Region::value_type  _GridAtom;
		_LJAtom ljatom1[vsz], ljatom2[atom_count];
		{
			unsigned n = vreal_t::size; // число компонент в sse векторе
			unsigned vsz = velement_count<vreal_t>(atom_count);
				// число элементов в массиве из векторных данных (для SSE оно кратно 4)

			unsigned sz1 = vsz * n;
			unsigned empty = sz1 - atom_count;

			const _GridAtom *atom__ = node->get(ATOM);
			for (unsigned i=empty, ndx=0; i<sz1; i++, ndx++)
			{
				_LJAtom &ljatom = ljatom1[i / n];
				unsigned m = i % n;
				const _GridAtom *atom = atom__ + ndx;
				*((int *)&ljatom.connect_data + m) = (int)atom->connect_data;
				*((int *)&ljatom.insert_data + m) = (int)atom->insert_data;
				*((real_t *)&ljatom.x + m) = (real_t)atom->X[0];
				*((real_t *)&ljatom.y + m) = (real_t)atom->X[1];
				*((real_t *)&ljatom.z + m) = (real_t)atom->X[2];
			}

			for (unsigned i=0; i<empty; i++)
			{
				_LJAtom &ljatom = ljatom1[i / n];
				unsigned m = i % n;
				*((int *)&ljatom.connect_data + m) = 0;
				*((int *)&ljatom.insert_data + m) = 0;
				*((real_t *)&ljatom.x + m) = *((real_t *)&ljatom.x + empty);
				*((real_t *)&ljatom.y + m) = *((real_t *)&ljatom.y + empty);
				*((real_t *)&ljatom.z + m) = *((real_t *)&ljatom.z + empty);
			}

			atom__ = node->get(ATOM);
			for (unsigned i=0, ndx=0; i<atom_count; i++, ndx++)
			{
				_LJAtom &ljatom = ljatom2[i];
				const _GridAtom *atom = atom__ + ndx;
				ljatom.connect_data = (int)atom->connect_data;
				ljatom.insert_data = (int)atom->insert_data;
				ljatom.x = (real_t)atom->X[0];
				ljatom.y = (real_t)atom->X[1];
				ljatom.z = (real_t)atom->X[2];
			}
		}

		vreal_t R[3];
		vector_t T = region_->translation_vector();
		vector_t _1T = vector_t(1./T[0], 1./T[1], 1./T[2]);
			// вектора трансляции используем от региона комплекса, а не текущего грида, в котором
			// они уменьшены из-за коррекции шага грида

		//---------------------------------------------------------------------------------------
		// Использована сложная схема суммирования, которая ищет сперва наиболее удаленные пары и
		// в последнюю очередь более близкие пары. Это позволяет несколько увеличить точность
		// расчетов. Однако суммирование выполняется в 3 цикла, вместо обычных 2-х.
		//---------------------------------------------------------------------------------------
		for (unsigned i=1; i<vsz; i++) // проход по диагоналям матрицы взаимодействий
		{
			int ii__ = atom_count - vreal_t::size * i;
			for (unsigned j=0,ii=0; j<i; j++,ii++,ii__+=vreal_t::size)
			{
				__LJAtom<vreal_t, vint_t> &atom1 = ljatom1[ii];
				for (unsigned k=0; k<vreal_t::size; k++)
				{
					__LJAtom<vreal_t, vint_t> &atom2 = ljatom2[ii__ + k];

					R[0] = atom2.x - atom1.x;
					R[1] = atom2.y - atom1.y;
					R[2] = atom2.z - atom1.z;
					R[0] -= round(R[0] * _1T[0]) * T[0];
					R[1] -= round(R[1] * _1T[1]) * T[1];
					R[2] -= round(R[2] * _1T[2]) * T[2];

					vbool_t is_interaction = sqr(R[0]) + sqr(R[1]) + sqr(R[2]) < rskin2;
					if (is_interaction == false) continue;
					CONTROL_CONTACT(is_interaction, atom1, atom2);

					for (unsigned i=0; i<vreal_t::size; i++)
					{
						if (is_interaction[i])
						{
							unsigned k = atom1.insert_data[i], k__ = atom2.insert_data[i];
							if (k > k__) std::swap(k, k__);
							verlet_table_[k].push_back(k__);
							pairs_count++;
						}
					}
				}
			}
		}

		for (unsigned ii=0; ii<vsz; ii++)
		{
			// Суммирование остатков, лежащих на диагонали матрицы (1, 2, 3, 4) -> (2, 3, 4, 1) и т.д.
			// делается в последнюю очередь, так как там находятся члены, даюшщие самые большие
			// результаты. Так пытаемся сохранить точность.

			__LJAtom<vreal_t, vint_t> &atom = ljatom1[ii]; real_t R[3];
			for (unsigned i=0; i<vreal_t::size; i++)
			{
				for (unsigned i__=i+1; i__<vreal_t::size; i__++)
				{
					R[0] = atom.x[i__] - atom.x[i];
					R[1] = atom.y[i__] - atom.y[i];
					R[2] = atom.z[i__] - atom.z[i];
					R[0] -= round(R[0] * _1T[0]) * T[0];
					R[1] -= round(R[1] * _1T[1]) * T[1];
					R[2] -= round(R[2] * _1T[2]) * T[2];

					vbool_t is_interaction = sqr(R[0]) + sqr(R[1]) + sqr(R[2]) < rskin2;
					if (is_interaction == false) continue;

					if (has_contact(atom.insert_data[i], atom.insert_data[i__],
						atom.connect_data[i], atom.connect_data[i__])) continue;

					if (atom.connect_data[i  ] == 0) continue;
					if (atom.connect_data[i__] == 0) continue;

					unsigned k = atom.insert_data[i], k__ = atom.insert_data[i__];
					if (k > k__) std::swap(k, k__);
					verlet_table_[k].push_back(k__);
					pairs_count++;
				}
			}
		}
		return pairs_count;
	}

	TEMPLATE_HEADER
	INLINE unsigned _Verlet_table<TEMPLATE_ARG>
	::update_(unsigned pos, unsigned pos__)
	{
		typedef __LJAtom<vreal_t, vint_t>     _LJAtom;
		typedef typename _Region::value_type  _GridAtom;
		typedef typename _Region::key_type    _Key;
		typedef typename _Region::index_type  _Index;
		typedef typename _Region::node_type   _Node;

		const _Node *node1 = region_->get(NODE, pos);
		const _Node *node2 = region_->get(NODE, pos__);
		unsigned atom_count1 = node1->count(ATOM);
		unsigned atom_count2 = node2->count(ATOM);
		if (atom_count1 == 0 || atom_count2 == 0) return 0;

		unsigned pairs_count = 0;
		real_t rskin2 = sqr(region_->get(EDGE_LENGTH));

		unsigned vsz = velement_count<vreal_t>(atom_count1);
			// число элементов в массиве из векторных данных (для SSE оно кратно 4)

		// Создаем 2 массива данных, причем в первом массиве начало заполняется пустыми элементами,
		// а во втором все данные дублируются во весь SSE вектор.	Такая техника позволяет избежать
		// проверок контроля последних элементов, при неполном заполнении SSE вектора данных (1 массив).
		__LJAtom<vreal_t, vint_t> ljatom1[vsz], ljatom2[atom_count2];

		vreal_t R[3];
		vector_t T = region_->translation_vector();
		vector_t _1T = vector_t(1./T[0], 1./T[1], 1./T[2]);
			// вектора трансляции используем от региона комплекса, а не текущего грида, в котором
			// они уменьшены из-за коррекции шага грида

		bool use_gonnet = !region_->is_edge(pos, pos__);
		if (use_gonnet)
		{
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
				while (sqr(ljatom2[kk__].hash[0] - hash1) < rskin2 && kk__ < atom_count2) kk__++;

			#endif

				// Счет в обратную сторону от самых дальних атомов, чтобы сохранить точность расчетов.
				for (int ii__=kk__ - 1; ii__>=0; ii__--)
				{
					__LJAtom<vreal_t, vint_t> &atom2 = ljatom2[ii__];
					R[0] = atom2.x - atom1.x;
					R[1] = atom2.y - atom1.y;
					R[2] = atom2.z - atom1.z;
					R[0] -= round(R[0] * _1T[0]) * T[0];
					R[1] -= round(R[1] * _1T[1]) * T[1];
					R[2] -= round(R[2] * _1T[2]) * T[2];

					vbool_t is_interaction = sqr(R[0]) + sqr(R[1]) + sqr(R[2]) < rskin2;
					if (is_interaction == false) continue;
					CONTROL_CONTACT(is_interaction, atom1, atom2);

					for (unsigned i=0; i<vreal_t::size; i++)
					{
						if (is_interaction[i])
						{
							unsigned k = atom1.insert_data[i], k__ = atom2.insert_data[i];
							if (k > k__) std::swap(k, k__);
							verlet_table_[k].push_back(k__);
							pairs_count++;
						}
					}
				}
			}
		}
		else // не используем Гоннет
		{
			piecewise_vector<const _GridAtom> atoms1(node1->get(ATOM), node1->count(ATOM));
			piecewise_vector<const _GridAtom> atoms2(node2->get(ATOM), node2->count(ATOM));

			_Index key = region_->index3(pos);
			_Index key__ = region_->index3(pos__);
			_Key direction(key__[0] - key[0], key__[1] - key[1], key__[2] - key[2]);

			region_->make_piecewise_vector(&atoms1, direction, pos);
			region_->make_piecewise_vector(&atoms2, direction, pos__);

			//=============================================================================================
			//                         Загрузка данных из ячеек региона.
			//=============================================================================================
			unsigned n = vreal_t::size; // число компонент в sse векторе
			unsigned sz1 = vsz * n;
			unsigned empty = sz1 - atom_count1;

			typedef typename piecewise_vector<const _GridAtom>::iterator  _Iterator;

			_Iterator it = atoms1.begin();
			for (unsigned i=empty, ndx=0; i<sz1; i++, ndx++)
			{
				_LJAtom &ljatom = ljatom1[i / n];
				unsigned m = i % n;
				const _GridAtom &atom = *it; ++it;

				*((int *)&ljatom.connect_data + m) = (int)atom.connect_data;
				*((int *)&ljatom.insert_data + m) = (int)atom.insert_data;
				*((real_t *)&ljatom.x + m) = (real_t)atom.X[0];
				*((real_t *)&ljatom.y + m) = (real_t)atom.X[1];
				*((real_t *)&ljatom.z + m) = (real_t)atom.X[2];
			}

			for (unsigned i=0; i<empty; i++)
			{
				_LJAtom &ljatom = ljatom1[i / n];
				unsigned m = i % n;
				*((int *)&ljatom.connect_data + m) = 0;
				*((int *)&ljatom.insert_data + m) = 0;
				*((real_t *)&ljatom.x + m) = *((real_t *)&ljatom.x + empty);
				*((real_t *)&ljatom.y + m) = *((real_t *)&ljatom.y + empty);
				*((real_t *)&ljatom.z + m) = *((real_t *)&ljatom.z + empty);
			}

			it = atoms2.begin();
			for (unsigned i=0, ndx=0; i<atom_count2; i++, ndx++)
			{
				_LJAtom &ljatom = ljatom2[i];
				const _GridAtom &atom = *it; ++it;
				ljatom.connect_data = (int)atom.connect_data;
				ljatom.insert_data = (int)atom.insert_data;
				ljatom.x = (real_t)atom.X[0];
				ljatom.y = (real_t)atom.X[1];
				ljatom.z = (real_t)atom.X[2];
			}

			// kk__ является индексом последнего атома, до которого имеет смысл считать взаимодействия.
			// Так как атомы упорядочены по расстояниям, то атомы за kk__ индексом слишком далеки.
			// Параметр настраивается перед счетом в Гоннет варианте через поиск вперед, стартуя от 0.
			unsigned kk__ = atom_count2;

			for (unsigned ii=0; ii<vsz; ii++)
			{
				_LJAtom &atom1 = ljatom1[ii];
				// Счет в обратную сторону от самых дальних атомов, чтобы сохранить точность расчетов.
				for (int ii__=kk__ - 1; ii__>=0; ii__--)
				{
					_LJAtom &atom2 = ljatom2[ii__];

					R[0] = atom2.x - atom1.x;
					R[1] = atom2.y - atom1.y;
					R[2] = atom2.z - atom1.z;
					R[0] -= round(R[0] * _1T[0]) * T[0];
					R[1] -= round(R[1] * _1T[1]) * T[1];
					R[2] -= round(R[2] * _1T[2]) * T[2];

					vbool_t is_interaction = sqr(R[0]) + sqr(R[1]) + sqr(R[2]) < rskin2;
					if (is_interaction == false) continue;
					CONTROL_CONTACT(is_interaction, atom1, atom2);

					for (unsigned i=0; i<vreal_t::size; i++)
					{
						if (is_interaction[i])
						{
							unsigned k = atom1.insert_data[i], k__ = atom2.insert_data[i];
							if (k > k__) std::swap(k, k__);
							verlet_table_[k].push_back(k__);
							pairs_count++;
						}
					}
				}
			}
		}
		return pairs_count;
	}

	TEMPLATE_HEADER
	INLINE _E(real_t) _Verlet_table<TEMPLATE_ARG>
	::dU__dX()
	{
		_E(real_t) energy = 0.; // накопители энергии должны иметь расширенный тип

	#if !defined(SKIP_COUL_ENERGY) || !defined(SKIP_VDW_ENERGY)

		update();

		unsigned atom_count = complex_->count(ATOM);
		unsigned thread_count = parall_.thread_count();

		typedef std::pair<_E(real_t), vector_t *>  wparam_t;
			// энергия и массив сил на запись для каждого потока

		wparam_t wparams[thread_count];

		// подготовим данные для параллельного счета
		for (unsigned i=0; i<thread_count; i++)
		{
			wparams[i].first = 0.;
			wparams[i].second = &verlet_forces_[i][0];
			for (unsigned k=0; k<atom_count; k++) wparams[i].second[k] = 0.;
		}

		unsigned row_count = atom_count;
#ifdef USE_VERLET_TABLE_GROUP
		row_count = verlet_groups_.size();
#endif

		unsigned rparams[row_count]; // идентификаторы рядов Верлет раблицы
		for (unsigned i=0; i<row_count; i++) rparams[i] = i;

		{
			DECLARE_AS_RESTORING_ROUND_MODE
			_MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);

			parall_.start_parallel_running(&_Verlet_table::parall_dU__dX_, (wparam_t *)wparams,
				row_count, (unsigned *)rparams);
			parall_.wait();
		}

		_Atom *atoms = complex_->get(ATOM);
		for (unsigned i=0; i<thread_count; i++)
		{
			energy += wparams[i].first;
			for (unsigned k=0; k<atom_count; k++) atoms[k].F += wparams[i].second[k];
		}
	#endif // !defined(SKIP_COUL_ENERGY) || !defined(SKIP_VDW_ENERGY)

		return energy;
	}

#ifdef USE_VERLET_TABLE_GROUP

	TEMPLATE_HEADER
	inline _E(real_t) _Verlet_table<TEMPLATE_ARG>
	::dU__dX_(vector_t *F__, unsigned igroup)
	{
		const Group_ &group = verlet_groups_[igroup];
		unsigned elem_count = group.elem_count_;

		vreal_t energy = (real_t) 0.f;
		real_t rmax2 = sqr(_Interaction::interaction_radius());

		__LJAtom<vreal_t, vint_t> ljatoms[elem_count];

		_Atom *atoms_ = complex_->get(ATOM);
		for (unsigned i=0; i<elem_count; i++)
		{
			_Atom &atom = atoms_[group.group_elems_[i]];
			ljatoms[i].x = atom.X[0];
			ljatoms[i].y = atom.X[1];
			ljatoms[i].z = atom.X[2];
			ljatoms[i].insert_data = atom.insert_data;
			ljatoms[i].connect_data = atom.connect_data;
			ljatoms[i].charge = atom.charge;
			ljatoms[i].sigma = atom.sigma;
			ljatoms[i].eps = atom.eps;
			ljatoms[i].fx = 0.f;
			ljatoms[i].fy = 0.f;
			ljatoms[i].fz = 0.f;
		}

		// расчет по разделяемым атомам
		const std::vector<unsigned> &row = group.share_atoms_;
		unsigned atom_count = row.size();

		if (atom_count != 0)
		{
			unsigned vsz = velement_count<vreal_t>(atom_count);
				// число элементов в массиве из векторных данных (для SSE оно кратно 4)

			// Загрузка данных из ячеек региона. Создание "псевдоатомов" в начале массива.
			__LJAtom<vreal_t, vint_t> ljatoms__[vsz];
			read_data(&ljatoms__[0], row);

			for (unsigned i=0; i<elem_count; i++)
			{
				__LJAtom<vreal_t, vint_t> &ljatom = ljatoms[i];
				for (unsigned ii=0; ii<vsz; ii++)
				{
					__LJAtom<vreal_t, vint_t> &ljatom__ = ljatoms__[ii];

					LJ_CALCULATE(is_interaction, energy__, fx, fy, fz, rmax2, ljatom, ljatom__);
					ljatom.fx += fx; ljatom.fy += fy; ljatom.fz += fz;
					ljatom__.fx -= fx; ljatom__.fy -= fy; ljatom__.fz -= fz;
					energy += energy__;
				}
			}

			save_data(F__, &ljatoms__[0], row);
		}

		// расчет по уникальным атомам
		for (unsigned i=0; i<elem_count; i++)
		{
			__LJAtom<vreal_t, vint_t> &ljatom = ljatoms[i];

			const std::vector<unsigned> &row = group.uniq_atoms_[i];
			unsigned atom_count = row.size();
			if (atom_count != 0)
			{
				unsigned vsz = velement_count<vreal_t>(atom_count);
					// число элементов в массиве из векторных данных (для SSE оно кратно 4)

				// Загрузка данных из ячеек региона. Создание "псевдоатомов" в начале массива.
				__LJAtom<vreal_t, vint_t> ljatoms__[vsz];
				read_data(&ljatoms__[0], row);

				for (unsigned ii=0; ii<vsz; ii++)
				{
					__LJAtom<vreal_t, vint_t> &ljatom__ = ljatoms__[ii];
					LJ_CALCULATE(is_interaction, energy__, fx, fy, fz, rmax2, ljatom, ljatom__);
					ljatom.fx += fx; ljatom.fy += fy; ljatom.fz += fz;
					ljatom__.fx -= fx; ljatom__.fy -= fy; ljatom__.fz -= fz;
					energy += energy__;
				}

				save_data(F__, &ljatoms__[0], row);
			}
		}

		for (unsigned i=0; i<elem_count; i++)
		{
			__LJAtom<vreal_t, vint_t> &ljatom = ljatoms[i];
			unsigned pos = group.group_elems_[i];
			_Atom &atom = atoms_[pos];
			atom.F[0] += summarize(ljatom.fx);
			atom.F[1] += summarize(ljatom.fy);
			atom.F[2] += summarize(ljatom.fz);
		}

		return summarize(energy);
	}

#else

	TEMPLATE_HEADER
	INLINE _E(real_t) _Verlet_table<TEMPLATE_ARG>
	::dU__dX_(vector_t *F__, unsigned irow)
	{
		std::vector<unsigned> &row = verlet_table_[irow];
		unsigned atom_count = row.size();
		if (atom_count == 0) return 0.;

		vreal_t energy = (real_t) 0.f;
		real_t rmax2 = sqr(_Interaction::interaction_radius());

		__LJAtom<vreal_t, vint_t> ljatom;

		_Atom *atoms_ = complex_->get(ATOM);
		_Atom &atom = atoms_[irow];
		ljatom.x = atom.X[0];
		ljatom.y = atom.X[1];
		ljatom.z = atom.X[2];
		ljatom.insert_data = atom.insert_data;
		ljatom.connect_data = atom.connect_data;

		ljatom.charge = atom.charge;
		ljatom.sigma = atom.sigma;
		ljatom.eps = atom.eps;
		ljatom.fx = 0.f;
		ljatom.fy = 0.f;
		ljatom.fz = 0.f;

		unsigned vsz = velement_count<vreal_t>(atom_count);
			// число элементов в массиве из векторных данных (для SSE оно кратно 4)

		// Загрузка данных из ячеек региона. Создание "псевдоатомов" в начале массива.
		__LJAtom<vreal_t, vint_t> ljatoms[vsz];
		read_data(&ljatoms[0], row);

		for (unsigned ii=0; ii<vsz; ii++)
		{
			__LJAtom<vreal_t, vint_t> &ljatom__ = ljatoms[ii];
			LJ_CALCULATE(is_interaction, energy__, fx, fy, fz, rmax2, ljatom, ljatom__);

			ljatom.fx += fx; ljatom.fy += fy; ljatom.fz += fz;
			ljatom__.fx -= fx; ljatom__.fy -= fy; ljatom__.fz -= fz;
			energy += energy__;
		}

		save_data(F__, &ljatoms[0], row);

		atom.F[0] += summarize(ljatom.fx);
		atom.F[1] += summarize(ljatom.fy);
		atom.F[2] += summarize(ljatom.fz);

		return summarize(energy);
	}
#endif

	TEMPLATE_HEADER
	template <typename T1, typename T2>
	INLINE void _Verlet_table<TEMPLATE_ARG>
	::read_data(__LJAtom<T1, T2> *ljatoms, const std::vector<unsigned> &row)
	{

		typedef typename T1::value_type real_t;
		typedef typename T2::value_type int_t;

		unsigned atom_count = row.size();

		unsigned n = T1::size; // число компонент в sse векторе
		unsigned vsz = velement_count<T1>(atom_count);
			// число элементов в массиве из векторных данных (для SSE оно кратно 4)

		unsigned sz = vsz * n;
		unsigned empty = sz - atom_count;

		const _Atom *atom__ = complex_->get(ATOM);
		for (unsigned i=empty; i<sz; i++)
		{
			__LJAtom<T1, T2> &ljatom = ljatoms[i / n];
			unsigned m = i % n;
			const _Atom *atom = atom__ + row[i-empty];
			*((int_t *)&ljatom.connect_data + m) = (int_t)atom->connect_data;
			*((int_t *)&ljatom.insert_data + m) = (int_t)atom->insert_data;
			*((real_t *)&ljatom.charge + m) = (real_t)atom->charge;
			*((real_t *)&ljatom.sigma + m) = (real_t)atom->sigma;
			*((real_t *)&ljatom.eps + m) = (real_t)atom->eps;
			*((real_t *)&ljatom.x + m) = (real_t)atom->X[0];
			*((real_t *)&ljatom.y + m) = (real_t)atom->X[1];
			*((real_t *)&ljatom.z + m) = (real_t)atom->X[2];
			*((real_t *)&ljatom.fx + m) = 0.f;
			*((real_t *)&ljatom.fy + m) = 0.f;
			*((real_t *)&ljatom.fz + m) = 0.f;
		}

		for (unsigned i=0; i<empty; i++)
		{
			__LJAtom<T1, T2> &ljatom = ljatoms[i / n];
			unsigned m = i % n;
			*((int_t *)&ljatom.connect_data + m) = 0;
			*((int_t *)&ljatom.insert_data + m) = 0;
			*((real_t *)&ljatom.charge + m) = 0.f;
			*((real_t *)&ljatom.sigma + m) = 1.f;
			*((real_t *)&ljatom.eps + m) = 0.f;
			*((real_t *)&ljatom.x + m) = *((real_t *)&ljatom.x + empty);
			*((real_t *)&ljatom.y + m) = *((real_t *)&ljatom.y + empty);
			*((real_t *)&ljatom.z + m) = *((real_t *)&ljatom.z + empty);
			*((real_t *)&ljatom.fx + m) = 0.f;
			*((real_t *)&ljatom.fy + m) = 0.f;
			*((real_t *)&ljatom.fz + m) = 0.f;
		}
	}

	TEMPLATE_HEADER
	template <typename T1, typename T2>
	INLINE void _Verlet_table<TEMPLATE_ARG>
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
	INLINE void _Verlet_table<TEMPLATE_ARG>
	::save_data(vector_t *F, const __LJAtom<T1, T2> *ljatoms, const std::vector<unsigned> &row)
	{
		typedef typename T1::value_type real_t;
		typedef typename T2::value_type int_t;

		unsigned atom_count = row.size();

		unsigned n = T1::size; // число компонент в sse векторе
		unsigned vsz = velement_count<T1>(atom_count);
			// число элементов в массиве из векторных данных (для SSE оно кратно 4)

		unsigned sz = vsz * n;
		unsigned empty = sz - atom_count;

		for (unsigned i=empty; i<sz; i++)
		{
			const __LJAtom<T1, T2> &ljatom = ljatoms[i / n];
			unsigned m = i % n;
			real_t *f = (real_t *)&(F + row[i-empty])[0];
			f[0] += *((const real_t *)&ljatom.fx + m);
			f[1] += *((const real_t *)&ljatom.fy + m);
			f[2] += *((const real_t *)&ljatom.fz + m);
		}
	}

#ifdef USE_VERLET_TABLE_GROUP

	TEMPLATE_HEADER
	INLINE void _Verlet_table<TEMPLATE_ARG>
	::compress_(bool print)
	{
		//function_timer_t<BUILD2_TIMER_> tm(0);

		vector_t T = region_->translation_vector();
		real_t len = pow(global_compress_factor / density_, 1./3.);

		_Index grid_sz
		(
			(unsigned)ceil(T[0] / len),
			(unsigned)ceil(T[1] / len),
			(unsigned)ceil(T[2] / len)
		);
		_Region grid(grid_sz, len);

		unsigned n = molkern::count(NEIGHBOR, rskin_, density_);
			// число соседей каждого атома

		unsigned group_size = grid.count(NODE);
		verlet_groups_.resize(group_size);

		for (unsigned i=0; i<group_size; i++)
		{
			Group_ &group = verlet_groups_[i];
			group.elem_count_ = 0;
			group.share_atoms_.reserve(n);
			for (unsigned k=0; k<MAX_GROUP_ELEM; k++) group.uniq_atoms_[k].reserve(n);
		}

		unsigned thread_count = parall_.thread_count();
		unsigned atom_count = complex_->count(ATOM);
		const _Atom *atoms_ = complex_->get(ATOM);

		share_.resize(thread_count);
		for (unsigned i=0; i<thread_count; i++) share_[i].resize(atom_count, 0);

		typename _LPComplex::ljatom_type ljatom;
		for (unsigned i=0; i<atom_count; i++)
		{
			ljatom.make(atoms_[i]);
			typename _Region::region_vector_type X = region_->make_insert_image(atoms_[i].X);
				// ограничиваемся областью основного региона, а не нового грида, который шире

			grid.insert(ljatom, X);
		}

		unsigned *wparams[thread_count];
		for (unsigned i=0; i<thread_count; i++) wparams[i] = &share_[i][0];

		unsigned node_count = grid.count(NODE);
		unsigned cells[node_count];
		for (unsigned i=0; i<node_count; i++) cells[i] = i;

		parall_.start_parallel_running(&_Verlet_table::parall_compress, (unsigned **)&wparams[0],
			node_count, (unsigned *)&cells[0], (void*)&grid);
		parall_.wait();

		//function_timer_t<BUILD2_TIMER_>::sinc();

		if (print)
		{
			unsigned c1 = 0, c2 = 0, full_group = 0;
			for (unsigned igroup=0; igroup<verlet_groups_.size(); igroup++)
			{
				const Group_ &group = verlet_groups_[igroup];
				unsigned elem_count = group.elem_count_;
				if (elem_count == 0) continue;

				full_group++;
				unsigned share_atom_count = group.share_atoms_.size();
				c1 += share_atom_count;
				c2 += share_atom_count * elem_count;

				for (unsigned m=0; m<elem_count; m++)
				{
					unsigned uniq_count = group.uniq_atoms_[m].size();
					c1 += uniq_count;
					c2 += uniq_count;
				}
			}
				_S msg = _S(" number of atoms in group is ") + make_string("%f", (float)atom_count / full_group)
					+ _S("  gain is ") + make_string("%f", (float)(c2 - c1) / c2);
				PRINT_MESSAGE(msg);
		}
	}

	TEMPLATE_HEADER
	inline void _Verlet_table<TEMPLATE_ARG>
	::parall_compress(void *wparam, void *rparam, void *params)
	{
		unsigned *share = *(unsigned **)wparam;
		unsigned igroup = *(unsigned *)rparam;
		const _Region *grid = (const _Region *)params;

		compress_(share, grid, igroup);
	}

	TEMPLATE_HEADER
	inline void _Verlet_table<TEMPLATE_ARG>
	::compress_(unsigned *share, const _Region *grid, unsigned igroup)
	{
		const _Node *node = grid->get(NODE);

		Group_ &group = verlet_groups_[igroup];
		group.elem_count_ = node[igroup].count(ATOM);

		if (group.elem_count_ == 0) return;
		if (group.elem_count_ > MAX_GROUP_ELEM)
		{
			std::cout << make_string("ERROR: number of group elements (%d) > MAX_GROUP_ELEM !!!",
				group.elem_count_) << std::endl;
			exit(0);
		}

		group.share_atoms_.resize(0);
		for (unsigned m=0; m<MAX_GROUP_ELEM; m++) group.uniq_atoms_[m].resize(0);

		typedef typename _Node::const_iterator _Iterator;
		_Iterator it = node[igroup].begin(), ite = node[igroup].end();
		for (unsigned m=0; it!=ite; ++it, m++)
		{
			const typename _LPComplex::ljatom_type &ljatom = *it;
			group.group_elems_[m] = ljatom.insert_data;
		}

		for (unsigned m=0; m<group.elem_count_; m++)
		{
			unsigned nverl = group.group_elems_[m];
			const std::vector<unsigned> &row = verlet_table_[nverl];
			unsigned sz = row.size();
			for (unsigned j=0; j<sz; j++)
			{
				unsigned l = row[j];
				share[l]++;
			}
		}

		// найдем разделяемые атомы группы
		{
			unsigned nverl = group.group_elems_[0];
			const std::vector<unsigned> &row = verlet_table_[nverl];
			unsigned sz = row.size();
			for (unsigned k=0; k<sz; k++)
			{
				unsigned elem = row[k];
				if (share[elem] == group.elem_count_)
				{
					group.share_atoms_.push_back(elem);
				}
			}
		}

		// уникальные элементы для каждого элемента группы
		for (unsigned k=0; k<group.elem_count_; k++)
		{
			unsigned n = group.group_elems_[k];

			const std::vector<unsigned> &row = verlet_table_[n];
			unsigned row_sz = row.size();
			unsigned cnt = row_sz - group.share_atoms_.size();
			group.uniq_atoms_[k].resize(cnt);

			if (cnt)
			{
				for (unsigned l=0, m=0; l<row_sz; l++)
				{
					unsigned elem = row[l];
					if (share[elem] != group.elem_count_)
						group.uniq_atoms_[k][m++] = elem;
				}
			}
		}

		for (unsigned k=0; k<group.elem_count_; k++)
		{
			const std::vector<unsigned> &row = group.uniq_atoms_[k];
			unsigned sz__ = row.size();
			for (unsigned m=0; m<sz__; m++) share[row[m]] = 0;
		}

		{
			const std::vector<unsigned> &row = group.share_atoms_;
			unsigned sz__ = row.size();
			for (unsigned m=0; m<sz__; m++) share[row[m]] = 0;
		}
	}

#endif

#undef TEMPLATE_HEADER
#undef TEMPLATE_ARG

}
#endif
