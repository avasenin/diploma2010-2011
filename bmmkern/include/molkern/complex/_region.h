#ifndef _REGION__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H
#define _REGION__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H

/*==============================================================================
 *                 ОПИСАНИЕ ОСОБЕННОСТЕЙ РЕАЛИЗАЦИИ ПРОЕКТА
 *
 *  Файл содержит описание объекта, который является  некоторой (прямоугольной)
 *  областью пространства. В отличие от ящика, регион управляет координатами
 *  объектов, поскольку знает свою форму и свойства циклической трансляции.
 *  Поскольку регион может быть определен не только в пространстве, но и на
 *  плоскости (задача морфинга), то у него появляется зависимость от размернсти
 *  пространства.
 =============================================================================*/

#include "molkern/__moldefs.h"
#include "molkern/complex/_parallel.h"

namespace molkern
{
	using namespace prgkern;

	template <typename T1, typename T2>
	struct __LJAtom
	{
		// connect_data : каждый N-бит числа определяет, есть ли связь (1-1 .. 1-4 )
		//   с атомом, смещенным вперед (!) на N-позиций от заданного атома
		// insert_data :
		//   0-30 биты (ref-bits) определяют номер атома в массиве атомов, они
		//   необходимы для пересылки результирующей силы в основное хранилище.
		//   31-бит (fix-bit) равняется 1 при фиксации атома в пространстве
		//   Такой порядок битов позволяет исключить преобразования поля
		//   при выключенной опции USE_FIXED_ATOMS.
		T2 connect_data; // таблица коннектов атома
		T2 insert_data; // тип вставки атома

		// Основные поля, необходимые для вычислений.
		T1 x, y, z; // текущие координаты атома
		T1 charge; // текущий заряд атома в единицах [a.u.q] * sqrt(ELECTRIC_FACTOR)
		T1 sigma; // параметр LJ 6-12 + 1-4 взаимодействия
		T1 eps; // параметр LJ 6-12 + 1-4 взаимодействия
		T1 fx, fy, fz; // текущая сила (coul + vdw), действующая на атом
		T1 hash; // хэш координат атома вдоль некоторого направления
	};

	template <typename T1, typename T2, typename S>
	INLINE void read_empty(__LJAtom<T1, T2> &ljatom, const S *atom) {}

	template <typename S>
	INLINE void read_empty(__LJAtom<vecreal_<4, real_t>, vecint_ <4, int> > &ljatom, const S *atom)
	{
		ljatom.connect_data = 0;
		ljatom.insert_data = 0;
		ljatom.charge = 0.f;
		ljatom.sigma = 1.f;
		ljatom.eps = 0.f;
		ljatom.x = *(const real_t *)&atom->X[0] + vreal_t(0.1f, 0.0f, 0.0f, 0.1f);
		ljatom.y = *(const real_t *)&atom->X[1] + vreal_t(0.0f, 0.1f, 0.0f, 0.1f);
		ljatom.z = *(const real_t *)&atom->X[2] + vreal_t(0.0f, 0.0f, 0.1f, 0.1f);
		ljatom.fx = 0.f;
		ljatom.fy = 0.f;
		ljatom.fz = 0.f;
	}

	template <typename T1, typename T2, typename S>
	INLINE void read_data(__LJAtom<T1, T2> &ljatom, const S *atom)
	{
		ljatom.connect_data = (int)atom->connect_data;
		ljatom.insert_data = (int)atom->insert_data;
		ljatom.charge = (real_t)atom->charge;
		ljatom.sigma = (real_t)atom->sigma;
		ljatom.eps = (real_t)atom->eps;
		ljatom.x = (real_t)atom->X[0];
		ljatom.y = (real_t)atom->X[1];
		ljatom.z = (real_t)atom->X[2];
		ljatom.fx = 0.f;
		ljatom.fy = 0.f;
		ljatom.fz = 0.f;
	}

	template <unsigned N, typename S>
	INLINE void read_data(unsigned m, __LJAtom<vecreal_<N, real_t>, vecint_ <N, int> > &ljatom,
		const S *atom)
	{
		assert(_LT(m, N));
		*((int *)&ljatom.connect_data + m) = (int)atom->connect_data;
		*((int *)&ljatom.insert_data + m) = (int)atom->insert_data;
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

#ifdef USE_VECTORIZATION
	INLINE void ror1(__LJAtom<vecreal_<4, real_t>, vecint_ <4, int> > &atom)
	{
		atom.connect_data = ror1(atom.connect_data);
		atom.insert_data = ror1(atom.insert_data);
		atom.x = ror1(atom.x);
		atom.y = ror1(atom.y);
		atom.z = ror1(atom.z);
		atom.charge = ror1(atom.charge);
		atom.sigma = ror1(atom.sigma);
		atom.eps = ror1(atom.eps);
		atom.fx = ror1(atom.fx);
		atom.fy = ror1(atom.fy);
		atom.fz = ror1(atom.fz);
	}
#else
	INLINE void ror1(__LJAtom<vecreal_<1, real_t>, vecint_ <1, int> > &atom)
	{
	}
#endif

	INLINE void ror1(__LJAtom<real_t, int> &atom)
	{
	}


	/**
	*  парсинг входящей строки, содержащей информацию о ящике
	*/
	inline index_<3, unsigned> parse_box_string(const std::string &s)
	{
	#define MAKE_CONTROL  \
		if (pos == (int)std::string::npos) { \
			std::string msg = _S("incorrect box format: ") + s; \
			PRINT_ERR(msg); \
		}

		index_<3, unsigned> sz;
		int pos = s.find_first_of("[");
		MAKE_CONTROL
		sz[0] = (unsigned)atoi(&s[pos + 1]);

		pos = s.find_first_of(",", pos + 1);
		MAKE_CONTROL
		sz[1] = (unsigned)atoi(&s[pos + 1]);

		pos = s.find_first_of(",", pos + 1);
		MAKE_CONTROL
		sz[2] = (unsigned)atoi(&s[pos + 1]);

	#undef MAKE_CONTROL
		return sz;
	}

	/**
	 * Восстановление оследовательного номера ячейки региона по ее 3d-индексу
	 * @param n 3d индекс ячейки
	 * @param s число ячеек по всем направлениям региона
	 * @return позиция ячейки в последовательном массиве ячеек
	 */
	INLINE unsigned make_pos(index_<3, unsigned> n, index_<3, unsigned> s)
	{
		// нормализуем индекс, чтобы не было компонент за пределами региона
		n = index_<3, unsigned>(mod(n[0], s[0]), mod(n[1], s[1]), mod(n[2], s[2]));

		unsigned pos = (unsigned)((n[0] * s[1] + n[1]) * s[2] + n[2]);
		assert(_LT(pos, (unsigned)(s[0] * s[1] * s[2])));

		return pos;
	}

	INLINE unsigned make_pos(index_<3, int> n, index_<3, unsigned> s)
	{
		// нормализуем индекс, чтобы не было компонент за пределами региона
		n = index_<3, int>(mod(n[0], (int)s[0]), mod(n[1], (int)s[1]), mod(n[2], (int)s[2]));

		unsigned pos = (unsigned)((n[0] * s[1] + n[1]) * s[2] + n[2]);
		assert(_LT(pos, (unsigned)(s[0] * s[1] * s[2])));

		return pos;
	}

	INLINE index_<3, unsigned> make_index(unsigned pos, index_<3, unsigned> s)
	{
		assert(_LT(pos, s[0] * s[1] * s[2]));

		unsigned p__ = pos / s[2];
		unsigned kx = p__ / s[1];
		unsigned ky = p__ - kx * s[1];
		unsigned kz = pos - p__* s[2];
		assert(_LT(kx, s[0])); assert(_GE(kx, (unsigned)0));
		assert(_LT(ky, s[1])); assert(_GE(ky, (unsigned)0));
		assert(_LT(kz, s[2])); assert(_GE(kz, (unsigned)0));

		return index_<3, unsigned>(kx, ky, kz);
	}

	#ifdef USE_GONNET
		const vecreal_<4, real_t> gonnet_coefs[14] =
		{
		#define _C(v) vreal_t::value_type(v)
			vecreal_<4, real_t>(_C(M_1_SQRT_3),  _C(M_1_SQRT_3),  _C(M_1_SQRT_3), _C(0.f)), /* 0 */
			vecreal_<4, real_t>(_C(M_1_SQRT_2),  _C(M_1_SQRT_2),     _C(0.f),     _C(0.f)), /* 1 */
			vecreal_<4, real_t>(_C(M_1_SQRT_3),  _C(M_1_SQRT_3), _C(-M_1_SQRT_3), _C(0.f)), /* 2 */
			vecreal_<4, real_t>(_C(M_1_SQRT_2),     _C(0.f),      _C(M_1_SQRT_2), _C(0.f)), /* 3 */
			vecreal_<4, real_t>(     _C(1.f),       _C(0.f),         _C(0.f),     _C(0.f)), /* 4 */
			vecreal_<4, real_t>(_C(M_1_SQRT_2),     _C(0.f),     _C(-M_1_SQRT_2), _C(0.f)), /* 5 */
			vecreal_<4, real_t>(_C(M_1_SQRT_3), _C(-M_1_SQRT_3),  _C(M_1_SQRT_3), _C(0.f)), /* 6 */
			vecreal_<4, real_t>(_C(M_1_SQRT_2), _C(-M_1_SQRT_2),     _C(0.f),     _C(0.f)), /* 7 */
			vecreal_<4, real_t>(_C(M_1_SQRT_3), _C(-M_1_SQRT_3), _C(-M_1_SQRT_3), _C(0.f)), /* 8 */
			vecreal_<4, real_t>(     _C(0.f),    _C(M_1_SQRT_2),  _C(M_1_SQRT_2), _C(0.f)), /* 9 */
			vecreal_<4, real_t>(     _C(0.f),       _C(1.f),         _C(0.f),     _C(0.f)), /* 10 */
			vecreal_<4, real_t>(     _C(0.f),    _C(M_1_SQRT_2), _C(-M_1_SQRT_2), _C(0.f)), /* 11 */
			vecreal_<4, real_t>(     _C(0.f),       _C(0.f),         _C(1.f)    , _C(0.f)), /* 12 */
			vecreal_<4, real_t>(     _C(0.f),       _C(0.f),         _C(0.f)    , _C(0.f))  /* 13 */
		#undef _C
		};
	#endif

	template <unsigned N, typename S> class Node_;

	template <> class Node_<3, void>
	{
		typedef index_<3, int>       _Key;
		typedef index_<3, unsigned>  _Index;

		_Index ndx_; // 3D индекс ячейки
		int neigbor_cells_[27]; // список номеров соседних ячеек (+ самовключение)

	public:

		typedef void value_type;

		Node_() {}

		/// Инициализация ячейки
		void init(_B2T<YES_PERIODIC_>, unsigned pos, const _Index &s)
		{
			_Index ndx = make_index(pos, s);
			neigbor_cells_[ 0] = make_pos(_Key(ndx[0] - 1, ndx[1] - 1, ndx[2] - 1), s);
			neigbor_cells_[ 1] = make_pos(_Key(ndx[0] - 1, ndx[1] - 1, ndx[2]    ), s);
			neigbor_cells_[ 2] = make_pos(_Key(ndx[0] - 1, ndx[1] - 1, ndx[2] + 1), s);
			neigbor_cells_[ 3] = make_pos(_Key(ndx[0] - 1, ndx[1]    , ndx[2] - 1), s);
			neigbor_cells_[ 4] = make_pos(_Key(ndx[0] - 1, ndx[1]    , ndx[2]    ), s);
			neigbor_cells_[ 5] = make_pos(_Key(ndx[0] - 1, ndx[1]    , ndx[2] + 1), s);
			neigbor_cells_[ 6] = make_pos(_Key(ndx[0] - 1, ndx[1] + 1, ndx[2] - 1), s);
			neigbor_cells_[ 7] = make_pos(_Key(ndx[0] - 1, ndx[1] + 1, ndx[2]    ), s);
			neigbor_cells_[ 8] = make_pos(_Key(ndx[0] - 1, ndx[1] + 1, ndx[2] + 1), s);
			neigbor_cells_[ 9] = make_pos(_Key(ndx[0]    , ndx[1] - 1, ndx[2] - 1), s);
			neigbor_cells_[10] = make_pos(_Key(ndx[0]    , ndx[1] - 1, ndx[2]    ), s);
			neigbor_cells_[11] = make_pos(_Key(ndx[0]    , ndx[1] - 1, ndx[2] + 1), s);
			neigbor_cells_[12] = make_pos(_Key(ndx[0]    , ndx[1]    , ndx[2] - 1), s);
			neigbor_cells_[13] = pos; // самовключение
			neigbor_cells_[14] = make_pos(_Key(ndx[0]    , ndx[1]    , ndx[2] + 1), s);
			neigbor_cells_[15] = make_pos(_Key(ndx[0]    , ndx[1] + 1, ndx[2] - 1), s);
			neigbor_cells_[16] = make_pos(_Key(ndx[0]    , ndx[1] + 1, ndx[2]    ), s);
			neigbor_cells_[17] = make_pos(_Key(ndx[0]    , ndx[1] + 1, ndx[2] + 1), s);
			neigbor_cells_[18] = make_pos(_Key(ndx[0] + 1, ndx[1] - 1, ndx[2] - 1), s);
			neigbor_cells_[19] = make_pos(_Key(ndx[0] + 1, ndx[1] - 1, ndx[2]    ), s);
			neigbor_cells_[20] = make_pos(_Key(ndx[0] + 1, ndx[1] - 1, ndx[2] + 1), s);
			neigbor_cells_[21] = make_pos(_Key(ndx[0] + 1, ndx[1]    , ndx[2] - 1), s);
			neigbor_cells_[22] = make_pos(_Key(ndx[0] + 1, ndx[1]    , ndx[2]    ), s);
			neigbor_cells_[23] = make_pos(_Key(ndx[0] + 1, ndx[1]    , ndx[2] + 1), s);
			neigbor_cells_[24] = make_pos(_Key(ndx[0] + 1, ndx[1] + 1, ndx[2] - 1), s);
			neigbor_cells_[25] = make_pos(_Key(ndx[0] + 1, ndx[1] + 1, ndx[2]    ), s);
			neigbor_cells_[26] = make_pos(_Key(ndx[0] + 1, ndx[1] + 1, ndx[2] + 1), s);

			ndx_ = ndx;
		}

		/// Инициализация ячейки
		void init(_B2T<NO_PERIODIC_>, unsigned pos, const _Index &s)
		{
			init(YES_PERIODIC, pos, s);
			_Index ndx = make_index(pos, s);

		#define OK(a, c, b) (a (c) && (c) b)

			unsigned k = 0;
			for (int x=-1; x<2; x++)
			for (int y=-1; y<2; y++)
			for (int z=-1; z<2; z++)
			{
				if ( OK(0 <=, ndx[0] + x, < s[0])
					&& OK(0 <=, ndx[1] + y, < s[1])
					&& OK(0 <=, ndx[2] + z, < s[2])
				)
				{
					continue;
				}
				else
				{
					// запишем недействительный индекс для несуществующей ячейки
					// чтобы позволять давать соседей для заданных направлений
					// -1 означает, что для заданного направления нет соседа
					neigbor_cells_[k] = -1;
				}
				k++;
			}

		#undef OK

			ndx_ = ndx;
		}

		/// Получение индекса ячейки
		_Index index() const { return ndx_; }

		/**
		 *  Получение соседа по заданному направлению
		 * @param direction
		 * @return -1, если нет соседа в заданном направлении
		 */
		int get_neighbor(unsigned direction) const
		{
			assert(_LT(direction, (unsigned)27));
			return neigbor_cells_[direction];
		}

		/**
		*  Формирует номера соседних ячеек.
		* @param neigs позиции ячеек, которые являются соседними для заданной
		* @param pos номер ячейки, для которой формируются соседи
		* @return число соседних ячеек
		*/
		unsigned get_neigbours(_B2T<YES_PERIODIC_>, unsigned *neigs) const
		{
			for (unsigned i=0; i<27; i++) neigs[i] = get_neighbor(i);
			return 27;
		}

		unsigned get_neigbours(_B2T< NO_PERIODIC_>, unsigned *neigs) const
		{
			unsigned cnt = 0; int nei;
			for (unsigned i=0; i<27; i++)
			{
				if ((nei = get_neighbor(i)) != -1) neigs[cnt++] = nei;
			}
			return cnt;
		}
	};


	//-------------------------------------------------------------------------------------------------
	//                                         Gonnet
	//-------------------------------------------------------------------------------------------------
	struct Gonnet_vector : public std::vector<std::pair<unsigned, real_t>  > // <индекс, хеш>
	{
		typedef std::pair<unsigned, real_t>  _Pair;

		/**
		* Формирует хэш-код (по гоннету) для заданного направления
		* @param direction направление в пространстве
		* @param x {x,y,z,0} вектор
		* @return хэш код вектора x
		*/
		static real_t make_hash(unsigned direction, vecreal_<4, real_t> x)
		{ return summarize(gonnet_coefs[direction] * x); }

		/// делает код направления по вектору направления
		static unsigned make_direction_code(const vector_t &d)
		{
			unsigned code = 13 + (sign(d[0]) * 9 + sign(d[1]) * 3 + sign(d[2]));
			if (code > 13) code = 26 - code;
			assert(_LT(code, (unsigned)13));
			return code;
		}

		/// рассчитывает смещения кодов Гоннета по вектору направления
		static real_t make_direction_shift(const vector_t &d)
		{
			return sqrt(sqr(d[0]) + sqr(d[1]) + sqr(d[2]));
		}

		/**
		 * Вставка хэш кода в массив. При вставке дополнительно вставляется номер по порядку,
		 * с целью дальнейшей сортировки кодов и возможности ссылок на номера вставляемых атомов.
		 * @param hash
		 */
		void insert(real_t hash) { push_back(_Pair(size(), hash)); }
		void make_ordering() { std::sort(begin(), end(), less_<_Pair, real_t, &_Pair::second>());	}
	};

	//-------------------------------------------------------------------------------------------------
	//                                         Node<3,S>
	//-------------------------------------------------------------------------------------------------
	template <typename S> class Node_<3, S> : public Node_<3, void>
	{
		typedef Node_<3, void>                _Base;
		typedef _Base::_Index                 _Index;

		std::vector<S> items_[8]; // массив элементов ячейки
		std::vector<vector_t> X_[8]; // координаты вставки элементов

	#ifdef USE_GONNET
		Gonnet_vector gv_[14];
	#endif

	public:

		typedef S value_type;
		typedef typename std::vector<S>::const_iterator  const_iterator;
		typedef typename std::vector<S>::iterator        iterator;

		using _Base::init;
		using _Base::index;
		using _Base::get_neighbor;

		/// Очистка ячейки от атомов.
		void reserve(unsigned n)
		{
			items_[0].reserve(n);
			X_[0].reserve(n);

		#ifdef USE_GONNET
			for (unsigned i=0; i<14; i++) gv_[i].reserve(n);
		#endif
		}

		/// Очистка ячейки от атомов.
		void clear()
		{
			for (unsigned i=0; i<8; i++)
			{
				items_[i].resize(0);
				X_[i].resize(0);
			}

		#ifdef USE_GONNET
			for (unsigned i=0; i<14; i++) gv_[i].resize(0);
		#endif
		}

		/// Вставка атома в ячейку
		void insert(const S &s, const vector_t &X, unsigned i=0)
		{
			items_[i].push_back(s);
			X_[i].push_back(X);
		}

		void insert(const S &s, unsigned i=0)
		{
			items_[i].push_back(s);
		}

		/// Получение числа элементов
		unsigned count(_I2T<ATOM_>, unsigned i=0) const { return items_[i].size(); }

	#ifdef USE_GONNET

		void make_gonnet_ordering()
		{
			unsigned count = items_[0].size();
			for (unsigned i=0; i<13; i++)
			{
				gv_[i].reserve(count);

				for (unsigned k=0; k<count; k++)
				{
					const vector_t &X = X_[0][k];
					real_t code = Gonnet_vector::make_hash(i, vecreal_<4, real_t>(X[0], X[1], X[2], real_t(0.)));
					gv_[i].insert(code);
				}
				// упорядочение индексов по гоннет-кодам
				gv_[i].make_ordering();
			}
			for (unsigned k=0; k<count; k++) gv_[13].insert(0.);
		}


		const std::pair<unsigned, real_t> *get_gonnet_pair(unsigned direction) const
		{
			assert(_LT(direction, (unsigned)13));
			return &gv_[direction][0];
		}

#endif

		/// Получение начала массива данных
		const S *get(_I2T<ATOM_>, unsigned i=0) const { return items_[i].size() ? &items_[i][0] : 0; }
		S *get(_I2T<ATOM_>, unsigned i=0) { return items_[i].size() ? &items_[i][0] : 0; }

		///  Возвращает итераторы по элементам ячейки.
		const_iterator begin(unsigned i=0) const { return items_[i].begin(); }
		const_iterator end(unsigned i=0) const { return items_[i].end(); }
		iterator begin(unsigned i=0) { return items_[i].begin(); }
		iterator end(unsigned i=0) { return items_[i].end(); }
	};


	template <unsigned N, typename S=unsigned> class Region_;

	/**
	*  Класс подобласти (грида), в котором можно построить решетку с/без
	*  циклических свойств. Грид выдает набор пар, попавших в область
	*  взаимодействия, отфильтровывая пары, оказавшиеся на большем расстоянии.
	*  При режиме вставки YES_PERIODIC данный грид может вставлять объекты только
	*  такого класса, которые содержат поле X. Это ограничение введено,
	*  поскольку для эффективности счета нужно приводить координаты к координатам,
	*  находящимся внутри грида. Если этого не делать, то для проверки расстояний
	*  требуется проверять расстояния x-x', x-x'-T, x-x'-2T, x-x'+T, x-x'+2T.
	*  Если не делать приведения координты, то обязательность наличия поля X
	*  исчезает.
	*
	*  Ящик всегда начинается в точке (0, 0, 0) и имеет положительную верхнюю
	*  точку (по построению и избегаются некоторые проверки для скорости)
	*
	*  Чтобы избежать ряда проблем (потери соседей на границах, неравномерной
	*  загрузки узлов кластера при разных размерах ячеек) мы принудительно
	*  масштабируем ящик под целое кратное радиуса взаимодействия.
	*/
	template <typename S> class Region_<3, S>
	{
	public:

		enum {
			dim = 3, // размерность решетки
			neighbours = 27 // число всех соседний ячеек (+ самовключение)
		};

	protected:

		typedef Box_<3, real_t>      _Box;
		typedef Node_<3, S>          _Node;
		typedef index_<3, int>       _Key;
		typedef index_<3, unsigned>  _Index;

		// совокупности независимых пар ячеек для каждого направления
		std::vector<std::vector<unsigned> > pairs_sets_;
		std::vector<unsigned> pairs_sets_direction_;

		std::vector<_Node> nodes_;

		_Index sz_; // число ячеек по каждой координате
		real_t h_, _1h_;  // прямой и обратный шаг ячейки = 1 / h
		vector_t T_, _1T_; // прямой и обратный вектор трансляции грида
		vecreal_<4, real_t> t_, _1t_;
		bool flag_;

		struct region_vector_t : public vecreal_<4, real_t>
		{
			region_vector_t(real_t x0=0.f, real_t x1=0.f, real_t x2=0.f, real_t x3=0.f)
			: vecreal_<4, real_t>(x0, x1, x2, x3) {}

			region_vector_t(const vecreal_<4, real_t> &v) : vecreal_<4, real_t>(v) {}

			region_vector_t(const vector_t &x) : vecreal_<4, real_t>(x[0], x[1], x[2], (real_t)0.f) {}

			operator vecreal_<4, real_t>() { return *this; }
		};

		struct region_key_t : public vecint_ <4, int>
		{
			region_key_t(int x0=0, int x1=0, int x2=0, int x3=0)
			: vecint_<4, int>(x0, x1, x2, x3) {}

			region_key_t(const vecint_<4, int> &v) : vecint_<4, int>(v) {}

			operator vecint_<4, int>() { return *this; }
		};

	public:

		typedef _Node                           node_type;
		typedef typename _Node::iterator        iterator;
		typedef typename _Node::const_iterator  const_iterator;
		typedef _Key                            key_type;
		typedef _Index                          index_type;
		typedef S                               value_type;
		typedef region_vector_t                 region_vector_type;
		typedef region_key_t                    region_key_type;

		Region_(const std::string &sbox, real_t interaction_radius, const vector_t *T=0)
		{ resize(parse_box_string(sbox), interaction_radius, T); }

		Region_(const index_type &sz, real_t interaction_radius, const vector_t *T=0)
		{ resize(sz, interaction_radius, T); }

		/**
		*  Функция увеличивает ящик таким образом, чтобы ящик покрыл заданые
		*  новые размеры. При этом новые параметры ящика удовлетворяют программным
		*  требованиям, чтобы эффективно генерились пары соседних ячеек и т.д.
		*  Размер ящика определяется как произведение числа ячеек(sz) на
		*  длину ячейки (interection_radius).
		* @param len новые размеры ящика
		*/
		void resize(index_type sz, real_t interaction_radius, const vector_t *T=0, bool is_print=NO_PRINT);

		/**
		*  Функция увеличивает ящик таким образом, чтобы ящик покрыл заданые
		*  новые размеры. При этом новые параметры ящика удовлетворяют программным
		*  требованиям, чтобы эффективно генерились пары соседних ячеек и т.д.
		* @param len новые размеры ящика
		*/
		void enlarge(const _Box &box)
		{
			vector_t len = box.length() + 0.5 * h_;
				// делаем чуть больше, чтобы при цикличности молекулы разделялись расстоянием,
				// не менее радиуса взаимодействия

			_Index sz
			(
				(unsigned) ceil(len[0] * _1h_),
				(unsigned) ceil(len[1] * _1h_),
				(unsigned) ceil(len[2] * _1h_)
			);
			resize(sz, h_, 0, YES_PRINT);
		}

		/// Очистка области от атомов.
		void clear() { for (unsigned i=0,sz=nodes_.size(); i<sz; i++) nodes_[i].clear(); }

		/// текущий ящик области
		_Box get(_I2T<BOX_>) const { return _Box(vector_t(0.,0.,0.), T_); }

		/// размерности текущего ящика области
		index_type box_size() const { return sz_; }

		bool is_edge(unsigned pos1, unsigned pos2) const
		{
			_Index key = nodes_[pos1].index();
			if ( (key[0] == (int)sz_[0] - 1) || (key[1] == (int)sz_[1] - 1) || (key[2] == (int)sz_[2] - 1)
		  ) return true;

			key = nodes_[pos2].index();
			if ( (key[0] == (int)sz_[0] - 1) || (key[1] == (int)sz_[1] - 1) || (key[2] == (int)sz_[2] - 1)
		  ) return true;

			return false;
		}

		/// трансляционный вектор всего региона
		vector_t translation_vector() const { return vector_t(T_[0], T_[1], T_[2]); }

		/// трансляционный вектор между ячеками pos1 & pos2
		vector_t translation_vector(unsigned pos1, unsigned pos2) const
		{
			_Index key1 = nodes_[pos1].index();
			_Index key2 = nodes_[pos2].index();
			_Key d(key2[0] - key1[0], key2[1] - key1[1], key2[2] - key1[2]);
				// смещения ячеек друг относительно друга в пространстве

			if (abs(d[0]) == (int)sz_[0] - 1) d[0] = -sign(1, d[0]);
			if (abs(d[1]) == (int)sz_[1] - 1) d[1] = -sign(1, d[1]);
			if (abs(d[2]) == (int)sz_[2] - 1) d[2] = -sign(1, d[2]);
			return vector_t(h_ * d[0], h_ * d[1], h_ * d[2]);
		}

		// получение координаты внутри региона
		region_vector_t make_insert_image(const vector_t &X) const
		{
			vecreal_<4, real_t> x(X[0], X[1], X[2], real_t(0.));
//			return region_vector_t(x - ifloor(x * _1t_) * t_);
			region_vector_t y(x - ifloor(x * _1t_) * t_);
			if (y[0] < 0.f) y[0] += T_[0];
			if (y[1] < 0.f) y[1] += T_[1];
			if (y[2] < 0.f) y[2] += T_[2];
			if (y[0] >= T_[0]) y[0] -= T_[0];
			if (y[1] >= T_[1]) y[1] -= T_[1];
			if (y[2] >= T_[2]) y[2] -= T_[2];
			return y;
		}

		vector_t make_insert_image_(const region_vector_t &x)
		{
			vecreal_<4, real_t> x__(x[0], x[1], x[2], real_t(0.));
			x__ -= ifloor(x__ * _1h_) * h_;
			return vector_t(x__[0], x__[1], x__[2]);
		}

		region_key_t make_insert_index3(const region_vector_t &X) const
		{
			return region_key_t(ifloor((vecreal_<4, real_t>&)X * _1h_));
		}

		unsigned make_insert_index(const region_key_t &key) const
		{
			return (key[0] * sz_[1] + key[1]) * sz_[2] + key[2];
		}

		unsigned make_insert_index(const region_vector_t &X) const
		{
			return make_insert_index(make_insert_index3(X));
		}

		/**
		 * Создает координаты вставки внутри циклически транслируемого региона.
		 * Исходная координата может быть вне региона. Функция ищет ее образ внутри региона.
		 * Поиск координаты и позиции вставки делаются синхронно для эффективности.
		 * @param y[out] координаты точки внутри региона (приведенные к ячейке)
		 * @param pos[out] позиция вставки
		 * @param X исходные координаты
		 * @return true, если координата X попадает внутрь региона
		 */
		_Index index3(unsigned pos) const
		{
			unsigned kz = pos % sz_[2];
			unsigned ky = (pos / sz_[2]) % sz_[1];
			unsigned kx = pos / (sz_[2] * sz_[1]);
			return _Index(kx, ky, kz);
		}

		template <typename T>
		void make_nearest_image_vector(T &dx, T &dy, T &dz) const
		{
			dx -= round(dx * _1T_[0]) * T_[0];
			dy -= round(dy * _1T_[1]) * T_[1];
			dz -= round(dz * _1T_[2]) * T_[2];
		}

		/**
		 * Вставка объекта в регион.
		 * @param pos номер ячейки для вставки
		 * @param s объект
		 * @param x приведенная (то есть внутри ячейки pos) координата вставки
		 */
		bool insert(_B2T<YES_PERIODIC_>, const S &s, const vector_t &X)
		{
			insert_(s, make_insert_image(X));
			return true;
		}

		bool insert(const S &s, const region_vector_t &X)
		{
			insert_(s, X);
			return true;
		}

		bool insert(_B2T<NO_PERIODIC_>, const S &s, const vector_t &X)
		{
			if ( !get(BOX).is_included(X) ) return false;
			insert_(s, region_vector_t(X));
			return true;
		}

		template <typename T>
		void make_piecewise_vector(piecewise_vector<T> *atoms, const _Key &key, unsigned pos) const
		{
			const _Node *node = get(NODE, pos);
			_Index ndx = index3(pos);

			if (flag_ && (ndx[0] == sz_[0] - 1 || ndx[1] == sz_[1] - 1 || ndx[2] == sz_[2] - 1))
			{
				unsigned code = 0;
				if (abs(key[0]) == 1) code += 4;
				if (abs(key[1]) == 1) code += 2;
				if (abs(key[2]) == 1) code += 1;

				// подгрузим атомы по коду
			#define ADD_ATOMS(__code) \
				if ((code & __code) == __code) { \
					atoms->insert(node->get(ATOM, __code), node->count(ATOM, __code)); \
				}

				ADD_ATOMS(0x7)
				ADD_ATOMS(0x6)
				ADD_ATOMS(0x5)
				ADD_ATOMS(0x4)
				ADD_ATOMS(0x3)
				ADD_ATOMS(0x2)
				ADD_ATOMS(0x1)

			#undef ADD_ATOMS
			}
		}

		real_t get(_I2T<EDGE_LENGTH_>) const { return h_; }

		/**
		*  Число ячеек в регионе.
		*/
		unsigned count(_I2T<NODE_>) const { return nodes_.size(); }
		const _Node *get(_I2T<NODE_>, unsigned i=0) const
		{
			assert(_LT(i, (unsigned) nodes_.size()));
			return &nodes_[i];
		}
		_Node *get(_I2T<NODE_>, unsigned i=0)
		{
			assert(_LT(i, (unsigned) nodes_.size()));
			return &nodes_[i];
		}

		/// Возвращает количество упорядочиваний ячеек в пары (26).
		unsigned count(_I2T<DIRECTION_>) const { return 13; }

		/// Возвращает количество упорядочиваний ячеек в пары (26).
		unsigned count(_I2T<SET_>) const { return pairs_sets_.size(); }

		/// Возвращает количество ячеек во множестве
		unsigned count(_I2T<SET_>, unsigned i) const { return pairs_sets_[i].size(); }

		/// номер направления, для заданного множества пар
		unsigned get(_I2T<DIRECTION_>, unsigned iset) const { return pairs_sets_direction_[iset]; }

		/// номер направления, когда ячейка указывает на себя (!) Связана со следующей функцией.
		unsigned get(_I2T<EMPTY_DIRECTION_>) const { return 13; }

		/// возращает указатель на тот массив пар ячеек для заданого множества
		const unsigned *get(_I2T<PAIR_>, unsigned nset) const { return &pairs_sets_[nset][0]; }

	protected:

		/// Генерит хранилище пар соседних ячеек
		unsigned generate_cell_pairs_();

		/// вставляет объект (и его образы, в случае T != sz * h) в заданные координаты
		void insert_(const S &s, const region_vector_t &X);

	};

	#define TEMPLATE_HEADER  template <typename S>
	#define TEMPLATE_ARG     3, S

	TEMPLATE_HEADER
	inline void Region_<TEMPLATE_ARG>
	::resize(index_type sz, real_t interaction_radius, const vector_t *T, bool is_print)
	{
		// Установим достаточный размер по всем направлениям
		if (sz[0] < 3) sz[0] = 3;
		if (sz[1] < 3) sz[1] = 3;
		if (sz[2] < 3) sz[2] = 3;

		if (sz_[0] == sz[0] && sz_[1] == sz[1] && sz_[2] == sz[2]) return;

		if (is_print)
		{
			_S msg = make_string("The region has been changed { %d, %d, %d} -> { %d, %d, %d } cells.",
				sz_[0], sz_[1], sz_[2], sz[0], sz[1], sz[2]);
			PRINT_MESSAGE(msg);
		}

		h_ = interaction_radius;
		_1h_ = (real_t) 1. / interaction_radius;

		flag_ = false;
		if (T) { T_ = *T; flag_ = true; }
		else T_ = vector_t
		(
			sz[0] * interaction_radius,
			sz[1] * interaction_radius,
			sz[2] * interaction_radius
		);

		_1T_[0] = (real_t) 1. / T_[0];
		_1T_[1] = (real_t) 1. / T_[1];
		_1T_[2] = (real_t) 1. / T_[2];

		t_ = vecreal_<4, real_t>(T_[0], T_[1], T_[2], real_t(0.));
		_1t_ = vecreal_<4, real_t>(_1T_[0],  _1T_[1],  _1T_[2], real_t(0.));
		sz_ = sz;

		const real_t DEFAULT_DISTANCE_BEETWEEN_ATOMS = (real_t) 2.1; // [Angstrom]
		//------------------------------------------------------------------------
		// резервируем память, чтобы после не перераспределять ее
		// с соответствующими накладными расходами
		//------------------------------------------------------------------------
		unsigned size = sz_[0] * sz_[1] * sz_[2];
		assert(_NE((unsigned)0, size));

		unsigned count = cube((unsigned) ceil(1. / (_1h_ * DEFAULT_DISTANCE_BEETWEEN_ATOMS)));

		nodes_.resize(size);
		for (unsigned pos=0; pos<size; pos++)
		{
			nodes_[pos].init(YES_PERIODIC, pos, sz_);
			nodes_[pos].reserve(count);
		}
		generate_cell_pairs_();
	}

	TEMPLATE_HEADER
	inline unsigned Region_<TEMPLATE_ARG>
	::generate_cell_pairs_()
	{
		unsigned nx = (unsigned)sz_[0];
		unsigned ny = (unsigned)sz_[1];
		unsigned nz = (unsigned)sz_[2];
		unsigned sz = nx * ny * nz;

		//-------------------------------------------------------------------------------
		//          подсчитаем возможное количество множеств независимых пар
		//-------------------------------------------------------------------------------
		unsigned odds = 0; // общее число нечетных размерностей
		odds += nx % 2 ? 1 : 0;
		odds += ny % 2 ? 1 : 0;
		odds += nz % 2 ? 1 : 0;

		unsigned set_count = 26;
		switch (odds)
		{
		case 1: set_count  += 1; break;
		case 2: set_count  += 4; break;
		case 3: set_count += 13; break;
		}
		pairs_sets_.resize(set_count + 1); // добавим множество одиночных ячеек
		pairs_sets_direction_.resize(set_count + 1); // добавим множество одиночных ячеек

		std::vector<bool> is_marked(sz);
		std::vector<unsigned> chains(sz);
		std::vector<unsigned> chain_end;

		unsigned cur_set = 0;
		for (unsigned direction=0; direction<13; direction++)
		{
			is_marked.assign(sz, false);
			chains.resize(0);
			chain_end.resize(0);

			// количество множеств в совокупности множеств пар для направления
			unsigned set_count__ = 2; // по умалчиванию минимально достаточное

			//------------------------------------------------------------------------------
			//                   построим множество направленных цепочек
			//------------------------------------------------------------------------------
			unsigned start_pos = 0; // первая неиспользовнный для построения пар ячейка
			do { // цикл по построению множества цепочек
				unsigned pos = start_pos;
				do { // цикл по построению одной циклической цепочки
					chains.push_back(pos);
					is_marked[pos] = true;
					pos = nodes_[pos].get_neighbor(direction);
				} while (pos != start_pos);
					// пользуемся тем, что BMMKERN работает только с циклическими ячейками

				unsigned count = chains.size();
				chain_end.push_back(count); // сохраним (недействительную) позицию конца цепочки
				if (count % 2) set_count__ = 3; // при нечетном числе элементов нужно 3 множества

				while (is_marked[++start_pos] == true) ;
					// ищем первый еще неиспользованный элемент для построения следующей цепочки

			} while (start_pos < sz);

			//------------------------------------------------------------------------------
			//                построим множество пар из направленных цепочек
			//------------------------------------------------------------------------------
			unsigned chain_count = chain_end.size();
			unsigned istart, iend = 0, iset;
			for (unsigned ichain=0; ichain<chain_count; ichain++)
			{
				istart = iend;
				iend = chain_end[ichain];
				for (unsigned i=istart; i<iend-1; i++)
				{
					iset = cur_set + i % set_count__;
					pairs_sets_[iset].push_back(chains[i + 1]);
					pairs_sets_[iset].push_back(chains[i]);
						// Делаем обмен ячеек, так как половина направлений по которым делется хеширование
						// в Гоннете убрана в силу симметрии. (!) После обмена направление от pos1 к pos2
						// теперь всегда продуцирует направление в диапозоне [0,13)
				}
				// вставка последней пары для циклического случая
				iset = cur_set + ((iend + 1 - istart) % set_count__ ?
					(iend - 1) % set_count__ : iend % set_count__);
				pairs_sets_[iset].push_back(chains[istart]);
				pairs_sets_[iset].push_back(chains[iend - 1]);
			}
			for (unsigned i=0; i<set_count__; i++) pairs_sets_direction_[cur_set + i] = direction;
			cur_set += set_count__;
		}

		// вставка множества одиночных ячеек
		for (unsigned i=0; i<sz; i++) pairs_sets_[cur_set].push_back(i);
		pairs_sets_direction_[cur_set] = 13;

		return set_count + 1;
	}

	TEMPLATE_HEADER
	INLINE void Region_<TEMPLATE_ARG>
	::insert_(const S &s, const region_vector_t &X__)
	{
		region_vector_t X(X__);
		region_key_t key = make_insert_index3(X);

		// вставляем дубли
		if (flag_ && (key[0] == 0 || key[1] == 0 || key[2] == 0))
		{
			unsigned code = 0;
			if (X[0] + T_[0] < sz_[0] * h_) code += 4;
			if (X[1] + T_[1] < sz_[1] * h_) code += 2;
			if (X[2] + T_[2] < sz_[2] * h_) code += 1;

		#define INSERT_MACRO(__code, __t0, __t1, __t2) \
			if ((code & __code) == __code) \
			{ \
				region_vector_t X__(X[0] + __t0, X[1] + __t1, X[2] + __t2); \
				unsigned pos__ = make_insert_index(X__); \
				nodes_[pos__].insert(s, __code); \
			}

			INSERT_MACRO(0x7, T_[0], T_[1], T_[2]);
			INSERT_MACRO(0x6, T_[0], T_[1],     0);
			INSERT_MACRO(0x5, T_[0],     0, T_[2]);
			INSERT_MACRO(0x4, T_[0],     0,     0);
			INSERT_MACRO(0x3,     0, T_[1], T_[2]);
			INSERT_MACRO(0x2,     0, T_[1],     0);
			INSERT_MACRO(0x1,     0,     0, T_[2]);
		#undef INSERT_MACRO
		}

		unsigned pos = make_insert_index(key);
#ifdef USE_GONNET
		nodes_[pos].insert(s, make_insert_image_(X));
#else
		nodes_[pos].insert(s);
#endif
	}

	#undef TEMPLATE_HEADER
	#undef TEMPLATE_ARG

	typedef Region_<3> Region;

}
#endif
