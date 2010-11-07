#ifndef _COMPLEX__09ED1116_C3FE_55a4_CA63_F84536C80D00__H
#define _COMPLEX__09ED1116_C3FE_55a4_CA63_F84536C80D00__H

#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_io.hpp"
#include "boost/thread/thread.hpp"
#include "boost/thread/barrier.hpp"

#include "molkern/__moldefs.h"

#include "molkern/complex/_region.h"
#include "molkern/forcefield/_interactions.h"
#include "molkern/complex/_archetype.h"
#include "molkern/complex/_geom_tool.h"
#include "molkern/complex/_molecule.h"
#include "molkern/complex/_parallel.h"
#include "molkern/complex/_verlet.h"
#include "molkern/complex/_linkcell.h"

namespace molkern
{
	using namespace prgkern;

	/**
	*  Комплекс взаимодействующих молекул.
	*
	*  Объект построен таким образом, что он является некой оберткой
	*  (или инструментом) для оперирования с молекулами, поскольку не хранит
	*  реальных данных. Все реальные данные (заряды, координаты, силы и т.д.)
	*  лежат в молекулах, а доступ к ним осуществляется через указатели.
	*
	*  Отсюда следуют ограничения использования данного класса. Класс не может
	*  оперировать с объектами, локализованными в различных адресных пространствах.
	*/
	template
	<
		int FORCEFIELD_TYPE = FORCEFIELD_AMBER_,
		int RESIDOME_TYPE = RESIDOME_AMBER_
	>
	class Complex_
	{
		typedef Forcefield_<FORCEFIELD_TYPE>                _Forcefield;
		typedef Residome_<RESIDOME_TYPE>                    _Residome;
		typedef Archetype_<FORCEFIELD_TYPE, RESIDOME_TYPE>  _Archetype;
		typedef Interaction                                 _Interaction;
		typedef typename _Archetype::atom_type              _Atomdata;
		typedef typename _Archetype::chain_type             _Chain;
		typedef typename _Archetype::rotamer_type           _Rotamer;
		typedef std::pair<unsigned, unsigned>               _Pair;
		typedef Box_<3, real_t>                             _Box;
		typedef Molecule_<_Archetype>                       _Molecule;
		typedef _Verlet_table<Complex_>                     _Verlet;
		typedef _Link_cell<Complex_>                        _LCell;

		/**
		 *  Базисный атомный тип для расчетов валентных и невалентных взаимодействий.
		 *  Из него вычленяются данные для пересылки в _LJAtom тип. Он содержит
		 *  существенно больше информации, хотя также выравнивается и педдится.
		 */
		struct _Atom
		{
			// Основные поля, необходимые для вычислений (заметим, что для вычислений
			// валентных взаимодействий поля charge, sigma, eps также нужны, поскольку
			// требуется рассчитывать 1-4 взаимодействия).

			vector_t X; ///< текущая координата атома
			real_t charge; ///< текущий заряд атома в единицах [a.u.q] * sqrt(ELECTRIC_FACTOR)
			real_t sigma; ///< параметр LJ 6-12 + 1-4 взаимодействия
			real_t eps; ///< параметр LJ 6-12 + 1-4 взаимодействия
			unsigned_t connect_data; ///< таблица коннектов атома
			unsigned_t insert_data; ///< тип вставки атома

			vector_t V; ///< скорость атома
				// поле V описывает точку в конфигурационном пространстве {X,V}
			vector_t F; ///< частичные силы от невалентных взаимодействий
				// используем расширенный тип (?), чтобы суммировать множество торсионных
				// вкладов без потери точности
				// поле F требуется для оптимизатора и для расчета вириала

			const _Atomdata *atomdata; ///< ссылка на основные параметры атома
					// заметим, что есть дублирование для полей sigma, eps и connect_data (они
					// явно в классе и неявно через Atomdata). Явное выделение сделано
					// для скорости копирования данных между хранилищами.

			const _Atomdata *operator->() const { return atomdata; }
					///< оператор доступа к интерфейсу _Atomdata

			/**
			 * По сути это конструтор.
			 * @param atomdata данные по атому
			 * @param atom_id - номер атома в массиве атомов (+ бит фиксации из atomdata)
			 */
			void make(const _Atomdata &atomdata, unsigned_t atom_id)
			{
				X = atomdata.X;
				charge = atomdata.charge * SQRT_ELECTRIC_FACTOR; // избегаем операции ELECTRIC_FACTOR
				sigma = atomdata.sigma;
				eps = (real_t) sqrt(atomdata.eps); // избегаем операции sqrt(epa1 * eps2) в E(VdW)
				connect_data = atomdata.connect_data;
				insert_data = atomdata.insert_data | atom_id;
				this->atomdata = &atomdata; // навешаем адрес, где лежат параметры атома
					// полагаемся на то, что массив atomdata сформирован и не будет перемещаться в памяти
			}
		};

		/**
		*  Атомный тип для расчетов невалентных LJ 6-12 взаимодействий. Тип
		*  предельно сокращается по памяти для минимизации пересылок, выравнивается,
		*  педдится (padding), чтобы эффективно использоваться на спецпроцесорах.
		*  Например, он не содержит поля F, которое предполагается размещается
		*  спецпроцессором в своей LS, чтобы не делать лишнюю пересылку.
		*  В таком сжатом виде он содержит 32 байта, что позволяет эффективно
		*  пересылать по 4 атома (128 байт).
		*/
		struct _LJAtom
		{
			// Основные поля, необходимые для вычислений.
			vector_t X; ///< текущая координата атома
			real_t charge; ///< текущий заряд атома в единицах [a.u.q] * sqrt(ELECTRIC_FACTOR)
			real_t sigma; ///< параметр LJ 6-12 + 1-4 взаимодействия
			real_t eps; ///< параметр LJ 6-12 + 1-4 взаимодействия

			// поля, содержащие информацию о 1-1 ... 1-4 контактах атома и о параметрах вставки атома
			// Формат записи следующий:
			// connect_data : каждый N-бит числа содержит информацию, есть ли связь с атомом,
			//   смещенным на N-позиций от заданного атома
			// insert_data : 0-30 биты (ref-bits) определяют номер атома в массиве атомов, они
			//   необходимы для пересылки результирующей силы в основное хранилище.
			//   31-бит (fix-bit) равняется 1 при фиксации атома в пространстве
			// Такой порядок битов позволяет исключить преобразования поля при отключении USE_FIXED_ATOMS
			unsigned_t connect_data; ///< таблица коннектов атома
			unsigned_t insert_data; ///< тип вставки атома

			vector_t F; ///< текущая сила (coul + vdw), действующая на атом
				// Силы вынесены из атомов и записываются в LJAtom, чтобы позволить алгоритмам расчета
				// невалентных и валентных взаимодействий работать независимо друг от друга. Только перед
				// переходом на следующую итерацию эти массивы сил должны сливаться.

			_LJAtom() : X(0.f), charge(0.f), sigma(1.f), eps(0.f),
				connect_data(0), insert_data(0), F(0.f) {}

			/**
			*  Это по сути конструктор LJ атома.
			*/
			void make(const _Atom &atom)
			{
				X = atom.X;
				charge = atom.charge;
				sigma = atom.sigma;
				eps = atom.eps;
				connect_data = atom.connect_data;
				insert_data = atom.insert_data;
				F = 0.0f;
			}
		};

		typedef Region_<3, _LJAtom>                         _Region;
		typedef typename _Region::node_type                 _Node;

	public:

		enum { dimension = 3 };

		typedef _Atom           atom_type;
		typedef _LJAtom         ljatom_type;
		typedef real_t          real_type;
		typedef _Archetype      archetype_type;
		typedef _Chain          chain_type;
		typedef _Rotamer        rotamer_type;
		typedef _Region         region_type;

		/**
		*  Конструирование пустого объекта. Параметры передаются как указатели на
		*  константные объекты, за исключением объекта region, размеры которого
		*  комплекс может менять.
		* @param box область взаимодействия
		* @param forcefield силовое поле
		* @param residome топология
		* @param atom_count приблизительное число атомов (для резервирования памяти)
		*/
		Complex_(const _Forcefield *forcefield, const _Residome *residome, _Region *region,
			int pH=NORMAL_PH_WATER, real_t density=NORMAL_WATER_DENSITY)
		: forcefield_(forcefield), residome_(residome),
			region_(region), pH_(pH), density_(density)
		{
#ifdef USE_VERLET_TABLE
			nr_integrator_ = new _Verlet(this);
#else
			nr_integrator_ = new _LCell(this);
#endif
		}

		~Complex_()
		{
			for (unsigned i=0,sz=archetypes_.size(); i<sz; i++) delete archetypes_[i];
			for (unsigned i=0,sz=molecules_.size(); i<sz; i++) delete molecules_[i];
			delete nr_integrator_;
		}

		/**
		*  Загрузить образец молекулы из файла.
		* @param filename полное имя файла (с путем)
		* @param freedom_type тип свободы молекулы
		* @param count количество молекул для заданного архетипа
		* @param altpos альтернативная позиция загрузки атомов (для *.pdb)
		* @return идентификатор загруженного образца молекулы или -1
		*/
		int load(_I2T<MOLECULE_>, const std::string &filename, unsigned freedom_type,
			unsigned count, char altpos='A');
		int load(_I2T<WATER_>, const std::string &solution, unsigned freedom_type);

		/**
		* Сохранить комплекс в файл. Формат файла определяется расширением имени файла.
		* @param filename полное с путем имя файла
		* @param prn_hydrogens печатать водороды (true) или нет
		* @param header комментарий, который вставляется в начало файла
		* @return количество сделанных записей в файл (для бинарных файлов)
		*/
		void save(const std::string &filename, bool prn_water=false,
			bool prn_hydrogens=true, const std::string &header=_S("")) const;

		/**
		*  Сохранение заголовка файла выполняется отдельно, чтобы позволить писать в
		*  файл непрерывный поток данных разных молекул, без прерывания потока
		*  какими-то заголовками.
		*/
		void write_header(_I2T<FORMAT_PDB_>, std::ofstream &file, const std::string &header=_S("")) const;
		void write_header(_I2T<FORMAT_HIN_>, std::ofstream &file, const std::string &header=_S("")) const;
		void write_header(_I2T<FORMAT_MOL2_>,std::ofstream &file, const std::string &header=_S("")) const;
		void write_header(_I2T<FORMAT_BMM_>, std::ofstream &file, const std::string &header=_S("")) const;

		template <int FORMAT> void save(_I2T<FORMAT>, std::ofstream &file, bool prn_water=false,
			bool prn_hydrogens=true, const std::string &header=_S("")) const
		{
			const _Atom *atoms__ = &atoms_[0];
			for (unsigned iarchetype=0,sz=archetypes_.size(); iarchetype<sz; iarchetype++)
			{
				const _Archetype *archetype = archetypes_[iarchetype];
				if (!prn_water && archetype->is_solution()) continue;

				unsigned molecule_count = count(MOLECULE, iarchetype);
				for (unsigned i=0; i<molecule_count; i++)
				{
					unsigned ndx = am_matrix_[iarchetype][i];
					const _Molecule *molecule = molecules_[ndx];
					molecule->save(_I2T<FORMAT>(), file, atoms__ + atom_start_[ndx], prn_hydrogens);
				}
			}
		}

		/**
		*  Размещение молекул и раствора, входящих в комплекс.
		* @param enable_clash разрешение иметь клеши при позиционировании молекул
		* @param iterations число итераций размещения каждой молекул при клешировании
		*/
		void build(bool enable_clash, unsigned iterations=DEFAULT_DISPOSE_ITERATIONS);

		/**
		* @brief energy & derivations of entire complex
		*/
		_E(real_t) U(bool make_print=YES_PRINT) const;
		_E(real_t) dU__dX();
		_E(real_t) dU__dX(_I2T<POSITION_>);

		/**
		* @brief extracts data
		*/
		DEFINE_VECTOR_ACCESS_FUNCTION(ARCHETYPE_, _Archetype, archetypes_);
		DEFINE_VECTOR_ACCESS_FUNCTION(MOLECULE_,  _Molecule,  molecules_ );
		DEFINE_VECTOR_ACCESS_FUNCTION(ATOM_,      _Atom,      atoms_     );

		DEFINE_VECTOR_ACCESS_FUNCTION(XPOSITION_, real_t, x_);
		DEFINE_VECTOR_ACCESS_FUNCTION(GPOSITION_, real_t, g_);

		_Box get(_I2T<BOX_>) const { return region_->get(BOX); }
		const _Region *get(_I2T<REGION_>) const { return region_; }
		_Region *get(_I2T<REGION_>) { return region_; }

		//----------------------------------------------------------------------------------------
		//                         функции для работы с ротамерами
		//----------------------------------------------------------------------------------------
		void read_md_position(const vector_t *x, const vector_t *v)
		{
			unsigned atom_count = atoms_.size();
			for (unsigned i=0; i<atom_count; i++)
			{
				_Atom &atom = atoms_[i];
				atom.X = x[i];
				atom.V = v[i];
			}
		}

		unsigned read(_I2T<POSITION_>, const real_t *x, _Atom *atoms)
		{
			unsigned cnt = 0;
			for (unsigned i=0,sz=molecules_.size(); i<sz; i++)
				cnt += molecules_[i]->read(POSITION, x + cnt, atoms + atom_start_[i]);
			return cnt;
		}

		unsigned write(_I2T<GRADIENT_>, real_t *g, const _Atom *atoms)
		{
			unsigned cnt = 0;
			for (unsigned i=0,sz=molecules_.size(); i<sz; i++)
				cnt += molecules_[i]->write(GRADIENT, g + cnt, atoms + atom_start_[i]);
			return cnt;
		}

		/// полное число степеней свободы системы (включая и движение как целое)
		unsigned count(_I2T<FREEDOM_>) const
		{
			unsigned count = 0;
			for (unsigned i=0,sz=molecules_.size(); i<sz; i++)
				count += molecules_[i]->count(FREEDOM);
			return count;
		}

		/// возвращает число свободных атомов комплекса
		unsigned count(_I2T<FREE_>) const
		{
		#ifdef USE_FIXED_ATOMS
			unsigned free_atom_count = 0;
			for (unsigned i=0,sz=molecules_.size(); i<sz; i++)
				free_atom_count += molecules_[i]->count(_I2T<FREE_>());
			return free_atom_count;
		#else
			return atoms_.size();
		#endif
		}

		/// число молекул для заданного архетипа
		unsigned count(_I2T<MOLECULE_>, unsigned iarchetype) const
		{ return am_matrix_[iarchetype].size();	}

		/**
		*  Конвертирует итератор по индексам в итератор по атомам молекулы.
		*  Позволяет пользователю, который создал группы атомов (через совокупность
		*  индексов) снаружи класса молекулы, получать соответствующие атомы группы.
		* @param it итератор (начала или конца) последовательности индексов
		* @return итератор по атомам молекулы
		*/
		template <typename _Iterator>
		array_iterator<_Atom, _Iterator> make_iterator(_Iterator it)
		{ return array_iterator<_Atom, _Iterator>(&atoms_[0], it); }

		template <typename _Iterator>
		const_array_iterator<_Atom, _Iterator> make_iterator(_Iterator it) const
		{ return const_array_iterator<_Atom, _Iterator>(&atoms_[0], it); }


//	#ifdef USE_VERLET_TABLE
//		unsigned count(_I2T<REBUILD_>) const { return vt_.count(REBUILD); }
//		real_t average(_I2T<RSKIN_>) const { return vt_.average(RSKIN); }
//	#endif

	protected:

		/**
		*   Заполнить все свободное пространство между молекулами в ящике раствором.
		* @return идентификатор молекулы
		*/
		void dispose_(_I2T<WATER_>, const _Box &box, unsigned iarchetype, real_t density);

		/**
		*  Позиционирование молекулы, задаваемой архетипом, в области пространства.
		* @param molecule молекула, чьи координаты используются для результата
		* @param iarchetype идентификатор архетипа молекулы
		* @param enable_clash разрешение иметь клеши в заданной области
		* @param iterations число попыток поиска безклешевой позиции
		* @return найдена или нет позиция без клеша
		*/
		boost::tuple<vector_t, vector_t> dispose_(_I2T<MOLECULE_>, const _Box &box, unsigned iarchetype,
			vector_t *X, bool enable_rotation, bool enable_clash, unsigned iterations);

		/**
		*  (Пере)заполнение грида атомами и прямой расчет энергий и сил пар.
		*  Функция оптимизирована только для параллельного выполнения.
		*/
		_E(real_t) U_(_I2T<PAIR_>, bool make_print=true) const;

	private:

		const _Forcefield *forcefield_; // указатель на силовое поле
		const _Residome *residome_; // указатель на топологию
		_Region *region_; // область взаимодействия
		int pH_; // pH раствора
		real_t density_; // плотность раствора (число молекул в 1 A**3)

		std::vector<_Archetype*> archetypes_; // указатели на типы молекул комплекса
		std::vector<_Molecule*> molecules_;  // указатели на молекулы комплекса
		std::vector<unsigned> atom_start_;  // стартовые номера атомов молекул в массиве атомов

		std::vector<std::vector<int> > am_matrix_;
			// матрица, хранящая идентификаторы молекул для каждого архетипа
			// она используется для группирования молекул, позволяет избежать печати
			// заданных типов молекул и т.д.

		std::vector<_Atom> atoms_; // все атомы комплекса
		std::vector<_Rotamer> rotamers_; // цепи всех входящих ротамеров
			// (!) цепи и ротамеры перенумерованны относительно комплекса

		// Нити инициализируются после того, как комплекс полностью построен.
		// Они используются только для расчета невалентных парных взаимодействий.
		// При деструкции комплекса, нити разрушаются.
		//_ParallWorker parall_worker_; // обеспечивает расчет LJ взаимодействий в параллельном режиме
		std::vector<real_t> x_, g_; // массивы обобщенных координат для оптимизатора

		near_range_integrator<Complex_> *nr_integrator_;
	};

	#define TEMPLATE_HEADER  template <int FORCEFIELD_TYPE, int RESIDOME_TYPE>
	#define TEMPLATE_ARG     FORCEFIELD_TYPE, RESIDOME_TYPE

	TEMPLATE_HEADER
	inline int Complex_<TEMPLATE_ARG>
	::load(_I2T<MOLECULE_>, const std::string &filename, unsigned freedom_type, unsigned count,
		char altpos)
	{
		std::ifstream file(filename.c_str());
		if (!file)
		{
			std::string msg = _S("[ERROR] can't open file ") + filename;
			PRINT_BREAK(msg);
		}

		unsigned archetype_id = archetypes_.size();

		_Archetype *archetype = new _Archetype(forcefield_, residome_, freedom_type);
		if (archetype->load(filename, altpos))
		{
			// Найдем максимальный ящик, в который помещаются все молекулы. Заметим, что такая стратегия
			// поиска (без проверки возможных вращений молекул) наиболее соответствует стандартным задачам.
			// Так например, она не вращает протеин, который может иметь существенно разные размеры
			// по разным осям, тем самым эта стратегия не увеличивает размеры задачи. Также для лигандов,
			// которые обычно маленькие, она не делает вращений, так как они в любом случае попадают в ящик.
			region_->enlarge(archetype->get(BOX));

			archetypes_.push_back(archetype);
			am_matrix_.resize(archetype_id + 1); // зафиксируем наличие нового архетипа
			am_matrix_[archetype_id].resize(count, -1); // забъем NULL указателями (передача числа молекул),
				// а инициализируем их только после удачного размещения атомов

			return archetype_id;
		}
		return -1;
	}

	TEMPLATE_HEADER
	INLINE int Complex_<TEMPLATE_ARG>
	::load(_I2T<WATER_>, const std::string &solution, unsigned freedom_type)
	{
		_Archetype *water = new _Archetype(forcefield_, residome_, freedom_type);
		if (water->load(_I2T<WATER_>(), solution))
		{
			unsigned archetype_id = archetypes_.size();
			archetypes_.push_back(water);
			am_matrix_.resize(archetype_id + 1); // зафиксируем наличие нового архетипа
			am_matrix_[archetype_id].resize(1, -1); // забъем NULL указателями (передача числа молекул),
				// а инициализируем их только после удачного размещения атомов

			return archetype_id;
		}
		return -1;
	}

	TEMPLATE_HEADER
	inline void Complex_<TEMPLATE_ARG>
	::save(const std::string &filename, bool prn_water, bool prn_hydrogens,
		const std::string &header) const
	{
		std::ofstream file(filename.c_str());
		if (!file)
		{
			std::string msg = _S("can't open file ") + filename;
			PRINT_ERR(msg);
		}
		std::string ext = extension(filename);

		if (ext == _S(".pdb") || ext == _S(".ent"))
		{
			write_header(_I2T<FORMAT_PDB_ >(), file, header);
			save(_I2T<FORMAT_PDB_ >(), file, prn_water, prn_hydrogens, header);
		}
		else if (ext == _S(".hin") )
		{
			write_header(_I2T<FORMAT_HIN_ >(), file, header);
			save(_I2T<FORMAT_HIN_ >(), file, prn_water, prn_hydrogens, header);
		}
		else if (ext == _S(".mol2"))
		{
			write_header(_I2T<FORMAT_MOL2_>(), file, header);
			save(_I2T<FORMAT_MOL2_>(), file, prn_water, prn_hydrogens, header);
		}
		else if (ext == _S(".bmm"))
		{
			write_header(_I2T<FORMAT_BMM_ >(), file, header);
			save(_I2T<FORMAT_BMM_ >(), file, prn_water, prn_hydrogens, header);
		}

		SAVED_OK_MESSAGE(filename);
	}

	TEMPLATE_HEADER
	INLINE void Complex_<TEMPLATE_ARG>
	::write_header(_I2T<FORMAT_PDB_>, std::ofstream &file, const std::string &header) const
	{
		if (header.length() != 0) file << "HEADER    " << header << std::endl;
	}

	TEMPLATE_HEADER
	INLINE void Complex_<TEMPLATE_ARG>
	::write_header(_I2T<FORMAT_HIN_>, std::ofstream &file, const std::string &header) const
	{
		if (header.length() != 0) file << "; " << header << std::endl;
	}

	TEMPLATE_HEADER
	INLINE void Complex_<TEMPLATE_ARG>
	::write_header(_I2T<FORMAT_MOL2_>, std::ofstream &file, const std::string &header) const
	{
		if (header.length() != 0) file << "@<TRIPOS>COMMENT " << header << std::endl;
	}

	TEMPLATE_HEADER
	INLINE void Complex_<TEMPLATE_ARG>
	::write_header(_I2T<FORMAT_BMM_>, std::ofstream &file, const std::string &header) const
	{
		// сохраним ящик в целочисленном представлении
		float radius = (float)_Interaction::interaction_radius();
		file.write((char*)&radius, sizeof(float));

		typename _Region::index_type boxsz = region_->box_size();
		unsigned sx = (unsigned)boxsz[0];
		unsigned sy = (unsigned)boxsz[1];
		unsigned sz = (unsigned)boxsz[2];
		file.write((char*)&sx, sizeof(unsigned));
		file.write((char*)&sy, sizeof(unsigned));
		file.write((char*)&sz, sizeof(unsigned));

		// сохраним число атомов комплекса
		unsigned atom_count = count(_I2T<ATOM_>());
		file.write((char*)&atom_count, sizeof(unsigned));
	}

	TEMPLATE_HEADER
	inline void Complex_<TEMPLATE_ARG>
	::build(bool enable_clash, unsigned iterations)
	{
		TIME_TESTING_START("", 1);
		_S msg = _S("\nBuilding of complex is started ...");
		PRINT_MESSAGE(msg);

		//--------------------------------------------------------------------------
		//         размещение молекул внутри ящика с возможным его ростом
		//--------------------------------------------------------------------------
		for (unsigned iarchetype=0,sz=archetypes_.size(); iarchetype<sz; ++iarchetype)
		{
			_Archetype *archetype = archetypes_[iarchetype];
			if (archetype->is_solution()) continue; // раствор размещаем в последнюю очередь
			unsigned molecule_count = count(MOLECULE, iarchetype);

			archetype->build(MOLECULE, pH_);

			unsigned nf = archetype->count(FREEDOM, archetype->get(FREEDOM_TYPE));
			unsigned prev_xsize = x_.size();
			x_.resize(prev_xsize + nf * molecule_count, 0.f);
			g_.resize(prev_xsize + nf * molecule_count, 0.f);

			unsigned atom_count = archetype->count(ATOMDATA);
			vector_t X[atom_count]; // позиция атомов вновь размещенной молекулы
			for (unsigned imolecule=0; imolecule<molecule_count; ++imolecule)
			{
				bool enable_rotation = (bool)imolecule; // нулевую молекулу не вращаем

				boost::tuple<vector_t, vector_t> tuple__ = dispose_(MOLECULE, region_->get(BOX),
					iarchetype, &X[0], enable_rotation, enable_clash, iterations);

				unsigned prev_atom_size = atoms_.size();
				atoms_.resize(prev_atom_size + atom_count);
				const _Atomdata *atomdata = archetype->get(ATOMDATA);
				for (unsigned i=0; i<atom_count; i++)
				{
					_Atom &atom = atoms_[prev_atom_size + i];
					atom.make(atomdata[i], prev_atom_size + i);
					atom.X = X[i];
				}

				unsigned prev_molecule_size = molecules_.size();

				_Molecule *molecule = new _Molecule(archetype, archetype->get(FREEDOM_TYPE));
				molecules_.push_back(molecule);
				atom_start_.push_back(prev_atom_size);
				am_matrix_[iarchetype][imolecule] = prev_molecule_size;

				vector_t angle = tuple__.get<0>();
				vector_t ranx = tuple__.get<1>();
				real_t *x__ = &x_[prev_xsize + imolecule * nf]; // адрес записи обобщенных координат
				unsigned freedom_type = archetype->get(FREEDOM_TYPE);
				archetype->write_position(MOLECULE, x__, freedom_type, &atoms_[prev_atom_size], angle, ranx);
			}
		}
		//--------------------------------------------------------------------------
		//                  заполнение промежутков атомами раствора
		//--------------------------------------------------------------------------
		for (unsigned iarchetype=0,sz=archetypes_.size(); iarchetype<sz; ++iarchetype)
		{
			_Archetype *archetype = archetypes_[iarchetype];
			if (!archetype->is_solution()) continue; // молекулы построены выше

			dispose_(WATER, region_->get(BOX), iarchetype, density_);

			archetype->build(WATER);
				// построение (связей, углов и т.д.) для всех вновь добавленных молекул воды
				// в отличие от молекулы, построение делается после размещения, поскольку
				// размещение изменяет число объектов

			const _Atomdata *atomdata = archetype->get(ATOMDATA);
			unsigned atomdata_count = archetype->count(ATOMDATA);
				// заново возьмем адрес, так как atomdata было перестроено в build(WATER)

			// перекопируем все координаты
			unsigned prev_atom_size = atoms_.size();
			atoms_.resize(prev_atom_size + atomdata_count);
			for (unsigned i=0; i<atomdata_count; i++)
			{
				_Atom &atom = atoms_[prev_atom_size + i];
				atom.make(atomdata[i], prev_atom_size + i);
				atom.X = atomdata[i].X;
			}

			unsigned prev_molecule_size = molecules_.size();
			_Molecule *molecule = new _Molecule(archetype, archetype->get(FREEDOM_TYPE));
			molecules_.push_back(molecule);
			atom_start_.push_back(prev_atom_size);
			am_matrix_[iarchetype][0] = prev_molecule_size;

			unsigned nf = archetype->count(FREEDOM, archetype->get(FREEDOM_TYPE));
			unsigned prev_xsize = x_.size();
			x_.resize(prev_xsize + nf, 0.f);
			g_.resize(prev_xsize + nf, 0.f);

			real_t *x__ = &x_[prev_xsize]; // адрес записи обобщенных координат
			unsigned freedom_type = archetype->get(FREEDOM_TYPE);
			archetype->write_position(WATER, x__, freedom_type, &atoms_[prev_atom_size]);
		}

		msg = _S("Building of complex is finished ...");
		PRINT_MESSAGE(msg);
		TIME_TESTING_FINISH;

		//--------------------------------------------------------------------------
		//                    дополнительный контроль данных
		//--------------------------------------------------------------------------
		real_t full_charge = 0.;
		for (unsigned i=0,sz=atoms_.size(); i<sz; i++) full_charge += atoms_[i].charge;

		full_charge /= SQRT_ELECTRIC_FACTOR; // коррекция к заряду в [a.e.q]
		if (fabs(full_charge) > 0.5)
		{
			std::string msg = make_string("[WARNING] Full charge of molecule is %10.3f \n",
				(float)full_charge);
			msg += _S("  You must include some ions to avoid problems with far coulomb calculations");
			PRINT_MESSAGE(msg);
		}

		nr_integrator_->resize(_Interaction::interaction_radius() + global_rskin_width, density_);
		nr_integrator_->update(YES_PRINT);
	}

	TEMPLATE_HEADER
	inline boost::tuple<vector_t, vector_t> Complex_<TEMPLATE_ARG>
	::dispose_(_I2T<MOLECULE_>, const _Box &box, unsigned iarchetype, vector_t *X,
		bool enable_rotation, bool enable_clash, unsigned iterations)
	{
		const _Archetype *archetype = archetypes_[iarchetype];

		//--------------------------------------------------------------------------
		// Функция использует простой метод O(N**2) по поиску клешей, что
		// обусловлено тем, что молекулы, которые вставляются в область являются
		// маленькими. Вплоть до размера в 400 атомов, такой алгоритм работает
		// быстрее алгоритма с поиском соседей
		//--------------------------------------------------------------------------
		std::vector<unsigned> atoms_inside_box;
			// массив атомов, которые попали в заданную область, и с которыми
			// нужно избежать клеша

		unsigned atoms_inside_box_count = 0; // число атомов, попавших в ящик
		if (!enable_clash)
		{
			atoms_inside_box.reserve((unsigned) (DEFAULT_ATOM_DENSITY * box.volume()));
				// резервируем памяти на приблизительное количество атомов в области

			for (unsigned i=0,sz=atoms_.size(); i<sz; i++)
			{
				if (box.is_included(atoms_[i].X)) atoms_inside_box.push_back(i);
					// при использовании подобной проверки предполагается,
					// что атомы, еще не зафиксированные в комплексе имеют
					// координаты вне любого разумного ящика X = {inf, inf, inf}
			}
			atoms_inside_box_count = atoms_inside_box.size();
		}

		unsigned atom_count = archetype->count(_I2T<ATOM_>());
		const _Atomdata *atomdata = archetype->get(ATOMDATA);

		vector_t geom_center = calculate(_I2T<GEOM_CENTER_>(), atomdata,
				range_iterator(0), range_iterator(atom_count));

		unsigned iteration = 0;
		while (iteration++ < iterations)
		{
			// Вращаем случайным образом молекулу около ее центра, записывая результат
			// на нужное место. Там и оставим, если все нормально. Ротация на 0-ой
			// угол нужна только для перезаписи данных на место молекулы.
			vector_t angle = ran0(vector_t(2*M_PI, M_PI, 2*M_PI));
			if (!enable_rotation) angle = vector_t(0., 0., 0.);

			// XYZ_ROTATOR для совместимости с построением обобщенных координат
			Rotator<XYZ_ROTATOR_, real_t> rotator(angle, geom_center);
			for (unsigned i=0; i<atom_count; i++) X[i] = rotator(atomdata[i].X);

			// найдем ее ящик после вращения (в зависимости от формы молекулы он
			// может меняться немного или сильно)
			Box_<3, real_t> mol_box(vector_t(INFINITY), vector_t(-INFINITY));
			for (unsigned i=0; i<atom_count; i++)
			{
				mol_box.bottom() = min(mol_box.bottom(), X[i]);
				mol_box.top()    = max(mol_box.top(),    X[i]);
			}
			{
				_S msg = make_string("molecule box is : %s", make_string(mol_box));
				PRINT_MESSAGE(msg);
			}
			// Позиционируем ящик случайным образом внутри заданного большого ящика,
			// так чтобы он не вышел за его границы. Для этого создаем безопасный
			// ящик, размер которого уменьшен на длину ящика молекулы.
			box_t safe_box = box;
			vector_t safe_len = mol_box.length();
			safe_len *= 0.5;
			safe_box.top() -= safe_len;
			safe_box.bottom() += safe_len;

			if (safe_box.volume() <= 0)
			{
				_S msg = make_string("target box is too small : %s", make_string(safe_box));
				PRINT_ERR(msg);
			}

			vector_t ranx = ran0(safe_box.bottom(), safe_box.top());
			ranx -= geom_center;
			for (unsigned i=0; i<atom_count; i++) X[i] += ranx;

			if (enable_clash || atoms_inside_box_count==0) goto LABEL_NO_CLASH;
				// при разрешении клешей, или отсутствия атомов в заданной области
				// не делаем проверок и новых попыток позиционирования

			//------------------------------------------------------------------------
			//                      проверка наличия клешей
			//------------------------------------------------------------------------
			for (unsigned i=0; i<atom_count; i++)
			{
				real_t sigma = atomdata[i].sigma;

				for (unsigned i__=0; i__<atoms_inside_box_count; i__++)
				{
					real_t sigma__ = atoms_[atoms_inside_box[i__]].sigma;
					vector_t X__ = atoms_[atoms_inside_box[i__]].X;
					if (distance2(X[i], X__) < sqr(sigma + sigma__)) goto LABEL_CLASH;
				}
			}

		LABEL_NO_CLASH: //--------------------- нет клешей -------------------------
			{
				_S msg = make_string("%s<%s> has been posed inside of target box in %s",
					archetype->name(_I2T<FILE_>()).c_str(),
					archetype->name(_I2T<MOLECULE_>()).c_str(),
					make_string(ranx).c_str()
				);
				PRINT_MESSAGE(msg);

				return boost::tuple<vector_t, vector_t>(angle, ranx);
			}
		LABEL_CLASH: ;
		}

		_S msg = make_string("have been done %d iterations to avoid clash without succeed", iteration);
		PRINT_ERR(msg);

		return boost::tuple<vector_t, vector_t>(0.f, 0.f);
	}

	TEMPLATE_HEADER
	inline void Complex_<TEMPLATE_ARG>
	::dispose_(_I2T<WATER_>, const _Box &box, unsigned iarchetype, real_t density)
	{
		_Archetype *water_archetype = archetypes_[iarchetype];
		if (water_archetype == 0) return;

		// Радиус клеширования определяет только размер ячейки, на который делится ящик при поиске
		// соседей. Большой радиус не изменяет число молекул отмеченных для клеширования, но
		// увеличивает время счета. Устанавливаем ящик с ячейками равными WATER_CLASH_RADIUS.

		vector_t T = region_->translation_vector();
		index_<3, unsigned> region_size
		(
			(unsigned)ceil(box.length(0) / WATER_CLASH_RADIUS),
			(unsigned)ceil(box.length(1) / WATER_CLASH_RADIUS),
			(unsigned)ceil(box.length(2) / WATER_CLASH_RADIUS)
		);

		_Region grid(region_size, WATER_CLASH_RADIUS, &T);
		_Box box__ = grid.get(BOX);

		typedef typename _Residome::residue_type      _Residue;
		typedef typename _Archetype::atom_type        _Atomdata;
		typedef typename _Region::region_vector_type  region_vector_t;
		typedef typename _Region::region_key_type     region_key_t;

		std::string molname = water_archetype->name(_I2T<MOLECULE_>());

		const _Residue *residue = residome_->get_data(molname);
		unsigned chain_count = residue->count(CHAIN);

		_Box residue_box = residue->box();
		vector_t cur_len = residue_box.length();

		if (chain_count == 1)
		{
			// Генерим новый ящик для растворов, у которых в AMBER *.lib файле нет ящиков.
			// Обычно это одиночные молекулы.
			unsigned atom_count = residue->count(ATOM);
			density /= atom_count; // плотность атомов кислорода, по которым мы отсекаем лишние
				// молекулы меньше на соответствующий фактор

			real_t len__ = 1. / pow(density, 1./3.);
			cur_len = vector_t(len__, len__, len__);
			residue_box = _Box(cur_len * (-0.5), cur_len * 0.5);
		}

		vector_t len = box__.length();
		vector_t X0 = box__.bottom() - residue_box.bottom();
			// смещение ящиков относительно друг друга

		unsigned nx = (unsigned) ceil(len[0] / cur_len[0]);
		unsigned ny = (unsigned) ceil(len[1] / cur_len[1]);
		unsigned nz = (unsigned) ceil(len[2] / cur_len[2]);
			// параметры дублирования по осям
		unsigned repeat_factor = nx * ny * nz;
			// коэффициент дублирования данных
		unsigned atomdata_count = water_archetype->count(ATOM);
		unsigned atom_count_per_chain = atomdata_count / chain_count;

		// функция построена так, что работает только с данными, загруженными
		// из lib файла. Функция полагается на допущения по размещению молекул,
		// принятые внутри класса комплекса, то есть по поводу хранения раствора
		// в нулевых элементах массивов archetypes_ & molecules_

		//--------------------------------------------------------------------------
		//                       вставим воду в ящик
		//--------------------------------------------------------------------------
		// Вставляем в регион информацию (X, sigma, chain=connect_data)
		// Используем такой порядок добавления атомов (работаем с атомами одной ячейки
		// и только после переходим к другой ячейке), чтобы выполнялось условие,
		// что близкие атомы в пространстве имеют близкие номера.
		// Возможно работать только с атомами != водородам, поскольку водороды
		// локализуются близко к атомам (1А), и только оттягивают на себя время.

		vector_t X__, Y; // текущий вектор трансляции
		const _Atomdata *atomdatas = water_archetype->get(ATOMDATA);

		// массив для хранения идентификаторов исключаемых цепей
		std::vector<bool> killdata(chain_count * repeat_factor, true);

		_LJAtom ljatom; unsigned count = 0, pos;
		for (unsigned ix=0; ix<nx; ix++)
		for (unsigned iy=0; iy<ny; iy++)
		for (unsigned iz=0; iz<nz; iz++)
		{
			// вектор трансляции атомов заданной ячейки
			X__[0] = X0[0] + cur_len[0] * ix;
			X__[1] = X0[1] + cur_len[1] * iy;
			X__[2] = X0[2] + cur_len[2] * iz;

			for (unsigned i=0; i<atomdata_count; i++, count++)
			{
				if (atomdatas[i].is_hydrogen()) continue; // игнорируем атомы водорода
				if (atomdatas[i].is_pseudo()) continue; // игнорируем атомы EPW

				ljatom.X = X__ + atomdatas[i].X;
				ljatom.sigma = atomdatas[i].sigma;
				ljatom.insert_data = count;

				if (grid.insert(NO_PERIODIC, ljatom, ljatom.X))
				{
					unsigned nchain = count / atom_count_per_chain;
					killdata[nchain] = false;
				}
			}
		}

		//--------------------------------------------------------------------------------------
		// Проверяем атомы молекул на клеши в их положении и удаляем соответствующие атомы воды.
		//--------------------------------------------------------------------------------------
		typedef typename _Region::key_type      _Key;
		typedef typename _Region::index_type      _Index;
		typedef typename _Node::const_iterator  _Iterator;
		unsigned neigbours[_Region::neighbours]; // соседи текущей ячейки региона

		const _Node *nodes = grid.get(NODE);
		for (unsigned i=0,sz=molecules_.size(); i<sz; i++)
		{
			if ((*molecules_[i])->is_solution()) continue;

			unsigned atom_count = (*molecules_[i])->count(ATOM);
			const _Atom *atom = &atoms_[0] + atom_start_[i];
			for (unsigned i=0; i<atom_count; i++)
			{
				if (atom[i]->is_hydrogen()) continue; // игнорируем атомы водорода

				ljatom.X = atom[i].X;
				ljatom.sigma = atom[i].sigma;

				region_vector_t X = grid.make_insert_image(ljatom.X);
				pos = grid.make_insert_index(X);
				{
					const _Node &node = nodes[pos];
					unsigned neig_count = node.get_neigbours(YES_PERIODIC, &neigbours[0]);
						// чтение данных методом YES_PERIODIC, чтобы избежать клешей на границах ячейки
					_Index key = grid.index3(pos);
					for (unsigned nei=0; nei<neig_count; nei++)
					{
						unsigned pos__ = neigbours[nei];
						const _Node &node__ = nodes[pos__];

						piecewise_vector<const _LJAtom> atoms__(node__.get(ATOM), node__.count(ATOM));
						_Index key__ = grid.index3(pos__);
						_Key direction(key__[0] - key[0], key__[1] - key[1], key__[2] - key[2]);
						grid.make_piecewise_vector(&atoms__, direction, pos__);

						typedef typename piecewise_vector<const _LJAtom>::iterator  _Iterator;
						_Iterator it = atoms__.begin(), ite = atoms__.end();

						for (; it!=ite; ++it)
						{
							int count = (int)((*it).insert_data);
							const _LJAtom &ljatom__ = *it;
							vector_t R
							(
								ljatom__.X[0] - ljatom.X[0],
								ljatom__.X[1] - ljatom.X[1],
								ljatom__.X[2] - ljatom.X[2]
							);
							grid.make_nearest_image_vector(R[0], R[1], R[2]);
							if (R.length2() > sqr(ljatom.sigma + ljatom__.sigma)) continue;

							unsigned nchain = count / atom_count_per_chain;
							killdata[nchain] = true;
						}
					}
				}
			}
		}

		//--------------------------------------------------------------------------
		//  перекопируем все непомеченные атомы в новый массив и переразметим цепи
		//--------------------------------------------------------------------------
		std::vector<_Atomdata> atomdatas__;
		atomdatas__.reserve(repeat_factor * atomdata_count); // разместили достаточно памяти

		int sid = 0; // упорядочиваем sid
		int res_seq = 0; // под этим номером заново упорядочиваются молекулы воды
		count = 0;

		unsigned nchain = 0;
		for (unsigned ix=0; ix<nx; ix++)
		for (unsigned iy=0; iy<ny; iy++)
		for (unsigned iz=0; iz<nz; iz++)
		{
			X__[0] = X0[0] + cur_len[0] * ix;
			X__[1] = X0[1] + cur_len[1] * iy;
			X__[2] = X0[2] + cur_len[2] * iz;
				// заново пересчитаем координаты вставляемой воды

			for (unsigned i=0; i<atomdata_count; i+=atom_count_per_chain, nchain++)
			{
				if (killdata[nchain]) continue;

				for (unsigned k=0; k<atom_count_per_chain; k++)
				{
					_Atomdata atomdata = atomdatas[i + k];
					atomdata.sid = sid++;
					atomdata.X += X__;
					atomdata.res_seq = res_seq;
					atomdatas__.push_back(atomdata);
				}
				res_seq++;
			}
		}

		atomdata_count = water_archetype->update(ATOMDATA, atomdatas__);
			// перенесем данные и уничтожим ненужные массивы

		PRINT_MESSAGE(_S("Have been rebuilt -> solution atomdata  (") + itoa(atomdata_count) + _S(")"));
	}

	TEMPLATE_HEADER
	inline _E(real_t) Complex_<TEMPLATE_ARG>
	::U(bool make_print) const
	{
		_E(real_t) energy = 0.;
		const _Atom *atoms__ = &atoms_[0];

		for (unsigned i=0,sz=molecules_.size(); i<sz; i++)
			energy += molecules_[i]->U(atoms__ + atom_start_[i], make_print);

		unsigned nf = count(FREEDOM);
		if (nf) energy += U_(_I2T<PAIR_>(), make_print);

		if (make_print)
		{
			// нормировка энергии на число степеней свободы при выводе
			if (nf) { PRINT_MESSAGE(make_string("complex : nfreedom = %d,  energy = %12.5e,  "
				"<energy/nfreedom> = %12.5e \n", nf, (float)energy, (float)(energy / nf))); }
			else { PRINT_MESSAGE(make_string("complex : nfreedom = %d,  energy = %12.5e\n",
				nf, (float)energy)); }
		}
		return energy;
	}

	TEMPLATE_HEADER
	INLINE _E(real_t) Complex_<TEMPLATE_ARG>
	::dU__dX()
	{
		_Atom *atoms__ = &atoms_[0];

		//-------------------------------------------------------------------------
		//            обнуление данных перед каждой итерацией
		//-------------------------------------------------------------------------
		_E(real_t) energy = 0.;
		for (unsigned i=0,sz=atoms_.size(); i<sz; i++) atoms__[i].F = (real_t)0.;

		//-------------------------------------------------------------------------
		//                             счет
		//-------------------------------------------------------------------------
		for (unsigned i=0,sz=molecules_.size(); i<sz; i++)
			energy += molecules_[i]->dU__dX(atoms__ + atom_start_[i]);

		energy += nr_integrator_->dU__dX();
		return energy;
	}

	TEMPLATE_HEADER
	inline _E(real_t) Complex_<TEMPLATE_ARG>
	::U_(_I2T<PAIR_>, bool make_print) const
	{
		_E(real_t) coul_energy = 0., vdw_energy = 0.;
			// накопители энергии должны иметь расширенный тип

	#if !defined(SKIP_COUL_ENERGY) || !defined(SKIP_VDW_ENERGY)

		unsigned pairs_count = 0;
		//--------------------------------------------------------------------------
		//          очистка, вставка атомов в грид, неявное формирование пар
		//--------------------------------------------------------------------------
		typedef typename _Region::region_vector_type  region_vector_t;
		typedef typename _Region::region_key_type  region_key_t;

		region_->clear(); // очистка от предыдущего расчета
		unsigned atom_count = atoms_.size();

		_LJAtom ljatom;
		for (unsigned i=0; i<atom_count; i++)
		{
			ljatom.make(atoms_[i]);
			region_->insert(YES_PERIODIC, ljatom, atoms_[i].X);
		}

		real_t rmax2 = sqr(_Interaction::interaction_radius());
		unsigned neigbours[_Region::neighbours]; // соседи текущей ячейки региона

		//--------------------------------------------------------------------------
		//                        поиск соседей
		//--------------------------------------------------------------------------
		_E(real_t) coul_energy__, max_coul_energy = 0.;
		_E(real_t) vdw_energy__, max_vdw_energy = 0.;
			// ищем максимум абсолютного значения, значит начальное значение = 0

		typedef typename _Node::const_iterator _Iterator;

		unsigned node_count = region_->count(NODE);
		const _Node *nodes = region_->get(NODE);
		for (unsigned pos=0; pos<node_count; pos++)
		{
			const _Node &node = nodes[pos];
			unsigned nei_count = node.get_neigbours(YES_PERIODIC, &neigbours[0]);
			for (unsigned nei=0; nei<nei_count; nei++)
			{
				unsigned pos__ = neigbours[nei];
				if (pos__ < pos) continue;
				const _Node &node__ = nodes[pos__];

				_Iterator it = node.begin(), ite = node.end();
				for (; it!=ite; ++it)
				{
					const _LJAtom &atom = *it;

					_Iterator it__ = node__.begin(), ite__ = node__.end();
					if (pos__ == pos) { it__ = it; if (it__!=ite__) ++it__; }
					for (; it__!=ite__; ++it__)
					{
						const _LJAtom &atom__ = *it__;

						vector_t R
						(
							atom__.X[0] - atom.X[0],
							atom__.X[1] - atom.X[1],
							atom__.X[2] - atom.X[2]
						);

						region_->make_nearest_image_vector(R[0], R[1], R[2]);
						real_t r2 = R.length2();

						if (r2 > rmax2) continue;
					#ifndef USE_RARE_GAS
						if (has_contact(atom.insert_data, atom__.insert_data, atom.connect_data, atom__.connect_data))
							continue;
					#endif

						real_t sigma2 = calculate_sigma(atom.sigma, atom__.sigma);

						if (r2 < sqr(DEBUG_PAIR_FACTOR) * sigma2 || r2 < sqr(MINIMAL_DISTANCE_BETWEEN_ATOMS))
						{
							_S msg = _S("[WARNING] It has been found extra small contact [" )
								+ make_string(sqrt(r2)) + _S(" vs ") + make_string(sqrt(sigma2))
								+ _S("] between atoms: \n")
								+ make_string(atoms_[atom.insert_data], *atoms_[atom.insert_data].atomdata) + _S("\n")
								+ make_string(atoms_[atom__.insert_data], *atoms_[atom__.insert_data].atomdata) + _S("\n");
							PRINT_MESSAGE(msg);
							msg = _S("ins:") + make_string(atom.insert_data) + _S(" ") + make_string(atom__.insert_data) + _S("\n");
							msg += _S("con:") + make_string(atom.connect_data) + _S(" ") + make_string(atom__.connect_data);
							PRINT_MESSAGE(msg);
						}

						real_t eps = atom.eps * atom__.eps;
						real_t charge = atom.charge * atom__.charge;
						real_t tau2 = r2 / sigma2;
						real_t r_1 = (real_t) (1. / sqrt(r2));

						coul_energy__ = _Interaction::U(_I2T<COUL_>(), sigma2, charge, r2, tau2, r_1);
						vdw_energy__  = _Interaction::U(_I2T<VDW_ >(), sigma2, eps, r2, tau2, r_1);

						coul_energy += coul_energy__;
						vdw_energy += vdw_energy__;

						if (make_print)
						{
							max_coul_energy = (fabs(max_coul_energy) < fabs(coul_energy__)) ?
								coul_energy__ : max_coul_energy;
							max_vdw_energy = (fabs(max_vdw_energy) < fabs(vdw_energy__)) ?
								vdw_energy__ : max_vdw_energy;
							pairs_count++;
						}
					}
				}
			}
		}

		if (make_print)
		{
		#ifndef SKIP_COUL_ENERGY
			{
				real_t coul_evarage = (real_t) (pairs_count ? coul_energy / pairs_count : 0.);
				std::string msg
					= make_string("  U/coul   / (%5d) : %12.5le", pairs_count, (float)coul_energy)
					+ make_string("   <U> : %12.5le", (float)coul_evarage)
					+ make_string("   <maxU> : %12.5le", (float)max_coul_energy);

				PRINT_MESSAGE(msg);
			}
		#endif

		#ifndef SKIP_VDW_ENERGY
			{
				real_t vdw_evarage = (real_t) (pairs_count ? vdw_energy / pairs_count : 0.);
				std::string msg
					= make_string("  U/vdw    / (%5d) : %12.5e", pairs_count, (float)vdw_energy)
					+ make_string("   <U> : %12.5e", (float)vdw_evarage)
					+ make_string("   <maxU> : %12.5e", (float)max_vdw_energy);
				PRINT_MESSAGE(msg);
			}
		#endif
		}
	#endif // !defined(SKIP_COUL_ENERGY) || !defined(SKIP_VDW_ENERGY)

		return coul_energy + vdw_energy;
	}

	#undef TEMPLATE_HEADER
	#undef TEMPLATE_ARG

	typedef Complex_<FORCEFIELD_AMBER_, RESIDOME_AMBER_> Complex;

}
#endif
