#ifndef MOLDEFS__F9ED1116_3790_5da1_93E4_D94353B00A01__H
#define MOLDEFS__F9ED1116_3790_5da1_93E4_D94353B00A01__H

/*==============================================================================
 *                 ОПИСАНИЕ ОСОБЕННОСТЕЙ РЕАЛИЗАЦИИ ПРОЕКТА
 *
 *  Везде, где это возможно в проекте используется минимальный по памяти выбор
 *  типа. Это делается не с целью экономии памяти (благо ее достаточно), а по
 *  приине того, что такие типы предельно оптимизированы по производительности
 *  и операции с такими типами могут выполняться векторным образом (SSE для
 *  x86 процессоров, SIMD для Cell процессора, и т.д.
 *
 *  real_t - стандарный float тип, который используется для хранения параметров
 *           баз данных и в любых вычислениях, которые не требуют накопления
 *  extended<real_t> - расширенный тип, который необходим для любых переменных,
 *           значения который получаются путем накопления по объему выборки
 *           не превышающей 10**6 элементов, например для расчета ближних сил,
 *           действующих на атом.
 *============================================================================*/
/// использвание одинарной точночти
#define USE_FLOAT

/// использование группирования данных в малые вектора (по 4 элемента) и счет группами
#define USE_VECTORIZATION

/// определяется при использовании SSE(,2,3,4) расширений
#define USE_SSE

/// определяются для пропуска счета тех или иных энергетических вкладов
//  #define SKIP_BONDS_ENERGY
//  #define SKIP_ANGLS_ENERGY
//  #define SKIP_TORAS_ENERGY
//  #define SKIP_COUL_ENERGY
//  #define SKIP_VDW_ENERGY

/// использование кода для одноатомных незаряженных молекул инертных газов
//#define USE_RARE_GAS

/// определяется в случае использования Верлет таблицы (малые rcutoff и менее 10^6 атомов)
//#define USE_VERLET_TABLE
//#define USE_VERLET_TABLE_TUNE
//#define USE_VERLET_TABLE_GROUP

/// используется эффективно работающий в параллельном коде метод связанных ячеек
/// с упорядочением атомов в ячейках по направлениям по Гоннету (USE_GONNET)
#define USE_GONNET

//#define USE_ERF

/// определяется при использовании библиотеки (необходима для расчета дальнего кулона)
//#define USE_FFTW_LIBRARY

/// определяются в случае тестирования энергетических вкладов
#ifdef USE_RARE_GAS
	#define SKIP_BONDS_ENERGY
	#define SKIP_ANGLS_ENERGY
	#define SKIP_TORAS_ENERGY
	#define SKIP_COUL_ENERGY
#endif

/// определяется при использовании фиксированных в пространстве атомов
//#define USE_FIXED_ATOMS

/// определяется при использовании обнуления несвязующих потенциалов в нуле
//#define USE_PSEUDO_POTENTIAL_IN_ZERO

/// Определяется, если мы хотим гармоническое, а не арифметическое усреднение
/// для параметров vdw потенциала. Гармоническое усреднение позволяет ввести
/// vdw потенциалы.
//#define USE_HARM_AVERAGE_OF_SIGMA

/// Определяется при использовании неявной воды, что не отменяет использование
/// молекул воды индивидуально, например, в активном центре.
//#define USE_IMPLICIT_WATER

/// определяется при тестировании этапов построения молекулы
// #define BUILD_DEBUG

/// определяется при тестировании соответствия требуемых молекулой параметров
/// силовых полей тем, которые имеются в базе силового поля
// #define FULL_TOPOLOGY_DEBUG

/// определяется при тестировании Верлет таблицы
 #define VERLET_DEBUG

/// (всегда) определяется для тестировании построения молекулы
#define STOP_AFTER_FIRST_TOPOLOGY_ERROR

/// определяется при тестировании времени расчета энергетических вкладов
// #define ENERGY_TIME_TESTING

/// определяется при тестировании времени построения внутренних структур данных
// #define BUILDING_TIME_TESTING

/// (всегда) определяется для устойчивости выполнения молекулярной динамики
#define USE_HEAVY_HYDROGENS

/// определяется контроля устойчивости выполнения молекулярной динамики
#define USE_STATISTICS_CONTROL

#include "prgkern/__prgkern.h"

namespace molkern
{
	using namespace prgkern;
	using prgkern::make_string;
		// использовать все ранние определения make_string в пространстве имен prgkern
		// наряду с новыми в пространстве имен molkern

#ifdef USE_SSE
	#define USE_FLOAT
	#define USE_VECTORIZATION
#endif

#if defined(USE_VERLET_TABLE_TUNE) || defined(USE_VERLET_TABLE_GROUP)
	#define USE_VERLET_TABLE
#endif

#ifdef USE_FLOAT
	typedef float real_t;
#else
	typedef double real_t;
#endif

#ifdef USE_VECTORIZATION
	typedef vecreal_<4, real_t>  vreal_t;
	typedef vecint_ <4, int>     vint_t;
	typedef vecbool_<4, int>     vbool_t;
#else
	typedef vecreal_<1, real_t>  vreal_t;
	typedef vecint_ <1, int>     vint_t;
	typedef vecbool_<1, int>     vbool_t;
#endif

	/// тип координаты 3D пространства
	typedef vdense_<3, real_t> vector_t;

	/// тип любого пространственного ящика в 3D пространстве
	typedef prgkern::Box_<3, real_t> box_t;


	/// тип любого атомного идентификатора
	typedef unsigned unsigned_t;

	typedef signed char     signed8_t;
	typedef signed short    signed16_t;
	typedef signed int      signed32_t;
	typedef unsigned char   unsigned8_t;
	typedef unsigned short  unsigned16_t;
	typedef unsigned int    unsigned32_t;

	const real_t ACCURACY = 10 * std::numeric_limits<real_t>::epsilon();
	const real_t ACCURACY2 = sqr(ACCURACY);

	/// using identifiers of types & objects
	enum IDENTS_OF_OBJECTS
	{
		NUCLEAR_,            // заряд ядра
		HYDROGEN_,           // водород
		ATOMDATA_,           // параметры атома
		ATOM_,               // атом
		BOND_,               // связь
		ANGLE_,              // валенный угол
		TORSION_,            // торсионный угол
		ROTAMER_,            // ротамер
		EQUI_,               // тип эквивалентности
		RESIDUE_,            // аминокислотный остаток
		CHAIN_,              // цепь свзанных атомов
		NTERM_,              // N терминал
		CTERM_,              // C терминал
		CONNECT_,            // контакт

		FORCEFIELD_,         // силовое поле
		RESIDOME_,           // база данных топологии
		ARCHETYPE_,          // архетип молекулы
		MOLECULE_,           // молекула
		WATER_,              // раствор

		PAIR_,               // парное взаимодействие
		COUL_,               // кулоновское взаимодействие
		VDW_,                // ван-дер-ваальсово взаимодействие
		PAIR14_,             // парное 1-4 взаимодействие
		COUL14_,             // кулоновское 1-4 взаимодействие
		VDW14_,              // ван-дер-ваальсово 1-4 взаимодействие
		FREE_,               // свободные атомы

		NODE_,               // узел данных (или копьютера)
		DIRECTION_,          // типы упорядочиваний по направлениям
		EMPTY_DIRECTION_,
		NODE_PAIR_,
		SET_,
		REBUILD_,

		EXT_ID_,             // внешний идентификатор атома
		FILE_,               // файл
		NAME_,
		EQTY_DISTANCE_,      // distance between types in table (we use the nearest)

		RESIDUE_NAME_,       // building internal residue names
		RESIDUE_CONTACT_,    // building internal residue contacts
		PDB_RESIDUE_NAME_,   // standard PDB residue name (HIS instead HIE or HID etc)

		SASA_,               // solvent accessible surface area
		MAX_VDW_RADIUS_,     // max forcefield atoms radius
		CHARGE1_,            // расчет изменений заряда

		MASS_CENTER_,        // center of mass
		GEOM_CENTER_,        // center of geometry
		INERTIA_TENSOR_DET_, // determinant of inertia tensor
		SOLID_TdS_,          // entropy of elimination of center mass & rotation as whole
		TdS_,                // полная энтропия
		A_DET_,

		EXT_CONTACT_PREV_,
		EXT_CONTACT_NEXT_,

                DCHARGE_,            // du__dq
		CHARGE_,             // заряд
		MASS_,               // масса
		TEMPERATURE_,        // температура ансамбля (K)
		PRESSURE_,           // давление
		VOLUME_,             // объем
		VELOCITY_,           // скорость

		FREEDOM_TYPE_,
		FREEDOM_FIXED_,      // нет свободы движения
		FREEDOM_CM_,         // свобода движения центра масс (все цепи молекулы синхронны)
		FREEDOM_CHAIN_,      // свобода движения подцепей как целого (цепи не синхронны)
		FREEDOM_ROTAMER_,    // свобода вращения ротамеров
		FREEDOM_ATOM_,        // свобода движения всех атомов
		FREEDOM_,

		START_POSITION_,
		POSITION_,
		XPOSITION_,
		GPOSITION_,

		BOX_,                // прямоугольная область пространства
		REGION_,             // область пространства модели

		OPTIMIZER_,
		THERMOSTATE_,

		/// matrix types for special purpouses
		ACCUMULATE_,         // ACCUMULATE_ matrix type

		STICK_,             // ось ротамера, по которой ротамеры связываются
		ROTAMER_EXT_ATOM_,  // атом на стике, но вне ротамера
		ROTAMER_INT_ATOM_,  // атом на стике, принадлежащий ротамеру
		ROTAMER_ADD_ATOM_,  // атом для расчетов относительных углов
		ROOT_ROTAMER_,      // ротамеры, являющиеся стартовыми в построении
		EDGE_,              // ротамерные оси
		ROTAMER_POSITION_,
		MOMENT_,
		GRADIENT_,

		FORCEFIELD_AMBER_,    // AMBER protein forcefield
		FORCEFIELD_GAFF_,     // AMBER inorganic ligand forcefield
		RESIDOME_AMBER_,      // AMBER protein residome
		FORMAT_UNKNOWN_,
		FORMAT_PDB_,          // PDB file format
		FORMAT_HIN_,          // HIN file format
		FORMAT_MOL2_,         // MOL2 file format
		FORMAT_BMM_,          // формат пакета MOLKERN
		FORMAT_VHIN_,          // формат пакета MOLKERN

		VERLET_PAIR_,
		RSKIN_,
		IMPULSE_,
		PARALL_,
		NEIGHBOR_,
		EDGE_LENGTH_,

		OBJECT_UNKNOWN_
	};

#define DEFINE_TAG_CONST(name)  const _I2T<name##_> name = _I2T<name##_>()

	//--------------------------------------------------------------------------
	//                   Различные теговые константы
	//--------------------------------------------------------------------------
	DEFINE_TAG_CONST(HYDROGEN);
	DEFINE_TAG_CONST(ATOMDATA);
	DEFINE_TAG_CONST(ATOM);
	DEFINE_TAG_CONST(BOND);
	DEFINE_TAG_CONST(ANGLE);
	DEFINE_TAG_CONST(TORSION);
	DEFINE_TAG_CONST(PAIR14);
	DEFINE_TAG_CONST(PAIR);
	DEFINE_TAG_CONST(CONNECT);
	DEFINE_TAG_CONST(ROTAMER);
	DEFINE_TAG_CONST(RESIDUE);
	DEFINE_TAG_CONST(CHAIN);
	DEFINE_TAG_CONST(EDGE);
	DEFINE_TAG_CONST(MOMENT);
	DEFINE_TAG_CONST(WATER);
	DEFINE_TAG_CONST(RESIDUE_NAME);
	DEFINE_TAG_CONST(RESIDUE_CONTACT);
	DEFINE_TAG_CONST(MOLECULE);

	DEFINE_TAG_CONST(COUL);
	DEFINE_TAG_CONST(VDW);

	DEFINE_TAG_CONST(FREEDOM_TYPE);
	DEFINE_TAG_CONST(FREEDOM_FIXED);
	DEFINE_TAG_CONST(FREEDOM_CM);
	DEFINE_TAG_CONST(FREEDOM_CHAIN);
	DEFINE_TAG_CONST(FREEDOM_ROTAMER);
	DEFINE_TAG_CONST(FREEDOM_ATOM);
	DEFINE_TAG_CONST(FREEDOM);

	DEFINE_TAG_CONST(REBUILD);
	DEFINE_TAG_CONST(BOX);
	DEFINE_TAG_CONST(REGION);
	DEFINE_TAG_CONST(NODE);
	DEFINE_TAG_CONST(SET);
	DEFINE_TAG_CONST(NODE_PAIR);
	DEFINE_TAG_CONST(DIRECTION);
	DEFINE_TAG_CONST(EMPTY_DIRECTION);
	DEFINE_TAG_CONST(MASS_CENTER);
	DEFINE_TAG_CONST(GEOM_CENTER);
	DEFINE_TAG_CONST(EULER_ROTATOR);
	DEFINE_TAG_CONST(XYZ_ROTATOR);
	DEFINE_TAG_CONST(MAX_VDW_RADIUS);
	DEFINE_TAG_CONST(STICK);
	DEFINE_TAG_CONST(ROTAMER_EXT_ATOM);
	DEFINE_TAG_CONST(ROTAMER_INT_ATOM);
	DEFINE_TAG_CONST(ROTAMER_ADD_ATOM);
	DEFINE_TAG_CONST(ROTAMER_POSITION);
	DEFINE_TAG_CONST(ROOT_ROTAMER);

	DEFINE_TAG_CONST(START_POSITION);
	DEFINE_TAG_CONST(CHARGE);
	DEFINE_TAG_CONST(DCHARGE);
	DEFINE_TAG_CONST(POSITION);
	DEFINE_TAG_CONST(XPOSITION);
	DEFINE_TAG_CONST(GPOSITION);
	DEFINE_TAG_CONST(GRADIENT);

	DEFINE_TAG_CONST(FORCEFIELD_AMBER);
	DEFINE_TAG_CONST(FORCEFIELD_GAFF);
	DEFINE_TAG_CONST(RESIDOME_AMBER);
	DEFINE_TAG_CONST(FORMAT_UNKNOWN);
	DEFINE_TAG_CONST(FORMAT_PDB);
	DEFINE_TAG_CONST(FORMAT_VHIN);
	DEFINE_TAG_CONST(FORMAT_HIN);
	DEFINE_TAG_CONST(FORMAT_MOL2);
	DEFINE_TAG_CONST(FORMAT_BMM);

	DEFINE_TAG_CONST(VERLET_PAIR);
	DEFINE_TAG_CONST(RSKIN);

	DEFINE_TAG_CONST(TEMPERATURE);
	DEFINE_TAG_CONST(PRESSURE);
	DEFINE_TAG_CONST(IMPULSE);
	DEFINE_TAG_CONST(PARALL);
	DEFINE_TAG_CONST(NEIGHBOR);
	DEFINE_TAG_CONST(EDGE_LENGTH);

#undef DEFINE_TAG_CONST

	const bool YES_HYDROGEN = true; // print hydrogen atoms to file
	const bool NO_HYDROGEN = false; // don't print hydrogen atoms to file
		// ignored for *.hin files

	const unsigned BUILD_NOTHING       = 0x00000000;
	const unsigned BUILD_HYDROGENS     = 0x00000001;
	const unsigned BUILD_BONDS         = 0x00000002;
	const unsigned BUILD_ANGLES        = 0x00000004;
	const unsigned BUILD_TORSIONS      = 0x00000008;
	const unsigned BUILD_PAIRS14       = 0x00000010;
	const unsigned BUILD_CONNECTS      = 0x00000020;
	const unsigned BUILD_ROTAMERS      = 0x00000040;
	const unsigned BUILD_CHAINS        = 0x00000080;

	const unsigned BUILD_ALL = BUILD_HYDROGENS
		| BUILD_BONDS | BUILD_ANGLES | BUILD_TORSIONS | BUILD_PAIRS14 | BUILD_CONNECTS
		| BUILD_ROTAMERS | BUILD_CHAINS;

	const unsigned NO_UNION_    = 0x00000000;
	const unsigned NO_CM_       = 0x00000000;
	const unsigned NO_ROTAMER_  = 0x00000000;
	const unsigned NO_ATOM_     = 0x00000000;

	const unsigned YES_UNION_    = 0x00000001; // бит объединения подцепей в движении
	const unsigned YES_CM_       = 0x00000002;
	const unsigned YES_ROTAMER_  = 0x00000004;
	const unsigned YES_ATOM_     = 0x00000008;

	const unsigned CALC_NOTHING_ = 0x00000000;
	const unsigned CALC_BOND_    = 0x00000001;
	const unsigned CALC_ANGLE_   = 0x00000002;
	const unsigned CALC_TORSION_ = 0x00000004;
	const unsigned CALC_PAIR14_  = 0x00000008;

	const bool YES_HETEROATOMS = true; // read heteroatoms from file
	const bool NO_HETEROATOMS = false; // don't read heteroatoms from file
		// ignored for *.hin files

	const unsigned_t USE_AS_FREE_  = 0x00000000; // нет битов
	const unsigned_t USE_AS_FIXED_ = 0x80000000; // бит в 31 позиции

	const bool YES_PRINT = true;
	const bool NO_PRINT = false;

	// длина заголовка в байтах в bmm формате
	const unsigned BMM_HEADER_LEN = 128;
	// длина одной записи в байтах в bmm формате
	const unsigned BMM_RECORD_LEN = 18;

	// точность данных в bmm формате (accuracy = 1A / BMM_GRID_STEP)
	const unsigned BMM_XGRID_STEP = 100;
	// точность хранения скоростей совпадает с точностью координат
	const unsigned BMM_VGRID_STEP = 100;
	// точность хранения заряда на 2 порядка выше точности хранения координат.
	const unsigned BMM_QGRID_STEP = 10000;

	/// error codes (file extractors, etc)
	enum code_error_
	{
		CODE_SUCCESS = 0, // no error
		CODE_ERROR,       // unknown error
		CODE_EOF          // end of file
	};

	const real_t INFINITY = infinity<real_t>();

	/// U() выводит все связи, длина которых превышает норму в DEBUG_BOND_FACTOR раз
	const real_t DEBUG_BOND_FACTOR = 2.;

	/// U() выводит все пары с расстоянием менее нормы в 1/DEBUG_PAIR_FACTOR раз
	const real_t DEBUG_PAIR_FACTOR = 0.5;

	const unsigned MAX_GROUP_ELEM = 32;

	/// lenght of atom force field type name
	const int MAX_ATOM_NAME   = 4; // 2 + reserve + '/n' (int)
		// must be fixed as 4 everytimes (used special code for fast)

	const int MAX_NAME       = 256; // lenght of name (molecule, ..)
	const int MAX_EQUIVALENT = 20; // max number of atom equivalent
		// if atom has more than 20 equivalents
		// than creates a new atom equivalence group with the same org

	/// max bond of atom
	const int MAX_ATOM_BOND = 6;

	/// количество атомов в среднем расчете
	const int DEFAULT_ATOM_COUNT = 4000;

	/// область контроля для удаление атомов воды [A]
	const real_t WATER_CLASH_RADIUS = 5.; // 3.3950 + 1.7683; (Cs + OW)
		// подгонка под популярный размер радиуса взаимодействия

	/// hydrogen - atom distance [angstrom]
	const real_t HYDROGEN_ATOM_DISTANCE = 1.09;

	const real_t MINIMAL_DISTANCE_BETWEEN_ATOMS = HYDROGEN_ATOM_DISTANCE - 0.2;

	/// коэффициент для увеличения массы водородов при USE_HEAVY_HYDROGENS
	const real_t HEAVY_MASS_FACTOR = 2;

	/// длина S-S связи (сульфидного мостика)
	const real_t SS_BOND_LENGTH = 2.2;

	/// default charge mesh step [angstrom]
	const real_t CHARGE_MESH_STEP = 1.;

	/// default charge radius of smoothing [angstrom]
	const real_t CHARGE_SMOOTH_RADIUS = 2.;

	/// число попыток заполнения вставки каждой молекулы в случайную позицию области
	const unsigned DEFAULT_DISPOSE_ITERATIONS = 10;

	/// среднее количество атомов в единице объема в 1 A**3 равно 0.0036
	/// получено через тестовые расчеты и чуть увеличено
	const real_t DEFAULT_ATOM_DENSITY = 0.05;

	/// default step of trajectory integration [ps]
	const unsigned DEFAULT_INTEGRATION_TIME_STEP = 1; // 1 fs
	const unsigned DEFAULT_SAMPLING_TIME_STEP = 10; // 10 fs

	/// coeffs of scaling 14-interactions
	const real_t AMBER_SCALE_COUL14 = 0.8333; // opls.pdf
	const real_t AMBER_SCALE_VDW14 = 0.5;

	const real_t DEFAULT_COUL_SPLIT_RADIUS = 8.0f;

	/*----------------------------------------------------------------------------
	*                     Используемые единицы измерения
	* ----------------------------------------------------------------------------
	* длина [Angstrom] 10**(-10) м
	* время [ps] 10**(-12) сек
	* масса [a.u.m.] 1/12 массы углерода
	* энергия [AUE/mol] атомная единица энергии (она ровно в 100 больше kJ/mol)
	* --------------------------------------------------------------------------*/
	/// функция печатает единицы измерения программы
	inline void print_user_information()
	{
		_S msg = _S("\n"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
			"!!                                                                     !!\n"
			"!!                           BMMKERN                                   !!\n"
			"!!                                                                     !!\n"
			"!!            Biological Molecular Modelling package (kern)            !!\n"
			"!!                                                                     !!\n"
			"!!                     version 2010.08                                 !!\n"
			"!!                                                                     !!\n"
			"!!                                                                     !!\n"
			"!!   Copyright 2006, Institute of Cytology and Genetics,               !!\n"
			"!!   Eduard S.Fomin  fomin@bionet.nsc.ru                               !!\n"
			"!!                                                                     !!\n"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
			"!!                                                                     !!\n"
			"!!  MEASURE UNITS        :  A(length), ps(time), a.m.u.(mass),         !!\n"
			"!!                          e(charge), K(temperature)                  !!\n"
			"!!  ENERGY UNIT (a.e.u)  :  a.m.u * (A / ps) ** 2 = 0.01 kJ /mol       !!\n"
			"!!                                                                     !!\n"
			"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		msg += _S("\n  Preprocessor directives used :");

	#ifdef USE_FLOAT
		msg += _S("\n    USE_FLOAT");
	#endif
	#ifdef USE_VECTORIZATION
		msg += _S("\n    USE_VECTORIZATION");
	#endif
	#ifdef USE_SSE
		msg += _S("\n    USE_SSE");
	#endif

	#ifdef USE_VERLET_TABLE
		msg += _S("\n    USE_VERLET_TABLE");
		#ifdef USE_VERLET_TABLE_TUNE
			msg += _S(" + USE_VERLET_TABLE_TUNE");
		#endif
		#ifdef USE_VERLET_TABLE_GROUP
			msg += _S(" + USE_VERLET_TABLE_GROUP");
		#endif
	#endif

	#ifdef USE_GONNET
		msg += _S("\n    USE_GONNET");
	#endif

	#ifdef USE_RARE_GAS
		msg += _S("\n    USE_RARE_GAS");
	#endif
	#ifdef SKIP_BONDS_ENERGY
		msg += _S("\n    SKIP_BONDS_ENERGY");
	#endif
	#ifdef SKIP_ANGLS_ENERGY
		msg += _S("\n    SKIP_ANGLS_ENERGY");
	#endif
	#ifdef SKIP_TORAS_ENERGY
		msg += _S("\n    SKIP_TORAS_ENERGY");
	#endif
	#ifdef SKIP_COUL_ENERGY
		msg += _S("\n    SKIP_COUL_ENERGY");
	#endif
	#ifdef SKIP_VDW_ENERGY
		msg += _S("\n    SKIP_VDW_ENERGY");
	#endif

	#ifdef USE_FIXED_ATOMS
		msg += _S("\n    USE_FIXED_ATOMS");
	#endif
	#ifdef USE_PSEUDO_POTENTIAL_IN_ZERO
		msg += _S("\n    USE_PSEUDO_POTENTIAL_IN_ZERO");
	#endif
	#ifdef USE_HARM_AVERAGE_OF_SIGMA
		msg += _S("\n    USE_HARM_AVERAGE_OF_SIGMA");
	#endif
	#ifdef USE_IMPLICIT_WATER
		msg += _S("\n    USE_IMPLICIT_WATER");
	#endif
	#ifdef USE_HEAVY_HYDROGENS
		msg += _S("\n    USE_HEAVY_HYDROGENS");
	#endif
	#ifdef USE_STATISTICS_CONTROL
		msg += _S("\n    USE_STATISTICS_CONTROL");
	#endif

		PRINT_MESSAGE(msg);
	}

	/// тип параметров, описывающих время модели
	struct model_time_t
	{
		unsigned tm_;
		model_time_t(unsigned t=0) : tm_(t) {}
		operator unsigned() { return tm_; }
		model_time_t &operator+=(unsigned dt) { tm_ += dt; return *this; }
	};

	INLINE model_time_t operator-(model_time_t t1, model_time_t t2)
	{ return model_time_t(t1.tm_ - t2.tm_); }

	extern model_time_t global_model_time; // текущее время для динамики (число fs)
		// все процедуры могут считывать это значение, например, для печати
		// но устаналивает время только одна функция (после завершения очередного
		// цикла динамики)
	extern time_t global_real_start_time; // физическое время начала динамики
	extern time_t global_real_prev_time;

	/// возвращает текущее время модели
	INLINE model_time_t current_model_time() { return global_model_time; }

	/// функция извлекает из глобального счетчика текущее время для динамики
	INLINE std::string make_string(model_time_t)
	{
		model_time_t tm = current_model_time();
		unsigned ns = tm / 1000000;
		unsigned ps = tm % 1000000 / 1000;
		unsigned fs = tm % 1000;
		return make_string("%3d,%03d'%03d ns", ns, ps, fs);
	}

	extern unsigned global_thread_count; //число потоков
	extern real_t global_rskin_width;
	extern real_t global_compress_factor;

	/// шаг изменения rskin при его оптимизации
	const double DEFAULT_RSKIN_STEP = 0.1;
	const double DEFAULT_RSKIN_WIDTH = 1.2;

	/// coefficient to transform cal -> J
	const real_t CAL2J = (real_t) 4.184; // US cal
	const real_t CAL2AUE = (real_t) 418.4; // US cal
		// AUE атомная единица энергии (а.е.м. * ангстрем**2 / пс**2)

	const real_t J2CAL = (real_t) (1. / CAL2J); // US cal

	const real_t AUE_BOLTZMAN = (real_t) 8.314510e-1; // AUE/(mol*K)

	template <typename T>
	INLINE T KT(T temperature) { return (T)(AUE_BOLTZMAN * temperature); }

	template <typename T>
	INLINE T Temperature(T energy) { return  (T)(energy / AUE_BOLTZMAN); }

	const bool NO_PERIODIC_ = false; // флаг отсутствия пространственной периодичности системы
	const bool YES_PERIODIC_ = true; // флаг пространственной периодичности системы

	const bool NO_CLASHES = false; // флаг запрета клеширования
	const bool YES_CLASHES = true; // флаг разрешения клеширования

	/// нормальная температура
	const real_t DEFAULT_TARGET_TEMPERATURE = (real_t) 300.; // 309.6 = 36.6 C

	/// максимально разрешенное смещение во время динамики
	const real_t DEFAULT_MAX_X = (real_t) 0.2; // [A]

	const real_t ELECTRIC_FACTOR = (real_t) 138935.485;
		// 1 / (4 pi e0) [AUE * A / (mol e^2)]
	const real_t SQRT_ELECTRIC_FACTOR = (real_t) sqrt(ELECTRIC_FACTOR);

	const real_t DIELECTRIC_WATER = (real_t) 78.5;
	const real_t DIELECTRIC_PROTEIN = (real_t) 4.0;

	/// коэффициент перевода из [kJ /mol A**3] в число атмосфер
	const real_t ATMOSPHERE_FACTOR = (real_t) 16388.5070058; //????

	/// the normal condition pH of water
	const int NORMAL_PH_WATER = 7;
	//const real_t NORMAL_WATER_DENSITY = 0.03345; // [молекул (!) в ангстрем **3]
	const real_t NORMAL_WATER_DENSITY = 0.10035; // [частиц (!) в ангстрем **3]

	/// standard radius of water molecule
	const real_t RADIUS_H2O = (real_t) 1.4; // Angstrom
		// http://nook.cs.ucdavis.edu:8080/~koehl/BioEbook/protsurf.html

	/// value of space integral for postevaluation calculation
	const real_t INTEGRAL__dXYZ = (real_t) 1660;

	/// value of euler integral for postevaluation calculation
	const real_t INTEGRAL__dOMEGA = (real_t) (8 * sqr(M_PI));

	/** @brief matrix for different purpouses
	* @param TYPE type of matrix
	*  matrix<ACCUMULATE_> is used to avoid dublications in topology building
	*/
	template <int TYPE> struct matrix;

	#define IMPLEMENT_MATRIX(NAME, ROW_TYPE) \
	template <> struct matrix<NAME> \
	: public std::vector<ROW_TYPE > \
	{ \
		typedef ROW_TYPE::iterator row_iterator; \
		typedef ROW_TYPE::const_iterator row_const_iterator; \
		typedef ROW_TYPE row_type; \
	};
	IMPLEMENT_MATRIX(ACCUMULATE_, std::set<int>);
	#undef IMPLEMENT_MATRIX

#define LOADED_OK_MESSAGE(filename) { \
		std::string msg = _S(filename) + _S(" has been loaded"); \
		PRINT_MESSAGE(msg); \
	}

	#define LOADED_ERR_MESSAGE(filename) { \
		std::string msg = _S(filename) + _S(" hasn't been loaded"); \
		PRINT_ERR(msg); \
	}

	#define SAVED_OK_MESSAGE(filename) { \
		std::string msg = _S(filename) + _S(" has been saved"); \
		PRINT_MESSAGE(msg); \
	}

	#define SKIP_LINE(src) { \
		src.ignore(MAX_STRING_LEN, '\n'); \
	}

	#define SKIP_UNTIL_EMPTY_LINE(src) { \
		char buffer___[MAX_STRING_LEN]; \
		do { \
			src.getline(buffer___, MAX_STRING_LEN); \
		} while (!isspace(buffer___)); \
	}

	#define NO_EQUIVALENTS(index__) { \
		std::string msg = _S("[ERROR] No equivalents for ") \
			+ make_string(index__); \
		PRINT_MESSAGE(msg); \
	}

	#ifdef FULL_TOPOLOGY_DEBUG
		#ifndef NDEBUG
		#define USES_EQUIVALENTS(ndx, index__) { \
			if (ndx != index__) { \
				std::string msg = _S("[WARNING] Database uses equivalent : ") \
					+ _S(make_string(ndx)) + _S(" -> ") \
					+ _S(make_string(index__)); \
				PRINT_MESSAGE(msg); \
			} \
		}
		#else
			#define USES_EQUIVALENTS(ndx, index__) {}
		#endif
	#else
		#define USES_EQUIVALENTS(ndx, index__)
	#endif

	#define DEFINE_VECTOR_ACCESS_FUNCTION(name, type, array) \
		unsigned count(_I2T<name>) const { return (unsigned)array.size(); } \
		\
		const type *get(_I2T<name>, unsigned i=0) const \
		{ assert(_LT((unsigned)i, array.size())); return &array[i]; } \
		\
		type *get(_I2T<name>, unsigned i=0) \
		{ assert(_LT((unsigned)i, array.size())); return &array[i]; } \
		\
		unsigned update(_I2T<name>, std::vector<type> &array__) \
		{ array.swap(array__); return array.size(); } \
		\
		void resize(_I2T<name>, unsigned n, type v=type()) { array.resize(n, v); } \
		\
		template <typename ITERATOR> void insert(_I2T<name>, ITERATOR from, ITERATOR to) \
		{ array.insert(array.end(), from, to); }

#define DEFINE_TAG_LOGICAL_CONST(name)  const _B2T<name##_> name = _B2T<name##_>()

	DEFINE_TAG_LOGICAL_CONST(YES_PERIODIC);
	DEFINE_TAG_LOGICAL_CONST(NO_PERIODIC);

#undef DEFINE_TAG_LOGICAL_CONST


	// таймеры для функций
	enum { READ_TIMER_=1000, SAVE_TIMER_, RUN_TIMER_, RUNV_TIMER_, NEAR_TIMER_, BUILD_TIMER_, CALC_TIMER_,
		GON_TIMER_, FILTER_TIMER_, BUILD2_TIMER_, TIMER_1, TIMER_2, TIMER_3, TIMER_4, TIMER_5, TIMER_6, TIMER_7, TIMER_8};

	extern unsigned wwww;
	extern long long unsigned global_pair;

	/**
	 * Объект позволяет унифицированно подключать любой расчетчик ближних взаимодействий
	 */
	template <typename _LPComplex>
	class near_range_integrator
	{
	public:

		virtual void resize(real_t rcut, real_t density) = 0;
		virtual void update(bool print=NO_PRINT) = 0;
		virtual _E(real_t) dU__dX() = 0;
	};
}

#endif
