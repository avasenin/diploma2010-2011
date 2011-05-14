#ifndef CONFIG__F9ED1116_3790_5da1_93E4_D94353B00A01__H
#define CONFIG__F9ED1116_3790_5da1_93E4_D94353B00A01__H

#include "prgkern/__prgkern.h"
#include "molkern/__moldefs.h"
#include <boost/program_options.hpp>

namespace molkern
{
	/*
	 * В программе используются регулярные выражения для анализа разных строк. Например, для
	 * извлечения данных из подстрок командной строки, типа -eNVT[300K]. Здесь дано некоторое
	 * описание синтаксиса некоторых конструкций регулярных выражений для быстрой ориентировки.
	 * Дело в том, что без непрерывного их использования, этот синтаксис стремительно забывается.
	 *
	 *    . любой символ
	 *    * любое число повторений предыдущего символа
	 *    *? "нежадный" поиск любого числа повторений предыдущего символа
	 *    \. символ точка
	 *    \\ вынужденный двойной слеш, по требованию языка С++
	 *    () подвыражение, которое запоминается
	 *    (?: ) незапоминаемое подвыражение, позволяет делать незапоминаемые конструкции выбора
	 *    ( / / ) незапоминаемая конструкция выбора
	 *    [] любой символ внутри скобок
	 *    [0-9] любая цифра внутри диапазона
	 *    \\[ квадратная скобка (двойной слеш по требованию С++)
	 *    ^ привязка к началу строки
	 *    $ привязка к концу строки
	 */

	/*
	 * (1) Формат строки "-sstatistics.bmm[0.1ns],wh"
	 */
	const RegEx FILE_NAME_REGEX     = RegEx(_S("(.*?\\.(?:pdb|ent|hin|mol2|bmm)).*"));
	const RegEx PRN_HYDROGEN_REGEX  = RegEx(_S(".*,.*?(h).*?$"));
	const RegEx PRN_WATER_REGEX     = RegEx(_S(".*,.*?(w).*?$"));
	const RegEx TIME_REGEX          = RegEx(_S(".*?\\[([0-9.]*)(sec|ms|mks|ns|ps|fs)\\].*"));

	/*
	 * (1) Формат строки -eNVT[300K], -eNVP[1.0atm], -eNPT[1atm,300K], -eNVE[3kJ/mol]
	 */
	const RegEx ENSEMBLE_REGEX      = RegEx(_S("^(N..).*"));
	const RegEx TEMPERATURE_REGEX   = RegEx(_S(".*?\\[.*?([0-9.eE+-]*?)(K|C).*?\\].*"));
	const RegEx PRESSURE_REGEX      = RegEx(_S(".*?\\[.*?([0-9.eE+-]*?)(atm).*?\\].*"));
	const RegEx ENERGY_REGEX        = RegEx(_S(".*?\\[.*?([0-9.eE+-]*?)(kJ/mol).*?\\].*"));

	/*
	 * (1) Формат строки "-mligand.hin[1.2nmol],hw{umrx}"
	 * (2) Формат строки "-mprotein.pdb[1.2nmol],Ahwl{umrx}"
	 * (3) Формат строки "-wSPCBOX[0.1atoms/A3],{umrx}"
	 */
	const RegEx WATER_REGEX         = RegEx(_S("^(.*?)(:?,|\\[)"));
	const RegEx CONCENTRATION_REGEX = RegEx(_S(".*?\\[.*?([0-9.eE+-]*?)"
		"(mol|mmol|mkmol|nmol|pmol|fmol).*?\\].*"));
		// Концентрация описывается в числе молекул лиганда на число молекул воды (молей лиганда на моль
		// воды), например, концентрация 1nmol соответствует 1 молекуле лиганда на 10**9 молекул воды.

	const RegEx DENSITY_REGEX	      = RegEx(_S(".*?\\[.*?([0-9.eE+-]*?)"
		"(molecule/A3|molecule/nm3).*?\\].*"));
		// плотность описывается как число молекул на заданный объем (A**3 или nm**3)

	const RegEx LIGAND_REGEX       = RegEx(_S(".*,.*?(l).*?$"));
		// вытаскивает из описателя файла молекулы символ 'l', который означает, что из файла
		// нужно загрузить не только основную молекулу, но и лиганды, если они есть. Это важно для
		// pdb файлов, которые могут включат дополнительные кофакторы, лиганды, атомы воды...

	const RegEx FREEDOM_REGEX      = RegEx(_S(".*,.*?([umrx]*).*?$"));
	const RegEx ALTPOS_REGEX       = RegEx(_S(".*,.*?([ABC]).*?$"));

	template <int TYPE> struct Interface_;

	/**
	 * @param time время в единицах, записанных в unit
	 * @param unit единица измерения времени
	 * @return время в [ps]
	 */
	INLINE real_t make_time(real_t time, const std::string &unit)
	{
		if      (unit == _S("fs" )) time *= 1e-3;
		else if (unit == _S("ps" )) time *= 1;
		else if (unit == _S("ns" )) time *= 1e3;
		else if (unit == _S("mks")) time *= 1e6;
		else if (unit == _S("ms" )) time *= 1e9;
		else if (unit == _S("sec")) time *= 1e12;
		return time;
	}

	/**
	 * Класс, хранящий все описания процесса симуляции. Объект этого класса является глобальным,
	 * и к нему имеют доступ все объекты программы.
	 */
	class Configure
	{
		boost::program_options::options_description desc_;
		boost::program_options::variables_map vm_;

	public:

		int threads; // число потоков выполнения
			// Число нитей должно согласовываться с числом ячеек в регионе (должно быть делителем
			// полного числа ячеек региона). Так как принято, что число ячеек в регионе должно быть
			// четным по всем направлениям, то максимальный делитель уже не менее 8.

		_S work_dir;   // директория, откуда забираются файлы молекул
		_S data_dir;   // директория, откуда забираются параметры силовых полей
		_S config;     // файл параметров конфигурации расчета

		_S box;        // Глобальный ящик, который задается строкой [nx,ny,nz]
			// в которой задается число ячеек ящика вдоль координат {x,y,z}.
			// Размер минимальной ячейки всегда равен радиусу взаимодействия.

		std::vector<_S> archetypes;
			// (1) Формат строки "-mligand.hin[N],Ahwl{umrx}", где :
			//    [N] - число копий молекул,
			//    A - альтернативная позиция атомов (работает только для *.pdb).
			// Наличие/отсутствие символов h и w означает использование/игнорирование
			// атомов водорода (h), молекул воды(w), лиганда(l) внутри файла
			//
			//  Тип движения молекул заданого архетипа.
			//    u   - объединие движения подцепей
			//    m   - движение центров масс
			//    r   - движение ротамеров
			//    x   - движение атомов
			//  Типы движения комбинируются произвольно, давая множество вариантов. Большие буквы (или
			//  отсутствие соответствующих малых букв) определяют запрет на движение заданных элементов.
			//    u     (пустой список) фиксация в пространстве
			//    rm    движение ротамеров и центров масс подцепей
			//    urm   движение ротамеров и центра масс всего архетипа
			//    x     движение атомов при фиксации центров масс подцепей
			//    xm    наиболее свободное движение
			//    xrm   движение атомов в рутовом ротамере + ротамеры + центр масс (гибридные координаты)
			//
			// При отутствии числа копий считается, что число молекул 1.
			// При отсутствии альтернативной позиций предполагается позиция A.

		_S outfile;
			// Формат строки "-oligand.hin,hw", где :
			// наличие/отсутствие символов h и w означает сохранение/игнорирование
			// атомов водорода (h) и молекул воды(w) при записи в файл
			// Заметим, что для HIN, MOL2 и BMM файлов молекулы водородов сохраняются
			// всегда.

		_S statfile;
			// Формат строки "-sligand.bmm[0.1ns],wh", где :
			// наличие/отсутствие символов w означает сохранение/игнорирование
			// молекул воды(w) и водородов(h) при записи в файл

		_S water; // Формат строки "-wSPCBOX,umcxr"
		int pH; // pH, к которому приводятся все молекулы комплекса, может быть без воды

		real_t cutoff_radius; // радиус обрезания потенциала
		real_t barrier; // высота барьера в 0 для 6-12 псевдопотенциала
		real_t mesh_fftw_step;
		real_t rskin_width;
		real_t rskin_compress_factor;

		/*
		 *________________________ ОПТИМИЗАЦИОННЫЕ ПАРАМЕТРЫ ______________________________
		 */
		_S optimizer;            // имя оптимизатора
		int iterations;          // число итераций оптимизации геометрии
		int steep;               // число начальных итераций, выполняемых steep методом
		real_t stpmin;           // ограничение по минимальному смещению по координатам
		real_t stpmax;           // ограничение по максимальному смещению по координатам
		real_t xtol;             // минимальное значение изменения координат, обрывающее оптимизацию
		real_t ftol;             // минимальное значение изменения энергии, обрывающее оптимизацию
		real_t gtol;             // минимальное значение изменения силы, обрывающее оптимизацию
		int maxfev;              // максимально разрешенное число оценок энергии, обрывающее оптимизацию
			// Данный параметр используется, чтобы избежать зацикливания расчета в силу любых причин
			// Такими причинами могут быть ошибки кода, ошибки эмпирических правил поиска минимуми,
			// расхождение по точности между переданными энергиями и ее производными и т.д.

		/*
		 * Параметры Вольфа используются, чтобы ускорить процесс нахождения минимума. Они основываются
		 * на опыте, что нахождение локального минимума вдоль некоторого направления с большой точностью
		 * является бессмысленной задачей, а расходы на это достаточно велики. По этой причине поиск
		 * локального минимума делается приблизительно, по некоторым эмпирическим правилам. Правила
		 * таковы, что если в некоторй промежуточной точке значение производной по абсолютному
		 * значению упало до указанной доли от исходной величины, то процесс можно прерывать.
		 * Что же касается самой функции, то данная доля используется на основе прогноза по значению
		 * функции в возможном минимуме.
		 */
		real_t wolfe1;           // параметр, принимающий промежуточную точку как допустимую (энергия)
		real_t wolfe2;           // параметр, принимающий промежуточную точку как допустимую (градиент)

		int m;                   // число запоминаемых предыдущих шагов в оптимизации
		int maxhalt;             // параметр подтверждающий конец оптимизации
			// Например, установив значение maxhalt = 3, то оптимизация прервется только, если сигнал
			// о завершении оптимизации произойдет 3 раза подряд на трех последовательных итерациях

		/*
		 *________________________ ТЕРМОДИНАМИЧЕСКИЕ ПАРАМЕТРЫ ______________________________
		 */
		_S ensemble;             // описатель ансамбля
			// формат строки -eNVT[300(C|K)], -eNVP[1.0atm], -eNPT[1atm,300K], -eNVE[3kJ/mol]

		real_t process_time;     // полное время динамики [ps]

		/*
		 * Coupling коэффициенты характеризуют насколько быстро передается тот или иной параметр
		 * термостата от него к системе. Коэффициент равен 1, если система достигает нужного значения
		 * за выбранное время релаксации. Если он мешьше, что система срелаксирует за большее время,
		 * и наоборот. Использование таких коэффициентов осмысленно для расчета систем с разными
		 * свойствами. Например, для систем с большей темплопроводностью этот коэффициент выше.
		 */
		real_t couplingT_coef;   // коэффициент для передачи температуры
		real_t couplingP_coef;   // коэффициент для передачи давления

		real_t integration_time; // шаг динамики [ps]
		real_t relaxation_time;  // время релаксации [ps]
		real_t print_time;       // время печати [ps]
		real_t average_time;     // время набора "мгновенной" статистики
			// под "мгновенной" статистикой понимается ее набор в некотором, пусть достаточно малом,
			// интервале времени (от 10 до 100 шагов динамики). Усреднение параметров в таком диапазоне
			// позволяет избежать резкого "раскачивания" системы, поскольку ее параметры резко меняются
			// от точки к точке (локальные удары атомов)

		real_t statistics_time;  // шаг сброса статистики в файл [ps]
			// Кадры текущих позиций атомов и скоростей системы сбрасываются в файл через заданное
			// время. Данное время должно быть очень длительным (0.1-10 ns), чтобы не перегружать
			// файловую систему компьютера и чтобы кадры действительно отличались друг от друга.
			// Вполне возможно, что нужно сбрасывать не сами кадры, а некоторое среднее значение по
			// набору кадров, усредненное по sampling_time. Тем самым, возможно, будут сохранены не
			// "импульсы" текущих положений, а более адекватная средняя их характеристика.

		real_t max_displacement; // максимально разрешенное смещение по координатам во время динамики
			// используется для ограничения скоростей, иначе за счет попадания в области с огромным
			// потенциалом, атомы получают огромные ускорения и следом огромные скорости, что
			// приводит к развалу системы

	public:

		Configure(int argc, char *argv[]) :
			desc_("Allowed options"), threads(1), iterations(0), process_time(0)
		{
			const char *env_variable = ::getenv("BMMKERN");
			if (!env_variable)
			{
				std::string msg = _S("you must set environment variable BMMKERN");
				PRINT_ERR(msg);
			}

			print_user_information();
				// информация для пользователя программы

			parse_(argc, argv);
			if (vm_.count("help"))
			{
				std::ostringstream oss_convert;
				oss_convert << desc_;
				_S msg = oss_convert.str();
				PRINT_MESSAGE(msg);
				exit(EXIT_SUCCESS);
			}

			data_dir = MAKE_DIRNAME(_S(env_variable));
			parse_(data_dir + config);
			PRINT_MESSAGE(this->str_());

			work_dir = MAKE_DIRNAME(_S(env_variable)) + work_dir;
		}

		/**
		*  Возвращает число элементов, соответствующих заданному ключу.
		*/
		unsigned count(const char *s) const { return (unsigned)vm_.count(s); }

		/**
		*  Возвращает число загруженных описателей молекул
		*/
		unsigned count(_I2T<ARCHETYPE_>) const { return archetypes.size(); }

		Interface_<WATER_>      get_interface(_I2T<WATER_>     ) const;
		Interface_<STATISTICS_> get_interface(_I2T<STATISTICS_>) const;
		Interface_<ENSEMBLE_>   get_interface(_I2T<ENSEMBLE_>  ) const;
		Interface_<OPTIMIZER_>  get_interface(_I2T<OPTIMIZER_> ) const;
		Interface_<ARCHETYPE_>  get_interface(_I2T<ARCHETYPE_>, unsigned i) const;

	protected:

		/**
		*  Делает парсинг командной строки. Формат командной строки должен
		*  соответствовать требованиям библиотеки boost_program_options.
		*/
		void parse_(int argc, char *argv[])
		{
		#define MAKE_STRING(reg_param, type, param, def_value, description) \
			(reg_param, boost::program_options::value<type >(&param)->default_value(def_value), description)

			desc_.add_options() ("help,h", "help message")
				MAKE_STRING("conf,c",          _S,       config,        "bmmkern.cfg", "configuration file")
				MAKE_STRING("work_dir,d",      _S,       work_dir,      "", "working directory")
				MAKE_STRING("water,w",         _S,       water,         "", "water model")
				MAKE_STRING("iterations,i",    int,      iterations,    0,  "number of the iterations")
				MAKE_STRING("simulations,t",   real_t,   process_time,  0., "the time of dynamics [ps]")
				MAKE_STRING("outfile,o",       _S,       outfile,       "", "output file")
				MAKE_STRING("statfile,s",      _S,       statfile,      "", "statistics file")
				MAKE_STRING("ensemble,e",      _S,       ensemble,      "", "ensemble")

				("box,b", boost::program_options::value<_S>(&box)->default_value("[3,3,3]"),
						"region of interaction")
				("molecule,m", boost::program_options::value<std::vector<_S> >(&archetypes),
						"molecule file[,number of archetypes | x,y,z,m], A")
			;
		#undef MAKE_STRING
			boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc_), vm_);
			boost::program_options::notify(vm_);
		}

		/**
		*  Делает парсинг конфигурационного файла. Формат этого файла должен
		*  соответствовать требованиям библиотеки boost_program_options.
		*/
		void parse_(const std::string &filename)
		{
			std::ifstream file(filename.c_str());
			if (!file)
			{
				std::string msg = _S("[ERROR] can't open file ") + filename;
				PRINT_BREAK(msg);
			}

		#define MAKE_STRING(reg_param, type, param, def_value, description) \
			(reg_param, boost::program_options::value<type >(&param)->default_value(def_value), description)

			desc_.add_options()
				MAKE_STRING("job.nthreads",              int,    threads,          DEFAULT_THREAD_NUMBER, "")
				MAKE_STRING("water.pH",                  int,    pH,               NORMAL_PH_WATER, "")
				MAKE_STRING("potentials.cutoff_radius",  real_t, cutoff_radius,    DEFAULT_CUTOFF,  "")
				MAKE_STRING("potentials.skin_radius",    real_t, rskin_width,      DEFAULT_SKIN_WIDTH, "")
				MAKE_STRING("potentials.compress_factor",real_t, rskin_compress_factor, DEFAULT_COMPRESS_COEF, "")
				MAKE_STRING("potentials.barrier",        real_t, barrier,          DEFAULT_BARRIER, "")
				MAKE_STRING("mesh.fftw_step",            real_t, mesh_fftw_step,   DEFAULT_FFTW_STEP, "")
				MAKE_STRING("optimization.stpmin",       real_t, stpmin,           DEFAULT_STPMIN, "")
				MAKE_STRING("optimization.stpmax",       real_t, stpmax,           DEFAULT_STPMAX, "")
				MAKE_STRING("optimization.xtol",         real_t, xtol,             DEFAULT_XTOL<real_t>(), "")
				MAKE_STRING("optimization.ftol",         real_t, ftol,             DEFAULT_FTOL<real_t>(), "")
				MAKE_STRING("optimization.gtol",         real_t, gtol,             DEFAULT_GTOL<real_t>(), "")
				MAKE_STRING("optimization.wolfe1",       real_t, wolfe1,           DEFAULT_WOLFE1, "")
				MAKE_STRING("optimization.wolfe2",       real_t, wolfe2,           DEFAULT_WOLFE2, "")
				MAKE_STRING("optimization.maxfev",       int,    maxfev,           DEFAULT_MAXFEV, "")
				MAKE_STRING("optimization.m",            int,    m,                DEFAULT_M, "")
				MAKE_STRING("optimization.maxhalt",      int,    maxhalt,          DEFAULT_MAX_HALTS, "")
				MAKE_STRING("optimization.steep",        int,    steep,            DEFAULT_STEEP_ITERATIONS, "")
				MAKE_STRING("dynamisc.integration_time", real_t, integration_time, DEFAULT_INTEGRATION_TIME, "")
				MAKE_STRING("dynamisc.relaxation_time",  real_t, relaxation_time,  DEFAULT_RELAXATION_TIME,   "")
				MAKE_STRING("dynamisc.print_time",       real_t, print_time,       DEFAULT_PRINT_TIME,   "")
				MAKE_STRING("dynamisc.average_time",     real_t, average_time,     DEFAULT_AVERAGE_TIME,    "")
				MAKE_STRING("dynamisc.couplingT_coef",   real_t, couplingT_coef,   DEFAULT_COUPLING_T_COEF,  "")
				MAKE_STRING("dynamisc.couplingP_coef",   real_t, couplingP_coef,   DEFAULT_COUPLING_P_COEF, "")
				MAKE_STRING("dynamisc.max_displacement", real_t, max_displacement, DEFAULT_MAX_DISPACEMENT, "")
			;

			boost::program_options::store(boost::program_options::parse_config_file(file, desc_), vm_);
			boost::program_options::notify(vm_);
		#undef MAKE_STRING

			global_thread_count = threads; // установили глобальное число потоков для программы
			global_rskin_width = rskin_width;
			global_compress_factor = rskin_compress_factor;
		}

		/**
		*  Дает текстовое представление текущего состояния объекта.
		*/
		std::string str_() const
		{
			const char *frmt = " %-35s"; _S msg = _S("");
			msg += _S("---------------------------  configuration  -------------------------------------\n");

		#define MAKE_STRING(prn_msg, param) \
			msg += prgkern::make_string(frmt, prn_msg) + prgkern::make_string(param) + _S("\n");

			MAKE_STRING("number of threads",           threads);
			MAKE_STRING("data directory",              data_dir);
			MAKE_STRING("working directory",           work_dir);
			MAKE_STRING("configuration file",          config);
			MAKE_STRING("output file",                 outfile);
			for (unsigned i=0,sz=archetypes.size(); i<sz; i++)
				MAKE_STRING("molecule", archetypes[i]);
			MAKE_STRING("water",                       water);
			MAKE_STRING("region",                      box);
			MAKE_STRING("environment pH",              pH);
			MAKE_STRING("cutoff_radius,A",             cutoff_radius);
			MAKE_STRING("barrier",                     barrier);
			MAKE_STRING("fftw mesh step,A",            mesh_fftw_step);
			msg += _S("-------------------------  optimization  -----------------------------------\n");
			MAKE_STRING("iterations",                  iterations);
			MAKE_STRING("min atom displacement,A",     stpmin);
			MAKE_STRING("max atom displacement,A",     stpmax);
			MAKE_STRING("min argument delta",          xtol);
			MAKE_STRING("min function delta",          ftol);
			MAKE_STRING("min gradient delta",          gtol);
			MAKE_STRING("1st wolf parameter",          wolfe1);
			MAKE_STRING("2nd wolf parameter",          wolfe2);
			MAKE_STRING("max function evaluations",    maxfev);
			MAKE_STRING("max previous steps usage",    m);
			MAKE_STRING("max trying to stop",          maxhalt);
			MAKE_STRING("max steep iterations",        steep);
			msg += _S("---------------------------  dynamics  -------------------------------------\n");
			MAKE_STRING("ensemble",                    ensemble);
			MAKE_STRING("simulation time,ps",          process_time);
			MAKE_STRING("integration time,ps",         integration_time);
			MAKE_STRING("relaxation time,ps",          relaxation_time);
			MAKE_STRING("console print time,ps",       print_time);
			MAKE_STRING("statistics average time,ps",  average_time);
			MAKE_STRING("statistics time,ps",          statistics_time);
			MAKE_STRING("coupling coefficient (T)",    couplingT_coef);
			MAKE_STRING("coupling coefficient (P)",    couplingP_coef);
			MAKE_STRING("statistics file",             statfile);
			msg += _S("============================================================================\n");

		#undef MAKE_STRING
			return msg;
		}
	};

	/*
	 * ------------------------------------------------------------------------------------------------
	 *                        Набор интерфейсов к Configure объекту.
	 * ------------------------------------------------------------------------------------------------
	 */
	template <> struct Interface_<ARCHETYPE_>
	{
		std::string    name;          // имя файла молекулы
		char           altpos;        // идентификатор альтернативной загрузки
		real_t         concentration; // концентрация загружаемых молекул
		bool           use_water;     // использовать молекулы кислорода (водные мостики) файла
		bool           use_hydrogens; // использовать водороды файла (может пригодиться для FORMAT_BMM_)
		unsigned       freedom_type;  // тип вставки (какие степени свободы разрешены)

		Interface_(const std::string &molecule) : altpos('A'), concentration(1.),
			use_water(false), use_hydrogens(false), freedom_type(0)
		{
			std::string s = FILE_NAME_REGEX.get_match(molecule);
			if (s != _S("")) name = s;

			s = ALTPOS_REGEX.get_match(molecule);
			if (s != _S("")) altpos = s[0];

			RegEx::match_type ss = CONCENTRATION_REGEX.get_match2(molecule);
			if (ss.first != _S(""))
			{
				concentration = make_value<real_t>(ss.first);
				if (ss.second == _S("mol")) concentration *= 1.;
				else if (ss.second == _S("mmol" )) concentration *= 1e-3;
				else if (ss.second == _S("mkmol")) concentration *= 1e-6;
				else if (ss.second == _S("nmol" )) concentration *= 1e-9;
				else if (ss.second == _S("pmol" )) concentration *= 1e-12;
				else if (ss.second == _S("fmol" )) concentration *= 1e-15;
			}

			use_water = PRN_WATER_REGEX.is_match(molecule);
			use_hydrogens = PRN_HYDROGEN_REGEX.is_match(molecule);

			s = FREEDOM_REGEX.get_match(molecule);
			for (unsigned i=0; i<s.length(); i++)
			{
				switch (s[i])
				{
				case 'u' : freedom_type |= YES_UNION_;   break;
				case 'm' : freedom_type |= YES_CM_;      break;
				case 'r' : freedom_type |= YES_ROTAMER_; break;
				case 'x' : freedom_type |= YES_ATOM_;    break;
				}
			}
		}
	};

	template <> struct Interface_<WATER_>
	{
		std::string name;        // имя воды
		real_t density;          // плотность воды
		unsigned freedom_type;   // тип вставки (какие степени свободы разрешены)

		Interface_(const Configure &config)
		: name         (_S("SPC"))
		, density      (NORMAL_WATER_DENSITY)
		, freedom_type (0)
		{
			std::string s = WATER_REGEX.get_match(config.water);
			if (s != _S("")) name = s;

			RegEx::match_type ss = DENSITY_REGEX.get_match2(config.water);
			if (ss.first != _S(""))
			{
				density = make_value<real_t>(ss.first);
				if (ss.second == _S("molecule/nm3")) density /= 1000;
			}

			s = FREEDOM_REGEX.get_match(config.water);
			for (unsigned i=0; i<s.length(); i++)
			{
				switch (s[i])
				{
				case 'm' : freedom_type |= YES_CM_;      break;
				case 'r' : freedom_type |= YES_ROTAMER_; break;
				case 'x' : freedom_type |= YES_ATOM_;    break;
				}
			}
		}
	};

	template <> struct Interface_<STATISTICS_>
	{
		_S filename;
		real_t statistics_time;
		bool prn_water;
		bool prn_hydrogens;

		Interface_(const Configure &config)
		: filename(_S("")), statistics_time(0.), prn_water(false), prn_hydrogens(false)
		{
			filename = FILE_NAME_REGEX.get_match(config.statfile);
			if (filename != _S(""))
			{
				filename = config.work_dir + filename;

				RegEx::match_type tmp = TIME_REGEX.get_match2(config.statfile);
				if (tmp.first != _S(""))
				{
					statistics_time = make_value<real_t>(tmp.first);
					statistics_time = make_time(statistics_time, tmp.second);

					prn_water = PRN_WATER_REGEX.is_match(config.statfile);
					prn_hydrogens = PRN_HYDROGEN_REGEX.is_match(config.statfile);
				}
			}
		}
	};

	template <> struct Interface_<ENSEMBLE_>
	{
		_S ensemble_type;
		real_t temperature;
		real_t pressure;
		real_t energy;

		Interface_(const Configure &config)
		: ensemble_type     ("NVT")
		, temperature       (DEFAULT_TEMPERATURE)
		, pressure          (DEFAULT_PRESSURE)
		, energy            (DEFAULT_ENERGY)
		{
			ensemble_type = ENSEMBLE_REGEX.get_match(config.ensemble);
			if (ensemble_type != _S(""))
			{
				RegEx::match_type sp = TEMPERATURE_REGEX.get_match2(config.ensemble);
				if (sp.first != _S(""))
				{
					temperature = make_value<real_t>(sp.first);
					if (sp.second == _S("C")) temperature -= DEFAULT_ABS_ZERO_TEMPERATURE;
				}

				sp = PRESSURE_REGEX.get_match2(config.ensemble);
				if (sp.first != _S(""))
				{
					pressure = make_value<real_t>(sp.first);
				}

				sp = ENERGY_REGEX.get_match2(config.ensemble);
				if (sp.first != _S(""))
				{
					energy = make_value<real_t>(sp.first);
				}
			}
		}
	};

	template <> struct Interface_<OPTIMIZER_>
	{
		_S optimizer;
		int maxiter;
		real_t stpmin;
		real_t stpmax;
		real_t xtol;
		real_t ftol;
		real_t gtol;
		real_t wolfe1;
		real_t wolfe2;
		int maxfev;
		int m;
		int maxhalt;
		int steep;

		Interface_<OPTIMIZER_>(const Configure &config)
		: optimizer(config.optimizer), maxiter(config.iterations),
		  stpmin(config.stpmin), stpmax(config.stpmax),
			xtol(config.xtol), ftol(config.ftol), gtol(config.gtol),
			wolfe1(config.wolfe1), wolfe2(config.wolfe2),
			maxfev(config.maxfev), m(config.m), maxhalt(config.maxhalt), steep(config.steep)
		{}
	};

	INLINE Interface_<WATER_> Configure::get_interface(_I2T<WATER_>) const
	{ return Interface_<WATER_>(*this); }

	INLINE Interface_<STATISTICS_> Configure::get_interface(_I2T<STATISTICS_>) const
	{ return Interface_<STATISTICS_>(*this); }

	INLINE Interface_<ENSEMBLE_> Configure::get_interface(_I2T<ENSEMBLE_>  ) const
	{ return Interface_<ENSEMBLE_>(*this); }

	INLINE Interface_<OPTIMIZER_>  Configure::get_interface(_I2T<OPTIMIZER_> ) const
	{ return Interface_<OPTIMIZER_>(*this); }

	INLINE Interface_<ARCHETYPE_>  Configure::get_interface(_I2T<ARCHETYPE_>, unsigned i) const
	{ return Interface_<ARCHETYPE_>(this->archetypes[i]); }

}
#endif

