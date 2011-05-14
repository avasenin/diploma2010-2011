#ifndef CONFIG__F9ED1116_3790_5da1_93E4_D94353B00A01__H
#define CONFIG__F9ED1116_3790_5da1_93E4_D94353B00A01__H

#include "prgkern/__prgkern.h"
#include "molkern/__moldefs.h"
#include <boost/program_options.hpp>

namespace molkern
{
	class Configure
	{
		boost::program_options::options_description desc_;
		boost::program_options::variables_map vm_;

	public:

		int threads; // число потоков выполнения
			// Число нитей должно согласовываться с числом ячеек в регионе (должно быть делителем
			// полного числа ячеек региона). Так как принято, что число ячеек в регионе должно быть
			// четным по всем направлениям, то максимальный делитель уже не менее 8.

		_S work_dir; // директория, откуда забираются файлы молекул
		_S data_dir; // директория, откуда забираются параметры силовых полей
		_S config; // файл параметров конфигурации расчета

		_S box; // Глобальный ящик, который задается строкой [nx,ny,nz]
			// в которой задается число ячеек ящика вдоль координат {x,y,z}.
			// Размер минимальной ячейки всегда равен радиусу взаимодействия.

		std::vector<_S> molecules;
			// (1) Формат строки "ligand.hin[N],Ahw{umrxUMRX}", где :
			//    [N] - число копий молекул,
			//    A - альтернативная позиция атомов (работает только для *.pdb).
			// Наличие/отсутствие символов h и w означает использование/игнорирование
			// атомов водорода (h) и молекул воды(w) файла
			//
			//  Тип движения молекул заданого архетипа.
			//    u   - объединие движения подцепей
			//    m   - движение центров масс
			//    r   - движение ротамеров
			//    x   - движение атомов
			//  Типы движения комбинируются произвольно, давая множество вариантов. Большие буквы (или
			//  отсутствие соответствующих малых букв) определяют запрет на движение заданных элементов.
			//    (u|U)MRX = (пустой список) фиксация в пространстве
			//    rm    движение ротамеров и центров масс подцепей
			//    urm   движение ротамеров и центра масс всего архетипа
			//    x     движение атомов при фиксации центров масс подцепей
			//    xm    наиболее свободное движение
			//    xrm   движение атомов в рутовом ротамере + ротамеры + центр масс (гибридные координаты)
			//
			// При отутствии числа копий считается, что число молекул 1.
			// При отсутствии альтернативной позиций предполагается позиция A.
			// При отсутствии скобок [] молекула размещается в координатах файла.
			// Наличие скобок [] всегда означает использование случайного размещения.
			// При размещении в координатах файла, первая молекула сдвигается так,
			// чтобы ее центр оказался в центре глобального ящика (лидер). Молекулы,
			// которые также не имеют скобок, и, которые должны также размещаться
			// в координатах файла, сдвигаются на вектор, задаваемый сдвигом лидера.
			// Таким образом, лидер (!) определяет сдвиг для всех остальных молекул.
			// Например, для протеин-лигандных комплексов протеин должен быть лидером.
			//
			// (2) Формат строки "ligand.hin[x,y,z,m],A", где :
			//   [x,y,z,m] - координата привязки локального ящика и его масштаб.
			// Формат используется для случайного размещения около указанной позиции
			// в ящике размером в радиус взаимодействия (при m=1). Коэффициенрт m
			// определяет размер локального ящика. Он показывает насколько нужно
			// умножить размеры ящика.
			// Координаты [x,y,z] всегда определяют координаты, относительно лидера.
			// Данные координаты корректируются согласно вектору сдвига лидера. При
			// отсутствии лидера текущая молекула вставляется в глобальный ящик
			// случайным образом.

		_S outfile_descriptor;
			// Формат строки "ligand.hin,hw", где :
			// наличие/отсутствие символов h и w означает сохранение/игнорирование
			// атомов водорода (h) и молекул воды(w) при записи в файл
			// Заметим, что для HIN, MOL2 и BMM файлов молекулы водородов сохраняются
			// всегда.

		_S statfile_descriptor;
			// Формат строки "ligand.bmm,w", где :
			// наличие/отсутствие символов w означает сохранение/игнорирование
			// молекул воды(w) при записи в файл

		_S water; // Формат строки "-wSPCBOX,{0|C|R|X}"
		int pH; // pH среды

		real_t cutoff_radius; // радиус обрезания потенциала
		real_t barrier; // высота барьера в 0 для 6-12 псевдопотенциала

		real_t mesh_fftw_step;
		real_t mesh_stat_step;
		real_t mesh_pots_step;
		real_t rskin_width;
		real_t rskin_compress_factor;

		_S optimizer; // имя оптимизатора
		int iterations; // число итераций оптимизации геометрии
		real_t stpmin; // минимальное смещение по координатам, обрывающее оптимизацию
		real_t stpmax; // ограничение по смещению по координатам
		real_t xtol; // минимальное значение изменения коорд, обрывающее оптимизацию
		real_t ftol; // минимальное значение изменения энергии, обрывающее оптимизацию
		real_t gtol; // минимальное значение изменения силы, обрывающее оптимизацию
		int maxfev; // максимально разрешенное число оценок энергии, обрывающее оптимизацию
		real_t wolfe1; // параметр, принимающий промежуточную точку как допустимую (энергия)
		real_t wolfe2; // параметр, принимающий промежуточную точку как допустимую (градиент)
		int m; // число запоминаемых последних результатов в опртимизации
		int maxhalt; // параметр подтверждающий конец оптимизации
		int steep; // число начальных итераций, выполняемых steep методом

		_S ensamble; // тип термодинамического ансамбля
		real_t temperature; // целевая температура
		real_t pressure; // целевое давление
		real_t density; // число частиц в 1 A**3
		unsigned process_time; // время динамики (fs), преобразовано из ввода (ps)
		unsigned integration_time; // шаг динамики (fs)
		unsigned sampling_time; // шаг набора статистики (fs)

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
		unsigned count(_I2T<MOLECULE_>) const { return molecules.size(); }

	protected:

		real_t process_time__; // переменная для преобразования вводимого времени

		/**
		*  Делает парсинг командной строки. Формат командной строки должен
		*  соответствовать требованиям библиотеки boost_program_options.
		*/
		void parse_(int argc, char *argv[])
		{
		#define MAKE_STRING(reg_param, type, param, def_value, description) \
			(reg_param, boost::program_options::value<type >(&param)->default_value(def_value), description)

			desc_.add_options()
				("help,h", "help message")
				MAKE_STRING("conf,c", _S, config, "bmmkern.cfg", "configuration file")

				("box,b", boost::program_options::value<_S>(&box)->default_value("[4,4,4]"), "region of interaction")

				MAKE_STRING("work_dir,d", _S, work_dir, "", "working directory")
				("molecule,m", boost::program_options::value<std::vector<_S> >(&molecules),
					"molecule file[,number of molecules | x,y,z,m], A")
				MAKE_STRING("water,w",  _S, water, "", "use of some water model")

				MAKE_STRING("iterations,i", int, iterations, 0, "number of the iterations")
				MAKE_STRING("md_time,t", real_t, process_time__, 0., "the time of dynamics [ps]")
				MAKE_STRING("outputfile,o", _S, outfile_descriptor, "", "output file")
				MAKE_STRING("statfile_descriptor,s", _S, statfile_descriptor, "", "statistics file")

				MAKE_STRING("ensamble,E",  _S, ensamble, "NVT", "ensamble type")
				MAKE_STRING("temperature,T",  real_t, temperature,  300., "temperature[K]")
				MAKE_STRING("pressure,P",  real_t, pressure, 1., "pressure [atm.]")
				MAKE_STRING("density,D",  real_t, density, NORMAL_WATER_DENSITY, "density [particles in A^3]")
			;
		#undef MAKE_STRING
			boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc_), vm_);
			boost::program_options::notify(vm_);

			process_time = unsigned(process_time__ / 0.001); // преобразуем к числу циклов
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
				MAKE_STRING("job.nthreads",               int,   threads, 1, "")
				MAKE_STRING("water.pH",                   int,   pH, NORMAL_PH_WATER, "")
				MAKE_STRING("potentials.cutoff_radius",   real_t, cutoff_radius, 10.,  "")
				MAKE_STRING("potentials.skin_radius",     real_t, rskin_width, 1., "")
				MAKE_STRING("potentials.compress_factor", real_t, rskin_compress_factor, 3., "")
				MAKE_STRING("potentials.barrier",         real_t, barrier,      1000., "")

				MAKE_STRING("mesh.fftw_step",           real_t, mesh_fftw_step,  1., "")
				MAKE_STRING("mesh.stat_step",           real_t, mesh_stat_step,  1., "")
				MAKE_STRING("mesh.pots_step",           real_t, mesh_pots_step,  1., "")

				MAKE_STRING("optimization.stpmin",      real_t, stpmin,  DEFAULT_STPMIN, "")
				MAKE_STRING("optimization.stpmax",      real_t, stpmax,  DEFAULT_STPMAX, "")
				MAKE_STRING("optimization.xtol",        real_t, xtol,    DEFAULT_XTOL<real_t>(), "")
				MAKE_STRING("optimization.ftol",        real_t, ftol,    DEFAULT_FTOL<real_t>(), "")
				MAKE_STRING("optimization.gtol",        real_t, gtol,    DEFAULT_GTOL<real_t>(), "")
				MAKE_STRING("optimization.wolfe1",      real_t, wolfe1,  DEFAULT_WOLFE1, "")
				MAKE_STRING("optimization.wolfe2",      real_t, wolfe2,  DEFAULT_WOLFE2, "")
				MAKE_STRING("optimization.maxfev",         int, maxfev,  DEFAULT_MAXFEV, "")
				MAKE_STRING("optimization.m",              int, m,       DEFAULT_M, "")
				MAKE_STRING("optimization.maxhalt",        int, maxhalt, DEFAULT_MAX_HALTS, "")
				MAKE_STRING("optimization.steep",          int, steep,   DEFAULT_STEEP_ITERATIONS, "")

				MAKE_STRING("dynamisc.integration_time", unsigned, integration_time,  1, "")
				MAKE_STRING("dynamisc.sampling_time",    unsigned,    sampling_time, 10, "")
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
			const char *frmt = " %-35s";
			_S msg = _S("\n============================ CONFIGURATION =================================\n");

		#define MAKE_STRING(reg_name, prn_msg, param) if (vm_.count(reg_name)) \
			msg += prgkern::make_string(frmt, prn_msg) + prgkern::make_string(param) + _S("\n");

			msg += prgkern::make_string(frmt, "data directory") + prgkern::make_string(data_dir) + _S("\n");
			MAKE_STRING("work_dir", "working directory",  work_dir);
			MAKE_STRING("conf",     "configuration file", config);
			MAKE_STRING("job.nthreads",  "number of threads", threads);
			MAKE_STRING("box",      "region",             box);

			for (unsigned i=0,sz=molecules.size(); i<sz; i++)
				msg += prgkern::make_string(frmt, "molecule ") + molecules[i] + _S("\n");
			MAKE_STRING("water", "water name", water);
			MAKE_STRING("water.pH", "environment pH", pH);

			MAKE_STRING("outputfile", "output file name",     outfile_descriptor);
			MAKE_STRING("statfile_descriptor",   "statistics file name", statfile_descriptor);

			MAKE_STRING("potentials.cutoff_radius", "cutoff_radius(A)", cutoff_radius);
			MAKE_STRING("potentials.barrier",       "barrier",       barrier);
			MAKE_STRING("mesh.fftw_step", "fftw mesh step(A)",          mesh_fftw_step);
			MAKE_STRING("mesh.stat_step", "statistic mesh step(A)",     mesh_stat_step);
			MAKE_STRING("mesh.pots_step", "potent mesh step(A)",        mesh_pots_step);

			msg += _S("-------------------------  optimization  -----------------------------------\n");
			MAKE_STRING("iterations",                "iterations",      iterations);
			MAKE_STRING("optimization.stpmin",       "min coordinate shift(A)", stpmin);
			MAKE_STRING("optimization.stpmax",       "max coordinate shift(A)", stpmax);
			MAKE_STRING("optimization.xtol",         "min coordinate shift(A)", xtol);
			MAKE_STRING("optimization.ftol",         "min energy shift", ftol);
			MAKE_STRING("optimization.gtol",         "min force shift", gtol);
			MAKE_STRING("optimization.wolfe1",       "1st wolf parameter", wolfe1);
			MAKE_STRING("optimization.wolfe2",       "2nd wolf parameter", wolfe2);
			MAKE_STRING("optimization.maxfev",       "max function evaluations", maxfev);
			MAKE_STRING("optimization.m",            "number of previous steps usage", m);
			MAKE_STRING("optimization.maxhalt",      "number of trying to stop", maxhalt);
			MAKE_STRING("optimization.steepit",      "number of steep iterations", steep);

			msg += _S("---------------------------  dynamics  -------------------------------------\n");
			MAKE_STRING("ensamble",                  "type of ensamble",     ensamble);
			MAKE_STRING("md_time",                   "full time(ps)",        process_time * 0.001);
			MAKE_STRING("temperature",               "temperature(K)",       temperature);
			MAKE_STRING("pressure",                  "pressure(atm)",        pressure);
			MAKE_STRING("density",                   "density(A**-3)",       density);
			MAKE_STRING("dynamisc.integration_time", "integration time(fs)", integration_time);
			MAKE_STRING("dynamisc.sampling_time",    "sample time(fs)",      sampling_time);
			msg += _S("============================================================================\n");

		#undef MAKE_STRING
			return msg;
		}
	};

	/**
	*  Наборы описателей тех или иных объектов. Описатели парсят текстовые строки,
	*  задаваемые в командной строке конфигуратора, и формируют для них
	*  набор параметров объектов.
	*/
	template <int TYPE> struct Descriptor_;

	/**
	*   Описатель молекулы и ее месторасположения в пространстве.
	*/
	template <> struct Descriptor_<MOLECULE_>
	{
		std::string name; // имя файла молекулы
		char altpos; // идентификатор альтернативной загрузки
		unsigned count; // количество загружаемых молекул
		bool is_random; // использование случайного позиционирования
		vector_t X; // точка привязки локального ящика (при count = 1)
		real_t box_scale; // масштаб локального ящика
		bool use_hydrogens; // использовать водороды файла (может пригодиться для FORMAT_BMM_)
		bool use_water; // использовать молекулы кислорода (водные мостики) файла
		unsigned freedom_type; // тип вставки (какие степени свободы разрешены)

		Descriptor_(std::string s) : altpos('A'), count(1), is_random(false),
			box_scale(1.), use_hydrogens(false), use_water(false), freedom_type(0)
			// установили значения по умалчиванию для случая координатного файла
			// и которые, если будет необходимость, переопределим
		{
		#define FORMAT_CONTROL(s, pos__) \
			if (pos__ == (int)std::string::npos) \
			{ \
				_S msg = _S("the molecule description string is : ") + s + _S("\n"); \
				msg += _S("the write format for description string is corrupted."); \
				PRINT_ERR(msg); \
			}

			int pos = s.find_first_of("[,");
			name = s.substr(0, pos);
				// независимо, от того, что вернет find_first_of()

			if (pos == (int)std::string::npos) return;
				// поскольку значения по умалчиванию установлены
			if (s[pos] == '[')
			{
				is_random = true;

				int pos__ = s.find_first_of("]", pos + 1);
				FORMAT_CONTROL(s, pos__)
				_S s__ = s.substr(pos + 1, pos__ - pos - 1);
					// сохранили содержимое [N] или [x,y,z,m]

				pos__ = s__.find_first_of(",");
				if (pos__ == (int)std::string::npos) // случай [N]
					count = atoi(s__.c_str());
				else
				{
					const char *p = s__.c_str();
					X[0] = (real_t)atof(p);
					X[1] = (real_t)atof(p + pos__ + 1);

					pos__ = s__.find_first_of(",", pos__ + 1);
					FORMAT_CONTROL(s, pos__)

					X[2] = (real_t)atof(p + pos__ + 1);
					pos__ = s__.find_first_of(",", pos__ + 1);
					FORMAT_CONTROL(s, pos__)

					box_scale = (real_t)atof(p + pos__ + 1);
				}
			}
		#undef FORMAT_CONTROL

			// разберем остаток строки "Ahw"
			s = s.substr(pos + 1); // остаток

			for (unsigned i=0; i<s.length(); i++)
			{
				switch (s[i])
				{
				case 'h' : use_hydrogens = true; break;
				case 'w' : use_water     = true; break;
				case 'A' : altpos = 'A'; break;
				case 'B' : altpos = 'B'; break;
				case 'C' : altpos = 'C'; break;

				// набираем побитово
				case 'u' : freedom_type |= YES_UNION_;   break;
				case 'm' : freedom_type |= YES_CM_;      break;
				case 'r' : freedom_type |= YES_ROTAMER_; break;
				case 'x' : freedom_type |= YES_ATOM_;    break;
				}
			}
		}
	};

	template <> struct Descriptor_<WATER_>
	{
		std::string name; // имя воды
		unsigned freedom_type; // тип вставки (какие степени свободы разрешены)

		Descriptor_(std::string s) : freedom_type(0)
			// установили значения по умалчиванию для случая координатного файла
			// и которые, если будет необходимость, переопределим
		{
		#define FORMAT_CONTROL(s, pos__) \
			if (pos__ == (int)std::string::npos) \
			{ \
				_S msg = _S("[ERROR] The water description string is : ") + s + _S("\n"); \
				msg += _S("[ERROR] The write format for description string is corrupted."); \
				PRINT_ERR(msg); \
			}

			int pos = s.find_first_of(",");
			name = s.substr(0, pos);
				// независимо, от того, что вернет find_first_of()
			if (name.length() == 0) name = _S("SPCBOX");

			if (pos == (int)std::string::npos) return;
				// поскольку значения по умалчиванию установлены

			// разберем остаток строки "Ahw"
			s = s.substr(pos + 1); // остаток

			for (unsigned i=0; i<s.length(); i++)
			{
				switch (s[i])
				{
				// набираем побитово
				case 'm' : freedom_type |= YES_CM_;      break;
				case 'r' : freedom_type |= YES_ROTAMER_; break;
				case 'x' : freedom_type |= YES_ATOM_;    break;
				}
			}
		}
	};

	/**
	*   Описатель выходных файлов.
	*/
	template <> struct Descriptor_<FILE_>
	{
		std::string name; // имя файла выводв для комплекса молекул
		bool use_hydrogens; // наличие водородов в выводе
		bool use_water; // наличие молкул воды в выводе

		Descriptor_(std::string s) : name(""), use_hydrogens(false), use_water(false)
			// установили значения по умалчиванию для случая координатного файла
			// и которые, если будет необходимость, переопределим
		{
			int pos = s.find_first_of(",");
			name = s.substr(0, pos);
				// независимо, от того, что вернет find_first_of()

			if (pos == (int)std::string::npos) return;
				// поскольку значения по умалчиванию установлены

			s = s.substr(pos + 1); // остаток строки
			std::set<char> symbols_in_line;
			for (unsigned i=0,sz=s.length(); i<sz; i++)
				symbols_in_line.insert(s[i]);

			std::set<char>::iterator it = symbols_in_line.begin(),
				ite = symbols_in_line.end();
			for (; it!=ite; ++it)
			{
				if (*it == 'h') use_hydrogens = true;
				if (*it == 'w') use_water = true;
			}
		}
	};

	/**
	*   Описатель параметров оптимизатора.
	*/
	template <> struct Descriptor_<OPTIMIZER_>
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

		Descriptor_<OPTIMIZER_>(const Configure &config)
		: optimizer(config.optimizer), maxiter(config.iterations),
		  stpmin(config.stpmin), stpmax(config.stpmax),
			xtol(config.xtol), ftol(config.ftol), gtol(config.gtol),
			wolfe1(config.wolfe1), wolfe2(config.wolfe2),
			maxfev(config.maxfev), m(config.m), maxhalt(config.maxhalt), steep(config.steep)
		{}
	};

	/**
	*   Описатель параметров термостата.
	*/
	template <> struct Descriptor_<THERMOSTATE_>
	{
		_S ensamble; // имя термодинамического ансамбля
		real_t temperature; // температура ансамбля
		real_t pressure; // давление
		unsigned integration_time; // шаг интегрирования
		unsigned sampling_time; // шаг считывания статистики

		Descriptor_<THERMOSTATE_>(const Configure &config) : ensamble(config.ensamble),
			temperature(config.temperature), pressure(config.pressure),
			integration_time(config.integration_time), sampling_time(config.sampling_time)
		{}
	};

}
#endif

