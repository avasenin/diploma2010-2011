/***************************************************************************
 *              Copyright (C) 2006-2009 by Eduard S. Fomin                 *
 *                       fomin@bionet.nsc.ru                               *
 ***************************************************************************/
#include "molkern/__molkern.h"

using namespace prgkern;
using namespace molkern;

int main(int argc, char *argv[])
{
	Configure config(argc, argv);

	//----------------------------------------------------------------------------
	//                          загрузка баз данных
	//----------------------------------------------------------------------------
#ifdef USE_ERF
	Interaction interaction();
#else
	Interaction interaction(config.cutoff_radius, config.barrier);
#endif
	ForceField forcefield(config.data_dir);
	Residome residome(config.data_dir);

	//----------------------------------------------------------------------------
	//                              построение
	//----------------------------------------------------------------------------
	Complex complex(&config, &forcefield, &residome);
		// создадим пустой комплекс

	// загрузим и построим все молекулы
	for (unsigned i=0,sz=config.count(ARCHETYPE); i<sz; i++)
	{
		Interface_<ARCHETYPE_> interface = config.get_interface(ARCHETYPE, i);
		_S filename = config.work_dir + interface.name;
		unsigned count = 1;
		complex.load(ARCHETYPE, filename, interface.freedom_type, count, interface.altpos);
	}

	// загрузим архетип воды
	if (config.water != _S(""))
	{
		Interface_<WATER_> interface = config.get_interface(WATER);
		complex.load(WATER, interface.name, interface.freedom_type);
	}

	// разместим молекулы и воду
	complex.build(YES_CLASHES);

	// распечатаем начальное состояние для контроля построения
	complex.U();

	//----------------------------------------------------------------------------
	//                              оптимизация
	//----------------------------------------------------------------------------
	if (config.iterations >= 0)
	{
		Interface_<OPTIMIZER_> interface = config.get_interface(OPTIMIZER);
		Optimizer_<LMBFGS_> optimizer; // инициализируем оптимизатор

		// выполним оптимизацию
		optimizer.run(&complex, interface);

		// распечатаем конечное состояние для контроля оптимизации
		complex.U();
	}

	//----------------------------------------------------------------------------
	//                              динамика
	//----------------------------------------------------------------------------
	if (config.process_time > 0.)
	{
		TIME_TESTING_START("Full cycle of dynamisc finished...", 1)

		Ensemble ensamble(&config, &complex);         // накопитель "мгновенной" статистики
		Statistics statistics(&config, &complex);     // файл глобальной статистики

		Thermostat thermostat(&config, &ensamble, &complex);
			// контактируем молекулу с термостатом и установим начальное состояние:
			// температуру, энергию, или давление

		Integrator integrator(&config, &ensamble, &complex);
			// контактируем молекулу с интегратором и установим его начальное состояние

		// зарегистрируем уведомляемые объекты у интегратора
		integrator.register_notified_object(&thermostat, 1);
		integrator.register_notified_object(&ensamble,   2);
		integrator.register_notified_object(&statistics   );

		integrator.run(&thermostat);
		TIME_TESTING_FINISH;

		// распечатаем конечное состояние для контроля динамики
		complex.U();
	}

	return 0;
}

