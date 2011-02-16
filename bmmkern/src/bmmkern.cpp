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
	Complex::region_type region(config.box, config.cutoff_radius);
	ForceField forcefield(config.data_dir);
	Residome residome(config.data_dir);

	//----------------------------------------------------------------------------
	//                              построение
	//----------------------------------------------------------------------------
	Complex complex(&forcefield, &residome, &region, config.pH, config.density);
		// создадим пустой комплекс

	// загрузим и построим все молекулы
	for (unsigned i=0,sz=config.count(MOLECULE); i<sz; i++)
	{
		Descriptor_<MOLECULE_> desc(config.molecules[i]);
			// конвертируем текстовую строку, описывающую молекулу, в набор параметров

		_S filename = config.work_dir + desc.name;

		// загрузим архетип молекулы в регион
		complex.load(MOLECULE, filename, desc.freedom_type, desc.count, desc.altpos);
	}

	// загрузим архетип воды
	if (config.water != _S(""))
	{
		Descriptor_<WATER_> desc(config.water);
		complex.load(WATER, desc.name, desc.freedom_type);
	}

	// разместим молекулы и воду
	complex.build(YES_CLASHES);

	//----------------------------------------------------------------------------
	//                              оптимизация
	//----------------------------------------------------------------------------
	if (config.iterations >= 0)
	{
		// распечатаем начальное состояние для контроля построения
		complex.U();

		Descriptor_<OPTIMIZER_> desc(config); // загрузим описание оптимизатора
		ChargeOptimizer_<LMBFGS_> optimizer; // инициализируем оптимизатор

		// выполним оптимизацию
		optimizer(&complex, desc);

		// распечатаем конечное состояние для контроля оптимизации
		complex.U();

		Descriptor_<FILE_> file(config.outfile_descriptor);
		if (file.name != _S(""))
		{
			_S filename = config.work_dir + file.name;
			complex.save(filename, file.use_water, file.use_hydrogens);
		}
	}

	//----------------------------------------------------------------------------
	//                              динамика
	//----------------------------------------------------------------------------
	if (config.process_time > 0.)
	{
		Descriptor_<THERMOSTATE_> desc(config); // загрузим описание термостата
		Thermostat thermostat(desc, &complex);
			// контактируем молекулу с термостатом и установим ее начальное состояние
			// (температуру, или полную энергию, или давление)

		Descriptor_<FILE_> statfile(config.statfile_descriptor);
		_S filename = config.work_dir + statfile.name;

		Integrator integrator(&complex, config.integration_time,
			config.sampling_time, filename, statfile.use_water);
			// контактируем молекулу с интегратором и установим его начальное состояние

		integrator.run(&thermostat, &complex, config.process_time);

		Descriptor_<FILE_> outfile(config.outfile_descriptor);
		if (outfile.name != _S(""))
		{
			_S filename = config.work_dir + outfile.name;
			complex.save(filename, outfile.use_water, outfile.use_hydrogens);
		}
	}

	// распечатаем конечное состояние для контроля оптимизации
	complex.U();

	return 0;
}

