#ifndef _ENSEMBLE__F9ED1116_EDB9_5e17_25FF_F745B15D0101__H
#define _ENSEMBLE__F9ED1116_EDB9_5e17_25FF_F745B15D0101__H

#include "molkern/__moldefs.h"
#include "molkern/__config.h"

namespace molkern
{
	using namespace prgkern;

	template <typename _System>
	class Ensemble_ : public notified_object_
	{
		mutable _System *system_;     // система, для которой накапливается статистика

		virtual void process_notification()
		{
			notified_object_::process_notification();

			sample_statistics(); // "мгновенная" статистика набирается на каждом шаге

			if (counter_ % notification_time_ == 0) // печать по заданому периоду времени
				print();
		}

	public:

		real_t N_; // число частиц
		Average_<real_t> U_; // статистика по потенциальной энергии системы
		Average_<real_t> V_; // статистика по объемy системы
		Average_<real_t> T_; // статистика по температуре в системе
		Average_<real_t> P_; // статистика по давлению в системе

		Ensemble_(const Configure *conf, _System *system)
		: notified_object_  (round(conf->print_time / conf->integration_time))
		, system_           (system)
		,	U_                (round(conf->average_time / conf->integration_time))
		, V_                (round(conf->average_time / conf->integration_time))
		, T_                (round(conf->average_time / conf->integration_time))
		, P_                (round(conf->average_time / conf->integration_time))
		{}

		void sample_statistics()
		{
			typedef typename _System::atom_type            _Atom;
			typedef array_iterator<_Atom, range_iterator>  _Iterator;

			unsigned N = system_->count(ATOM);
			real_t volume = system_->get(VOLUME);

			_E(real_t) kenergy = 0.;
			_E(real_t) virial = 0.;

			_Iterator it = system_->make_iterator(range_iterator(0));
			_Iterator ite = system_->make_iterator(range_iterator(N));
			for (; it!=ite; ++it)
			{
				const _Atom &atom = *it;
				kenergy += atom->mass * scalar_product(atom.V, atom.V);
				virial += scalar_product(atom.X, atom.F);
			}

			N_ = N;
			U_.push(system_->get(POTENT_ENERGY));
			V_.push(volume);
			T_.push((real_t)Temperature(kenergy / (3 * N)));
			P_.push((real_t)Pressure((kenergy + virial) / (3 * N)));
		}

		/**
		*  Печать статистики текущего состояния.
		*/
		std::string print_statistics() const
		{
			float P = P_.average();
			float V = V_.average();
			float T = T_.average();
//			float P = P_.top();
//			float V = V_.top();
//			float T = T_.top();

			float K =  1.5 * N_ * KT(T);
			float U = U_.average();
			float E = K + U;
			return make_string("PVT={%10.3e %9.3e %4.0f}  KUE={%10.3e %10.3e %12.5e}",
				P, V, T, K, U, E);
		}

		void print() const
		{
			_S msg = make_string(global_model_time) + _S(" ") + print_statistics()
				+ make_string(current_time() - global_start_time);
			PRINT_MESSAGE(msg);
		}

	};
	typedef Ensemble_<Complex>   Ensemble;


	template <typename _System>
	class Statistics_ : public notified_object_
	{
		_System *system_;             // объект, для которого набирается глобальная статистика
		real_t integration_time_;     // время шага динамики (требуется для печати статистики)
		std::ofstream file_;          // файловый объект для сброса статистики
		std::string filename_;        // имя файла для сброса статистики
		bool file_opened_;            // открыт/закрыт файл?
		bool prn_water_;              // печатать ли в файл молекулы воды
		bool prn_hydrogens_;          // печатать ли в файл водороды

		virtual void process_notification()
		{
			notified_object_::process_notification();

			if (counter_ % notification_time_ == 0)
			{
				_S time = make_string(global_model_time);
				save(make_string("FRAME %s", time.c_str()));
			}
		}

	public:

		Statistics_() : system_(0), file_opened_(false) {}

		Statistics_(const Configure *conf, _System *system)
		: notified_object_(1), system_(system), filename_("")
		, file_opened_(false), prn_water_(false), prn_hydrogens_(false)
		{
			Interface_<STATISTICS_> interface = conf->get_interface(STATISTICS);

			filename_ = interface.filename;
			real_t statistics_time = interface.statistics_time;
			prn_water_ = interface.prn_water;
			prn_hydrogens_ = interface.prn_hydrogens;

			if (filename_ != _S(""))
			{
				set_notification_period(round(statistics_time / conf->integration_time));
				open(filename_, system);
			}
			else
			{
				_S msg = _S("\n[WARNING] The statistics would not be saved !\n");
				PRINT_MESSAGE(msg);
			}
		}

		~Statistics_() { if (file_opened_) close(); }

		void open(std::string &filename, _System *system)
		{
			std::string ext = extension(filename);
			if (ext == _S(".pdb" ) || ext == _S(".ent" ) || ext == _S(".hin" ) || ext == _S(".mol2")
				|| ext == _S(".bmm" ))
			{
				if (ext != _S(".bmm" )) file_.open(filename.c_str(), std::ios_base::out);
				else file_.open(filename.c_str(), std::ios_base::out | std::ios_base::binary);

				if (!file_)
				{
					std::string msg = _S("can't open statistics file ") + filename;
					PRINT_ERR(msg);
				}
				filename_ = filename;

				// сохраним ящик и число атомов
				if (ext == _S(".bmm" )) system->write_header(_I2T<FORMAT_BMM_>(), file_);
				file_opened_ = true;
			}
		}

		void save(const std::string &header)
		{
			if (file_opened_)
			{
				std::string ext = extension(filename_);

				if      (ext == _S(".pdb" )) system_->save(FORMAT_PDB,  file_, prn_water_, prn_hydrogens_, header);
				else if (ext == _S(".ent" )) system_->save(FORMAT_PDB,  file_, prn_water_, prn_hydrogens_, header);
				else if (ext == _S(".hin" )) system_->save(FORMAT_HIN,  file_, prn_water_, prn_hydrogens_, header);
				else if (ext == _S(".mol2")) system_->save(FORMAT_MOL2, file_, prn_water_, prn_hydrogens_, header);
				else if (ext == _S(".bmm" )) system_->save(FORMAT_BMM,  file_, prn_water_, prn_hydrogens_, header);

				file_.flush(); // реальная запись, чтобы можно было оборвать выполнение и сохранить файл
			}
		}

		void close()
		{
			if (file_opened_)
			{
				file_.close();
				SAVED_OK_MESSAGE(filename_);
			}
		}
	};
	typedef Statistics_<Complex>  Statistics;

}
#endif
