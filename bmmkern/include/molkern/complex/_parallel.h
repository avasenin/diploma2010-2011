#ifndef _PARALLEL__09ED1116_C3FE_55a4_CA63_F84536C80D11__H
#define _PARALLEL__09ED1116_C3FE_55a4_CA63_F84536C80D11__H

#include "boost/thread/thread.hpp"
#include "boost/thread/barrier.hpp"

#include "molkern/__moldefs.h"

namespace molkern
{
	using namespace prgkern;

	/**
	 * Менеджер параллельного выполнения задач для заданного объекта класса.
	 * Например, для обеспечения счета LJ взаимодействий в параллельном режиме.
	 * Класс написан таким образом, что может подключаться на выполнение любого
	 * количества любых функций, но одну за раз. Подключемые функции должны
	 * иметь заданный интерфейс. Класс "не знает", что за функции он выполняет,
	 * потому он независим ни от выходных, ни от входных параметров этих функций.
	 */
	template <typename Object>
	class ParallWorker_
	{
		std::vector<boost::thread *> thread_; // потоки
		boost::barrier *pre__barrier_; // барьер синхронизации перед стартом работы
		boost::barrier *post_barrier_; // барьер синхронизации после выполнения работы
		unsigned thread_count_; // число потоков для обработки

		typedef void (Object::*_Fn)(unsigned, void *, void *);
			// 1 параметр - номер нити
			// 2 параметр - возвращаемое значение
			// 3 параметр - параметры функции

		Object *object_; // комплекс для которого выполняется работа
		_Fn fn_; // функция объекта, выполняемая в параллельном режиме
		void *retval_; // начало массива для складывания результатов (return)
		void *params_; // параметры функции объекта

		struct thread_data_
		{
			ParallWorker_ *worker_; // ссылка на нахождение данных
			unsigned nth_; // номер нити, чтобы нить знала откуда брать данные

			thread_data_(ParallWorker_ *worker, unsigned nth) : worker_(worker), nth_(nth) {}
			void operator()() // оператор, необходимый для boost::thread
			{
				while(true)
				{
					worker_->pre__barrier_->wait(); // подождем установки работы
					if (worker_->fn_ == 0) break; // отсутствие функции - есть сигнал завершения
					(worker_->object_->*worker_->fn_)(nth_, worker_->retval_, worker_->params_);
					worker_->post_barrier_->wait(); // стоп, чтобы внешняя нить могла собрать результаты
				}
			}
		};

	public:

		~ParallWorker_()
		{
			fn_ = 0;
			if (pre__barrier_) pre__barrier_->wait(); // сигнал завершения для нитей
			for (unsigned i=0; i<thread_count_; i++) if (thread_[i]) thread_[i]->join();

			for (unsigned i=0; i<thread_count_; i++) delete thread_[i];
			delete pre__barrier_;
			delete post_barrier_;
		}

		/// инициализация планировщика параллельных задач
		void init(Object *object, unsigned thread_count=global_thread_count)
		{
			thread_count_ = thread_count;
			object_ = object;

			pre__barrier_ = new boost::barrier(thread_count_ + 1);
			post_barrier_ = new boost::barrier(thread_count_ + 1);
				// расчетные нити + основная нить потока

			thread_.resize(thread_count_);
			for (unsigned i=0; i<thread_count_; i++)
				thread_[i] = new boost::thread(thread_data_(this, i));
		}

		/// число нитей выполнения
		unsigned thread_count() const { return thread_count_; }

		/// запуск задачи в параллельном режиме
		void start_parallel_running(_Fn fn, void *retval=0, void *params=0)
		{
			fn_ = fn;
			retval_ = retval;
			params_ = params;
			pre__barrier_->wait();
		}

		/// ожидание завершения выполнения всех нитей
		void wait() { post_barrier_->wait(); }

	};

	/**
	 * Менеджер параллельного выполнения задач для заданного объекта класса.
	 * Класс написан таким образом, что может подключаться на выполнение любого
	 * количества любых функций, но одну за раз. Подключемые функции должны
	 * иметь заданный интерфейс. Класс "не знает", что за функции он выполняет,
	 * потому он независим ни от выходных, ни от входных параметров этих функций.
	 * Менеджер сам распределяет задания по своим потокам. Вызывающий поток не
	 * обязан знать, сколько нитей у менеджера.
	 */
	template <typename Object>
	class Parall_
	{
		std::vector<boost::thread *> thread_; // потоки
		boost::barrier *pre__barrier_; // барьер синхронизации перед стартом работы
		boost::barrier *post_barrier_; // барьер синхронизации после выполнения работы
		unsigned thread_count_; // число потоков для обработки

		typedef void (Object::*_Fn)(void *, void *, void *);
			// 1 параметр - возвращаемое значение
			// 2 параметр - параметры функции

		Object *object_; // комплекс для которого выполняется работа
		_Fn fn_; // функция объекта, выполняемая в параллельном режиме
		void *params_;

		void *retval_; // адрес начала памяти для сохранения результата
		void *parval_; // адрес начала памяти для чтения данных
		unsigned retlen_; // длина одной порции результата (в байтах)
		unsigned parlen_; // длина одной порции данных для счета (в байтах)
		unsigned count_; // число порций данных (число вызовов внешней функции,
			// по одной для каждой порции данных)

		struct thread_data_
		{
			Parall_ *worker_; // ссылка на нахождение данных
			unsigned nth_; // номер нити, чтобы нить знала откуда брать данные

			thread_data_(Parall_ *worker, unsigned nth) : worker_(worker), nth_(nth) {}
			void operator()() // оператор, необходимый для boost::thread
			{
				while (true)
				{
					worker_->pre__barrier_->wait(); // подождем установки работы
					if (worker_->fn_ == 0) break; // отсутствие функции - есть сигнал завершения

					// определим границы подмассивов и число элементов для данной нити
					unsigned odd = worker_->count_ % worker_->thread_count_? 1 : 0;
					unsigned count = (unsigned)(worker_->count_ / worker_->thread_count_) + odd;

					char *rval = (char *)worker_->retval_ + nth_ * worker_->retlen_;
					char *pval = (char *)worker_->parval_ + nth_ * count * worker_->parlen_;
						// обработка выхода начала за границы массива делается параметром count

					if (odd) // корректируем число элементов для последних нитей
					{
						unsigned s = nth_ * count;
						if (s + count > worker_->count_) count = worker_->count_ - s;
						if (s > worker_->count_) count = 0;
					}

					// запускаем нить на счет по всем элементам подмассива
					for (unsigned i=0; i<count; i++, pval+=worker_->parlen_)
						(worker_->object_->*worker_->fn_)((void *)rval, (void *)pval, worker_->params_);

					worker_->post_barrier_->wait(); // стоп, чтобы внешняя нить могла собрать результаты
				}
			}
		};

	public:

		~Parall_()
		{
			fn_ = 0;
			if (pre__barrier_) pre__barrier_->wait(); // сигнал завершения для нитей
			for (unsigned i=0; i<thread_count_; i++) if (thread_[i]) thread_[i]->join();

			for (unsigned i=0; i<thread_count_; i++) delete thread_[i];
			delete pre__barrier_;
			delete post_barrier_;
		}

		/// инициализация планировщика параллельных задач
		void init(Object *object, unsigned thread_count=global_thread_count)
		{
			thread_count_ = thread_count;
			object_ = object;

			pre__barrier_ = new boost::barrier(thread_count_ + 1);
			post_barrier_ = new boost::barrier(thread_count_ + 1);
				// расчетные нити + основная нить потока

			thread_.resize(thread_count_);
			for (unsigned i=0; i<thread_count_; i++)
				thread_[i] = new boost::thread(thread_data_(this, i));
		}

		/// число нитей выполнения
		unsigned thread_count() const { return thread_count_; }

		/**
		 * Запуск задачи в параллельном режиме.
		 * @param fn функция класса Object, выполняемая в параллельном режиме
		 * @param rarr массив для записи выходных значений (число элементов = thread_count)
		 * @param n полное число элементов массива parr
		 * @param parr массив для чтения входных значений (число элементов = n)
		 * @param params возможный объект для счета
		 */
		template <typename R, typename T>
		void start_parallel_running(void (Object::*fn)(void *, void *, void *),
			R *rarr, unsigned n, const T *parr, void *params=0)
		{
			fn_ = fn;
			retval_ = (void*)rarr; retlen_ = sizeof(R);
			parval_ = (void*)parr; parlen_ = sizeof(T);
			count_ = n;
			params_ = params;

			pre__barrier_->wait();
		}

		/// ожидание завершения выполнения всех нитей
		void wait() { post_barrier_->wait(); }

	};
}
#endif

