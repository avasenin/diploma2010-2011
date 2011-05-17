#include "molkern/__molkern.h"
/**
*  Определения всех переменных, которые вынуждено должны быть в *.cpp.
*/

namespace molkern
{

#ifdef USE_PSEUDO_POTENTIAL_IN_ZERO
	const real_t Interaction_<C612>::c1__ = (real_t)( 0.1132524796725085633035822270098921L);
	const real_t Interaction_<C612>::c2__ = (real_t)(-0.1009131954055306831946948492286285L);
	const real_t Interaction_<C612>::sqrt_pi = (real_t)( 1.7724538509055158819194275565678254L);

	real_t Interaction_<C612>::alpha_;
	real_t Interaction_<C612>::tau2_;
	real_t Interaction_<C612>::gamma_;
#endif

	real_t Interaction_<C612>::shift_factor_;
	real_t Interaction_<C612>::rcutoff_;
	real_t Interaction_<E612>::alpha_;

	model_timer_t global_model_time;     // счетчик времени динамики
	system_time_t global_start_time;     // физическое время старта динамики [sec]

	real_t global_rskin_width;
	real_t global_compress_factor;

	unsigned global_thread_count; // число запущенных процессов
	RappleGoddardParams* RappleGoddardParams::s_instance = NULL;
	CoulombParams* CoulombParams::s_instance = NULL;

};
