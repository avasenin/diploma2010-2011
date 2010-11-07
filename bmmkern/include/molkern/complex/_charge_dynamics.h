#ifndef _MD_TOOL__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H
#define _MD_TOOL__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H

#include "molkern/__moldefs.h"

namespace molkern
{
	using namespace prgkern;

	template <typename _Real>
	/**
	*  Рассчитывает вектор изменений зарядов на атомах
	*   Функция не контролирует правильность размеров массивов Q1, D, I и выход
	*   за их границы.
	* @param Q1 вектор производных по времени от зарядов
	* @param incq смещения между элементами Q1
	* @param D вектор проводимостей связей
	* @param incd смещения между элементами D
	* @param nv число связей
	* @param V вектор связей
	* @param I вектор токов, связанных с атомным центром
	* @param inci смещения между элементами I
	*/
	INLINE void calculate(_I2T<CHARGE1_>, _Real *Q1, unsigned incq,
		_Real *D, unsigned incd, unsigned nv, index_<2, unsigned> *V, _Real *I, unsigned inci)
	{
		unsigned i, j; _Real d;
		for (unsigned k=0; k<nv; k++)
		{
			i = V[k][0];
			j = V[k][1];
			d = *(D + k * incd);
			*(Q1 + i * incq) -= *(I + j * inci) * d;
			*(Q1 + j * incq) -= *(I + i * inci) * d;
		}
	}

}
#endif
