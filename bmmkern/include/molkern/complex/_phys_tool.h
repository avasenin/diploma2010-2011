#ifndef _PHYS_TOOL__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H
#define _PHYS_TOOL__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H

#include "molkern/__moldefs.h"

namespace molkern
{
	using namespace prgkern;

	/**
	* @brief calculate mass of complex molecule
	* @return the mass of molecule
	*/
	template <typename Molecule, typename Iterator>
	INLINE typename Molecule::real_type
	calculate(_I2T<MASS_>, const Molecule *molecule, Iterator it, Iterator ite)
	{
		typedef typename Molecule::real_type     _Real;
		typedef typename Molecule::atom_type     _Atom;
		typedef array_iterator<_Atom, Iterator>  _Iterator;

		_Iterator it__ = molecule->make_iterator(it);
		_Iterator ite__ = molecule->make_iterator(ite);

		_Real mass = 0.;
		for (; it__!=ite__; ++it__) mass += (*it__).mass;
		return mass;
	}

	template <typename Molecule>
	INLINE typename Molecule::real_type
	calculate(_I2T<MASS_>, const Molecule *molecule)
	{
		unsigned count = molecule->count(_I2T<ATOM_>());
		return calculate(_I2T<MASS_>(), molecule,
			range_iterator(0), range_iterator(count));
	}

	/**
	* @brief calculate mass center of molecule
	* @return the mass center of molecule
	*/
//	template <typename Molecule, typename Iterator>
//	INLINE const vdense_<Molecule::dimension, real_t>
//	calculate(_I2T<MASS_CENTER_>, const Molecule *molecule,
//		Iterator it, Iterator ite)
//	{
//		typedef vdense_<Molecule::dimension, real_t>   _Point;
//		typedef typename Molecule::atom_type           _Atom;
//		typedef const_array_iterator<_Atom, Iterator>  _Iterator;
//
//		_Point cm = 0; unsigned count = 0; real_t mass = 0.;
//
//		_Iterator it__ = molecule->make_iterator(it);
//		_Iterator ite__ = molecule->make_iterator(ite);
//		for (; it__!=ite__; ++it__)
//		{
//			real_t mass__ = (*it__).mass;
//			cm += (*it__).X * mass__;
//			mass += mass__;
//		}
//		cm *= (1. / mass);
//		return cm;
//	}
//
//	template <typename Molecule>
//	INLINE const vdense_<Molecule::dimension, real_t>
//	calculate(_I2T<MASS_CENTER_>, const Molecule *molecule)
//	{
//		unsigned count = molecule->count(_I2T<ATOM_>());
//		return calculate(_I2T<MASS_CENTER_>(), molecule,
//			range_iterator(0), range_iterator(count));
//	}

	/**
	* @brief calculate determinant of inertia center
	* @return the determinant
	*/
	template <typename Molecule, typename Iterator>
	INLINE typename Molecule::real_type
	calculate(_I2T<INERTIA_TENSOR_DET_>, const Molecule *molecule,
		Iterator it, Iterator ite)
	{
		const unsigned dim = Molecule::dimension;
		typedef typename Molecule::real_type           _Real;
		typedef vdense_<dim, _Real>                    _Point;
		typedef typename Molecule::atom_type           _Atom;
		typedef const_array_iterator<_Atom, Iterator>  _Iterator;
		typedef mdense_<dim, dim, _Real>               _m3x3dense;

//		_Point cm = calculate(_I2T<MASS_CENTER_>(), molecule, it, ite);
		_Point cm = calculate(MASS_CENTER, molecule->get(ATOM), it, ite);
		_m3x3dense m; _Point X; _Real mass;

		_Iterator it__ = molecule->make_iterator(it);
		_Iterator ite__ = molecule->make_iterator(ite);
		for (; it__!=ite__; ++it__)
		{
			X = (*it__).X - cm;
			mass = (*it__).mass;
			m[0][0] += mass * (sqr(X[1]) + sqr(X[2]));
			m[1][1] += mass * (sqr(X[2]) + sqr(X[0]));
			m[2][2] += mass * (sqr(X[0]) + sqr(X[1]));
			m[0][1] += mass * (-X[0] * X[1]);
			m[0][2] += mass * (-X[0] * X[2]);
			m[1][2] += mass * (-X[1] * X[2]);
		}
		m[1][0] = m[0][1]; m[2][0] = m[0][2]; m[2][1] = m[1][2];

		_Real det = determinant(m);
		assert(_GT(det, 0.));
		return det;
	}

	template <typename Molecule>
	INLINE typename Molecule::real_type
	calculate(_I2T<INERTIA_TENSOR_DET_>, const Molecule *molecule)
	{
		unsigned count = molecule->count(_I2T<ATOM_>());
		return calculate(_I2T<INERTIA_TENSOR_DET_>(), molecule,
			range_iterator(0), range_iterator(count));
	}

	/**
	* @brief calculate full charge of molecule
	* @return the charge of molecule
	*/
	template <typename Molecule, typename Iterator>
	INLINE typename Molecule::real_type
	calculate(_I2T<CHARGE_>, const Molecule *molecule,
		Iterator it, Iterator ite)
	{
		const unsigned dim = Molecule::dimension;
		typedef typename Molecule::real_type           _Real;
		typedef vdense_<dim, _Real>                    _Point;
		typedef typename Molecule::atom_type           _Atom;
		typedef const_array_iterator<_Atom, Iterator>  _Iterator;

		_Real charge = 0;
		_Iterator it__ = molecule->make_iterator(it);
		_Iterator ite__ = molecule->make_iterator(ite);
		for (; it__!=ite__; ++it__) charge += (*it__).charge;
		return charge / SQRT_ELECTRIC_FACTOR;
	}

	template <typename Molecule>
	INLINE typename Molecule::real_type
	calculate(_I2T<CHARGE_>, const Molecule *molecule)
	{
		unsigned count = molecule->count(_I2T<ATOM_>());
		return calculate(_I2T<CHARGE_>(), molecule,
			range_iterator(0), range_iterator(count));
	}

#ifdef USE_FFTW_LIBRARY

	/**
	*  Рассчитывает потенциал дальнего кулона на точках сетки PPPM методом,
	*  аппроксимирует силы на атомах молекулы с помощью ближайших точек сетки,
	*  рассчитывает кулоновскую энергию.
	* @param molecule[in,out] молекула (комплекс молекул)
	* @param h - шаг сетки
	* @param smooth_radius - радиус сглаживания заряда атомов
	* @return энергия дальнего кулона (дополнительно корректирует силы на атомах)
	*/
// 	template <class _COMPLEX>
// 	INLINE _Real calculate(_I2T<COULFAR_>, _COMPLEX *molecule,
// 		_Real h=CHARGE_MESH_STEP,
// 		_Real smooth_radius=CHARGE_SMOOTH_RADIUS)
// 	{
// 		typedef typename _COMPLEX::atom_iterator iterator;
// 		typedef std::complex<_Real> complex_;
// 		typedef vdense_<4, _Real> v4dense_;
// 		typedef vdense_<4, complex_> c4dense_;
// 		typedef Mesh_<3>::index_type _Index;
//
// 		typedef Dfft<3, _Real, complex_ > Dfft_forward;
// 			// при прямом фурье преобразовании рассчитываем только фурье образ заряда
// 		typedef Dfft<3, c4dense_, v4dense_ > Dfft_backward;
// 			// при обратном фурье преобразовании рассчитываем как энергию, так и силы
//
// 		Box box = calculate(_I2T<BOX_>(), molecule);
// 		const Stencil_<3, 0> stencil(h, smooth_radius);
// 		Mesh_<3> mesh(box, h, smooth_radius);
// 		unsigned meshsz = mesh.size(); // размер всех сеток
//
// 		std::vector<_Real> netg(meshsz); // сетка для заряда
// 	#ifdef PPPM_DEBUG
// 		std::vector<_Real> netg__(meshsz); // сетка для заряда (тестовый вариант)
// 		// Удаление из системы излишнего заряда.
// 		// Фокус только для теста, так как в расчетах дает неверный результат.
// 		_Real q = 0.; unsigned n = 0;
// 		for (iterator it=molecule->begin(_I2T<ATOMS_>()),
// 			ite=molecule->end(_I2T<ATOMS_>()); it!=ite; ++it)
// 		{
// 			_Real charge = (*it).charge;
// 			q += charge; n++;
// 		}
// 		q /= n;
// 	#endif
//
// 		_Real Q = 0.;
// 	TIME_TESTING_START("mesh", 1)
// 		for (iterator it=molecule->begin(_I2T<ATOMS_>()),
// 			ite=molecule->end(_I2T<ATOMS_>()); it!=ite; ++it)
// 		{
// 			_Real charge = (*it).charge // закрытие после #ifdef
// 	#ifdef PPPM_DEBUG
// 			- q;
// 			mesh.insert(&netg__[0], stencil, (*it).X, charge);
// 	#endif
// 			; // закрытие оператора перед #ifdef PPPM_DEBUG
// 			mesh.insert(&netg[0], stencil, (*it).X, charge);
// 			Q += charge;
// 		}
// 	TIME_TESTING_FINISH
//
// // 		if (std::abs(Q) > M_SAFETY)
// // 		{
// // 			_S msg = _S("\n[USING ERROR] Full charge of system is ")
// // 				+ make_string(Q / SQRT_ELECTRIC_FACTOR) + _S("\n")
// // 				+ _S("   You must provide full Q == 0 to use PPPM method\n");
// // 			PRINT_MESSAGE(msg);
// // 		}
//
// 		Mesh_<3>::index_type meshdim = mesh.dimension();
// 		Dfft_forward  dfft_forward(meshdim);
// 		unsigned csz = dfft_forward.osize();
// 		unsigned rsz = dfft_forward.isize();
//
// 		std::vector<complex_> netgk(csz);
// 		// рассчитаем фурье образ заряда
// 	TIME_TESTING_START("forward", 1)
// 		dfft_forward(&netgk[0], &netg[0]);
// 	TIME_TESTING_FINISH
//
// 		Dfft_backward dfft_backward(mesh.dimension());
// 		std::vector<c4dense_> netkfp(csz);
// 		std::vector<v4dense_> netfp(rsz);
//
// 		netgk[0].real() = 0.;
// 		netgk[0].imag() = 0.;
// 			// Данный метод определен только для заряда = 0,
// 			// а точка G[0] и есть пространственный интеграл по плотности заряда
//
// 		Potential<COUL_> coulomb; _Point k; complex_ u;
// 		_Real klen; _Real coef = 1. / h;
// 	TIME_TESTING_START("extract", 1)
// 		for (unsigned i=1; i<csz; i++) // исключая K=0
// 		{
// 			dfft_forward.extract_k(k, i);
// 			klen = coef * sqrt(scalar_product(k, k));  // ускорить через k**2
// 			u = netgk[i] * coulomb.fourier_far3(klen);
// 			netkfp[i][0] = u * complex_(0., k[0]);
// 			netkfp[i][1] = u * complex_(0., k[1]);
// 			netkfp[i][2] = u * complex_(0., k[2]);
// 			netkfp[i][3] = u;
// 		}
// 	TIME_TESTING_FINISH
//
// 	TIME_TESTING_START("backward", 1)
// 		dfft_backward(&netkfp[0], &netfp[0]);
// 	TIME_TESTING_FINISH
//
// 		_Real norm = sqr(dfft_forward.norm()) / (cube(h) * SQRT_ELECTRIC_FACTOR);
// 			// нормировка включает 1/h**3 множитель от отображения фурье образа кулона
// 			// на дискретную сетку с шагом h
// 	TIME_TESTING_START("norm", 1)
// 		for (unsigned i=0; i<rsz; i++)
// 		{
// 			netfp[i][0] *= norm;
// 			netfp[i][1] *= norm;
// 			netfp[i][2] *= norm;
// 			netfp[i][3] *= norm;
// 		}
// 	TIME_TESTING_FINISH
//
// 		_Point X__; // координата, относительно начала той ячейки, с которой связывается X
// 		unsigned neig[8]; // линейные индексы узлов ячейки, куда попадает X
// 		v4dense_ F__; // силы [0..2] и потенциал [3]
// 		_Real energy = 0.;
// 	TIME_TESTING_START("linear", 1)
// 		for (iterator it=molecule->begin(_I2T<ATOMS_>()),
// 			ite=molecule->end(_I2T<ATOMS_>()); it!=ite; ++it)
// 		{
// 			_Real charge = (*it).charge / SQRT_ELECTRIC_FACTOR;
// 			mesh.get_neighbours((unsigned *)neig, X__, (*it).X);
//
// 			linear_interpolation(4, &F__[0], X__[0], X__[1], X__[2],
// 				&netfp[neig[0]][0], &netfp[neig[1]][0],
// 				&netfp[neig[2]][0], &netfp[neig[3]][0],
// 				&netfp[neig[4]][0], &netfp[neig[5]][0],
// 				&netfp[neig[6]][0], &netfp[neig[7]][0]);
//
// 			_Point &F = (*it).F;
// 			F[0] += F__[0] * charge;
// 			F[1] += F__[1] * charge;
// 			F[2] += F__[2] * charge;
// 				// меняю знак только для сил на атомах, согласно F = -dU/dX
// 				// само поле mesh_fp - есть поле производных, оно не нужно далее,
// 				// потому игнорируем смену знака
//
// 		#ifdef PPPM_DEBUG
// 			F[0] = F__[0] * charge;
// 			F[1] = F__[1] * charge;
// 			F[2] = F__[2] * charge;
// 		#endif
// 			energy += charge * F__[3];
// 		}
// 	TIME_TESTING_FINISH
//
// 	#ifdef PPPM_DEBUG
// 	//----------------------------------------------------------------------------
// 	// Эта часть функции используется только для тестирования, в нормальном режиме
// 	// она автоматически отключена (PPPM_DEBUG не определена)
// 	//----------------------------------------------------------------------------
// 		char buf[120]; _Real r, p; _Point rr;
//
// 		PRINT_LINE(97,'-');
// 		sprintf(buf, "node      direct dV/dX [kJ/mol A]   &   V[kJ/mol]   vs.   fourie dV/dX [kJ/mol A]   &   V[kJ/mol]");
// 		std::cout << buf << std::endl;
// 		PRINT_LINE(97,'-');
//
// 		for (int i=0; i<meshdim[0]; i++)
// 		{
// 			// Выведем данные только по главной диагонали, что дает разумный компромисс
// 			// между скоростью теста и представлением результатов. Данные можно посмотреть
// 			// в любом редакторе, который отрисовывает одномерные графики (например, Excel)
// 			Mesh_<3>::index_type ndx(i, i, i);
// 			unsigned n = mesh[ndx];
// 			v4dense_ v = v4dense_(0., 0., 0., 0.);
// 			for (int i__=0; i__<meshdim[0]; i__++)
// 			for (int j__=0; j__<meshdim[1]; j__++)
// 			for (int k__=0; k__<meshdim[2]; k__++)
// 			{
// 				Mesh_<3>::index_type ndx__(i__, j__, k__);
// 				_Real charge__ = netg__[mesh[ndx__]];
// 				if (charge__ == 0.) continue;
//
// 			#ifdef PPPM_BOUNDS_DEBUG
// 				for (int ii=-1; ii<2; ii++)
// 				for (int jj=-1; jj<2; jj++)
// 				for (int kk=-1; kk<2; kk++)
// 				{
// 					// Сделаем суммирование от точек в соседних ячейках, чтобы снять
// 					// влияние границ на сопоставление результатов
// 					rr[0] = h * (i - i__ - meshdim[0] * ii);
// 					rr[1] = h * (i - j__ - meshdim[1] * jj);
// 					rr[2] = h * (i - k__ - meshdim[2] * kk);
// 					r = h * sqrt(
// 							sqr(i - i__ - meshdim[0] * ii)
// 						+ sqr(i - j__ - meshdim[1] * jj)
// 						+ sqr(i - k__ - meshdim[2] * kk)
// 					);
// 					v[3] += charge__ * coulomb.far(r); // необходимо добавлять потенциал от r=0
//
// 					if (r == 0.) continue;
// 					p = coulomb.far1(r) / r;
// 					v[0] += charge__ * p * rr[0];
// 					v[1] += charge__ * p * rr[1];
// 					v[2] += charge__ * p * rr[2];
// 				}
// 			#else
// 				rr[0] = h * (i - i__);
// 				rr[1] = h * (i - j__);
// 				rr[2] = h * (i - k__);
// 				r = h * sqrt( sqr(i - i__) + sqr(i - j__) + sqr(i - k__));
// 				v[3] += charge__ * coulomb.far(r);
//
// 				if (r == 0.) continue;
// 				p = coulomb.far1(r) / r;
// 				v[0] += charge__ * p * rr[0];
// 				v[1] += charge__ * p * rr[1];
// 				v[2] += charge__ * p * rr[2];
// 			#endif
// 			}
// 			v[0] /= SQRT_ELECTRIC_FACTOR;
// 			v[1] /= SQRT_ELECTRIC_FACTOR;
// 			v[2] /= SQRT_ELECTRIC_FACTOR;
// 			v[3] /= SQRT_ELECTRIC_FACTOR;
// 			sprintf(buf, "%3d [%10.3e %10.3e %10.3e] %10.3e   [%10.3e %10.3e %10.3e] %10.3e",
// 				i, v[0], v[1], v[2], v[3], netfp[n][0], netfp[n][1], netfp[n][2], netfp[n][3]);
// 			std::cout << buf << std::endl;
// 		}
// 		PRINT_LINE(97,'-');
// 		sprintf(buf, "atoms     direct dV/dX [kJ/mol A]   vs.   fourie dV/dX [kJ/mol A]");
// 		std::cout << buf << std::endl;
// 		PRINT_LINE(97,'-');
//
// 		// Так как размер ребра ячейки и произведение h на размерность ячейки не совпадают,
// 		// в силу того, что h не корректируется под размерность ячейки, то неизбежны
// 		// расхождения между тестами на границах диапазона. Возможно отсутствие совпадения
// 		// значений fftw на границах.
//
// 		unsigned i = 0;
// 		for (iterator it=molecule->begin(_I2T<ATOMS_>()),
// 			ite=molecule->end(_I2T<ATOMS_>()); it!=ite; ++it, i++)
// 		{
// 			const _Point &F = (*it).F;
// 			const _Point &X = (*it).X;
// 			_Real charge = (*it).charge / ELECTRIC_FACTOR;
// 			_Point testF = 0.;
// 			for (iterator it__=molecule->begin(_I2T<ATOMS_>()); it__!=ite; ++it__)
// 			{
// 				const _Point &X__ = (*it__).X;
// 				_Point R = X - X__;
// 				r = distance1(X, X__);
// 				if (r == 0.) continue;
//
// 				_Real charge__ = (*it__).charge;
// 				p = charge__ * charge * coulomb.far1(r) / r;
// 				testF[0] += p * R[0];
// 				testF[1] += p * R[1];
// 				testF[2] += p * R[2];
// 			}
// 			sprintf(buf, "%3d [%10.3e %10.3e %10.3e]    [%10.3e %10.3e %10.3e]",
// 				i, testF[0], testF[1], testF[2], F[0], F[1], F[2]);
// 			std::cout << buf << std::endl;
// 		}
// 		PRINT_LINE(97,'-');
// 	//----------------------------------------------------------------------------
// 	//               Конец тестирующей части (PPPM_DEBUG)
// 	//----------------------------------------------------------------------------
// 	#endif
//
// 		return energy;
// 	}

#endif

	/**
	* @brief calculate full charge of molecule
	* @return the charge of molecule
	*/
	template <typename Molecule>
	INLINE typename Molecule::real_type
	calculate(_I2T<SOLID_TdS_>, const Molecule *molecule,
		typename Molecule::real_type kT=KT(300.0))
	{
		typedef typename Molecule::real_type  _Real;

		_Real mass = calculate(_I2T<MASS_>(), molecule);
		_Real TdS = 0.;
		if (molecule->count(_I2T<ATOM_>()) == 1)
		{
			TdS = kT * (1.5 * log(mult2(M_PI) * kT * mass)
				+ log(INTEGRAL__dXYZ) );
		}
		else
		{
			_Real detI = calculate(_I2T<INERTIA_TENSOR_DET_>(), molecule);
			TdS = kT * ( log(cube(mult2(M_PI) * kT))
				+ log(INTEGRAL__dXYZ * INTEGRAL__dOMEGA)
				+ 0.5 * log(cube(mass) * detI) );
		}
		return TdS;
	}

	/**
	* @brief calculates the molecule effecitve electrostatic size
	* @note "Sigalov2.pdf" Eq. (23)
	* @return determinant of effecitve electrostatic size
	*/
	template <typename Molecule>
	INLINE typename Molecule::real_type
	calculate(_I2T<A_DET_>, const Molecule *molecule)
	{
		const unsigned dim = Molecule::dimension;
		typedef typename Molecule::real_type                 _Real;
		typedef vdense_<dim, _Real>                          _Point;
		typedef typename Molecule::atom_type                 _Atom;
		typedef const_array_iterator<_Atom, range_iterator>  _Iterator;

		unsigned count = molecule->count(_I2T<ATOM_>());

		// calculates Sigalov center mass & Sigalov full mass
		_Point scm = 0.; _Real mass = 0.;

		_Iterator start = molecule->make_iterator(range_iterator(0));
		_Iterator end = molecule->make_iterator(range_iterator(count));
		for (_Iterator it=start; it!=end; ++it)
		{
			_Real mass__ = cube((*it).radius);
			scm += (*it).X * mass__;
			mass += mass__;
		}
		scm *= (1. / mass);

		_Real I11 = 0;
		_Real I12 = 0;
		_Real I13 = 0;
		_Real I22 = 0;
		_Real I23 = 0;
		_Real I33 = 0;
			// Inertia tensor components

		for (_Iterator it=start; it!=end; ++it)
		{
			const _Atom &atom = *it;
			_Real x = atom.X[0] - scm[0];
			_Real y = atom.X[1] - scm[1];
			_Real z = atom.X[2] - scm[2];
			_Real radius = atom.radius;
			_Real radius2__ = sqr(radius);
			_Real m = radius2__ * radius;
			_Real x2 = x * x;
			_Real y2 = y * y;
			_Real z2 = z * z;
			_Real a__ = (2./5.) * radius2__;
			I11 += m * (y2 + z2 + a__);
			I22 += m * (x2 + z2 + a__);
			I33 += m * (x2 + y2 + a__);
			I12 -= m * x * y;
			I13 -= m * x * z;
			I23 -= m * y * z;
		}

		_Real det = I11*I22*I33 + 2*I12*I23*I13 - I11*I23*I23 - I22*I13*I13 - I33*I12*I12;
		return sqrt((5./2.) / mass) * pow(det, 1./6.);
	}

}
#endif
