#ifndef _CHARGE_OPTIMIZE__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H
#define _CHARGE_OPTIMIZE__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H

#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_io.hpp"

#include "molkern/__moldefs.h"

namespace molkern
{
  using namespace prgkern;

  /**
  * @brief geometry optimizer for molecule (or complex)
  */
  template <int OPTIMIZER_TYPE>
  class Optimizer_ : protected Minimizer_<OPTIMIZER_TYPE>
  {
private:
std::vector<double> m_aCharges, m_aDerivatives;
  public:

    template <typename LPComplex, typename Param>
    real_t operator()(LPComplex *complex, const Param &param)
    {
      real_t energy = 0.;
      unsigned number_of_atoms = complex->count(ATOM);
      m_aCharges.resize(number_of_atoms);
      m_aDerivatives.resize(number_of_atoms);
      energy = Minimizer_<OPTIMIZER_TYPE>::operator()(
        &dU__dQ<LPComplex>, (void*)complex,
        number_of_atoms, &m_aCharges[0], &m_aDerivatives[0], param.maxiter,
        param.stpmin, param.stpmax,
        param.maxfev, param.maxhalt,
        param.wolfe1, param.wolfe2,
        param.xtol, param.ftol, param.gtol,
        param.m, param.steep
      );
      _Atom *atoms = complex->get(ATOM);
      for (int i=0; i < number_of_atoms; i++)
      {
        PRINT_MESSAGE("i=%d CHARGE OPTIMIZE=%.5\n",i,atoms[i].charge);
      }
      return energy;
    }

    template <typename LPComplex>
    static real_t dU__dQ(unsigned n, const real_t *charges, real_t *derivatives, void *param)
    {
      typedef typename LPComplex::atom_type _Atom;
      LPComplex* complex = (LPComplex*)param;
      _Atom *atoms = complex->get(ATOM);

      complex->read(POSITION, x, atoms);
      real_t energy = complex->dU__dX();
      complex->write(GRADIENT, g, atoms);

      return energy;
    }
  };

}
#endif
