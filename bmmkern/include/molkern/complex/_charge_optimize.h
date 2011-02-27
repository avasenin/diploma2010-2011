#ifndef _CHARGE_OPTIMIZE__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H
#define _CHARGE_OPTIMIZE__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H

#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_io.hpp"

#include "molkern/__moldefs.h"
#include "molkern/complex/_coulomb_params.h"
#include <vector>
namespace molkern
{
  using namespace prgkern;
  /**
  * @brief geometry optimizer for molecule (or complex)
  */
  template <int OPTIMIZER_TYPE>
  class ChargeOptimizer_ : protected Minimizer_<OPTIMIZER_TYPE>
  {
    struct ChargeOptimizerParams
    {
      void* complex;
    };
  private:
    std::vector<real_t> m_aCharges, m_aDerivatives;
    mdense_<UNLIMITED_,UNLIMITED_,real_t> m_coulomb;
    Basis_ m_basis;
    public:
    template <typename LPComplex, typename Param>
    real_t operator()(LPComplex *complex, const Param &param)
    {
      typedef typename LPComplex::atom_type _Atom;
      real_t energy = 0.;
      unsigned number_of_atoms = complex->count(ATOM);
      m_aCharges.resize(number_of_atoms);
      m_aDerivatives.resize(number_of_atoms);
      _Atom *atoms = complex->get(ATOM);
      for (unsigned i=0; i < number_of_atoms; i++)
      {
        std::string atom_symbol = make_string(atoms[i].atomdata->name);
        double gaussian_exponent = RappleGoddardParams::instance()->find(atom_symbol).gaussian_exponent;
        m_basis.insert(Basis_::value_type(atom_symbol, gaussian_exponent));
      }
      CoulombParams::build(m_basis);
      ChargeOptimizerParams params = {(void*) complex};
      energy = Minimizer_<OPTIMIZER_TYPE>::operator()(
        &dU__dQ<LPComplex>, (void*)&params,
        number_of_atoms - 1, &m_aCharges[0], &m_aDerivatives[0], param.maxiter,
        param.stpmin, param.stpmax,
        param.maxfev, param.maxhalt,
        param.wolfe1, param.wolfe2,
        param.xtol, param.ftol, param.gtol,
        param.m, param.steep
      );
      for (unsigned i=0; i < number_of_atoms; i++)
      {
        PRINT_MESSAGE(make_string("i=%d CHARGE OPTIMIZE=%.5e\n",i,atoms[i].charge));
      }
      return energy;
    }

    template <typename LPComplex>
    static real_t dU__dQ(unsigned n, const real_t *charges, real_t *derivatives, void *params)
    {
      typedef typename LPComplex::atom_type _Atom;
      ChargeOptimizerParams* casted_params = (ChargeOptimizerParams*) params;
      LPComplex* complex = (LPComplex*)(casted_params->complex);
      _Atom *atoms = complex->get(ATOM);

      complex->read(CHARGE, charges, atoms);
      real_t energy = complex->dU__dQ();
      complex->write(DCHARGE, derivatives, atoms);

      return energy;
    }
  };

}
#endif
