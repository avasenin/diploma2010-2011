#ifndef _COULOMB_PARAMS__0077A726_CBE6_5c54_73A2_F9434C9C0000__H
#define _COULOMB_PARAMS__0077A726_CBE6_5c54_73A2_F9434C9C0000__H

#pragma once
namespace molkern
{
  class CoulombParams
  {

    /**
    * @brief Predefined Rappe and Goddard elements && gaussians
    */

    static CoulombParams* s_instance;
    Basis_ m_basis;
  public:
    CoulombParams(Basis_ basis)
    {
      m_basis = basis;
    }
    static void build(Basis_ basis)
    {
        s_instance = new CoulombParams(basis);
    }
    static CoulombParams* instance()
    {
      return s_instance;
    }
    TEMPLATE_HEADER
    template <typename _Atom>
    inline _E(real_t) get(_Atom atomA, _Atom atomB)
    {
        const sGTO &sgto_from_first_atom= m_basis[make_string(atomA.atomdata->name)];
        const sGTO &sgto_from_second_atom= m_basis[make_string(atomB.atomdata->name)];
        double r = distance1(atomA.atomdata->X, atomB.atomdata->X);
        return m_basis.coulomb_integral(sgto_from_first_atom, sgto_from_second_atom, r);
    }
  };
}
#endif
