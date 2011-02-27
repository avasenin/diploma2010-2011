#ifndef _COULOMB_PARAMS__0077A726_CBE6_5c54_73A2_F9434C9C0000__H
#define _COULOMB_PARAMS__0077A726_CBE6_5c54_73A2_F9434C9C0000__H

#pragma once
namespace molkern
{
  struct sGTO
  {
    sGTO(double exponent=0.) : zeta(exponent) {}
    double zeta;
  };
  /// Threshold for calculating Coulomb integrals
  const double COULOMB_INTEGRAL_THRESHOLD = 1.0e-9;
  const double LOG_COULOMB_INTEGRAL_THRESHOLD = log(COULOMB_INTEGRAL_THRESHOLD);

  /// Threshold for calculating overlap integrals
  const double OVERLAP_INTEGRAL_THRESHOLD = 1.0e-9;
  const double SQR_OVERLAP_INTEGRAL_THRESHOLD = sqr(OVERLAP_INTEGRAL_THRESHOLD);

  class Basis_ : public std::map<std::string, sGTO>
  {
    typedef std::map<std::string, sGTO> _Base;

    double Coulomn_Integral_Distance_Threshold;
    double Overlap_Integral_Distance_Threshold;
    double smallest_gaussian_exponent;

  public:
    typedef std::map<std::string, sGTO>::value_type  value_type;

    Basis_() : smallest_gaussian_exponent(1.e20) {}

    void insert(const value_type &val)
    {
      _Base::insert(val);
      double zeta = val.second.zeta;
      if (smallest_gaussian_exponent > zeta)
      {
        smallest_gaussian_exponent = zeta;
        Coulomn_Integral_Distance_Threshold =
          2 * sqrt(-LOG_COULOMB_INTEGRAL_THRESHOLD / smallest_gaussian_exponent);
        Overlap_Integral_Distance_Threshold =
          sqrt(
            log( (M_PI / cube(2 * smallest_gaussian_exponent))
                / SQR_OVERLAP_INTEGRAL_THRESHOLD
            ) / smallest_gaussian_exponent
          );
      }
    }

    /**
    * @brief Computes the 2-center Coulomb integral over Gaussian-type s-orbitals analytically.
    * @note Reference: T. Helgaker, P. Jorgensen, J. Olsen, Molecular Electronic Structure Theory
    *   Wiley, NY, 2000, Equations (9.7.21) and (9.8.23)
    * @param a,b atom names
    * @param R internuclear distance in atomic units (bohr)
    * @return the value of the Coulomb potential energy integral
    */
    double coulomb_integral(sGTO a, sGTO b, double R)
    {
      if (R > Coulomn_Integral_Distance_Threshold) return 1. / R;
      double p = sqrt(a.zeta * b.zeta / (a.zeta + b.zeta));
      return erf(p * R) / R;
    }

    /**
    * @brief Computes overlap integral analytically over s-type GTOs
    * @note Reference: T. Helgaker, P. Jorgensen, J. Olsen, Molecular Electronic Structure Theory
    *    Wiley, NY, 2000, Equation (9.2.41)
    * @param a Gaussian exponent of first atom in atomic units (inverse squared Bohr)
    * @param b Gaussian exponent of second atom in atomic units (inverse squared Bohr)
    * @param R internuclear distance in atomic units (bohr)
    */
    double overlap_integral(sGTO a, sGTO b, double R)
    {
      if (R > Overlap_Integral_Distance_Threshold) return 0.;
      double p = a.zeta + b.zeta;
      double q = a.zeta * b.zeta / p;
      return pow(4 * q / p, 0.75) * exp(-q * R * R);
    }
  };

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
    template <typename _Atom>
    inline _E(real_t) get(_Atom atomA, _Atom atomB)
    {
        double r = distance1(atomA.atomdata->X, atomB.atomdata->X);
        if (0 == r)
        {
          return RappleGoddardParams::instance()->find(make_string(atomA.atomdata->name)).hardness;
        }
        const sGTO &sgto_from_first_atom= m_basis[make_string(atomA.atomdata->name)];
        const sGTO &sgto_from_second_atom= m_basis[make_string(atomB.atomdata->name)];
        return m_basis.coulomb_integral(sgto_from_first_atom, sgto_from_second_atom, r);
    }
  };
}
#endif
