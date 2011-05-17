#ifndef _RAPPLE_GODDARD_PARAMS__0077A726_CBE6_5c54_73A2_F9434C9C0000__H
#define _RAPPLE_GODDARD_PARAMS__0077A726_CBE6_5c54_73A2_F9434C9C0000__H

#pragma once
namespace molkern
{
  //using namespace prgkern;
  const double EV_2_HARTREE          = 3.67493245e-2;
  struct Nucleus
  {
    std::string symbol;
    int nuclear_charge;
    int formal_charge;
    double electronegativity;
    double hardness;
    double gaussian_exponent;
  };
  const Nucleus RAPPLE_GODDARD_PARAMS[] =
  {
    { "UNKNOW",  0, 0, 0, 0, 0.},
    { "H",  1, 0,  4.528 * EV_2_HARTREE, 13.890 * EV_2_HARTREE, 0.534337523756312 },
    { "LI", 3, 0,  3.006 * EV_2_HARTREE,  4.772 * EV_2_HARTREE, 0.166838519142176 },
    { "C",  6, 0,  5.343 * EV_2_HARTREE, 10.126 * EV_2_HARTREE, 0.206883838259186 },
    { "N",  7, 0,  7.139 * EV_2_HARTREE, 12.844 * EV_2_HARTREE, 0.221439796025873 },
    { "O",  8, 0,  8.741 * EV_2_HARTREE, 13.364 * EV_2_HARTREE, 0.223967308625516 },
    { "F",  9, 0, 10.874 * EV_2_HARTREE, 14.948 * EV_2_HARTREE, 0.231257590182828 },
    { "NA", 11, 0,  2.843 * EV_2_HARTREE,  4.592 * EV_2_HARTREE, 0.095892938712585 },
    { "SI", 14, 0,  4.168 * EV_2_HARTREE,  6.974 * EV_2_HARTREE, 0.105219608142377 },
    { "P",  15, 0,  5.463 * EV_2_HARTREE,  8.000 * EV_2_HARTREE, 0.108476721661715 },
    { "S",  16, 0,  6.084 * EV_2_HARTREE, 10.660 * EV_2_HARTREE, 0.115618357843499 },
    { "CL", 17, 0,  8.564 * EV_2_HARTREE,  9.892 * EV_2_HARTREE, 0.113714050615107 },
    { "K",  19, 0,  2.421 * EV_2_HARTREE,  3.840 * EV_2_HARTREE, 0.060223294377778 },
    { "BR", 35, 0,  7.790 * EV_2_HARTREE,  8.850 * EV_2_HARTREE, 0.070087547802259 },
    { "RB", 37, 0,  2.331 * EV_2_HARTREE,  3.692 * EV_2_HARTREE, 0.041999054745368 },
    { "I",  53, 0,  6.822 * EV_2_HARTREE,  7.524 * EV_2_HARTREE, 0.068562697575073 },
    { "CS", 55, 0,  2.183 * EV_2_HARTREE,  3.422 * EV_2_HARTREE, 0.030719481189777 }
  };
  class RappleGoddardParams
  {

    /**
    * @brief Predefined Rappe and Goddard elements && gaussians
    */

    static RappleGoddardParams* s_instance;
  public:
    static RappleGoddardParams* instance()
    {
      if (!s_instance)
      {
        s_instance = new RappleGoddardParams();
      }
      return s_instance;
    }
    inline const Nucleus find(const std::string &name)
    {
      unsigned number_of_nucleus_in_base = sizeof(RAPPLE_GODDARD_PARAMS) / sizeof(Nucleus);
      for (unsigned i=0; i<=number_of_nucleus_in_base; i++)
          if (trim(make_string(RAPPLE_GODDARD_PARAMS[i].symbol)) == trim(name)) return RAPPLE_GODDARD_PARAMS[i];
      return RAPPLE_GODDARD_PARAMS[0];
    }
  };
}
#endif
