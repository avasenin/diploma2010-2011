#ifndef _CHARGE_OPTIMIZE__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H
#define _CHARGE_OPTIMIZE__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H

#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_io.hpp"

#include "molkern/__moldefs.h"
#include "molkern/complex/_coulomb_params.h"
#include "molkern/forcefield/_charge_optimize_params.h"
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
  private:
    real_t *m_variables, *m_derivatives;
    public:
    template <typename LPComplex, typename Param>
    real_t operator()(LPComplex *complex, const Param &param)
    {
      typedef typename LPComplex::atom_type _Atom;
      real_t energy = 0.;
      unsigned number_of_atoms = complex->count(ATOM);
      m_variables = new real_t[2 * number_of_atoms];
      memset(m_variables, 0, 2 * number_of_atoms * sizeof(real_t));
      m_derivatives = new real_t[2 * number_of_atoms];
      memset(m_derivatives, 0, 2 * number_of_atoms * sizeof(real_t));
      _Atom *atoms = complex->get(ATOM);
      Basis_ basis;
      for (unsigned i=0; i < number_of_atoms; i++)
      {
        std::string atom_symbol = make_string(atoms[i].atomdata->name);
        double gaussian_exponent = RappleGoddardParams::instance()->find(atom_symbol).gaussian_exponent;
        basis.insert(Basis_::value_type(atom_symbol, gaussian_exponent));
      }
      //init atoms charges
      /*real_t charges[2*number_of_atoms];
      memset(charges, 0, 2 * number_of_atoms * sizeof(real_t));
*/    for (int i = 0; i < number_of_atoms; i++)
      {
        //m_variables[i] = atoms[i].atomdata->charge;
	/*m_variables[0] = -0.896;
	m_variables[1] = 0.418;
	m_variables[2] = 0.478;
     */     //RappleGoddardParams::instance()->find(make_string(atoms[i].atomdata->name)).formal_charge;
      }
      CoulombParams::build(basis);
      ChargeOptimizeParams<LPComplex>* params = new ChargeOptimizeParams<LPComplex>(complex);
      //collapseEnumeration(charges, m_variables, params);
    //  for (unsigned i=0; i < number_of_atoms; i++)
     // {
      //  PRINT_MESSAGE(make_string("i=%d CHARGE OPTIMIZE=%.5e\n",i,atoms[i].charge));
     // }
      /*
      energy = Minimizer_<OPTIMIZER_TYPE>::operator()(
        &dU__dQ<LPComplex>, params,
        2 * number_of_atoms - 1, m_variables, m_derivatives, param.maxiter,
        param.stpmin, param.stpmax,
        param.maxfev, param.maxhalt,
        param.wolfe1, param.wolfe2,
        param.xtol, param.ftol, param.gtol,
        param.m, param.steep
      ); */

      for (int k=0; k < 200; k++)
      {
        for (int i=0; i < 1000; i++)
        {
          modeling(m_variables, 0.0001, params);
        }
        PRINT_MESSAGE(make_string("%.4f %.4f %.4f", m_variables[0], m_variables[1], m_variables[2]));
        //modeling(m_variables, 0.0001, params);
      }

      modeling(m_variables, 0.00001, params);
      modeling(m_variables, 0.00001, params);
      modeling(m_variables, 0.00001, params);
      delete params;
      delete m_variables;
      delete m_derivatives;
      for (unsigned i=0; i < number_of_atoms; i++)
      {
        PRINT_MESSAGE(make_string("i=%d CHARGE OPTIMIZE=%.5e\n",i,atoms[i].charge));
      }
      return energy;
    }

    template <typename TComplex>
    static void collapseEnumeration(const real_t *source, real_t *dest, ChargeOptimizeParams<TComplex>* params)
    {
        typedef typename TComplex::molecule_type _Molecule;
        typedef typename _Molecule::archetype_type _Archetype;
        TComplex* complex = params->getComplex();
        std::vector<_Molecule*> molecules = complex->get(MOLECULE);
        unsigned numberOfMolecules = molecules.size();
        int destIndex = 0;
        int sourceIndex = 0;
        real_t* q = dest;
        for (unsigned moleculeIndex = 0; moleculeIndex < numberOfMolecules; moleculeIndex++)
        {
          const _Archetype* currentArchetype = molecules[moleculeIndex]->archetype();
          unsigned numberOfAtomsInMolecule = currentArchetype->count(ATOM);
          for (unsigned atomIndex = 0; atomIndex < numberOfAtomsInMolecule - 1; atomIndex++)
          {
            dest[destIndex++] = source[sourceIndex++];
          }
          sourceIndex++;
        }
        for (int i = 0; i < params->getNumberOfAtoms(); i++)
        {
          dest[destIndex++] = source[sourceIndex++];
        }
    }

    template <typename TComplex>
    static void expandEnumeration(const real_t *source, real_t *dest, ChargeOptimizeParams<TComplex>* params)
    {
      typedef typename TComplex::molecule_type _Molecule;
      typedef typename _Molecule::archetype_type _Archetype;
      TComplex* complex = params->getComplex();
      std::vector<_Molecule*> molecules = complex->get(MOLECULE);
      unsigned numberOfMolecules = molecules.size();
      int atomOffset = 0;
      int variableIndex = 0;
      real_t* q = dest;
      for (unsigned moleculeIndex = 0; moleculeIndex < numberOfMolecules; moleculeIndex++)
      {
        if (0 != moleculeIndex)
        {
          const _Archetype* prevArchetype = molecules[moleculeIndex - 1]->archetype();
          atomOffset += prevArchetype->count(ATOM);
        }
        real_t charge = 0.;
        const _Archetype* currentArchetype = molecules[moleculeIndex]->archetype();
        unsigned numberOfAtomsInMolecule = currentArchetype->count(ATOM);
        for (unsigned atomIndex = 0; atomIndex < numberOfAtomsInMolecule - 1; atomIndex++, variableIndex++)
        {
          q[atomIndex + atomOffset] = source[variableIndex];
          charge += source[variableIndex];
        }
        q[numberOfAtomsInMolecule + atomOffset - 1] = -charge;
      }
      real_t* p = dest + params->getNumberOfAtoms();
      memcpy(p, source +params->getNumberOfAtoms()-numberOfMolecules, params->getNumberOfAtoms() * sizeof(real_t));
    }

    template <typename TComplex>
    static real_t dU__dQ(unsigned n, const real_t *variablesCompressed, real_t *derivativesCompressed, void *uncastedParams)
    {
      typedef typename TComplex::atom_type _Atom;
      typedef typename TComplex::molecule_type _Molecule;
      typedef typename _Molecule::archetype_type _Archetype;

      ChargeOptimizeParams<TComplex>* params = (ChargeOptimizeParams<TComplex>*) uncastedParams;
      TComplex* complex = params->getComplex();
      real_t variables[2 * params->getNumberOfAtoms()];
      memset(variables, 0, 2 * params->getNumberOfAtoms() * sizeof(real_t));
      real_t* q = variables;
      real_t* p = variables+ params->getNumberOfAtoms();
      expandEnumeration(variablesCompressed, variables, params);

      real_t derivatives[2 * params->getNumberOfAtoms()];
      memset(derivatives, 0, 2 * params->getNumberOfAtoms() * sizeof(real_t));
      real_t* dq = derivatives;
      real_t* dp = derivatives + params->getNumberOfAtoms();

      const mdense_<UNLIMITED_, UNLIMITED_, real_t>& matrixJ = params->getMatrixJ();
      const mdense_<UNLIMITED_, UNLIMITED_, real_t>& matrixAKA = params->getMatrixAKA();
      _Atom *atoms = complex->get(ATOM);
      complex->read(CHARGE, variables, atoms);
      typedef typename _Molecule::archetype_type _Archetype;
      std::vector<_Molecule*> molecules = complex->get(MOLECULE);
      unsigned numberOfMolecules = molecules.size();
//calculating dh_dp and dh_dq
      for (unsigned i=0; i < params->getNumberOfAtoms(); i++)
      {
        dp[i] = 0;
        double electronegativity = RappleGoddardParams::instance()->find(make_string(atoms[i].atomdata->name)).electronegativity;
        dq[i] = electronegativity;
        for (unsigned j=0; j < params->getNumberOfAtoms(); j++)
        {
          dp[i] += matrixAKA(i,j)*p[j];
          dq[i] += matrixJ(i,j)*q[j];
        }
      }
//normalization
      unsigned atomOffset = 0;
      for (unsigned moleculeIndex = 0; moleculeIndex < numberOfMolecules; moleculeIndex++)
      {
        if (0 != moleculeIndex)
        {
          const _Archetype* prevArchetype = molecules[moleculeIndex - 1]->archetype();
          atomOffset += prevArchetype->count(ATOM);
        }
        const _Archetype* currentArchetype = molecules[moleculeIndex]->archetype();
        unsigned numberOfAtomsInMolecule = currentArchetype->count(ATOM);
        for (unsigned atomIndex = 0; atomIndex < numberOfAtomsInMolecule - 1; atomIndex++)
        {
          dq[atomIndex + atomOffset] -= dq[numberOfAtomsInMolecule + atomOffset - 1];
        }
      }

//energy
      real_t energy = 0.;
      //   1/2 P A K-1 AT P
      for (unsigned i=0; i < params->getNumberOfAtoms(); i++)
      {
        for (unsigned j=0; j < params->getNumberOfAtoms(); j++)
        {
          energy += (1. / 2) * p[i]*matrixAKA(i,j)*p[j];
        }
      }
      // (x,Q)
      for (unsigned i=0; i < params->getNumberOfAtoms(); i++)
      {
        double electronegativity = RappleGoddardParams::instance()->find(make_string(atoms[i].atomdata->name)).electronegativity;
        energy += electronegativity * q[i];
      }
      //  1/2 Q J Q
      for (unsigned i=0; i < params->getNumberOfAtoms(); i++)
      {
        for (unsigned j=0; j < params->getNumberOfAtoms(); j++)
        {
          energy += (1. / 2) * q[i]*matrixJ(i,j)*q[j];
        }
      }
      collapseEnumeration(derivatives, derivativesCompressed, params);
      return energy;
    }
//---------------------
    static void dpCalculator(unsigned atomIndexA, unsigned atomIndexB, unsigned atomOffset, unsigned numberOfAtomsInMolecule, const mdense_<UNLIMITED_, UNLIMITED_, real_t>& matrixAKA, real_t* dp_dt, real_t* numerator, real_t* denomirator)
    {
      if (atomIndexA == numberOfAtomsInMolecule)
        *numerator += matrixAKA(atomIndexA + atomOffset, atomIndexB + atomOffset) * dp_dt[atomIndexA];
      else
        *denomirator += matrixAKA(atomIndexA + atomOffset, atomIndexB + atomOffset);
    }
    template <typename TComplex>
    static void modeling(real_t *variables, real_t time, ChargeOptimizeParams<TComplex>* params)
    {
      typedef typename TComplex::atom_type _Atom;
      typedef typename TComplex::molecule_type _Molecule;
      typedef typename _Molecule::archetype_type _Archetype;
      typedef typename TComplex::bond_type _Bond;

      TComplex* complex = params->getComplex();
      real_t* q = variables;
      real_t* p = variables+ params->getNumberOfAtoms();

      real_t derivatives[2*params->getNumberOfAtoms()];
      real_t* dq_dt = derivatives;
      real_t* dp_dt = derivatives + params->getNumberOfAtoms();

      const mdense_<UNLIMITED_, UNLIMITED_, real_t>& matrixJ = params->getMatrixJ();
      const mdense_<UNLIMITED_, UNLIMITED_, real_t>& matrixAKA = params->getMatrixAKA();
      _Atom *atoms = complex->get(ATOM);
      complex->read(CHARGE, variables, atoms);
      typedef typename _Molecule::archetype_type _Archetype;
      std::vector<_Molecule*> molecules = complex->get(MOLECULE);
      unsigned numberOfMolecules = molecules.size();
//we should keep total molecule charge as const
      int atomOffset = 0;
      for (unsigned moleculeIndex = 0; moleculeIndex < numberOfMolecules; moleculeIndex++)
      {
        if (0 != moleculeIndex)
        {
          const _Archetype* prevArchetype = molecules[moleculeIndex - 1]->archetype();
          atomOffset += prevArchetype->count(ATOM);
        }
        real_t charge = 0.;
        const _Archetype* currentArchetype = molecules[moleculeIndex]->archetype();
        unsigned numberOfAtomsInMolecule = currentArchetype->count(ATOM);
        for (unsigned atomIndex = 0; atomIndex < numberOfAtomsInMolecule - 1; atomIndex++)
        {
          charge += q[atomIndex+atomOffset];
        }
        q[numberOfAtomsInMolecule + atomOffset - 1] = -charge;
      }
//calculating dh_dp and dh_dq
      for (unsigned i=0; i < params->getNumberOfAtoms(); i++)
      {
        real_t speedQ = 0;
        real_t speedP = RappleGoddardParams::instance()->find(make_string(atoms[i].atomdata->name)).electronegativity;
        for (unsigned j=0; j < params->getNumberOfAtoms(); j++)
        {
          speedQ += matrixAKA(i,j)*p[j];
          speedP += matrixJ(i,j)*q[j];
        }
        dp_dt[i] = -speedP;
        dq_dt[i] = speedQ;
      }
//normalization
      atomOffset = 0;
      for (unsigned moleculeIndex = 0; moleculeIndex < numberOfMolecules; moleculeIndex++)
      {
        if (0 != moleculeIndex)
        {
          const _Archetype* prevArchetype = molecules[moleculeIndex - 1]->archetype();
          atomOffset += prevArchetype->count(ATOM);
        }
        const _Archetype* currentArchetype = molecules[moleculeIndex]->archetype();
        unsigned numberOfAtomsInMolecule = currentArchetype->count(ATOM);
        real_t numerator = 0.;
        real_t denomirator = 0.;
        for (unsigned atomIndex = 0; atomIndex < numberOfAtomsInMolecule - 1; atomIndex++)
        {
          dp_dt[atomIndex + atomOffset] -= dp_dt[numberOfAtomsInMolecule + atomOffset - 1];
          dpCalculator(atomIndex, atomIndex, atomOffset, numberOfAtomsInMolecule, matrixAKA, dp_dt, &numerator, &denomirator);
        }
        dpCalculator(numberOfAtomsInMolecule - 1, numberOfAtomsInMolecule - 1, atomOffset, numberOfAtomsInMolecule, matrixAKA, dp_dt, &numerator, &denomirator);
        const _Bond* bonds = currentArchetype->get(BOND);
        int numberOfBonds = currentArchetype->count(BOND);
        for (unsigned bondIndex = 0; bondIndex < numberOfBonds; bondIndex++)
        {
          unsigned atomIndexA = bonds[bondIndex].ndx[0];
          unsigned atomIndexB = bonds[bondIndex].ndx[1];
          dpCalculator(atomIndexA, atomIndexB, atomOffset, numberOfAtomsInMolecule, matrixAKA, dp_dt, &numerator, &denomirator);
          dpCalculator(atomIndexB, atomIndexA, atomOffset, numberOfAtomsInMolecule, matrixAKA, dp_dt, &numerator, &denomirator);
        }
        dp_dt[numberOfAtomsInMolecule -1] = 0.;

      }
      for (unsigned i=0; i < params->getNumberOfAtoms(); i++)
      {
        q[i] += dq_dt[i]*time;
        p[i] += dp_dt[i]*time;
      }
    }
  };
}
#endif
