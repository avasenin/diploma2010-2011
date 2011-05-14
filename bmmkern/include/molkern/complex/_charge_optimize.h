#ifndef _CHARGE_OPTIMIZE__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H
#define _CHARGE_OPTIMIZE__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H

#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_io.hpp"

#include "molkern/__moldefs.h"
#include "molkern/complex/_coulomb_params.h"
#include "molkern/forcefield/_charge_optimize_params.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
namespace molkern
{
  using namespace prgkern;
  using namespace std;
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

      CoulombParams::build(basis);
      ChargeOptimizeParams<LPComplex>* params = new ChargeOptimizeParams<LPComplex>(complex);

      //collapseEnumeration(charges, m_variables, params);

      energy = Minimizer_<OPTIMIZER_TYPE>::operator()(
        &dU__dQ<LPComplex>, params,
        2 * number_of_atoms - 1, m_variables, m_derivatives, param.maxiter,
        param.stpmin * (real_t) 0.001, param.stpmax,
        param.maxfev, param.maxhalt,
        param.wolfe1 * (real_t) 0.001, param.wolfe2 * (real_t) 0.001,
        param.xtol * (real_t) 0.001, param.ftol * (real_t) 0.001, param.gtol * (real_t) 0.001,
        param.m * 1000, param.steep * 1000
      );
      ofstream fqeq("./qeq_.txt");
      for (unsigned i=0; i < number_of_atoms; i++)
      {
    	  fqeq << setiosflags(ios::fixed) << setprecision(5) << atoms[i].charge << endl;
      }
      fqeq.close();
      memset(m_variables, 0, 2 * number_of_atoms * sizeof(real_t));
      for (unsigned i=0; i < number_of_atoms; i++)
      {
    	m_variables[i] = atoms[i].charge;
      }
      //m_variables[0] = atoms[0].charge;
      //m_variables[1] = 0.41886;
      //m_variables[2] = 0.47857;
      ofstream own("./own_.txt");
      own << setiosflags(ios::fixed) << setprecision(5) << "QEq "<< atoms[2].charge << endl;
      real_t* avarage = new real_t[number_of_atoms];
      memset(avarage, 0, number_of_atoms * sizeof(real_t));
      real_t last = 0.;
      int period = -1;
      for (int k=0; k < 2000; k++)
      {
        for (int i=0; i < 1000; i++)
        {
          modeling(m_variables, 0.0001, params);
        }
        avarage[2] = (avarage[2] * k + m_variables[2]) / (k + 1);
        own << setiosflags(ios::fixed) << setprecision(5) <<  m_variables[2] << " " << avarage[2] << endl;
        if (-1 != period && 0 != k && (0 <= max(last, m_variables[5])) &&  (min(last, m_variables[5]) < 0))
        {
        	period = k;
        }
        last = m_variables[5];
      }
      own.close();
      delete params;
      delete avarage;
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
      memset(derivatives, 0, 2 * params->getNumberOfAtoms()*sizeof(real_t));

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
//calculating dq_dt
      atomOffset = 0;
      int bondOffset = 0;
      for (unsigned moleculeIndex = 0; moleculeIndex < numberOfMolecules; moleculeIndex++)
      {
        if (0 != moleculeIndex)
        {
          const _Archetype* prevArchetype = molecules[moleculeIndex - 1]->archetype();
          atomOffset += prevArchetype->count(ATOM);
          bondOffset += prevArchetype->count(BOND);
        }
        const _Archetype* currentArchetype = molecules[moleculeIndex]->archetype();
        for (unsigned bondIndex = 0; bondIndex < currentArchetype->count(BOND); bondIndex++)
        {
          unsigned i = currentArchetype->get(BOND)[bondIndex].ndx[0];
          unsigned j = currentArchetype->get(BOND)[bondIndex].ndx[1];
          dq_dt[i] += matrixAKA(i,j)*p[j];
          dq_dt[j] += matrixAKA(j,i)*p[i];
        }
      }
      for (unsigned i = 0; i < params->getNumberOfAtoms(); i++)
      {
        dq_dt[i] += matrixAKA(i,i)*p[i];
      }
//calculation dp_dt
      const vector< pair<int, int> >& notZeroMatrixJIndexes = params->getNotZeroMatrixJIndexes();
      for (unsigned i=0; i < params->getNumberOfAtoms(); i++)
      {
        dp_dt[i] = -RappleGoddardParams::instance()->find(make_string(atoms[i].atomdata->name)).electronegativity;
      }
      for (unsigned k=0; k < notZeroMatrixJIndexes.size(); k++)
      {
        int i = notZeroMatrixJIndexes[k].first;
        int j = notZeroMatrixJIndexes[k].second;
        dp_dt[i] -= matrixJ(i,j)*q[j];
        if (i != j)
        {
          dp_dt[j] -= matrixJ(j,i)*q[i];
        }
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
        for (unsigned atomIndex = 0; atomIndex < numberOfAtomsInMolecule - 1; atomIndex++)
        {
          dp_dt[atomIndex + atomOffset] -= dp_dt[numberOfAtomsInMolecule + atomOffset - 1];
        }
        const _Bond* bonds = currentArchetype->get(BOND);
        int numberOfBonds = currentArchetype->count(BOND);
        for (unsigned bondIndex = 0; bondIndex < numberOfBonds; bondIndex++)
        {
          unsigned atomIndexA = bonds[bondIndex].ndx[0];
          unsigned atomIndexB = bonds[bondIndex].ndx[1];
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
