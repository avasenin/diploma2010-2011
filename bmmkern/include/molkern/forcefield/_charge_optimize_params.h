#ifndef _CHARGE_OPTIMIZE_PARAMS__0077A726_CBE6_5c54_73A2_F9434C9C0000__H
#define _CHARGE_OPTIMIZE_PARAMS__0077A726_CBE6_5c54_73A2_F9434C9C0000__H

#pragma once

using namespace std;
namespace molkern
{
  template<typename TComplex>
  class ChargeOptimizeParams
  {
private:
    mdense_<UNLIMITED_, UNLIMITED_, real_t> m_A;
    mdense_<UNLIMITED_, UNLIMITED_, real_t> m_AKA;
    mdense_<UNLIMITED_, UNLIMITED_, real_t> m_K;
    mdense_<UNLIMITED_, UNLIMITED_, real_t> m_J;


    unsigned m_numberOfBonds;
    unsigned m_numberOfAtoms;
    TComplex* m_complex;
public:

    unsigned getNumberOfAtoms() { return m_numberOfAtoms; }
    unsigned getNumberOfBonds() { return m_numberOfBonds; }

    const mdense_<UNLIMITED_, UNLIMITED_, real_t>& getMatrixJ() { return m_J; };
    const mdense_<UNLIMITED_, UNLIMITED_, real_t>& getMatrixAKA() { return m_AKA; };

    template<typename T>
    void printMatrix(mdense_<UNLIMITED_, UNLIMITED_, T> &U, const char* header)
     {
       PRINT('\n');
       PRINT(header);
       PRINT('\n');
       for (unsigned i=0; i < U.dimension()[0]; i++)
       {
         for (unsigned j=0; j < U.dimension()[1]; j++)
         {
           PRINT(U(i,j));
           PRINT(" ");
         }
         PRINT('\n');
       }
     };
    TComplex* getComplex() const
    {
      return m_complex;
    }
    ChargeOptimizeParams(TComplex* complex)
    {
      typedef typename TComplex::molecule_type _Molecule;
      m_numberOfAtoms = complex->count(ATOM);
      m_numberOfBonds = 0;
      std::vector<_Molecule*> molecules = complex->get(MOLECULE);
      unsigned numberOfMolecules = molecules.size();
      for (unsigned moleculeIndex = 0; moleculeIndex < numberOfMolecules; moleculeIndex++)
      {
         m_numberOfBonds += molecules[moleculeIndex]->archetype()->count(BOND);
      }
      m_complex = complex;
      buildMatrixA();
      buildMatrixK();
      buildMatrixAKA();
      buildMatrixJ();
    };

    ~ChargeOptimizeParams()
    {
    };
    //A*K^-1*At
    void buildMatrixAKA()
    {
      typedef typename TComplex::molecule_type _Molecule;
      typedef typename _Molecule::archetype_type _Archetype;
      std::vector<_Molecule*> molecules = m_complex->get(MOLECULE);
      unsigned numberOfMolecules = molecules.size();
      m_AKA.resize(m_numberOfAtoms, m_numberOfAtoms);
      int atomOffset = 0;
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
          unsigned k = bondIndex + bondOffset;
          unsigned i = currentArchetype->get(BOND)[bondIndex].ndx[0];
          unsigned j = currentArchetype->get(BOND)[bondIndex].ndx[1];
          m_AKA(i,j) = -m_K(k,k);
          m_AKA(j,i) = -m_K(k,k);
          m_AKA(i,i) += m_K(k,k);
          m_AKA(j,j) += m_K(k,k);
        }
      }
  //    printMatrix(m_AKA, "Matrix AKA");
    };

    void buildMatrixK()
    {
      typedef typename TComplex::atom_type _Atom;
      m_K.resize(m_numberOfBonds, m_numberOfBonds);
      _Atom* atoms = m_complex->get(ATOM);
      for (unsigned bondIndex = 0; bondIndex < m_numberOfBonds; bondIndex++)
      {
        unsigned firstAtomIndexInCurrentBond = -1;
        unsigned secondAtomIndexInCurrentBond = -1;
        for (unsigned atomIndex = 0; atomIndex < m_numberOfAtoms; atomIndex++)
        {
          switch ((int)m_A(atomIndex, bondIndex))
          {
            case -1:
              firstAtomIndexInCurrentBond = atomIndex;
              break;
            case 1:
              secondAtomIndexInCurrentBond = atomIndex;
              break;
          }
        }
        m_K(bondIndex, bondIndex) = CoulombParams::instance()->get(atoms[firstAtomIndexInCurrentBond], atoms[secondAtomIndexInCurrentBond]);
      }
      //printMatrix(m_K, "Matrix K");
    }

    vector< pair<int,int> > m_effectiveJValues;
    const vector< pair<int,int> >& getNotZeroMatrixJIndexes() const {return m_effectiveJValues;}
#define EPS 0.0001
    void buildMatrixJ()
    {
      typedef typename TComplex::atom_type _Atom;
      _Atom* atoms = m_complex->get(ATOM);
      m_J.resize(m_numberOfAtoms, m_numberOfAtoms);
      for (unsigned i=0; i < m_numberOfAtoms; i++)
      {
        for (unsigned j=i; j < m_numberOfAtoms; j++)
        {
          double value =  CoulombParams::instance()->get(atoms[i], atoms[j]);
          m_J(i, j) = value;
          m_J(j, i) = value;
          if (EPS < abs(value)) {
            m_effectiveJValues.push_back(pair<int, int>(i,j));
          }
        }
      }
      printMatrix(m_J, "Matrix J");
    }

    void buildMatrixA()
    {
      typedef typename TComplex::molecule_type _Molecule;
      typedef typename _Molecule::archetype_type _Archetype;
      std::vector<_Molecule*> molecules = m_complex->get(MOLECULE);
      unsigned numberOfMolecules = molecules.size();
      m_A.resize(m_numberOfAtoms, m_numberOfBonds);

      int atomOffset = 0;
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
          m_A(currentArchetype->get(BOND)[bondIndex].ndx[0] + atomOffset, bondIndex+bondOffset) = 1;
          m_A(currentArchetype->get(BOND)[bondIndex].ndx[1] + atomOffset, bondIndex+bondOffset) = -1;
        }

      }
    //  printMatrix(m_A, "Matrix A");
    }
  };

}


#endif
