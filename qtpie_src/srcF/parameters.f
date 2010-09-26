!>
!! Stores parameters for our charge models
!<
      module Parameters
        use AtomicUnits
        implicit none
        save

        double precision, parameter :: pi =  3.141592653589793d0

!>
!! Parameters for a s-type Slater type orbital (STO) basis function
!!
!! \param n : principal quantum number
!! \param zeta : zeta exponent with dimensions of inverse length in atomic units
!<
      type sSTO
        double precision, dimension(1:3) :: Position
        integer :: n
        double precision :: zeta
      end type

!> Parameters for a s-type Gaussian type orbital (GTO) basis function
!!
!! \param Position : an array of three double precisions describing Cartesian coordinates
!! \param zeta: exponent with dimensions of inverse square length in atomic units
!<
      type sGTO
        double precision, dimension(1:3) :: Position
        double precision :: zeta
      end type

!>
!! Atomic parameters
!! \param Symbol            : elemental symbol
!! \param Z                 : atomic number
!! \param FormalCharge      : formal charge, integers only
!! \param Electronegativity : Mulliken electronegativity in atomic units
!! \param Hardness          : Parr-Pearson chemical hardness in atomic units
!<
      type AtomData
        character (len = 2) :: Symbol
        integer             :: Z, FormalCharge
        double precision    :: Electronegativity, Hardness
      end type AtomData

!>
!! Describes an atom in a molecule
!!
!! \param  Element : type(AtomData) containing atomic parameters
!! \param Basis   : A basis function associated with the atom
!! \param Position: double precision(3) vector of Cartesian coordinates describing spatial location
!! \param Charge  : double precision, result of charge model calculation
!<
      type Atom
        type (AtomData) :: Element
        type (sSTO) :: Basis
        double precision, dimension(1:3) :: Position
        double precision :: Charge
      end type Atom
      
!>
!! Describes a molecular system
!!
!! \param Description: a text label of 132 characters
!! \param    NumAtoms: number of atoms (integer)
!! \param TotalCharge: total charge of system (double precision)
!! \param       Atoms: array of atoms
!! \param     Overlap: overlap matrix
!! \param      OvNorm: overlap norm vector (useful temporary variable)
!! \param     Coulomb: Coulomb matrix
!! \param      Energy: QTPIE contribution to the potential energy
!! \param   EGradient: Energy gradients
        type Molecule
          character (len=132) :: Description
          integer :: NumAtoms
          double precision :: TotalCharge
          Type(Atom), dimension(:), allocatable :: Atoms
          double precision, dimension(:,:), allocatable :: Overlap
          double precision, dimension(:), allocatable :: OvNorm
          double precision, dimension(:,:), allocatable :: Coulomb
          double precision :: Energy
          double precision, dimension(:,:), allocatable :: EGradient
        end type Molecule

C       Here are a bunch of predefined elements
C       As parameterized by Rappe and Goddard for QEq

        type(AtomData), parameter :: Hydrogen    = 
     &       AtomData( "H",  1, 0, 4.528*eV, 13.890*eV)
        type(AtomData), parameter :: Lithium     =
     &       AtomData("Li",  3, 0, 3.006*eV,  4.772*eV)
        type(AtomData), parameter :: Carbon      =
     &       AtomData( "C",  6, 0, 5.343*eV, 10.126*eV)
        type(AtomData), parameter :: Nitrogen   =
     &       AtomData( "N",  7, 0, 7.139*eV, 12.844*eV)
        type(AtomData), parameter :: Oxygen     =
     &       AtomData( "O",  8, 0, 8.741*eV, 13.364*eV)
        type(AtomData), parameter :: Fluorine   =
     &       AtomData( "F",  9, 0,10.874*eV, 14.948*eV)
        type(AtomData), parameter :: Sodium     = 
     &       AtomData("Na", 11, 0, 2.843*eV,  4.592*eV)
        type(AtomData), parameter :: Silicon    = 
     &       AtomData("Si", 14, 0, 4.168*eV,  6.974*eV)
        type(AtomData), parameter :: Phosphorus = 
     &       AtomData( "P", 15, 0, 5.463*eV,  8.000*eV)
        type(AtomData), parameter :: Sulphur    =
     &       AtomData( "S", 16, 0, 6.084*eV, 10.660*eV)
        type(AtomData), parameter :: Chlorine   =
     &       AtomData("Cl", 17, 0, 8.564*eV,  9.892*eV)
        type(AtomData), parameter :: Potassium  =
     &       AtomData( "K", 19, 0, 2.421*eV,  3.840*eV)
        type(AtomData), parameter :: Bromine    =
     &       AtomData("Br", 35, 0, 7.790*eV,  8.850*eV)
        type(AtomData), parameter :: Rubidium   =
     &       AtomData("Rb", 37, 0, 2.331*eV,  3.692*eV)
        type(AtomData), parameter :: Iodine =
     &       AtomData( "I", 53, 0, 6.822*eV,  7.524*eV)
        type(AtomData), parameter :: Cesium =
     &       AtomData("Cs", 55, 0, 2.183*eV,  3.422*eV)

        integer, parameter :: numParameterizedAtoms = 16 !< Number of defined atomic parameters
        type(AtomData), parameter :: 
     &   ParameterizedAtoms(numParameterizedAtoms) = 
     &   (/Hydrogen, Lithium, Carbon, Nitrogen, Oxygen, Fluorine,
     &     Sodium, Silicon, Phosphorus, Sulphur, Chlorine,
     &     Potassium, Bromine, Rubidium, Iodine, Cesium /) !< Array of defined atomic parameters

C       Parameters for cations. All experimental values!
!>
!!      Sodium cation
!<
        type(AtomData), parameter :: SodiumCation = 
     &       AtomData("Na",11,+1,4562*kJ_mol, 5.13908*eV)



c       Data for newly parameterized Gaussian basis set
        double precision, parameter, dimension(numParameterizedAtoms) ::
     &       GaussianExponent =
     &   (/  0.534337523756312, 0.166838519142176, 0.206883838259186,
     &       0.221439796025873, 0.223967308625516, 0.231257590182828,
     &       0.095892938712585, 0.105219608142377, 0.108476721661715,
     &       0.115618357843499, 0.113714050615107, 0.060223294377778,
     &       0.070087547802259, 0.041999054745368, 0.068562697575073,
     &       0.030719481189777 /)

!>
!overlape!      Threshold for calculating overlap integrals
!<
        double precision, parameter :: OvIntThreshold = 1.0d-9

        double precision :: SmallestGaussianExponentInSystem = 1.0d40
!>
!!      Store pre-calculated thresholds for prescreening
!<
        double precision :: OvIntMaxR
!>
!!      Threshold for calculating Coulomb integrals
!<
        double precision, parameter :: CoulIntThreshold = 1.0d-9
!>
!!      Store pre-calculated thresholds for prescreening
!<
        double precision :: CoulIntMaxR
        contains

!>
!! Computes the expectation value of the radial distance over s-type STOs
!! \param Basis: s-type STO basis function
!! \return the expectation value of the radial distance over s-type STOs
!<
          double precision function ExpectR(Basis)
            implicit none
            type(sSTO), intent(in) :: Basis
            ExpectR = (Basis%n + 0.5) / Basis%zeta
          end function ExpectR

!>
!! Computes the principal quantum number of an atom given its atomic number
!! \param theAtom Atom to determine principle quantum number for
!! \return the principal quantum number
!<
      integer function pqn(theAtom)
        implicit none
        integer :: j
        integer, parameter :: maxelectrons(7)=
     &    (/ 2, 10, 18, 36, 54, 86, 118 /) !< Lookup table for max number of electrons for that quantum number
        type(AtomData), intent(in) :: theAtom
        pqn=1
C       work through each shell
        do j=1,7
          if (theAtom%Z.gt.maxelectrons(j)) then
             pqn=pqn+1
          end if
        end do
      end function pqn

      end module Parameters
