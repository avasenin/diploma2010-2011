!>
!! Assigns a Gaussian-type orbital to the atom
!! 
!! \note See research notes dated 2008-03-14
!<
      subroutine AssignsGTOBasis(theAtom)
        use Parameters
        implicit none
        type(Atom) :: theAtom
        double precision, external :: sGTOCoulInt
        integer :: i

!       Assign position
        theAtom%Basis%Position = theAtom%Position

!       Assign Gaussian orbital exponent
        do i=1,numParameterizedAtoms
           if (theAtom%Element%Z .eq. ParameterizedAtoms(i)%Z) then
              theAtom%Basis%zeta = GaussianExponent(i)
           end if
        end do

!       Scaling - does not work well with QEq parameters
!        theAtom%Basis%zeta = (theAtom%Element%Hardness/
!     &       sGTOCoulInt(1.0d0, 1.0d0, 0.0d0))**2
      end subroutine AssignsGTOBasis

!>
!! Calculates best-fit GTO exponent given best-fit STO exponent
!! \param n: principal quantum number
!! \param zeta: exponent for s-type Slater orbital
!! \return the best-fit exponent for the s-type Gaussian orbital
!! \note See research notes dated 2007-08-31
!! \deprecated
!<
      double precision  function sSTO2sGTO(n, zeta)
        implicit none
        integer, intent(in) :: n
        double precision, intent(in) :: zeta
        double precision, parameter :: conversion(1:7) = (/
     &       0.2709498089, 0.2527430925, 0.2097635701,
     &       0.1760307725, 0.1507985107, 0.1315902101,
     &       0.1165917484 /)
        sSTO2sGTO = conversion(n) * zeta * zeta
      end function sSTO2sGTO

!>
!! Calculates a best-fit Gaussian-type orbital (STO-1G) to
!! the Slater-type orbital defined from the hardness parameters
!! \param Hardness: chemical hardness in atomic units
!! \param n: principal quantum number
!! \note See research notes dated 2007-08-30
!! \deprecated
!<
      subroutine AssignsSTO1GBasis(theAtom)
        use Parameters
        implicit none
        type(Atom) :: theAtom
        double precision, external :: sSTOCoulInt
        double precision, external :: sSTO2sGTO
        integer :: n
        double precision :: zeta

C       Approximate the exact value of the constant of proportionality
C       by its value at a very small distance epsilon
C       since the exact R = 0 case has not be programmed into STOIntegrals
        double precision :: epsilon = 1.0d-6

C       Assign position
        theAtom%Basis%Position = theAtom%Position

C       Assign principal quantum number
        n = pqn(theAtom%Element)
        theAtom%Basis%n = n

C       Assign orbital exponent
        zeta = (sSTOCoulInt(1.0d0, 1.0d0, n, n, epsilon)
     &      /theAtom%Element%Hardness)**(-1.0d0/(3.0d0 + 2.0d0*n))

C       Rewrite it with best-fit Gaussian
        theAtom%Basis%zeta = sSTO2sGTO(n, zeta)
      end subroutine AssignsSTO1GBasis

!>
!! Computes Coulomb integral analytically over s-type GTOs
!!
!! Computes the two-center Coulomb integral over Gaussian-type orbitals
!! of s symmetry.
!!
!! \param a: Gaussian exponent of first atom in atomic units (inverse squared Bohr)
!! \param b: Gaussian exponent of second atom in atomic units (inverse squared Bohr)
!! \param R: internuclear distance in atomic units (bohr)
!! \return the value of the Coulomb potential energy integral
!! \note Reference: T. Helgaker, P. Jorgensen, J. Olsen, Molecular Electronic Structure Theory
!!                  Wiley, NY, 2000, Equations (9.7.21) and (9.8.23)
!<
      double precision  function sGTOCoulInt(a, b, R)
        implicit none
        double precision, intent(in) :: a,b,R
        intrinsic :: erf 
        double precision :: p
          
        p = sqrt(a * b / (a + b))
        sGTOCoulInt = erf(p * R) / R 
      end function sGTOCoulInt

!> Computes overlap integral analytically over s-type GTOs
!!
!! Computes the overlap integral over two Gaussian-type orbitals of s symmetry.
!! \param a: Gaussian exponent of first atom in atomic units (inverse squared Bohr)
!! \param b: Gaussian exponent of second atom in atomic units (inverse squared Bohr)
!! \param R: internuclear distance in atomic units (bohr)
!! \note Reference: T. Helgaker, P. Jorgensen, J. Olsen, Molecular Electronic Structure Theory
!!                  Wiley, NY, 2000, Equation (9.2.41)
!! \note With normalization constants added, calculates
!! \f[
!! S = \left(\frac{4\alpha\beta}{(\alpha + \beta)^2}\right)^\frac{3}{4}
!!     \exp\left(-\frac{\alpha\beta}{\alpha+\beta} R^2 \right)
!! \f]
!<
      double precision  function sGTOOvInt(a, b, R)
        implicit none
        double precision, intent(in) :: a,b,R 
        double precision :: p, q

        p = a + b
        q = a * b / p
        sGTOOvInt = (4*q/p)**0.75d0 * exp(-q*R*R)

!c       Sanity check
!        if (sGTOOvInt .ge. 1.0d0 .or. sGTOOvInt .lt. 0.0d0) then
!           print *, "Error: Overlap integral exceeds bounds: ",sGTOOvInt
!           print *, a, b, R
!           stop
!        end if
      end function sGTOOvInt

!>
!! Computes derivative of Coulomb integral wrt R
!! \param a: Gaussian exponent of first atom in atomic units (inverse squared Bohr)
!! \param b: Gaussian exponent of second atom in atomic units (inverse squared Bohr)
!! \param R: internuclear distance in atomic units (bohr)
!> \return the derivative of the Coulomb potential energy integral
!<
      double precision  function sGTOCoulIntGrad(a, b, R)
        implicit none
        double precision , intent(in) :: a,b,R
        double precision, parameter :: pi =  3.141592653589793d0
        double precision, external :: sGTOCoulInt
        double precision :: p
        
        if (abs(R) .eq. 0) then
           print *, "FATAL ERROR: R = 0 in sGTOCoulIntGrad"
           stop
        end if

        p = sqrt(a * b / (a + b))
        sGTOCoulIntGrad = 2.0d0 * p / (R * sqrt(pi)) * exp(-(p*R)**2)
     &       - sGTOCoulInt(a,b,R) / R
      end function sGTOCoulIntGrad

!>
!! Computes gradient of overlap integral wrt R
!!
!! Computes the derivative of the overlap integral over two Gaussian-type orbitals of s symmetry.
!! \param a: Gaussian exponent of first atom in atomic units (inverse squared Bohr)
!! \param b: Gaussian exponent of second atom in atomic units (inverse squared Bohr)
!! \param R: internuclear distance in atomic units (bohr)
!> \return the derivative of the sGTOOvInt integral
!<
      double precision function sGTOOvIntGrad(a,b,R)
        implicit none
        double precision, intent(in) :: a,b,R
        double precision, external :: sGTOOvInt

        sGTOOvIntGrad = -2 * (a*b)/(a+b)* R * sGTOOvInt(a,b,R)
      end function sGTOOvIntGrad
