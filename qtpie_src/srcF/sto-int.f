!>
!! Calculates a Slater-type orbital exponent
!! based on the hardness parameters
!! \param Hardness: chemical hardness in atomic units
!! \param        n: principal quantum number
!! \note See research notes dated 2007-08-30
!<
      subroutine AssignsSTOBasis(theAtom)
        use Parameters
        implicit none
        type(Atom) :: theAtom
        double precision, external :: sSTOCoulInt
        integer :: n

C       Approximate the exact value of the constant of proportionality
C       by its value at a very small distance epsilon
C       since the exact R = 0 case has not be programmed
        double precision :: epsilon = 1.0d-8

C       Assign position
        theAtom%Basis%Position = theAtom%Position

C       Assign principal quantum number
        n = pqn(theAtom%Element)
        theAtom%Basis%n = n

C       Assign orbital exponent
        theAtom%Basis%zeta = (sSTOCoulInt(1.0d0, 1.0d0, n, n, epsilon)
     &      /theAtom%Element%Hardness)**(-1.0d0/(3.0d0 + 2.0d0*n))
      end subroutine AssignsSTOBasis

!>
!! Computes Rosen's Guillimer-Zener function A
!!
!! Computes Rosen's A integral, an auxiliary quantity needed to
!! compute integrals involving Slater-type orbitals of s symmetry.
!! \f[
!! A_n(\alpha) = \int_1^\infty x^n e^{-\alpha x}dx
!! = \frac{n! e^{-\alpha}}{\alpha^{n+1}}\sum_{\nu=0}^n
!! \frac{\alpha^\nu}{\nu!}
!! \f]
!! \param n - principal quantum number
!! \param alpha - Slater exponent 
!! \return the value of Rosen's A integral
!! \note N. Rosen, Phys. Rev., 38 (1931), 255
!<
      double precision function RosenA(n,a)
        implicit none
        integer, intent(in) :: n
        double precision, intent(in) :: a
        double precision :: Term
        integer :: nu
        if (a.ne.0.0d0) then
          Term = 1.0d0
          RosenA = Term
          do nu = 1,n
            Term = a/nu*Term
            RosenA = RosenA + Term
          end do
          RosenA=RosenA/Term*exp(-a)/a
        else
          RosenA=0.0d0
        end if
      end function RosenA

!>
!! Computes Rosen's Guillimer-Zener function B
!!
!! Computes Rosen's B integral, an auxiliary quantity needed to
!! compute integrals involving Slater-type orbitals of s symmetry.
!! \f[
!! B_n(\alpha) = \int_{-1}^1 x^n e^{-\alpha x} dx
!!             = \frac{n!}{\alpha^{n+1}} 
!!               \sum_{\nu=0}^n \frac{e^\alpha(-\alpha)^\nu
!!                 - e^{-\alpha} \alpha^\nu}{\nu!}
!! \f]
!! \param n - principal quantum number
!! \param alpha - Slater exponent 
!! \return the value of Rosen's B integral
!! \note N. Rosen, Phys. Rev., 38 (1931), 255
!<
      double precision  function RosenB(n,alpha)
        implicit none
        integer, intent(in) :: n
        double precision , intent(in) :: alpha
        double precision :: TheSum, Term
        double precision :: PSinhRosenA, PCoshRosenA, PHyperRosenA
        integer :: nu
        logical :: IsPositive
        if (alpha.ne.0.0d0) then
          Term = 1.0d0
          TheSum = 1.0d0
          IsPositive = .True.
C         These two expressions are (up to constant factors) equivalent
C         to computing the hyperbolic sine and cosine of a respectively
C         The series consists of adding up these terms in an
C         alternating fashion
          PSinhRosenA =  exp(alpha) - exp(-alpha)
          PCoshRosenA = -exp(alpha) - exp(-alpha)
          TheSum=PSinhRosenA
          do nu = 1,n
            if (isPositive) then
              PHyperRosenA = PCoshRosenA
              isPositive = .False.
            else !term to add should be negative
              PHyperRosenA = PSinhRosenA
              isPositive = .True.
            end if
            Term=alpha/(1.0d0*nu)*Term
            TheSum=TheSum+Term*PHyperRosenA
          end do
          RosenB=TheSum/(alpha*Term)
        else
C         pathological case of a=0
          print *, "WARNING, a = 0 in RosenB"
          RosenB=(1.0d0-(-1.0d0)**n)/(n+1.0d0)
        end if
      end function RosenB

!>
!! Computes Rosen's D combinatorial factor  
!!
!! Computes Rosen's D factor, an auxiliary quantity needed to
!! compute integrals involving Slater-type orbitals of s symmetry.
!! \f[
!! RosenD^{mn}_p = \sum_k (-1)^k \frac{m! n!}
!!                 {(p-k)!(m-p+k)!(n-k)!k!}
!! \f]
!! \return the value of Rosen's D factor
!! \note N. Rosen, Phys. Rev., 38 (1931), 255
!<
      integer function RosenD(m,n,p)
        use Factorial
        implicit none
        integer, intent(in) :: m,n,p
        integer k
        RosenD = 0

        if (m+n+p.gt.maxFact) then
          print *, "Error, arguments exceed maximum factorial computed"
     &      , m+n+p,">",maxFact
          stop
        end if
        do k=max(p-m,0),min(n,p)
           if (mod(k,2).eq.0) then
           RosenD = RosenD + fact(m) / ( fact(p-k) *
     &       fact(m-p+k)) * fact(n) / (fact(n-k) * fact(k))
           else
           RosenD = RosenD - fact(m) / ( fact(p-k) *
     &       fact(m-p+k)) * fact(n) / (fact(n-k) * fact(k))
           end if
        end do
      end function RosenD

!>
!! Computes Coulomb integral analytically over s-type STOs
!!
!! Computes the two-center Coulomb integral over Slater-type
!! orbitals of s symmetry.
!! \param a: Slater zeta exponent of first atom in inverse Bohr (au)
!! \param b: Slater zeta exponent of second atom in inverse Bohr (au)
!! \param m: principal quantum number of first atom
!! \param n: principal quantum number of second atom
!! \param R: internuclear distance in atomic units (bohr)
!! \return value of the Coulomb potential energy integral
!! \note N. Rosen, Phys. Rev., 38 (1931), 255
!! \note In Rosen's paper, this integral is known as K2.
!<
      double precision  function sSTOCoulInt(a, b, m, n, R)
        use Factorial
        implicit none
        integer, intent(in) :: m,n
        double precision , intent(in) :: a,b,R
        double precision , external :: RosenA, RosenB
        integer, external :: RosenD
        integer :: nu, p
        double precision :: x, K2
        double precision :: Factor1, Factor2, Term, OneElectronTerm
        double precision :: eps, epsi

C       To speed up calculation, we terminate loop once contributions
C       to integral fall below the bound, epsilon
        double precision, parameter :: epsilon = 0.0d0

C       x is the argument of the auxiliary RosenA and RosenB functions
        x=2.0*a*R
C       First compute the two-electron component
        sSTOCoulInt = 0.0d0
        if (x.eq.0) then
C         Pathological case
          if ((a.eq.b).and.(m.eq.n)) then
            do nu = 0,2*n-1
              K2 = 0.0d0
              do p = 0, 2*n+m
                K2 = K2 + 1.0d0 / fact(p)
              end do
              sSTOCoulInt = sSTOCoulInt + fact(2*n+m)/fact(m)*K2
            end do  
            sSTOCoulInt = 2*a/(n*fact(2*n))*sSTOCoulInt
          else
C          Not implemented
            print *, "ERROR, sSTOCoulInt cannot compute from arguments"
            print *, "a = ",a,"b = ",b,"m =",m,"n = ",n,"R = ",R
            stop
          end if
        else
          OneElectronTerm = 1.0d0/R + x**(2*m)/(fact(2*m)*R)*
     &              ((x-2*m)*RosenA(2*m-1,x)-exp(-x)) + sSTOCoulInt
          eps = epsilon / OneElectronTerm
          if (a.eq.b) then
C           Apply Rosen (48)
            Factor1 = -a*(a*R)**(2*m)/(n*fact(2*m))
            do nu=0,2*n-1
              Factor2 = (2.0d0*n-nu)/fact(nu)*(a*R)**nu
              epsi = eps / abs(Factor1 * Factor2)
              K2=0.0d0
              do p=0,m+(nu-1)/2
                Term = RosenD(2*m-1,nu,2*p)/(2.0d0*p+1.0d0)
     &               *RosenA(2*m+nu-1-2*p,x)
                K2=K2 + Term
                if ((Term.gt.0).and.(Term.lt.epsi)) then
                  goto 1
                end if
              end do
              sSTOCoulInt=sSTOCoulInt+K2*Factor2
            end do
 1          sSTOCoulInt=sSTOCoulInt*Factor1
          else
            Factor1 = -a*(a*R)**(2*m)/(2.0d0*n*fact(2*m))
            epsi = eps/abs(Factor1)
            if (b.eq.0.0d0) then
              print *, "WARNING: b = 0 in sSTOCoulInt"
            else  
C           Apply Rosen (54)
            do nu=0,2*n-1
              K2=0
              do p=0,2*m+nu-1
                K2=K2+RosenD(2*m-1,nu,p)*RosenB(p,R*(a-b))
     &                 *RosenA(2*m+nu-1-p,R*(a+b))
              end do
              Term = K2*(2*n-nu)/fact(nu)*(b*R)**nu
              sSTOCoulInt=sSTOCoulInt+Term
              if (abs(Term) .lt. epsi) then
                goto 2
              end if
            end do
 2          sSTOCoulInt=sSTOCoulInt*Factor1
            end if
          end if  
C         Now add the one-electron term from Rosen (47) = Rosen (53)
          sSTOCoulInt=sSTOCoulInt + OneElectronTerm
        end if
      end function sSTOCoulInt

!>
!! Computes overlap integral analytically over s-type STOs
!!
!!       Computes the overlap integral over two
!!      Slater-type orbitals of s symmetry.
!! \param a: Slater zeta exponent of first atom in inverse Bohr (au)
!! \param b: Slater zeta exponent of second atom in inverse Bohr (au)
!! \param m: principal quantum number of first atom
!! \param n: principal quantum number of second atom
!! \param R: internuclear distance in atomic units (bohr)
!! \return the value of the sSTOOvInt integral
!! \note N. Rosen, Phys. Rev., 38 (1931), 255
!! \note In the Rosen paper, this integral is known as I.
!<
      double precision  function sSTOOvInt(a,b,m,n,R)
        use Factorial
        implicit none
        integer, intent(in) :: m,n
        double precision , intent(in) :: a,b,R
        double precision , external :: RosenA, RosenB
        integer, external :: RosenD
        integer :: q

        double precision :: Factor, Term, eps

C       To speed up calculation, we terminate loop once contributions
C       to integral fall below the bound, epsilon
        double precision, parameter :: epsilon = 0.0d0

        sSTOOvInt=0.0d0

        if (a.eq.b) then
          Factor = (a*R)**(m+n+1)/sqrt(fact(2*m)*fact(2*n))
          eps = epsilon / abs(Factor)
          do q=0,(m+n)/2
            Term = RosenD(m,n,2*q)/(2.0d0*q+1.0d0)*RosenA(m+n-2*q,a*R)
            sSTOOvInt=sSTOOvInt+Term
            if (abs(Term).lt.eps) then
              exit
            end if
          end do
          sSTOOvInt=sSTOOvInt*Factor
        else
          Factor = 0.5d0*(a*R)**(m+0.5d0)*(b*R)**(n+0.5d0)
     &         /sqrt(fact(2*m)*fact(2*n))
          eps = epsilon / abs(Factor)
          do q=0,m+n
            Term = RosenD(m,n,q)*RosenB(q,R/2.0d0*(a-b))
     &            *RosenA(m+n-q,R/2.0d0*(a+b))
            sSTOOvInt=sSTOOvInt+Term
            if (abs(Term) .lt. eps) then
              exit
            end if
          end do
          sSTOOvInt=sSTOOvInt*Factor
        end if
      end function sSTOOvInt

!>
!! Computes kinetic energy integral analytically over s-type STOs
!!
!! Computes the overlap integral over two Slater-type orbitals of s symmetry.
!! \param a: Slater zeta exponent of first atom in inverse Bohr (au)
!! \param b: Slater zeta exponent of second atom in inverse Bohr (au)
!! \param m: principal quantum number of first atom
!! \param n: principal quantum number of second atom
!! \param R: internuclear distance in atomic units (bohr)
!! \return the value of the kinetic energy integral
!! \note N. Rosen, Phys. Rev., 38 (1931), 255
!! \note untested
!<
      double precision  function KinInt(a,b,m,n,R)
        implicit none
        integer, intent(in) :: m,n
        double precision , intent(in) :: a,b,R
        double precision , external :: sSTOOvInt
        KinInt=-0.5*b*b*sSTOOvInt(a,b,m,n,R)
        if (n.gt.0) then
          KinInt=KinInt+b*b*(2*b/(2*b-1))**0.5*sSTOOvInt(a,b,m,n-1,R)
          if (n.gt.1) then
            KinInt=KinInt+(n*(n-1)/((n-0.5)*(n-1.5)))**0.5
     &                       *sSTOOvInt(a,b,m,n-2,R)
          end if
        end if
      end function

!>
!! Computes derivative of Coulomb integral with respect to the interatomic distance
!!
!! Computes the two-center Coulomb integral over Slater-type orbitals of s symmetry.
!! \param a: Slater zeta exponent of first atom in inverse Bohr (au)
!! \param b: Slater zeta exponent of second atom in inverse Bohr (au)
!! \param m: principal quantum number of first atom
!! \param n: principal quantum number of second atom
!! \param R: internuclear distance in atomic units (bohr)
!! \return the derivative of the Coulomb potential energy integral
!! \note Derived in QTPIE research notes, May 15 2007
!<
      double precision  function sSTOCoulIntGrad(a, b, m, n, R)
        use Factorial
        implicit none
        integer, intent(in) :: m,n
        double precision , intent(in) :: a,b,R
        double precision , external :: RosenA, RosenB, sSTOCoulInt
        integer, external :: RosenD
C       loop counters
        integer :: nu, p
C       temporary quantities
        double precision :: x, y, z, K2, TheSum
C       x is the argument of the auxiliary RosenA and RosenB functions
        x=2.0*a*R
C       First compute the two-electron component
        sSTOCoulIntGrad = 0.0d0
        if (x.eq.0) then
C         Pathological case
          print *, "WARNING: argument given to sSTOCoulIntGrad is 0"
          print *, "a = ", a, "R = ", R
        else
          if (a.eq.b) then
            TheSum = 0.0d0
            do nu=0,2*(n-1)
              K2 = 0.0d0
              do p=0,(m+nu)/2
                K2 = K2 + RosenD(2*m-1, nu+1, 2*p)/(2*p + 1.0d0)
     &                  * RosenA(2*m+nu-1-2*p, x)
              end do
              TheSum = TheSum + (2*n-nu-1)/fact(nu)*(a*R)**(nu) * K2
            end do
            sSTOCoulIntGrad = -a**(2*m+2)*R**(2*m)
     &                      /(n*fact(2*m))*TheSum
            TheSum = 0.0d0
            do nu=0,2*n-1
              K2 = 0.0d0
              do p=0,(m+nu-1)/2
                K2 = K2 + RosenD(2*m-1, nu, 2*p)/(2*p + 1.0d0)
     &                  * RosenA(2*m+nu-2*p, x)
              end do
              TheSum = TheSum + (2*n-nu)/fact(nu)*(a*R)**nu * K2
            end do
            sSTOCoulIntGrad = sSTOCoulIntGrad + 2*a**(2*m+2)*R**(2*m)
     &                                  /(n*fact(2*m))*TheSum
          else
C           Slater exponents are different
C           First calculate some useful arguments
            y = R*(a+b)
            z = R*(a-b)
            TheSum = 0.0d0
            do nu=0,2*n-1              
              K2 = 0.0d0
              do p=0,2*m+nu
                K2 = K2 + RosenD(2*m-1, nu+1, p)
     &                  * RosenB(p,z)*RosenA(2*m+nu-p, y)
              end do
              TheSum = TheSum + (2*n-nu-1)/fact(nu)*(b*R)**nu * K2
            end do
            sSTOCoulIntGrad = -b*a**(2*m+1)*R**(2*m)/
     &                       (2*n*fact(2*m))*TheSum
            TheSum = 0.0d0
            do nu=0,2*n
              K2 = 0.0d0
              do p=0,2*m-1+nu
                K2 = K2 + RosenD(2*m-1, nu, p)
     &             * ((a-b)*RosenB(p+1,z)*RosenA(2*m+nu-p-1, y)
     &               +(a+b)*RosenB(p  ,z)*RosenA(2*m+nu-p  , y))
              end do
              TheSum = TheSum + (2*n-nu)/fact(nu)*(b*R)**nu * K2
            end do
            sSTOCoulIntGrad = sSTOCoulIntGrad
     &           + a**(2*m+1)*R**(2*m)/(2*n*fact(2*m))*TheSum
          end if
C         Now add one-electron terms and common term
          sSTOCoulIntGrad = sSTOCoulIntGrad - (2.0d0*m+1.0d0)/R**2
     &                  + 2.0d0*m/R * sSTOCoulInt(a,b,m,n,R)
     &         +x**(2*m)/(fact(2*m)*R**2) * ((2.0d0*m+1.0d0)*exp(-x)
     &          +2.0d0*m*(1.0d0+2.0d0*m-x)*RosenA(2*m-1,x))
        end if
      end function sSTOCoulIntGrad

!> Computes gradient of overlap integral with respect to the interatomic diatance
!!
!! Computes the derivative of the overlap integral over two Slater-type orbitals of s symmetry.
!! \param a: Slater zeta exponent of first atom in inverse Bohr (au)
!! \param b: Slater zeta exponent of second atom in inverse Bohr (au)
!! \param m: principal quantum number of first atom
!! \param n: principal quantum number of second atom
!! \param R: internuclear distance in atomic units (bohr)
!! \return the derivative of the sSTOOvInt integral
!! \note Derived in QTPIE research notes, May 15 2007
!<
      double precision  function sSTOOvIntGrad(a,b,m,n,R)
        use Factorial
        implicit none
        integer, intent(in) :: m,n
        double precision, intent(in) :: a,b,R
        double precision, external :: RosenA, RosenB
        integer, external :: RosenD
        double precision, external :: sSTOOvInt
C       Useful temporary quantities
        double precision :: w, x, y, z, TheSum
C       Loop variable
        integer :: q

C       Calculate first term
        sSTOOvIntGrad=(m+n+1.0d0)/R * sSTOOvInt(a,b,m,n,R)
C       Calculate remaining terms; answers depend on exponents 
        TheSum = 0.0d0
        x = a * R
        if (a.eq.b) then
          do q = 0,(m+n)/2
            TheSum = TheSum + RosenD(m,n,2*q) / (2*q + 1.0d0)
     &                * RosenA(m+n-2*q+1, x)
          end do
          sSTOOvIntGrad = sSTOOvIntGrad - a*x**(m+n+1)/
     &                      sqrt(fact(2*m)*fact(2*n))*TheSum
        else
C         Useful arguments
          w = b*R
          y = 0.5d0*R*(a+b)
          z = 0.5d0*R*(a-b)
          do q = 0,m+n
            TheSum = TheSum + RosenD(m,n,q) *
     &            ((a-b)*RosenB(q+1,z)*RosenA(m+n-q  ,y)
     &            +(a+b)*RosenB(q  ,z)*RosenA(m+n-q+1,y))
          end do
          sSTOOvIntGrad = sSTOOvIntGrad - 0.25d0*sqrt((x**(2*m+1)
     &      *w**(2*n+1))/(fact(2*m)*fact(2*n)))*TheSum

        end if
      end function sSTOOvIntGrad
