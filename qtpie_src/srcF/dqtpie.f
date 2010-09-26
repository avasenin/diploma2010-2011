!>
!! Populates atomic charges according to the QTPIE charge model
!! \param Mol : of the Molecule data type
!! Mol%Atoms(i)%Charge are computed
!! \note The model is described in the paper below:
!!       J. Chen and T. J. Martinez, Chem. Phys. Lett., 438 (4-6), 2007, 315-320
!!       doi:10.1016/j.cplett.2007.02.065
!<
      subroutine dQTPIE(Mol)
        use Parameters
        implicit none
        type(Molecule) :: Mol

        double precision, dimension(:), allocatable :: Voltage, Charge

        double precision :: ThisCharge !< Temporary atomic charge variable
        double precision :: Overlap !< Temporary overlap integral
        double precision :: R !< Temporary distance vairiable

C       i1-i2 loop over atoms
        integer :: stat, i1, i2

        integer :: N !< size of matrix problem

        double precision, external :: Distance, sGTOOvInt
        external :: dsolver

C       Define size of matrix problem
        N = Mol%NumAtoms

C       Allocate scratch memory
        allocate(Voltage(N+1), Charge(N+1), STAT=stat)
        if (stat .ne. 0) then
           print *,"QTPIE: ERROR: could not allocate sufficient memory"
           stop
        end if

*       Construct problem

C       Calculate integral pre-screening threshold
        do i1 = 1,Mol%NumAtoms
          SmallestGaussianExponentInSystem = min(
     &       SmallestGaussianExponentInSystem,
     &       Mol%Atoms(i1)%Basis%zeta)
        end do

        OvIntMaxR = sqrt(
     &      log( (pi/(2*SmallestGaussianExponentInSystem)**3)
     &           / OvIntThreshold**2)
     &       /SmallestGaussianExponentInSystem)

        CoulIntMaxR = 4.32*2/sqrt(SmallestGaussianExponentInSystem) !<Hard coded to threshold of 1d-9

C       Construct voltages
c$omp   parallel do private(R, Overlap, ThisCharge) schedule(static, 32)
        do i1 = 1,N
           ThisCharge = ZERO
           do i2 = 1,N
              if (i1.ne.i2) then
                 R = Distance(Mol%Atoms(i1)%Basis%Position,
     &                Mol%Atoms(i2)%Basis%Position)

                 if (R .lt. OvIntMaxR) then
                    Overlap = sGTOOvInt(
     &                   Mol%Atoms(i1)%Basis%zeta,
     &                   Mol%Atoms(i2)%Basis%zeta, R)

                 ThisCharge = ThisCharge - Overlap
     &                *(Mol%Atoms(i1)%Element%Electronegativity
     &                -Mol%Atoms(i2)%Element%Electronegativity)
                 end if
              end if
           end do
           Voltage(i1) = ThisCharge/N
        end do
c$omp   end parallel do

c       Put in charge constraint
        Voltage(N+1) = Mol%TotalCharge

c       Use internal conjugate gradients routine
        call dsolver(Mol, Voltage, N+1)

C       Copy out solution from work into charges
        Mol%Atoms(1:N)%Charge = Charge(1:N)
c       print *, "Chemical potential =", Charge(N+1)

C       Calculate energy
        Mol%Energy = 0.0d0
        do i1=1,N
           ThisCharge = Mol%Atoms(i1)%Charge
           Mol%Energy = Mol%Energy - ThisCharge * Voltage(i1)
           do i2=1,i1-1
c             Calculate the contribution to the electrostatic energy. If
c             we are interfacing with TINKER, remember to turn off
c             corresponding calculation in TINKER to avoid double
c             counting
c           if (i1.ne.i2) then
              Mol%Energy = Mol%Energy + ThisCharge
     &             * Mol%Atoms(i2)%Charge * Mol%Coulomb(i1, i2)
c           end if
           end do
           Mol%Energy = Mol%Energy + 0.5d0 * ThisCharge
     &         * ThisCharge * Mol%Coulomb(i1, i1)
        end do

C       Deallocate scratch memory
        deallocate(Voltage, Charge, STAT=stat)
        if (stat .ne. 0) then
          print *, "QTPIE: ERROR: could not deallocate memory"
          stop
        end if
      end subroutine dQTPIE

      subroutine dsolver(Mol, b, N, x)
      use Parameters
      implicit none

      integer, intent(in) :: N
      type(Molecule), intent(in) :: Mol
      real*8, dimension(N), intent(in) :: b
      real*8, dimension(N), intent(out) :: x

      integer :: max_k = 10 !< maximum number of iterations
      real*8, parameter :: tol = 1.0d-8 !< Convergence tolerance

c     Iteration and row loop counters
      integer :: k,i

c     Residual vector, p, q, z, tmp
      real*8, dimension(N) :: r, p, q, z, tmp
      real*8 :: alpha, norm, critical_norm, gamma, gamma0
c     Use the BLAS routines dcopy, ddot, dnrm2, daxpy
      real*8, external :: ddot, dnrm2
      external :: dcopy, daxpy

      external :: MultiplyByA

      logical, parameter :: Verbose = .True.

      max_k = max(10, N+1)

*     Termination criterion norm
      critical_norm = tol * dnrm2(N,b,1)

*     Calculate initial guess from diagonal part
      do i=1,N-1
            x(i) = b(i) / Mol%Atoms(i)%Element%Hardness
      end do
            x(N) = ONE

*     Calculate residual r = b - Ax
      call MultiplyByA(Mol, x, N, tmp)
      r = b - tmp

c     Calculate norm
      norm = dnrm2(N,r,1)
      if (Verbose) then
         print *, "Iteration",0,":",norm, norm/critical_norm
      end if

      do k=1,max_k
*        Generate preconditioned P z = r
         do i=1,N-1
            z(i) = r(i) / Mol%Atoms(i)%Element%Hardness
         end do
         z(N) = ONE

*        Propagate old vectors
         gamma0 = gamma
c        gamma  = r . z
         gamma  = ddot(N,r,1,z,1)

         if (k.ne.1) then
c           p = z + gamma/gamma0 * p
c           With BLAS, first overwrite z,then copy result from z to p
c           z = z + gamma/gamma0 * p
            call daxpy(N,gamma/gamma0,p,1,z,1)
         end if
c        p = z
         call dcopy(N,z,1,p,1)
*        Form matrix-vector product
c        q = A p
         call MultiplyByA(Mol, p, N, q)

*        Calculate step size
c        alpha = gamma / p.q
         alpha = gamma / ddot(N,p,1,q,1)

*        Propagate by step size
c        x = x + alpha * p
         call daxpy(N, alpha,p,1,x,1)
c        r = r - alpha * q
         call daxpy(N,-alpha,q,1,r,1)

*        Calculate new norm of residual
         norm = dnrm2(N,r,1)

*        If requested, print convergence information
         if (Verbose) then
            print *, "Iteration",k,":",norm, norm/critical_norm
         end if
*        Check termination criterion
c        Done if || r || < tol || b ||
         if (norm.lt.critical_norm) then
            goto 1
         end if
      end do

c     Oops, reached maximum iterations without convergence
      print *, "dcg: maximum iterations reached."
      print *, "WARNING: Solution may not be converged."

c     Finally, return the answer
 1    if (Verbose) then
         print *, "dcg: solution found:"
c         print *, x
         print *, "with residual", norm
      end if
      end subroutine dsolver

!>
!! Calculates the matrix to be solved in QTPIE
!! \f[
!! \left(\begin{array}{cc}
!! \mathbf{J} & \vec{1}\\
!! \vec{1}^{T} & 0\end{array}\right)\left(\begin{array}{c}
!! \vec{x}\\
!! y\end{array}\right)=\left(\begin{array}{c}
!! \mathbf{J}\vec{x}+y\\
!! \vec{1}\cdot\vec{x}\end{array}\right)
!! \f]
!! where \f$J\f$ is the classical Coulomb matrix.
!! To speed this up, employ prescreening
!<
      subroutine MultiplyByA(Mol, X, N, V)
      use Parameters
      implicit none  
      integer, intent(in) :: N
      type(Molecule), intent(in) :: Mol
      double precision, dimension(N), intent(in) :: X
      double precision, dimension(N), intent(out) :: V
      double precision :: R, MatrixElement
      double precision, external :: Distance, sGTOCoulInt
      integer :: i, j

      V = 0.0d0

!     (A 1) (x) = (Ax + y)
!     (1 0) (y) = ( 1.x  )

      do i=1,N-1
         do j=1,N-1
            if (i.eq.j) then
               MatrixElement = Mol%Atoms(i)%Element%Hardness
            else
               R = Distance(Mol%Atoms(i)%Basis%Position,
     &                      Mol%Atoms(j)%Basis%Position)

               if (R .gt. CoulIntMaxR) then
                  MatrixElement = 1.0d0 / R
               else
                  MatrixElement = sGTOCoulInt(
     &                 Mol%Atoms(i)%Basis%zeta,
     &                 Mol%Atoms(j)%Basis%zeta, R)
              end if
            end if

            V(i) = V(i) + MatrixElement * X(j)
         end do
         V(i) = V(i) + X(N)
      end do

      do i=1,N-1
         V(N) = V(N) + X(i)
      end do

      end subroutine MultiplyByA



