!>
!! Populates integral matrices in Mol data type
!!
!! Mol%Coulomb and Mol%Overlap are initialized
!! \param Mol : of the Molecule data type
!<
      subroutine DosGTOIntegrals(Mol)
        use Parameters
        implicit none
        double precision, external :: sGTOCoulInt, sGTOOvInt, Distance
        type(Molecule) :: Mol
        integer :: i1, i2, stat, N

!       Temporary variables for caching
        double precision :: R !< Temporary distance
        double precision :: Integral, Norm !< Temporary integrals
        double precision, dimension(3) :: Pos

        N = Mol%NumAtoms

C       Check if memory for integral matrices have been allocated
        if (.not. allocated(Mol%Coulomb)) then
           allocate(Mol%Coulomb(N+1, N+1), STAT = stat)
           if (stat .ne. 0) then
             print *, "DosGTOIntegrals: Error allocating Coulomb matrix"
             print *, "Error code =",stat
             stop
           end if
        end if
        if (.not. allocated(Mol%Overlap)) then
           allocate(Mol%Overlap(N, N), STAT=stat)
           if (stat .ne. 0) then
             print *, "DosGTOIntegrals: Error allocating Overlap matrix"
             print *, "Error code =",stat
             stop
           end if
        end if
        if (.not. allocated(Mol%OvNorm)) then
           allocate(Mol%OvNorm(N), STAT=stat)
           if (stat .ne. 0) then
             print *, "DosGTOIntegrals: Error allocating OvNorm array"
             print *, "Error code =",stat
             stop
           end if
        end if
        
C       Calculate integral pre-screening thresholds
        do i1 = 1,N
          SmallestGaussianExponentInSystem = min(
     &       SmallestGaussianExponentInSystem,
     &       Mol%Atoms(i1)%Basis%zeta)
        end do

        OvIntMaxR = sqrt(
     &      log( (pi/(2*SmallestGaussianExponentInSystem)**3)
     &           / OvIntThreshold**2)
     &       /SmallestGaussianExponentInSystem)

!       An asymptotic expansion of erfc-1(x) gives this formula
        CoulIntMaxR = 2 * sqrt(-log(CoulIntThreshold)/
     &                    SmallestGaussianExponentInSystem)

C       Populate integral matrices

C       Note: Only the lower (i2<i1) triangle of Mol%Overlap is populated
C             because the rest of the code is written to take advantage of
C             its symmetry. This avoids a significant portion of memory
C             access costs.

c        Mol%Overlap = ZERO
c$omp   parallel do private(R, Integral, Pos) schedule(dynamic, 32)
        do i1 = 1, N
           Pos  = Mol%Atoms(i1)%Basis%Position
           do i2 = 1, i1-1
              R = Distance(Pos, Mol%Atoms(i2)%Basis%Position)

              if (R .gt. CoulIntMaxR) then
                 Integral = 1.0d0 / R
              else
                 Integral = sGTOCoulInt(
     &            Mol%Atoms(i1)%Basis%zeta, Mol%Atoms(i2)%Basis%zeta,
     &            R)
              end if

              Mol%Coulomb(i1, i2) = Integral
              Mol%Coulomb(i2, i1) = Integral

              if (R .gt. OvIntMaxR) then
                 Integral = 0.0d0
              else
                 Integral = sGTOOvInt(
     &            Mol%Atoms(i1)%Basis%zeta, Mol%Atoms(i2)%Basis%zeta,
     &            R )
              end if
              Mol%Overlap(i1, i2) = Integral
c              Mol%Overlap(i2, i1) = Integral

           end do
C          For the diagonal elements, use hardness
           Mol%Coulomb(i1, i1) = Mol%Atoms(i1)%Element%Hardness
c           Mol%Overlap(i1, i1) = ONE
        end do
c$omp   end parallel do

c       Calculate due normalization
c$omp   parallel do private(Norm) schedule(static, 32)
        do i1 = 1, N
           Norm = ONE
           do i2 = 1, i1 - 1
              Norm = Norm + Mol%Overlap(i1, i2)
           end do
           do i2 = i1 + 1, N
              Norm = Norm + Mol%Overlap(i2, i1)
           end do
           Mol%OvNorm(i1) = N / Norm
        end do
c$omp   end parallel do
      end subroutine DosGTOIntegrals

!>
!! Populates integral matrices in Mol data type
!!
!! Mol%Coulomb and Mol%Overlap are initialized
!! \param Mol : of the Molecule data type
!<
      subroutine DosSTOIntegrals(Mol)
        use Parameters
        implicit none
        double precision, external :: sSTOCoulInt, sSTOOvInt, Distance
        type(Molecule) :: Mol
        integer :: i1, i2, stat

        double precision, dimension(:,:), allocatable :: RefOverlap

C       Check if memory for integral matrices have been allocated
        if (.not. allocated(Mol%Coulomb)) then
           allocate(Mol%Coulomb(Mol%NumAtoms, Mol%NumAtoms), STAT=stat)
           if (stat .ne. 0) then
             print *, "DosSTOIntegrals: Error allocating Coulomb matrix"
             stop
           end if
        end if
        if (.not. allocated(Mol%Overlap)) then
           allocate(Mol%Overlap(Mol%NumAtoms, Mol%NumAtoms), STAT=stat)
           if (stat .ne. 0) then
             print *, "DosSTOIntegrals: Error allocating Overlap matrix"
             stop
           end if
        end if
        if (.not. allocated(Mol%OvNorm)) then
           allocate(Mol%OvNorm(Mol%NumAtoms), STAT=stat)
           if (stat .ne. 0) then
             print *, "DosSTOIntegrals: Error allocating OvNorm array"
             stop
           end if
        end if
C       Allocate memory for reference overlaps
        allocate(RefOverlap(Mol%NumAtoms, Mol%NumAtoms), STAT=stat)
        if (stat .ne. 0) then
          print *, "DosSTOIntegrals: Error allocating RefOverlap matrix"
          stop
        end if

C       Now compute Coulomb matrix
        do i1 = 1,Mol%NumAtoms
           do i2 = 1, i1-1
                 Mol%Coulomb(i1, i2) = sSTOCoulInt(
     &               Mol%Atoms(i1)%Basis%zeta, Mol%Atoms(i2)%Basis%zeta,
     &               Mol%Atoms(i1)%Basis%n , Mol%Atoms(i2)%Basis%n,
     &               Distance(Mol%Atoms(i1)%Basis%Position,
     &                        Mol%Atoms(i2)%Basis%Position) )
c              print *, "Co", i1, i2, Mol%Coulomb(i1,i2)
C             Fill in the other triangle
              Mol%Coulomb(i2, i1) = Mol%Coulomb(i1,i2)
           end do
C          For the diagonal elements, use hardness
              Mol%Coulomb(i1, i1) = Mol%Atoms(i1)%Element%Hardness
        end do
        
C       Now compute Overlap and RefOverlap matrices
        do i1 = 1,Mol%NumAtoms
           do i2 = 1, i1-1
              Mol%Overlap(i1, i2) = sSTOOvInt(
     &            Mol%Atoms(i1)%Basis%zeta, Mol%Atoms(i2)%Basis%zeta,
     &            Mol%Atoms(i1)%Basis%n , Mol%Atoms(i2)%Basis%n,
     &            Distance(Mol%Atoms(i1)%Basis%Position,
     &                     Mol%Atoms(i2)%Basis%Position) )

C             Calculate the same quantity but referenced to an intrinsic
C             length scale
              RefOverlap(i1, i2) = sSTOOvInt(
     &            Mol%Atoms(i1)%Basis%zeta, Mol%Atoms(i2)%Basis%zeta,
     &            Mol%Atoms(i1)%Basis%n , Mol%Atoms(i2)%Basis%n,
     &            ExpectR(Mol%Atoms(i1)%Basis) 
     &           +ExpectR(Mol%Atoms(i2)%Basis) )

C             Fill in the other triangle
c              print *, "Ov", i1, i2, Mol%Overlap(i1,i2)
              Mol%Overlap(i2, i1) = Mol%Overlap(i1, i2)
              RefOverlap(i2, i1) = RefOverlap(i1, i2)
           end do
C          For the diagonal elements, the overlap is just the orbital normalization
           Mol%Overlap(i1, i1) = 1.0d0
           RefOverlap(i1, i1) = 1.0d0
        end do

C     Now compute normalization of Attenuation (overlap) matrix
      do i1 = 1,Mol%NumAtoms
        Mol%OvNorm(i1) = 0.0d0
        do i2 = 1,Mol%NumAtoms
          Mol%OvNorm(i1) = Mol%OvNorm(i1) + RefOverlap(i1, i2)
        end do
        Mol%OvNorm(i1) = Mol%OvNorm(i1) / Mol%NumAtoms
      end do

C     Multiply in Norm
      do i1 = 1,Mol%NumAtoms
        do i2 = 1,Mol%NumAtoms
          Mol%Overlap(i1,i2) = Mol%Overlap(i1,i2) / Mol%OvNorm(i1)
        end do
      end do

C     Deallocate temporary variables
      deallocate(RefOverlap, STAT=stat)
      end subroutine DosSTOIntegrals

!>
!! Populates atomic charges according to the QEq(-H) charge model
!! \param Mol : of the Molecule data type
!! Mol%Atoms(i)%Charge are computed
!> \note The model is described in the seminal paper below:
!!       "Charge equilibration for Molecular dynamics simulations"
!!       A. K. Rappe and W. A. Goddard, J. Phys. Chem., 1991, 95(8), 3358-3363
!!       doi:10.1021/j100161a070
!> \note This implementation does not do the additional procedure for H atoms
!!       nor does it check for overly large charges that exceed the principal
!!      quantum number of the given atom.
!<
      subroutine QEq(Mol)
        use Parameters
        implicit none
        type(Molecule) :: Mol
        double precision, dimension(:,:), allocatable :: Capacitance
        double precision, dimension(:), allocatable :: Voltage, Charge
        integer :: i1, i2, stat

        integer :: N !< size of problem

C       Wrapper for linear algebra solver
        external :: solver

        N = Mol%NumAtoms

*       Allocate scratch memory
        allocate(Capacitance(N+1,N+1),Voltage(N+1), Charge(N+1),
     &           STAT=stat)
        if (stat .ne. 0) then
           print *,"QEq: ERROR: could not allocate sufficient memory"
           stop
        end if

*       Calculate problem
C       Construct the capacitance matrix
        Capacitance(1:N,1:N) = Mol%Coulomb
        Capacitance = Mol%Coulomb
c       Add charge conservation constraints
        Capacitance(1:N,N+1) = ONE
        Capacitance(N+1,1:N) = ONE
        Capacitance(N+1,N+1) = ZERO

c       Construct voltage matrix
        Voltage(1:N) = -Mol%Atoms(1:N)%Element%Electronegativity
c       The last row expresses charge conservation
        Voltage(N+1) = Mol%TotalCharge

*       Solve for Charge, where Capacitance * Charge = Voltage
c       i.e. an A*x = b problem where A is symmetric and positive nonnegative

c       Call linear algebra solver
        call solver(N+1, Capacitance, Voltage, Charge)

c       Copy charges into the Mol structure
        Mol%Atoms(1:N)%Charge = Charge(1:N)
c       print *, "Chemical potential =",Charge(N+1)

*       Calculate energy
        Mol%Energy = 0.0d0
        do i1=1,N
           Mol%Energy = Mol%Energy + Mol%Atoms(i1)%Charge
     &          * Mol%Atoms(i1)%Element%Electronegativity
           do i2=1,N
c             Calculate the contribution to the electrostatic energy. If
c             we are interfacing with TINKER, remember to turn off
c             corresponding calculation in TINKER to avoid double
c             counting
              Mol%Energy = Mol%Energy + 0.5d0 * Mol%Atoms(i1)%Charge
     &             * Mol%Atoms(i2)%Charge * Mol%Coulomb(i1, i2)
           end do
        end do

*       Done. Cleanup.
c       Deallocate scratch memory
        deallocate(Capacitance, Voltage, Charge, STAT=stat)
        if (stat .ne. 0) then
          print *, "QEq: ERROR: could not deallocate memory"
          stop
        end if 

      end subroutine QEq

!>
!! Populates atomic charges according to the QTPIE charge model
!! \param Mol : of the Molecule data type
!! Mol%Atoms(i)%Charge are computed
!! \note The model is described in the paper below:
!!       J. Chen and T. J. Martinez, Chem. Phys. Lett., 438 (4-6), 2007, 315-320
!!       doi:10.1016/j.cplett.2007.02.065
!<
      subroutine QTPIE(Mol)
        use Parameters
        implicit none
        type(Molecule) :: Mol

        double precision, dimension(:), allocatable :: Voltage, Charge
        save Charge

        double precision :: ThisCharge !< Temporary atomic charge variable
        double precision :: Overlap    !< Temporary integral
        double precision :: VoltageDifference, Norm
C       i1-i2 loop over atoms
        integer :: stat, i1, i2

        integer :: N !< size of matrix problem

C       Wrapper for linear algebra solver
        external :: solver

C       Define size of matrix problem
        N = Mol%NumAtoms

C       Allocate scratch memory
        if (.not. allocated(Voltage)) then
           allocate(Voltage(N+1), STAT=stat)
           if (stat .ne. 0) then
              print *,"QTPIE: ERROR: could not allocate memory"
              stop
           end if
        end if

        if (.not. allocated(Charge)) then
           allocate(Charge(N+1), STAT=stat)
           Charge=ZERO
           if (stat .ne. 0) then
              print *,"QTPIE: ERROR: could not allocate memory"
              stop
           end if
        end if

*       Construct problem

c       Add charge conservation constraints
        Mol%Coulomb(1:N,N+1) = ONE
        Mol%Coulomb(N+1,1:N) = ONE
        Mol%Coulomb(N+1,N+1) = ZERO
        
C       Construct voltages
c       The code here is a little convoluted but knows that Overlap
c       is sparse and symmetric

c$omp   parallel do private(ThisCharge, Overlap, Norm,
c$omp&     VoltageDifference, i2) schedule(dynamic, 32)
        do i1 = 1,N
           ThisCharge = ZERO
           Norm = Mol%OvNorm(i1)
           do i2 = 1,i1-1
              Overlap = Mol%Overlap(i1,i2)
              if (Overlap.ne.ZERO) then
                 VoltageDifference = 
     &                ( Mol%Atoms(i1)%Element%Electronegativity
     &                - Mol%Atoms(i2)%Element%Electronegativity)
                 if (VoltageDifference.ne.ZERO) then
                   ThisCharge = ThisCharge - Overlap * VoltageDifference
     &                   * Norm
                 end if
              end if
           end do
           do i2 = i1+1,N
              Overlap = Mol%Overlap(i2,i1)
              if (Overlap.ne.ZERO) then
                 VoltageDifference = 
     &                ( Mol%Atoms(i1)%Element%Electronegativity
     &                - Mol%Atoms(i2)%Element%Electronegativity)
                 if (VoltageDifference.ne.ZERO) then
                   ThisCharge = ThisCharge - Overlap * VoltageDifference
     &                   * Norm
                 end if
              end if
           end do
           Voltage(i1) = ThisCharge/N
        end do
c$omp   end parallel do

c       Put in charge constraint
        Voltage(N+1) = Mol%TotalCharge

c       Use internal conjugate gradients routine
        call solver(N+1, Mol%Coulomb, Voltage, Charge)

C       Copy out solution from work into charges
        Mol%Atoms(1:N)%Charge = Charge(1:N)
c       print *, "Chemical potential =", Charge(N+1)

C       Calculate energy
        Mol%Energy = 0.0d0
        do i1=1,N
           ThisCharge = Mol%Atoms(i1)%Charge
           Mol%Energy = Mol%Energy - ThisCharge * Voltage(i1)
           do i2=1,i1-1
              Mol%Energy = Mol%Energy + ThisCharge
     &             * Mol%Atoms(i2)%Charge * Mol%Coulomb(i1, i2)
           end do
           Mol%Energy = Mol%Energy + 0.5d0 * ThisCharge
     &         * ThisCharge * Mol%Coulomb(i1, i1)
        end do

C       Deallocate scratch memory
!        deallocate(Voltage, Charge, STAT=stat)
!        if (stat .ne. 0) then
!          print *, "QTPIE: ERROR: could not deallocate memory"
!          stop
!        end if 

      end subroutine QTPIE

!>
!! Computes energy gradients numerically
!!
!! Calculates energy gradients using the method of finite differences
!! using forward gradients
!! As you can imagine, this is pretty slow
!! You should not use this routine!
!! \param Mol : of the Molecule data type
!! Mol%EGradient is calculated
!<
      subroutine DoGradientsByFiniteDifference(Mol)
        use Parameters
        implicit none
        type(Molecule) :: Mol
        integer :: i1, i2, stat
        double precision :: OriginalEnergy
        double precision, parameter :: Eps = 1.0d-3

C       Check if memory for gradient matrix has been allocated
        if (.not. allocated(Mol%EGradient)) then
           allocate(Mol%EGradient(Mol%NumAtoms, 3), STAT=stat)
           if (stat .ne. 0) then
              print *, "DoGradientsByFiniteDifference: allocation error"
              stop
           end if
        end if

C       Save current energy
        OriginalEnergy = Mol%Energy
C       Calculate energy gradients
        do i1=1,Mol%NumAtoms
          do i2=1,3
C           Perturb Geometry
            Mol%Atoms(i1)%Position(i2) =
     &      Mol%Atoms(i1)%Position(i2) + Eps
            Mol%Atoms(i1)%Basis%Position(i2) =
     &      Mol%Atoms(i1)%Basis%Position(i2) + Eps
C           Redo QTPIE
            call DosGTOIntegrals(Mol)
            call QTPIE(Mol)
C           Calculate gradient
            Mol%EGradient(i1, i2) =
     &           (Mol%Energy - OriginalEnergy)  / ( Eps) 
C           Perturb Geometry
            Mol%Atoms(i1)%Position(i2) =
     &      Mol%Atoms(i1)%Position(i2) - Eps
            Mol%Atoms(i1)%Basis%Position(i2) =
     &      Mol%Atoms(i1)%Basis%Position(i2) - Eps
          end do
        end do
C       Redo integrals
        call DosGTOIntegrals(Mol)
      end subroutine DoGradientsByFiniteDifference

!>
!! Computes energy gradients analytically
!!
!! Calculates energy gradients using analytic derivatives
!! \param Mol : of the Molecule data type
!! Mol%EGradient is calculated
!<
      subroutine DoGradientsAnalytically(Mol)
        use Parameters
        implicit none
        type(Molecule) :: Mol
        double precision :: a,b,R, Force
        double precision, external :: Distance
        double precision, external :: sGTOOvIntGrad, sGTOCoulIntGrad
        integer :: i1,i2,i3,m,n,stat

C       Check if memory for gradient matrix has been allocated
        if (.not. allocated(Mol%EGradient)) then
           allocate(Mol%EGradient(Mol%NumAtoms, 3), STAT=stat)
           if (stat .ne. 0) then
              print *, "cannot allocate memory for energy gradients"
              stop
           end if
        end if

c       Initialize gradients
        Mol%EGradient=0.0d0

        do i1=1,Mol%NumAtoms
c          Obtain basis set parameters for atom i1
           a = Mol%Atoms(i1)%Basis%zeta
           m = Mol%Atoms(i1)%Basis%n

c          Calculate contribution to gradient from voltage term
           do i2=1,Mol%NumAtoms
c          Diagonal part has no contribution to gradient
           if (i1.ne.i2) then 
c             Obtain basis set parameters for atom i2
              b = Mol%Atoms(i2)%Basis%zeta
              n = Mol%Atoms(i2)%Basis%n
c             Calculate pairwise distance
              R = Distance(Mol%Atoms(i1)%Basis%Position,
     &                     Mol%Atoms(i2)%Basis%Position)
              Force = 2 * Mol%Atoms(i1)%Charge / Mol%NumAtoms
     &          * ( Mol%Atoms(i1)%Element%Electronegativity
     &            - Mol%Atoms(i2)%Element%Electronegativity )
     &          * sGTOOvIntGrad(a,b,R) / Mol%OvNorm(i1)
              Force = Force - Mol%Atoms(i1)%Charge / Mol%NumAtoms
     &          * ( Mol%Atoms(i1)%Element%Electronegativity
     &            - Mol%Atoms(i2)%Element%Electronegativity )
     &          * Mol%Overlap(i1,i2) / (Mol%OvNorm(i1) ** 2)
     &          * sGTOOvIntGrad(a,b,R)
c     &          * (1/Mol%OvNorm(i1) - 1/Mol%OvNorm(i2)) 
             Force = Force + Mol%Atoms(i1)%Charge * Mol%Atoms(i2)%Charge
     &                 * sGTOCoulIntGrad(a,b,R)
c             Calculates projection onto direction vector
c             $Temp*\frac{\partial R_{i1,i2}}{\partial R_{k,i3}}
c             * (\delta_{i1,k} - \delta_{i2,k})$
              Force = Force / R
              do i3=1,3
                 Mol%EGradient(i1,i3)=Mol%EGradient(i1,i3) +
     &                (Mol%Atoms(i1)%Basis%Position(i3)-
     &                 Mol%Atoms(i2)%Basis%Position(i3))*Force
              end do
           end if
           end do
        end do
      end subroutine DoGradientsAnalytically
