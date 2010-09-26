!>
!! Runs QTPIE for a single XYZ geometry
!<
      program onexyz 
        use Parameters
        implicit none
        type(Molecule) :: Mol
        type(Molecule), external :: loadXYZ
        integer :: NumArgs
        intrinsic :: iargc
        character (len = 50) :: fileName

C        intrinsic :: etime
        real elapsed(2)
        real total, old
        total = 0.0
        print *, "Single geometry mode"
        NumArgs = iargc()
        if (NumArgs.ge.1) then
           call getarg(1, fileName)
           print *, "Reading in file ", fileName
        else
           print *, "Reading default file name qtpie.xyz"
           fileName = "qtpie.xyz"
        end if

        old = 0.0d0
        Mol =  loadXYZ(fileName)
C        total = etime(elapsed)
        print *, "Read file", total-old

        old = total
        call DosGTOIntegrals(Mol)
C        total = etime(elapsed)
        print *, "Integrals", total-old

        old = total
        call QTPIE(Mol)
C        total = etime(elapsed)
        print *, "Calculated", total-old

        print *, "QTPIE Energy is", Mol%Energy

!        old = total
!        call dQTPIE(Mol)
!        total = etime(elapsed)
!        print *, "Calculated (direct)", total-old
!       print *, "QTPIE Energy is", Mol%Energy

        call WriteLog(Mol, "qtpie.log")
        print *, "Calculated charges written to qtpie.log"

      stop
      print *, "Numerical forces"
      call DoGradientsByFiniteDifference(Mol)
      print *, Mol%EGradient
      print *, "Analytic forces"
      call DoGradientsAnalytically(Mol)
      print *, Mol%EGradient
      
      end program onexyz
