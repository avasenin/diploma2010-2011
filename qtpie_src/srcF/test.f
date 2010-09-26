!>
!! A simple program to test the functions implemented with some test values
!<
      program test
        use Parameters
        use factorial
        implicit none
        double precision, external :: RosenA, RosenB
        integer, external :: RosenD
        double precision, external :: sSTOCoulInt, sSTOOvInt
        double precision, external :: sSTOCoulIntGrad, sSTOOvIntGrad
        double precision, external :: sGTOCoulInt, sGTOOvInt
        double precision, external :: sGTOCoulIntGrad, sGTOOvIntGrad
        type(Molecule) :: Mol1, Mol2
        type(Molecule), external :: loadXYZ
        
        real*8, parameter :: epsilon = 1.0d-6

        print *, "Testing mode"

        if (fact(6).eq.720) then
          print *, "Factorial correct"
        else
          print *, "FATAL ERROR: Factorials incorrectly computed"
          print *, "6! = ", fact(6), ", expected 720"
          print *, "There is an error in factorial.f"
          stop
        end if  

        if (fact(10).eq.3628800) then
          print *, "Factorial correct"
        else
          print *, "FATAL ERROR: Factorials incorrectly computed"
          print *, "10! = ", fact(10), ", expected 3628800"
          print *, "There is an error in factorial.f"
          stop
        end if  
        
        if (abs(RosenA(4,3.0d0)-8.05198d-2).lt.epsilon) then
          print *, "Rosen A integral correct"
        else
          print *, "FATAL ERROR: Rosen A integral incorrectly computed"
          print *, "RosenA(4,3.0) = ", RosenA(4,3.0d0),
     &         "expected 0.0805198"
          print *, "There is an error in sto-int.f"
          stop
        end if  

        if (abs(RosenA(12,10.0d0)-3.79157d-5).lt.epsilon) then
          print *, "Rosen A integral correct"
        else
          print *, "FATAL ERROR: Rosen A integral incorrectly computed"
          print *, "RosenA(12,10.0) = ", RosenA(12,10.0d0), 
     &             ", expected", 3.79157d-5
          print *, "There is an error in sto-int.f"
          stop
        end if  

        if (abs(RosenB(4,3.0d0)-2.6471457).lt.epsilon) then
          print *, "Rosen B integral correct"
        else
          print *, "FATAL ERROR: Rosen B integral incorrectly computed"
          print *, "RosenB(4,3.0) = ", RosenB(4,3.0d0), ", expected",
     &      2.6471457
          print *, "There is an error in sto-int.f"
          stop
        end if  

        if (abs(RosenB(8,3.0d0)-1.715602).lt.epsilon) then
          print *, "Rosen B integral correct"
        else
          print *, "FATAL ERROR: Rosen B integral incorrectly computed"
          print *, "RosenB(8,3.0) = ", RosenB(8,3.0d0), ", expected",
     &      3.75628-2.04068
          print *, "There is an error in sto-int.f"
          stop
        end if  

        if (abs(RosenB(12,10.0d0)-9.759958896510301d2).lt.epsilon) then
          print *, "Rosen B integral correct"
        else
          print *, "FATAL ERROR: Rosen B integral incorrectly computed"
          print *, "RosenB(12,10.0) = ", RosenB(12,10.0d0), 
     &         ", expected", 9.75996d2-3.79157d-5
          print *, "There is an error in sto-int.f"
          stop
        end if  

        if (RosenD(1,2,3).eq.1) then
          print *, "Rosen D factor correct"
        else
          print *, "FATAL ERROR: Rosen D factor incorrectly computed"
          print *, "RosenD(1,2,3) = ", RosenD(1,2,3), ", expected 1"
          print *, "There is an error in sto-int.f"
          stop
        end if  

        if (RosenD(5,4,3).eq.-4) then
          print *, "Rosen D factor correct"
        else
          print *, "FATAL ERROR: Rosen D factor incorrectly computed"
          print *, "RosenD(5,4,3) = ", RosenD(5,4,3), ", expected -4"
          print *, "There is an error in sto-int.f"
          stop
        end if  

        if (RosenD(5,3,8).eq.-1) then
          print *, "Rosen D factor correct"
        else
          print *, "FATAL ERROR: Rosen D factor incorrectly computed"
          print *, "RosenD(5,3,8) = ", RosenD(5,3,8), ", expected -1"
          print *, "There is an error in sto-int.f"
          stop
        end if  

        if (abs(sSTOCoulInt(1.0d0,2.0d0,3,4,5.0d0)-0.1903871)
     &       .lt.epsilon) then
          print *, "Coulomb integral correct"
        else
          print *, "FATAL ERROR: Coulomb integral incorrectly computed"
          print *, "Coulomb(1,2,3,4,5) = ",
     &         sSTOCoulInt(1.0d0,2.0d0,3,4,5.0d0),
     &      ", expected", 0.1903871
          print *, "There is an error in sto-int.f"
          stop
        end if  

        if (abs(sSTOCoulInt(1.0d0,1.0d0,2,3,5.0d0)-0.1879457)
     &       .lt.epsilon) then
          print *, "Coulomb integral correct"
        else
          print *, "FATAL ERROR: Coulomb integral incorrectly computed"
          print *, "Coulomb(1,1,2,3,5) = ",
     &         sSTOCoulInt(1.0d0,1.0d0,2,3,5.0d0),
     &         ", expected", 0.1879457
          print *, "There is an error in sto-int.f"
          stop
        end if  

        if (abs(sSTOCoulInt(5.0d0,4.0d0,3,2,1.0d0)-0.9135013)
     &       .lt.epsilon) then
          print *, "Coulomb integral correct"
        else
          print *, "FATAL ERROR: Coulomb integral incorrectly computed"
          print *, "Coulomb(5,4,3,2,1) = ",
     &         sSTOCoulInt(5.0d0,4.0d0,3,2,1.0d0),
     &         ", expected", 0.9135013
          print *, "There is an error in sto-int.f"
          stop
        end if  

        if (abs(sSTOOvInt(1.0d0,2.0d0,3,4,5.0d0)-0.3145446)
     &       .lt.epsilon) then
          print *, "Overlap integral correct"
        else
          print *, "FATAL ERROR: Overlap integral incorrectly computed"
          print *, "Overlap(1,2,3,4,5) = ",
     &         sSTOOvInt(1.0d0,2.0d0,3,4,5.0d0),
     &      ", expected", 0.3145446
          print *, "There is an error in sto-int.f"
          stop
        end if  

        if (abs(sSTOOvInt(1.0d0,1.0d0,2,3,5.0d0)-0.3991235)
     &       .lt.epsilon) then
          print *, "Overlap integral correct"
        else
          print *, "FATAL ERROR: Overlap integral incorrectly computed"
          print *, "Overlap(1,1,2,3,5) = ",
     &         sSTOOvInt(1.0d0,1.0d0,2,3,5.0d0),
     &      ", expected", 0.3991235
          print *, "There is an error in sto-int.f"
          stop
        end if  

        Mol1 =  loadXYZ("../test/nacl.xyz")
        print *, "Load XYZ successful"

c        call DosSTOIntegrals(Mol1)
        call DosGTOIntegrals(Mol1)
        call QEq(Mol1)
 
        if (abs(Mol1%Atoms(1)%Charge - 1.3895802392931).lt.epsilon) then
          print *, "QEq Charge calculation for sodium chloride correct"
        else
          print *, "FATAL ERROR: QEq charges incorrectly computed"
          print *, "Charge(1) = ", Mol1%Atoms(1)%Charge,
     &         " expected", 1.3895802392931755
          print *, "There is an error in qtpie.f"
c          stop
        end if

        call QTPIE(Mol1)
        if (abs(Mol1%Atoms(1)%Charge - 0.7252290067905).lt.epsilon) then
          print *, "QTPIE Charge calculation for ",
     &          "sodium chloride correct"
        else
          print *, "FATAL ERROR: QTPIE charges incorrectly computed"
          print *, "Charge(1) = ", Mol1%Atoms(1)%Charge,
     &         " expected", 0.72522900679059155
          print *, "There is an error in qtpie.f"
c          stop
        end if

        Mol2 =  loadXYZ("../test/h2o.xyz")
        print *, "Load XYZ successful"

c        call DosSTOIntegrals(Mol2)
        call DosGTOIntegrals(Mol2)
        call QEq(Mol2)

        if ((abs(Mol2%Atoms(1)%Charge + 0.98965172663).lt.epsilon) .and.
     &      (abs(Mol2%Atoms(2)%Charge - 0.4943811799).lt. epsilon)) then
          print *, "QEq Charge calculation for water correct"
        else
          print *, "FATAL ERROR: QEq charges incorrectly computed"
          print *, "Charge(1) = ", Mol2%Atoms(1)%Charge,
     &         " expected", -0.98965172663781498
          print *, "Charge(2) = ", Mol2%Atoms(2)%Charge,
     &         " expected",  0.49438117990925540
          print *, "There is an error in qtpie.f"
c          stop
        end if

        call QTPIE(Mol2)

        if ((abs(Mol2%Atoms(1)%Charge + 0.81213640965).lt.epsilon) .and.
     &      (abs(Mol2%Atoms(2)%Charge - 0.4055667959).lt. epsilon)) then
          print *, "QTPIE Charge calculation for water correct"
        else
          print *, "FATAL ERROR: QTPIE charges incorrectly computed"
          print *, "Charge(1) = ", Mol2%Atoms(1)%Charge,
     &         " expected", -0.40556679590604666
          print *, "Charge(2) = ", Mol2%Atoms(2)%Charge,
     &         " expected",  0.40556679590604666
          print *, "There is an error in qtpie.f"
c          stop
        end if

*     Now compare numerical and analytic gradients

      print *, "Testing gradients for Slater orbitals"

      if (abs(sSTOOvIntGrad(1.0d0,1.0d0,3,4,5.0d0) -
     &     (sSTOOvInt(1.0d0,1.0d0,3,4,5.0d0+epsilon)
     &     -sSTOOvInt(1.0d0,1.0d0,3,4,5.0d0-epsilon))
     &     /(2*epsilon)).lt.epsilon) then 
         print *, "Overlap gradient correct"
      else
         print *, "FATAL ERROR: Overlap gradients don't match"
         print *, "sSTOOvIntGrad(1.0, 1.0, 3, 4, 5.0)"
         print *, "Analytic = ", sSTOOvIntGrad(1.0d0,1.0d0,3,4,5.0d0)
         print *, "Numerical= ",
     &         (sSTOOvInt(1.0d0,1.0d0,3,4,5.0d0+epsilon)
     &         -sSTOOvInt(1.0d0,1.0d0,3,4,5.0d0-epsilon))/(2*epsilon) 
          print *, "There is an error in sto-int.f"
      end if

      if (abs(sSTOOvIntGrad(1.0d0,2.0d0,3,4,5.0d0) -
     &     (sSTOOvInt(1.0d0,2.0d0,3,4,5.0d0+epsilon)
     &     -sSTOOvInt(1.0d0,2.0d0,3,4,5.0d0-epsilon))
     &     /(2*epsilon)).lt.epsilon) then 
         print *, "Overlap gradient correct"
      else
         print *, "FATAL ERROR: Overlap gradients don't match"
         print *, "sSTOOvIntGrad(1.0, 2.0, 3, 4, 5.0)"
         print *, "Analytic = ", sSTOOvIntGrad(1.0d0,2.0d0,3,4,5.0d0)
         print *, "Numerical= ",
     &         (sSTOOvInt(1.0d0,2.0d0,3,4,5.0d0+epsilon)
     &         -sSTOOvInt(1.0d0,2.0d0,3,4,5.0d0-epsilon))/(2*epsilon) 
          print *, "There is an error in sto-int.f"
      end if

      if (abs(sSTOCoulIntGrad(4.0d0,4.0d0,2,2,5.0d0) -
     &     (sSTOCoulInt(4.0d0,4.0d0,2,2,5.0d0+epsilon)
     &     -sSTOCoulInt(4.0d0,4.0d0,2,2,5.0d0-epsilon))
     &     /(2*epsilon)).lt.epsilon) then 
         print *, "Coulomb gradient correct"
      else
         print *, "FATAL ERROR: Coulomb gradients don't match"
         print *, "sSTOCoulIntGrad(4.0, 4.0, 2, 2, 5.0)"
         print *, "Analytic = ", sSTOCoulIntGrad(4.0d0,4.0d0,2,2,5.0d0)
         print *, "Numerical= ",
     &         (sSTOCoulInt(4.0d0,4.0d0,2,2,5.0d0+epsilon)
     &         -sSTOCoulInt(4.0d0,4.0d0,2,2,5.0d0-epsilon))/(2*epsilon) 
          print *, "There is an error in sto-int.f"
      end if

      if (abs(sSTOCoulIntGrad(4.0d0,4.0d0,3,4,5.0d0) -
     &     (sSTOCoulInt(4.0d0,4.0d0,3,4,5.0d0+epsilon)
     &     -sSTOCoulInt(4.0d0,4.0d0,3,4,5.0d0-epsilon))
     &     /(2*epsilon)).lt.epsilon) then 
         print *, "Coulomb gradient correct"
      else
         print *, "FATAL ERROR: Coulomb gradients don't match"
         print *, "sSTOCoulIntGrad(4.0, 4.0, 3, 4, 5.0)"
         print *, "Analytic = ", sSTOCoulIntGrad(4.0d0,4.0d0,3,4,5.0d0)
         print *, "Numerical= ",
     &         (sSTOCoulInt(4.0d0,4.0d0,3,4,5.0d0+epsilon)
     &         -sSTOCoulInt(4.0d0,4.0d0,3,4,5.0d0-epsilon))/(2*epsilon) 
          print *, "There is an error in sto-int.f"
      end if

      if (abs(sSTOCoulIntGrad(1.0d0,2.0d0,3,4,5.0d0) -
     &     (sSTOCoulInt(1.0d0,2.0d0,3,4,5.0d0+epsilon)
     &     -sSTOCoulInt(1.0d0,2.0d0,3,4,5.0d0-epsilon))
     &     /(2*epsilon)).lt.epsilon) then 
         print *, "Coulomb gradient correct"
      else
         print *, "FATAL ERROR: Coulomb gradients don't match"
         print *, "sSTOCoulIntGrad(1.0, 2.0, 3, 4, 5.0)"
         print *, "Analytic = ", sSTOCoulIntGrad(1.0d0,2.0d0,3,4,5.0d0)
         print *, "Numerical= ",
     &         (sSTOCoulInt(1.0d0,2.0d0,3,4,5.0d0+epsilon)
     &         -sSTOCoulInt(1.0d0,2.0d0,3,4,5.0d0-epsilon))/(2*epsilon) 
          print *, "There is an error in sto-int.f"
      end if

      print *, "Testing gradients for Gaussian orbitals"

      if (abs(sGTOOvIntGrad(1.0d0,1.0d0,2.0d0) -
     &     (sGTOOvInt(1.0d0,1.0d0,2.0d0+epsilon)
     &     -sGTOOvInt(1.0d0,1.0d0,2.0d0-epsilon))
     &     /(2*epsilon)).lt.epsilon) then 
         print *, "Overlap gradient correct"
      else
         print *, "FATAL ERROR: Overlap gradients don't match"
         stop
      end if

      if (abs(sGTOCoulIntGrad(1.0d0,1.0d0,2.0d0) -
     &     (sGTOCoulInt(1.0d0,1.0d0,2.0d0+epsilon)
     &     -sGTOCoulInt(1.0d0,1.0d0,2.0d0-epsilon))
     &     /(2*epsilon)).lt.epsilon) then 
         print *, "Coulomb gradient correct"
      else
         print *, "FATAL ERROR: Coulomb gradients don't match"
         print *, sGTOCoulIntGrad(1.0d0, 1.0d0, 2.0d0)
         print *, (sGTOCoulInt(1.0d0,1.0d0,2.0d0+10*epsilon)
     &     -sGTOCoulInt(1.0d0,1.0d0,2.0d0-10*epsilon))
     &     /(20*epsilon)
         stop
      end if

      print *, "forces for NaCl"
      print *, "Numerical forces"
      call DoGradientsByFiniteDifference(Mol1)
      print *, Mol1%EGradient
      print *, "Analytic forces"
      call DoGradientsAnalytically(Mol1)
      print *, Mol1%EGradient
      print *, "forces for water"
      print *, "Numerical forces"
      call DoGradientsByFiniteDifference(Mol2)
      print *, Mol2%EGradient
      print *, "Analytic forces"
      call DoGradientsAnalytically(Mol2)
      print *, Mol2%EGradient
      end program test
