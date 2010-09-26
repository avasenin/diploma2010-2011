!>
!! Computes the dipole moment
!! \param Centroid Point about which to define dipole moment
!! \return the dipole moment vector
!<
      function dipmom(Mol)
        use Parameters
        type(Molecule) :: Mol
        double precision, dimension(3) :: Centroid
        double precision, dimension(3) :: dipmom
        integer i,j

        do i=1,Mol%NumAtoms
           do j=1,3
              dipmom(j) = dipmom(j) + Mol%Atoms(i)%Charge *
     &          (Mol%Atoms(i)%Position(j) - Centroid(j))
           end do
        end do
      end function dipmom
      
!>
!! Computes the dipole polarizability tensor
!! \note currently returns 0
!<
      double precision function polarizability()
        double precision, dimension(3:3) :: polarizability
        polarizability = 0
      end function polarizability
