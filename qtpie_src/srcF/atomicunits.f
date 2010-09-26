!>
!! A Fortran module storing conversion factors
!!
!!
!!       Our QTPIE charge model works exclusively in atomic units
!!       The values stored here are conversion factors to convert
!!       from that unit into atomic units
!!       There is no dimensional checking implemented!
!<
      module AtomicUnits
        implicit none
        save
        double precision, parameter :: ONE = 1.0d0
        double precision, parameter :: ZERO = 0.0d0
        double precision, parameter :: eV = 3.67493245d-2 !< electron volt to Hartree
        double precision, parameter :: kJ_mol = 6.6744644952d-3 !<kilojoule per mole to Hartree
        double precision, parameter :: kcal_mol = 1.5952353d-3 !<kilocalorie per mole to Hartree
        double precision, parameter :: invAngstrom = 455.6335252760d0 !<inverse Ångstrom to Hartree
        double precision, parameter :: Debye = 0.3934302014076827d0 !<Debye to atomic unit of dipole moment
        double precision, parameter :: Angstrom = 1.0d0/0.529177249d0 !<Ångstrom to bohr
      end module AtomicUnits
