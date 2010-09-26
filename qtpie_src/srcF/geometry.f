!>
!! Computes pairwise distances from Cartesian coordinates
!!
!! \param Point1, Point2: 3-vectors of double precisions
!! \return Cartesian distance in atomic units
!<
      double precision function Distance(Point1, Point2)
        implicit none
        double precision, dimension(3), intent(in) :: Point1, Point2
        double precision :: x, y, z
        x = Point2(1) - Point1(1)
        y = Point2(2) - Point1(2)
        z = Point2(3) - Point1(3)
        Distance = sqrt(x*x + y*y + z*z)
      end function Distance
