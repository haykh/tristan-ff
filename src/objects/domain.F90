#include "../defs.F90"

!--- DOMAIN ----------------------------------------------------!
! To store all domain related variables and constructions
!   - all grid meshblocks are stored as `mesh` type
!...............................................................!

module m_domain
  use m_globalnamespace
  implicit none

  type :: meshptr
    type(mesh), pointer :: ptr
  end type meshptr

  type :: mesh
    integer :: rnk          ! rank of the cpu that takes care of the current meshblock
    integer :: x0, y0, z0   ! global coordinates (in cells) of the corner
    integer :: sx, sy, sz   ! # of cells in each dimension
    ! pointers to the neighboring meshblocks
    type(meshptr), dimension(-1:1, -1:1, -1:1) :: neighbor
  end type mesh

  type(meshptr) :: this_meshblock  ! pointer to current (rank) meshblock
  type(mesh) :: global_mesh     ! global mesh parameters
  type(mesh), allocatable, target :: meshblocks(:)   ! meshblocks for all cpus

  ! boundary conditions for all dimensions
  !   - boundary = 1: periodic
  !   - boundary = 0: open
  integer :: boundary_x, boundary_y, boundary_z
  integer :: absorb_buff
  integer :: sendrecv_neighbors
end module m_domain
