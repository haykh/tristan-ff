#include "../defs.F90"

!--- FIELDS ----------------------------------------------------!
! To store all the field related quantities
!   - field indices go
!       from `-NGHOST`
!       to `sx/sy/sz - 1 + NGHOST` inclusively
!...............................................................!

module m_fields
  use m_globalnamespace
  implicit none

  real, allocatable   :: ex(:,:,:), ey(:,:,:), ez(:,:,:),&
                       & bx(:,:,:), by(:,:,:), bz(:,:,:)
  real, allocatable   :: dex(:,:,:), dey(:,:,:), dez(:,:,:),&
                       & dbx(:,:,:), dby(:,:,:), dbz(:,:,:)
  real, allocatable   :: b0x(:,:,:), b0y(:,:,:), b0z(:,:,:),&
                       & enx(:,:,:), eny(:,:,:), enz(:,:,:),&
                       & bnx(:,:,:), bny(:,:,:), bnz(:,:,:)
  real, allocatable   :: rho(:,:,:)
  real, allocatable   :: recv_fld(:), send_fld(:)
  integer             :: sendrecv_buffsz, sendrecv_offsetsz
end module m_fields
