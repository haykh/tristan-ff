#include "../src/defs.F90"

module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_fields
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real    :: L0, r0, beta_twist, xc, yc, zc, bnorm, t0
  private :: L0, r0, beta_twist, xc, yc, zc, bnorm, t0
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userInitFields
  !...............................................................!
contains
  subroutine userInitialize()
    implicit none
    call userReadInput()
    call userInitFields()
  end subroutine userInitialize

  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'L0', L0)
    call getInput('problem', 'r_L0', r0)
    call getInput('problem', 't0', t0)
    call getInput('problem', 'beta_twist', beta_twist)

    xc = 0.5 * REAL(global_mesh%sx)
    yc = 0.5 * REAL(global_mesh%sy)
    zc = 0.0

    r0 = r0 * L0

    bnorm = 1.0
  end subroutine userReadInput

  subroutine userInitFields()
    implicit none
    integer :: i, j, k, i1, i2, j1, j2, k1, k2
    real    :: x_glob, y_glob, z_glob, rx, ry, rr

    ex(:,:,:) = 0; ey(:,:,:) = 0; ez(:,:,:) = 0
    bx(:,:,:) = 0; by(:,:,:) = 0

    i1 = 0; i2 = this_meshblock%ptr%sx - 1
    j1 = 0; j2 = this_meshblock%ptr%sy - 1
    k1 = 0; k2 = this_meshblock%ptr%sz - 1

    do i = i1 - NGHOST, i2 + NGHOST
      x_glob = REAL(i + this_meshblock%ptr%x0)
      do j = j1 - NGHOST, j2 + NGHOST
        y_glob = REAL(j + this_meshblock%ptr%y0)
        do k = k1 - NGHOST, k2 + NGHOST
          rx = (x_glob + 0.5 - xc)
          ry = (y_glob + 0.5 - yc)
          rr = sqrt(rx**2 + ry**2 + TINY)
          bz(i, j, k) = bnorm / (1.0 + (rr / r0))
        end do
      end do
    end do

    b0x(:,:,:) = bx(:,:,:); b0y(:,:,:) = by(:,:,:); b0z(:,:,:) = bz(:,:,:)
  end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!

  !............................................................!

  !--- boundaries ---------------------------------------------!
  subroutine userFieldBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
    integer :: i, j, k
    real    :: x_glob, y_glob, z_glob
    integer :: i1, i2, j1, j2, k1, k2
    real    :: rx, ry, rr, omega, ee0
    real    :: time
    time = REAL(step)

    i1 = 0; i2 = this_meshblock%ptr%sx - 1
    j1 = 0; j2 = this_meshblock%ptr%sy - 1
    k1 = 0; k2 = this_meshblock%ptr%sz - 1

    do i = i1 - NGHOST, i2 + NGHOST
      x_glob = REAL(i + this_meshblock%ptr%x0)
      do j = j1 - NGHOST, j2 + NGHOST
        y_glob = REAL(j + this_meshblock%ptr%y0)
        do k = k1 - NGHOST, k2 + NGHOST
          z_glob = REAL(k + this_meshblock%ptr%z0)

          ! if (z_glob .le. 0) then
          !   ex(i, j, k) = bnorm
          !   ey(i, j, k) = 0.0
          ! end if

          ! if ((z_glob .gt. 100) .and. (z_glob .lt. 119)) then
          !   by(i, j, k) = 0.0
          !   bx(i, j, k) = 0.0
          ! end if

          ! if (rr .lt. r0) then
          !   omega = beta_twist * CC / r0
          ! else
          !   omega = 0.0
          ! end if

          if (z_glob .le. 2) then
            rx = (x_glob + 0.5 - xc)
            ry = (y_glob - yc)
            rr = sqrt(rx**2 + ry**2 + TINY)
            omega = exp(-rr / r0) * beta_twist * CC / r0
            ee0 = -omega * bnorm * tanh(time / t0)
            ex(i, j, k) = rx * ee0

            rx = (x_glob - xc)
            ry = (y_glob + 0.5 - yc)
            rr = sqrt(rx**2 + ry**2 + TINY)
            omega = exp(-rr / r0) * beta_twist * CC / r0
            ee0 = -omega * bnorm * tanh(time / t0)
            ey(i, j, k) = ry * ee0

            rx = (x_glob + 0.5 - xc)
            ry = (y_glob + 0.5 - yc)
            rr = sqrt(rx**2 + ry**2 + TINY)
            bz(i, j, k) = bnorm / (1.0 + (rr / r0))
          end if

          if (z_glob .ge. REAL(global_mesh%sz) - 2) then
            rx = (x_glob + 0.5 - xc)
            ry = (y_glob + 0.5 - yc)
            rr = sqrt(rx**2 + ry**2 + TINY)
            bz(i, j, k) = bnorm / (1.0 + (rr / r0))
            ex(i, j, k) = 0.0
            ey(i, j, k) = 0.0
          end if
        end do
      end do
    end do

  end subroutine userFieldBoundaryConditions
  !............................................................!

  !--- custom output ------------------------------------------!
  subroutine userOutput(var, temp_, i, j, k)
    implicit none
    integer, intent(in)         :: var
    integer(kind=2), intent(in) :: i, j, k
    real, intent(out)           :: temp_
    real                :: dummy1_, dummy2_, dummy3_
    if (var .eq. 1) then
      temp_ = rho(i, j, k)
    else if (var .eq. 2) then
      ! ExB / (r |B|^2)
      temp_ = sqrt((-bz(i, j, k) * ex(i, j, k) + bx(i, j, k) * ez(i, j, k))**2 +&
                  & (bz(i, j, k) * ey(i, j, k) - by(i, j, k) * ez(i, j, k))**2)
      dummy1_ = REAL(this_meshblock%ptr%x0 + i, 4)
      dummy1_ = (dummy1_ - 0.5 * REAL(global_mesh%sx))
      dummy2_ = REAL(this_meshblock%ptr%y0 + j, 4)
      dummy2_ = (dummy2_ - 0.5 * REAL(global_mesh%sy))
      dummy3_ = sqrt(dummy1_**2 + dummy2_**2 + TINY)
      temp_ = temp_ / dummy3_
      dummy3_ = bx(i, j, k)**2 + by(i, j, k)**2 + bz(i, j, k)**2
      temp_ = temp_ / dummy3_
    else
      temp_ = 0.0
    end if
  end subroutine userOutput
  !............................................................!
end module m_userfile
