#include "../src/defs.F90"

module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_fields
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real :: a0, L0, h0, r0, t0, beta_twist, xc, yc1, yc2, zc, bnorm
  private :: a0, L0, h0, r0, t0, beta_twist, xc, yc1, yc2, zc, bnorm
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
    call getInput('problem', 'a_L0', a0)
    call getInput('problem', 'h_L0', h0)
    call getInput('problem', 'r_L0', r0)
    call getInput('problem', 't0', t0)
    call getInput('problem', 'beta_twist', beta_twist)

    xc = 0.5 * REAL(global_mesh % sx)
    yc1 = 0.5 * REAL(global_mesh % sy) - a0 * L0
    yc2 = 0.5 * REAL(global_mesh % sy) + a0 * L0
    zc = -h0 * L0

    bnorm = 1000.0 * h0**2 * L0**2
  end subroutine userReadInput

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    real :: x_glob, y_glob, z_glob
    real :: rx, ry1, ry2, rz, rr1, rr2
    integer :: i1, i2, j1, j2, k1, k2

    i1 = -NGHOST; i2 = this_meshblock % ptr % sx - 1 + NGHOST
    j1 = -NGHOST; j2 = this_meshblock % ptr % sy - 1 + NGHOST
    k1 = -NGHOST; k2 = this_meshblock % ptr % sz - 1 + NGHOST

    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = 0; by(:, :, :) = 0; bz(:, :, :) = 0

    do i = i1, i2
      x_glob = REAL(i + this_meshblock % ptr % x0)
      do j = j1, j2
        y_glob = REAL(j + this_meshblock % ptr % y0)
        do k = k1, k2
          z_glob = REAL(k + this_meshblock % ptr % z0)

          rx = (x_glob - xc)
          ry1 = (y_glob + 0.5 - yc1)
          ry2 = (y_glob + 0.5 - yc2)
          rz = (z_glob + 0.5 - zc)
          rr1 = sqrt(rx * rx + ry1 * ry1 + rz * rz)
          rr2 = sqrt(rx * rx + ry2 * ry2 + rz * rz)
          bx(i, j, k) = bnorm * (rx / rr1**3 - rx / rr2**3)

          rx = (x_glob + 0.5 - xc)
          ry1 = (y_glob - yc1)
          ry2 = (y_glob - yc2)
          rz = (z_glob + 0.5 - zc)
          rr1 = sqrt(rx * rx + ry1 * ry1 + rz * rz)
          rr2 = sqrt(rx * rx + ry2 * ry2 + rz * rz)
          by(i, j, k) = bnorm * (ry1 / rr1**3 - ry2 / rr2**3)

          rx = (x_glob + 0.5 - xc)
          ry1 = (y_glob + 0.5 - yc1)
          ry2 = (y_glob + 0.5 - yc2)
          rz = (z_glob - zc)
          rr1 = sqrt(rx * rx + ry1 * ry1 + rz * rz)
          rr2 = sqrt(rx * rx + ry2 * ry2 + rz * rz)
          bz(i, j, k) = bnorm * (rz / rr1**3 - rz / rr2**3)
        end do
      end do
    end do

    b0x(:, :, :) = bx(:, :, :); b0y(:, :, :) = by(:, :, :); b0z(:, :, :) = bz(:, :, :)
  end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!

  !............................................................!

  !--- boundaries ---------------------------------------------!
  subroutine userFieldBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
    integer :: i, j, k
    real :: x_glob, y_glob, z_glob
    integer :: i1, i2, j1, j2, k1, k2
    real :: rx, ry1, ry2, rz, rr1, rr2, omega, ee0, bz_loc
    real :: time
    time = REAL(step)

    i1 = -NGHOST; i2 = this_meshblock % ptr % sx - 1 + NGHOST
    j1 = -NGHOST; j2 = this_meshblock % ptr % sy - 1 + NGHOST
    k1 = -NGHOST; k2 = this_meshblock % ptr % sz - 1 + NGHOST

    do i = i1, i2
      x_glob = REAL(i + this_meshblock % ptr % x0)
      do j = j1, j2
        y_glob = REAL(j + this_meshblock % ptr % y0)
        do k = k1, k2
          z_glob = REAL(k + this_meshblock % ptr % z0)
          if (z_glob .lt. 0) then
            ! rx = (x_glob + 0.5 - xc)
            ! ry1 = (y_glob - yc1)
            ! ry2 = (y_glob - yc2)
            ! rz = (z_glob - zc)
            ! rr1 = sqrt(rx * rx + ry1 * ry1 + rz * rz)
            ! rr2 = sqrt(rx * rx + ry2 * ry2 + rz * rz)
            ! bz_loc = bnorm * (rz / rr1**3 - rz / rr2**3)
            ! omega = exp(-rr1 / r0) * beta_twist * CC / r0
            ! ee0 = -omega * bz_loc! * tanh(time / t0)
            ! ex(i, j, k) = rx * ee0
            !
            ! rx = (x_glob - xc)
            ! ry1 = (y_glob + 0.5 - yc1)
            ! ry2 = (y_glob + 0.5 - yc2)
            ! rz = (z_glob - zc)
            ! rr1 = sqrt(rx * rx + ry1 * ry1 + rz * rz)
            ! rr2 = sqrt(rx * rx + ry2 * ry2 + rz * rz)
            ! bz_loc = bnorm * (rz / rr1**3 - rz / rr2**3)
            ! omega = exp(-rr1 / r0) * beta_twist * CC / r0
            ! ee0 = -omega * bz_loc! * tanh(time / t0)
            ! ey(i, j, k) = ry1 * ee0
            !
            ! rx = (x_glob + 0.5 - xc)
            ! ry1 = (y_glob + 0.5 - yc1)
            ! ry2 = (y_glob + 0.5 - yc2)
            ! rz = (z_glob - zc)
            ! rr1 = sqrt(rx * rx + ry1 * ry1 + rz * rz)
            ! rr2 = sqrt(rx * rx + ry2 * ry2 + rz * rz)
            ! bz(i, j, k) = bnorm * (rz / rr1**3 - rz / rr2**3)
            ex(i, j, k) = ex(i, j, 0)
            ey(i, j, k) = ey(i, j, 0)
            bz(i, j, k) = bz(i, j, 0)
          end if

          if (z_glob .ge. REAL(global_mesh % sz)) then
            ex(i, j, k) = ex(i, j, this_meshblock % ptr % sz - 1)
            ey(i, j, k) = ey(i, j, this_meshblock % ptr % sz - 1)
            bz(i, j, k) = bz(i, j, this_meshblock % ptr % sz - 1)
          end if
        end do
      end do
    end do
  end subroutine userFieldBoundaryConditions
  !............................................................!

  !--- custom output ------------------------------------------!
  subroutine userOutput(var, temp_, i, j, k)
    implicit none
    integer, intent(in) :: var, i, j, k
    real, intent(out) :: temp_
    if (var .eq. 1) then
      temp_ = rho(i, j, k)
    else if (var .eq. 2) then
      temp_ = 0.0
    else
      temp_ = 0.0
    end if
  end subroutine userOutput
  !............................................................!
end module m_userfile
