#include "../defs.F90"

module m_fldsolver
  use m_globalnamespace
  use m_aux
  use m_domain
  use m_fields
  use m_helpers
  implicit none

  !--- PRIVATE functions -----------------------------------------!

  !...............................................................!
contains
  subroutine initRKstep()
    implicit none
    ! copy fields
    enx(:,:,:) = ex(:,:,:); eny(:,:,:) = ey(:,:,:); enz(:,:,:) = ez(:,:,:)
    bnx(:,:,:) = bx(:,:,:); bny(:,:,:) = by(:,:,:); bnz(:,:,:) = bz(:,:,:)
    ! init dE, dB
    dex(:,:,:) = -1000.0; dey(:,:,:) = -1000.0; dez(:,:,:) = -1000.0
    dbx(:,:,:) = -1000.0; dby(:,:,:) = -1000.0; dbz(:,:,:) = -1000.0
    ! init rho
    rho(:,:,:) = -1000.0
  end subroutine initRKstep

  subroutine rk3Step(rk_c1, rk_c2, rk_c3)
    implicit none
    real, intent(in) :: rk_c1, rk_c2, rk_c3
    call computeRho()
    call rk3Push()
    call rk3Update(rk_c1, rk_c2, rk_c3)
    ! call checkEB()
  end subroutine rk3Step

  subroutine computeRho()
    implicit none
    integer :: i1, i2, j1, j2, k1, k2
    integer :: i, j, k

    i1 = -1; i2 = this_meshblock%ptr%sx + 1
    j1 = -1; j2 = this_meshblock%ptr%sy + 1
    k1 = -1; k2 = this_meshblock%ptr%sz + 1

    do i = i1, i2
      do j = j1, j2
        do k = k1, k2
          rho(i, j, k) = (ex(i, j, k) - ex(i - 1, j, k)) +&
                       & (ey(i, j, k) - ey(i, j - 1, k)) +&
                       & (ez(i, j, k) - ez(i, j, k - 1))
        end do
      end do
    end do
  end subroutine computeRho

  subroutine rk3Push()
    implicit none
    integer :: i1, i2, j1, j2, k1, k2
    integer :: i, j, k
    real    :: jx, jy, jz
    real    :: corr
    real    :: intex, intey, intez, intbx, intby, intbz, bmag, emag, intrho

    i1 = -1; i2 = this_meshblock%ptr%sx
    j1 = -1; j2 = this_meshblock%ptr%sy
    k1 = -1; k2 = this_meshblock%ptr%sz

    do i = i1, i2
      do j = j1, j2
        do k = k1, k2
          dbx(i, j, k) = CC * ((ey(i, j, k + 1) - ey(i, j, k)) - (ez(i, j + 1, k) - ez(i, j, k)));
          dby(i, j, k) = CC * ((ez(i + 1, j, k) - ez(i, j, k)) - (ex(i, j, k + 1) - ex(i, j, k)));
          dbz(i, j, k) = CC * ((ex(i, j + 1, k) - ex(i, j, k)) - (ey(i + 1, j, k) - ey(i, j, k)));

          dex(i, j, k) = -CC * (((by(i, j, k) - by(i, j, k - 1)) - (bz(i, j, k) - bz(i, j - 1, k))) -&
                              & ((b0y(i, j, k) - b0y(i, j, k - 1)) - (b0z(i, j, k) - b0z(i, j - 1, k))));
          dey(i, j, k) = -CC * (((bz(i, j, k) - bz(i - 1, j, k)) - (bx(i, j, k) - bx(i, j, k - 1))) -&
                              & ((b0z(i, j, k) - b0z(i - 1, j, k)) - (b0x(i, j, k) - b0x(i, j, k - 1))));
          dez(i, j, k) = -CC * (((bx(i, j, k) - bx(i, j - 1, k)) - (by(i, j, k) - by(i - 1, j, k))) -&
                              & ((b0x(i, j, k) - b0x(i, j - 1, k)) - (b0y(i, j, k) - b0y(i - 1, j, k))));

          !--------X current-------------------------------------
          intrho = 0.5 * (rho(i + 1, j, k) + rho(i, j, k))

          ! call interpFromEdges(0.5, 0.0, 0.0, i, j, k, &
          !                    & ex, ey, ez, &
          !                    & intex, intey, intez)
          ! call interpFromFaces(0.5, 0.0, 0.0, i, j, k, &
          !                    & bx, by, bz, &
          !                    & intbx, intby, intbz)
          intex = ex(i, j, k)
          intey = 0.25 * (ey(i, j, k) + ey(i + 1, j, k) + ey(i, j - 1, k) + ey(i + 1, j - 1, k))
          intez = 0.25 * (ez(i, j, k) + ez(i + 1, j, k) + ez(i, j, k - 1) + ez(i + 1, j, k - 1))

          intbx = 0.125 * (bx(i, j, k) + bx(i, j - 1, k) + bx(i + 1, j - 1, k) + bx(i + 1, j, k) +&
                         & bx(i, j, k - 1) + bx(i, j - 1, k - 1) + bx(i + 1, j - 1, k - 1) + bx(i + 1, j, k - 1))
          intby = 0.5 * (by(i, j, k) + by(i, j, k - 1))
          intbz = 0.5 * (bz(i, j, k) + bz(i, j - 1, k))
          !---------------------------------------------------

          bmag = intbx**2 + intby**2 + intbz**2
          emag = intex**2 + intey**2 + intez**2
          if (emag .gt. bmag) then
            corr = sqrt(bmag / emag)
            intey = intey * corr
            intez = intez * corr
          endif

          jx = CC / (bmag + TINY) * (intrho * (intey * intbz - intby * intez))

          !--------Y current--------------------------------------
          intrho = 0.5 * (rho(i, j + 1, k) + rho(i, j, k))

          ! call interpFromEdges(0.0, 0.5, 0.0, i, j, k, &
          !                    & ex, ey, ez, &
          !                    & intex, intey, intez)
          ! call interpFromFaces(0.0, 0.5, 0.0, i, j, k, &
          !                    & bx, by, bz, &
          !                    & intbx, intby, intbz)
          intex = 0.25 * (ex(i, j, k) + ex(i - 1, j, k) + ex(i, j + 1, k) + ex(i - 1, j + 1, k))
          intey = ey(i, j, k)
          intez = 0.25 * (ez(i, j, k) + ez(i, j + 1, k) + ez(i, j, k - 1) + ez(i, j + 1, k - 1))

          intbx = 0.5 * (bx(i, j, k) + bx(i, j, k - 1))
          intby = 0.125 * (by(i, j, k) + by(i - 1, j, k) + by(i - 1, j + 1, k) + by(i, j + 1, k) +&
                         & by(i, j, k - 1) + by(i - 1, j, k - 1) + by(i - 1, j + 1, k - 1) + by(i, j + 1, k - 1))
          intbz = 0.5 * (bz(i, j, k) + bz(i - 1, j, k))
          !---------------------------------------------------

          bmag = intbx**2 + intby**2 + intbz**2
          emag = intex**2 + intey**2 + intez**2
          if (emag .gt. bmag) then
            corr = sqrt(bmag / emag)
            intex = intex * corr
            intez = intez * corr
          endif

          jy = CC / (bmag + TINY) * (intrho * (intez * intbx - intex * intbz))

          !--------Z current--------------------------------------
          intrho = 0.5 * (rho(i, j, k) + rho(i, j, k + 1))

          ! call interpFromEdges(0.0, 0.0, 0.5, i, j, k, &
          !                    & ex, ey, ez, &
          !                    & intex, intey, intez)
          ! call interpFromFaces(0.0, 0.0, 0.5, i, j, k, &
          !                    & bx, by, bz, &
          !                    & intbx, intby, intbz)
          intex = 0.25 * (ex(i, j, k) + ex(i - 1, j, k) + ex(i, j, k + 1) + ex(i - 1, j, k + 1))
          intey = 0.25 * (ey(i, j, k) + ey(i, j - 1, k) + ey(i, j, k + 1) + ey(i, j - 1, k + 1))
          intez = ez(i, j, k)

          intbx = 0.5 * (bx(i, j, k) + bx(i, j - 1, k))
          intby = 0.5 * (by(i, j, k) + by(i - 1, j, k))
          intbz = 0.125 * (bz(i, j, k) + bz(i - 1, j, k) + bz(i - 1, j - 1, k) + bz(i, j - 1, k) +&
                         & bz(i, j, k + 1) + bz(i - 1, j, k + 1) + bz(i - 1, j - 1, k + 1) + bz(i, j - 1, k + 1) )
          !---------------------------------------------------

          bmag = intbx**2 + intby**2 + intbz**2
          emag = intex**2 + intey**2 + intez**2
          if (emag .gt. bmag) then
            corr = sqrt(bmag / emag)
            intex = intex * corr
            intey = intey * corr
          endif

          jz = CC / (bmag + TINY) * (intrho * (intex * intby - intbx * intey))
          !---------------------------------------------------

          dex(i, j, k) = dex(i, j, k) - jx
          dey(i, j, k) = dey(i, j, k) - jy
          dez(i, j, k) = dez(i, j, k) - jz
        end do
      end do
    end do
  end subroutine rk3Push

  subroutine rk3Update(rk_c1, rk_c2, rk_c3)
    implicit none
    real, intent(in) :: rk_c1, rk_c2, rk_c3
    integer :: i, j, k
    integer :: i1, i2, j1, j2, k1, k2

    i1 = -1; i2 = this_meshblock%ptr%sx
    j1 = -1; j2 = this_meshblock%ptr%sy
    k1 = -1; k2 = this_meshblock%ptr%sz

    do i = i1, i2
      do j = j1, j2
        do k = k1, k2
          ! update E-field
          ex(i, j, k) = rk_c1 * enx(i, j, k) + rk_c2 * ex(i, j, k) + rk_c3 * dex(i, j, k)
          ey(i, j, k) = rk_c1 * eny(i, j, k) + rk_c2 * ey(i, j, k) + rk_c3 * dey(i, j, k)
          ez(i, j, k) = rk_c1 * enz(i, j, k) + rk_c2 * ez(i, j, k) + rk_c3 * dez(i, j, k)
          ! save E-field
          dex(i, j, k) = ex(i, j, k)
          dey(i, j, k) = ey(i, j, k)
          dez(i, j, k) = ez(i, j, k)
          ! update B-field
          bx(i, j, k) = rk_c1 * bnx(i, j, k) + rk_c2 * bx(i, j, k) + rk_c3 * dbx(i, j, k)
          by(i, j, k) = rk_c1 * bny(i, j, k) + rk_c2 * by(i, j, k) + rk_c3 * dby(i, j, k)
          bz(i, j, k) = rk_c1 * bnz(i, j, k) + rk_c2 * bz(i, j, k) + rk_c3 * dbz(i, j, k)
        end do
      end do
    end do
  end subroutine rk3Update

  subroutine cleanEpar()
    implicit none
    integer :: i, j, k
    real    :: jx, jy, jz
    real    :: intex, intey, intez, intbx, intby, intbz, bmag, epar
    integer :: i1, i2, j1, j2, k1, k2

    i1 = 0; i2 = this_meshblock%ptr%sx - 1
    j1 = 0; j2 = this_meshblock%ptr%sy - 1
    k1 = 0; k2 = this_meshblock%ptr%sz - 1

    do i = i1, i2
      do j = j1, j2
        do k = k1, k2
          ! -------interpolate E x-----------------------------
          intex = ex(i, j, k)
          intey = 0.25 * (ey(i, j, k) + ey(i + 1, j, k) + ey(i, j - 1, k) + ey(i + 1, j - 1, k))
          intez = 0.25 * (ez(i, j, k) + ez(i + 1, j, k) + ez(i, j, k - 1) + ez(i + 1, j, k - 1))
          !-------interpolate B x-----------------------------
          intbx = 0.125 * (bx(i, j, k) + bx(i, j - 1, k) + bx(i + 1, j - 1, k) + bx(i + 1, j, k) +&
                         & bx(i, j, k - 1) + bx(i, j - 1, k - 1) + bx(i + 1, j - 1, k - 1) + bx(i + 1, j, k - 1))
          intby = 0.5 * (by(i, j, k) + by(i, j, k - 1))
          intbz = 0.5 * (bz(i, j, k) + bz(i, j - 1, k))
          !---------------------------------------------------
          ! call interpFromEdges(0.5, 0.0, 0.0, i, j, k, &
          !                    & ex, ey, ez, &
          !                    & intex, intey, intez)
          ! call interpFromFaces(0.5, 0.0, 0.0, i, j, k, &
          !                    & bx, by, bz, &
          !                    & intbx, intby, intbz)
          bmag = intbx**2 + intby**2 + intbz**2
          epar = (intex * intbx + intey * intby + intez * intbz)
          dex(i, j, k) = ex(i, j, k) - epar * intbx / bmag

          !-------interpolate E y-----------------------------
          intex = 0.25 * (ex(i, j, k) + ex(i - 1, j, k) + ex(i, j + 1, k) + ex(i - 1, j + 1, k))
          intey = ey(i, j, k)
          intez = 0.25 * (ez(i, j, k) + ez(i, j + 1, k) + ez(i, j, k - 1) + ez(i, j + 1, k - 1))
          !-------interpolate B y-----------------------------
          intbx = 0.5 * (bx(i, j, k) + bx(i, j, k - 1))
          intby = 0.125 * (by(i, j, k) + by(i - 1, j, k) + by(i - 1, j + 1, k) + by(i, j + 1, k) +&
                         & by(i, j, k - 1) + by(i - 1, j, k - 1) + by(i - 1, j + 1, k - 1) + by(i, j + 1, k - 1))
          intbz = 0.5 * (bz(i, j, k) + bz(i - 1, j, k))
          !---------------------------------------------------
          ! call interpFromEdges(0.0, 0.5, 0.0, i, j, k, &
          !                    & ex, ey, ez, &
          !                    & intex, intey, intez)
          ! call interpFromFaces(0.0, 0.5, 0.0, i, j, k, &
          !                    & bx, by, bz, &
          !                    & intbx, intby, intbz)
          bmag = intbx**2 + intby**2 + intbz**2
          epar = (intex * intbx + intey * intby + intez * intbz)
          dey(i, j, k) = ey(i, j, k) - epar * intby / bmag

          !-------interpolate E z-----------------------------
          intex = 0.25 * (ex(i, j, k) + ex(i - 1, j, k) + ex(i, j, k + 1) + ex(i - 1, j, k + 1))
          intey = 0.25 * (ey(i, j, k) + ey(i, j - 1, k) + ey(i, j, k + 1) + ey(i, j - 1, k + 1))
          intez = ez(i, j, k)
          !-------interpolate B z-----------------------------
          intbx = 0.5 * (bx(i, j, k) + bx(i, j - 1, k))
          intby = 0.5 * (by(i, j, k) + by(i - 1, j, k))
          intbz = 0.125 * (bz(i, j, k) + bz(i - 1, j, k) + bz(i - 1, j - 1, k) + bz(i, j - 1, k) +&
                         & bz(i, j, k + 1) + bz(i - 1, j, k + 1) + bz(i - 1, j - 1, k + 1) + bz(i, j - 1, k + 1) )
          !---------------------------------------------------
          ! call interpFromEdges(0.0, 0.0, 0.5, i, j, k, &
          !                    & ex, ey, ez, &
          !                    & intex, intey, intez)
          ! call interpFromFaces(0.0, 0.0, 0.5, i, j, k, &
          !                    & bx, by, bz, &
          !                    & intbx, intby, intbz)
          bmag = intbx**2 + intby**2 + intbz**2
          epar = (intex * intbx + intey * intby + intez * intbz)
          dez(i, j, k) = ez(i, j, k) - epar * intbz / bmag
        end do
      end do
    end do
    ex(:,:,:) = dex(:,:,:); ey(:,:,:) = dey(:,:,:); ez(:,:,:) = dez(:,:,:)
  end subroutine cleanEpar

  subroutine checkEB()
    implicit none
    integer :: i, j, k
    real    :: intex, intey, intez, intbx, intby, intbz, bmag, emag, corr
    integer :: i1, i2, j1, j2, k1, k2

    i1 = -1; i2 = this_meshblock%ptr%sx + 1
    j1 = -1; j2 = this_meshblock%ptr%sy + 1
    k1 = -1; k2 = this_meshblock%ptr%sz + 1

    do i = i1, i2
      do j = j1, j2
        do k = k1, k2
          !-------interpolate on Ex---------------------------
          intex = ex(i, j, k)
          intey = 0.25 * (ey(i, j, k) + ey(i + 1, j, k) + ey(i, j - 1, k) + ey(i + 1, j - 1, k))
          intez = 0.25 * (ez(i, j, k) + ez(i + 1, j, k) + ez(i, j, k - 1) + ez(i + 1, j, k - 1))

          intbx = 0.125 * (bx(i, j, k) + bx(i, j - 1, k) + bx(i + 1, j - 1, k) + bx(i + 1, j, k) +&
                         & bx(i, j, k - 1) + bx(i, j - 1, k - 1) + bx(i + 1, j - 1, k - 1) + bx(i + 1, j, k - 1))
          intby = 0.5 * (by(i, j, k) + by(i, j, k - 1))
          intbz = 0.5 * (bz(i, j, k) + bz(i, j - 1, k))
          ! call interpFromEdges(0.5, 0.0, 0.0, i, j, k, &
          !                    & ex, ey, ez, &
          !                    & intex, intey, intez)
          ! call interpFromFaces(0.5, 0.0, 0.0, i, j, k, &
          !                    & bx, by, bz, &
          !                    & intbx, intby, intbz)

          bmag = intbx**2 + intby**2 + intbz**2
          emag = intex**2 + intey**2 + intez**2
          if (emag .gt. bmag) then
            corr = sqrt(bmag / emag)
            intey = intey * corr
            intez = intez * corr
          endif
          ex(i, j, k) = ex(i, j, k) * corr
          !---------------------------------------------------

          !-------interpolate on Ey---------------------------
          ! call interpFromEdges(0.0, 0.5, 0.0, i, j, k, &
          !                    & ex, ey, ez, &
          !                    & intex, intey, intez)
          ! call interpFromFaces(0.0, 0.5, 0.0, i, j, k, &
          !                    & bx, by, bz, &
          !                    & intbx, intby, intbz)
          intex = 0.25 * (ex(i, j, k) + ex(i - 1, j, k) + ex(i, j + 1, k) + ex(i - 1, j + 1, k))
          intey = ey(i, j, k)
          intez = 0.25 * (ez(i, j, k) + ez(i, j + 1, k) + ez(i, j, k - 1) + ez(i, j + 1, k - 1))

          intbx = 0.5 * (bx(i, j, k) + bx(i, j, k - 1))
          intby = 0.125 * (by(i, j, k) + by(i - 1, j, k) + by(i - 1, j + 1, k) + by(i, j + 1, k) +&
                         & by(i, j, k - 1) + by(i - 1, j, k - 1) + by(i - 1, j + 1, k - 1) + by(i, j + 1, k - 1))
          intbz = 0.5 * (bz(i, j, k) + bz(i - 1, j, k))

          bmag = intbx**2 + intby**2 + intbz**2
          emag = intex**2 + intey**2 + intez**2
          if (emag .gt. bmag) then
            corr = sqrt(bmag / emag)
          endif
          ey(i, j, k) = ey(i, j, k) * corr
          !---------------------------------------------------

          !-------interpolate on Ez---------------------------
          intex = 0.25 * (ex(i, j, k) + ex(i - 1, j, k) + ex(i, j, k + 1) + ex(i - 1, j, k + 1))
          intey = 0.25 * (ey(i, j, k) + ey(i, j - 1, k) + ey(i, j, k + 1) + ey(i, j - 1, k + 1))
          intez = ez(i, j, k)

          intbx = 0.5 * (bx(i, j, k) + bx(i, j - 1, k))
          intby = 0.5 * (by(i, j, k) + by(i - 1, j, k))
          intbz = 0.125 * (bz(i, j, k) + bz(i - 1, j, k) + bz(i - 1, j - 1, k) + bz(i, j - 1, k) +&
                         & bz(i, j, k + 1) + bz(i - 1, j, k + 1) + bz(i - 1, j - 1, k + 1) + bz(i, j - 1, k + 1))
          ! call interpFromEdges(0.0, 0.0, 0.5, i, j, k, &
          !                    & ex, ey, ez, &
          !                    & intex, intey, intez)
          ! call interpFromFaces(0.0, 0.0, 0.5, i, j, k, &
          !                    & bx, by, bz, &
          !                    & intbx, intby, intbz)

          bmag = intbx**2 + intby**2 + intbz**2
          emag = intex**2 + intey**2 + intez**2
          if (emag .gt. bmag) then
            corr = sqrt(bmag / emag)
          endif
          ez(i, j, k) = ez(i, j, k) * corr
          !---------------------------------------------------
        end do
      end do
    end do

  end subroutine checkEB

  end module m_fldsolver
