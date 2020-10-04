#include "../defs.F90"

module m_helpers
  use m_globalnamespace
  use m_domain
  use m_fields
  implicit none
contains
  function rnkToInd(rnk)
    implicit none
    integer, intent(in)   :: rnk
    integer, dimension(3) :: rnkToInd
    if ((rnk .lt. 0) .or. (rnk .ge. mpi_size)) then
      rnkToInd = (/-1, -1, -1/)
    else
      rnkToInd(3) = rnk / (sizex * sizey)
      rnkToInd(2) = (rnk - sizex * sizey * rnkToInd(3)) / sizex
      rnkToInd(1) = rnk - sizex * sizey * rnkToInd(3) - sizex * rnkToInd(2)
    end if
  end function rnkToInd

  function indToRnk(ind)
    implicit none
    integer, intent(in)                 :: ind(3)
    integer                             :: ind_t(3), indToRnk
    ind_t = ind
    if (boundary_Px .eq. 1) then
      if (ind_t(1) .lt. 0) ind_t(1) = ind_t(1) + sizex
      ind_t(1) = modulo(ind_t(1), sizex)
    end if
    if (boundary_Py .eq. 1) then
      if (ind_t(2) .lt. 0) ind_t(2) = ind_t(2) + sizey
      ind_t(2) = modulo(ind_t(2), sizey)
    end if
    if (boundary_Pz .eq. 1) then
      if (ind_t(3) .lt. 0) ind_t(3) = ind_t(3) + sizez
      ind_t(3) = modulo(ind_t(3), sizez)
    end if
    if ((ind_t(1) .lt. 0) .or. (ind_t(1) .ge. sizex) .or.&
      & (ind_t(2) .lt. 0) .or. (ind_t(2) .ge. sizey) .or.&
      & (ind_t(3) .lt. 0) .or. (ind_t(3) .ge. sizez)) then
      indToRnk = -1
    else
      indToRnk = ind_t(3) * sizex * sizey + ind_t(2) * sizex + ind_t(1)
    end if
  end function indToRnk

  subroutine reassignNeighborsForAll()
    implicit none
    integer   :: rnk
    integer   :: ind1, ind2, ind3
    do rnk = 0, mpi_size - 1
      do ind1 = -1, 1
        do ind2 = -1, 1
          do ind3 = -1, 1
            call assignNeighbor(rnk, (/ ind1, ind2, ind3/))
          end do
        end do
      end do
    end do
    call computeNumberOfNeighbors()
  end subroutine reassignNeighborsForAll

  subroutine assignNeighbor(rnk, inds1)
    implicit none
    integer, intent(in)       :: rnk, inds1(3)
    integer                   :: rnk2, inds0(3)
    inds0 = rnkToInd(rnk)
    rnk2 = indToRnk([inds0(1) + inds1(1), inds0(2) + inds1(2), inds0(3) + inds1(3)])
    if (rnk2 .eq. -1) then
      meshblocks(rnk + 1)%neighbor(inds1(1), inds1(2), inds1(3))%ptr => null()
    else
      meshblocks(rnk + 1)%neighbor(inds1(1), inds1(2), inds1(3))%ptr => meshblocks(rnk2 + 1)
    end if
  end subroutine assignNeighbor

  subroutine computeNumberOfNeighbors()
    integer   :: ind1, ind2, ind3
    integer   :: cntr
    cntr = 0
    do ind1 = -1, 1
      do ind2 = -1, 1
        do ind3 = -1, 1
          if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
          if (.not. associated(this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr)) cycle
          cntr = cntr + 1
        end do
      end do
    end do
    sendrecv_neighbors = cntr
  end subroutine computeNumberOfNeighbors

  subroutine interpFromEdges(dx, dy, dz, i, j, k, &
                           & fx, fy, fz, &
                           & intfx, intfy, intfz)
    !$omp declare simd(interpFromEdges)
    implicit none
    integer, intent(in)           :: i, j, k
    real, intent(in)              :: dx, dy, dz
    real, intent(in)              :: fx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST)
    real, intent(in)              :: fy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST)
    real, intent(in)              :: fz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST)
    real, intent(out)             :: intfx, intfy, intfz
    real                          :: c000, c100, c001, c101, c010, c110, c011, c111,&
                                   & c00, c01, c10, c11, c0, c1
    ! f_x
    c000 = 0.5 * (fx(     i,      j,      k) + fx(  i - 1,      j,      k))
    c100 = 0.5 * (fx(     i,      j,      k) + fx(  i + 1,      j,      k))
    c010 = 0.5 * (fx(     i,  j + 1,      k) + fx(  i - 1,  j + 1,      k))
    c110 = 0.5 * (fx(     i,  j + 1,      k) + fx(  i + 1,  j + 1,      k))
    c00 = c000 * (1 - dx) + c100 * dx
    c10 = c010 * (1 - dx) + c110 * dx
    c0 = c00 * (1 - dy) + c10 * dy
    c001 = 0.5 * (fx(     i,      j,  k + 1) + fx(  i - 1,      j,  k + 1))
    c101 = 0.5 * (fx(     i,      j,  k + 1) + fx(  i + 1,      j,  k + 1))
    c011 = 0.5 * (fx(     i,  j + 1,  k + 1) + fx(  i - 1,  j + 1,  k + 1))
    c111 = 0.5 * (fx(     i,  j + 1,  k + 1) + fx(  i + 1,  j + 1,  k + 1))
    c01 = c001 * (1 - dx) + c101 * dx
    c11 = c011 * (1 - dx) + c111 * dx
    c1 = c01 * (1 - dy) + c11 * dy
    intfx = c0 * (1 - dz) + c1 * dz

    ! f_y
    c000 = 0.5 * (fy(     i,      j,      k) + fy(      i,  j - 1,      k))
    c100 = 0.5 * (fy( i + 1,      j,      k) + fy(  i + 1,  j - 1,      k))
    c010 = 0.5 * (fy(     i,      j,      k) + fy(      i,  j + 1,      k))
    c110 = 0.5 * (fy( i + 1,      j,      k) + fy(  i + 1,  j + 1,      k))
    c00 = c000 * (1 - dx) + c100 * dx
    c10 = c010 * (1 - dx) + c110 * dx
    c0 = c00 * (1 - dy) + c10 * dy
    c001 = 0.5 * (fy(     i,      j,  k + 1) + fy(      i,  j - 1,  k + 1))
    c101 = 0.5 * (fy( i + 1,      j,  k + 1) + fy(  i + 1,  j - 1,  k + 1))
    c011 = 0.5 * (fy(     i,      j,  k + 1) + fy(      i,  j + 1,  k + 1))
    c111 = 0.5 * (fy( i + 1,      j,  k + 1) + fy(  i + 1,  j + 1,  k + 1))
    c01 = c001 * (1 - dx) + c101 * dx
    c11 = c011 * (1 - dx) + c111 * dx
    c1 = c01 * (1 - dy) + c11 * dy
    intfy = c0 * (1 - dz) + c1 * dz

    ! f_z
    c000 = 0.5 * (fz(     i,      j,      k) + fz(      i,      j,  k - 1))
    c100 = 0.5 * (fz( i + 1,      j,      k) + fz(  i + 1,      j,  k - 1))
    c010 = 0.5 * (fz(     i,  j + 1,      k) + fz(      i,  j + 1,  k - 1))
    c110 = 0.5 * (fz( i + 1,  j + 1,      k) + fz(  i + 1,  j + 1,  k - 1))
    c001 = 0.5 * (fz(     i,      j,      k) + fz(      i,      j,  k + 1))
    c101 = 0.5 * (fz( i + 1,      j,      k) + fz(  i + 1,      j,  k + 1))
    c011 = 0.5 * (fz(     i,  j + 1,      k) + fz(      i,  j + 1,  k + 1))
    c111 = 0.5 * (fz( i + 1,  j + 1,      k) + fz(  i + 1,  j + 1,  k + 1))
    c00 = c000 * (1 - dx) + c100 * dx
    c01 = c001 * (1 - dx) + c101 * dx
    c10 = c010 * (1 - dx) + c110 * dx
    c11 = c011 * (1 - dx) + c111 * dx
    c0 = c00 * (1 - dy) + c10 * dy
    c1 = c01 * (1 - dy) + c11 * dy
    intfz = c0 * (1 - dz) + c1 * dz
  end subroutine interpFromEdges

  subroutine interpFromFaces(dx, dy, dz, i, j, k, &
                           & fx, fy, fz, &
                           & intfx, intfy, intfz)
    !$omp declare simd(interpFromFaces)
    implicit none
    integer, intent(in)           :: i, j, k
    real, intent(in)              :: dx, dy, dz
    real, intent(in)              :: fx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST)
    real, intent(in)              :: fy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST)
    real, intent(in)              :: fz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST)
    real, intent(out)             :: intfx, intfy, intfz
    real                          :: c000, c100, c001, c101, c010, c110, c011, c111,&
                                   & c00, c01, c10, c11, c0, c1
    ! f_x
    c000 = 0.25 * (fx(      i,      j,      k) + fx(      i,  j - 1,      k) +&
                 & fx(      i,      j,  k - 1) + fx(      i,  j - 1,  k - 1))
    c100 = 0.25 * (fx(  i + 1,      j,      k) + fx(  i + 1,  j - 1,      k) +&
                 & fx(  i + 1,      j,  k - 1) + fx(  i + 1,  j - 1,  k - 1))
    c001 = 0.25 * (fx(      i,      j,      k) + fx(      i,      j,  k + 1) +&
                 & fx(      i,  j - 1,      k) + fx(      i,  j - 1,  k + 1))
    c101 = 0.25 * (fx(  i + 1,      j,      k) + fx(  i + 1,      j,  k + 1) +&
                 & fx(  i + 1,  j - 1,      k) + fx(  i + 1,  j - 1,  k + 1))
    c010 = 0.25 * (fx(      i,      j,      k) + fx(      i,  j + 1,      k) +&
                 & fx(      i,      j,  k - 1) + fx(      i,  j + 1,  k - 1))
    c110 = 0.25 * (fx(  i + 1,      j,      k) + fx(  i + 1,      j,  k - 1) +&
                 & fx(  i + 1,  j + 1,  k - 1) + fx(  i + 1,  j + 1,      k))
    c011 = 0.25 * (fx(      i,      j,      k) + fx(      i,  j + 1,      k) +&
                 & fx(      i,  j + 1,  k + 1) + fx(      i,      j,  k + 1))
    c111 = 0.25 * (fx(  i + 1,      j,      k) + fx(  i + 1,  j + 1,      k) +&
                 & fx(  i + 1,  j + 1,  k + 1) + fx(  i + 1,      j,  k + 1))
    c00 = c000 * (1 - dx) + c100 * dx
    c01 = c001 * (1 - dx) + c101 * dx
    c10 = c010 * (1 - dx) + c110 * dx
    c11 = c011 * (1 - dx) + c111 * dx
    c0 = c00 * (1 - dy) + c10 * dy
    c1 = c01 * (1 - dy) + c11 * dy
    intfx = c0 * (1 - dz) + c1 * dz

    ! f_y
    c000 = 0.25 * (fy(  i - 1,      j,  k - 1) + fy(  i - 1,      j,      k) +&
                 & fy(      i,      j,  k - 1) + fy(      i,      j,      k))
    c100 = 0.25 * (fy(      i,      j,  k - 1) + fy(      i,      j,      k) +&
                 & fy(  i + 1,      j,  k - 1) + fy(  i + 1,      j,      k))
    c001 = 0.25 * (fy(  i - 1,      j,      k) + fy(  i - 1,      j,  k + 1) +&
                 & fy(      i,      j,      k) + fy(      i,      j,  k + 1))
    c101 = 0.25 * (fy(      i,      j,      k) + fy(      i,      j,  k + 1) +&
                 & fy(  i + 1,      j,      k) + fy(  i + 1,      j,  k + 1))
    c010 = 0.25 * (fy(  i - 1,  j + 1,  k - 1) + fy(  i - 1,  j + 1,      k) +&
                 & fy(      i,  j + 1,  k - 1) + fy(      i,  j + 1,      k))
    c110 = 0.25 * (fy(      i,  j + 1,  k - 1) + fy(      i,  j + 1,      k) +&
                 & fy(  i + 1,  j + 1,  k - 1) + fy(  i + 1,  j + 1,      k))
    c011 = 0.25 * (fy(  i - 1,  j + 1,      k) + fy(  i - 1,  j + 1,  k + 1) +&
                 & fy(      i,  j + 1,      k) + fy(      i,  j + 1,  k + 1))
    c111 = 0.25 * (fy(      i,  j + 1,      k) + fy(      i,  j + 1,  k + 1) +&
                 & fy(  i + 1,  j + 1,      k) + fy(  i + 1,  j + 1,  k + 1))
    c00 = c000 * (1 - dx) + c100 * dx
    c01 = c001 * (1 - dx) + c101 * dx
    c10 = c010 * (1 - dx) + c110 * dx
    c11 = c011 * (1 - dx) + c111 * dx
    c0 = c00 * (1 - dy) + c10 * dy
    c1 = c01 * (1 - dy) + c11 * dy
    intfy = c0 * (1 - dz) + c1 * dz

    ! f_z
    c000 = 0.25 * (fz(  i - 1,  j - 1,      k) + fz(  i - 1,      j,      k) +&
                 & fz(      i,  j - 1,      k) + fz(      i,      j,      k))
    c100 = 0.25 * (fz(      i,  j - 1,      k) + fz(      i,      j,      k) +&
                 & fz(  i + 1,  j - 1,      k) + fz(  i + 1,      j,      k))
    c001 = 0.25 * (fz(  i - 1,  j - 1,  k + 1) + fz(  i - 1,      j,  k + 1) +&
                 & fz(      i,  j - 1,  k + 1) + fz(      i,      j,  k + 1))
    c101 = 0.25 * (fz(      i,  j - 1,  k + 1) + fz(      i,      j,  k + 1) +&
                 & fz(  i + 1,  j - 1,  k + 1) + fz(  i + 1,      j,  k + 1))
    c010 = 0.25 * (fz(  i - 1,      j,      k) + fz(  i - 1,  j + 1,      k) +&
                 & fz(      i,      j,      k) + fz(      i,  j + 1,      k))
    c110 = 0.25 * (fz(      i,      j,      k) + fz(      i,  j + 1,      k) +&
                 & fz(  i + 1,      j,      k) + fz(  i + 1,  j + 1,      k))
    c011 = 0.25 * (fz(  i - 1,      j,  k + 1) + fz(  i - 1,  j + 1,  k + 1) +&
                 & fz(      i,      j,  k + 1) + fz(      i,  j + 1,  k + 1))
    c111 = 0.25 * (fz(      i,      j,  k + 1) + fz(      i,  j + 1,  k + 1) +&
                 & fz(  i + 1,      j,  k + 1) + fz(  i + 1,  j + 1,  k + 1))
    c00 = c000 * (1 - dx) + c100 * dx
    c01 = c001 * (1 - dx) + c101 * dx
    c10 = c010 * (1 - dx) + c110 * dx
    c11 = c011 * (1 - dx) + c111 * dx
    c0 = c00 * (1 - dy) + c10 * dy
    c1 = c01 * (1 - dy) + c11 * dy
    intfz = c0 * (1 - dz) + c1 * dz
  end subroutine interpFromFaces

  subroutine computeForceFreeCurrent(jx0, jy0, jz0, i, j, k)
    implicit none
    real, intent(out)   :: jx0, jy0, jz0
    integer, intent(in) :: i, j, k
    real                :: divE, ex0, ey0, ez0, bx0, by0, bz0, bb2
    real                :: dx1, dx2, dy1, dy2, dz1, dz2, curls
    real                :: curlBx, curlBy, curlBz
    real                :: curlEx, curlEy, curlEz
    real                :: EcrossBx, EcrossBy, EcrossBz
    ! `div E` on nodes
    divE = (ex(i, j, k) - ex(i - 1, j, k)) +&
         & (ey(i, j, k) - ey(i, j - 1, k)) +&
         & (ez(i, j, k) - ez(i, j, k - 1))

    ! `E` and `B` on nodes
    call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
    call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
    bb2 = bx0**2 + by0**2 + bz0**2

    ! `curl B` on nodes
    dx1 = (bz(    i,    j,    k) - bz(    i,j - 1,    k)) - (by(    i,    j,    k) - by(    i,    j,k - 1))
    dx2 = (bz(i - 1,    j,    k) - bz(i - 1,j - 1,    k)) - (by(i - 1,    j,    k) - by(i - 1,    j,k - 1))
    curlBx = 0.5 * (dx1 + dx2)
    dy1 = (bx(    i,    j,    k) - bx(    i,    j,k - 1)) - (bz(    i,    j,    k) - bz(i - 1,    j,    k))
    dy2 = (bx(    i,j - 1,    k) - bx(    i,j - 1,k - 1)) - (bz(    i,j - 1,    k) - bz(i - 1,j - 1,    k))
    curlBy = 0.5 * (dy1 + dy2)
    dz1 = (by(    i,    j,    k) - by(i - 1,    j,    k)) - (bx(    i,    j,    k) - bx(    i,j - 1,    k))
    dz2 = (by(    i,    j,k - 1) - by(i - 1,    j,k - 1)) - (bx(    i,    j,k - 1) - bx(    i,j - 1,k - 1))
    curlBz = 0.5 * (dz1 + dz2)

    ! `curl E` on nodes
    curlEx = 0.5 * (0.5 * (ez(i, j + 1, k) - ez(i, j - 1, k)) + 0.5 * (ez(i, j + 1, k - 1) - ez(i, j - 1, k - 1))) -&
           & 0.5 * (0.5 * (ey(i, j, k + 1) - ey(i, j, k - 1)) + 0.5 * (ey(i, j - 1, k + 1) - ey(i, j - 1, k - 1)))
    curlEy = 0.5 * (0.5 * (ex(i, j, k + 1) - ex(i, j, k - 1)) + 0.5 * (ex(i - 1, j, k + 1) - ex(i - 1, j, k - 1))) -&
           & 0.5 * (0.5 * (ez(i + 1, j, k) - ez(i - 1, j, k)) + 0.5 * (ez(i + 1, j, k - 1) - ez(i - 1, j, k - 1)))
    curlEz = 0.5 * (0.5 * (ey(i + 1, j, k) - ey(i - 1, j, k)) + 0.5 * (ey(i + 1, j - 1, k) - ey(i - 1, j - 1, k))) -&
           & 0.5 * (0.5 * (ex(i, j + 1, k) - ex(i, j - 1, k)) + 0.5 * (ex(i - 1, j + 1, k) - ex(i - 1, j - 1, k)))

    curls = ((bx0 * curlBx + by0 * curlBy + bz0 * curlBz) - (ex0 * curlEx + ey0 * curlEy + ez0 * curlEz))
    EcrossBx = (ey0 * bz0 - ez0 * by0)
    EcrossBy = -(ex0 * bz0 - ez0 * bx0)
    EcrossBz = (ex0 * by0 - ey0 * bx0)

    jx0 = (divE * EcrossBx - bx0 * curls) / bb2
    jy0 = (divE * EcrossBy - by0 * curls) / bb2
    jz0 = (divE * EcrossBz - bz0 * curls) / bb2
  end subroutine computeForceFreeCurrent
end module m_helpers
