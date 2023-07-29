#include "../defs.F90"

module m_aux
  use m_globalnamespace
  implicit none
  real(dprec) :: dseed

  abstract interface
    function spatialDistribution(x_glob, y_glob, z_glob,&
                               & dummy1, dummy2, dummy3)
      real :: spatialDistribution
      real, intent(in), optional :: x_glob, y_glob, z_glob
      real, intent(in), optional :: dummy1, dummy2, dummy3
    end function spatialDistribution
  end interface

  interface STR
    module procedure intToStr
    module procedure realToStr
  end interface STR

  type :: generic_var
    integer :: value_int
    real :: value_real
    logical :: value_bool
  end type generic_var

  type generic_string
    character(len=STR_MAX), allocatable :: str
  end type generic_string

  type :: simulation_params
    integer :: count
    integer, allocatable :: param_type(:) ! 1 = int, 2 = float, 3 = bool
    type(generic_string), allocatable :: param_group(:)
    type(generic_string), allocatable :: param_name(:)
    type(generic_var), allocatable :: param_value(:)
  end type simulation_params

  type(simulation_params) :: sim_params

  !--- PRIVATE functions -----------------------------------------!
  private :: intToStr, realToStr
  !...............................................................!
contains
  subroutine printDiag(bool, msg, prepend)
    implicit none
    character(len=*), intent(in) :: msg
    logical, intent(in) :: bool
    logical, optional, intent(in) :: prepend
    character(len=STR_MAX) :: dummy
    integer :: sz, i, ierr
#ifdef DEBUG
    if (bool) then
      sz = len(trim(msg))
      if (present(prepend)) then
        if (prepend) then
          dummy = '...'
          sz = sz + 3
        end if
      else
        dummy = ''
      end if
      dummy = trim(dummy)//trim(msg)
      do i = 1, (38 - sz)
        dummy = trim(dummy)//'.'
      end do
      dummy = trim(dummy)//'[OK]'
      print *, trim(dummy)
    end if
#endif
  end subroutine printDiag

  function getFMTForReal(value, w) result(FMT)
    implicit none
    real, intent(in) :: value
    character(len=STR_MAX) :: FMT
    integer, intent(in), optional :: w
    integer :: w_
    character(len=10) :: dummy
    if (.not. present(w)) then
      w_ = 10
    else
      w_ = w
    end if
    write (dummy, '(I10)') w_

    if ((abs(value) .ge. 100000) .or.&
      & ((abs(value) .lt. 1e-4) .and.&
        & (abs(value) .ne. 0.0))) then
      FMT = 'ES'//trim(dummy)//'.2'
    else
      FMT = 'F'//trim(dummy)//'.2'
    end if
  end function getFMTForReal

  subroutine printReport(bool, msg, prepend)
    implicit none
    character(len=*), intent(in) :: msg
    logical, intent(in) :: bool
    logical, optional, intent(in) :: prepend
    character(len=STR_MAX) :: dummy
    integer :: sz, i, ierr
    if (bool) then
      sz = len(trim(msg))
      if (present(prepend)) then
        if (prepend) then
          dummy = '...'
          sz = sz + 3
        end if
      else
        dummy = ''
      end if
      dummy = trim(dummy)//trim(msg)
      do i = 1, (38 - sz)
        dummy = trim(dummy)//'.'
      end do
      dummy = trim(dummy)//'[OK]'
      print *, trim(dummy)
    end if
  end subroutine printReport

  subroutine printTimeHeader(tstep)
    implicit none
    integer, intent(in) :: tstep
    character(len=STR_MAX) :: dummy
    integer :: sz, i

    ! printing divider
    do i = 72, 72
      dummy(i:i) = ' '
    end do
    do i = 1, 71
      dummy(i:i) = '-'
    end do
    print *, dummy(1:72)

    ! printing timestep
    sz = len(trim("Timestep: "//STR(tstep)))
    do i = 1, 71
      dummy(i:i) = '.'
    end do
    dummy(1:sz) = trim("Timestep: "//STR(tstep))
    dummy(66:71) = '[DONE]'
    print *, dummy(1:72)

    ! printing header
    do i = 1, 72
      dummy(i:i) = ' '
    end do
    dummy(1:71) = '[ROUTINE]          [TIME, ms]      [MIN  /  MAX, ms]      [FRACTION, %]'
    print *, dummy(1:72)
  end subroutine printTimeHeader

  subroutine printTimeFooter()
    implicit none
    character(len=STR_MAX) :: dummy
    integer :: i

    do i = 72, 72
      dummy(i:i) = ' '
    end do
    do i = 1, 71
      dummy(i:i) = '.'
    end do
    print *, dummy(1:72)
  end subroutine printTimeFooter

  subroutine printTime(dt_arr, msg, fullstep)
    implicit none
    character(len=*), intent(in) :: msg
    character(len=STR_MAX) :: dummy, dummy1, FMT
    real(kind=8), intent(in) :: dt_arr(:)
    real, optional, intent(in) :: fullstep
    real :: dt_mean, dt_max, dt_min
    integer :: sz, sz1, i
    dt_mean = SUM(dt_arr) * 1000 / mpi_size
    dt_max = MAXVAL(dt_arr) * 1000
    dt_min = MINVAL(dt_arr) * 1000
    if (present(fullstep)) then
      if (dt_mean / fullstep .lt. 1e-4) then
        dt_mean = 0; dt_min = 0; dt_max = 0
      end if
    end if

    do i = 1, 72
      dummy(i:i) = ' '
    end do

    sz = len(msg)
    dummy(1:sz) = msg

    FMT = "("//trim(getFMTForReal(dt_mean))//")"
    write (dummy1, FMT) dt_mean
    sz = len_trim(dummy1)
    dummy(20:20 + sz - 1) = trim(dummy1)

    FMT = "("//trim(getFMTForReal(dt_min))//")"
    write (dummy1, FMT) dt_min
    sz = len_trim(dummy1)
    dummy(32:32 + sz - 1) = trim(dummy1)

    FMT = "("//trim(getFMTForReal(dt_max))//")"
    write (dummy1, FMT) dt_max
    sz = len_trim(dummy1)
    dummy(43:43 + sz - 1) = trim(dummy1)
    if (present(fullstep)) then
      FMT = "("//trim(getFMTForReal(dt_mean * 100 / fullstep))//")"
      write (dummy1, FMT) dt_mean * 100 / fullstep
      sz1 = len_trim(dummy1)
      dummy(62:62 + sz1 - 1) = trim(dummy1)
    end if

    print *, dummy(1:72)
  end subroutine printTime

  function intToStr(my_int) result(string)
    implicit none
    integer, intent(in) :: my_int
    character(:), allocatable :: string
    character(len=STR_MAX) :: temp
    write (temp, '(i0)') my_int
    string = trim(temp)
  end function intToStr

  function realToStr(my_real) result(string)
    implicit none
    real, intent(in) :: my_real
    character(:), allocatable :: string
    character(len=STR_MAX) :: temp
    if ((my_real .ge. 1000) .or. ((my_real .lt. 1e-2) .and. (my_real .ne. 0.0))) then
      write (temp, '(ES10.2)') my_real
    else
      write (temp, '(F10.2)') my_real
    end if
    string = trim(temp)
  end function realToStr

  function STRtoINT(my_str) result(my_int)
    implicit none
    character(len=*), intent(in) :: my_str
    integer :: my_int
    read (my_str, *) my_int
  end function STRtoINT

  real(dprec) function randomNum(DSEED)
    implicit none
    real(dprec) :: DSEED
    integer :: I
    real(dprec) :: S2P31, S2P31M, SEED
    DATA S2P31M/2147483647.D0/, S2P31/2147483648.D0/
    SEED = DSEED
    SEED = DMOD(16807.D0 * SEED, S2P31M)
    randomNum = SEED / S2P31
    DSEED = SEED
    return
  end function randomNum

  real function random(DSEED)
    implicit none
    real(dprec) :: DSEED
    real :: rnd
    rnd = 1.0
    do while (rnd .eq. 1.0)
      rnd = randomNum(DSEED)
    end do
    random = rnd
    return
  end function random

  integer function randomInt(DSEED, amin, amax)
    implicit none
    real(dprec) :: DSEED
    integer, intent(in) :: amin, amax
    randomInt = amin + INT((amax - amin) * random(dseed))
    return
  end function randomInt

  real function poisson(num)
    implicit none
    real, intent(in) :: num
    real(kind=8) :: Lps, pps
    real :: kps, ups
    Lps = EXP(-REAL(num, 8))
    kps = 0
    pps = 1
    do while (pps .ge. Lps)
      kps = kps + 1
      ups = random(dseed)
      pps = pps * ups
    end do
    poisson = kps - 1
    return
  end function poisson

  subroutine initializeRandomSeed(rank)
    implicit none
    integer, intent(in) :: rank
    dseed = 123457.D0
    dseed = dseed + rank
  end subroutine initializeRandomSeed
end module m_aux
