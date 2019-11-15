#include "defs.F90"

module m_mainloop
  use m_globalnamespace
  use m_helpers
  use m_aux
  use m_writeoutput
  use m_fldsolver
  use m_exchangefields
  use m_userfile
  use m_errors
  implicit none

  integer       :: timestep

  real(kind=8)  :: t_fullstep, t_outputstep,&
                 & t_fldexchstep, t_fldslvrstep,&
                 & t_usrfuncs


  !--- PRIVATE functions -----------------------------------------!
  private :: makeReport
  !...............................................................!

  !--- PRIVATE variables -----------------------------------------!
  private :: t_fullstep, t_outputstep,&
           & t_fldexchstep, t_fldslvrstep,&
           & t_usrfuncs
  !...............................................................!
contains
  subroutine mainloop()
    implicit none
    integer       :: ierr, i
    integer       :: s, ti, tj, tk, p

    ! ADD needs to be changed for restart
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call printReport((mpi_rank .eq. 0), "Starting mainloop()")

    t_fullstep = 0;       t_outputstep = 0
    t_fldexchstep = 0;    t_fldslvrstep = 0
    t_usrfuncs = 0

    ! filling ghost zones
    call exchangeFields(exchangeB0 = .true.)

    do timestep = 0, final_timestep
        t_fullstep = MPI_WTIME()

      ! MAINLOOP >
      !-------------------------------------------------
        t_usrfuncs = MPI_WTIME()
      call userFieldBoundaryConditions(timestep)
        t_usrfuncs = MPI_WTIME() - t_usrfuncs

        t_fldexchstep = MPI_WTIME() - t_fldexchstep
      call exchangeFields()
        t_fldexchstep = MPI_WTIME() - t_fldexchstep

      !-------------------------------------------------
      ! Output
      t_outputstep = 0
      if ((modulo(timestep, output_interval) .eq. 0) .and.&
        & (timestep .ge. output_start)) then
        t_outputstep = MPI_WTIME()
        call computeRho()
        call writeOutput(timestep)
        t_outputstep = MPI_WTIME() - t_outputstep
      end if
      !.................................................

        t_fldslvrstep = MPI_WTIME()
      call initRKstep()
      call rk3Step(1.0, 0.0, 1.0)
      call cleanEpar()
        t_fldslvrstep = MPI_WTIME() - t_fldslvrstep

        t_usrfuncs = MPI_WTIME()
      call userFieldBoundaryConditions(timestep)
        t_usrfuncs = MPI_WTIME() - t_usrfuncs

        t_fldexchstep = MPI_WTIME()
      call exchangeFields()
        t_fldexchstep = MPI_WTIME() - t_fldexchstep

        t_fldslvrstep = MPI_WTIME() - t_fldslvrstep
      call rk3Step(3.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0)
      call cleanEpar()
        t_fldslvrstep = MPI_WTIME() - t_fldslvrstep

        t_usrfuncs = MPI_WTIME()
      call userFieldBoundaryConditions(timestep)
        t_usrfuncs = MPI_WTIME() - t_usrfuncs

        t_fldexchstep = MPI_WTIME() - t_fldexchstep
      call exchangeFields()
        t_fldexchstep = MPI_WTIME() - t_fldexchstep

        t_fldslvrstep = MPI_WTIME() - t_fldslvrstep
      call rk3Step(1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0)
      call cleanEpar()
        t_fldslvrstep = MPI_WTIME() - t_fldslvrstep

        t_fldexchstep = MPI_WTIME() - t_fldexchstep
      call exchangeFields()
        t_fldexchstep = MPI_WTIME() - t_fldexchstep
      !.................................................

      ! </ MAINLOOP
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        t_fullstep = MPI_WTIME() - t_fullstep
      if (ierr .eq. MPI_SUCCESS) then
        call makeReport(timestep)
      end if
    end do
  end subroutine mainloop

  subroutine makeReport(tstep)
    implicit none
    integer, intent(in)           :: tstep
    integer                       :: ierr
    real                          :: fullstep
    real(kind=8), allocatable     :: dt_fullstep(:), dt_outputstep(:),&
                                   & dt_fldexchstep(:), dt_fldslvrstep(:),&
                                   & dt_usrfuncs(:)

    allocate(dt_fullstep(mpi_size), dt_outputstep(mpi_size))
    allocate(dt_fldexchstep(mpi_size), dt_fldslvrstep(mpi_size))
    allocate(dt_usrfuncs(mpi_size))

    call MPI_GATHER(t_fullstep, 1, MPI_REAL8,&
                  & dt_fullstep, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER(t_outputstep, 1, MPI_REAL8,&
                  & dt_outputstep, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER(t_fldexchstep, 1, MPI_REAL8,&
                  & dt_fldexchstep, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER(t_fldslvrstep, 1, MPI_REAL8,&
                  & dt_fldslvrstep, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER(t_usrfuncs, 1, MPI_REAL8,&
                  & dt_usrfuncs, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)

    if (mpi_rank .eq. 0) then
      fullstep = SUM(dt_fullstep) * 1000 / mpi_size
      call printReport(.true., "timestep: " // STR(tstep))
      call printTime(dt_fullstep, "Full_step: ", is_first_row = .true.)
      call printTime(dt_fldexchstep, "  fld_exchange: ", fullstep)
      call printTime(dt_fldslvrstep, "  fld_solver: ", fullstep)
      call printTime(dt_usrfuncs, "  usr_funcs: ", fullstep)
      call printTime(dt_outputstep, "  output_step: ", fullstep)
      print *, ""
    end if

    deallocate(dt_fullstep, dt_outputstep)
    deallocate(dt_fldexchstep, dt_fldslvrstep)
    deallocate(dt_usrfuncs)
  end subroutine makeReport

end module m_mainloop
