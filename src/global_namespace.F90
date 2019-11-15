#include "defs.F90"

!--- GLOBAL_NAMESPACE ------------------------------------------!
! To store predefined(!) parameters and functions/interfaces
!   and share them between modules
!...............................................................!

module m_globalnamespace
  #ifdef MPI08
    use mpi_f08
  #endif

  #ifdef MPI
    include "mpif.h"
  #endif

  integer, parameter     :: dprec = kind(1.0d0)
  integer, parameter     :: sprec = kind(1.0e0)
  integer, parameter     :: UNIT_input = 10, UNIT_output = 20, UNIT_history = 30

  ! simulation parameters
  integer                :: final_timestep, output_index = 0
  character(len=STR_MAX) :: input_file_name = 'input',&
                          & output_dir_name = 'output',&
                          & restart_dir_name = 'restart'

  ! mpi variables
  integer       :: mpi_rank, mpi_size, mpi_statsize
  integer       :: sizex, sizey, sizez
  #ifdef HDF5
    integer, parameter  :: UNIT_xdmf = 40
    integer             :: h5comm, h5info
  #endif
end module m_globalnamespace