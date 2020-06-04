#include "../defs.F90"

module m_initialize
  #ifdef IFPORT
      use ifport, only : makedirqq
    #endif
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_writeoutput
  use m_domain
  use m_fields
  use m_userfile
  use m_helpers
  use m_errors
  implicit none

  !--- PRIVATE functions -----------------------------------------!
  private :: initializeCommunications, initializeOutput,&
           & firstRankInitialize, initializeFields,&
           & distributeMeshblocks, initializeDomain,&
           & initializeSimulation, checkEverything,&
           & allocateFieldArray, printParams
  !...............................................................!
contains
  ! initialize all the necessary arrays and variables
  subroutine initializeAll()
    implicit none
    call readCommandlineArgs()

    ! initializing the simulation parameters class ...
    ! ... which stores all the input values for the simulation
    sim_params%count = 0
    allocate(sim_params%param_type(1000))
    allocate(sim_params%param_group(1000))
    allocate(sim_params%param_name(1000))
    allocate(sim_params%param_value(1000))

    call initializeDomain()

    call initializeCommunications()
      call printDiag((mpi_rank .eq. 0), "initializeCommunications()", .true.)

    call distributeMeshblocks()
      call printDiag((mpi_rank .eq. 0), "distributeMeshblocks()", .true.)

    call initializeSimulation()
      call printDiag((mpi_rank .eq. 0), "initializeSimulation()", .true.)

    call initializeOutput()
      call printDiag((mpi_rank .eq. 0), "initializeOutput()", .true.)

    call initializeFields()
      call printDiag((mpi_rank .eq. 0), "initializeFields()", .true.)

    call initializeRandomSeed(mpi_rank)
      call printDiag((mpi_rank .eq. 0), "initializeRandomSeed()", .true.)

    if (mpi_rank .eq. 0) then
      call firstRankInitialize()
      call printDiag(.true., "firstRankInitialize()", .true.)
    end if

    call userInitialize()
      call printDiag((mpi_rank .eq. 0), "userInitialize()", .true.)

    call checkEverything()
      call printDiag((mpi_rank .eq. 0), "checkEverything()", .true.)

    call printParams()

    call printReport((mpi_rank .eq. 0), "InitializeAll()")
  end subroutine initializeAll

  subroutine printParams()
    implicit none
    integer                 :: n
    character(len=STR_MAX)  :: FMT
    ! printing simulation parameters in the report

    if (mpi_rank .eq. 0) then
      FMT = '== Simulation parameters ==============================================='
      write(*, '(A)') trim(FMT)
      do n = 1, sim_params%count
        if (sim_params%param_type(n) .eq. 1) then
          FMT = '(A30,A1,A20,A1,I10)'
          write (*, FMT) trim(sim_params%param_group(n)%str), ':',&
                       & trim(sim_params%param_name(n)%str), ':',&
                       & sim_params%param_value(n)%value_int
        else if (sim_params%param_type(n) .eq. 2) then
          FMT = getFMTForReal(sim_params%param_value(n)%value_real)
          FMT = '(A30,A1,A20,A1,' // trim(FMT) // ')'
          write (*, FMT) trim(sim_params%param_group(n)%str), ':',&
                       & trim(sim_params%param_name(n)%str), ':',&
                       & sim_params%param_value(n)%value_real
        else if (sim_params%param_type(n) .eq. 3) then
          FMT = '(A30,A1,A20,A1,L10)'
          write (*, FMT) trim(sim_params%param_group(n)%str), ':',&
                       & trim(sim_params%param_name(n)%str), ':',&
                       & sim_params%param_value(n)%value_bool
        else
          call throwError('ERROR. Unknown `param_type` in `saveAllParameters`.')
        end if
      end do
      FMT = '........................................................................'
      write(*, '(A)') trim(FMT)
    end if
  end subroutine printParams


  subroutine initializeCommunications()
    implicit none
    integer :: ierr
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)
    mpi_statsize = MPI_STATUS_SIZE
    if (mpi_size .ne. sizex * sizey * sizez) then
      call throwError('ERROR: # of processors is not equal to the number of processors from input')
    end if
  end subroutine initializeCommunications

  subroutine initializeDomain()
    implicit none
    call getInput('node_configuration', 'sizex', sizex)
    call getInput('node_configuration', 'sizey', sizey)
    call getInput('node_configuration', 'sizez', sizez)
    global_mesh%x0 = 0
    global_mesh%y0 = 0
    global_mesh%z0 = 0
    call getInput('grid', 'mx0', global_mesh%sx)
    call getInput('grid', 'my0', global_mesh%sy)
    call getInput('grid', 'mz0', global_mesh%sz)

    if ((modulo(global_mesh%sx, sizex) .ne. 0) .or.&
      & (modulo(global_mesh%sy, sizey) .ne. 0) .or.&
      & (modulo(global_mesh%sz, sizez) .ne. 0)) then
      call throwError('ERROR: grid size is not evenly divisible by the number of cores')
    end if

    call getInput('grid', 'boundary_Px', boundary_Px)
    call getInput('grid', 'boundary_Mx', boundary_Mx)
    if ((boundary_Px .eq. 1) .or. (boundary_Mx .eq. 1)) then
      boundary_Px = 1; boundary_Mx = 1
    end if
    call getInput('grid', 'boundary_Py', boundary_Py)
    call getInput('grid', 'boundary_My', boundary_My)
    if ((boundary_Py .eq. 1) .or. (boundary_My .eq. 1)) then
      boundary_Py = 1; boundary_My = 1
    end if
    call getInput('grid', 'boundary_Pz', boundary_Pz)
    call getInput('grid', 'boundary_Mz', boundary_Mz)
    if ((boundary_Pz .eq. 1) .or. (boundary_Mz .eq. 1)) then
      boundary_Pz = 1; boundary_Mz = 1
    end if
    call getInput('grid', 'ds_abs', ds_abs, 5.0)
  end subroutine initializeDomain

  subroutine distributeMeshblocks()
    implicit none
    integer, dimension(3) :: ind, m
    integer               :: rnk
    m(1) = global_mesh%sx / sizex
    m(2) = global_mesh%sy / sizey
    m(3) = global_mesh%sz / sizez
    allocate(meshblocks(mpi_size))
    this_meshblock%ptr => meshblocks(mpi_rank + 1)
    do rnk = 0, mpi_size - 1
      ind = rnkToInd(rnk)
      meshblocks(rnk + 1)%rnk = rnk
      ! find sizes and corner coords
      meshblocks(rnk + 1)%sx = m(1)
      meshblocks(rnk + 1)%sy = m(2)
      meshblocks(rnk + 1)%sz = m(3)
      meshblocks(rnk + 1)%x0 = ind(1) * m(1) + global_mesh%x0
      meshblocks(rnk + 1)%y0 = ind(2) * m(2) + global_mesh%y0
      meshblocks(rnk + 1)%z0 = ind(3) * m(3) + global_mesh%z0
    end do
    ! assign all neighbors
    call reassignNeighborsForAll()
  end subroutine distributeMeshblocks

  subroutine initializeOutput()
    implicit none
    call getInput('output', 'start', output_start, 0)
    call getInput('output', 'interval', output_interval, 10)
    call getInput('output', 'istep', output_istep, 4)

    #ifdef HDF5

    #ifdef MPI08
      h5comm = MPI_COMM_WORLD%MPI_VAL
      h5info = MPI_INFO_NULL%MPI_VAL
    #endif

    #ifdef MPI
      h5comm = MPI_COMM_WORLD
      h5info = MPI_INFO_NULL
    #endif

    #endif
  end subroutine initializeOutput

  subroutine initializeSimulation()
    implicit none
    call getInput('time', 'last', final_timestep, 1000)
    call getInput('algorithm', 'c', CC, 0.45)
  end subroutine initializeSimulation

  subroutine allocateFieldArray(fld)
    implicit none
    real, allocatable, intent(inout) :: fld(:,:,:)
    if (allocated(fld)) deallocate(fld)
    allocate(fld(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
               & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
               & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
  end subroutine allocateFieldArray

  subroutine initializeFields()
    implicit none
    call allocateFieldArray(ex)
    call allocateFieldArray(ey)
    call allocateFieldArray(ez)
    call allocateFieldArray(bx)
    call allocateFieldArray(by)
    call allocateFieldArray(bz)

    call allocateFieldArray(b0x)
    call allocateFieldArray(b0y)
    call allocateFieldArray(b0z)

    call allocateFieldArray(dex)
    call allocateFieldArray(dey)
    call allocateFieldArray(dez)
    call allocateFieldArray(dbx)
    call allocateFieldArray(dby)
    call allocateFieldArray(dbz)

    call allocateFieldArray(enx)
    call allocateFieldArray(eny)
    call allocateFieldArray(enz)
    call allocateFieldArray(bnx)
    call allocateFieldArray(bny)
    call allocateFieldArray(bnz)

    call allocateFieldArray(rho)

    ! exchange fields
    ! 20 = max # of fields sent in each direction
    sendrecv_offsetsz = MAX0(this_meshblock%ptr%sx, this_meshblock%ptr%sy, this_meshblock%ptr%sz)**2 * NGHOST * 20
    sendrecv_buffsz = sendrecv_offsetsz * 30
    ! 26 (~30) directions to send/recv in 3D

    if (allocated(send_fld)) deallocate(send_fld)
    allocate(send_fld(sendrecv_buffsz))
    if (allocated(recv_fld)) deallocate(recv_fld)
    allocate(recv_fld(sendrecv_offsetsz))
  end subroutine initializeFields

  subroutine firstRankInitialize()
    ! create output/restart directories
    !   if does not already exist
    !     note: some compilers may not support IFPORT
    #ifdef IFPORT
      logical :: result
      result = makedirqq(trim(output_dir_name))
      ! result = makedirqq(trim(restart_dir_name))
    #else
      call system('mkdir -p ' // trim(output_dir_name))
      ! call system('mkdir -p ' // trim(restart_dir_name))
    #endif
  end subroutine firstRankInitialize

  subroutine checkEverything()
    implicit none
    ! check that the domain size is larger than the number of ghost zones
    if ((this_meshblock%ptr%sx .lt. NGHOST) .or.&
      & (this_meshblock%ptr%sy .lt. NGHOST) .or.&
      & (this_meshblock%ptr%sz .lt. NGHOST)) then
      call throwError('ERROR: ghost zones overflow the domain size in ' // trim(STR(mpi_rank)))
    end if
  end subroutine checkEverything
end module m_initialize
