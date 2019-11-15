#include "../defs.F90"

module m_finalize
  use m_globalnamespace
  use m_aux
  use m_domain
  use m_fields
  implicit none

  !--- PRIVATE functions -----------------------------------------!
  private :: finalizeCommunications, deallocateArrays
  !...............................................................!
contains
  subroutine finalizeAll()
    implicit none
    integer :: ierr
    call deallocateArrays()
      call printDiag((mpi_rank .eq. 0), "deallocateArrays()", .true.)
    call finalizeCommunications()
      call printDiag((mpi_rank .eq. 0), "finalizeCommunications()", .true.)

    call printReport((mpi_rank .eq. 0), "FinalizeAll()")
  end subroutine finalizeAll

  subroutine deallocateArrays()
    implicit none
    integer :: i, j, k, ierr

    ! dealloc meshblocks
    if (allocated(meshblocks)) deallocate(meshblocks)
    nullify(this_meshblock%ptr)

    ! dealloc field arrays
    if (allocated(ex)) deallocate(ex)
    if (allocated(ey)) deallocate(ey)
    if (allocated(ez)) deallocate(ez)
    if (allocated(bx)) deallocate(bx)
    if (allocated(by)) deallocate(by)
    if (allocated(bz)) deallocate(bz)

    if (allocated(b0x)) deallocate(b0x)
    if (allocated(b0y)) deallocate(b0y)
    if (allocated(b0z)) deallocate(b0z)

    if (allocated(enx)) deallocate(enx)
    if (allocated(eny)) deallocate(eny)
    if (allocated(enz)) deallocate(enz)
    if (allocated(bnx)) deallocate(bnx)
    if (allocated(bny)) deallocate(bny)
    if (allocated(bnz)) deallocate(bnz)

    if (allocated(dex)) deallocate(dex)
    if (allocated(dey)) deallocate(dey)
    if (allocated(dez)) deallocate(dez)
    if (allocated(dbx)) deallocate(dbx)
    if (allocated(dby)) deallocate(dby)
    if (allocated(dbz)) deallocate(dbz)

    if (allocated(rho)) deallocate(rho)

    ! dealloc field exchange
    if (allocated(send_fld)) deallocate(send_fld)
    if (allocated(recv_fld)) deallocate(recv_fld)
  end subroutine deallocateArrays

  subroutine finalizeCommunications()
    implicit none
    integer :: ierr
    call MPI_FINALIZE(ierr)
  end subroutine finalizeCommunications
end module m_finalize
