#include "defs.F90"

program tristan
  use mpi_f08
  use hdf5
  use m_initialize
  use m_mainloop
  use m_finalize
  implicit none
  !----- main code --------------------------

  call initializeAll()
  call mainloop()
  call finalizeAll()

  !..... main code ..........................
end program tristan
