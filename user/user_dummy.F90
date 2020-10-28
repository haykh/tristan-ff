#include "../src/defs.F90"

module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_fields
  implicit none

  !--- PRIVATE variables -----------------------------------------!

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
  end subroutine userReadInput


  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    bx(:,:,:) = -1.0; by(:,:,:) = -1.0; bz(:,:,:) = -1.0
    do i = 0, this_meshblock%ptr%sx - 1
      do j = 0, this_meshblock%ptr%sy - 1
        do k = 0, this_meshblock%ptr%sz - 1
          ex(i, j, k) = 0.0; ey(i, j, k) = 0.0; ez(i, j, k) = 0.0
          bx(i, j, k) = 0.0; by(i, j, k) = 0.0; bz(i, j, k) = 0.0
          b0x(i, j, k) = 0.0; b0y(i, j, k) = 0.0; b0z(i, j, k) = 1.0
        end do
      end do
    end do
  end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!

  !............................................................!

  !--- boundaries ---------------------------------------------!
  subroutine userFieldBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userFieldBoundaryConditions
  !............................................................!

  !--- custom output ------------------------------------------!
  subroutine userOutput(var, temp_, i, j, k)
    implicit none
    integer, intent(in)         :: var
    integer(kind=2), intent(in) :: i, j, k
    real, intent(out)           :: temp_
    if (var .eq. 1) then
      temp_ = 1.0
    else if (var .eq. 2) then
      temp_ = 0.0
    else
      temp_ = 0.0
    end if
  end subroutine userOutput
  !............................................................!

end module m_userfile
