#include "../defs.F90"

module m_writeoutput
  #ifdef HDF5
    use hdf5
  #endif
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_fields
  use m_helpers
  use m_userfile
  implicit none

  integer                 :: output_start, output_interval, output_istep
  integer                 :: n_fld_vars
  character(len=STR_MAX)  :: fld_vars(100)

  !--- PRIVATE functions -----------------------------------------!
  #ifdef HDF5
    private :: writeFields_hdf5, writeDomain_hdf5, writeXDMF_hdf5
  #endif
  private :: initializeOutput
  !...............................................................!

  !--- PRIVATE variables -----------------------------------------!
  private :: n_fld_vars, fld_vars, n_dom_vars, dom_vars
  !...............................................................!
contains
  subroutine writeOutput(time)
    implicit none
    integer, intent(in)        :: time
    integer                    :: step, ierr
    call defineFieldVarsToOutput()
    step = output_index
    #ifdef HDF5
      call writeParams_hdf5(step, time)
        call printDiag((mpi_rank .eq. 0), "...writeParams_hdf5()", .true.)
      call writeFields_hdf5(step, time)
        call printDiag((mpi_rank .eq. 0), "...writeFields_hdf5()", .true.)
    #endif
    call printDiag((mpi_rank .eq. 0), "output()", .true.)
    output_index = output_index + 1
  end subroutine writeOutput

  subroutine defineFieldVarsToOutput()
    implicit none
    ! initialize field variables
    !   total number of fields
    n_fld_vars = 0
    fld_vars(n_fld_vars + 1 : n_fld_vars + 1 + 12) =&
                                 & (/'ex   ', 'ey   ', 'ez   ',&
                                   & 'bx   ', 'by   ', 'bz   ',&
                                   & 'jx   ', 'jy   ', 'jz   ',&
                                   & 'xx   ', 'yy   ', 'zz   '/)
    n_fld_vars = n_fld_vars + 12

    fld_vars(n_fld_vars + 1 : n_fld_vars + 1 + 4) = (/'curlBx', 'curlBy', 'curlBz', 'divE'/)
    n_fld_vars = n_fld_vars + 4

    fld_vars(n_fld_vars + 1 : n_fld_vars + 1 + 3) = (/'var1', 'var2', 'var3'/)
    n_fld_vars = n_fld_vars + 3
  end subroutine defineFieldVarsToOutput

  ! writes a field specified by `fld_var` from gridcell `i,j,k` ...
  ! ... to `sm_arr(i1, j1, k1)` with proper interpolation etc for the output
  subroutine selectFieldForOutput(fld_var, i1, j1, k1, i, j, k)
    implicit none
    character(len=STR_MAX), intent(in)  :: fld_var
    integer(kind=2), intent(in)         :: i1, j1, k1, i, j, k
    real                                :: ex0, ey0, ez0, bx0, by0, bz0, jx0, jy0, jz0
    real                                :: dx1, dx2, dy1, dy2, dz1, dz2, divE
    select case (trim(fld_var))
    case('ex')
      #ifndef DEBUG
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
      #else
        ex0 = ex(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = ex0
    case('ey')
      #ifndef DEBUG
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
      #else
        ey0 = ey(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = ey0
    case('ez')
      #ifndef DEBUG
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
      #else
        ez0 = ez(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = ez0
    case('bx')
      #ifndef DEBUG
        call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
      #else
        bx0 = bx(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = bx0
    case('by')
      #ifndef DEBUG
        call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
      #else
        by0 = by(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = by0
    case('bz')
      #ifndef DEBUG
        call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
      #else
        bz0 = bz(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = bz0
    case('jx')
      call computeForceFreeCurrent(jx0, jy0, jz0, i, j, k)
      sm_arr(i1, j1, k1) = jx0
    case('jy')
      call computeForceFreeCurrent(jx0, jy0, jz0, i, j, k)
      sm_arr(i1, j1, k1) = jy0
    case('jz')
      call computeForceFreeCurrent(jx0, jy0, jz0, i, j, k)
      sm_arr(i1, j1, k1) = jz0
    case('curlBx')
      dx1 = (bz(    i,    j,    k) - bz(    i,j - 1,    k)) - (by(    i,    j,    k) - by(    i,    j,k - 1))
      dx2 = (bz(i - 1,    j,    k) - bz(i - 1,j - 1,    k)) - (by(i - 1,    j,    k) - by(i - 1,    j,k - 1))
      sm_arr(i1, j1, k1) = 0.5 * (dx1 + dx2)
    case('curlBy')
      dy1 = (bx(    i,    j,    k) - bx(    i,    j,k - 1)) - (bz(    i,    j,    k) - bz(i - 1,    j,    k))
      dy2 = (bx(    i,j - 1,    k) - bx(    i,j - 1,k - 1)) - (bz(    i,j - 1,    k) - bz(i - 1,j - 1,    k))
      sm_arr(i1, j1, k1) = 0.5 * (dy1 + dy2)
    case('curlBz')
      dz1 = (by(    i,    j,    k) - by(i - 1,    j,    k)) - (bx(    i,    j,    k) - bx(    i,j - 1,    k))
      dz2 = (by(    i,    j,k - 1) - by(i - 1,    j,k - 1)) - (bx(    i,    j,k - 1) - bx(    i,j - 1,k - 1))
      sm_arr(i1, j1, k1) = 0.5 * (dz1 + dz2)
    case('divE')
      divE = (ex(i, j, k) - ex(i - 1, j, k)) +&
           & (ey(i, j, k) - ey(i, j - 1, k)) +&
           & (ez(i, j, k) - ez(i, j, k - 1))
      sm_arr(i1, j1, k1) = divE
    case('xx')
      sm_arr(i1, j1, k1) = REAL(this_meshblock%ptr%x0 + i, 4)
    case('yy')
      sm_arr(i1, j1, k1) = REAL(this_meshblock%ptr%y0 + j, 4)
    case('zz')
      sm_arr(i1, j1, k1) = REAL(this_meshblock%ptr%z0 + k, 4)
    case('var1')
      call userOutput(1, dummy1_, i, j, k)
      sm_arr(i1, j1, k1) = dummy1_
    case('var2')
      call userOutput(2, dummy1_, i, j, k)
      sm_arr(i1, j1, k1) = dummy1_
    case('var3')
      call userOutput(3, dummy1_, i, j, k)
      sm_arr(i1, j1, k1) = dummy1_
    case default
      call throwError("ERROR: unrecognized `fldname`")
    end select
  end subroutine selectFieldForOutput

  #ifdef HDF5

  subroutine writeParams_hdf5(step, time)
    implicit none
    integer, intent(in)               :: step, time
    character(len=STR_MAX)            :: stepchar, filename
    integer                           :: n, error, datarank
    integer(HID_T)                    :: file_id, dspace_id, dset_id
    integer(HSIZE_T), dimension(1)    :: data_dims
    integer, allocatable              :: data_int(:)
    real, allocatable                 :: data_real(:)
    character(len=STR_MAX)            :: dsetname

    if (mpi_rank .eq. 0) then
      datarank = 1
      data_dims(1) = 1
      allocate(data_real(1))
      allocate(data_int(1))

      write(stepchar, "(i5.5)") step
      filename = trim(output_dir_name) // '/params.' // trim(stepchar)

      dsetname = 'timestep'
      call h5open_f(error)
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      call h5screate_simple_f(datarank, data_dims, dspace_id, error)
      call h5dcreate_f(file_id, trim(dsetname), H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
      data_int(1) = time
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data_int, data_dims, error)
      call h5dclose_f(dset_id, error)
      call h5sclose_f(dspace_id, error)

      do n = 1, sim_params%count
        dsetname = trim(sim_params%param_group(n)%str) // ':' // trim(sim_params%param_name(n)%str)
        call h5screate_simple_f(datarank, data_dims, dspace_id, error)
        if (sim_params%param_type(n) .eq. 1) then
          call h5dcreate_f(file_id, trim(dsetname), H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
          data_int(1) = sim_params%param_value(n)%value_int
          call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data_int, data_dims, error)
        else if (sim_params%param_type(n) .eq. 2) then
          call h5dcreate_f(file_id, trim(dsetname), H5T_NATIVE_REAL, dspace_id, dset_id, error)
          data_real(1) = sim_params%param_value(n)%value_real
          call h5dwrite_f(dset_id, H5T_NATIVE_REAL, data_real, data_dims, error)
        else if (sim_params%param_type(n) .eq. 3) then
          call h5dcreate_f(file_id, trim(dsetname), H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
          data_int(1) = sim_params%param_value(n)%value_bool
          call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data_int, data_dims, error)
        else
          call throwError('ERROR. Unknown `param_type` in `saveAllParameters`.')
        end if
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)
      end do
      call h5fclose_f(file_id, error)
      call h5close_f(error)
    end if
  end subroutine writeParams_hdf5

  subroutine writeXDMF_hdf5(step, time, ni, nj, nk)
    implicit none
    integer, intent(in)               :: step, time, ni, nj, nk
    character(len=STR_MAX)            :: stepchar, filename
    integer                           :: var

    write(stepchar, "(i5.5)") step
    filename = trim(output_dir_name) // '/flds.tot.' // trim(stepchar) // '.xdmf'

    open (UNIT_xdmf, file=filename, status="replace", access="stream", form="formatted")
    write (UNIT_xdmf, "(A)")&
                           & '<?xml version="1.0" ?>'
    write (UNIT_xdmf, "(A)")&
                           & '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
    write (UNIT_xdmf, "(A)")&
                           & '<Xdmf Version="2.0">'
    write (UNIT_xdmf, "(A)")&
                           & '  <Domain>'
    write (UNIT_xdmf, "(A)")&
                           & '    <Grid Name="domain" GridType="Uniform">'
    write (UNIT_xdmf, "(A,A,A,I10,I10,I10,A)")&
                           & '      <Topology TopologyType=',&
                              & '"3DCoRectMesh"', ' Dimensions="', &
                              & nk, nj, ni, '"/>'
    write (UNIT_xdmf, "(A)")&
                           & '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
    write (UNIT_xdmf, "(A,A)")&
                           & '        <DataItem Format="XML" Dimensions="3"',&
                              & ' NumberType="Float" Precision="4">'
    write (UNIT_xdmf, "(A)")&
                           & '          0.0 0.0 0.0'
    write (UNIT_xdmf, "(A)")&
                           & '        </DataItem>'
    write (UNIT_xdmf, "(A,A)")&
                           & '        <DataItem Format="XML" Dimensions="3"',&
                              & ' NumberType="Float" Precision="4">'
    write (UNIT_xdmf, "(A)")&
                           & '          1.0 1.0 1.0'
    write (UNIT_xdmf, "(A)")&
                           & '        </DataItem>'
    write (UNIT_xdmf, "(A)")&
                           & '      </Geometry>'

    do var = 1, n_fld_vars
      write (UNIT_xdmf, "(A)")&
                           & '      <Attribute Name="' // trim(fld_vars(var)) // '" Center="Node">'
      write (UNIT_xdmf, "(A,I10,I10,I10,A)")&
                           & '        <DataItem Format="HDF" Dimensions="',&
                           & nk, nj, ni,&
                           & '" NumberType="Float" Precision="4">'
      write (UNIT_xdmf, "(A,A)")&
                           & '          flds.tot.' // trim(stepchar) // ':/',&
                           & trim(fld_vars(var))
      write (UNIT_xdmf, "(A)")&
                           & '        </DataItem>'
      write (UNIT_xdmf, "(A)")&
                           & '      </Attribute>'
    end do

    write (UNIT_xdmf, "(A)")&
                           & '    </Grid>'
    write (UNIT_xdmf, "(A)")&
                           & '  </Domain>'
    write (UNIT_xdmf, "(A)")&
                           & '</Xdmf>'
    close (UNIT_xdmf)
  end subroutine writeXDMF_hdf5

  subroutine writeFields_hdf5(step, time)
    implicit none
    integer, intent(in)               :: step, time
    character(len=STR_MAX)            :: stepchar, filename
    integer(HID_T)                    :: file_id, dset_id(40), filespace(40), memspace, plist_id
    integer                           :: error, f, s
    integer                           :: i, j, k
    integer                           :: dataset_rank = 3
    integer(HSSIZE_T), dimension(3)   :: offsets
    integer(HSIZE_T), dimension(3)    :: global_dims, blocks
    real                              :: ex0, ey0, ez0, bx0, by0, bz0, dummy1_
    real                              :: jx0, jy0, jz0
    real, allocatable                 :: sm_arr(:,:,:)

    ! downsampling variables
    integer :: this_x0, this_y0, this_z0, this_sx, this_sy, this_sz
    integer :: i_start, i_end, j_start, j_end, k_start, k_end
    integer :: offset_i, offset_j, offset_k, i1, j1, k1
    integer :: n_i, n_j, n_k, glob_n_i, glob_n_j, glob_n_k

    ! for convenience
    this_x0 = this_meshblock%ptr%x0
    this_y0 = this_meshblock%ptr%y0
    this_z0 = this_meshblock%ptr%z0
    this_sx = this_meshblock%ptr%sx
    this_sy = this_meshblock%ptr%sy
    this_sz = this_meshblock%ptr%sz

    write(stepchar, "(i5.5)") step
    filename = trim(output_dir_name) // '/flds.tot.' // trim(stepchar)

    ! assuming `global_mesh%{x0,y0,z0} .eq. 0`
    ! if debug enabled write also the outermost ghost zones
    if (output_istep .eq. 1) then
      #ifndef DEBUG
        offset_i = this_x0;   offset_j = this_y0;   offset_k = this_z0
        n_i = this_sx - 1;    n_j = this_sy - 1;    n_k = this_sz - 1
        glob_n_i = global_mesh%sx
        glob_n_j = global_mesh%sy
        glob_n_k = global_mesh%sz

        i_start = 0; j_start = 0; k_start = 0
      #else
        ! also save ghosts
        if (this_x0 .eq. 0) then ! -x end
          offset_i = 0
          n_i = this_sx - 1 + NGHOST
          i_start = -NGHOST
        else if (this_x0 + this_sx .eq. global_mesh%sx) then ! +x end
          offset_i = this_x0 + NGHOST
          n_i = this_sx - 1 + NGHOST
          i_start = 0
        else
          offset_i = this_x0 + NGHOST
          n_i = this_sx - 1
          i_start = 0
        end if

        if (this_y0 .eq. 0) then ! -y end
          offset_j = 0
          n_j = this_sy - 1 + NGHOST
          j_start = -NGHOST
        else if (this_y0 + this_sy .eq. global_mesh%sy) then ! +y end
          offset_j = this_y0 + NGHOST
          n_j = this_sy - 1 + NGHOST
          j_start = 0
        else
          offset_j = this_y0 + NGHOST
          n_j = this_sy - 1
          j_start = 0
        end if

        if (this_z0 .eq. 0) then ! -z end
          offset_k = 0
          n_k = this_sz - 1 + NGHOST
          k_start = -NGHOST
        else if (this_z0 + this_sz .eq. global_mesh%sz) then ! +z end
          offset_k = this_z0 + NGHOST
          n_k = this_sz - 1 + NGHOST
          k_start = 0
        else
          offset_k = this_z0 + NGHOST
          n_k = this_sz - 1
          k_start = 0
        end if

        glob_n_i = global_mesh%sx + 2 * NGHOST
        glob_n_j = global_mesh%sy + 2 * NGHOST
        glob_n_k = global_mesh%sz + 2 * NGHOST
      #endif
    else
      offset_i = CEILING(REAL(this_x0) / REAL(output_istep))
      offset_j = CEILING(REAL(this_y0) / REAL(output_istep))

      i_start = CEILING(REAL(this_x0) / REAL(output_istep)) * output_istep - this_x0
      i_end = (CEILING(REAL(this_x0 + this_sx) / REAL(output_istep)) - 1) * output_istep - this_x0

      j_start = CEILING(REAL(this_y0) / REAL(output_istep)) * output_istep - this_y0
      j_end = (CEILING(REAL(this_y0 + this_sy) / REAL(output_istep)) - 1) * output_istep - this_y0

      n_i = (i_end - i_start) / output_istep
      n_j = (j_end - j_start) / output_istep

      glob_n_i = CEILING(REAL(global_mesh%sx) / REAL(output_istep))
      glob_n_i = MAX(1, glob_n_i)

      glob_n_j = CEILING(REAL(global_mesh%sy) / REAL(output_istep))
      glob_n_j = MAX(1, glob_n_j)

      offset_k = CEILING(REAL(this_z0) / REAL(output_istep))
      k_start = CEILING(REAL(this_z0) / REAL(output_istep)) * output_istep - this_z0
      k_end = (CEILING(REAL(this_z0 + this_sz) / REAL(output_istep)) - 1) * output_istep - this_z0
      n_k = (k_end - k_start) / output_istep
      glob_n_k = CEILING(REAL(global_mesh%sz) / REAL(output_istep))
      glob_n_k = MAX(1, glob_n_k)
    end if

    if (mpi_rank .eq. 0) then
      call writeXDMF_hdf5(step, time, glob_n_i, glob_n_j, glob_n_k)
    end if

    allocate(sm_arr(0:n_i, 0:n_j, 0:n_k))

    offsets(1) = offset_i
    offsets(2) = offset_j
    offsets(3) = offset_k
    blocks(1) = n_i + 1
    blocks(2) = n_j + 1
    blocks(3) = n_k + 1
    global_dims(1) = glob_n_i
    global_dims(2) = glob_n_j
    global_dims(3) = glob_n_k

    call h5open_f(error)
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, h5comm, h5info, error)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id, error)

    do f = 1, n_fld_vars
      call h5screate_simple_f(dataset_rank, global_dims, filespace(f), error)
    end do

    do f = 1, n_fld_vars
      call h5dcreate_f(file_id, fld_vars(f), H5T_NATIVE_REAL, filespace(f), &
                     & dset_id(f), error)
      call h5sclose_f(filespace(f), error)
    end do

    call h5screate_simple_f(dataset_rank, blocks, memspace, error)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

    do f = 1, n_fld_vars
      call h5dget_space_f(dset_id(f), filespace(f), error)
      call h5sselect_hyperslab_f(filespace(f), H5S_SELECT_SET_F, offsets, blocks, error)

      ! Create dataset by interpolating fields
      do i1 = 0, n_i
        do j1 = 0, n_j
          do k1 = 0, n_k
            i = i_start + i1 * output_istep
            j = j_start + j1 * output_istep
            k = k_start + k1 * output_istep
            call selectFieldForOutput(fld_vars(f), i1, j1, k1, i, j, k)
          end do
        end do
      end do

      ! Write the dataset collectively
      call h5dwrite_f(dset_id(f), H5T_NATIVE_REAL, sm_arr(:,:,:), global_dims, error, &
                    & file_space_id = filespace(f), mem_space_id = memspace, xfer_prp = plist_id)
    end do

    do f = 1, n_fld_vars
      call h5sclose_f(filespace(f), error)
    end do
    call h5sclose_f(memspace, error)

    do f = 1, n_fld_vars
      call h5dclose_f(dset_id(f), error)
    end do

    call h5pclose_f(plist_id, error)
    call h5fclose_f(file_id, error)
    call h5close_f(error)

    deallocate(sm_arr)
  end subroutine writeFields_hdf5

  #endif

end module m_writeoutput
