#include "../defs.F90"

module m_writeoutput
  use hdf5
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_fields
  use m_helpers
  use m_userfile
  implicit none

  integer :: output_start, output_interval, output_istep
  integer :: n_fld_vars, n_dom_vars
  character(len=STR_MAX) :: fld_vars(100), dom_vars(100)

  !--- PRIVATE functions -----------------------------------------!
  private :: writeFields_hdf5, writeDomain_hdf5, writeXDMF_hdf5
  private :: initializeOutput
  !...............................................................!

  !--- PRIVATE variables -----------------------------------------!
  private :: n_fld_vars, fld_vars, n_dom_vars, dom_vars
  !...............................................................!
contains
  subroutine writeOutput(time)
    implicit none
    integer, intent(in) :: time
    integer :: step, ierr

    call initializeOutput()

    step = output_index
    call writeParams_hdf5(step, time)
    call printDiag((mpi_rank .eq. 0), "...writeParams_hdf5()", .true.)
    call writeFields_hdf5(step, time)
    call printDiag((mpi_rank .eq. 0), "...writeFields_hdf5()", .true.)
    call writeDomain_hdf5(step, time)
    call printDiag((mpi_rank .eq. 0), "...writeDomain_hdf5()", .true.)
    call printDiag((mpi_rank .eq. 0), "output()", .true.)
    output_index = output_index + 1
  end subroutine writeOutput

  subroutine initializeOutput()
    implicit none
    ! initialize field variables
    n_fld_vars = 12
    fld_vars(1:n_fld_vars) = (/'ex   ', 'ey   ', 'ez   ',&
                               & 'bx   ', 'by   ', 'bz   ',&
                               & 'var1 ', 'var2 ', 'var3 ',&
                               & 'xx   ', 'yy   ', 'zz   '/)

    ! initialize domain output variables
    n_dom_vars = 6
    dom_vars(1:6) = (/'x0   ', 'y0   ', 'z0   ',&
                      & 'sx   ', 'sy   ', 'sz   '/)
  end subroutine initializeOutput

  subroutine writeParams_hdf5(step, time)
    implicit none
    integer, intent(in) :: step, time
    character(len=STR_MAX) :: stepchar, filename
    integer :: n, error, datarank
    integer(HID_T) :: file_id, dspace_id, dset_id
    integer(HSIZE_T), dimension(1) :: data_dims
    integer, allocatable :: data_int(:)
    real, allocatable :: data_real(:)
    character(len=STR_MAX) :: dsetname

    if (mpi_rank .eq. 0) then
      datarank = 1
      data_dims(1) = 1
      allocate (data_real(1))
      allocate (data_int(1))

      write (stepchar, "(i5.5)") step
      filename = trim(output_dir_name)//'/params.'//trim(stepchar)

      dsetname = 'timestep'
      call h5open_f(error)
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      call h5screate_simple_f(datarank, data_dims, dspace_id, error)
      call h5dcreate_f(file_id, trim(dsetname), H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
      data_int(1) = time
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data_int, data_dims, error)
      call h5dclose_f(dset_id, error)
      call h5sclose_f(dspace_id, error)

      do n = 1, sim_params % count
        dsetname = trim(sim_params % param_group(n) % str)//':'//trim(sim_params % param_name(n) % str)
        call h5screate_simple_f(datarank, data_dims, dspace_id, error)
        if (sim_params % param_type(n) .eq. 1) then
          call h5dcreate_f(file_id, trim(dsetname), H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
          data_int(1) = sim_params % param_value(n) % value_int
          call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data_int, data_dims, error)
        else if (sim_params % param_type(n) .eq. 2) then
          call h5dcreate_f(file_id, trim(dsetname), H5T_NATIVE_REAL, dspace_id, dset_id, error)
          data_real(1) = sim_params % param_value(n) % value_real
          call h5dwrite_f(dset_id, H5T_NATIVE_REAL, data_real, data_dims, error)
        else if (sim_params % param_type(n) .eq. 3) then
          call h5dcreate_f(file_id, trim(dsetname), H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
          data_int(1) = sim_params % param_value(n) % value_bool
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
    integer, intent(in) :: step, time, ni, nj, nk
    character(len=STR_MAX) :: stepchar, filename
    integer :: var

    write (stepchar, "(i5.5)") step
    filename = trim(output_dir_name)//'/flds.tot.'//trim(stepchar)//'.xdmf'

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
                           & '      <Attribute Name="'//trim(fld_vars(var))//'" Center="Node">'
      write (UNIT_xdmf, "(A,I10,I10,I10,A)")&
                           & '        <DataItem Format="HDF" Dimensions="',&
                           & nk, nj, ni,&
                           & '" NumberType="Float" Precision="4">'
      write (UNIT_xdmf, "(A,A)")&
                           & '          flds.tot.'//trim(stepchar)//':/',&
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
    integer, intent(in) :: step, time
    character(len=STR_MAX) :: stepchar, filename
    integer(HID_T) :: file_id, dset_id(40), filespace(40), memspace, plist_id
    integer :: error, f, s
    integer :: i, j, k
    integer :: dataset_rank = 3
    integer(HSSIZE_T), dimension(3) :: offsets
    integer(HSIZE_T), dimension(3) :: global_dims, blocks
    real :: ex0, ey0, ez0, bx0, by0, bz0, dummy1_
    real :: jx0, jy0, jz0
    real, allocatable :: sm_arr(:, :, :)

    ! downsampling variables
    integer :: this_x0, this_y0, this_z0, this_sx, this_sy, this_sz
    integer :: i_start, i_end, j_start, j_end, k_start, k_end
    integer :: offset_i, offset_j, offset_k, i1, j1, k1
    integer :: n_i, n_j, n_k, glob_n_i, glob_n_j, glob_n_k

    ! for convenience
    this_x0 = this_meshblock % ptr % x0
    this_y0 = this_meshblock % ptr % y0
    this_z0 = this_meshblock % ptr % z0
    this_sx = this_meshblock % ptr % sx
    this_sy = this_meshblock % ptr % sy
    this_sz = this_meshblock % ptr % sz

    write (stepchar, "(i5.5)") step
    filename = trim(output_dir_name)//'/flds.tot.'//trim(stepchar)

    ! assuming `global_mesh%{x0,y0,z0} .eq. 0`
    ! if debug enabled write also the outermost ghost zones
    if (output_istep .eq. 1) then
      offset_i = this_x0; offset_j = this_y0; offset_k = this_z0
      n_i = this_sx - 1; n_j = this_sy - 1; n_k = this_sz - 1
      glob_n_i = global_mesh % sx
      glob_n_j = global_mesh % sy
      glob_n_k = global_mesh % sz

      i_start = 0; j_start = 0; k_start = 0
      ! also save ghosts
      if (this_x0 .eq. 0) then ! -x end
        offset_i = 0
        n_i = this_sx - 1 + NGHOST
        i_start = -NGHOST
      else if (this_x0 + this_sx .eq. global_mesh % sx) then ! +x end
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
      else if (this_y0 + this_sy .eq. global_mesh % sy) then ! +y end
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
      else if (this_z0 + this_sz .eq. global_mesh % sz) then ! +z end
        offset_k = this_z0 + NGHOST
        n_k = this_sz - 1 + NGHOST
        k_start = 0
      else
        offset_k = this_z0 + NGHOST
        n_k = this_sz - 1
        k_start = 0
      end if

      glob_n_i = global_mesh % sx + 2 * NGHOST
      glob_n_j = global_mesh % sy + 2 * NGHOST
      glob_n_k = global_mesh % sz + 2 * NGHOST
    else
      offset_i = CEILING(REAL(this_x0) / REAL(output_istep))
      offset_j = CEILING(REAL(this_y0) / REAL(output_istep))

      i_start = CEILING(REAL(this_x0) / REAL(output_istep)) * output_istep - this_x0
      i_end = (CEILING(REAL(this_x0 + this_sx) / REAL(output_istep)) - 1) * output_istep - this_x0

      j_start = CEILING(REAL(this_y0) / REAL(output_istep)) * output_istep - this_y0
      j_end = (CEILING(REAL(this_y0 + this_sy) / REAL(output_istep)) - 1) * output_istep - this_y0

      n_i = (i_end - i_start) / output_istep
      n_j = (j_end - j_start) / output_istep

      glob_n_i = CEILING(REAL(global_mesh % sx) / REAL(output_istep))
      glob_n_i = MAX(1, glob_n_i)

      glob_n_j = CEILING(REAL(global_mesh % sy) / REAL(output_istep))
      glob_n_j = MAX(1, glob_n_j)

      offset_k = CEILING(REAL(this_z0) / REAL(output_istep))
      k_start = CEILING(REAL(this_z0) / REAL(output_istep)) * output_istep - this_z0
      k_end = (CEILING(REAL(this_z0 + this_sz) / REAL(output_istep)) - 1) * output_istep - this_z0
      n_k = (k_end - k_start) / output_istep
      glob_n_k = CEILING(REAL(global_mesh % sz) / REAL(output_istep))
      glob_n_k = MAX(1, glob_n_k)
    end if

    if (mpi_rank .eq. 0) then
      call writeXDMF_hdf5(step, time, glob_n_i, glob_n_j, glob_n_k)
    end if

    allocate (sm_arr(0:n_i, 0:n_j, 0:n_k))

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
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp=plist_id)
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
            select case (trim(fld_vars(f)))
            case ('ex')
#ifndef DEBUG
              call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
#else
              ex0 = ex(i, j, k)
#endif
              sm_arr(i1, j1, k1) = ex0
            case ('ey')
#ifndef DEBUG
              call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
#else
              ey0 = ey(i, j, k)
#endif
              sm_arr(i1, j1, k1) = ey0
            case ('ez')
#ifndef DEBUG
              call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
#else
              ez0 = ez(i, j, k)
#endif
              sm_arr(i1, j1, k1) = ez0
            case ('bx')
#ifndef DEBUG
              call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
#else
              bx0 = bx(i, j, k)
#endif
              sm_arr(i1, j1, k1) = bx0
            case ('by')
#ifndef DEBUG
              call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
#else
              by0 = by(i, j, k)
#endif
              sm_arr(i1, j1, k1) = by0
            case ('bz')
#ifndef DEBUG
              call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
#else
              bz0 = bz(i, j, k)
#endif
              sm_arr(i1, j1, k1) = bz0
            case ('var1')
              call userOutput(1, dummy1_, i, j, k)
              sm_arr(i1, j1, k1) = dummy1_
            case ('var2')
              call userOutput(2, dummy1_, i, j, k)
              sm_arr(i1, j1, k1) = dummy1_
            case ('var3')
              call userOutput(3, dummy1_, i, j, k)
              sm_arr(i1, j1, k1) = dummy1_
            case ('xx')
              sm_arr(i1, j1, k1) = REAL(this_meshblock % ptr % x0 + i, 4)
            case ('yy')
              sm_arr(i1, j1, k1) = REAL(this_meshblock % ptr % y0 + j, 4)
            case ('zz')
              sm_arr(i1, j1, k1) = REAL(this_meshblock % ptr % z0 + k, 4)
            case default
              call throwError("ERROR: unrecognized `fld_vars(f)`")
            end select
          end do
        end do
      end do

      ! Write the dataset collectively
      call h5dwrite_f(dset_id(f), H5T_NATIVE_REAL, sm_arr(:, :, :), global_dims, error, &
                    & file_space_id=filespace(f), mem_space_id=memspace, xfer_prp=plist_id)
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

    deallocate (sm_arr)
  end subroutine writeFields_hdf5

  subroutine writeDomain_hdf5(step, time)
    implicit none
    integer, intent(in) :: step, time
    character(len=STR_MAX) :: stepchar, filename
    integer :: error, s, i, datarank, d, rnk
    integer(HID_T) :: file_id, dset_id, dspace_id
    integer(HSIZE_T), dimension(1) :: data_dims
    integer, allocatable, dimension(:) :: domain_data

    datarank = 1
    data_dims(1) = mpi_size

    ! only root rank writes the domain file
    if (mpi_rank .eq. 0) then
      write (stepchar, "(i5.5)") step
      filename = trim(output_dir_name)//'/domain.'//trim(stepchar)

      ! Initialize FORTRAN interface
      call h5open_f(error)
      ! Create a new file using default properties
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      allocate (domain_data(mpi_size))

      do d = 1, n_dom_vars
        select case (trim(dom_vars(d)))
        case ('x0')
          do rnk = 0, mpi_size - 1
            domain_data(rnk + 1) = meshblocks(rnk + 1) % x0
          end do
        case ('y0')
          do rnk = 0, mpi_size - 1
            domain_data(rnk + 1) = meshblocks(rnk + 1) % y0
          end do
        case ('z0')
          do rnk = 0, mpi_size - 1
            domain_data(rnk + 1) = meshblocks(rnk + 1) % z0
          end do
        case ('sx')
          do rnk = 0, mpi_size - 1
            domain_data(rnk + 1) = meshblocks(rnk + 1) % sx
          end do
        case ('sy')
          do rnk = 0, mpi_size - 1
            domain_data(rnk + 1) = meshblocks(rnk + 1) % sy
          end do
        case ('sz')
          do rnk = 0, mpi_size - 1
            domain_data(rnk + 1) = meshblocks(rnk + 1) % sz
          end do
        case default
          call throwError('ERROR: unrecognized `dom_vars`: `'//trim(dom_vars(d))//'`')
        end select

        call h5screate_simple_f(datarank, data_dims, dspace_id, error)
        call h5dcreate_f(file_id, dom_vars(d), H5T_NATIVE_INTEGER, dspace_id, &
                       & dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, domain_data, data_dims, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)
      end do

      ! Close the file
      call h5fclose_f(file_id, error)
      ! Close FORTRAN interface
      call h5close_f(error)

      if (allocated(domain_data)) deallocate (domain_data)
    end if
  end subroutine writeDomain_hdf5

end module m_writeoutput
