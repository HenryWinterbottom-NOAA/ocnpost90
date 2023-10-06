program ocnpost

  use utils_mod
  use ocnvars
  use netcdf

  implicit none

  character(len=240) :: filesrc, filedst, wgtsfile, fout
  ! # TODO: Hardcoded to use pre-computed weights; CDO does not
  ! # require this and the weights are computed on the fly; further,
  ! # the ESMF weights do not account for bathymetry, only the
  ! # geo-coordinate map and there for the sub-surface ocean variables
  ! # will be incorrect as the bathymetry changes in the vertical.
  character(len=120) :: wgtsdir = '/scratch1/NCEPDEV/climate/climpara/S2S/FIX/fix_UFSp6/fix_reg2grb2/'

  ! # TODO: Hardcorded to assumed a MOM6 or CICE tripolar grid; also
  ! hardcoded to assume 40-vertical levels; the CDO method cares
  ! nothing about the `nxt` and `nyt` nor does it care about the
  ! number of vertical levels; if this were to be moved to a namelist
  ! type file it would still require more code than a single command
  ! provided to CDO.  source grid, tripole 1/4 deg, 40 vertical levels
  integer, parameter :: nxt = 1440, nyt = 1080, nlevs = 40

  ! destination grids
  integer, parameter :: ndest = 3
  ! # TODO: Hardcoded to support only specific grid configurations;
  ! # CDO is agnostic.
  integer, parameter, dimension(ndest) :: nxrs = (/1440, 720, 360/)
  integer, parameter, dimension(ndest) :: nyrs = (/ 721, 361, 181/)
  character(len=4), dimension(ndest)   :: dstgrds = (/'0p25', '0p5 ', '1p0 '/)

  ! packed arrays on source grid
  real, allocatable, dimension(:,:)   :: bilin2d, consd2d
  real, allocatable, dimension(:,:,:) :: bilin3d
  ! variable types
  type(vardefs), allocatable :: b2d(:)
  type(vardefs), allocatable :: c2d(:)
  type(vardefs), allocatable :: b3d(:)

  ! source grid fields
  real, dimension(nxt,nyt,nlevs) :: tmp3d
  real, dimension(nxt*nyt,nlevs) :: mask3d
  real, dimension(nxt*nyt)       :: cosrot, sinrot
  real, dimension(nlevs)         :: z_l
  real, dimension(0:nlevs)       :: z_i

  ! destination grid fields, vectors
  real, allocatable, dimension(:,:)  :: dstlon, dstlat
  real, allocatable, dimension(:,:)  :: rgmask3d
  ! output fields (x,y,z)
  real, allocatable, dimension(:)     :: out1d
  real, allocatable, dimension(:,:)   :: out2d
  real, allocatable, dimension(:,:,:) :: out3d
  ! packed remapped fields
  real, allocatable, dimension(:,:)   :: rgb2d, rgc2d
  real, allocatable, dimension(:,:,:) :: rgb3d

  real(kind=8) :: timestamp
  character(len= 40) :: timeunit, timecal
  character(len= 20) :: vname, vunit
  character(len=120) :: vlong
  character(len=4)   :: dstgrid

  integer :: nd, nxr, nyr
  integer :: i,j,k,n,nn,nvalid
  integer :: rc,ncid,varid,dimid
  integer :: nbilin2d,nbilin3d,nconsd2d
  real    :: vfill
  integer :: idimid,jdimid,kdimid,edimid,timid
  integer :: idx1,idx2,idx3

  logical :: debug = .true.

  ! initialize the variable type
  call ovars_typedefine

  ! --------------------------------------------------------
  ! read the source file and obtain the units and long name,
  ! rotation angles, vertical grid and time axis
  ! --------------------------------------------------------

  filesrc = '/scratch1/NCEPDEV/climate/Jiande.Wang/For-others/For-Eric/ocn2013100106.01.2013100100.nc'
  !filesrc='ocn_2013_04_01_03.nc'

  nvalid = 0
  rc = nf90_open(trim(filesrc), nf90_nowrite, ncid)
  do i = 1,maxvars
     if (len_trim(ovars(i)%input_var_name) > 0 ) then
        rc = nf90_inq_varid(ncid, trim(ovars(i)%input_var_name), varid)
        rc = nf90_get_att(ncid, varid, 'long_name', ovars(i)%long_name)
        rc = nf90_get_att(ncid, varid, 'units', ovars(i)%units)
        rc = nf90_get_att(ncid, varid, '_FillValue', ovars(i)%var_fillvalue)
        nvalid = nvalid+1
        if (trim(ovars(i)%input_var_name) == 'temp')vfill = ovars(i)%var_fillvalue
     end if
  end do

  ! timestamp
  rc = nf90_inq_varid(ncid, 'time', varid)
  rc = nf90_get_var(ncid, varid, timestamp)
  rc = nf90_get_att(ncid, varid, 'units', timeunit)
  rc = nf90_get_att(ncid, varid, 'calendar', timecal)
  ! cell centers
  rc = nf90_inq_varid(ncid, 'z_l', varid)
  rc = nf90_get_var(ncid, varid, z_l)
  ! cell edges
  rc = nf90_inq_varid(ncid, 'z_i', varid)
  rc = nf90_get_var(ncid, varid, z_i)
  rc = nf90_close(ncid)
  ! rotation angles

  ! # TODO: What happens if the rotation angle variable names are
  ! # different between input files?; for example, the MOM6/CICE
  ! # tripolar grids use `angle` as their rotation which results in a
  ! # different computation for the rotation of the vectors (the CDO
  ! # script takes care of this) while the forecast output assumes
  ! # `cos_rot` and `sin_rot` (the CDO script knows how to
  ! # differentiate and is once again agnostic to the input file.
  call getfield(trim(filesrc), 'cos_rot', dims=(/nxt,nyt/), field=cosrot)
  call getfield(trim(filesrc), 'sin_rot', dims=(/nxt,nyt/), field=sinrot)

  ! --------------------------------------------------------
  ! mask3d contain 1's on land and 0's at valid points.
  ! when remapped, any mask value > 0 identifies land values that
  ! have crept into the field. remapped model fields are then
  ! masked with this interpolation mask
  ! --------------------------------------------------------

  rc = nf90_open(trim(filesrc), nf90_nowrite, ncid)
  ! 3D temp to use as mask, obtain directly from file to preserve vfill
  ! # TODO: What happens if the variable name changes?
  rc = nf90_inq_varid(ncid, 'temp', varid)
  rc = nf90_get_var(ncid, varid, tmp3d)
  rc = nf90_close(ncid)

  mask3d = reshape(tmp3d, (/nxt*nyt,nlevs/))
  ! set mask to 0 on ocean, 1 on land on source grid
  where(mask3d .eq. vfill)mask3d = 1.0
  where(mask3d .ne.   1.0)mask3d = 0.0
  if (debug) print '(a,2g14.7)',' mask min/max on source grid ',minval(mask3d),maxval(mask3d)
  if (debug) call dumpnc('mask3d.nc', 'mask3d', dims=(/nxt,nyt,nlevs/), field=mask3d)

  ! --------------------------------------------------------
  ! count numbers of fields to remapped for each
  ! mapping type; these can be remapped as packed arrays
  ! --------------------------------------------------------

  nbilin2d = 0; nbilin3d = 0; nconsd2d = 0
  do n = 1,nvalid
     if (trim(ovars(n)%var_remapmethod)  == 'bilinear') then
        if (ovars(n)%var_dimen == 2) nbilin2d = nbilin2d + 1
        if (ovars(n)%var_dimen == 3) nbilin3d = nbilin3d + 1
     end if
     if (trim(ovars(n)%var_remapmethod)  == 'conserve')nconsd2d = nconsd2d + 1  !no 3d variables w/ conservative mapping
  end do
  if (debug) print '(3(a,i4))','bilin 2d ',nbilin2d,' bilin 3d ',nbilin3d,' conserv 2d ',nconsd2d

  ! initialization required when compile in sinit_arrays=nan
  allocate(bilin2d(nxt*nyt,nbilin2d)); bilin2d = 0.0
  allocate(consd2d(nxt*nyt,nconsd2d)); consd2d = 0.0
  allocate(bilin3d(nxt*nyt,nlevs,nbilin3d)); bilin3d = 0.0

  allocate(b2d(1:nbilin2d))
  allocate(c2d(1:nconsd2d))
  allocate(b3d(1:nbilin3d))

  ! --------------------------------------------------------
  ! create types for each packed array
  ! --------------------------------------------------------

  ! # TODO: Only supports two types of remapping; CDO supports all of
  ! # these and more (e.g., nearest-neighbor) and computes the weights
  ! # for the respective method on the fly.
  i = 0; j = 0; k = 0
  do n = 1,nvalid
     if (trim(ovars(n)%var_remapmethod)  == 'bilinear') then
        if (ovars(n)%var_dimen == 2) then
           i = i+1; b2d(i) = ovars(n)
        end if
        if (ovars(n)%var_dimen == 3) then
           j = j+1; b3d(j) = ovars(n)
        end if
     end if
     if (trim(ovars(n)%var_remapmethod)  == 'conserve') then
        k = k+1; c2d(k) = ovars(n)
     end if
  end do

  if (debug) then
     print *,'2D fields mapped bilin'
     do n = 1,nbilin2d
        print '(i6,4(a,a))',n,'  ',trim(b2d(n)%input_var_name),'  ',trim(b2d(n)%var_grid), &
             '  ',trim(b2d(n)%var_pair),'  ', trim(b2d(n)%var_pair_grid)
     end do
     print *,'3D fields mapped bilin'
     do n = 1,nbilin3d
        print '(i6,4(a,a))',n,'  ',trim(b3d(n)%input_var_name),'  ',trim(b3d(n)%var_grid), &
             '  ',trim(b3d(n)%var_pair),'  ', trim(b3d(n)%var_pair_grid)
     end do
     print *,'2D fields mapped conserv'
     do n = 1,nconsd2d
        print '(i6,4(a,a))',n,'  ',trim(c2d(n)%input_var_name),'  ',trim(c2d(n)%var_grid), &
             '  ',trim(c2d(n)%var_pair),'  ', trim(c2d(n)%var_pair_grid)
     end do
  end if

  ! --------------------------------------------------------
  ! create packed arrays for mapping
  ! --------------------------------------------------------

  ! 2D bilin
  call packarrays(trim(filesrc), trim(wgtsdir), cosrot, sinrot, b2d, dims=(/nxt,nyt/), nflds=nbilin2d, fields=bilin2d)
  if (debug) call dumpnc('bilin2d.nc', 'bilin2d', dims=(/nxt,nyt/), nflds=nbilin2d, field=bilin2d)

  ! 2D conserv
  call packarrays(trim(filesrc), trim(wgtsdir), cosrot, sinrot, c2d, dims=(/nxt,nyt/), nflds=nconsd2d, fields=consd2d)
  if (debug) call dumpnc('consd2d.nc', 'consd2d', dims=(/nxt,nyt/), nflds=nconsd2d, field=consd2d)

  ! 3D bilin
  call packarrays(trim(filesrc), trim(wgtsdir), cosrot, sinrot, b3d, dims=(/nxt,nyt,nlevs/), nflds=nbilin3d, fields=bilin3d)
  if (debug) call dumpnc('bilin3d.nc', 'bilin3d', dims=(/nxt,nyt,nlevs/), nk=nlevs, nflds=nbilin3d, field=bilin3d)

  ! --------------------------------------------------------
  ! remap packed arrays to each destination grid
  ! --------------------------------------------------------

  do nd = 1,ndest
     dstgrid =  trim(dstgrds(nd))
     nxr = nxrs(nd); nyr = nyrs(nd)

     allocate(rgb2d(nxr*nyr,nbilin2d)); rgb2d = 0.0
     allocate(rgc2d(nxr*nyr,nconsd2d)); rgc2d = 0.0
     allocate(rgb3d(nxr*nyr,nlevs,nbilin3d)); rgb3d = 0.0
     allocate(dstlon(nxr,nyr)); dstlon = 0.0
     allocate(dstlat(nxr,nyr)); dstlat = 0.0
     allocate(out1d(nxr*nyr)); out1d = 0.0
     allocate(out2d(nxr,nyr)); out2d = 0.0
     allocate(out3d(nxr,nyr,nlevs)); out3d = 0.0
     allocate(rgmask3d(nxr*nyr,nlevs)); rgmask3d = 0.0

     ! lat,lon of destination grid can be obtained from xc_b,yc_b in wgtsfile
     ! # TODO: What happens if the variable name changes?
     wgtsfile = trim(wgtsdir)//'tripole.mx025.Ct.to.rect.'//trim(dstgrid)//'.bilinear.nc'
     rc = nf90_open(trim(wgtsfile), nf90_nowrite, ncid)
     rc = nf90_inq_varid(ncid, 'xc_b', varid)
     rc = nf90_get_var(ncid,    varid, out1d)
     dstlon = reshape(out1d,(/nxr,nyr/))
     rc = nf90_inq_varid(ncid, 'yc_b', varid)
     rc = nf90_get_var(ncid,    varid, out1d)
     dstlat = reshape(out1d,(/nxr,nyr/))
     rc = nf90_close(ncid)

     ! # TODO: Hardcoded paths; the path to the ESMF remapping
     ! # coefficients is not required by CDO.
     wgtsfile = trim(wgtsdir)//'tripole.mx025.Ct.to.rect.'//trim(dstgrid)//'.bilinear.nc'
     if (debug) print '(a)','remapping 2D fields bilinear with '//trim(wgtsfile)
     call remap(trim(wgtsfile), dim2=nbilin2d, src_field=bilin2d, dst_field=rgb2d)
     if (debug) call dumpnc('rgbilin2d.'//dstgrid//'.nc', 'rgbilin2d', dims=(/nxr,nyr/), nflds=nbilin2d, field=rgb2d)

     ! # TODO: Hardcoded paths; the path to the ESMF remapping
     ! # coefficients is not required by CDO.
     wgtsfile = trim(wgtsdir)//'tripole.mx025.Ct.to.rect.'//trim(dstgrid)//'.conserve.nc'
     if (debug) print '(a)','remapping 2D fields conserv with '//trim(wgtsfile)
     call remap(trim(wgtsfile), dim2=nconsd2d, src_field=consd2d, dst_field=rgc2d)
     if (debug) call dumpnc('rgconsd2d.'//dstgrid//'.nc', 'rgconsd2d', dims=(/nxr,nyr/), nflds=nconsd2d, field=rgc2d)

     ! # TODO: Hardcoded paths; the path to the ESMF remapping
     ! # coefficients is not required by CDO.
     wgtsfile = trim(wgtsdir)//'tripole.mx025.Ct.to.rect.'//trim(dstgrid)//'.bilinear.nc'
     if (debug) print '(a)','remapping 3D fields bilinear with '//trim(wgtsfile)
     call remap(trim(wgtsfile), nk=nlevs, nflds=nbilin3d, src_field=bilin3d, dst_field=rgb3d)
     if (debug) call dumpnc('rgbilin3d.'//dstgrid//'.nc', 'rgbilin3d', dims=(/nxr,nyr,nlevs/), nk=nlevs, nflds=nbilin3d, field=rgb3d)

     ! --------------------------------------------------------
     ! remap the source grid 3D mask to obtain the interpolation mask.
     ! --------------------------------------------------------

     ! # TODO: Hardcoded paths; the path to the ESMF remapping
     ! # coefficients is not required by CDO.
     wgtsfile = trim(wgtsdir)//'tripole.mx025.Ct.to.rect.'//trim(dstgrid)//'.bilinear.nc'
     call remap(trim(wgtsfile), dim2=nlevs, src_field=mask3d, dst_field=rgmask3d)
     ! set interpolation mask missing on land, 1.0 on ocean on destination grids
     where(rgmask3d > 0.0)rgmask3d = vfill
     where(rgmask3d /= vfill)rgmask3d = 1.0
     if (debug) print '(a,2g14.7)',' mask min/max on destination grid ',minval(rgmask3d),maxval(rgmask3d)
     if (debug) call dumpnc('rgmask3d.nc', 'rgmask3d', dims=(/nxr,nyr,nlevs/), field=rgmask3d)

     ! --------------------------------------------------------
     ! mask the mapped fields
     ! --------------------------------------------------------

     do n = 1,nbilin2d
        where(rgmask3d(:,1) .eq. vfill)rgb2d(:,n) = vfill
     end do
     do n = 1,nconsd2d
        where(rgmask3d(:,1) .eq. vfill)rgc2d(:,n) = vfill
     end do
     do n = 1,nbilin3d
        where(rgmask3d(:,:) .eq. vfill)rgb3d(:,:,n) = vfill
     end do

     ! --------------------------------------------------------
     ! replace model native speed field with a value calculated
     ! from remapped ssu,ssv
     ! --------------------------------------------------------

     ! # TODO: What happens if the variable names change?
     do n = 1,nbilin2d
        if (trim(b2d(n)%output_var_name) == 'speed')idx1 = n
        if (trim(b2d(n)%output_var_name) ==   'SSU')idx2 = n
        if (trim(b2d(n)%output_var_name) ==   'SSV')idx3 = n
     enddo
     where(rgb2d(:,idx1) .ne. vfill)rgb2d(:,idx1) = &
          sqrt(rgb2d(:,idx2)**2 + rgb2d(:,idx3)**2)

     ! # TODO: This entire next block of FORTRAN code is done via a
     ! single CDO command.
     ! -------------------------------------------------------- write
     ! the mapped fields
     ! --------------------------------------------------------

     fout = 'test.'//dstgrid//'.nc'

     ! # TODO: What happens if the variable names change?
     rc = nf90_create(trim(fout), nf90_clobber, ncid)
     rc = nf90_def_dim(ncid, 'longitude', nxr, idimid)
     rc = nf90_def_dim(ncid, 'latitude',  nyr, jdimid)
     rc = nf90_def_dim(ncid,  'z_l',    nlevs, kdimid)
     rc = nf90_def_dim(ncid,  'z_i',  nlevs+1, edimid)
     rc = nf90_def_dim(ncid, 'time', nf90_unlimited, timid)

     ! define the time variable
     rc = nf90_def_var(ncid, 'time', nf90_double, (/timid/), varid)
     rc = nf90_put_att(ncid, varid, 'units', trim(timeunit))
     rc= nf90_put_att(ncid, varid, 'calendar', trim(timecal))
     ! spatial grid
     rc = nf90_def_var(ncid, 'longitude', nf90_float,  (/idimid/), varid)
     rc = nf90_put_att(ncid, varid, 'units', 'degrees_east')
     rc = nf90_def_var(ncid, 'latitude', nf90_float,  (/jdimid/), varid)
     rc = nf90_put_att(ncid, varid, 'units', 'degrees_north')
     ! vertical grid
     rc = nf90_def_var(ncid, 'z_l', nf90_float,  (/kdimid/), varid)
     rc = nf90_put_att(ncid, varid, 'units', 'm')
     rc = nf90_put_att(ncid, varid, 'positive', 'down')
     rc = nf90_def_var(ncid, 'z_i', nf90_float,  (/edimid/), varid)
     rc = nf90_put_att(ncid, varid, 'units', 'm')
     rc = nf90_put_att(ncid, varid, 'positive', 'down')

     do n = 1,nbilin2d
        vname = trim(b2d(n)%output_var_name)
        vunit = trim(b2d(n)%units)
        vlong = trim(b2d(n)%long_name)
        vfill = b2d(n)%var_fillvalue
        rc = nf90_def_var(ncid, vname, nf90_float, (/idimid,jdimid,timid/), varid)
        rc = nf90_put_att(ncid, varid,     'units',  vunit)
        rc = nf90_put_att(ncid, varid, 'long_name',  vlong)
        rc = nf90_put_att(ncid, varid, '_FillValue', vfill)
     enddo
     do n = 1,nconsd2d
        vname = trim(c2d(n)%output_var_name)
        vunit = trim(c2d(n)%units)
        vlong = trim(c2d(n)%long_name)
        vfill = c2d(n)%var_fillvalue
        rc = nf90_def_var(ncid, vname, nf90_float, (/idimid,jdimid,timid/), varid)
        rc = nf90_put_att(ncid, varid,     'units',  vunit)
        rc = nf90_put_att(ncid, varid, 'long_name',  vlong)
        rc = nf90_put_att(ncid, varid, '_FillValue', vfill)
     enddo
     do n = 1,nbilin3d
        vname = trim(b3d(n)%output_var_name)
        vunit = trim(b3d(n)%units)
        vlong = trim(b3d(n)%long_name)
        vfill = b3d(n)%var_fillvalue
        rc = nf90_def_var(ncid, vname, nf90_float, (/idimid,jdimid,kdimid,timid/), varid)
        rc = nf90_put_att(ncid, varid,     'units',  vunit)
        rc = nf90_put_att(ncid, varid, 'long_name',  vlong)
        rc = nf90_put_att(ncid, varid, '_FillValue', vfill)
     enddo
     rc = nf90_enddef(ncid)

     ! dimensions
     rc = nf90_inq_varid(ncid, 'z_l', varid)
     rc = nf90_put_var(ncid,   varid,   z_l)
     rc = nf90_inq_varid(ncid, 'z_i', varid)
     rc = nf90_put_var(ncid,   varid,   z_i)
     rc = nf90_inq_varid(ncid, 'longitude', varid)
     rc = nf90_put_var(ncid,   varid, dstlon(:,1))
     rc = nf90_inq_varid(ncid,  'latitude', varid)
     rc = nf90_put_var(ncid,   varid, dstlat(1,:))
     ! time
     rc = nf90_inq_varid(ncid, 'time', varid)
     rc = nf90_put_var(ncid, varid, timestamp)

     do n = 1,nbilin2d
        out2d(:,:) = reshape(rgb2d(:,n), (/nxr,nyr/))
        vname = trim(b2d(n)%output_var_name)
        rc = nf90_inq_varid(ncid, vname, varid)
        rc = nf90_put_var(ncid,   varid, out2d)
     end do
     do n = 1,nconsd2d
        out2d(:,:) = reshape(rgc2d(:,n), (/nxr,nyr/))
        vname = trim(c2d(n)%output_var_name)
        rc = nf90_inq_varid(ncid, vname, varid)
        rc = nf90_put_var(ncid,   varid, out2d)
     end do
     do n = 1,nbilin3d
        out3d(:,:,:) = reshape(rgb3d(:,:,n), (/nxr,nyr,nlevs/))
        vname = trim(b3d(n)%output_var_name)
        rc = nf90_inq_varid(ncid, vname, varid)
        rc = nf90_put_var(ncid,   varid, out3d)
     end do
     rc = nf90_close(ncid)

     deallocate(rgb2d, rgc2d, rgb3d)
     deallocate(out1d, out2d, out3d, dstlon, dstlat, rgmask3d)
  end do !nd
  print *,'all done!'

end program ocnpost
