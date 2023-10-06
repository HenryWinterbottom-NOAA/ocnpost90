module ocnvars

  implicit none

  integer, parameter :: maxvars = 40          !< The maximum number of variables written to a file

  type :: vardefs
     character(len= 10)   :: input_var_name    !< A variable's  input variable name
     character(len= 10)   :: output_var_name   !< A variable's output variable name
     character(len=120)   :: long_name         !< A variable's long name
     character(len= 20)   :: units             !< A variable's unit
     character(len= 10)   :: var_remapmethod   !< A variable's mapping method
     integer              :: var_dimen         !< A variable's dimensionality
     character(len=  2)   :: var_grid          !< A variable's input grid location; all output locations are on cell centers
     character(len= 10)   :: var_pair          !< A variable's pair
     character(len=  2)   :: var_pair_grid     !< A pair variable grid
     real                 :: var_fillvalue     !< A variable's fillvalue
  end type vardefs

  type(vardefs), public :: ovars(maxvars)   !< Attribute definitions for the variables

contains

  ! # TODO: This entire module will need to be converted to a
  ! # namelist; the CDO scripts does not assume any variable names;
  ! # further the variables to write, and the corresponding
  ! # interpolation method, is defined in a single `.csv` file and
  ! # passed to the script; this section is not acceptable as it is
  ! # hardcoded and assumes a single application.
  subroutine ovars_typedefine

    ! local variables
    integer :: ii = 0

    !set defaults
    ovars(:)%input_var_name=''
    ovars(:)%var_grid  = 'Ct'
    ovars(:)%var_remapmethod  = 'bilinear'
    ovars(:)%var_dimen = 2
    ovars(:)%var_pair = ''
    ovars(:)%var_pair_grid = ''
    ovars(:)%long_name = ''          ! obtained from input file
    ovars(:)%units = ''              ! obtained from input file
    ovars(:)%var_fillvalue = -1.0    ! obtained from input file

    ! 2D states with native grid location on cell centers; remapped bilinearly
    ii = ii + 1; ovars(ii)%input_var_name  = 'SSH'
    ii = ii + 1; ovars(ii)%input_var_name  = 'SST'
    ii = ii + 1; ovars(ii)%input_var_name  = 'SSS'
    ii = ii + 1; ovars(ii)%input_var_name  = 'speed'
    !ii = ii + 1; ovars(ii)%input_var_name  = 'mld'
    ii = ii + 1; ovars(ii)%input_var_name  = 'ePBL'
    ii = ii + 1; ovars(ii)%input_var_name  = 'MLD_003'
    ii = ii + 1; ovars(ii)%input_var_name  = 'MLD_0125'

    ! 2D fluxes with native grid location on cell centers; remapped conservatively
    ii = ii + 1; ovars(ii)%input_var_name  = 'latent'    ; ovars(ii)%var_remapmethod  = 'conserve'
    ii = ii + 1; ovars(ii)%input_var_name  = 'sensible'  ; ovars(ii)%var_remapmethod  = 'conserve'
    ii = ii + 1; ovars(ii)%input_var_name  = 'SW'        ; ovars(ii)%var_remapmethod  = 'conserve'
    ii = ii + 1; ovars(ii)%input_var_name  = 'LW'        ; ovars(ii)%var_remapmethod  = 'conserve'
    ii = ii + 1; ovars(ii)%input_var_name  = 'evap'      ; ovars(ii)%var_remapmethod  = 'conserve'
    ii = ii + 1; ovars(ii)%input_var_name  = 'lprec'     ; ovars(ii)%var_remapmethod  = 'conserve'
    ii = ii + 1; ovars(ii)%input_var_name  = 'fprec'     ; ovars(ii)%var_remapmethod  = 'conserve'
    ii = ii + 1; ovars(ii)%input_var_name  = 'LwLatSens' ; ovars(ii)%var_remapmethod  = 'conserve'
    ii = ii + 1; ovars(ii)%input_var_name  = 'Heat_PmE'  ; ovars(ii)%var_remapmethod  = 'conserve'

    ! # TODO: These are concessions for a single grid type; CDO does
    ! this in a single command.  2D vector states on stagger
    ! locations; remapped bilinearly
    ii = ii + 1; ovars(ii)%input_var_name  = 'SSU'
                 ovars(ii)%var_grid = 'Cu'
                 ovars(ii)%var_pair = 'SSV'
                 ovars(ii)%var_pair_grid = 'Cv'

    ii = ii + 1; ovars(ii)%input_var_name  = 'SSV'
                 ovars(ii)%var_grid = 'Cv'
                 ovars(ii)%var_pair = 'SSU'
                 ovars(ii)%var_pair_grid = 'Cu'

    ! 2D vector fluxes on stagger locations; remapped conservatively
    ii = ii + 1; ovars(ii)%input_var_name = 'taux'
                 ovars(ii)%var_grid = 'Cu'
                 ovars(ii)%var_pair = 'tauy'
                 ovars(ii)%var_pair_grid = 'Cv'
                 ovars(ii)%var_remapmethod = 'conserve'

    ii = ii + 1; ovars(ii)%input_var_name = 'tauy'
                 ovars(ii)%var_grid = 'Cv'
                 ovars(ii)%var_pair = 'taux'
                 ovars(ii)%var_pair_grid = 'Cu'
                 ovars(ii)%var_remapmethod = 'conserve'

    ! 3D scalars with native grid location on cell centers; remapped bilinearly
    ii = ii + 1; ovars(ii)%input_var_name  = 'temp' ; ovars(ii)%var_dimen = 3
    ii = ii + 1; ovars(ii)%input_var_name  = 'so'   ; ovars(ii)%var_dimen = 3

    ! 3D vectors on stagger locations; remapped bilinearly
    ii = ii + 1; ovars(ii)%input_var_name  = 'uo'
                 ovars(ii)%var_grid = 'Cu'
                 ovars(ii)%var_pair = 'vo'
                 ovars(ii)%var_pair_grid = 'Cv'
                 ovars(ii)%var_dimen = 3

    ii = ii + 1; ovars(ii)%input_var_name  = 'vo'
                 ovars(ii)%var_grid = 'Cv'
                 ovars(ii)%var_pair = 'uo'
                 ovars(ii)%var_pair_grid = 'Cu'
                 ovars(ii)%var_dimen = 3

    ! set default output name
    ovars(:)%output_var_name = ovars(:)%input_var_name

  end subroutine ovars_typedefine
end module ocnvars
