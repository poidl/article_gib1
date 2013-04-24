/*******+*********+*********+*********+*********+*********+*********+*
 *    This header file determies many of the adjustable parameters    *
 *  for the Hallberg Isopycnal Model (HIM).  Where appropriate, MKS   *
 *  units are used.  The program param_suggest.c will analyze this    *
 *  header to suggest appropriate values for DT, DTBT, and AH.        *
 *                                                                    *
 * @  Values marked with @ at the start of the describing comment can *
 *  be reset to different values at run time via the file given to    *
 *  HIM_parser.c.                                                     *
 *                                                                    *
 ********+*********+*********+*********+*********+*********+*********+*/

/*  Specify the physical domain.                                      */
#define OMEGA 7.2921e-5          /*    The rotation rate of the earth   */
                               /*  in s-1.                            */
#define RE 6.378e6             /*    The radius of the earth in m.    */

#define GFS 9.8                /*    The reduced gravity at the free  */
                               /*  surface, in m s-2.                 */
#define MAXIMUM_DEPTH 1500.0   /*    The maximum depth of the ocean,  */
                               /*  in m.  This is required for valid  */
                               /*  advice from param_suggest.c.       */
#define MINIMUM_DEPTH 1.0      /*    The minimum ocean depth, in m.   */
                               /*  Anything shallower than this depth */
                               /*  is assumed to be on land, and all  */
                               /*  appropriate fluxes are masked out. */
#define GINT 0.0190            /*  The reduced gravity of the       */
                               /*  internal interfaces, in m s-2.     */
                               /*  GINT may be superseded in          */
                               /*  initialize.c, but GINT should be   */
                               /*  the average of all interface re-   */
                               /*  duced gravities for param_suggest.c*/
                               /*  to give appropriate advice.        */

#define LENLAT 170000.0         /*    The length of the domain in units*/
#define LENLON 1033970.72      /*  that are defined by AXIS_UNITS; by */
                               /*  default the units are degrees of   */
                               /*  latitude and longitude.            */
#define LOWLAT -84000.0        /*  The domain's southern latitude.    */
#define WESTLON -92997.366     /*  The domain's western longitude.    */
#define  AXIS_UNITS 'm'         /*    AXIS_UNITS should be defined as  */
                               /*  'k' for the axis units to be km or */
                               /*  'm' for the axis units to be m.    */
                               /*  Other units may be defined later,  */
                               /*  but the default is degrees of      */
                               /*  latitude and longitude.  Except on */
                               /*  a Cartesian grid, only the default */
                               /*  is currently implemented.          */

#undef REENTRANT_X             /*    If the domain is zonally         */
                               /*  reentrant this should be '#define  */
                               /*  REENTRANT_X'.  Otherwise, use      */
                               /*  '#undef REENTRANT_X'.              */
#undef REENTRANT_Y             /*    If the domain is meridionally    */
                               /*  reentrant this should be '#define  */
                               /*  REENTRANT_Y'.  Otherwise, use      */
                               /*  '#undef REENTRANT_Y'.              */
#undef  TRIPOLAR_N             /*    Use tripolar connectivity at the */
                               /*  northern edge of the domain.  With */
                               /*  TRIPOLAR_N, NXTOT must be even.    */

#define  CARTESIAN             /*    Use a uniform Cartesian grid if  */
                               /*  CARTESIAN is defined, and a spher- */
                               /*  ical coordinate grid with uniform  */
                               /*  spacing in latitude and longitude  */
                               /*  otherwise.                         */
#undef ISOTROPIC               /*    If ISOTROPIC is defined, an      */
                               /*  isotropic grid on a sphere (also   */
                               /*  known as a Mercator grid) is used. */
                               /*  With an isotropic grid, the merid- */
                               /*  ional extent of the domain         */
                               /*  (LENLAT), the zonal extent         */
                               /*  (LENLON), and the number of grid   */
                               /*  points in each direction are _not_ */
                               /*  independent.  Here the meridional  */
                               /*  extent will be determined to fit   */
                               /*  the zonal extent and the number of */
                               /*  grid points.  The grid is          */
                               /*  perfectly isotropic.               */
#undef  EQUATOR_REFERENCE      /*    If EQUATOR_REFERENCE is defined, */
                               /*  the grid is defined to have the    */
                               /*  equator at the nearest q or h grid */
                               /*  point to (-LOWLAT*NYTOT/LENLAT).   */
#undef LAT_EQ_ENHANCE 0.0      /*   The latitude (north and south) to */
                               /*  which the resolution is enhanced.  */
#undef LAT_ENHANCE_FACTOR 1.0  /*   The amount by which the meridional*/
                               /*  resolution is enhanced within      */
                               /*  LAT_EQ_ENHANCE of the equator.     */
#undef  XMETRIC_J              /*    Define XMETRIC_J if the x-direc- */
#undef  XMETRIC_I              /*  tion metrics vary in the y- (or j^)*/
#undef  YMETRIC_J              /*  direction.  Otherwise undefine     */
#undef  YMETRIC_I              /*  XMETRIC_J.  XMETRIC_I, YMETRIC_J,  */
                               /*  and YMETRIC_I are used similarly.  */
                               /*  CARTESIAN overrides all of these   */
                               /*  choices.                           */
                               /*    For example, on a regular lat-   */
                               /*  itude longitude grid, define       */
                               /*  XMETRIC_J and undefine the rest.   */
                               /*  For a Mercator grid, define        */
                               /*  only XMETRIC_J and YMETRIC_J.      */

/*  Specify the numerical domain.  */
#define NXTOT 517               /*    NXTOT and NYTOT are the number   */
#define NYTOT 85               /*  of thickness grid points in the    */
                               /*  zonal and meridional directions of */
                               /*  the physical domain.               */
#define NZ 2                   /*    The number of layers.            */

#ifndef SERIAL_IO_CODE         /*    SERIAL_IO_CODE is specified on   */
                               /*  the compile line of postprocessing */
                               /*  code.                              */
# define MPI_PARALLEL          /*    If MPI_PARALLEL is defined, the  */
                               /*  parallelization is done using the  */
                               /*  MPI subroutines.  Otherwise SHMEM  */
                               /*  subroutines are used.  This option */
                               /*  affects only file HIM_parallel.c.  */
# define  PARALLEL_X            /*    Define PARALLEL_X if there is to */
                               /*  be domain decomposition in the X-  */
                               /*  direction.                         */
# define  PARALLEL_Y            /*    Define PARALLEL_Y if there is to */
                               /*  be domain decomposition in the Y-  */
                               /*  direction.                         */
# undef  PARALLEL_IO           /*  With PARALLEL_IO and NETCDF_OUTPUT */
                               /*  defined, each processor writes out */
                               /*  its own NetCDF file.  These files  */
                               /*  can be combined using the utility  */
                               /*  mppnccombine.                      */
#endif

#ifdef PARALLEL_X
# define NXPROC 4              /*    NXPROC is the minimum number of  */
                               /*  processors in the x-direction.     */
#endif
#ifdef PARALLEL_Y
# define NYPROC 2              /*    NYPROC is the minimum number of  */
                               /*  processors in the y-direction.     */
                               /*  The minimum total number of        */
                               /*  processors is NXPROC*NYPROC.       */
#endif
#define MAXPROC 8              /*    MAXPROC is the maximum number of */
                               /*  processors that might be used      */
                               /*  without recompiling.  MAXPROC must */
                               /*  exceed NXPROC*NYPROC or NXPROC or  */
                               /*  NYPROC, or be undefined.           */

#undef  CHECKPARALLEL          /*    If CHECKPARALLEL is defined, it  */
                               /*  causes the model to run simultan-  */
                               /*  eously on one and several process- */
                               /*  ors.  The subroutine check_field   */
                               /*  can then be called to compare the  */
                               /*  fields in the serial and parallel  */
                               /*  simulations, reporting any diff-   */
                               /*  erences.                           */
#undef  CHECK_TRIPOLAR         /*    If CHECK_TRIPOLAR and TRIPOLAR_N */
                               /*  are defined, 2 1-PE copies of the  */
                               /*  model are run, one on a tripolar   */
                               /*  grid, the other on a non-tripolar  */
                               /*  grid that is half as wide in X and */
                               /*  twice as long in Y. CHECK_TRIPOLAR */
                               /*  needs CHECKPARALLEL to be defined, */
                               /*  and probably won't work unless     */
                               /*  there is a strip of land at the    */
                               /*  eastern and western edges of the   */
                               /*  domain.                            */

/*  Specify the time integration scheme.                              */
#define SPLIT                  /*    Use "#define SPLIT" to use the   */
                               /*  split time stepping.               */
#define DT 75.0               /* @  The time step, in s.             */
#define DTBT 8.0                /* @  The barotropic time step, in s.  */
                               /*  DTBT is only used with the split   */
                               /*  explicit time stepping.            */
#define BE 0.7                 /* @  BE determines whether the neutral*/
                               /*  baroclinic time stepping scheme    */
                               /*  (0.5) or a backward Euler scheme   */
                               /*  (1) is used.  BE may be from 0.5   */
                               /*  to 1, but instability may occur    */
                               /*  near 0.5.                          */
#define BEBT 0.2               /* @  BEBT determines whether the baro-*/
                               /*  tropic time stepping uses the for- */
                               /*  ward-backward timestepping scheme  */
                               /*  or a backward Euler scheme.  BEBT  */
                               /*  is valid in the range from 0 (for  */
                               /*  a forward-backward treatment of    */
                               /*  nonrotating gravity waves) to 1    */
                               /*  (for a backward Euler treatment).  */
                               /*  In practice, BEBT must be greater  */
                               /*  than about 0.05.                   */

/*  Specify properties of the I/O and length of the integration.      */
#define TIMEUNIT 86400.0       /*    The time unit in seconds for the */
                               /*  following fields.                  */
#define DAYMAX 1080.0           /* @  The final day of the simulation. */
#define DAYSAVE 0.0           /* @  The first day to save output.    */
#define SAVEINT 3.0            /* @  The number of days between saves.*/
#define SAVES_PER_FILE 1000     /* @  The number time levels to be     */
                               /*  put into each file.                */
#define ENERGYSAVEDAYS 3.0     /* @  The number of days between saves */
                               /*  of the energies of the run.        */
#define ENERGYFILE "timestats" /* @  The file to use to save the      */
                               /*  energies of the run.               */
#define RESTARTFILE "RESTART"  /* @  The name of the restart file.    */
#define RESTINT 2000.0           /* @  The number of days between saves */
                               /*  of the restart file.  Use a value  */
                               /*  that is larger than DAYMAX not to  */
                               /*  save incremental restart files     */
                               /*  within a run.  Use 0 not to save   */
                               /*  restart files at all.              */
#define MAXCPU -1.0         /* @  The maximum amount of cpu time   */
                               /*  per processor for which HIM should */
                               /*  run before saving a restart file   */
                               /*  and quiting with a return value    */
                               /*  that indicates that a further exe- */
                               /*  cution is required to complete the */
                               /*  simulation.  If automatic restarts */
                               /*  are not desired, use a negative    */
                               /*  value for MAXCPU.  MAXCPU has units*/
                               /*  of wall-clock seconds (i.e. CPU    */
                               /*  time limit is larger by a factor of*/
                               /*  the number of processors used.     */
#define SAVEAVERAGES           /* @  #define SAVEAVERAGES to save any */
                               /*  time averaged fields.  This only   */
                               /*  affects the flow in USER_set_output*/
                               /*  in initialize_output.c.            */
#define SAVE_DOUBLE            /*    If SAVE_DOUBLE is defined, all   */
                               /*  fields are written in double pre-  */
                               /*  cision.  Otherwise fields may be   */
                               /*  written as floating point numbers. */
#define INPUTDIR "/home/stefan/github_final/run181/data"
                               /* @  INPUTDIR is a directory in which */
                               /*  NetCDF input files might be found. */
#define MAX_FILES 50           /*    The maximum permitted number     */
                               /*  (each) of output and restart files.*/
#define MAX_FIELDS 100         /*    The maximum permitted number     */
                               /*  (each) of output and restart       */
                               /*  variables, summed across all the   */
                               /*  files of a type.                   */

/*  Specify the horizontal (along-isopycnal) viscosity.               */
#undef LAPLACIAN              /*    LAPLACIAN is defined to use a    */
                               /*  Laplacian horizontal viscosity.    */
#define BIHARMONIC              /*    BIHARMONIC is defined to use a   */
                               /*  biharmonic horizontal viscosity.   */
                               /*  BIHARMONIC may be used with        */
                               /*  LAPLACIAN, and it is automatically */
                               /*  defined if LAPLACIAN is undefined. */
#undef BOUND_KH               /*    If BOUND_KH is defined, the      */
                               /*  Laplacian coefficient is locally   */
                               /*  limited to guarantee stability.    */
#undef BOUND_AH               /*    If BOUND_AH is defined, the bi-  */
                               /*  harmonic coefficient is locally    */
                               /*  limited to guarantee stability.    */
#define KH 0.0               /* @  KH is the Laplacian horizontal   */
                               /*  viscosity, in m2 s-1.  KH is only  */
                               /*  used if LAPLACIAN is defined.      */
#define AH 4e8                /* @  AH is the biharmonic horizontal  */
                               /*  viscosity, in m4 s-1.  AH is only  */
                               /*  used if BIHARMONIC is defined.     */
#undef KH_VEL_SCALE 0.01      /* @  The velocity scale which is mult-*/
                               /*  tiplied by the grid spacing to     */
                               /*  calculate the Laplacian viscosity  */
                               /*  if LAPLACIAN is defined, in m s-1. */
                               /*  The final viscosity is the largest */
                               /*  of this scaled viscosity, the Smag-*/
                               /*  orinsky viscosity and KH.          */
#undef AH_VEL_SCALE 0.003     /* @  The velocity scale which is mult-*/
                               /*  tiplied by the cube of the grid    */
                               /*  spacing to calculate the biharmonic*/
                               /*  viscosity if BIHARMONIC is defined,*/
                               /*  in units of m s-1. The final vis-  */
                               /*  cosity is the largest of this      */
                               /*  scaled viscosity, the Smagorinsky  */
                               /*  viscosity and AH.                  */
#undef SMAGORINSKY_KH         /*    Use Smagorinsky's nonlinear eddy */
                               /*  viscosity if SMAGORINSKY_KH is     */
                               /*  defined.  KH is the background.    */
#define SMAG_LAP_CONST 0.15    /* @  The nondimensional Laplacian     */
                               /*  Smagorinsky constant.  Often 0.15. */
#undef SMAGORINSKY_AH         /*    Use a biharmonic form of Smag-   */
                               /*  orinsky's nonlinear eddy viscosity */
                               /*  if SMAGORINSKY_AH is defined.      */
#define SMAG_BI_CONST 0.032    /* @  The nondimensional biharmonic    */
                               /*  Smagorinsky constant.  Often 0.015.*/
#undef NOSLIP                  /*    This should be '#define NOSLIP'  */
                               /*  for no slip boundary conditions    */
                               /*  or '#undef NOSLIP' for free slip   */
                               /*  boundary conditions (the default). */
                               /*    The implementation of the free   */
                               /*  slip boundary conditions on a C-   */
                               /*  grid is much cleaner than the      */
                               /*  no slip boundary conditions.  The  */
                               /*  use of free slip b.c.s is strongly */
                               /*  encouraged.  The no slip b.c.s are */
                               /*  not implemented with the biharmonic*/
                               /*  viscosity.                         */

/*  Specify the horizontal interface depth diffusion.                 */
#undef  THICKNESSDIFFUSE       /*    If THICKNESSDIFFUSE is defined,  */
                               /*  interfaces are diffused with a     */
                               /*  coefficient of KHTH.               */
#define KHTH 500.0             /* @  KHTH is the interface depth      */
                               /*  diffusivity, in m2 s-1.            */

/*  Specify the properties of the passive tracers.                    */
#define MAXTR 3                /*    The maximum number tracers (such */
                               /*  as T and S) that might be used.    */
#define KHTR 0.0               /* @  KHTR is the along-isopycnal      */
                               /*  tracer diffusivity, in m2 s-1. No  */
                               /*  KHTR is needed for numerical       */
                               /*  stability.  If this diffusion is   */
                               /*  not needed, memory is saved if KHTR*/
                               /*  is undefined.                      */

/*  Specify the scheme for the continuity equation.                   */
#define IORD 3                 /* @  If Smolarkiewicz's MPDATA is     */
                               /*  used to step the continuity        */
                               /*  equation, IORD is the number of    */
                               /*  iterations. IORD=1 is the donor    */
                               /*  cell scheme.  Results improve up to*/
                               /*  IORD ~ 4. IORD 2 or 3 are typical. */

/*  Specify the scheme for the Coriolis and momentum advection terms. */
#undef SADOURNY                /*    If SADOURNY is defined, the Cor- */
                               /*  iolis terms are parameterized with */
                               /*  Sadourny's energy conserving       */
                               /*  scheme, otherwise Arakawa & Hsu's  */
                               /*  scheme is used.  If the deform-    */
                               /*  ation radius is not resolved, Sad- */
                               /*  ourny's scheme should probably be  */
                               /*  used.                              */
#undef BOUND_CORIOLIS          /*    If BOUND_CORIOLIS is defined, the*/
                               /*  Coriolis terms at u points are     */
                               /*  bounded by the four estimates of   */
                               /*  (f+rv)v from the four neighboring v*/
                               /*  points, and similarly at v points. */
                               /*  This option would have no effect on*/
                               /*  the SADOURNY scheme if it were     */
                               /*  possible to use centered difference*/
                               /*  thickness fluxes.  In addition, if */
                               /*  SMAGORINSKY_AH is used, the        */
                               /*  biharmonic viscosity is modified to*/
                               /*  include a term that scales quad-   */
                               /*  ratically with the velocity shears.*/

/*  Specify the properties of the active tracers and Eqn of state.    */
#undef  TEMPERATURE            /*    Temperature and salinity are used*/
                               /*  as state variables if TEMPERATURE  */
                               /*  is defined.                        */
#undef  FRAZIL                 /*    If FRAZIL is defined, water      */
                               /*  freezes if it gets too cold, and   */
                               /*  the accumulated heat deficit is    */
                               /*  returned in the surface state.     */
#undef  NONLINEAR_EOS          /*    If NONLINEAR_EOS is defined,     */
                               /*  density is calculated with a full  */
                               /*  nonlinear equation of state.       */
                               /*  TEMPERATURE must be defined for    */
                               /*  NONLINEAR_EOS to work.             */
#undef P_REF 1e5            /*    P_REF is the pressure that is    */
                               /*  used for calculating the coordinate*/
                               /*  density, in Pa (=1e4 dbar).        */
#define RHO_0 1035.0           /*    RHO_0 is used in the Boussinesq  */
                               /*  approximation to calculations of   */
                               /*  pressure and pressure gradients,   */
                               /*  in units of kg m-3.                */
#define C_P 3925.0             /*    C_P is the heat capacity of sea  */
                               /*  water in J kg-1 K-1, approximated  */
                               /*  as a constant.                     */
#define G_EARTH 9.80           /*    G_EARTH is the Earth's gravi-    */
                               /*  tational acceleration, in m s-2.   */
#undef  CORRECT_DENSITY        /*    If CORRECT_DENSITY is defined,   */
                               /*  the layer densities are restored   */
                               /*  toward their target variables by   */
                               /*  the diapycnal mixing.              */

/*  Specify the properties of the diapycnal viscosity and diffusion.  */
#undef ADIABATIC              /*    There are no diapycnal mass      */
                               /*  fluxes if ADIABATIC is defined.    */
                               /*  This assumes that KD = KDML = 0.0  */
                               /*  and that there is no buoyancy      */
                               /*  forcing, but makes the model faster*/
                               /*  by elminating subroutine calls.    */
#define NTSTEP 1               /* @  The number of time steps between */
                               /*  diabatic forcing or tracer updates.*/

#undef  BULKMIXEDLAYER         /*    Use a bulk mixed layer parameter-*/
                               /*  ization if this is defined.  Layers*/
                               /*  0 through NML+1 have variable dens-*/
                               /*  ities, and there must be at least  */
                               /*  NML+2 layers if this is defined.   */
/* The following parameters only apply when BULKMIXEDLAYER is defined.*/
#define NML 1                  /*    NML is the number of sublayers   */
                               /*  within the mixed layer.  The buffer*/
                               /*  layer is not included in NML.      */
#define MSTAR 1.25             /* @  MSTAR is a non-dimensional con-  */
                               /*  stant of proportionality between   */
                               /*  the turbulent kinetic energy input */
                               /*  from the surface and the cube of   */
                               /*  the surface friction velocity.  If */
                               /*  MSTAR is not defined in this file  */
                               /*  a default value of 1.25 is used.   */
#define NSTAR 1.0              /* @  NSTAR is the portion of the buoy-*/
                               /*  ant potential energy imparted by   */
                               /*  surface fluxes that is available   */
                               /*  to drive entrainment at the base of*/
                               /*  mixed layer when that energy is    */
                               /*  positive.                          */
#define NSTAR2 0.0             /* @  NSTAR2 is the portion of any     */
                               /*  potential energy released by con-  */
                               /*  vective adjustment that is avail-  */
                               /*  able to drive entrainment at the   */
                               /*  base of the mixed layer.           */
#define PEN_SW_FRAC 0.42       /* @  PEN_SW_FRAC is the fraction of   */
                               /*  the shortwave radiation that pen-  */
                               /*  etrates below the surface.         */
#define PEN_SW_SCALE 15.0      /* @  PEN_SW_SCALE is the vertical     */
                               /*  absorption e-folding depth of the  */
                               /*  penetrating shortwave radiation,   */
                               /*  in m.                              */
#define TKE_DECAY 2.5          /* @  TKE_DECAY relates the vertical   */
                               /*  rate of decay of the TKE available */
                               /*  for mechanical entrainment to the  */
                               /*  natural Ekman depth.  Nondim.      */
#define CONV_DECAY 0.5         /* @  CONV_DECAY relates the vertical  */
                               /*  rate of decay of the convectively  */
                               /*  released TKE available for penetra-*/
                               /*  ting entrainment to the natural    */
                               /*  Ekman length.  Nondimensional.     */
/*  End of the BULKMIXEDLAYER parameters.                             */

#define HMIX 0.0               /* @  The depth of the assumed mixed   */
                               /*  layer for distribution of wind     */
                               /*  forcing, in m.  If BULKMIXEDLAYER  */
                               /*  is defined, the buoyancy fluxes    */
                               /*  are scaled away when the total     */
                               /*  depth is less than HMIX/2.         */
#define KVML 0.0               /* @  The kinematic viscosity in the   */
                               /*  mixed layer, in m2 s-1.  A typical */
                               /*  value is ~1e-2 m2 s-1. KVML is not */
                               /*  used if BULKMIXEDLAYER is defined. */
#define KDML 0.0               /* @  The diapycnal diffusivity of     */
                               /*  density in the mixed layer, in     */
                               /*  m2 s-1.  This value may be 0.0.    */
                               /*  KDML is not used if BULKMIXEDLAYER */
                               /*  is defined.                        */
#undef DIRECT_STRESS           /*    If DIRECT_STRESS is defined, the */
                               /*  wind stress is distributed over the*/
                               /*  topmost HMIX of fluid, and KVML    */
                               /*  may be set to a very small value.  */

#define KV 1.00e-6            /* @  The kinematic viscosity below    */
                               /*  the mixed layer, in m2 s-1.  The   */
                               /*  molecular value, ~1e-6 m2 s-1,     */
                               /*  might be used.                     */
#define KD 0.0                 /* @  The diapycnal diffusivity of     */
                               /*  density below the mixed layer,     */
                               /*  in m2 s-1.  Zero or the molecular  */
                               /*  value, ~1e-7 m2 s-1, may be used.  */
#define MAX_ENT_IT 5           /* @  The maximum number of iterations */
                               /*  that may be used to calculate the  */
                               /*  interior diapycnal entrainment.    */
#undef  RINOMIX                /*    Use Richardson number dependent  */
                               /*  mixing.  The mixing rate is pro-   */
                               /*  portional to the velocity shears   */
                               /*  when the shear Richardson number   */
                               /*  drops below RINO_CRIT.             */
#define RINO_CRIT 0.2          /* @  The critical shear Richardson    */
                               /*  number for shear-driven entrainment*/
                               /*  The theoretical value is 1/4, not 1*/
                               /*  as in Hallberg (MWR 2000).         */
#define MAX_RINO_IT 2          /* @  The maximum number of iterations */
                               /*  that may be used to estimate the   */
                               /*  Richardson number driven mixing.   */

#define HBBL 10.0              /* @  The thickness of a bottom        */
                               /*  boundary layer with a viscosity of */
                               /*  KVBBL if BOTTOMDRAGLAW is not      */
                               /*  defined, or the thickness over     */
                               /*  which near-bottom velocities are   */
                               /*  averaged for the drag law if       */
                               /*  BOTTOMDRAGLAW is defined but       */
                               /*  LINEAR_DRAG is not. HBBL is in m.  */
#define KVBBL 1e-3          /* @ The kinematic viscosity in the   */
                               /*  benthic boundary layer, in m2 s-1. */
                               /*  A typical value is ~1e-3 m2 s-1.   */
                               /*  KVBBL is not used with if          */
                               /*  BOTTOMDRAGLAW is defined.          */
#undef BOTTOMDRAGLAW          /*    If BOTTOMDRAGLAW is defined, the */
                               /*  bottom stress is calculated with a */
                               /*  drag law c_drag*|u|*u. The velocity*/
                               /*  magnitude may be an assumed value  */
                               /*  or it may be based on the actual   */
                               /*  velocity in the bottommost HBBL,   */
                               /*  depending on LINEAR_DRAG.          */
#define CDRAG 0.002            /* @  CDRAG is the drag coefficient    */
                               /*  relating the magnitude of the velo-*/
                               /*  city field to the bottom stress.   */
                               /*  CDRAG is only used if BOTTOMDRAGLAW*/
                               /*  is defined.                        */
#undef LINEAR_DRAG            /*    If LINEAR_DRAG and BOTTOMDRAGLAW */
                               /*  are defined, the drag law is       */
                               /*  cdrag*DRAG_BG_VEL*u.               */
#define DRAG_BG_VEL 0.10       /* @  DRAG_BG_VEL is either the assumed*/
                               /*  bottom velocity (with LINEAR_DRAG) */
                               /*  or an unresolved velocity that is  */
                               /*  combined with the resolved velocity*/
                               /*  to estimate the velocity magnitude,*/
                               /*  in m s-1.  DRAG_BG_VEL is only used*/
                               /*  when BOTTOMDRAGLAW is defined.     */
#define BBL_THICK_MIN 0.1      /* @  The minimum bottom boundary layer*/
                               /*  thickness, in m, that can be used  */
                               /*  with BOTTOMDRAGLAW.  This might be */
                               /*  Kv / (cdrag * drag_bg_vel) to give */
                               /*  Kv as the minimum near-bottom      */
                               /*  viscosity.                         */

/*  Specify the properties of surface forcing.                        */
#undef VARIABLE_WINDS          /*    If the wind stresses vary with   */
                               /*  time, define VARIABLE_WINDS, which */
                               /*  will cause wind_forcing to be      */
                               /*  called.  (wind_forcing will most   */
                               /*  likely have to be modified.)       */
#undef RESTOREBUOY             /*    If RESTOREBUOY is defined, the   */
                               /*  buoyancy fluxes drive the model    */
                               /*  back toward some observed state    */
                               /*  with a timescale determined from   */
                               /*  FLUXCONST.  Otherwise a specified  */
                               /*  buoyancy flux is used.             */
#undef FLUXCONST (50.0/2.592e6) /* A constant that relates the       */
                               /* surface fluxes to the mixed layer   */
                               /* property anomalies that is used if  */
                               /* RESTOREBUOY is defined. In m s-1.   */
#undef VARIABLE_BUOYFORCE      /*    If VARIABLE_BUOYFORCE is defined */
                               /*  a subroutine (buoyancy_forcing)    */
                               /*  that calculates temporally varying */
                               /*  surface fluxes of buoyancy, heat,  */
                               /*  and fresh water is called.         */
#ifdef RESTOREBUOY             /*    VARIABLE_BUOYFORCE must be       */
# define VARIABLE_BUOYFORCE    /*  defined if RESTOREBUOY is defined. */
#endif
#ifdef VARIABLE_BUOYFORCE      /*    If VARIABLE_BUOYFORCE is defined,*/
# undef  ADIABATIC             /*  ADIABATIC must be undefined.       */
#endif

/*  Specify whether sponges are used.                                 */
/*    It is possible to use the model in robust diagnostic mode by    */
/*  applying sponges that span the entire domain.                     */
#define SPONGE                 /*    If SPONGE is defined, sponges    */
                               /*  may be applied anywhere in the     */
                               /*  domain. The exact location and     */
                               /*  properties of those sponges are    */
                               /*  specified from initialize.c.       */

/* Specify a few miscellaneous limits.                                */
#define MAXVEL 5.            /* @  This is the maximum velocity     */
                               /*  allowed before the velocity is     */
                               /*  truncated, in units of m s-1.      */
#define MAXTRUNC 10            /* @  The run will be stopped, and the */
                               /*  day set to 9.9e9 if the velocity   */
                               /*  is truncated more than MAXTRUNC    */
                               /*  times between energy saves.  Set   */
                               /*  MAXTRUNC to 0 to stop if there is  */
                               /*  any truncation of velocities.      */
#define EPSILON 1.0e-10        /*   The minimum layer thickness, in m.*/



/*  The following lines should not be modified.                       */
#define X1 2                   /*    X1 and Y1 are the lowest indices */
#define Y1 2                   /*  of the physical domain on each PE  */
                               /*  and the size of the memory halos.  */
#define X0 (X1-1)              /*    X0 and Y0 are the offsets between*/
#define Y0 (Y1-1)              /*  the physical and memory domains.   */
#if defined(PARALLEL_X) && !defined(CHECKPARALLEL)
# define NXMEM (((NXTOT-1)/NXPROC)+1+2*X1) /* NXMEM is the size in    */
#else                          /*  memory of 2- and 3-d arrays in the */
# define NXMEM (NXTOT+2*X1)    /*  the x-direction on each processor. */
#endif

#if defined(PARALLEL_Y) && !defined(CHECKPARALLEL)
# define NYMEM (((NYTOT-1)/NYPROC)+1+2*Y1) /* NYMEM is the size in    */
#else                          /*  memory of 2- and 3-d arrays in the */
# define NYMEM (NYTOT+2*Y1)    /*  y-direction on each processor.     */
#endif

#define GLOBAL_2D_SIZE ((NXTOT+1)*(NYTOT+1))

#ifndef NXPROC
# define NXPROC 1
#endif
#ifndef NYPROC
# define NYPROC 1
#endif
#ifndef MAXPROC
# define MAXPROC ((NXPROC)*(NYPROC))
#endif

#if defined(CHECK_TRIPOLAR) && defined(TRIPOLAR_N)
# ifndef CHECKPARALLEL
#  error "CHECKPARALLEL must be defined to use CHECK_TRIPOLAR"
# endif
# undef NYMEM
# define NYMEM (2*NYTOT+2*Y1)
# undef GLOBAL_2D_SIZE
# define GLOBAL_2D_SIZE (2*(NXTOT+1)*(NYTOT+1))
#endif
