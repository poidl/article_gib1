/**********************************************************************
 *                   GNU General Public License                       *
 * This file is a part of HIM.                                        *
 *                                                                    *
 * HIM is free software; you can redistribute it and/or modify it and *
 * are expected to follow the terms of the GNU General Public License *
 * as published by the Free Software Foundation; either version 2 of  *
 * the License, or (at your option) any later version.                *
 *                                                                    *
 * HIM is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public   *
 * License for more details.                                          *
 *                                                                    *
 * For the full text of the GNU General Public License,               *
 * write to: Free Software Foundation, Inc.,                          *
 *           675 Mass Ave, Cambridge, MA 02139, USA.                  *
 * or see:   http://www.gnu.org/licenses/gpl.html                     *
 **********************************************************************/

/********+*********+*********+*********+*********+*********+*********+*
 *                                                                    *
 *  By Robert Hallberg, April 1994 - June 2002                        *
 *                                                                    *
 *    This subroutine initializes the fields for the simulations.     *
 *  The one argument passed to initialize, *day, is set to the        *
 *  current day of the simulation.  The fields which are initialized  *
 *  here are:                                                         *
 *    u - Zonal velocity in m s-1.                                    *
 *    v - Meridional velocity in m s-1.                               *
 *    h - Layer thickness in m.  (Must be positive.)                  *
 *    D - Basin depth in m.  (Must be positive.)                      *
 *    f - The Coriolis parameter, in s-1.                             *
 *    g - The reduced gravity at each interface, in m s-2.            *
 *    Rlay - Layer potential density (coordinate variable) in kg m-3. *
 *  If TEMPERATURE is defined:                                        *
 *    T - Temperature in C.                                           *
 *    S - Salinity in psu.                                            *
 *  If BULKMIXEDLAYER is defined:                                     *
 *    Rml - Mixed layer and buffer layer potential densities in       *
 *          units of kg m-3.                                          *
 *  If SPONGE is defined:                                             *
 *    A series of subroutine calls are made to set up the damping     *
 *    rates and reference profiles for all variables that are damped  *
 *    in the sponge.                                                  *
 *  Any user provided tracer code is also first linked through this   *
 *  subroutine.                                                       *
 *                                                                    *
 *    Forcing-related fields (taux, tauy, buoy, ustar, etc.) are set  *
 *  in surface_forcing.c.                                             *
 *                                                                    *
 *    These variables are all set in the set of subroutines (in this  *
 *  file) USER_initialize_bottom_depth, USER_initialize_thickness,    *
 *  USER_initialize_velocity,  USER_initialize_temperature_salinity,  *
 *  USER_initialize_mixed_layer_density, USER_initialize_sponges,     *
 *  USER_set_coordinate, and USER_set_ref_profile.                    *
 *                                                                    *
 *    The names of these subroutines should be self-explanatory. They *
 *  start with "USER_" to indicate that they will likely have to be   *
 *  modified for each simulation to set the initial conditions and    *
 *  boundary conditions.  Most of these take two arguments: an integer*
 *  argument specifying whether the fields are to be calculated       *
 *  internally or read from a NetCDF file; and a string giving the    *
 *  path to that file.  If the field is initialized internally, the   *
 *  path is ignored.                                                  *
 *                                                                    *
 *  Variables written all in capital letters are defined in init.h    *
 *                                                                    *
 *     A small fragment of the grid is shown below:                   *
 *                                                                    *
 *    j+1  x ^ x ^ x   At x:  q, f                                    *
 *    j+1  > o > o >   At ^:  v, tauy                                 *
 *    j    x ^ x ^ x   At >:  u, taux                                 *
 *    j    > o > o >   At o:  h, D, buoy, tr, T, S, Rml, ustar        *
 *    j-1  x ^ x ^ x                                                  *
 *        i-1  i  i+1  At x & ^:                                      *
 *           i  i+1    At > & o:                                      *
 *                                                                    *
 *  The boundaries always run through q grid points (x).              *
 *                                                                    *
 ********+*********+*********+*********+*********+*********+*********+*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <init.h>
#include <HIM_io.h>
#include <HIM.h>
#include <metrics.h>

extern double u[3][NZ][NYMEM][NXMEM];/* Zonal velocity, in m s-1.     */
extern double v[3][NZ][NYMEM][NXMEM];/* Meridional velocity, in m s-1.*/
extern double D[NYMEM][NXMEM];       /* Basin depth, in m.            */
extern double h[3][NZ][NYMEM][NXMEM];/* Layer thickness, in m.        */

extern struct params HIM_params;   /* HIM_params is a structure that  */
                                   /* contains a number of parameters */
                                   /* of the simulation which can be  */
                                   /* modified from a namelist file   */
                                   /* at run time.                    */

extern time_type Start_time;       /* The simulation's starting time. */

#ifdef BULKMIXEDLAYER
extern double Rml[NML+1][NYMEM][NXMEM]; /* The mixed and buffer layer */
                                     /* potential densities in kg m-2.*/
#endif
#ifdef TEMPERATURE
extern double T[NZ][NYMEM][NXMEM];   /* Potential temperature in C.   */
extern double S[NZ][NYMEM][NXMEM];   /* Salinity in PSU.              */
#endif


extern double f[NYMEM][NXMEM];       /* The Coriolis parameter in s-1.*/
extern double geolath[NYMEM][NXMEM]; /* The geographic latitude at h, */
extern double geolatq[NYMEM][NXMEM]; /* q, u, or v points in degrees  */
extern double geolatu[NYMEM][NXMEM]; /* of latitude or m.             */
extern double geolatv[NYMEM][NXMEM];
extern double geolonh[NYMEM][NXMEM]; /* The geographic longitude at h,*/
extern double geolonq[NYMEM][NXMEM]; /* q, u, or v points in degrees  */
extern double geolonu[NYMEM][NXMEM]; /* of longitude or m.            */
extern double geolonv[NYMEM][NXMEM];

extern double g[NZ];                 /* g[k] is the reduced gravity   */
                                     /* at interface k, in m s-2.     */
extern double Rlay[NZ+1];            /* The target potential density  */
                                     /* in each layer in kg m-3.      */

extern char directory[];             /* The directory to use to save  */
                                     /* the energy and restart files. */

extern int pe_here;                  /* The current processor's label.*/
extern int X1off, Y1off;             /* The offset in the global      */
                                     /* indices of X1 and Y1 on the   */
                                     /* current processor, relative   */
                                     /* to X1 and Y1 on PE 0.         */
extern int nx, ny;                   /* The largest index value in the*/
                                     /* x- and y- directions of the   */
                                     /* local computational domain.   */

static void USER_set_coordinate(int from_file, char filename[]);
static void USER_initialize_bottom_depth(int from_file, char filename[]);
static void USER_initialize_thickness(int from_file, char filename[]);
static void USER_initialize_velocity(int from_file, char filename[]);
#ifdef TEMPERATURE
static void USER_initialize_temperature_salinity(int from_file, char filename[]);
#endif
#ifdef BULKMIXEDLAYER
static void USER_initialize_mixed_layer_density(void);
#endif
#ifdef SPONGE
static void USER_initialize_sponges(int from_file, char filename[],
                                    char damp_file[]);
#endif
#ifdef NONLINEAR_EOS
static void USER_set_ref_profile(void);
#endif

#define FROM_FILE 1
#define NOT_FROM_FILE 0
#define FAIL 1

void initialize(time_type *day) {
/* Arguments: day - Time of the start of the run.  day is intent out. */

  char input_filename[200];             /* The name of the save file. */
  char filename[200];                   /* The name of an input file. */
  int i, j, new_sim;
  double phi0;
/*    See if the run is to be started from saved conditions, and get  */
/*  the names of the I/O directories and initialization file.  This   */
/*  subroutine also calls the subroutine that allows run-time changes */
/*  in parameters.  GetInputLines is in initialize_output.c           */
  GetInputLines(input_filename,sizeof(input_filename));

/*   At this point, add the prototypes and call to any user-supplied  */
/*  tracer initialization subroutines.                                */
  {
  /*  extern void USER_register_tracer_example(void); */
  /*  USER_register_tracer_example();                 */
  }

/*====================================================================*/
/*    Initialize fields that are time invariant - metrics, topography,*/
/*  masks, vertical coordinate, Coriolis parameter, a reference T & S */
/*  profile for compensating for compressibility.                     */
  {
/*    Set up the parameters of the physical domain, potentially       */
/*  by reading a grid file.                                           */
    set_metrics();
    Start_time = double_to_time(0.0,TIMEUNIT);

/*    This call calculates the layer densities and reduced gravities. */
/*  The two arguments here are (1) a flag indicating whether the      */
/*  layer densities are to be read from a file and (2) the path to    */
/*  that file.                                                        */
    strcpy(filename,HIM_params.inputdir);
    strcat(filename,"/init.nc");
    USER_set_coordinate(FROM_FILE,filename);

/*  The two arguments here are (1) a flag indicating whether the      */
/*  bottom depths are to be read from a file and (2) the path to      */
/*  that file.                                                        */
    strcpy(filename,HIM_params.inputdir);
    strcat(filename,"/grd.nc");
    USER_initialize_bottom_depth(FROM_FILE,filename);

/*    This call sets masks that prohibit flow over any point with     */
/*  a bottom that is shallower than the argument.                     */
    initialize_masks(MINIMUM_DEPTH);

/*    Calculate the value of the Coriolis parameter at the latitude   */
/*  of the q grid points, in s-1. */
    phi0 = (36. /180.)*M_PI;                             
    for (j=Y1-1;j<=ny+1;j++) for (i=X1-1;i<=nx+1;i++) {
      f[j][i] = 2.0*OMEGA*sin(phi0) + 2.0*OMEGA*(cos(phi0)/RE)*geolatq[j][i];
    }
    printf("geolatq: %f \n",geolatq[9][30]);
    printf("geolatq: %f \n",geolatq[8][30]);
    printf("f(1,30): %f \n",f[1][30]);
    printf("f(19,30): %f \n",f[19][30]);

#ifdef NONLINEAR_EOS
    USER_set_ref_profile();
#endif
  }

/*====================================================================*/
/*    Initialize temporally evolving fields, either as initial        */
/*  conditions or by reading them from a restart (or saves) file.     */
  if ((input_filename[0]=='n') && ((int) strlen(input_filename) == 1))
    new_sim = 1;
  else new_sim = 0;

  if (new_sim == 1)  {
/*  This block initializes all of the fields internally.              */
    if (pe_here == 0) printf("Run initialized internally.\n");

/*  The two arguments here are (1) a flag indicating whether the      */
/*  initial thicknesses are to be read from a file and (2) the        */
/*  path to that file.                                                */
    strcpy(filename,HIM_params.inputdir);
    strcat(filename,"/init.nc");
    USER_initialize_thickness(FROM_FILE,filename);

/*  The two arguments here are (1) a flag indicating whether the      */
/*  initial velocities are to be read from a file and (2) the         */
/*  path to that file.                                                */
    strcpy(filename,HIM_params.inputdir);
    strcat(filename,"/init.nc");
    USER_initialize_velocity(FROM_FILE,filename);

#ifdef TEMPERATURE
/*  The two arguments here are (1) a flag indicating whether the      */
/*  initial temperature and salinity fields are to be read from a     */
/*  file and (2) the path to that file.                               */
    strcpy(filename,HIM_params.inputdir);
    strcat(filename,"/MESO_Init_2.nc");
    USER_initialize_temperature_salinity(FROM_FILE,filename);
#endif

#ifdef BULKMIXEDLAYER
/*  The mixed layer density is typically determined from T & S.       */
    USER_initialize_mixed_layer_density();
#endif

    *day = Start_time;

/*  write_grid_file is in initialize_output.c.                        */
    write_grid_file();

  } else {
/*    This line calls a subroutine that reads the initial conditions  */
/*  from a previously generated file.                                 */
    restore_state(input_filename,directory,day);
  }

#ifdef SPONGE
/* The 3 arguments here are (1) a flag indicating whether the sponge  */
/* values are to be read from files, (2) the name of a file containing*/
/* the state toward which the model is damped, and (3) the file in    */
/* which the 2-D damping rate field can be found.                     */
  {
    char filename2[200];              /* The name of an input files.  */
    strcpy(filename,HIM_params.inputdir);
    strcat(filename,"/bdy.nc");
    strcpy(filename2,HIM_params.inputdir);
    strcat(filename2,"/sponge.nc");
    USER_initialize_sponges(FROM_FILE,filename,filename2);
  }
#endif

/* USER_set_diagnostics and USER_set_output are in initialize_output.c*/
  USER_set_diagnostics();

  USER_set_output(*day);

/* User-supplied tracer initialization routines are called here.      */
  call_tracer_init_fns(!new_sim, *day);

}


static void USER_set_coordinate(int from_file, char filename[]) {
/* This subroutine sets the layer densities (Rlay) and the interface  */
/* reduced gravities (g).                                             */

/* Arguments: from_file - 1 if the variables that are set here are to */
/*                        be read from a file; 0 to be set internally.*/
/*  (in)      filename - The name of the file to read.                */

  int k;

  if (from_file) {
    int err, cdfid;
    err = open_input_file(filename,&cdfid,NULL);
    if (err != 0) {
      printf("PE %d: Unable to open %s.\n",pe_here,filename); quit(-10);
    }
    Read_Field(cdfid,"LAYER",'1',NZ, 0, FAIL, &Rlay[0]);
    close_file(&cdfid);

    g[0] = GFS;
    for (k=1;k<NZ;k++) g[k] = (G_EARTH/RHO_0) * (Rlay[k] - Rlay[k-1]);
    //stef
    printf("g[0]: %f\n g[1]: %f\n Rlay[0]: %f\n Rlay[1]: %f\n", g[0],g[1],
       Rlay[0],Rlay[1]);
  } else {
    int set_g=1; /* If 0, a profile of T and S is specified, and g    */
                 /* are set to match.  Otherwise, g and the uppermost */
                 /* layer's density (or T & S) is specified, and the  */
                 /* rest of the layers' densities are set accordingly.*/
    if (set_g == 0) {
/*     Option 1 - set a target profile of T and S, and calculate the  */
/*       layer densities from these profiles, or (alternately) set the*/
/*       layer densities directly.  The reduced gravities across the  */
/*       interfaces are set to match the layer densities.             */
#ifdef TEMPERATURE
/*
      double T_tar[NZ] = {22.014, 20.024, 16.590, ..., -0.320};
      double S_tar[NZ] = {37.03, 36.89, 36.26, ..., 34.90};
      double pres[NZ];
      for (k=0;k<NZ;k++) pres[k] = P_REF;

      calculate_density(T_tar, S_tar, pres, Rlay, 0, NZ);
      for (k=1;k<NZ;k++) if (Rlay[k] <= Rlay[k-1])
        printf("WARNING: Inverted target densities R[%d] = %g, R[%d] = %g!\n",
                k-1,Rlay[k-1],k,Rlay[k]);
*/
#else
/*    The following lines could be replaced by something like ...     */
/*      double R_vals[NZ] = {1033.0, 1033.2, ... };         */
/*      for (k=0;k<NZ;k++) Rlay[k] = R_vals[k];             */
/*    and the next few lines would then be commented out.             */
      Rlay[0] = RHO_0;
      for (k=1;k<NZ;k++) Rlay[k] = Rlay[k-1] + GINT*RHO_0/G_EARTH;
#endif

/*    These statements set the interface reduced gravities.           */
      g[0] = GFS;
      for (k=1;k<NZ;k++) g[k] = (G_EARTH/RHO_0) * (Rlay[k] - Rlay[k-1]);

    } else {
/*    Option 2 - set the reduced gravities across the interfaces and  */
/*      the density of the top layer, and determine the deeper layer  */
/*      densities accordingly.                                        */

/*    These statements set the interface reduced gravities.           */
  /*  double g_val[NZ] = {GFS,   0.005, 0.005, 0.003, ..., 0.001 }; */
  /*  for (k=0;k<NZ;k++) g[k] = g_val[k];                           */
      g[0] = GFS;
      for (k=1;k<NZ;k++) g[k] = GINT;

/*    The uppermost layer's density is set here.  Subsequent layers'  */
/*  densities are determined from this value and the g values.        */
#ifdef TEMPERATURE
      {
        double T0 = 28.228, S0 = 34.5848, pref = P_REF;
        calculate_density(&T0, &S0, &pref, Rlay, 0, 1);
      }
#else
      Rlay[0] = RHO_0;
#endif

/*    These statements set the layer densities.                       */
      for (k=1;k<NZ;k++) Rlay[k] = Rlay[k-1] + g[k]*RHO_0/G_EARTH;
    }
  }
}


static void USER_initialize_bottom_depth(int from_file, char filename[]) {
/*  This subroutine places the bottom depth in m into D[][].          */

/* Arguments: from_file - 1 if the variables that are set here are to */
/*                        be read from a file; 0 to be set internally.*/
/*  (in)      filename - The name of the file to read.                */

  if (from_file==FROM_FILE) {
/*  Read the depth from a netcdf file.                                */
    int err, cdfid;
    err = open_input_file(filename,&cdfid,NULL);
    if (err != 0) {
      printf("PE %d: Unable to open %s.\n",pe_here,filename); quit(-10);
    }

    Read_Field(cdfid,"D",'h',1, 0, FAIL, &D[0][0]);

    close_file(&cdfid);
  } else {
/*  Calculate the depth of the bottom.                                */
    double D0;                  /* A constant to make the maximum     */
                                /* basin depth MAXIMUM_DEPTH.         */
    double expdecay = 400000.0; /* A decay scale of associated with   */
                                /* the sloping boundaries, in m.      */
    double Dedge = 100.0;       /* The depth in m at the basin edge.  */
    int i, j;

    D0 = (MAXIMUM_DEPTH - Dedge) /
           ((1.0 - exp(-0.5*LENLAT*RE*M_PI/(180.0 *expdecay))) *
            (1.0 - exp(-0.5*LENLAT*RE*M_PI/(180.0 *expdecay))));

    for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {
/*  This sets a bowl shaped (sort of) bottom topography, with a       */
/*  maximum depth of MAXIMUM_DEPTH.                                   */
      D[j][i] =  Dedge + D0 *
         (sin(M_PI * geolonh[j][i] / LENLON) *
         ((1.0 - exp(-(geolath[j][i] - (LOWLAT))*RE*M_PI/(180.0*expdecay))) *
         (1.0 - exp((geolath[j][i] - (LOWLAT+LENLAT))*RE*M_PI/(180.0*expdecay)))));

/*  This sets the bottom to be flat.                                  */
/*      D[j][i] = MAXIMUM_DEPTH;  */

      if (D[j][i] > MAXIMUM_DEPTH) D[j][i] = MAXIMUM_DEPTH;

      if (D[j][i] < MINIMUM_DEPTH) D[j][i] = 0.5*MINIMUM_DEPTH;
    }
  }
}

static void USER_initialize_thickness(int from_file, char filename[]) {
/*  This function puts the initial layer thicknesses into h[0][][][], */
/* using D[][] to constraint the thicknesses.                         */

/* Arguments: from_file - 1 if the variables that are set here are to */
/*                        be read from a file; 0 to be set internally.*/
/*  (in)      filename - The name of the file to read.                */

  int i, j, k;

  if (from_file==FROM_FILE) {
/*  Read the thicknesses from a netcdf file.                          */
    int err, cdfid;
    err = open_input_file(filename,&cdfid,NULL);
    if (err != 0) {
      printf("PE %d: Unable to open %s.\n",pe_here,filename); quit(-10);
    }

/*    If the file contains thicknesses, use the next line and comment */
/*  out the subsequent block.  If the file contains interface heights,*/
/*  use the extended block instead.                                   */

  /*  Read_Field(cdfid,"h",'h',NZ,0,FAIL,&h[0][0][0][0]);  */
    {
      double eta[NZ+1][NYMEM][NXMEM];
      int inconsistent = 0;
      Read_Field(cdfid,"ETA",'h',NZ+1,0,FAIL,&eta[0][0][0]);

      for (k=NZ-1;k>=0;k--) for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {
        if (eta[k][j][i] < (eta[k+1][j][i] + EPSILON)) {
          eta[k][j][i] = eta[k+1][j][i] + EPSILON;
          h[0][k][j][i] = EPSILON;
        } else h[0][k][j][i] = eta[k][j][i] - eta[k+1][j][i];
      }

/* Check for consistency between the interface heights and topography.*/
      for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++)
        if (fabs(eta[NZ][j][i] + D[j][i]) > 1.0) inconsistent++;
      sum_int_fields(&inconsistent,1);

      if ((inconsistent > 0) && (pe_here == 0))
        printf("WARNING: Thickness initial conditions inconsistent with "
            "topography in %d places!\n",inconsistent);
    }

    close_file(&cdfid);

  } else {
    double e0[NZ];    /* The resting interface heights, in m, usually */
                      /* negative because it is positive upward.      */
    double e_pert[NZ];/* Interface height perturbations, positive     */
                      /* upward, in m.                                */
    double eta[NZ+1]; /* Interface height relative to the sea surface */
                      /* positive upward, in m.                       */

    for (k=0;k<NZ;k++) e_pert[k] = 0.0;

    for (k=0;k<NZ;k++) e0[k] = -MAXIMUM_DEPTH * (double) k / (double) NZ;

    for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {

/*    This sets e_pert, the upward interface height perturbations     */
/*  relative to e0[k].                                                */
  /*    for (k=0;k<NZ;k++) e_pert[k] = e_here[k] - e0[k];  */

/*  The remainder of this subroutine should not be changed.           */

/*    This sets the initial thickness (in m) of the layers.  The      */
/*  thicknesses are set to insure that: 1.  each layer is at least    */
/*  EPSILON thick, and 2.  the interfaces are where they should be    */
/*  based on the resting depths and interface height perturbations,   */
/*  as long at this doesn't interfere with 1.                         */
      eta[NZ] = -D[j][i];
      for (k=NZ-1;k>=0;k--) {
        eta[k] = e0[k] + e_pert[k];
        if (eta[k] < (eta[k+1] + EPSILON)) {
          eta[k] = eta[k+1] + EPSILON;
          h[0][k][j][i] = EPSILON;
        } else h[0][k][j][i] = eta[k] - eta[k+1];
      }
    }
  }
}

static void USER_initialize_velocity(int from_file, char filename[]) {
/*   This function puts the initial layer velocities into u[0][][][]  */
/* and v[0][][][].                                                    */

/* Arguments: from_file - 1 if the variables that are set here are to */
/*                        be read from a file; 0 to be set internally.*/
/*  (in)      filename - The name of the file to read.                */

  if (from_file==FROM_FILE) {
/*  Read the velocities from a netcdf file.                           */
    int err, cdfid;
    err = open_input_file(filename,&cdfid,NULL);
    if (err != 0) {
      printf("PE %d: Unable to open %s.\n",pe_here,filename); quit(-10);
    }

    Read_Field(cdfid,"u", 'h', NZ, 0, FAIL, &u[0][0][0][0]);
    Read_Field(cdfid,"v", 'h', NZ, 0, FAIL, &v[0][0][0][0]);

    close_file(&cdfid);
  }
  else {
    int i, j, k;
    for (k=0;k<=NZ-1;k++) for (j=Y1-1;j<=ny;j++) for (i=X1-1;i<=nx;i++) {
      u[0][k][j][i] = 0.0; v[0][k][j][i] = 0.0;
    }
  }
}

#ifdef TEMPERATURE
static void USER_initialize_temperature_salinity(int from_file,
                  char filename[]) {
/*  This function puts the initial layer temperatures and salinities  */
/* into T[][][] and S[][][].                                          */

/* Arguments: from_file - 1 if the variables that are set here are to */
/*                        be read from a file; 0 to be set internally.*/
/*  (in)      filename - The name of the file to read.                */

  if (from_file==FROM_FILE) {
/* Read the temperatures and salinities from a netcdf file.           */
    int err, cdfid;

    err = open_input_file(filename,&cdfid,NULL);
    if (err != 0) {
      printf("PE %d: Unable to open %s.\n",pe_here,filename); quit(-10);
    }

    Read_Field(cdfid,"PTEMP", 'h', NZ, 0, FAIL, &T[0][0][0]);
    Read_Field(cdfid,"SALT", 'h', NZ, 0, FAIL, &S[0][0][0]);

    close_file(&cdfid);

  } else {
    double T0[NZ], S0[NZ];
    double pres[NZ];     /* Reference pressure in kg m-3.             */
    double drho_dT[NZ];  /* Derivative of density with temperature in */
                         /* kg m-3 K-1.                               */
    double drho_dS[NZ];  /* Derivative of density with salinity in    */
                         /* kg m-3 PSU-1.                             */
    double rho_guess[NZ];/* Potential density at T0 & S0 in kg m-3.   */
    int i,j,k,itt;

    for (k=0;k<NZ;k++) pres[k] = P_REF;
    for (k=0;k<NZ;k++) S0[k] = 35.0;
    T0[0] = 25.0;

    calculate_density(T0,S0,pres,rho_guess,0,1);
    calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,0,1);

/* A first guess of the layers' temperatures.                         */
    for (k=NZ-1;k>=0;k--)
      T0[k] = T0[0] + (Rlay[k] - rho_guess[0]) / drho_dT[0];

/* Refine the guesses for each layer.                                 */
    for (itt = 1; itt <= 6; itt++) {
      calculate_density(T0,S0,pres,rho_guess,0,NZ);
      calculate_density_derivs(T0,S0,pres,drho_dT,drho_dS,0,NZ);
      for (k=NZ-1;k>=0;k--)
        T0[k] = T0[k] + (Rlay[k] - rho_guess[k]) / drho_dT[k];
    }

    for (k=0;k<=NZ-1;k++) for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {
      T[k][j][i] = T0[k]; S[k][j][i] = S0[k];
    }

  }
}
#endif

#ifdef BULKMIXEDLAYER
static void USER_initialize_mixed_layer_density(void) {
/* Set the initial mixed layer and buffer layer densities.            */
/*
  if (from_file==FROM_FILE) {
  } else {
*/
    double pres[NXMEM];  /* Reference pressure, in Pa.                */
    int i, j, k;
    for (i=X1-1;i<=nx;i++) pres[i] = P_REF;

# ifdef TEMPERATURE
    for (k=0;k<=NML;k++) {
      for (j=Y1-1;j<=ny;j++) {
        calculate_density(T[k][j],S[k][j],pres,Rml[k][j],X1-1,nx-X1+2);
      }
    }
# else
    for (k=0;k<=NML;k++) for (j=Y1-1;j<=ny;j++) for (i=X1-1;i<=nx;i++)
      Rml[k][j][i] = RHO_0+(g[1]+g[2])*RHO_0/G_EARTH;
# endif
/*}*/
}
#endif

#ifdef SPONGE
static void USER_initialize_sponges(int from_file, char filename[],
                                    char damp_file[]) {
/*   This subroutine sets the inverse restoration time (Idamp), and   */
/* the values toward which the interface heights and an arbitrary     */
/* number of tracers should be restored within each sponge. The       */
/* interface height is always subject to damping, and must always be  */
/* the first registered field.                                        */

/* Arguments: from_file - 1 if the variables that are used here are to*/
/*                        be read from a file; 0 to be set internally.*/
/*  (in)      filename - The name of the file to read for all fields  */
/*                       except the inverse damping rate.             */
/*  (in)      damp_file - The name of the file from which to read the */
/*                        inverse damping rate.                       */

  double temp[NZ][NYMEM][NXMEM]; /* A temporary array.                */
  double Idamp[NYMEM][NXMEM];    /* The inverse damping rate, in s-1. */

  int i, j, k;

  if (from_file==FROM_FILE) {
    int err, cdfid;

    err = open_input_file(damp_file,&cdfid,NULL);
    if (err != 0) {
      printf("PE %d: Unable to open %s.\n",pe_here,damp_file); quit(-10);
    }

    Read_Field(cdfid,"Idamp", 'h', 1, 0, FAIL, &Idamp[0][0]);
    /*stef*/
    printf("Idamp: %f \n", Idamp[1][1]);
    initialize_sponge(Idamp);

    close_file(&cdfid);

/* Now register all of the fields which are damped in the sponge.     */
/* By default, momentum is advected vertically within the sponge, but */
/* momentum is typically not damped within the sponge.                */

    err = open_input_file(filename,&cdfid,NULL);
    if (err != 0) {
      printf("PE %d: Unable to open %s.\n",pe_here,filename); quit(-10);
    }

/*  The first call to set_up_sponge_field is for the interface height.*/
    Read_Field(cdfid,"ETA", 'h', NZ+1, 0, FAIL, &temp[0][0][0]);
    for (k=NZ-1;k>=0;k--) for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {
      if (temp[k][j][i] < (temp[k+1][j][i] + EPSILON))
        temp[k][j][i] = temp[k+1][j][i] + EPSILON;
    }
    set_up_sponge_field(temp,NULL,NZ);

# ifdef BULKMIXEDLAYER
/*   The second call to set_up_sponge_field must be for Rml if        */
/* BULKMIXEDLAYER is defined. The remaining calls can be in any order.*/
/* Only a single layer's worth of the Rml reference is used, even if  */
/* there are multiple parts of the mixed layer (i.e. NML>1).          */
# ifdef TEMPERATURE
    {
      double pres[NXMEM];  /* Reference pressure, in Pa.              */
      int i, j;
      for (i=X1-1;i<=nx;i++) pres[i] = P_REF;

      Read_Field(cdfid,"PTEMP", 'h', 1, 0, FAIL, &temp[1][0][0]);
      Read_Field(cdfid,"SALT", 'h', 1, 0, FAIL, &temp[2][0][0]);

      for (j=Y1;j<=ny;j++)
        calculate_density(temp[1][j],temp[2][j],pres,temp[0][j],X1,nx-X1+1);
    }
#  else
    Read_Field(cdfid,"Rml", 'h', 1, 0, FAIL, &temp[0][0][0]);
#  endif

    set_up_sponge_field(temp,&Rml,1);
# endif

/*  The remaining calls to set_up_sponge_field can be in any order.   */
# ifdef TEMPERATURE
    Read_Field(cdfid,"PTEMP", 'h', NZ, 0, FAIL, &temp[0][0][0]);
    set_up_sponge_field(temp,&T,NZ);
    Read_Field(cdfid,"SALT", 'h', NZ, 0, FAIL, &temp[0][0][0]);
    set_up_sponge_field(temp,&S,NZ);
# endif
//stef
    Read_Field(cdfid,"u", 'u', NZ, 0, FAIL, &temp[0][0][0]);
    set_up_sponge_field(temp,&u,NZ);
    Read_Field(cdfid,"v", 'v', NZ, 0, FAIL, &temp[0][0][0]);
    set_up_sponge_field(temp,&v,NZ);
    close_file(&cdfid);
  } else {
/*  Here the inverse damping time, in s-1, is set. Set Idamp to 0     */
/*  wherever there is no sponge, and the subroutines that are called  */
/*  will automatically set up the sponges only where Idamp is positive*/
/*  and hmask is 1.                                                   */

/*  This particular example sets sponges at the eastern and western   */
/*  sides of the domain.                                              */
    double Idamp_all[NXTOT+1];
    int sp_width_w = 9;
    int sp_width_e = 9;
    for (i=1;i<=NXTOT;i++) Idamp_all[i] = 0.0;

    Idamp_all[sp_width_w] = 1.0/(5*86400.0);
    Idamp_all[sp_width_w-1] = 1.0/(4*86400.0);
    Idamp_all[sp_width_w-2] = 1.0/(3*86400.0);
    for (i=sp_width_w-3;i>=1;i--) Idamp_all[i] = 2.0*Idamp_all[i+1];

    Idamp_all[NXTOT+1-sp_width_e] = 1.0/(5.0*86400.0);
    Idamp_all[NXTOT+2-sp_width_e] = 1.0/(4.0*86400.0);
    Idamp_all[NXTOT+3-sp_width_e] = 1.0/(3.0*86400.0);
    for (i=NXTOT+4-sp_width_e;i<=NXTOT;i++) Idamp_all[i] = 2.0*Idamp_all[i-1];

    for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {
      Idamp[j][i] = Idamp_all[i+X1off-X1+1];
    }
    initialize_sponge(Idamp);

/* Now register all of the fields which are damped in the sponge.     */
/* By default, momentum is advected vertically within the sponge, but */
/* momentum is typically not damped within the sponge.                */

/*  The first call to set_up_sponge_field is for the interface height.*/
    {
/*
      double e_Med[NZ] = {   0.0,    0.0,  -10.0,  -36.0,  -71.0,
                          -109.0, -139.0, -158.0, -168.0, -174.0,
                          -176.0, -178.0, -183.0, -188.0, -200.0, -240.0};
      double e_Atl[NZ] = {   0.0,   -80.0,  -170.0,  -380.0,  -620.0,
                          -900.0, -1300.0, -1800.0, -2400.0, -3000.0,
                         -3800.0, -6000.0, -6000.0, -6000.0, -6000.0, -6000.0};
      for (i=X1;i<=nx;i++) for (j=Y1;j<=ny;j++) for (k=0;k<NZ;k++) {
        temp[k][j][i] = (geolonh[j][i] > LENLON + (WESTLON) -2.5) ? e_Med[k] : e_Atl[k];
      }
*/
      set_up_sponge_field(temp,NULL,NZ);
    }

# ifdef BULKMIXEDLAYER
/*   The second call to set_up_sponge_field must be for Rml, if       */
/* BULKMIXEDLAYER is defined. The remaining calls can be in any order.*/
/* Only a single layer's worth of the Rml reference is used, even if  */
/* there are multiple parts of the mixed layer (i.e. NML>1).          */
    {
      double Rml_Med, Rml_Atl = 0.0;
      Rml_Med = Rlay[2];
      for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {
        temp[0][j][i] = (geolonh[j][i] > LENLON + (WESTLON) -2.5) ? Rml_Med : Rml_Atl;
      }
      set_up_sponge_field(temp,&Rml,1);
    }
# endif
/*  The remaining calls to set_up_sponge_field can be in any order.   */
# ifdef TEMPERATURE
    {
      /*
      double T_Atl[NZ] = {22.014, 20.024, 16.590, 12.916, 9.786,
                           7.563,  6.716,  5.361,  4.844, 3.981,
                           3.450,  3.159,  2.685,  2.182, 1.361, -0.320};
      double T_Med[NZ] = {22.035, 20.399, 18.700, 17.192, 15.819,
                          14.912, 14.388, 14.183, 14.112, 14.049,
                          14.034, 14.013, 13.984, 13.886, 13.787, 13.624};
      for (k=0;k<NZ;k++) for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {
        temp[k][j][i] = (geolonh[j][i] > LENLON + (WESTLON) -2.5) ? T_Med[k] : T_Atl[k];
      }
      */
      set_up_sponge_field(temp,&T,NZ);
    }
    {
      /*
      double S_Med[NZ] = {37.04, 37.04, 37.04, 37.14, 37.24,
                          37.39, 37.54, 37.67, 37.75, 37.77,
                          37.79, 37.81, 37.84, 37.90, 38.01, 38.22};
      double S_Atl[NZ] = {37.03, 36.89, 36.26, 35.68, 35.34,
                          35.21, 35.20, 35.21, 35.20, 35.06,
                          34.98, 34.95, 34.90, 34.90, 34.90, 34.90};
      for (k=0;k<NZ;k++) for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {
        temp[k][j][i] = (geolonh[j][i] > LENLON + (WESTLON) -2.5) ? S_Med[k] : S_Atl[k];
      }
      */
      set_up_sponge_field(temp,&S,NZ);
    }
# endif
  }
}
#endif                                                      /* SPONGE */

#ifdef NONLINEAR_EOS
static void USER_set_ref_profile(void) {
/* This subroutine specifies a reference profile to use to compensate */
/* for the compressibility of seawater.  This profile should be both  */
/* representative and smooth.  A fairly uniform (with depth) distrib- */
/* ution of points should be selected.                                */

/* This particular example is for a reference profile in the eastern  */
/* North Atlantic, somewhere near 20 W, 30 N.                         */
/*
  double T_ref[27] = {22.01, 19.26, 13.39, 10.05, 8.36,
                       7.35,  6.93,  6.42,  5.82, 5.36,
                       5.19,  5.02,  4.84,  4.50, 4.22,
                       3.98,  3.80,  3.65,  3.54, 3.45,
                       3.39,  3.34,  3.29,  3.24, 3.19,
                       3.17,  3.16};
  double S_ref[27] = {37.03, 36.75, 35.76, 35.37, 35.26,
                      35.21, 35.21, 35.20, 35.20, 35.20,
                      35.20, 35.20, 35.18, 35.14, 35.10,
                      35.06, 35.04, 35.02, 35.00, 34.98,
                      34.97, 34.96, 34.96, 34.95, 34.95,
                      34.95, 34.95};
  double e_ref[27] = {   0.0,  -100.0,   -350.0,  -600.0,  -800.0,
                      -1000.0, -1200.0, -1400.0, -1600.0, -1800.0,
                      -2000.0, -2200.0, -2400.0, -2600.0, -2800.0,
                      -3000.0, -3200.0, -3400.0, -3600.0, -3800.0,
                      -4000.0, -4200.0, -4400.0, -4600.0, -4800.0,
                      -5000.0, -5200.0};
  register_compress(T_ref, S_ref, e_ref, 27);
*/
/* This particular example is vaguely representative of the ACC south */
/* of 55s.                                                            */
  double T_ref[13] = {2.00, 1.93, 1.73, 1.63, 0.9,
                      0.1, -0.5, -0.5, -0.5, -0.5,
                      -0.5, -0.5, -0.5};
  double S_ref[13] = {34.65, 34.65, 34.65, 34.65, 34.65, 34.65, 34.65,
                      34.65, 34.65, 34.65, 34.65, 34.65, 34.65};
  double e_ref[13] = {    0.0,  -100.0,  -300.0,  -500.0, -1000.0,
                      -1500.0, -2000.0, -2500.0, -3000.0, -3500.0,
                      -4000.0, -4500.0, -5000.0};
  register_compress(T_ref, S_ref, e_ref, 13);
  
}
#endif
