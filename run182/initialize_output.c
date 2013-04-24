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
 *    Subroutines in this file specify the diagnostic fields that are *
 *  calculated (USER_set_diagnostics), and which fields are written   *
 *  as snapshots or as time means (USER_set_output).  An additional   *
 *  subroutine, GetInputLines, prompts the user for runtime input     *
 *  which includes setting directory, the path to the directory for   *
 *  all of the subsequent saves.  write_grid_file writes out one file *
 *  containing grid information, and another file containing the      *
 *  topography and several other fields that do not change through    *
 *  the course of the run.  The order of the subroutines in this      *
 *  file reflects the likelihood that users will want to modify them, *
 *  not the order in which they will probably be called.              *
 *                                                                    *
 *  Variables written all in capital letters are defined in init.h    *
 *                                                                    *
 *     A small fragment of the grid is shown below:                   *
 *                                                                    *
 *    j+1  x ^ x ^ x   At x:  q, f                                    *
 *    j+1  > o > o >   At ^:  v                                       *
 *    j    x ^ x ^ x   At >:  u                                       *
 *    j    > o > o >   At o:  h, D, T, S                              *
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

extern struct diag_fld diag;       /* A structure containing pointers */
                                   /* to diagnostic fields which might*/
                                   /* be calculated.                  */

#ifdef BULKMIXEDLAYER
extern double Rml[NML+1][NYMEM][NXMEM]; /* The mixed and buffer layer */
                                     /* potential densities in kg m-2.*/
#endif
#ifdef TEMPERATURE
extern double T[NZ][NYMEM][NXMEM];   /* Potential temperature, in C.  */
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

char directory[80];                  /* The directory to use to save  */
                                     /* the energy and restart files. */
char directory2[80];                 /* The directory to use to save  */
                                     /* output files.                 */

extern int pe_here;                  /* The current processor's label.*/
extern int nx, ny;                   /* The largest index value in the*/
                                     /* x- and y- directions of the   */
                                     /* local computational domain.   */

static double *allocate_dbl(int size);

void USER_set_output(time_type day) {
/*   This subroutine sets up the fields to be output, either as       */
/* snapshots or time averages.  A save file is set up by first calling*/
/* specify_saves_file to register the name-root, and path of the file,*/
/* along with arguments specifying whether it is intended to contain  */
/* time-averaged or instantaneous quantities, and time_type arguments */
/* indicating the current time, the start time, stop time, and save   */
/* interval for a file, and the number of time levels that can go into*/
/* each file.  A number of calls to set_up_save_field are then made   */
/* register each of the fields to be saved to that file.              */

/* Argument: day - Time of the start of the run.                      */

/*   The following structures describes the output and metadata for   */
/* each field.                                                        */
/*   vardesc is a structure defined in HIM_io.h.   The elements of    */
/* this structure, in order, are: (1) the variable name for the NetCDF*/
/* file; (2) the variable's long name; (3) a character indicating the */
/* horizontal grid, which may be '1' (column), 'h', 'q', 'u', or 'v', */
/* for the corresponding C-grid variable; (4) a character indicating  */
/* the vertical grid, which may be 'L' (layer), 'i' (interface),      */
/* '2' (mixed-layers), or '1' (no vertical location); (5) a character */
/* indicating the time levels of the field, which may be 's' (snap-   */
/* shot), 'a' (average between snapshots), 'm' (monthly average), or  */
/* '1' (no time variation); (6) the variable's units; and (7) a       */
/* character indicating the size in memory to write, which may be     */
/* 'd' (8-byte) or 'f' (4-byte).                                      */
  vardesc var_desc_u = {"u","Zonal velocity",'u','L','s',"meter second-1",'f'};
  vardesc var_desc_v = {"v","Meridional velocity",'v','L','s',"meter second-1",'f'};
  vardesc var_desc_h = {"h","Layer Thickness",'h','L','s',"meter",'d'};
#ifdef BULKMIXEDLAYER
  vardesc var_desc_Rml = {"Rml","Mixed Layer Potential Density",
                          'h','2','s',"kg meter-3", 'f'};
#endif
#ifdef TEMPERATURE
  vardesc var_desc_T = {"Temp","Potential Temperature",'h','L','s',"degC", 'f'};
  vardesc var_desc_S = {"Salt","Salinity",'h','L','s',"PSU", 'f'};
#endif
  vardesc var_desc_e = {"e","Interface Height",'h','i','s',"meter",'d'};
  vardesc var_desc_e_D = {"emD","Interface Height above Bottom",'h','i','s',
                          "meter",'d'};

  specify_saves_file("save", directory2, 's', day, HIM_params.daysave,
        HIM_params.daymax, HIM_params.saveint, HIM_params.saves_per_file);

/* To change the variables in the output files, comment, uncomment,   */
/* or add lines calling set_up_save_field.                            */

  set_up_save_field(&u[0][0],&u[1][0],var_desc_u);
  set_up_save_field(&v[0][0],&v[1][0],var_desc_v);
  set_up_save_field(&h[0][0],&h[1][0],var_desc_h);
  if (diag.e != NULL)
    set_up_save_field(&diag.e[0],NULL,var_desc_e);
  if (diag.e_D != NULL)
    set_up_save_field(&diag.e_D[0],NULL,var_desc_e_D);
#ifdef BULKMIXEDLAYER
  set_up_save_field(&Rml[0],NULL,var_desc_Rml);
#endif
#ifdef TEMPERATURE
  set_up_save_field(&T[0],NULL,var_desc_T);
  set_up_save_field(&S[0],NULL,var_desc_S);
#endif


/*  Set up the fields to be averaged.                                 */
  if (HIM_params.save_aves > 0) {
    extern double uh[NZ][NYMEM][NXMEM];/* uh = u * h * dy at u points.*/
    extern double vh[NZ][NYMEM][NXMEM];/* vh = v * h * dx at v points.*/

/*  These structures describe the output and metadata for each field. */
/*
    vardesc var_desc_uh = {"uhtm","Time Mean Zonal Thickness Flux",
                           'u','L','a',"meter3 second-1",'f'};
    vardesc var_desc_vh = {"vhtm","Time Mean Meridional Thickness Flux",
                           'v','L','a',"meter3 second-1",'f'};
    vardesc var_desc_h = {"htm","Time Mean Layer Thickness",
                          'h','L','a',"meter",'d'};
    vardesc var_desc_u = {"utm","Time Mean Zonal Velocity",
                          'u','L','a',"meter second-1",'f'};
    vardesc var_desc_v = {"vtm","Time Mean Meridional Velocity",
                          'v','L','a',"meter second-1",'f'};
*/
    vardesc var_desc_dia = {"wd","Time Mean Diapycnal Velocity",
                            'h','i','a',"meter second-1",'f'};
    vardesc var_desc_e = {"etm","Time Mean Interface Height",
                           'h','i','a',"meter",'d'};
    vardesc var_desc_e_D = {"emDtm","Time Mean Interface Height above Bottom",
                            'h','i','a',"meter",'d'};
    vardesc var_desc_rv = {"rvtm","Time Mean Relative Vorticity",
                           'q','L','a',"second-1",'f'};
    vardesc var_desc_q = {"q","Time Mean layer Potential Vorticity",
                          'q','L','a',"meter-1 second-1",'f'};


/*   specify_saves_file provides certain information about the file   */
/* and indicates the file that will be used for subsequently          */
/* registered fields.                                                 */
    specify_saves_file("avfld", directory2, 'a', day, HIM_params.daysave,
         HIM_params.daymax, HIM_params.saveint, HIM_params.saves_per_file);

/*   set_up_save_field registers a field to be saved or averaged.     */
/* To change the variables in the average files, comment, uncomment,  */
/* or add lines calling set_up_save_field.                            */

/*
    set_up_save_field(&uh[0],NULL,var_desc_uh);
    set_up_save_field(&vh[0],NULL,var_desc_vh);
    set_up_save_field(&h[0][0],&h[1][0],var_desc_h);
    set_up_save_field(&u[0][0],&u[1][0],var_desc_u);
    set_up_save_field(&v[0][0],&v[1][0],var_desc_v);
*/
    if (diag.e != NULL)
      set_up_save_field(&diag.e[0],NULL,var_desc_e);
    if (diag.e_D != NULL)
      set_up_save_field(&diag.e_D[0],NULL,var_desc_e_D);

    if (diag.diapyc_vel != NULL)
      set_up_save_field(&diag.diapyc_vel[0],NULL,var_desc_dia);
    if (diag.rv != NULL)
      set_up_save_field(&diag.rv[0],NULL,var_desc_rv);
    if (diag.q != NULL)
      set_up_save_field(&diag.q[0],NULL,var_desc_q);

/* Following are a variety of terms in the momentum equations.        */
    {
      vardesc var_desc_PFu = {"PFu","Zonal Pressure Force Acceleration",
        'u','L','a',"meter second-2",'f'};
      vardesc var_desc_PFv = {"PFv","Meridional Pressure Force Acceleration",
        'v','L','a',"meter second-2",'f'};
      vardesc var_desc_du_dt = {"dudt","Zonal Acceleration",'u','L','a',"meter second-2",'f'};
      vardesc var_desc_dv_dt = {"dvdt","Meridional Acceleration",'v','L','a',"meter second-2",'f'};
      vardesc var_desc_CAu = {"CAu","Zonal Coriolis and Advective Acceleration",
        'u','L','a',"meter second-2",'f'};
      vardesc var_desc_CAv = {"CAv","Meridional Coriolis and Advective Acceleration",
        'v','L','a',"meter second-2",'f'};
      vardesc var_desc_diffu = {"diffu","Zonal Acceleration from Horizontal Viscosity",
        'u','L','a',"meter second-2",'f'};
      vardesc var_desc_diffv = {"diffv","Meridional Acceleration"
        "from Horizontal Viscosity",'v','L','a',"meter second-2",'f'};
      vardesc var_desc_vertvu = {"vertu","Zonal Acceleration"
        "from Vertical Viscosity",'u','L','a',"meter second-2",'f'};
      vardesc var_desc_vertvv = {"vertv","Meridional Acceleration"
        "from Vertical Viscosity",'v','L','a',"meter second-2",'f'};
      vardesc var_desc_dia_u = {"diau","Zonal Acceleration"
        "from Diapycnal Mixing",'u','L','a',"meter second-2",'f'};
      vardesc var_desc_dia_v = {"diav","Meridional Acceleration"
        "from Diapycnal Mixing",'v','L','a',"meter second-2",'f'};
      vardesc var_desc_gKEu = {"gKEu","Zonal Acceleration"
        "from Grad. Kinetic Energy",'u','L','a',"meter second-2",'f'};
      vardesc var_desc_gKEv = {"gKEv","Meridional Acceleration"
        "from Grad. Kinetic Energy",'v','L','a',"meter second-2",'f'};
      vardesc var_desc_rvxv = {"rvxv","Zonal Acceleration"
        "from Relative Vorticity",'u','L','a',"meter second-2",'f'};
      vardesc var_desc_rvxu = {"rvxu","Meridional Acceleration"
        "from Relative Vorticity",'v','L','a',"meter second-2",'f'};
      vardesc var_desc_PFubc = {"PFubc","Density Gradient Zonal"
        "Pressure Force Accel.",'u','2','a',"meter second-2",'f'};
      vardesc var_desc_PFvbc = {"PFvbc","Density Gradient Meridional"
        "Pressure Force Accel.",'v','2','a',"meter second-2",'f'};

/*    To select any of the following variables,  allocate the memory  */
/*  for it in USER_set_diagnostics.                                   */
      if (diag.PFu != NULL)
        set_up_save_field(&diag.PFu[0],NULL,var_desc_PFu);
      if (diag.PFv != NULL)
        set_up_save_field(&diag.PFv[0],NULL,var_desc_PFv);
      if (diag.CAu != NULL)
        set_up_save_field(&diag.CAu[0],NULL,var_desc_CAu);
      if (diag.CAv != NULL)
        set_up_save_field(&diag.CAv[0],NULL,var_desc_CAv);
      if (diag.diffu != NULL)
        set_up_save_field(&diag.diffu[0],NULL,var_desc_diffu);
      if (diag.diffv != NULL)
        set_up_save_field(&diag.diffv[0],NULL,var_desc_diffv);
      if (diag.du_dt != NULL)
        set_up_save_field(&diag.du_dt[0],NULL,var_desc_du_dt);
      if (diag.dv_dt != NULL)
        set_up_save_field(&diag.dv_dt[0],NULL,var_desc_dv_dt);
      if (diag.du_dt_visc != NULL)
        set_up_save_field(&diag.du_dt_visc[0],NULL,var_desc_vertvu);
      if (diag.dv_dt_visc != NULL)
        set_up_save_field(&diag.dv_dt_visc[0],NULL,var_desc_vertvv);
      if (diag.du_dt_dia != NULL)
        set_up_save_field(&diag.du_dt_dia[0],NULL,var_desc_dia_u);
      if (diag.dv_dt_dia != NULL)
        set_up_save_field(&diag.dv_dt_dia[0],NULL,var_desc_dia_v);
      if (diag.gradKEu != NULL)
        set_up_save_field(&diag.gradKEu[0],NULL,var_desc_gKEu);
      if (diag.gradKEv != NULL)
        set_up_save_field(&diag.gradKEv[0],NULL,var_desc_gKEv);
      if (diag.rv_x_v != NULL)
        set_up_save_field(&diag.rv_x_v[0],NULL,var_desc_rvxv);
      if (diag.rv_x_u != NULL)
        set_up_save_field(&diag.rv_x_u[0],NULL,var_desc_rvxu);
      if (diag.PFu_bc != NULL)
        set_up_save_field(&diag.PFu_bc[0],NULL,var_desc_PFubc);
      if (diag.PFv_bc != NULL)
        set_up_save_field(&diag.PFv_bc[0],NULL,var_desc_PFvbc);
    }

#ifdef SPLIT
    {
      vardesc var_desc_PFubt = {"PFuBT","Zonal Barotropic Pressure Force"
        "Acceleration",'u','1','a',"meter second-2",'f'};
      vardesc var_desc_PFvbt = {"PFvBT","Meridional Barotropic Pressure Force"
        "Acceleration",'v','1','a',"meter second-2",'f'};
      vardesc var_desc_Corubt = {"CoruBT","Zonal Barotropic Coriolis"
        "Acceleration",'u','1','a',"meter second-2",'f'};
      vardesc var_desc_Corvbt = {"CorvBT","Meridional Barotropic Coriolis"
        "Acceleration",'v','1','a',"meter second-2",'f'};
      vardesc var_desc_Nonlnu = {"NluBT","Zonal Barotropic"
        "Nonlinear Acceleration",'u','1','a',"meter second-2",'f'};
      vardesc var_desc_Nonlnv = {"NlvBT","Meridional Barotropic"
        "Nonlinear Acceleration",'v','1','a',"meter second-2",'f'};
      vardesc var_desc_uhbt = {"uhBT","Zonal Barotropic"
        "Mass Flux",'u','1','a',"meter3 second-1",'f'};
      vardesc var_desc_vhbt = {"vhBT","Zonal Barotropic"
        "Mass Flux",'v','1','a',"meter3 second-1",'f'};

      if (diag.PFu_bt != NULL)
        set_up_save_field((double (*)[NYMEM][NXMEM]) diag.PFu_bt,NULL,var_desc_PFubt);
      if (diag.PFv_bt != NULL)
        set_up_save_field((double (*)[NYMEM][NXMEM]) diag.PFv_bt,NULL,var_desc_PFvbt);
      if (diag.Coru_bt != NULL)
        set_up_save_field((double (*)[NYMEM][NXMEM]) diag.Coru_bt,NULL,var_desc_Corubt);
      if (diag.Corv_bt != NULL)
        set_up_save_field((double (*)[NYMEM][NXMEM]) diag.Corv_bt,NULL,var_desc_Corvbt);
      if (diag.Nonlnu_bt != NULL)
        set_up_save_field((double (*)[NYMEM][NXMEM]) diag.Nonlnu_bt,NULL,var_desc_Nonlnu);
      if (diag.Nonlnv_bt != NULL)
        set_up_save_field((double (*)[NYMEM][NXMEM]) diag.Nonlnv_bt,NULL,var_desc_Nonlnv);
      if (diag.ubt_flux != NULL)
        set_up_save_field((double (*)[NYMEM][NXMEM]) diag.ubt_flux,NULL,var_desc_uhbt);
      if (diag.vbt_flux != NULL)
        set_up_save_field((double (*)[NYMEM][NXMEM]) diag.vbt_flux,NULL,var_desc_vhbt);
    }

#endif

/* Following are a variety of incremental quantities.                 */
    {
      vardesc var_desc_au = {"au_visc","Zonal Viscous Vertical Coupling"
          "Coefficient",'u','L','a',"meter second-1",'f'};
      vardesc var_desc_av = {"av_visc","Meridional Viscous Vertical Coupling"
          "Coefficient",'v','L','a',"meter second-1",'f'};
      vardesc var_desc_Ahh = {"Ahh","Biharmonic Horizontal Viscosity at h Points",
          'h','L','a',"meter4 second-1",'f'};
      vardesc var_desc_Ahq = {"Ahq","Biharmonic Horizontal Viscosity at q Points",
          'q','L','a',"meter4 second-1",'f'};
      vardesc var_desc_Khh = {"Khh","Laplacian Horizontal Viscosity at h Points",
          'h','L','a',"meter2 second-1",'f'};
      vardesc var_desc_Khq = {"Khq","Laplacian Horizontal Viscosity at q Points",
          'q','L','a',"meter2 second-1",'f'};
/*    To select any of the following variables,  allocate the memory  */
/*  for it in USER_set_diagnostics.                                   */
      if (diag.au_visc != NULL)
        set_up_save_field(&diag.au_visc[0],NULL,var_desc_au);
      if (diag.av_visc != NULL)
        set_up_save_field(&diag.av_visc[0],NULL,var_desc_av);
      if (diag.Ah_h != NULL)
        set_up_save_field(&diag.Ah_h[0],NULL,var_desc_Ahh);
      if (diag.Ah_q != NULL)
        set_up_save_field(&diag.Ah_q[0],NULL,var_desc_Ahq);
      if (diag.Kh_h != NULL)
        set_up_save_field(&diag.Kh_h[0],NULL,var_desc_Khh);
      if (diag.Kh_q != NULL)
        set_up_save_field(&diag.Kh_q[0],NULL,var_desc_Khq);
    }

  }

}

void USER_set_diagnostics(void) {
/*   This subroutine is used to specify which diagnostic quantities   */
/* will be calculated.  Those fields in diag that contain pointers to */
/* an appropriate location in memory will be calculated, while those  */
/* fields that are NULL will not be calculated.                       */

  extern double CAu[NZ][NYMEM][NXMEM];/* CAu = fv - u.grad(u), m s-2. */
  extern double CAv[NZ][NYMEM][NXMEM];/* CAv = -fu - u.grad(v), m s-2.*/
  extern double PFu[NZ][NYMEM][NXMEM];/* PFu = -dM/dx, in m s-2.      */
  extern double PFv[NZ][NYMEM][NXMEM];/* PFv = -dM/dy, in m s-2.      */
  extern double diffu[NZ][NYMEM][NXMEM];/* Horizontal diffusive       */
  extern double diffv[NZ][NYMEM][NXMEM];/* accelerations, in m s-2.   */

/*   Set every element of diag to NULL.  These will be replaced with  */
/* pointers to the fields to be calculated by uncommenting lines later*/
/* in this subroutine.                                                */
  diag.e = NULL;  diag.e_D = NULL; diag.diapyc_vel = NULL;

  diag.du_dt_visc = NULL; diag.dv_dt_visc = NULL;
  diag.du_dt_dia = NULL; diag.dv_dt_dia = NULL;
  diag.diffu = NULL; diag.diffv = NULL;
  diag.CAu = NULL; diag.CAv = NULL;
  diag.PFu = NULL; diag.PFv = NULL;
  diag.gradKEu = NULL; diag.gradKEv = NULL;
  diag.rv_x_v = NULL; diag.rv_x_u = NULL;
  diag.rv = NULL; diag.q = NULL;
  diag.du_dt = NULL; diag.dv_dt = NULL;

  diag.au_visc = NULL; diag.av_visc = NULL;
  diag.Ah_h = NULL; diag.Ah_q = NULL;
  diag.Kh_h = NULL; diag.Kh_q = NULL;

  diag.PFu_bc = NULL; diag.PFv_bc = NULL;

  diag.PFu_tot = NULL; diag.PFv_tot = NULL;
  diag.CAu_tot = NULL; diag.CAv_tot = NULL;
  diag.PFu_bt = NULL; diag.PFv_bt = NULL;
  diag.Coru_bt = NULL; diag.Corv_bt = NULL;
  diag.Nonlnu_bt = NULL; diag.Nonlnv_bt = NULL;
  diag.ubt_flux = NULL; diag.vbt_flux = NULL;

  diag.u_trunc_diag_file = NULL; diag.v_trunc_diag_file = NULL;

/*   To use a diagnostic field, just allocate the memory for the      */
/* space to holds that field.  To do that, just uncomment the appro-  */
/* priate line below.  See the documentation of the diag structure    */
/* in HIM.h for a description of each of these possible fields.       */

  diag.e = (double (*)[NYMEM][NXMEM]) allocate_dbl((NZ+1)*NYMEM*NXMEM);
//  diag.e_D = (double (*)[NYMEM][NXMEM]) allocate_dbl((NZ+1)*NYMEM*NXMEM);
  diag.diapyc_vel = (double (*)[NYMEM][NXMEM]) allocate_dbl((NZ+1)*NYMEM*NXMEM);

//  diag.rv = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
//  diag.q = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
//  diag.du_dt_visc = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
//  diag.dv_dt_visc = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
  diag.du_dt_dia = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
  diag.dv_dt_dia = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
//  diag.diffu = &diffu[0];
//  diag.diffv = &diffv[0];
//  diag.gradKEu = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
//  diag.gradKEv = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
//  diag.rv_x_v = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
//  diag.rv_x_u = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
//  diag.du_dt = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
//  register_time_deriv(&u[0][0],&u[1][0],diag.du_dt,NZ);
//  diag.dv_dt = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
//  register_time_deriv(&v[0][0],&v[1][0],diag.dv_dt,NZ);

//  diag.au_visc = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
//  diag.av_visc = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
#ifdef BIHARMONIC
//  diag.Ah_h = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
//  diag.Ah_q = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
#endif
#ifdef LAPLACIAN
//  diag.Kh_h = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
//  diag.Kh_q = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
#endif

#ifdef BULKMIXEDLAYER
//  diag.PFu_bc = (double (*)[NYMEM][NXMEM]) allocate_dbl((NML+1)*NYMEM*NXMEM);
//  diag.PFv_bc = (double (*)[NYMEM][NXMEM]) allocate_dbl((NML+1)*NYMEM*NXMEM);
#endif

#ifdef SPLIT
//  diag.PFu_bt = (double (*)[NXMEM]) allocate_dbl(NYMEM*NXMEM);
//  diag.PFv_bt = (double (*)[NXMEM]) allocate_dbl(NYMEM*NXMEM);
//  diag.Coru_bt = (double (*)[NXMEM]) allocate_dbl(NYMEM*NXMEM);
//  diag.Corv_bt = (double (*)[NXMEM]) allocate_dbl(NYMEM*NXMEM);
//  diag.Nonlnu_bt = (double (*)[NXMEM]) allocate_dbl(NYMEM*NXMEM);
//  diag.Nonlnv_bt = (double (*)[NXMEM]) allocate_dbl(NYMEM*NXMEM);
  diag.ubt_flux = (double (*)[NXMEM]) allocate_dbl(NYMEM*NXMEM);
  diag.vbt_flux = (double (*)[NXMEM]) allocate_dbl(NYMEM*NXMEM);
//  diag.CAu_tot = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
//  diag.CAv_tot = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
//  diag.PFu_tot = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);
//  diag.PFv_tot = (double (*)[NYMEM][NXMEM]) allocate_dbl(NZ*NYMEM*NXMEM);

/* The lines in the following block should not be commented out.  To  */
/* write CAu, CAv, etc., uncomment CAu_tot, CAv_tot, etc., above.     */
  {
    if ((diag.PFu_bt == NULL) && ((diag.PFu_tot != NULL) || (diag.CAu_tot != NULL)))
      diag.PFu_bt = (double (*)[NXMEM]) allocate_dbl(NYMEM*NXMEM);
    if ((diag.PFv_bt == NULL) && ((diag.PFv_tot != NULL) || (diag.CAv_tot != NULL)))
      diag.PFv_bt = (double (*)[NXMEM]) allocate_dbl(NYMEM*NXMEM);
    diag.CAu = diag.CAu_tot;
    diag.CAv = diag.CAv_tot;
    diag.PFu = diag.PFu_tot;
    diag.PFv = diag.PFv_tot;
  }
#else
//  diag.CAu = &CAu[0];
//  diag.CAv = &CAv[0];
//  diag.PFu = &PFu[0];
//  diag.PFv = &PFv[0];
#endif

/*   The following 6 lines specify files to use to report the accel-  */
/* erations in columns with velocities larger than HIM_params.maxvel. */
/* This reporting can be disabled by commenting out these lines.      */

/*
  diag.u_trunc_diag_file = (char *) calloc(120,sizeof(char));
  strcpy(diag.u_trunc_diag_file,directory2);
  strcat(diag.u_trunc_diag_file,"U_diag_report");
  diag.v_trunc_diag_file = (char *) calloc(120,sizeof(char));
  strcpy(diag.v_trunc_diag_file,directory2);
  strcat(diag.v_trunc_diag_file,"V_diag_report");
*/

}

void GetInputLines(char filename[], int namelength) {
/*  This subroutine reads a few input fields to determine how the     */
/* model is to be initialized.                                        */

/* Arguments: filename - A string that determines how the model state */
/*                       is to be initialized; filename is intent out.*/
/*  (in)      namelength - The available length of filename.          */
  char nl_filename[200];
  int k;

  if (pe_here == 0) {
    printf("  Enter the directory (including final /) for initial conditions.\n");
    fgets(directory,sizeof(directory),stdin);
    directory[(int) strlen(directory)-1] = '\0';
    while (directory[(int) strlen(directory)-1] == ' ')
      directory[(int) strlen(directory)-1] = '\0';
    if ((directory[(int) strlen(directory)-1] != '/') &&
        ((int) strlen(directory) > 0)) {
      int endpoint;
      endpoint = strlen(directory);
      directory[endpoint] = '/'; directory[endpoint+1] = '\0';
      /* printf("Will use directory %s.\n",directory); */
    }

    printf("\n  Enter 'n' to start a new run, 'r' to restart from the restart file(s),\n");
    printf("or the names of the input files separated by spaces (without directory).\n");
    fgets(filename,namelength,stdin);
    filename[(int) strlen(filename)-1] = '\0';
    while (filename[(int) strlen(filename)-1] == ' ')
      filename[(int) strlen(filename)-1] = '\0';

    printf("\n  Enter the directory to use to write output files,\n");
    printf("or nothing to use the same directory as used for input.\n");
    fgets(directory2,sizeof(directory2),stdin);
    if ((int) strlen(directory2) == 0) directory2[0] = '\0';
    else directory2[(int) strlen(directory2)-1] = '\0';
    while (directory2[(int) strlen(directory2)-1] == ' ')
      directory2[(int) strlen(directory2)-1] = '\0';
    if ((int) strlen(directory2) == 0)
      for (k=0;k<sizeof(directory);k++) directory2[k] = directory[k];
    else {
      if ((directory2[(int) strlen(directory2)-1] != '/') &&
          ((int) strlen(directory2) > 0)) {
        int endpoint;
        endpoint = strlen(directory2);
        directory2[endpoint] = '/'; directory2[endpoint+1] = '\0';
        /* printf("Will use directory %s.\n",directory2); */
      }
    }

    printf("\n  Enter the name of the namelist file to use to modify some of the\n"
           "parameters that are set in init.h at run time.  Alternately, enter\n"
           "nothing just to use the values set in init.h\n");
    printf("  The format for each line of the namelist file is the same as in init.h,\n"
           "for example \"#  define DT 100.0 /* COMMENTS MAY FOLLOW */\".\n"
           "Any values which are changed will be reported.  Not all variables\n"
           "in init.h can be changed in this way.\n");

    fgets(nl_filename,sizeof(nl_filename),stdin);
    if ((int) strlen(nl_filename) == 0) nl_filename[0] = '\0';
    else nl_filename[(int) strlen(nl_filename)-1] = '\0';
    while (nl_filename[(int) strlen(nl_filename)-1] == ' ')
      nl_filename[(int) strlen(nl_filename)-1] = '\0';
  }

/*  Broadcast filename, directory, and directory2 from pe 0.          */
  spread_string(filename,namelength);
  spread_string(directory,sizeof(directory));
  spread_string(directory2,sizeof(directory2));
  spread_string(nl_filename,sizeof(nl_filename));

#ifdef CHECK_RESTART
  if (pe_here > 0) {
    strcpy(filename,"r");
    strcpy(directory,"saves/");
    strcpy(directory2,"saves/");
    strcpy(nl_filename,"input");
  }
#endif

/*   Call the subroutine that will change parameters based on the     */
/* namelist file.                                                     */
  HIM_parser(nl_filename);
}

void write_grid_file(void) {
/*    The following subroutine saves the depth, reduced gravity of    */
/*  each interface, the potential density of each layer, and the      */
/*  Coriolis parameter.  A variety of metric terms are written to a   */
/*  separate file.                                                    */
  char filename[40];           /* The file name of the save file.     */
  char filepath[120];          /* The path (dir/file) to the file.    */
  int i, j, cdfid, timeid;
  size_t err = 1;

  struct varcdfinfo varinfo[10];
/*   vardesc is a structure defined in HIM_io.h.   The elements of    */
/* this structure, in order, are: (1) the variable name for the NetCDF*/
/* file; (2) the variable's long name; (3) a character indicating the */
/* horizontal grid, which may be '1' (column), 'h', 'q', 'u', or 'v', */
/* for the corresponding C-grid variable; (4) a character indicating  */
/* the vertical grid, which may be 'L' (layer), 'i' (interface),      */
/* '2' (mixed-layers), or '1' (no vertical location); (5) a character */
/* indicating the time levels of the field, which may be 's' (snap-   */
/* shot), 'a' (average between snapshots), 'm' (monthly average), or  */
/* '1' (no time variation); (6) the variable's units; and (7) a       */
/* character indicating the size in memory to write, which may be     */
/* 'd' (8-byte) or 'f' (4-byte).                                      */
  vardesc vars[4] = {
    {"D","Basin Depth",'h','1','1',"meter", 'd'},
    {"g","Reduced gravity",'1','L','1',"meter second-2", 'd'},
    {"R","Target Potential Density",'1','L','1',"kilogram meter-3", 'd'},
    {"f","Coriolis Parameter",'q','1','1',"second-1", 'd'}
  };

  sprintf(filename,"D.%d.%d.%d",NXTOT,NYTOT,NZ);
  strcpy(filepath, directory);
  strcat(filepath, filename);

  create_file(filepath, vars, 4, &cdfid, &timeid, varinfo);
  err *= write_field(cdfid, vars[0], varinfo[0], 0, D[0]);
  err *= write_field(cdfid, vars[1], varinfo[1], 0, g);
  err *= write_field(cdfid, vars[2], varinfo[2], 0, Rlay);
  err *= write_field(cdfid, vars[3], varinfo[3], 0, f[0]);

  close_file(&cdfid);
  if (err == 0)
    printf("Problems saving general parameters.\n");

  {
/* ### REVISIT THIS PART OF THE ROUTINE AFTER METRICS... ### */
    double out[NYMEM][NXMEM];     /* An array for output.  */
    extern double hmask[NYMEM][NXMEM];

    vardesc vars2[10]={
      {"geolatb","latitude at q points",'q','1','1',"degree",'d'},
      {"geolonb","longitude at q points",'q','1','1',"degree",'d'},
      {"wet", "land or ocean?", 'h','1','1',"none",'d'},
      {"geolat", "latitude at h points", 'h','1','1',"degree",'d'},
      {"geolon","longitude at h points",'h','1','1',"degree",'d'},
      {"dxh","Zonal grid spacing at h points",'h','1','1',"m",'d'},
      {"dxq","Zonal grid spacing at q points",'q','1','1',"m",'d'},
      {"dyh","Meridional grid spacing at h points",'h','1','1',"m",'d'},
      {"dyq","Meridional grid spacing at q points",'q','1','1',"m",'d'},
      {"Ah","Area of h cells",'h','1','1',"m2",'d'},
    };

    sprintf(filename,"grid.%d.%d",NXTOT,NYTOT);
    strcpy(filepath, directory);
    strcat(filepath, filename);

    create_file(filepath, vars2, 10, &cdfid, &timeid, varinfo);
    err *= write_field(cdfid, vars2[0], varinfo[0], 0, geolatq[0]);
    err *= write_field(cdfid, vars2[1], varinfo[1], 0, geolonq[0]);
    err *= write_field(cdfid, vars2[2], varinfo[2], 0, hmask[0]);
    err *= write_field(cdfid, vars2[3], varinfo[3], 0, geolath[0]);
    err *= write_field(cdfid, vars2[4], varinfo[4], 0, geolonh[0]);

    for (j=0;j<NYMEM;j++) for (i=0;i<NXMEM;i++) out[j][i]=DXh(j,i);
    err *= write_field(cdfid, vars2[5], varinfo[5], 0, out[0]);
    for (j=0;j<NYMEM;j++) for (i=0;i<NXMEM;i++) out[j][i]=DXq(j,i);
    err *= write_field(cdfid, vars2[6], varinfo[6], 0, out[0]);
    for (j=0;j<NYMEM;j++) for (i=0;i<NXMEM;i++) out[j][i]=DYh(j,i);
    err *= write_field(cdfid, vars2[7], varinfo[7], 0, out[0]);
    for (j=0;j<NYMEM;j++) for (i=0;i<NXMEM;i++) out[j][i]=DYq(j,i);
    err *= write_field(cdfid, vars2[8], varinfo[8], 0, out[0]);
    for (j=0;j<NYMEM;j++) for (i=0;i<NXMEM;i++) out[j][i]=DXDYh(j,i);
    err *= write_field(cdfid, vars2[9], varinfo[9], 0, out[0]);

    close_file(&cdfid);

    if (err == 0)
      printf("Problems saving latitude and longitude and wet.\n");
  }
}

static double *allocate_dbl(int size) {
/* This subroutine allocates space for a double array of size size,   */
/* initializes it to 0, and returns a pointer to it.                  */
  double *ptr;
  int i;
  ptr = (double *) calloc(size,sizeof(double));
  if (ptr == NULL)
    printf("Unable to allocate %d words for a diagnostic field.\n",size);
  else for (i=0;i<size;i++) ptr[i] = 0.0;
  return ptr;
}
