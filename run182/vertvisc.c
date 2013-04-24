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
 *  By Robert Hallberg, April 1994 - December 2000                    *
 *  Quadratic Bottom Drag by James Stephens and R. Hallberg.          *
 *                                                                    *
 *    This file contains the subroutine that implements vertical      *
 *  viscosity (vertvisc) and the subroutine that calculates the       *
 *  viscosity and thickness of the BBL (set_viscous_BBL).             *
 *                                                                    *
 *    The vertical diffusion of momentum is fully implicit.  This is  *
 *  necessary to allow for vanishingly small layers.  The coupling    *
 *  is based on the distance between the centers of adjacent layers,  *
 *  except where a layer is close to the bottom compared with a       *
 *  bottom boundary layer thickness when a bottom drag law is used.   *
 *  A stress top b.c. and a no slip bottom  b.c. are used.  There     *
 *  is no limit on the time step for vertvisc.                        *
 *                                                                    *
 *    Near the bottom, the horizontal thickness interpolation scheme  *
 *  changes to an upwind biased estimate to control the effect of     *
 *  spurious Montgomery potential gradients at the bottom where       *
 *  nearly massless layers layers ride over the topography.  Within a *
 *  few boundary layer depths of the bottom, the harmonic mean        *
 *  thickness (i.e. (2 h+ h-) / (h+ + h-) ) is used if the velocity   *
 *  is from the thinner side and the arithmetic mean thickness        *
 *  (i.e. (h+ + h-)/2) is used if the velocity is from the thicker    *
 *  side.  Both of these thickness estimates are second order         *
 *  accurate.  Above this the arithmetic mean thickness is used.      *
 *                                                                    *
 *    In addition, vertvisc truncates any velocity component that     *
 *  exceeds maxvel to truncvel. This basically keeps instabilities    *
 *  spatially localized.  The number of times the velocity is         *
 *  truncated is reported each time the energies are saved, and if    *
 *  exceeds HIM_params.Maxtrunc the model will stop itself and change *
 *  the day to 9.9e9.  This has proven very useful in (1) diagnosing  *
 *  model failures and (2) letting the model settle down to a         *
 *  meaningful integration from a poorly specified initial condition. *
 *                                                                    *
 *    If the split time stepping is used, the time average velocities *
 *  of each layer (uav & vav) are determined by separating the        *
 *  vertical average of the instantaneous velocities, and adding the  *
 *  time average barotropic velocities.  Both the instantaneous       *
 *  velocities and the time average velocities are subjected to both  *
 *  the vertical viscosity and the possible truncation.  The same     *
 *  thicknesses are used for the vertical diffusion of both the       *
 *  instantaneous and time average velocities.                        *
 *                                                                    *
 *    The same code is used for the two velocity components, by       *
 *  indirectly referencing the velocities and defining a handful of   *
 *  direction-specific defined variables.                             *
 *                                                                    *
 *  Variables written all in capital letters are defined in init.h    *
 *                                                                    *
 *     A small fragment of the grid is shown below:                   *
 *                                                                    *
 *    j+1  x ^ x ^ x   At x:  q                                       *
 *    j+1  > o > o >   At ^:  v, frhatv, tauy                         *
 *    j    x ^ x ^ x   At >:  u, frhatu, taux                         *
 *    j    > o > o >   At o:  h                                       *
 *    j-1  x ^ x ^ x                                                  *
 *        i-1  i  i+1  At x & ^:                                      *
 *           i  i+1    At > & o:                                      *
 *                                                                    *
 *  The boundaries always run through q grid points (x).              *
 *                                                                    *
 ********+*********+*********+*********+*********+*********+*********+*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <init.h>
#include <HIM.h>
#include <HIM_io.h>

extern double umask[NYMEM][NXMEM];   /* _mask are 1 over ocean and 0  */
extern double vmask[NYMEM][NXMEM];   /* over land on the u & v grids. */

extern struct params HIM_params;   /* HIM_params is a structure that  */
                                   /* contains a number of parameters */
                                   /* of the simulation which can be  */
                                   /* modified from a namelist file   */
                                   /* at run time.                    */

extern struct diag_fld diag;       /* A structure containing pointers */
                                   /* to diagnostic fields which might*/
                                   /* be calculated.                  */

extern long ntrunc;                  /* ntrunc is the number of times */
                                     /* the velocity has been trunc-  */
                                     /* ated since the last call to   */
                                     /* write_energy.                 */

#ifdef BOTTOMDRAGLAW
double bbl_thick[2][NYMEM][NXMEM];   /* The bottom boundary layer     */
                                     /* thickness at the 2 velocity   */
                                     /* points in m.                  */
double kv_bbl[2][NYMEM][NXMEM];      /* The bottom boundary layer     */
                                     /* viscosity at the 2 velocity   */
                                     /* points, in m2 s-1.            */
extern double Rlay[NZ+1];            /* The target potential density  */
                                     /* in each layer in kg m-3.      */
extern double f[NYMEM][NXMEM];       /* The Coriolis parameter in s-1.*/
#endif

extern int nx, ny;                   /* The largest index value in the*/
                                     /* x- and y- directions of the   */
                                     /* local computational domain.   */

void vertvisc(double u[NZ][NYMEM][NXMEM], double v[NZ][NYMEM][NXMEM],
        double h[NZ][NYMEM][NXMEM], struct forcing fluxes,
        struct bt_vars_ptrs bt, double dt,
        double uav[NZ][NYMEM][NXMEM], double vav[NZ][NYMEM][NXMEM]) {
/*    This subroutine does a fully implicit vertical diffusion        */
/*  of momentum.  Stress top and bottom b.c.s are used.               */

/* Arguments: u - Zonal velocity, in m s-1.  Intent in/out.           */
/*  (in/out)  v - Meridional velocity, in m s-1.                      */
/*  (in)      h - Layer thickness, in m.                              */
/*  (in)      fluxes - A structure containing pointers to any possible*/
/*                     forcing fields.  Unused fields have NULL ptrs. */
/*  (in)      bt - pointers to a collection of variables related to   */
/*                 the split time stepping scheme, all set elsewhere. */
/*  (in)      dt - Time increment in s.                               */
/*  (out)     uav - u with its vertical mean replaced by ubtav, m s-1.*/
/*  (out)     vav - v with its vertical mean replaced by vbtav, m s-1.*/

/* Fields from fluxes used in this subroutine:                        */
/*   taux: Zonal wind stress in Pa.                                   */
/*   tauy: Meridional wind stress in Pa.                              */
/*   ustar: the friction velocity in m s-1, used here as the mixing   */
/*     velocity in the mixed layer if NML > 1 in a bulk mixed layer.  */

/* Fields from bt used in this subroutine:                            */
/*  ubtav, vbtav: The barotropic zonal and meridional velocities      */
/*    averaged over a baroclinic time step, in m s-1.                 */
/*  IDatu, IDatv: Inverse of the basin depth at u or v points, in m-1.*/
/*  frhatu, frhatv:  Fraction of the total column thickness interpol- */
/*    ated to u and v grid points in each layer, nondimensional.      */

  double hvel[NZ][NXMEM];  /* hvel is the thickness at a velocity     */
                           /* grid point that is actually used, in m. */

  double a[NZ+1][NXMEM];   /* a is the drag coefficient across an     */
                           /* interface time integrated over dt, in m.*/
                           /* a times the velocity difference gives   */
                           /* the stress across an interface.         */
  double b1[NXMEM];        /* b1 and c1 are variables used by the     */
  double c1[NZ][NXMEM];    /* tridiagonal solver.                     */

#if !defined(BULKMIXEDLAYER) || defined(DIRECT_STRESS)
  double I_Hmix;           /* The inverse of the mixed layer          */
                           /* thickness, in m-1.                      */
#endif
  double I_Hbbl[NXMEM];    /* The inverse of the bottom boundary      */
                           /* layer thickness, in m-1.                */
  double Idt;              /* The inverse of the time step, in s-1.   */
  double dt_Rho0;          /* The time step divided by the mean       */
                           /* density, in s m3 kg-1.                  */

  double C2dtKv;          /* C2dtKv = 2*dt*Kv, in m2.                 */
#ifndef BULKMIXEDLAYER
  double C2dtKvml;         /* C2dtKvml = 2*dt*Kvml, in m2.            */
#endif
#ifndef BOTTOMDRAGLAW
  double C2dtKvbbl;        /* C2dtKvbbl = 2*dt*Kvbbl, in m2.          */
  double CdtKvbbl;         /* CdtKvbbl = dt*Kvbbl, in m2.             */
  double C2dtKvbblmKv;     /* C2dtKvbblmKv = 2*dt*(Kvbbl-Kv), in m2.  */
#endif

  double maxvel;           /* Velocities components greater than      */
  double truncvel;         /* maxvel are truncated to truncvel, both  */
                           /* in m s-1.                               */

  double stress[NXMEM];    /*   The masked surface stress times the   */
                           /* time step, divided by the density, in   */
                           /* units of m2 s-1.                        */
  double surface_stress[NXMEM];/* The same as stress, unless the wind */
                           /* stress is applied as a body force, in   */
                           /* units of m2 s-1.                        */
  int do_i[NXMEM];
  int i, j, k, m;
  int dowrite[2][NXMEM], maywrite[2];
  int do_ave;
  static int first_calc = 0;

  double (*U[2])[NZ][NYMEM][NXMEM]; /* Pointers to the velocities.    */
#ifdef SPLIT
  double (*Uav[2])[NZ][NYMEM][NXMEM]; /* Pointers to the velocities   */
                           /* with the altered vertical means.        */
  double velbtcor[NXMEM];  /* The difference between the time mean    */
                           /* barotropic velocity and the depth       */
                           /* averaged velocity in m s-1.             */
#endif

  do_ave = query_averaging_enabled();

  maywrite[0] = (diag.u_trunc_diag_file == NULL) ? 0 : 1;
  maywrite[1] = (diag.v_trunc_diag_file == NULL) ? 0 : 1;

  if (first_calc == 0) { first_calc = 1; }

  U[0] = (double (*)[NZ][NYMEM][NXMEM]) u;
  U[1] = (double (*)[NZ][NYMEM][NXMEM]) v;
#ifdef SPLIT
  Uav[0] = (double (*)[NZ][NYMEM][NXMEM]) uav;
  Uav[1] = (double (*)[NZ][NYMEM][NXMEM]) vav;
#endif

  dt_Rho0 = dt/RHO_0;
  Idt = 1.0 / dt;

#if !defined(BULKMIXEDLAYER) || defined(DIRECT_STRESS)
  I_Hmix = 1.0/HIM_params.Hmix;
#endif

  maxvel = HIM_params.maxvel;  truncvel = 0.9*maxvel;

  C2dtKv = 2.0*HIM_params.Kv*dt;
#ifndef BULKMIXEDLAYER
  C2dtKvml = 2.0*HIM_params.Kvml*dt;
#endif

#ifndef BOTTOMDRAGLAW
  for (i=X1;i<=nx;i++) I_Hbbl[i] = 1.0/HIM_params.Hbbl;
  CdtKvbbl = HIM_params.Kvbbl*dt; C2dtKvbbl = 2.0*CdtKvbbl;
  C2dtKvbblmKv = C2dtKvbbl-C2dtKv;
#endif                                           /* End BOTTOMDRAGLAW */

  for (j=Y1;j<=ny;j++) {

    if (diag.du_dt_visc != NULL) for (k=0;k<NZ;k++) for (i=X1;i<=nx;i++)
      diag.du_dt_visc[k][j][i] = u[k][j][i];
    if (diag.dv_dt_visc != NULL) for (k=0;k<NZ;k++) for (i=X1;i<=nx;i++)
      diag.dv_dt_visc[k][j][i] = v[k][j][i];

    for (m=0; m<=1; m++) {
      if (m==0) for (i=X1;i<=nx;i++) do_i[i] = (umask[j][i]>0.0) ? 1 : 0;
      else for (i=X1;i<=nx;i++) do_i[i] = (vmask[j][i]>0.0) ? 1 : 0;

      for (i=X1;i<=nx;i++) if (do_i[i]) {
        stress[i] = dt_Rho0 * ((m==0) ? (umask[j][i]*(*(fluxes.taux))[j][i]) :
                                        (vmask[j][i]*(*(fluxes.tauy))[j][i]));
        surface_stress[i] = stress[i];
#ifdef DIRECT_STRESS
/*   One option is to have the wind stress applied as a body force    */
/* over the topmost HIM_params.Hmix fluid.  If DIRECT_STRESS is not   */
/* defined, the wind stress is applied as a stress boundary condition.*/
        {
          double z, hfr;
          surface_stress[i] = 0.0;
          z = 0.0;
          for (k=0;k<NZ;k++) {
            if (z >= HIM_params.Hmix) break;
            else {
              double h_a;
              h_a = 0.5 * (h[k][j][i] + ((m==0)?h[k][j][i+1]:h[k][j+1][i]));
              if ((z+h_a) <= HIM_params.Hmix) hfr = I_Hmix;
              else hfr = (1.0 - z*I_Hmix) / h_a;
              (*U[m])[k][j][i] += hfr * stress[i];
              z += h_a;
            }
          }
        }
#endif                                               /* DIRECT_STRESS */
      }

#ifdef SPLIT
      if (m==0) {
        for (i=X1;i<=nx;i++) if (do_i[i])
          velbtcor[i] = u[0][j][i]*bt.frhatu[0][j][i] - bt.ubtav[j][i]
                      + stress[i]*bt.IDatu[j][i];
        for (k=1;k<NZ;k++) for (i=X1;i<=nx;i++) if (do_i[i])
          velbtcor[i] += u[k][j][i] * bt.frhatu[k][j][i];
      } else {
        for (i=X1;i<=nx;i++) if (do_i[i])
          velbtcor[i] = v[0][j][i]*bt.frhatv[0][j][i] - bt.vbtav[j][i]
                      + stress[i]*bt.IDatv[j][i];;
        for (k=1;k<NZ;k++) for (i=X1;i<=nx;i++) if (do_i[i])
          velbtcor[i] += v[k][j][i] * bt.frhatv[k][j][i];
      }
#endif

/**********************************************************************/
/*    The following loop calculates (1) the thicknesses at velocity   */
/*  grid points for the vertical viscosity (hvel[k]) and (2) the      */
/*  'coupling coefficient' (a[k]) at the interfaces. If BOTTOMDRAGLAW */
/*  is defined, the minimum of Hbbl and the half the adjacent layer   */
/*  thicknesses are used to calculate a[k] near the bottom.  Near the */
/*  bottom an upwind biased thickness is used to control the effect   */
/*  of spurious Montgomery potential gradients at the bottom where    */
/*  nearly massless layers layers ride over the topography.  The terms*/
/*  including 1e-10 in the denominators of the final expressions for  */
/*  a[k] are here to avoid truncation error problems in the           */
/*  tridiagonal solver.  Effectively, this sets the maximum coupling  */
/*  coefficient at 1e10 m.  'botfn' determines when a point is within */
/*  the influence of the bottom boundary layer, going from 1 at the   */
/*  bottom to 0 in the interior.                                      */
      {
/*    The following loop calculates the vertical average velocity and */
/*  surface mixed layer contributions to the vertical viscosity.      */
        double h_harm[NZ][NXMEM]; /* Harmonic mean of the thicknesses */
                             /* around a velocity grid point, defined */
                             /* by 2*(h+ * h-)/(h+ + h-), in m.       */
        double h_pl[NZ][NXMEM]; /* Thickness at i+1 or j+1 for u or v.*/
        double z[NXMEM];     /* The distance from the top or bottom,  */
                             /* normalized by Hmix or Hbbl, nondim.   */
        double r;            /* A thickness to compare with Hbbl in m.*/

        for (i=X1;i<=nx;i++) if (do_i[i]) {
          h_pl[0][i] = (m==0) ? h[0][j][i+1] : h[0][j+1][i];
          h_harm[0][i] = 2.0*h[0][j][i]*h_pl[0][i] / (h[0][j][i]+h_pl[0][i]);
#ifndef BULKMIXEDLAYER
          z[i] = 0.0;
#endif
        }

        for (k=1;k<NZ;k++) {
          for (i=X1;i<=nx;i++) if (do_i[i]) {
            h_pl[k][i] = (m==0) ? h[k][j][i+1] : h[k][j+1][i];
            h_harm[k][i] = 2.0*h[k][j][i]*h_pl[k][i] / (h[k][j][i]+h_pl[k][i]);

#ifdef BULKMIXEDLAYER
            a[k][i] = -C2dtKv;
#else
            z[i] += h_harm[k-1][i]*I_Hmix;
/*        a[k][i] = -C2dtKv - C2dtKvml / (1.0 + 0.09*z[i]*z[i]*z[i]*z[i]*z[i]*z[i]);  */
            a[k][i] = -C2dtKv - C2dtKvml / ((z[i]*z[i]) *
                           (1.0 + 0.09*z[i]*z[i]*z[i]*z[i]*z[i]*z[i]));
/*           (           (1.0 + 0.05*z[i]*z[i]*z[i]*z[i]*z[i]*z[i]*z[i]*z[i]));          */
#endif
          }
        }

        for (i=X1;i<=nx;i++) if (do_i[i]) {
          hvel[NZ-1][i] = ((*U[m])[NZ-1][j][i] * (h_pl[NZ-1][i]-h[NZ-1][j][i]) >= 0) ?
                        h_harm[NZ-1][i] : 0.5*(h_pl[NZ-1][i]+h[NZ-1][j][i]);
          z[i] = 0.0;
#ifdef BOTTOMDRAGLAW
          r = hvel[NZ-1][i]*0.5;
          a[NZ][i] = -dt*kv_bbl[m][j][i] / (1e-10*dt*kv_bbl[m][j][i] +
                (((r < bbl_thick[m][j][i]) ? r : bbl_thick[m][j][i])));
          I_Hbbl[i] = 1.0 / bbl_thick[m][j][i];
#else
          a[NZ][i] = -C2dtKvbbl / (hvel[NZ-1][i] + 1.0e-10*C2dtKvbbl);
#endif
        }

        for (k=NZ-1;k>=1;k--) {
          for (i=X1;i<=nx;i++) if (do_i[i]) {
            double botfn;    /* A function which goes from 1 at the   */
                             /* bottom to 0 much more than Hbbl into  */
                             /* the interior.                         */
            double h_shear;/* The distance over which shears occur, m.*/
            z[i] += h_harm[k][i]*I_Hbbl[i];
/*        botfn = 1.0 / (1.0 + 0.05*z[i]*z[i]*z[i]*z[i]*z[i]*z[i]*z[i]*z[i]); */
//stef
            botfn = (z[i] <= 1.0) ? 1.0 : ((z[i] >= 2.0) ? 0.0 : 2.0-z[i]); 
//            botfn = 1.0 / (1.0 + 0.09*z[i]*z[i]*z[i]*z[i]*z[i]*z[i]);

            hvel[k-1][i] = ((*U[m])[k-1][j][i] * (h_pl[k-1][i]-h[k-1][j][i]) < 0) ?
                (1.0-botfn)*h_harm[k-1][i] + botfn*0.5*(h_pl[k-1][i]+h[k-1][j][i]) :
                h_harm[k-1][i];

#ifdef BOTTOMDRAGLAW
            a[k][i] += -2.0*dt*(kv_bbl[m][j][i]-HIM_params.Kv)*botfn;
            r = (hvel[k][i]+hvel[k-1][i]);
            h_shear = (r > 2.0*bbl_thick[m][j][i]) ?
                      ((1.0 - botfn) * r + botfn*2.0*bbl_thick[m][j][i]) : r;
#else
            a[k][i] += -C2dtKvbblmKv*botfn;
            h_shear = (hvel[k][i]+hvel[k-1][i]);
#endif
            a[k][i] /= (h_shear - 1.0e-10*a[k][i]);
          }

        }

#ifdef BULKMIXEDLAYER
        for (i=X1;i<=nx;i++) if (do_i[i]) {
          double dtustar;
          dtustar = 0.5*dt * ( (*(fluxes.ustar))[j][i] + ((m==0) ?
              (*(fluxes.ustar))[j][i+1] : (*(fluxes.ustar))[j+1][i]) );
          for (k=1;k<NML;k++) a[k][i] = ((-a[k][i])<dtustar) ? -dtustar : a[k][i];
        }
#endif

      }

/* Store the coupling coefficients if the space is allocated.         */
      if (do_ave) {
        if ((m==0) && (diag.au_visc != NULL)) {
          for (i=X1;i<=nx;i++) diag.au_visc[0][j][i] = 0.0;
          for (k=1;k<NZ;k++) for (i=X1;i<=nx;i++)
            diag.au_visc[k][j][i] = (do_i[i] > 0) ? a[k][i]*Idt : 0.0;
        } else if (diag.av_visc != NULL) {
          for (i=X1;i<=nx;i++) diag.av_visc[0][j][i] = 0.0;
          for (k=1;k<NZ;k++) for (i=X1;i<=nx;i++)
            diag.av_visc[k][j][i] = (do_i[i] > 0) ? a[k][i]*Idt : 0.0;
        }
      }

/*  This is a modification of a standard tridiagonal solver.          */
      for (i=X1;i<=nx;i++) if (do_i[i]) {
        b1[i]=1.0/(hvel[0][i]-a[1][i]);
#ifdef SPLIT
        (*Uav[m])[0][j][i] = b1[i] * (hvel[0][i] * ((*U[m])[0][j][i]-velbtcor[i]) +
                                    surface_stress[i]);
#endif
        (*U[m])[0][j][i] = b1[i] * (hvel[0][i] * (*U[m])[0][j][i] + surface_stress[i]);
        dowrite[m][i] = 0;
      }
      for (k=1;k<=NZ-1;k++) {
        for (i=X1;i<=nx;i++) if (do_i[i]) {
          c1[k][i]=a[k][i]*b1[i];
          b1[i]=1.0/(hvel[k][i] - a[k+1][i] - a[k][i]*(1.0+c1[k][i]));
#ifdef SPLIT
          (*Uav[m])[k][j][i] = (hvel[k][i] * ((*U[m])[k][j][i] - velbtcor[i]) -
                         a[k][i] * (*Uav[m])[k-1][j][i]) * b1[i];
#endif
          (*U[m])[k][j][i] = (hvel[k][i] * (*U[m])[k][j][i] -
                          a[k][i] * (*U[m])[k-1][j][i]) * b1[i];
        }
      }
      for (k=(NZ-2);k>=0;k--) {
        for (i=X1;i<=nx;i++) if (do_i[i]) {
          (*U[m])[k][j][i] -= c1[k+1][i] * (*U[m])[k+1][j][i];
#ifdef SPLIT
          (*Uav[m])[k][j][i] -= c1[k+1][i] * (*Uav[m])[k+1][j][i];
#endif
/*  The tridiagonal solver is complete at this point.                 */

/*   At this point, velocity components which exceed a threshold for  */
/* physically reasonable values are truncated, or the column is       */
/* column is flagged to be sent to a diagnostic reporting subroutine. */
          if (fabs((*U[m])[k+1][j][i]) > maxvel) {
            if (maywrite[m] > 0) dowrite[m][i] = 1;
            else (*U[m])[k+1][j][i] = (((*U[m])[k+1][j][i] >= 0) ?
                                       truncvel : -truncvel);
            ntrunc++;
          }
#ifdef SPLIT
          if (fabs((*Uav[m])[k+1][j][i]) > maxvel) {
            if (maywrite[m] > 0) dowrite[m][i] = 1;
            else (*Uav[m])[k+1][j][i] = (((*Uav[m])[k+1][j][i] >= 0) ?
                                         truncvel : -truncvel);
            ntrunc++;
          }
#endif
        }
      }
      for (i=X1;i<=nx;i++) if (do_i[i]) {
        if (fabs((*U[m])[0][j][i]) > maxvel) {
          if (maywrite[m] > 0) dowrite[m][i] = 1;
          else (*U[m])[0][j][i] = (((*U[m])[0][j][i] >= 0) ?
                                   truncvel : -truncvel);
          ntrunc++;
        }
#ifdef SPLIT
        if (fabs((*Uav[m])[0][j][i]) > maxvel) {
          if (maywrite[m] > 0) dowrite[m][i] = 1;
          else (*Uav[m])[0][j][i] = (((*Uav[m])[0][j][i] >= 0) ?
                                     truncvel : -truncvel);
          ntrunc++;
        }
#endif

/*   Here the diagnostic reporting subroutines may be called if       */
/* unphysically large values were found.                              */
        if ((m==0) && (dowrite[0][i]==1)) write_u_acc(i,j,u,h,a,hvel,stress[i]);
        if ((m==1) && (dowrite[1][i]==1)) write_v_acc(i,j,v,h,a,hvel,stress[i]);

      }

    }

    if (diag.du_dt_visc != NULL) for (k=0;k<NZ;k++) for (i=X1;i<=nx;i++)
      diag.du_dt_visc[k][j][i] = (u[k][j][i] - diag.du_dt_visc[k][j][i])*Idt;
    if (diag.dv_dt_visc != NULL) for (k=0;k<NZ;k++) for (i=X1;i<=nx;i++)
      diag.dv_dt_visc[k][j][i] = (v[k][j][i] - diag.dv_dt_visc[k][j][i])*Idt;

  }

/* Offer diagnostic fields for averaging.                             */
  if (query_averaging_enabled() > 0) {
    static int ref_du_dt = -1, ref_dv_dt = -1, ref_au = -1, ref_av = -1;

    if (ref_du_dt < 0) ref_du_dt = get_field_ref(&diag.du_dt_visc[0]);
    if (ref_dv_dt < 0) ref_dv_dt = get_field_ref(&diag.dv_dt_visc[0]);
    if (ref_au < 0) ref_au = get_field_ref(&diag.au_visc[0]);
    if (ref_av < 0) ref_av = get_field_ref(&diag.av_visc[0]);

    if (ref_du_dt > 0) average_field(&diag.du_dt_visc[0], ref_du_dt);
    if (ref_dv_dt > 0) average_field(&diag.dv_dt_visc[0], ref_dv_dt);
    if (ref_au > 0) average_field(&diag.au_visc[0], ref_au);
    if (ref_av > 0) average_field(&diag.av_visc[0], ref_av);
  }
}


void set_viscous_BBL(double u[NZ][NYMEM][NXMEM], double v[NZ][NYMEM][NXMEM],
         double h[NZ][NYMEM][NXMEM], struct thermo_var_ptrs tv) {
/*   The following subroutine calculates the thickness of the bottom  */
/* boundary layer and the viscosity within that layer.  A drag law is */
/* used, either linearized about an assumed bottom velocity or using  */
/* the actual near-bottom velocities combined with an assumed         */
/* unresolved velocity.  The bottom boundary layer thickness is       */
/* limited by a combination of stratification and rotation, as in the */
/* paper of Killworth and Edwards, JPO 1999.  It is not necessary to  */
/* calculate the thickness and viscosity every time step; instead     */
/* previous values may be used.                                       */

/* Arguments: u - Zonal velocity, in m s-1.                           */
/*  (in)      v - Meridional velocity, in m s-1.                      */
/*  (in)      h - Layer thickness, in m.                              */
/*  (in)      tv - A structure containing pointers to any available   */
/*                 thermodynamic fields. Absent fields have NULL ptrs.*/
#ifdef BOTTOMDRAGLAW

  double htot;             /* Sum of the layer thicknesses up to some */
                           /* point, in m.                            */
  double Rhtot;            /* Running sum of thicknesses times the    */
                           /* layer potential densities in kg m-2.    */
  double h_at_vel[NZ][NXMEM];/* Layer thickness at a velocity point,  */
                           /* using an upwind-biased second order     */
                           /* accurate estimate based on the previous */
                           /* velocity direction, in m.               */
  double ustar;            /* The bottom friction velocity, in m s-1. */
  double ustarsq;          /* 400 times the square of ustar, times    */
                           /* RHO_0 divided by G_EARTH, in kg m-2.    */
  int do_i[NXMEM];
  int i, j, k, m;

  double cdrag_sqrt;       /* Square root of the drag coefficient, nd.*/
# ifndef LINEAR_DRAG
  double U_bg_sq;          /* The square of an assumed background     */
                           /* velocity, for calculating the mean      */
                           /* magnitude near the bottom for use in the*/
                           /* quadratic bottom drag, in m2.           */

  U_bg_sq = HIM_params.drag_bg_vel * HIM_params.drag_bg_vel;
# endif
  cdrag_sqrt=sqrt(HIM_params.cdrag);

# ifdef LINEAR_DRAG
/*  With a linear drag law, the friction velocity is already known.   */
  ustar = cdrag_sqrt*HIM_params.drag_bg_vel;
# endif

  for (j=Y1;j<=ny;j++) {
    for (m=0; m<=1; m++) {

      if (m==0) for (i=X1;i<=nx;i++) do_i[i] = ((umask[j][i]>0) > 0.0) ? 1 : 0;
      else for (i=X1;i<=nx;i++) do_i[i] = ((vmask[j][i]>0) > 0.0) ? 1 : 0;

      for (k=0;k<NZ;k++) {
        if (m==0) {
          for (i=X1;i<=nx;i++) if (do_i[i]) {
            h_at_vel[k][i] = ((u[k][j][i] * (h[k][j][i+1] - h[k][j][i]) >= 0) ?
                2.0 * h[k][j][i] * h[k][j][i+1] / (h[k][j][i] + h[k][j][i+1]) :
                0.5 * (h[k][j][i] + h[k][j][i+1]));
          }
        } else {
          for (i=X1;i<=nx;i++) if (do_i[i]) {
            h_at_vel[k][i] = ((v[k][j][i] * (h[k][j+1][i] - h[k][j][i]) >= 0) ?
                2.0 * h[k][j][i] * h[k][j+1][i] / (h[k][j][i] + h[k][j+1][i]) :
                0.5 * (h[k][j][i] + h[k][j+1][i]));
          }
        }
      }

      for (i=X1;i<=nx;i++) if (do_i[i]) {
# ifndef LINEAR_DRAG
/*   This block of code calculates the mean velocity magnitude over   */
/* the bottommost HIM_params.Hbbl of the water column for determining */
/* the quadratic bottom drag.                                         */
       {
        double hwtot;      /* Sum of the thicknesses used to calculate*/
                           /* the near-bottom velocity magnitude, m.  */
        double hutot;      /* Running sum of thicknesses times the    */
                           /* velocity magnitudes, in m2 s-1.         */
        htot  = 0.0; hwtot = 0.0; hutot = 0.0;
        for (k=NZ-1;k>=0;k--) {
          double hweight;          /* The thickness of a layer that   */
                                   /* is within Hbbl of the bottom.   */
          double v_at_u, u_at_v;   /* v at a u point or vice versa.   */

          if (htot>=HIM_params.Hbbl) break;
          else if (htot + h_at_vel[k][i] > HIM_params.Hbbl)
            hweight = HIM_params.Hbbl - htot;
          else hweight = h_at_vel[k][i];
          if (hweight < 1.5*EPSILON) continue;

          htot  += h_at_vel[k][i];  hwtot += hweight;

          if (m==0) {
            v_at_u = (hweight == 0.0) ? 0.0 :
               ((v[k][j][i] * (h[k][j][i]+h[k][j+1][i]) +
                 v[k][j][i+1] * (h[k][j][i+1]+h[k][j+1][i+1]) +
                 v[k][j-1][i+1] * (h[k][j-1][i+1]+h[k][j][i+1]) +
                 v[k][j-1][i] * (h[k][j-1][i]+h[k][j][i])) /
                (h[k][j-1][i] + 2.0*h[k][j][i] + h[k][j+1][i] +
                 h[k][j-1][i+1] + 2.0*h[k][j][i+1] + h[k][j+1][i+1]));

            hutot += hweight * sqrt(u[k][j][i]*u[k][j][i] +
	                            v_at_u*v_at_u + U_bg_sq);
          } else {
            u_at_v = (hweight == 0.0) ? 0.0 :
               ((u[k][j][i] * (h[k][j][i]+h[k][j][i+1]) +
                 u[k][j][i-1] * (h[k][j][i-1]+h[k][j][i]) +
                 u[k][j+1][i-1] * (h[k][j+1][i-1]+h[k][j+1][i]) +
                 u[k][j+1][i] * (h[k][j+1][i]+h[k][j+1][i+1])) /
                (h[k][j][i-1] + 2.0*h[k][j][i] + h[k][j][i+1] +
                 h[k][j+1][i-1] + 2.0*h[k][j+1][i] + h[k][j+1][i+1]));

            hutot += hweight * sqrt(v[k][j][i]*v[k][j][i] +
                                    u_at_v*u_at_v + U_bg_sq);
          }
        }
        ustar = hutot*cdrag_sqrt/hwtot;
       }
# endif

/*  The 400.0 in this expression is the square of a constant proposed */
/*  by Killworth and Edwards, 1999, in equation (2.20).               */
        ustarsq = 400.0*ustar*ustar*RHO_0/G_EARTH;
        Rhtot = 0.0; htot = 0.0;

/*   This block of code calculates the thickness of a stratification  */
/* limited bottom boundary layer, using the prescription from         */
/* Killworth and Edwards, 1999, as described in Stephens and Hallberg */
/* 2000.                                                              */
# ifdef BULKMIXEDLAYER
        for (k=NZ-1;k>=NML+1;k--)
# else
        for (k=NZ-1;k>=1;k--)
# endif
        {
          double oldfn;    /* The integrated energy required to       */
                           /* entrain up to the bottom of the layer,  */
                           /* divided by G_EARTH, in kg m-2.          */
          double Dfn;      /* The increment in oldfn for entraining   */
                           /* the layer, in kg m-1.                   */
          double Dh;       /* The increment in layer thickness from   */
                           /* the present layer, in m.                */

          oldfn = Rhtot - Rlay[k]*htot;
          Dfn = (Rlay[k] - Rlay[k-1])*(h_at_vel[k][i]+htot);

          if (oldfn >= ustarsq) break;
          else if ((oldfn + Dfn) <= ustarsq) Dh = h_at_vel[k][i];
          else Dh = h_at_vel[k][i] * sqrt((ustarsq-oldfn)/Dfn);

          htot += Dh;
          Rhtot += Rlay[k]*Dh;
        }
# ifndef BULKMIXEDLAYER              /* BULKMIXEDLAYER _Not_ defined. */
        if (Rhtot - Rlay[0]*htot < ustarsq) htot += h_at_vel[0][i];
# else                                     /* BULKMIXEDLAYER defined. */
        {
          double Rl, Rla;
          Rla = 0.5*(tv.Rml[NML][j][i] +
                     ((m==0) ? tv.Rml[NML][j][i+1] : tv.Rml[NML][j+1][i]));
          for (k=NML;k>=1;k--) {
            double oldfn, Dfn, Dh;  /* These are the same as above.   */
            Rl = Rla;
            Rla = 0.5*(tv.Rml[k-1][j][i] +
                       ((m==0) ? tv.Rml[k-1][j][i+1] : tv.Rml[k-1][j+1][i]));

            oldfn = Rhtot - Rl*htot;  Dfn = (Rl - Rla)*(h_at_vel[k][i]+htot);

            if (oldfn >= ustarsq) break;
            else if ((oldfn + Dfn) <= ustarsq) Dh = h_at_vel[k][i];
            else Dh = h_at_vel[k][i] * sqrt((ustarsq-oldfn)/Dfn);

            htot += Dh;   Rhtot += Rl*Dh;
          }
          if (Rhtot - Rla*htot < ustarsq) htot += h_at_vel[0][i];
        }
# endif                                        /* end BULKMIXEDLAYER. */

/* The Coriolis limit is 0.5*ustar/f. The buoyancy limit here is htot.*/
/* The  bottom boundary layer thickness is found by solving the same  */
/* equation as in Killworth and Edwards:    (h/h_f)^2 + h/h_N = 1.    */
        {
          double C2f;             /* C2f = 2*f at velocity points.    */
          C2f = (m==0) ? (f[j-1][i]+f[j][i]) : (f[j][i-1]+f[j][i]);
          bbl_thick[m][j][i] = htot /
                 (0.5 + sqrt(0.25 + htot*htot*C2f*C2f/(ustar*ustar)));
        }
        if (bbl_thick[m][j][i] < HIM_params.BBL_thick_min)
            bbl_thick[m][j][i] = HIM_params.BBL_thick_min;
# ifdef RINOMIX
/* If there is Richardson number dependent mixing, that determines    */
/* the vertical extent of the bottom boundary layer, and there is no  */
/* need to set that scale here.  In fact, viscously reducing the      */
/* shears over an excessively large region reduces the efficacy of    */
/* the Richardson number dependent mixing.                            */
        if (bbl_thick[m][j][i] > 0.5*HIM_params.Hbbl)
            bbl_thick[m][j][i] = 0.5*HIM_params.Hbbl;
# endif

/*   Here the near-bottom viscosity is set to a value which will give */
/* the correct stress when the shear occurs over bbl_thick.           */
        kv_bbl[m][j][i] = cdrag_sqrt*ustar*bbl_thick[m][j][i];
      }
    }
  }
#endif                                               /* BOTTOMDRAGLAW */

}

