
!=============================================================================
!
! $Id: smv2chem_diag.F90,v 1.4 2013-09-03 16:16:43 jkouatch Exp $
!
! CODE DEVELOPER
!   Original code from Peter Connell, LLNL
!   Gmimod modifications:  John Tannahill
!                          jrt@llnl.gov
!
! FILE
!   smv2chem_diag.F
!
! ROUTINES
!   Do_Smv2_Diag
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Smv2_Diag
!
! DESCRIPTION
!   This routine collects the Smvgear II chemical diagnostics.  It accumulates
!   averages for species and reaction rates.
!
! ARGUMENTS
!   jlooplo      : low ntloop grid-cell - 1 in a grid-block
!   ktloop       : number of grid-cells     in a grid-block
!   pr_nc_period : NetCDF output period
!   tdt          : model time step (s)
!   told         : stores last value of xelaps in case current step fails
!   do_cell_chem : do chemistry for a particular cell?
!   jreorder     : gives original grid-cell from re-ordered grid-cell
!   inewold      : original spc # of each new jnew spc
!   denair       : density of air (molec/cm^3)
!   cnew         : init (and final) spc conc
!                  (# cm^-3-air or moles l^-1-h2o (?))
!   xtimestep    : xelaps - told
!
!-----------------------------------------------------------------------------

      subroutine Do_Smv2_Diag  &
     &  (savedVars, jlooplo, ktloop, pr_nc_period, tdt, told, do_cell_chem,  &
     &   jreorder, inewold, denair, cnew, xtimestep, &
     &   yda, qqkda, qqjda, arateOrg, prateOrg, &
     &   ilong, ilat, ivert, itloop, &
     &   i1, i2, ju1, j2, k1, k2, &
     &   num_qks, num_qjs, num_active)

      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

      IMPLICIT NONE

#     include "smv2chem_par.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: ilong, ilat, ivert, itloop
      integer, intent(in) :: num_qks, num_qjs, num_active
      real*8 , intent(inout) :: qqjda(i1:i2, ju1:j2, k1:k2, num_qjs)
      real*8 , intent(inout) :: qqkda(i1:i2, ju1:j2, k1:k2, num_qks)
      real*8 , intent(inout) :: yda  (i1:i2, ju1:j2, k1:k2, num_active)
      real*8 , intent(in) :: arateOrg (ktloop, num_qks)
      real*8 , intent(in) :: prateOrg (ktloop, num_qjs)


      integer, intent(in)  :: jlooplo
      integer, intent(in)  :: ktloop
      real*8,  intent(in)  :: pr_nc_period
      real*8,  intent(in)  :: tdt
      real*8,  intent(in)  :: told
      logical, intent(in)  :: do_cell_chem(ilong, ilat, ivert)
      integer, intent(in)  :: jreorder(itloop)
      integer, intent(in)  :: inewold (MXGSAER, ICS)
      real*8,  intent(in)  :: denair  (ktloop)
      real*8,  intent(in)  :: cnew    (KBLOOP, MXGSAER)

      real*8,  intent(inout) :: xtimestep

      type(t_Smv2Saved), intent(inOut) :: savedVars

!     ----------------------
!     Variable declarations.
!     ----------------------

      logical :: end_chem_tstep

      integer :: ic, ix
      integer :: il, ij, ik
      integer :: ind
      integer :: jloop, kloop
      integer :: smv2di

      real*8  :: rsteps_per_period
      real*8  :: tacc_tdt
      real*8  :: tfrac

      real*8  :: qjblock (ktloop, num_qjs)
      real*8  :: qqjblock(ktloop, num_qjs)

      real*8  :: qkblock (ktloop, num_qks)
      real*8  :: qqkblock(ktloop, num_qks)

      real*8  :: yblock  (ktloop, num_active)


!     ----------------
!     Begin execution.
!     ----------------

!     ===============
      if (savedVars%end_period) then
!     ===============

        savedVars%end_period = .false.

!       -------------------------------------
!       Initialize arrays at start of period.
!       -------------------------------------

        savedVars%isteps (:,:,:) = 0

        savedVars%taccum (:,:,:) = 0.0d0

        qqjda(:,:,:,:) = 0.0d0
        qqkda(:,:,:,:) = 0.0d0
        yda  (:,:,:,:) = 0.0d0

        savedVars%qqjts(:,:,:,:) = 0.0d0
        savedVars%qqkts(:,:,:,:) = 0.0d0
        savedVars%yts  (:,:,:,:) = 0.0d0

      end if

        qjblock(:,:) = prateOrg(:,:)
        qkblock(:,:) = arateOrg(:,:)

      do kloop = 1, ktloop

        jloop  = jlooplo + kloop
        smv2di = jreorder(jloop)

        il = savedVars%lonloop(smv2di)
        ij = savedVars%latloop(smv2di)
        ik = savedVars%altloop(smv2di)

!        qjblock(kloop,:) = qjgmi(il,ij,ik,:)
!        qkblock(kloop,:) = qkgmi(il,ij,ik,:)

        do ic = 1, num_active
          ix = inewold(ic,1)
          yblock(kloop,ix) = cnew(kloop,ic)
        end do

      end do


!     =====================
      call Calc_Rate_Setkin  &
!     =====================
     &  (ktloop, num_qjs, num_qks, num_active, qkblock, qjblock,  &
     &   yblock, qqkblock, qqjblock)


!     -------------------------------------------------
!     Test for complete restart or time step overshoot.
!     -------------------------------------------------

      jloop  = jlooplo + 1
      smv2di = jreorder(jloop)

      il = savedVars%lonloop(smv2di) - i1  + 1
      ij = savedVars%latloop(smv2di) - ju1 + 1
      ik = savedVars%altloop(smv2di)

      tacc_tdt = savedVars%taccum(il,ij,ik)


      if ((told == 0.0d0) .and. (tacc_tdt /= 0.0d0)) then

!       ------------------------------------------------
!       SmvgearII has restarted timestep for this block.
!       ------------------------------------------------

        tacc_tdt = 0.0d0

        do kloop = 1, ktloop

          jloop  = jlooplo + kloop
          smv2di = jreorder(jloop)

          il = savedVars%lonloop(smv2di) - i1  + 1
          ij = savedVars%latloop(smv2di) - ju1 + 1
          ik = savedVars%altloop(smv2di)

          savedVars%taccum(il,ij,ik)   = 0.0d0

          savedVars%qqjts (il,ij,ik,:) = 0.0d0
          savedVars%qqkts (il,ij,ik,:) = 0.0d0
          savedVars%yts   (il,ij,ik,:) = 0.0d0

        end do

      end if


      if ((Abs ((tacc_tdt + xtimestep) - tdt) <= 1.0d-06) .or.  &
     &    ((tacc_tdt + xtimestep) > tdt )) then

!       -----------------------------------------------
!       SmvgearII has overshot timestep for this block.
!       -----------------------------------------------

        xtimestep = tdt - tacc_tdt

      end if

      tfrac = xtimestep / pr_nc_period


!     ------------------------------------
!     Accumulate rates and concentrations.
!     ------------------------------------

      do kloop = 1, ktloop

        jloop  = jlooplo + kloop
        smv2di = jreorder(jloop)

        il = savedVars%lonloop(smv2di) - i1  + 1
        ij = savedVars%latloop(smv2di) - ju1 + 1
        ik = savedVars%altloop(smv2di)

        where (qqjblock(kloop,:) > 0.0d0) savedVars%qqjts(il,ij,ik,:)  =  &
     &    savedVars%qqjts(il,ij,ik,:) +  &
     &    ((qqjblock(kloop,:)  / denair(kloop)) * tfrac)

        where (qqkblock(kloop,:) > 0.0d0) savedVars%qqkts(il,ij,ik,:)  =  &
     &    savedVars%qqkts(il,ij,ik,:) +  &
     &    ((qqkblock(kloop,:)  / denair(kloop)) * tfrac)

        savedVars%yts(il,ij,ik,:) =  &
     &    savedVars%yts(il,ij,ik,:) +  &
     &    (yblock(kloop,:) * tfrac)

        savedVars%taccum(il,ij,ik) = savedVars%taccum(il,ij,ik) + xtimestep

      end do


      if (Any (do_cell_chem(:,:,:) .and.  &
     &         (tdt - savedVars%taccum(:,:,:) > 1.0d-06))) then
        end_chem_tstep = .false.
      else
        end_chem_tstep = .true.
      end if


!     ===================
      if (end_chem_tstep) then
!     ===================

        savedVars%nsteps = savedVars%nsteps + 1

        savedVars%taccum(:,:,:) = 0.0d0

!       ----------------------------------------------------------------
!       Keep count of number of time steps cells have calculated values.
!       ----------------------------------------------------------------

        where (do_cell_chem(:,:,:)) savedVars%isteps(:,:,:) = savedVars%isteps(:,:,:) + 1


        qqjda(i1:i2,ju1:j2,k1:k2,:) =  &
     &    qqjda(i1:i2,ju1:j2,k1:k2,:) + savedVars%qqjts(1:ilong,1:ilat,1:ivert,:)

        qqkda(i1:i2,ju1:j2,k1:k2,:) =  &
     &    qqkda(i1:i2,ju1:j2,k1:k2,:) + savedVars%qqkts(1:ilong,1:ilat,1:ivert,:)

        yda  (i1:i2,ju1:j2,k1:k2,:) =  &
     &    yda  (i1:i2,ju1:j2,k1:k2,:) + savedVars%yts  (1:ilong,1:ilat,1:ivert,:)


        savedVars%qqjts(:,:,:,:) = 0.0d0
        savedVars%qqkts(:,:,:,:) = 0.0d0
        savedVars%yts  (:,:,:,:) = 0.0d0

      end if


!     ================================
      if (savedVars%nsteps == savedVars%nsteps_per_period) then
!     ================================

        savedVars%end_period = .true.

        savedVars%nsteps = 0

        rsteps_per_period = savedVars%nsteps_per_period


        do ic = 1, num_qjs
          where (do_cell_chem(1:ilong,1:ilat,1:ivert))

            qqjda(i1:i2,ju1:j2,k1:k2,ic) =  &
     &        qqjda(i1:i2,ju1:j2,k1:k2,ic) *  &
     &        (rsteps_per_period / savedVars%isteps(1:ilong,1:ilat,1:ivert))

          end where
        end do

        do ic = 1, num_qks
          where (do_cell_chem(1:ilong,1:ilat,1:ivert))

            qqkda(i1:i2,ju1:j2,k1:k2,ic) =  &
     &        qqkda(i1:i2,ju1:j2,k1:k2,ic) *  &
     &        (rsteps_per_period / savedVars%isteps(1:ilong,1:ilat,1:ivert))

          end where
        end do

        do ic = 1, num_active
          where (do_cell_chem(1:ilong,1:ilat,1:ivert))

            yda(i1:i2,ju1:j2,k1:k2,ic) =  &
     &        yda(i1:i2,ju1:j2,k1:k2,ic) *  &
     &        (rsteps_per_period / savedVars%isteps(1:ilong,1:ilat,1:ivert))

!c?       elsewhere

!c          yda(:,:,:,ic) = yinit(:,:,:,ic)

          end where
        end do

      end if


      return

      end

