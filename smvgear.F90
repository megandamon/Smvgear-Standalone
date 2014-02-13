
!=============================================================================
!
! $Id: smvgear.F90,v 1.5 2013-09-03 16:16:43 jkouatch Exp $
!
! CODE DEVELOPER
!   Original code from Mark Z. Jacobson ((C) COPYRIGHT, 1993).
!   LLNL modifications:  John Tannahill
!                        jrt@llnl.gov
!
! FILE
!   smvgear.F
!
! ROUTINES
!   Backsub
!   Smvgear
!   Decomp
!   Update
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Smvgear
!
! DESCRIPTION
!   This routine is the driver for the Smvgear (Sparse Matrix Vector Gear code)
!   chemistry solver.  It uses a Gear-type integrator that solves first order
!   ordinary differential equations with initial value boundary conditions.
!   Smvgear differs from an original Gear code in that it uses sparse matrix
!   and vectorization techniques to improve its computational speed.
!
!   This version is Smvgear II, 9/96.  It has been modified to include
!   grid-cell reordering prior to each time interval and different chemistry
!   for different atmospheric regions.  The purpose of the reordering is to
!   group cells with stiff equations together and those with non-stiff
!   equations together.  This reordering can save signifcant computer time
!   (e.g., speed the code by a factor of two or more), depending on the
!   variation in stiffness throughout the grid-domain.  When the stiffness is
!   the same throughout the grid-domain (e.g., if all concentrations and rates
!   are the same), then reordering is unnecessary and will not speed solutions.
!
!   This version includes a variable absolute error tolerance.  The absolute
!   tolerance is recalculated every few Gear time steps.  This version also
!   contains different sets of chemistry for different regions of the
!   atmosphere.  Thus, urban, free tropospheric, and stratospheric chemistry
!   can be solved during the same model run.
!
!   References =>
!   ----------
!
!     Jacobson M. Z. (1997) Improvement in Smvgear II through Absolute
!     Error Tolerance Control; in submission.
!
!     Jacobson M. Z. (1995) Computation of Global Photochemistry with Smvgear
!     II, Atmos. Environ., 29a, 2541-2546.
!
!     Jacobson M. Z. (1994) Developing, Coupling, and Applying a Gas, Aerosol,
!     Transport, and Radiation Model to Studying Urban and Regional Air
!     Pollution, PhD thesis, University of California, Los Angeles.
!
!     Jacobson M. Z. and Turco R. P. (1994) Smvgear: A Sparse Matrix,
!     Vectorized Gear Code for Atmospheric Models, Atmos. Environ. 28a,
!     273-284.
!
!     The origins of the Gear integrator used in Smvgear are found in:
!       Gear C. W. (1971) Numerical Initial Value Problems in Ordinary
!       Differential Equations, Prentice-Hall, NJ, pp. 158-166.
!
!     Finally, in subroutine Smvgear, the following ideas originated from
!     Lsodes, the Livermore solver for ordinary differential with sparse
!     matrices (Hindmarsh A. C. and Sherman A. H.):
!       (a) predicting the first time-step;
!       (b) determining corrector convergence differently than in Gear's
!           original code (goc);
!       (c) determining error differently than in goc;
!       (d) summing up the pascal matrix differently than in goc.
!
!     References for the 1987 Lsodes version include:
!
!       Sherman A. H. and Hindmarsh A. C. (1980) Gears: A Package for the
!       Solution of Sparse, Stiff Ordinary Differential Equations, Lawrence
!       Livermore National Laboratory Report UCRL-84102.
!
!       Hindmarsh A. C. (1983) Odepack, A Systematized Collection of ODE
!       Solvers, Scientific Computing, R.S. Stepleman et. al., eds.,
!       North-Holland, Amsterdam, pp. 55-74.
!
! ARGUMENTS
!   do_qqjk_inchem   : if pr_qqjk is on, should qqj's & qqk's be determined
!                      inside the chemistry solver, or outside?
!   do_semiss_inchem : do surface emissions inside the chemistry solver, or
!                      outside?
!   pr_qqjk  : should the periodic qqjk output file be written?
!   pr_smv2  : should the SmvgearII     output file be written
!              (non-parallel mode only)?
!   ifsun    : identifies whether sun is up (=1) or down (=2)
!   ilat     : # of latitudes
!   ilong    : # of longitudes
!   ivert    : # of vertical layers
!   ireord   : 1 => reorder grid-cells and blocks for chemistry
!              2 => solve chemistry
!   itloop   : # of zones (ilong * ilat * ivert)
!   jlooplo  : low ntloop grid-cell - 1 in a grid-block
!   ktloop   : # of grid-cells in a grid-block
!   lunsmv   : logical unit number to write to when pr_smv2 is true
!   numActiveReactants    : # of active rxns
!   ncs      : identifies gas chemistry type (1.NCSGAS)
!   nfdh2    : nfdh3 + # of rxns with two   active reactants
!   nfdh3    :         # of rxns with three active reactants
!   nfdl1    : nfdh2 + 1
!   nfdl2    : nfdh3 + 1
!   nfdrep   : nfdh3 + # of rxns with two active reactants that are not
!              followed by a rxn with the same reactants
!   nfdrep1  : nfdrep + 1
!   fracdec  : fraction time step is decreased in Smvgear if convergence
!              test fails
!   hmaxnit  : max time step for night all chem (s)
!   pr_nc_period     : NetCDF output period
!   tdt      : model time step (s)
!   do_cell_chem     : do chemistry for a particular cell?
!   irma,b,c : spc # of each reactant; locates reordered active spc #s
!   jreorder : gives original grid-cell from re-ordered grid-cell
!   jphotrat : tbd
!   ntspec   : # of active + inactive gases
!   inewold  : original spc # of each new jnew spc
!   denair   : density of air (molec/cm^3)
!   corig    : original gas-phase concs used to restart Smvgear if a failure
!              occurs (molec/cm^3)
!   pratk1   : tbd
!   yemis    : surface emissions (units?)
!   smvdm    : amount added to each spc at each grid-cell, for mass balance
!              accounting (# cm^-3 for gas chemistry (?))
!   nfdh1    : nfdh2 + # of rxns with one   active reactant
!   errmx2   : measure of stiffness/nearness to convergence of each block
!              sum ydot/y for all species (MRD per Kareem Sorathia)
!   cc2      : array holding values of decomposed matrix
!   cnew     : stores cnewDerivatives (y (estimated)) (molec/cm^3)
!   gloss    : value of first derivatives on output from velocity; right-side
!              of eqn on input to Backsub; error term (solution from Backsub)
!              on output from Backsub
!   vdiag    : 1 / current diagonal term of the decomposed matrix
!   rrate    : rate constants
!   trate    : rxn rate (moles l^-1-h2o s^-1 or # cm^-3 s^-1 (?)) !REMOVED!
!
!-----------------------------------------------------------------------------

      subroutine Smvgear  &
     &  (savedVars, do_qqjk_inchem, do_semiss_inchem, pr_qqjk, pr_smv2, ifsun,  &
     &   ilat, ilong, ivert, ireord, itloop, jlooplo, ktloop, lunsmv,  &
     &   numActiveReactants, ncs, nfdh2, nfdh3, nfdl1, nfdl2, nfdrep, nfdrep1,  &
     &   fracdec, hmaxnit, pr_nc_period, tdt, do_cell_chem, irma, irmb,  &
     &   irmc, jreorder, jphotrat, ntspec, inewold, denair, corig,  &
     &   pratk1, yemis, smvdm, nfdh1, errmx2, cc2, cnew, gloss, vdiag,  &
     &   rrate, &
     &   yda, qqkda, qqjda, &
     &   CTMi1, CTMi2, CTMju1, CTMj2, CTMk1, CTMk2, &
     &   num_qks, num_qjs, num_active)

      use GmiPrintError_mod, only : GmiPrintError
      use GmiMechanism_mod
      use GmiManager_mod
      use GmiSparseMatrix_mod
      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

      implicit none

#     include "smv2chem_par.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: CTMi1, CTMi2, CTMju1, CTMj2, CTMk1, CTMk2
      integer, intent(in) :: num_qks, num_qjs, num_active
      real*8 , intent(inout) :: qqjda(CTMi1:CTMi2, CTMju1:CTMj2, CTMk1:CTMk2, num_qjs)
      real*8 , intent(inout) :: qqkda(CTMi1:CTMi2, CTMju1:CTMj2, CTMk1:CTMk2, num_qks)
      real*8 , intent(inout) :: yda  (CTMi1:CTMi2, CTMju1:CTMj2, CTMk1:CTMk2, num_active)
      logical, intent(in)  :: do_qqjk_inchem
      logical, intent(in)  :: do_semiss_inchem
      logical, intent(in)  :: pr_qqjk
      logical, intent(in)  :: pr_smv2
      integer, intent(in)  :: ifsun
      integer, intent(in)  :: ilat, ilong, ivert
      integer, intent(in)  :: ireord
      integer, intent(in)  :: itloop
      integer, intent(in)  :: jlooplo
      integer, intent(in)  :: ktloop
      integer, intent(in)  :: lunsmv
      integer, intent(in)  :: numActiveReactants
      integer, intent(in)  :: ncs
      integer, intent(in)  :: nfdh2,  nfdh3
      integer, intent(in)  :: nfdl1,  nfdl2
      integer, intent(in)  :: nfdrep, nfdrep1
      real*8,  intent(in)  :: fracdec
      real*8,  intent(in)  :: hmaxnit
      real*8,  intent(in)  :: pr_nc_period
      real*8,  intent(in)  :: tdt
      logical, intent(in)  :: do_cell_chem(ilong, ilat, ivert)
      !K: irm[a,b,c] can be removed as they are now in GenChem and don't vary with block
      integer, intent(in)  :: irma    (NMTRATE)
      integer, intent(in)  :: irmb    (NMTRATE)
      integer, intent(in)  :: irmc    (NMTRATE)
      integer, intent(in)  :: jreorder(itloop)
      integer, intent(in)  :: jphotrat(ICS)
      integer, intent(in)  :: ntspec  (ICS)
      integer, intent(in)  :: inewold (MXGSAER, ICS)
      real*8,  intent(in)  :: denair  (ktloop)
      real*8,  intent(in)  :: corig   (KBLOOP, MXGSAER)
      real*8,  intent(in)  :: pratk1  (KBLOOP, IPHOT)
      real*8,  intent(in)  :: yemis   (ilat*ilong, IGAS)

      real*8,  intent(inout) :: errmx2(itloop)
      !K: Why is cc2 passed in/out?  Seems silly
      !K: Same question for everything but cnew
      real*8,  intent(inout) :: cc2   (KBLOOP, 0:MXARRAY)
      real*8,  intent(inout) :: cnew  (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: gloss (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: smvdm (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: vdiag (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: rrate (KBLOOP, NMTRATE)

      integer, intent(out) :: nfdh1

      type(t_Smv2Saved), intent(inOut) :: savedVars

!     ----------------------
!     Variable declarations.
!     ----------------------

!     ------------------------------------------------------------------------
!
!     jeval     :  1 => call Pderiv the next time through the corrector steps;
!                  0 => last step successful and do not need to call Pderiv;
!                 -1 => Pderiv just called, and do not need to call again
!                  until jeval switched to 1
!     ------------------------------------------------------------------------

      integer :: i, j, k
      integer :: i1, i2
      integer :: jb
      integer :: jeval
      integer :: jg1
      integer :: jgas
      integer :: jnew
      integer :: jspc
      integer :: k1, k2, k3, k4, k5
      integer :: kloop
      integer :: kstepisc
      integer :: l3
      integer :: nact
      integer :: ncsp  ! ncs       => for daytime   gas chemistry
                       ! ncs + ICS => for nighttime gas chemistry
      integer :: nqisc
      integer :: ibcb(IGAS)

!     ------------------------------------------------
!     concAboveAbtolCount : counts # of concs above abtol(i), i = 1..
!     ------------------------------------------------
      integer :: concAboveAbtolCount(KBLOOP, 5)

!     ------------------------------------------------------------------------
!     delt      : current time step (s)
!     MAX_REL_CHANGE     : max relative change in delt*aset(1) before Pderiv is called
!     order     : floating point value of num1stOEqnsSolve, the order of # of ODEs
!     ------------------------------------------------------------------------

      real*8  :: cnewylow
      real*8  :: cnw
      real*8  :: consmult
      real*8  :: delt
      real*8  :: dtasn1
      real*8  :: der1max, der3max
      real*8  :: errymax
      real*8  :: r1delt, rdelta
      real*8  :: real_kstep
      real*8  :: rmsErrorPrevious, rmsrat
      real*8  :: xtimestep

      real*8, parameter  :: MAX_REL_CHANGE = 0.3d0

      real*8 :: eup ! pertst^2*order for one order higher than current order
      real*8  :: edwn ! pertst^2*order for one order lower  than current order


!     -------------------------------------------------------------------------
!     dely   : tbd
!     yabst  : absolute error tolerance (molec/cm^-3 for gases)
!     cest   : stores value of dtlos when idoub = 1
!     explic : tbd
!     cnewDerivatives   : an array of length num1stOEqnsSolve*(MAXORD+1) that carries the
!              derivatives of cnew, scaled by delt^j/factorial(j), where j is
!              the jth derivative; j varies from 1 to nqq; e.g., cnewDerivatives(jspc,2)
!              stores delt*y' (estimated)
!     -------------------------------------------------------------------------

      real*8  :: dely  (KBLOOP)
      real*8  :: yabst (KBLOOP)
      real*8  :: cest  (KBLOOP, MXGSAER)
      real*8  :: explic(KBLOOP, MXGSAER)
      real*8  :: cnewDerivatives  (KBLOOP, MXGSAER*7)

      type (Mechanism_type) :: mechanismObject
      type (Manager_type) :: managerObject
      integer :: nondiag     ! # of final matrix positions, excluding diagonal

      call initializeMechanism (mechanismObject, ktloop, irma, &
                              &  irmb, irmc, nfdh2, nfdh3, nfdrep, rrate)

      nact = nnact

!     =======================
#     include "setkin_ibcb.h"
!     =======================

      call resetGear (managerObject, ncsp, ncs, ifsun, hmaxnit, savedVars)

!     ----------------------------------------------------
!     Start time interval or re-enter after total failure.
!     ----------------------------------------------------

!     ========
 100  continue
!     ========

      call startTimeInterval (managerObject, ncs, savedVars)
      call initConcentrationArray(ktloop, cnew, corig, managerObject)

!     --------------------------------------------------------------------
!     Re-enter here if total failure or if restarting with new cell block.
!     --------------------------------------------------------------------

!     ========
 150  continue
!     ========

      call resetBeforeUpdate (managerObject)

!!DIR$ INLINE
      call updatePhotoDissRates  (mechanismObject, ktloop, numActiveReactants, &
         &  ncs, ncsp, jphotrat, pratk1, savedVars)
!!DIR$ NOINLINE

      call velocity (mechanismObject, managerObject%num1stOEqnsSolve, ncsp, cnew, gloss, nfdh1, savedVars)

      managerObject%numCallsVelocity = managerObject%numCallsVelocity + 1

      call setBoundaryConditions (mechanismObject, itloop, &
         & jreorder, jlooplo, ilat, &
         & ilong, ntspec, ncs, inewold, &
         & do_semiss_inchem, gloss, yemis, ibcb)


! MRD: Take the reordering and error/tolerance calculations and keep them in the solver for now
!     -------------------------------------------
!     Determine initial absolute error tolerance.
!     -------------------------------------------

      !do kloop = 1, ktloop
      !  dely(kloop) = 0.0d0
      !end do

! get rid of magic number 1
! refactor this into routine(s) smvgear until we deterine their resting place.
! probably goes in the manager, possibly in gear
!     ==========================
      IREORDIF: if (ireord /= 1) then
!     ==========================

      call calcNewAbsoluteErrorTolerance (managerObject, cnew, concAboveAbtolCount, &
         & ktloop, yabst, ncs, savedVars)

      !MRD: see manager routine "calculateErrorTolerances"
        do kloop = 1, ktloop
          dely(kloop) = 0.0d0
          do jspc = 1, managerObject%num1stOEqnsSolve
            cnewylow    = cnew (kloop,jspc) + (yabst(kloop) * managerObject%reltol1)
            errymax     = gloss(kloop,jspc) / cnewylow
            dely(kloop) = dely (kloop) + (errymax * errymax) ! this is an error (not a tolerance)
          end do
        end do

!     ====
      else
!     ====

      call calculateErrorTolerances (managerObject, ktloop, jlooplo, itloop, cnew, gloss, dely, errmx2)
      return

!     ===============
      end if IREORDIF
!     ===============


      call calcInitialTimeStepSize (managerObject, ktloop, dely, &
         delt, ncs, savedVars)

!     -----------------------
!     Set initial order to 1.
!     -----------------------

      managerObject%nqqold = 0
      managerObject%nqq    = 1
      jeval  = 1
      managerObject%rdelt  = 1.0d0


!     --------------------------------------------------------------
!     Store initial concentration and first derivatives x time step.
!     --------------------------------------------------------------

      do jspc = 1, managerObject%num1stOEqnsSolve
         j = jspc + managerObject%num1stOEqnsSolve

        do kloop = 1, ktloop
          cnewDerivatives(kloop,jspc) = cnew(kloop,jspc)
          cnewDerivatives(kloop,j)    = delt * gloss(kloop,jspc)
        end do

      end do


!     ========
 200  continue
!     ========

      if (managerObject%nqq /= managerObject%nqqold) call updateCoefficients (managerObject, savedVars)
      call calculateTimeStep (managerObject, delt, jeval, MAX_REL_CHANGE)


      if (delt < HMIN) then
        call tightenErrorTolerance (managerObject, pr_smv2, &
               lunsmv, ncs, delt, savedVars)
        !     ========================================
        go to 100 ! routine start startTimeInterval?
        !     ========================================
      end if


!     -------------------------------------------------------------------
!     If the delt is different than during the last step (if rdelt /= 1),
!     then scale the derivatives.
!     -------------------------------------------------------------------

      if (managerObject%rdelt /= 1.0d0) then
         call scaleDerivatives (managerObject, ktloop, cnewDerivatives)
      end if


!     --------------------------------------------------------------
!     If the last step was successful, reset rdelmax = 10 and update
!     the chold array with current values of cnew.
!     --------------------------------------------------------------

!     ================================
      IFSUCCESSIF: if (managerObject%ifsuccess == 1) then
!     ================================

        managerObject%rdelmax = 10.0d0

!       ---------------------------------------
!       Determine new absolute error tolerance.
!       ---------------------------------------

        if (Mod (managerObject%numSuccessTdt, 3) == 2) then

            call calcNewAbsoluteErrorTolerance (managerObject, cnew, concAboveAbtolCount, &
               &  ktloop, yabst, ncs, savedVars)

        end if

         call updateChold (managerObject, ktloop, cnew, yabst)

!     ==================
      end if IFSUCCESSIF
!     ==================

      call predictConcAndDerivatives (managerObject, cnewDerivatives, explic, ktloop)


!     -------------------------------------------------------------------
!     Correction loop.
!
!     Take up to 3 corrector iterations.  Test convergence by requiring
!     that changes be less than the rms norm weighted by chold.
!     Accumulate the correction in the array dtlos.  It equals the
!     jth derivative of concentration multiplied by delt^kstep /
!     (factorial(kstep-1) * aset(kstep)); thus, it is proportional to the
!     actual errors to the lowest power of delt present (delt^kstep).
!     -------------------------------------------------------------------


!     ========
 250  continue
!     ========

      l3 = 0
      do jspc = 1, managerObject%num1stOEqnsSolve
        do kloop = 1, ktloop
          cnew (kloop,jspc) = cnewDerivatives(kloop,jspc)
          managerObject%dtlos(kloop,jspc) = 0.0d0
        end do
      end do


!     ------------------------------------------------------------------
!     If jeval = 1, re-evaluate predictor matrix P = I - H * aset(1) * J
!     before starting the corrector iteration.  After calling Pderiv,
!     set jeval = -1 to prevent recalling Pderiv unless necessary later.
!     Call Decomp to decompose the matrix.
!     ------------------------------------------------------------------

      if (jeval == 1) then

         r1delt = -managerObject%asn1 * delt
         nondiag  = savedVars%iarray(ncsp) - managerObject%num1stOEqnsSolve ! iarray is in common block
         mechanismObject%numRxns1 = nfdh2 + savedVars%ioner(ncsp)

         !K: Need to send whole mech object to get rrate in predictor, also need cnew
         call calculatePredictor (nondiag, savedVars%iarray(ncsp), mechanismObject, cnew, &
            &  savedVars%npdhi(ncsp), savedVars%npdlo(ncsp), r1delt, cc2, savedVars)
         managerObject%numCallsPredict  = managerObject%numCallsPredict + 1

         ! MRD: End block of code that was in Pderiv

         !K: Consider un-inlining this.
!!DIR$   INLINE
!       ===========
        call Decomp  &
!       ===========
     &    (savedVars, managerObject%num1stOEqnsSolve, ktloop, ncsp, cc2, vdiag)
!DIR$   NOINLINE

        jeval  = -1
        managerObject%hratio = 1.0d0
        managerObject%nslp   = managerObject%numSuccessTdt + MBETWEEN
        managerObject%drate  = 0.7d0

      end if


!     -------------------------------------------------------------
!     Evaluate the first derivative using corrected values of cnew.
!     -------------------------------------------------------------

!     ========
 300  continue
!     ========

!   print*, "in 300, evaluating first derivative"
     call velocity (mechanismObject, managerObject%num1stOEqnsSolve, ncsp, cnew, gloss, nfdh1, savedVars)

     managerObject%numCallsVelocity = managerObject%numCallsVelocity + 1

      call setBoundaryConditions (mechanismObject, itloop, &
         & jreorder, jlooplo, ilat, &
         & ilong, ntspec, ncs, inewold, &
         & do_semiss_inchem, gloss, yemis, ibcb)



!     ---------------------------------------------------------------
!     In the case of the chord method, compute error (gloss) from the
!     corrected calculation of the first derivative.
!     ---------------------------------------------------------------

      do jspc = 1, managerObject%num1stOEqnsSolve

        j = jspc + managerObject%num1stOEqnsSolve

        do kloop = 1, ktloop
          gloss(kloop,jspc) = (delt * gloss(kloop,jspc)) -  &
     &                        (cnewDerivatives(kloop,j) + managerObject%dtlos(kloop,jspc))
        end do

      end do


!     --------------------------------------------------------------
!     Solve the linear system of equations with the corrector error;
!     Backsub solves backsubstitution over matrix of partial derivs.
!     --------------------------------------------------------------
! MRD: stays in gear b/c manager doesn't know about backsubstitution
! As part of gear's timestep it calls backsub
!     ============
      call Backsub  &
!     ============
     &  (savedVars, managerObject%num1stOEqnsSolve, ktloop, ncsp, cc2, vdiag, gloss)


!     ----------------------------------------------------------------
!     Sum up the accumulated error, correct the concentration with the
!     error, and begin to calculate the rmsnorm of the error relative
!     to chold.
!     ----------------------------------------------------------------

      do kloop = 1, ktloop
        dely(kloop) = 0.0d0
      end do

      ! MRD: removed an optimization for the case of asn1 = 1  (saves a multiplication per loop)
      do i = 1, managerObject%num1stOEqnsSolve !*
         do kloop = 1, ktloop
            managerObject%dtlos(kloop,i) = managerObject%dtlos(kloop,i) + gloss(kloop,i) !*
            cnew(kloop,i)  = cnewDerivatives(kloop,i)  + (managerObject%asn1 * managerObject%dtlos(kloop,i))
            errymax        = gloss(kloop,i) * managerObject%chold(kloop,i) !*
            dely(kloop)    = dely(kloop)    + (errymax * errymax) !*
         end do
      end do


      call calculateNewRmsError (managerObject, ktloop, dely, l3, savedVars)


!     --------------------------------------------------------
!     If convergence occurs, go on to check accumulated error.
!     --------------------------------------------------------

      if (managerObject%dcon > 1.0d0) then

!       -------------------------------------------------------------------
!       If nonconvergence after one step, re-evaluate first derivative with
!       new values of cnew.
!       -------------------------------------------------------------------

        if (l3 == 1) then

!         =========
          go to 300
!         =========

!         ----------------------------------------------------------------
!         The corrector iteration failed to converge.
!
!         If the Jacobian matrix is more than one step old, update the
!         Jacobian and try convergence again.  If the Jacobian is current,
!         then reduce the time step, reset the accumulated derivatives to
!         their values before the failed step, and retry with the smaller
!         step.
!         ----------------------------------------------------------------

        else if (jeval == 0) then

          managerObject%numFailOldJacobian = managerObject%numFailOldJacobian + 1
          jeval = 1

!         =========
          go to 250
!         =========

        end if

        managerObject%numFailAfterPredict     = managerObject%numFailAfterPredict + 1
        managerObject%rdelmax   = 2.0d0
        jeval     = 1
        managerObject%ifsuccess = 0
        managerObject%xelaps    = managerObject%told
        managerObject%rdelt     = fracdec

        call resetCnewDerivatives(managerObject, cnewDerivatives, ktloop)

!       =========
        go to 200
!       =========

      end if


!     -------------------------------------------------------------------
!     The corrector iteration converged.
!
!     Set jeval = 0, so that it does not need to be called the next step.
!     If all else goes well.  Next, test the accumulated error from the
!     convergence process above.
!     -------------------------------------------------------------------

      jeval = 0
      if (l3 > 1) then
         call testAccumulatedError (managerObject, ktloop, dely)
      end if


!     ----------------------------------------------------------------
!     The accumulated error test failed.
!
!     In all cases, reset the derivatives to their values before the
!     last time step.  Next:
!       (a) re-estimate a time step at the same or one lower order and
!           retry the step;
!       (b) if the first attempts fail, retry the step at fracdec the
!           the prior step;
!       (c) iF this fails, reset the order to 1 and go back to the
!           beginning, at order = 1, because errors of the wrong order
!           have accumulated.
!     ----------------------------------------------------------------

!     ==============================
      DER2MAXIF: if (managerObject%der2max > managerObject%enqq) then
!     ==============================

        managerObject%xelaps = managerObject%told
        managerObject%numFailErrorTest  = managerObject%numFailErrorTest + 1
        managerObject%jfail  = managerObject%jfail  + 1

        call resetCnewDerivatives(managerObject, cnewDerivatives, ktloop)

        managerObject%rdelmax = 2.0d0

        ! MRD: magic numbers
        ! enums
        ! prefer strings to integers
        if (managerObject%jFail <= 6) then

          managerObject%ifsuccess = 0
          managerObject%rdeltup   = 0.0d0

!         =========
          go to 400
!         =========

        else if (managerObject%jfail <= 20) then

          managerObject%ifsuccess = 0
          managerObject%rdelt     = fracdec

!         =========
          go to 200
!         =========

        else

          delt    = delt * 0.1d0
          managerObject%rdelt   = 1.0d0
          managerObject%jfail   = 0
          managerObject%jrestar = managerObject%jrestar + 1
          managerObject%idoub   = 5

          do jspc = 1, managerObject%num1stOEqnsSolve
            do kloop = 1, ktloop
              cnew(kloop,jspc) = cnewDerivatives(kloop,jspc)
            end do
          end do

          if (pr_smv2) then
            Write (lunsmv,970) delt, managerObject%xelaps
          end if

 970      format ('delt dec to ', e13.5, ' at time ', e13.5,  &
     &            ' because of excessive errors.')

          if (managerObject%jrestar == 100) then

            if (pr_smv2) then
              Write (lunsmv,980)
            end if

 980        format ('Smvgear:  Stopping because of excessive errors.')
            print*, "Excessive errors"
            call GmiPrintError ('Problem in Smvgear', .true., 0, 0, 0, 0, 0.0d0, 0.0d0)

          end if

!         =========
          go to 150
!         =========

        end if

!     ====
      else
!     ====

!       -------------------------------------------------------------
!       All successful steps come through here.
!
!       After a successful step, update the concentration and all
!       derivatives, reset told, set ifsuccess = 1, increment numSuccessTdt,
!       and reset jfail = 0.
!       -------------------------------------------------------------


        if (pr_qqjk .and. do_qqjk_inchem) then
          xtimestep = managerObject%xelaps - managerObject%told

!      print*, "about to call Do_Smv2_Diag"
!         =================
          call Do_Smv2_Diag  &
!         =================
     &      (savedVars, jlooplo, ktloop, pr_nc_period, tdt, managerObject%told, &
             do_cell_chem, jreorder, inewold, denair, cnew, xtimestep, &
     &       yda, qqkda, qqjda, rrate, pratk1, &
     &       ilong, ilat, ivert, itloop, &
     &       CTMi1, CTMi2, CTMju1, CTMj2, CTMk1, CTMk2, &
     &       num_qks, num_qjs, num_active)
        end if


        managerObject%jfail     = 0
        managerObject%ifsuccess = 1
        managerObject%numSuccessTdt    = managerObject%numSuccessTdt + 1
        managerObject%told      = managerObject%xelaps

        call updateDerivatives(managerObject, cnewDerivatives, ktloop, savedVars)

         ! TODO: Megan FIX this routine. Does not give 0 diff
        !call updateChemistryMassBalance (ktloop, cnewDerivatives, explic, smvdm, &
         !& prDiag, managerObject)

        if (managerObject%asn1 == 1.0d0) then

          do i = 1, managerObject%num1stOEqnsSolve
            do kloop = 1, ktloop
              smvdm(kloop,i) =  &
     &          smvdm(kloop,i) + managerObject%dtlos(kloop,i) + explic(kloop,i)
              cnewDerivatives(kloop,i) = cnewDerivatives(kloop,i) + managerObject%dtlos(kloop,i)
            end do
          end do

        else

          do i = 1, managerObject%num1stOEqnsSolve
            do kloop = 1, ktloop
              dtasn1         = managerObject%asn1 * managerObject%dtlos(kloop,i)
              smvdm(kloop,i) = smvdm(kloop,i) + dtasn1 + explic(kloop,i)
              cnewDerivatives (kloop,i) = cnewDerivatives (kloop,i) + dtasn1
            end do
          end do

        end if

!       ---------------------------------------------------
!       Exit smvgear if a time interval has been completed.
!       ---------------------------------------------------

        managerObject%timeremain = managerObject%chemTimeInterval - managerObject%xelaps
        !   print*, "checking time interval"

        if (managerObject%timeremain <= 1.0d-06) return

!       -------------------------------------------------------------------
!       idoub counts the number of successful steps before re-testing the
!       step-size and order:
!         if idoub > 1, decrease idoub and go on to the next time step with
!                       the current step-size and order;
!         if idoub = 1, store the value of the error (dtlos) for the time
!                       step prediction, which will occur when idoub = 0,
!                       but go on to the next step with the current step
!                       size and order;
!         if idoub = 0, test the time step and order for a change.
!       -------------------------------------------------------------------

        if (managerObject%idoub > 1) then

          managerObject%idoub = managerObject%idoub - 1

          if (managerObject%idoub == 1) then

            do jspc = 1, managerObject%num1stOEqnsSolve, 2

              jg1 = jspc + 1

              do kloop = 1, ktloop
                cest(kloop,jspc) = managerObject%dtlos(kloop,jspc)
                cest(kloop,jg1)  = managerObject%dtlos(kloop,jg1)
              end do

            end do

          end if

          managerObject%rdelt = 1.0d0

!         =========
          go to 200
!         =========

        end if

!     ================
      end if DER2MAXIF
!     ================


!     ------------------------------------------------------------------
!     Test whether to change the step-size and order.
!
!     Determine the time step at (a) one order lower than, (b) the same
!     order as, and (c) one order higher than the current order.  In the
!     case of multiple grid-cells in a grid-block, find the minimum
!     step size among all the cells for each of the orders.  Then, in
!     all cases, choose the longest time step among the three steps
!     paired with orders, and choose the order allowing this longest
!     step.
!     ------------------------------------------------------------------

!     ---------------------------------------------------------------
!     Estimate the time step ratio (rdeltup) at one order higher than
!     the current order.  If nqq >= MAXORD, then we do not allow the
!     order to increase.
!     ---------------------------------------------------------------

      if (managerObject%nqq < MAXORD) then

        do kloop = 1, ktloop
          dely(kloop) = 0.0d0
        end do

        do jspc = 1, managerObject%num1stOEqnsSolve
          do kloop = 1, ktloop
            errymax     = (managerObject%dtlos(kloop,jspc) - cest(kloop,jspc)) *  &
     &                    managerObject%chold(kloop,jspc)
            dely(kloop) = dely(kloop) + (errymax * errymax)
          end do
        end do

        der3max = 0.0d0

        do kloop = 1, ktloop

          if (dely(kloop) > der3max) then
            der3max = dely(kloop)
          end if

        end do

        managerObject%rdeltup = 1.0d0 / ((managerObject%conp3 * der3max**savedVars%enqq3(managerObject%nqq)) + 1.4d-6)

      else

        managerObject%rdeltup = 0.0d0

      end if


!     ========
 400  continue
!     ========

      call estimateTimeStepRatio (managerObject, ktloop, dely, cnewDerivatives, savedVars)


!     ---------------------------------------------------------------
!     If the last step was successful and rdelt is small, keep the
!     current step and order, and allow three successful steps before
!     re-checking the time step and order.
!     ---------------------------------------------------------------

      if ((managerObject%rdelt < 1.1d0) .and. (managerObject%ifsuccess == 1)) then

        managerObject%idoub = 3

!       =========
        go to 200
!       =========

!       --------------------------------------------------------------
!       If the maximum time step ratio is that of one order lower than
!       the current order, decrease the order.  Do not minimize rdelt
!       to <= 1, when ifsuccess = 0 since this is less efficient.
!       --------------------------------------------------------------

      else if (managerObject%rdelt == managerObject%rdeltdn) then


        managerObject%nqq = managerObject%nqq - 1

!       ---------------------------------------------------------------
!       If the maximum time step ratio is that of one order higher than
!       the current order, increase the order and add a derivative term
!       for the higher order.
!       ---------------------------------------------------------------


      else if (managerObject%rdelt == managerObject%rdeltup) then

        real_kstep = managerObject%kstep
        consmult   = savedVars%aset(managerObject%nqq,managerObject%kstep) / real_kstep
        managerObject%nqq        = managerObject%kstep
        nqisc      = managerObject%nqq * managerObject%num1stOEqnsSolve

        do jspc = 1, managerObject%num1stOEqnsSolve, 2

          jg1 = jspc + 1
          i1  = jspc + nqisc
          i2  = jg1  + nqisc

          do kloop = 1, ktloop
            cnewDerivatives(kloop,i1) = managerObject%dtlos(kloop,jspc) * consmult
            cnewDerivatives(kloop,i2) = managerObject%dtlos(kloop,jg1)  * consmult
          end do

        end do

      end if

!     ----------------------------------------------------------------
!     If the last two steps have failed, re-set idoub to the current
!     order + 1.  Do not minimize rdelt if jfail >= 2 since tests show
!     that this merely leads to additional computations.
!     ----------------------------------------------------------------

      managerObject%idoub = managerObject%nqq + 1

!     =========
      go to 200
!     =========

      return

      end subroutine Smvgear

      subroutine initConcentrationArray(ktloop, concentrationsNew, concentrationsOld, managerObject)

         use GmiManager_mod
         implicit none

#     include "smv2chem_par.h"

         integer, intent(in) :: ktloop
         real*8, intent(out) :: concentrationsNew(KBLOOP, MXGSAER)
         real*8,  intent(in)  :: concentrationsOld(KBLOOP, MXGSAER)
         type (Manager_type) :: managerObject

         integer :: jnew, kloop

         do jnew = 1, managerObject%num1stOEqnsSolve
           do kloop = 1, ktloop
             concentrationsNew(kloop, jnew) = concentrationsOld(kloop, jnew) ! why save this?
           end do
         end do

      end subroutine


      !   smvdm    : amount added to each spc at each grid-cell (# cm^-3 for gas chemistry (?))
      subroutine updateChemistryMassBalance (ktloop, cnewDerivatives, explic, smvdm, &
         & managerObject)

         use GmiManager_mod
         implicit none

#     include "smv2chem_par.h"

         integer, intent(in) :: ktloop
         real*8, intent(inout) :: cnewDerivatives(KBLOOP, MXGSAER*7)
         real*8, intent(in) :: explic(KBLOOP, MXGSAER)
         real*8, intent(inout) :: smvdm(KBLOOP, MXGSAER)
         type (Manager_type) :: managerObject

         real :: dtasn1
         integer :: i,kloop

         ! double check, but it looks like we can eliminate the first part

       if (managerObject%asn1 == 1.0d0) then
          do i = 1, managerObject%num1stOEqnsSolve
            do kloop = 1, ktloop
              smvdm(kloop,i) = smvdm(kloop,i) + managerObject%dtlos(kloop,i) + explic(kloop,i)
              cnewDerivatives(kloop,i) = cnewDerivatives(kloop,i) + managerObject%dtlos(kloop,i)
            end do
          end do

        else

          do i = 1, managerObject%num1stOEqnsSolve
            do kloop = 1, ktloop
              dtasn1         = managerObject%asn1 * managerObject%dtlos(kloop,i)
              smvdm(kloop,i) = smvdm(kloop,i) + dtasn1 + explic(kloop,i)
              cnewDerivatives (kloop,i) = cnewDerivatives (kloop,i) + dtasn1
            end do
          end do

        end if

      end subroutine



!-----------------------------------------------------------------------------
!
! ROUTINE
!   Backsub
!
! DESCRIPTION
!   This routine performs back-substitutions on the decomposed matrix.  It
!   solves the linear set of equations Ax = B FOR x, the correction vector,
!   where "A" is the L-U decompostion of the original matrix =>
!
!     P = I - H x Bo x J
!
!   I = identity matrix, H = time step, Bo = a coefficient that depends on
!   the order of the integration method, and J is the matrix of partial
!   derivatives.  B is sent from Smvgear as a corrected value of the first
!   derivatives of the ordinary differential equations.  Decomp solved for
!   "A", the decomposed matrix.  See Press, et. al. (1992), Numerical
!   Recipes, Cambridge University Press, for a better description of the
!   back-substitution process.
!
!   This back-substitution process uses sparse matrix techniques,
!   vectorizes around the grid-cell dimension, and uses no partial
!   pivoting.  Tests by Sherman & Hindmarsh (1980), Lawrence Livermore
!   Livermore Laboratory, Rep. UCRL-84102, and by us have confirmed that
!   the removal of partial pivoting has little effect on results.
!
!   Backsub loop # 1 =>
!     First, adjust right side of Ax = B using lower triangular matrix.
!     Sum 1,2,3,4, or 5 terms at a time to improve vectorization.
!
! ARGUMENTS
!   num1stOEqnsSolve : # of first-order eqns to solve, = # of spc = order of original
!            matrix; num1stOEqnsSolve has a different value for day and night, and for
!            gas- and aqueous-phase chemistry;
!            # spc with prod or loss terms in Smvgear (?)
!   ktloop : # of grid-cells in a grid-block
!   ncsp   : ncs       => for daytime   gas chemistry
!            ncs + ICS => for nighttime gas chemistry
!   cc2    : array holding values of decomposed matrix
!   vdiag  : 1 / current diagonal term of the decomposed matrix
!   gloss  : first derivative = sum of prod minus loss rates for a spc
!
!-----------------------------------------------------------------------------

      subroutine Backsub  &
     &  (savedVars, num1stOEqnsSolve, ktloop, ncsp, cc2, vdiag, gloss)

      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved
      implicit none

#     include "smv2chem_par.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: num1stOEqnsSolve
      integer, intent(in)  :: ktloop
      integer, intent(in)  :: ncsp
      real*8,  intent(in)  :: cc2  (KBLOOP, 0:MXARRAY)
      real*8,  intent(in)  :: vdiag(KBLOOP, MXGSAER)

      real*8,  intent(inout) :: gloss(KBLOOP, MXGSAER)

      type(t_Smv2Saved), intent(inOut) :: savedVars

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: i, ij
      integer :: ij0, ij1, ij2, ij3, ij4

      integer :: j0, j1, j2, j3, j4

      integer :: k, kc, kzt
      integer :: kh1, kh2, kh3, kh4, kh5
      integer :: kl1, kl2, kl3, kl4, kl5

      integer :: mc, mzt
      integer :: mh1, mh2, mh3, mh4, mh5
      integer :: ml1, ml2, ml3, ml4, ml5


!     ----------------
!     Begin execution.
!     ----------------

!c    Write (6,*) 'Backsub called.'


      ij = 1


!     ==========================================
      KZTLOOP: do kzt = savedVars%kztlo(ncsp), savedVars%kzthi(ncsp)
!     ==========================================

        i = savedVars%ikztot(kzt)

        kl5 = savedVars%kbl5(kzt)
        kh5 = savedVars%kbh5(kzt)
        kl4 = savedVars%kbl4(kzt)
        kh4 = savedVars%kbh4(kzt)
        kl3 = savedVars%kbl3(kzt)
        kh3 = savedVars%kbh3(kzt)
        kl2 = savedVars%kbl2(kzt)
        kh2 = savedVars%kbh2(kzt)
        kl1 = savedVars%kbl1(kzt)
        kh1 = savedVars%kbh1(kzt)

!       -- Sum 5 terms at a time. --

        do kc = kl5, kh5

          ij0 = ij
          ij1 = ij + 1
          ij2 = ij + 2
          ij3 = ij + 3
          ij4 = ij + 4
          ij  = ij + 5

          j0  = savedVars%kzeroa(kc)
          j1  = savedVars%kzerob(kc)
          j2  = savedVars%kzeroc(kc)
          j3  = savedVars%kzerod(kc)
          j4  = savedVars%kzeroe(kc)

          !K: Hot loop
          do k = 1, ktloop
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (cc2(k,ij0) * gloss(k,j0)) -  &
     &        (cc2(k,ij1) * gloss(k,j1)) -  &
     &        (cc2(k,ij2) * gloss(k,j2)) -  &
     &        (cc2(k,ij3) * gloss(k,j3)) -  &
     &        (cc2(k,ij4) * gloss(k,j4))
          end do

        end do

!       -- Sum 4 terms at a time. --

        do kc = kl4, kh4

          ij0 = ij
          ij1 = ij + 1
          ij2 = ij + 2
          ij3 = ij + 3
          ij  = ij + 4

          j0  = savedVars%kzeroa(kc)
          j1  = savedVars%kzerob(kc)
          j2  = savedVars%kzeroc(kc)
          j3  = savedVars%kzerod(kc)

          do k = 1, ktloop
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (cc2(k,ij0) * gloss(k,j0)) -  &
     &        (cc2(k,ij1) * gloss(k,j1)) -  &
     &        (cc2(k,ij2) * gloss(k,j2)) -  &
     &        (cc2(k,ij3) * gloss(k,j3))
          end do

        end do

!       -- Sum 3 terms at a time. --

        do kc = kl3, kh3

          ij0 = ij
          ij1 = ij + 1
          ij2 = ij + 2
          ij  = ij + 3

          j0  = savedVars%kzeroa(kc)
          j1  = savedVars%kzerob(kc)
          j2  = savedVars%kzeroc(kc)

          do k = 1, ktloop
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (cc2(k,ij0) * gloss(k,j0)) -  &
     &        (cc2(k,ij1) * gloss(k,j1)) -  &
     &        (cc2(k,ij2) * gloss(k,j2))
          end do

        end do

!       -- Sum 2 terms at a time. --

        do kc = kl2, kh2

          ij0 = ij
          ij1 = ij + 1
          ij  = ij + 2

          j0  = savedVars%kzeroa(kc)
          j1  = savedVars%kzerob(kc)

          do k = 1, ktloop
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (cc2(k,ij0) * gloss(k,j0)) -  &
     &        (cc2(k,ij1) * gloss(k,j1))
          end do

        end do

!       -- Sum 1 term at a time. --

        do kc = kl1, kh1

          ij0 = ij
          ij  = ij + 1

          j0  = savedVars%kzeroa(kc)

          do k = 1, ktloop
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (cc2(k,ij0) * gloss(k,j0))
          end do

        end do

!     ==============
      end do KZTLOOP
!     ==============


!     ---------------------------------------------------------------
!     Backsub loop # 2.
!
!     Backsubstite with upper triangular matrix to find solution.
!     Again, sum up several terms at a time to improve vectorization.
!     ---------------------------------------------------------------

!     ===========================
      ILOOP: do i = num1stOEqnsSolve, 1, -1
!     ===========================

        mzt = savedVars%imztot(i,ncsp)

!       ===================
        MZTIF: if (mzt > 0) then
!       ===================

          ml5 = savedVars%mbl5(mzt)
          mh5 = savedVars%mbh5(mzt)
          ml4 = savedVars%mbl4(mzt)
          mh4 = savedVars%mbh4(mzt)
          ml3 = savedVars%mbl3(mzt)
          mh3 = savedVars%mbh3(mzt)
          ml2 = savedVars%mbl2(mzt)
          mh2 = savedVars%mbh2(mzt)
          ml1 = savedVars%mbl1(mzt)
          mh1 = savedVars%mbh1(mzt)

!         -- Sum 5 terms at a time. --

          do mc = ml5, mh5

            ij0 = ij
            ij1 = ij + 1
            ij2 = ij + 2
            ij3 = ij + 3
            ij4 = ij + 4
            ij  = ij + 5

            j0  = savedVars%mzeroa(mc)
            j1  = savedVars%mzerob(mc)
            j2  = savedVars%mzeroc(mc)
            j3  = savedVars%mzerod(mc)
            j4  = savedVars%mzeroe(mc)

            do k = 1, ktloop
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (cc2(k,ij0) * gloss(k,j0)) -  &
     &          (cc2(k,ij1) * gloss(k,j1)) -  &
     &          (cc2(k,ij2) * gloss(k,j2)) -  &
     &          (cc2(k,ij3) * gloss(k,j3)) -  &
     &          (cc2(k,ij4) * gloss(k,j4))
            end do

          end do

!         -- Sum 4 terms at a time. --

          do mc = ml4, mh4

            ij0 = ij
            ij1 = ij + 1
            ij2 = ij + 2
            ij3 = ij + 3
            ij  = ij + 4

            j0  = savedVars%mzeroa(mc)
            j1  = savedVars%mzerob(mc)
            j2  = savedVars%mzeroc(mc)
            j3  = savedVars%mzerod(mc)

            do k = 1, ktloop
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (cc2(k,ij0) * gloss(k,j0)) -  &
     &          (cc2(k,ij1) * gloss(k,j1)) -  &
     &          (cc2(k,ij2) * gloss(k,j2)) -  &
     &          (cc2(k,ij3) * gloss(k,j3))
            end do

          end do

!         -- Sum 3 terms at a time. --

          do mc = ml3, mh3

            ij0 = ij
            ij1 = ij + 1
            ij2 = ij + 2
            ij  = ij + 3

            j0  = savedVars%mzeroa(mc)
            j1  = savedVars%mzerob(mc)
            j2  = savedVars%mzeroc(mc)

            do k = 1, ktloop
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (cc2(k,ij0) * gloss(k,j0)) -  &
     &          (cc2(k,ij1) * gloss(k,j1)) -  &
     &          (cc2(k,ij2) * gloss(k,j2))
            end do

          end do

!         -- Sum 2 terms at a time. --

          do mc = ml2, mh2

            ij0 = ij
            ij1 = ij + 1
            ij  = ij + 2

            j0  = savedVars%mzeroa(mc)
            j1  = savedVars%mzerob(mc)

            do k = 1, ktloop
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (cc2(k,ij0) * gloss(k,j0)) -  &
     &          (cc2(k,ij1) * gloss(k,j1))
            end do

          end do

!         -- Sum 1 term at a time. --

          do mc = ml1, mh1

            ij0 = ij
            ij  = ij + 1

            j0  = savedVars%mzeroa(mc)

            do k = 1, ktloop
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (cc2(k,ij0) * gloss(k,j0))
            end do

          end do

!       ============
        end if MZTIF
!       ============

!       -- Adjust gloss with diagonal element. --

        do k = 1, ktloop
          gloss(k,i) = gloss(k,i) * vdiag(k,i)
        end do

!     ============
      end do ILOOP
!     ============


      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Decomp
!
! DESCRIPTION
!   This routine decomposes the sparse matrix "P" into the matrix "A" in
!   order to solve the linear set of equations Ax = B for x, which is a
!   correction vector.  Ax = B is solved in Backsub, the original matrix
!   "P" is =>
!
!     P = I - H x Bo x J
!
!   where I = identity matrix, H = time step, Bo = a coefficient that
!   depends on the order of the integration method, and J is the matrix of
!   partial derivatives.  See Press, et. al. (1992), Numerical Recipes,
!   Cambridge University Press, for a better description of the L-U
!   decompostion process.
!
!   This L-U decompostion process uses sparse matrix techniques, vectorizes
!   around the grid-cell dimension, and uses no partial pivoting.  Tests by
!   Sherman & Hindmarsh (1980), Lawrence Livermore National Laboratory,
!   Rep. UCRL-84102, and by us have confirmed that the removal of partial
!   pivoting has little effect on results.
!
! ARGUMENTS
!   num1stOEqnsSolve : # of first-order eqns to solve, = # of spc = order of original
!            matrix; num1stOEqnsSolve has a different value for day and night, and for
!            gas- and aqueous-phase chemistry;
!            # spc with prod or loss terms in Smvgear (?)
!   ktloop : # of grid-cells in a grid-block
!   ncsp   : ncs       => for daytime   gas chemistry
!            ncs + ICS => for nighttime gas chemistry
!   cc2    : array of iarray units holding values of each matrix
!            position actually used; originally,
!            cc2 = P = I - delt * aset(nqq,1) * partial_derivatives;
!            however, cc2 is decomposed here
!   vdiag  : 1 / current diagonal term of the decomposed matrix
!
!-----------------------------------------------------------------------------
! MRD: LU Decomp  - should go into sparseMatrix module
      subroutine Decomp  &
     &  (savedVars, num1stOEqnsSolve, ktloop, ncsp, cc2, vdiag)

      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved
      implicit none

#     include "smv2chem_par.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: num1stOEqnsSolve
      integer, intent(in)  :: ktloop
      integer, intent(in)  :: ncsp

      real*8,  intent(inout) :: cc2  (KBLOOP, 0:MXARRAY)
      real*8,  intent(inout) :: vdiag(KBLOOP, MXGSAER)

      type(t_Smv2Saved), intent(inOut) :: savedVars

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: iar, ic
      integer :: ih1, ih2, ih3, ih4, ih5
      integer :: ij, ija, ijt
      integer :: ik0, ik1, ik2, ik3, ik4
      integer :: il1, il2, il3, il4, il5
      integer :: j, jc, jh, jl
      integer :: k
      integer :: kj0, kj1, kj2, kj3, kj4


!     ----------------
!     Begin execution.
!     ----------------

!c    Write (6,*) 'Decomp called.'


!     -----------------------------------------------------------
!     First loop of L-U decompostion.
!
!     Sum 1,2,3,4, OR 5 terms at a time to improve vectorization.
!     -----------------------------------------------------------
!     =======================
      JLOOP: do j = 1, num1stOEqnsSolve
!     =======================

!       ==============================================
        IJTLOOP: do ijt =  savedVars%ijtlo(j,ncsp),  savedVars%ijthi(j,ncsp)
!       ==============================================

         !MRD: all things with 5 terms
         ! should be part of sparse matrix type
          ij  =  savedVars%ijval(ijt)
          il5 =  savedVars%idl5 (ijt)
          ih5 =  savedVars%idh5 (ijt)
          il4 =  savedVars%idl4 (ijt)
          ih4 =  savedVars%idh4 (ijt)
          il3 =  savedVars%idl3 (ijt)
          ih3 =  savedVars%idh3 (ijt)
          il2 =  savedVars%idl2 (ijt)
          ih2 =  savedVars%idh2 (ijt)
          il1 =  savedVars%idl1 (ijt)
          ih1 =  savedVars%idh1 (ijt)

!         -- Sum 5 terms at a time. --
          ! MRD: does this unrolling really help...?
          !K: Unclear, I'd rather remove it for clarity if nothing else
          !K: But this involves unifying ikdeca,etc

          ! should be part of sparse matrix type
          do ic = il5, ih5

            ik0 =  savedVars%ikdeca(ic)
            ik1 =  savedVars%ikdecb(ic)
            ik2 =  savedVars%ikdecc(ic)
            ik3 =  savedVars%ikdecd(ic)
            ik4 =  savedVars%ikdece(ic)

            kj0 =  savedVars%kjdeca(ic)
            kj1 =  savedVars%kjdecb(ic)
            kj2 =  savedVars%kjdecc(ic)
            kj3 =  savedVars%kjdecd(ic)
            kj4 =  savedVars%kjdece(ic)

            do k = 1, ktloop
              cc2(k,ij) =  &
     &          cc2(k,ij) -  &
     &          (cc2(k,ik0) * cc2(k,kj0)) -  &
     &          (cc2(k,ik1) * cc2(k,kj1)) -  &
     &          (cc2(k,ik2) * cc2(k,kj2)) -  &
     &          (cc2(k,ik3) * cc2(k,kj3)) -  &
     &          (cc2(k,ik4) * cc2(k,kj4))
            end do

          end do

!         -- Sum 4 terms at a time. --

          do ic = il4, ih4

            ik0 =  savedVars%ikdeca(ic)
            ik1 =  savedVars%ikdecb(ic)
            ik2 =  savedVars%ikdecc(ic)
            ik3 =  savedVars%ikdecd(ic)

            kj0 =  savedVars%kjdeca(ic)
            kj1 =  savedVars%kjdecb(ic)
            kj2 =  savedVars%kjdecc(ic)
            kj3 =  savedVars%kjdecd(ic)

            do k = 1, ktloop
              cc2(k,ij) =  &
     &          cc2(k,ij) -  &
     &          (cc2(k,ik0) * cc2(k,kj0)) -  &
     &          (cc2(k,ik1) * cc2(k,kj1)) -  &
     &          (cc2(k,ik2) * cc2(k,kj2)) -  &
     &          (cc2(k,ik3) * cc2(k,kj3))
            end do

          end do

!         -- Sum 3 terms at a time. --

          do ic = il3, ih3

            ik0 =  savedVars%ikdeca(ic)
            ik1 =  savedVars%ikdecb(ic)
            ik2 =  savedVars%ikdecc(ic)

            kj0 =  savedVars%kjdeca(ic)
            kj1 =  savedVars%kjdecb(ic)
            kj2 =  savedVars%kjdecc(ic)

            do k = 1, ktloop
              cc2(k,ij) =  &
     &          cc2(k,ij) -  &
     &          (cc2(k,ik0) * cc2(k,kj0)) -  &
     &          (cc2(k,ik1) * cc2(k,kj1)) -  &
     &          (cc2(k,ik2) * cc2(k,kj2))
            end do

          end do

!         -- Sum 2 terms at a time. --

          do ic = il2, ih2

            ik0 =  savedVars%ikdeca(ic)
            ik1 =  savedVars%ikdecb(ic)

            kj0 =  savedVars%kjdeca(ic)
            kj1 =  savedVars%kjdecb(ic)

            do k = 1, ktloop
              cc2(k,ij) =  &
     &          cc2(k,ij) -  &
     &          (cc2(k,ik0) * cc2(k,kj0)) -  &
     &          (cc2(k,ik1) * cc2(k,kj1))
            end do

          end do

!         -- Sum 1 term  at a time. --

          do ic = il1, ih1

            ik0 =  savedVars%ikdeca(ic)

            kj0 =  savedVars%kjdeca(ic)

            do k = 1, ktloop
              cc2(k,ij) =  &
     &          cc2(k,ij) -  &
     &          (cc2(k,ik0) * cc2(k,kj0))
            end do

          end do

!       ==============
        end do IJTLOOP
!       ==============

        iar =  savedVars%jarrdiag(j,ncsp)

        do k = 1, ktloop
          vdiag(k,j) = 1.0d0 / cc2(k,iar)
        end do

!       ----------------------------
!       Second loop of decompostion.
!       ----------------------------

        jl =  savedVars%jloz1(j,ncsp)
        jh =  savedVars%jhiz1(j,ncsp)

        do jc = jl, jh

          ija =  savedVars%jzeroa(jc)

          do k = 1, ktloop
            cc2(k,ija) = cc2(k,ija) * vdiag(k,j)
          end do

        end do

!     ============
      end do JLOOP
!     ============


      return

      end subroutine Decomp






