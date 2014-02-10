
!=============================================================================
!
! $Id: physproc.F90,v 1.9 2013-09-03 16:16:43 jkouatch Exp $
!
! CODE DEVELOPER
!   Original code from Mark Z. Jacobson ((C) COPYRIGHT, 1993).
!   LLNL modifications:  John Tannahill
!                        jrt@llnl.gov
!
! FILE
!   physproc.F
!
! ROUTINES
!   Physproc
!   Deter_Block_Size
!   Solve_Block
!   Calcrate
!   Reorder_Grid_Cells
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   physProc
!
! DESCRIPTION
!   This routine solves gas-phase chemical equations.  It divides the grid-
!   domain into grid-blocks, and the code vectorizes around the number of
!   grid-cells in each block;
!     grid-block  => a group of grid-cells in the grid-domain,
!     grid-cell   => a single grid-box,
!     grid-domain => all grid-cells in the model.
!   It sets compact arrays and calls Smvgear (Sparse Matrix Vectorized Gear
!   code).
!
! ARGUMENTS
!   doQqjkInChem   : if prQqjk is on, should qqj's & qqk's be determined
!                      inside the chemistry solver, or outside?
!   doSurfEmissInChem : do surface emissions inside the chemistry solver, or
!                      outside?
!   prQqjk   : should the periodic qqjk output file be written?
!   prSmv2   : should the SmvgearII     output file be written
!               (non-parallel mode only)?
!   numLats      : # of latitudes
!   numLongs     : # of longitudes
!   numVert     : # of vertical layers
!   doReOrder   : if 1, then reorder grid-cells by stiffness
!   airDensityIndex     : array index for air density
!   nitrogenSpNum : identifies spc # of nitrogen gas
!   oxygenSpNum   : identifies spc # of oxygen   gas
!   numZones    : # of zones (numLongs * numLats * numVert)
!   numGridCellsInBlock    : intended # of grid-cells in a grid-block
!   smvUnitNumber    : logical unit number to write to when prSmv2 is true
!   gasChemistryType       : identifies gas chemistry type (1.NCSGAS)
!   fractionDecrease   : fraction time step is decreased in Smvgear if convergence
!               test fails
!   maxTimeStepNight   : max time step for night all chem (s)
!   ncOutPeriod     : NetCDF output period
!   modelTimeStep       : model time step (s)
!   doCellChem     : do chemistry for a particular cell?
!   jPhotRate  : tbd
!   numKineticRxns    : # of kinetic rxns (non-photo)
!   ntLoopNcs : tbd
!   numActInActGases    : # of active + inactive gases
!   origSpcNumber   : original spc # of each new jnew spc
!   npphotrat : tbd
!   arate     : thermal    rate constants (units vary)
!   prate     : photolysis rate constants (s^-1)
!   yemis     : surface emissions (molecules/cm^3/sec)
!   jreorder  : gives original grid-cell from re-ordered grid-cell
!   lreorder  : gives original grid-cell from re-ordered cell, except when
!               cell is a virtual boundary cell, then it gives original
!               edge cell from re-ordered virtual boundary cell
!   csuma     : tbd
!   csumc     : tbd
!   errmx2    : tbd
!   cx        : spc conc (molec/cm^3)
!
! HISTORY
!   - July 1, 2004 - Jules Kouatchou
!       o Use the pre-processing option MSG_OPTION to remove MPI
!         statements when MSG_OPTION is set to MSG_NONE.
!-----------------------------------------------------------------------------

      subroutine Physproc  &
     &  (savedVars, doQqjkInChem, doSurfEmissInChem, prQqjk, prSmv2, numLats,  &
     &   numLongs, numVert, doReOrder, airDensityIndex, nitrogenSpNum, oxygenSpNum, numZones,  &
     &   numGridCellsInBlock, smvUnitNumber, gasChemistryType, fractionDecrease, maxTimeStepNight, ncOutPeriod, modelTimeStep,  &
     &   doCellChem, jPhotRate, numKineticRxns, ntLoopNcs, numActInActGases, origSpcNumber,  &
     &   npphotrat, arate, prate, yemis, jreorder, lreorder, csuma,  &
     &   csumc, errmx2, cx, yda, qqkda, qqjda, &
     &   i1, i2, ju1, j2, k1, k2, &
     &   num_qks, num_qjs, num_active, commuWorld)

      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

      implicit none

#     include "smv2chem_par.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: commuWorld
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: num_qks, num_qjs, num_active
      real*8 , intent(inout) :: qqjda(i1:i2, ju1:j2, k1:k2, num_qjs)
      real*8 , intent(inout) :: qqkda(i1:i2, ju1:j2, k1:k2, num_qks)
      real*8 , intent(inout) :: yda  (i1:i2, ju1:j2, k1:k2, num_active)
      logical, intent(in)  :: doQqjkInChem
      logical, intent(in)  :: doSurfEmissInChem
      logical, intent(in)  :: prQqjk
      logical, intent(in)  :: prSmv2
      integer, intent(in)  :: numLats, numLongs, numVert
      integer, intent(in)  :: doReOrder
      integer, intent(in)  :: airDensityIndex
      integer, intent(in)  :: nitrogenSpNum
      integer, intent(in)  :: oxygenSpNum
      integer, intent(in)  :: numZones
      integer, intent(in)  :: numGridCellsInBlock
      integer, intent(in)  :: smvUnitNumber
      integer, intent(in)  :: gasChemistryType
      real*8,  intent(in)  :: fractionDecrease
      real*8,  intent(in)  :: maxTimeStepNight
      real*8,  intent(in)  :: ncOutPeriod
      real*8,  intent(in)  :: modelTimeStep
      logical, intent(in)  :: doCellChem(numZones)
      integer, intent(in)  :: jPhotRate (ICS)
      integer, intent(in)  :: numKineticRxns   (ICS)
      integer, intent(in)  :: ntLoopNcs(ICS)
      integer, intent(in)  :: numActInActGases   (ICS)
      integer, intent(in)  :: origSpcNumber  (MXGSAER, ICS)
      integer, intent(in)  :: npphotrat(IPHOT,   ICS)
      real*8,  intent(in)  :: arate    (numZones,  ITHERM)
      real*8,  intent(in)  :: prate    (numZones,  IPHOT)
      real*8,  intent(in)  :: yemis    (numLats*numLongs, IGAS)

      integer, intent(inout) :: jreorder(numZones)
      integer, intent(inout) :: lreorder(numZones)
      real*8,  intent(inout) :: csuma   (numZones)
      real*8,  intent(inout) :: csumc   (numZones)
      real*8,  intent(inout) :: errmx2  (numZones)
      real*8,  intent(inout) :: cx      (numZones, IGAS)

      type(t_Smv2Saved), intent(inOut) :: savedVars

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: idaynt
      integer :: ifSun      ! 1 => for daytime, 2 => for nighttime
      integer :: ireord     ! 1 => reorder grid-cells and blocks for chemistry
                            ! 2 => solve chemistry
      integer :: iday
      integer :: jloop, jloopn
      integer :: loreord    ! 1 => if reordering, 2 => if no reordering
      integer :: nblockuse  ! (ireord == loreord) => # of original  blocks
                            ! (ireord /= loreord) => # of reordered blocks
      integer :: nreblock
      integer :: ntloopuse

      integer :: jlowvar(MXBLOCK)
      integer :: ktlpvar(MXBLOCK)

      real*8, parameter :: PHOTOLYSIS_THRESHOLD = 1.0d-80
      logical :: isDaytimeInFirstCell, isDaytime, dayAndNight
      integer ::  isDaytimeInFirstCell_I, isDaytime_I, dayAndNight_I

!     ----------------
!     Begin execution.
!     ----------------

#ifndef nonZeroInd_tracers
!     --------------------------------------------------------
!     Can the domain be separated into day and night sections?
!     --------------------------------------------------------

      if (prate(1,1) >= PHOTOLYSIS_THRESHOLD) then
        ifSun = 1
      else
        ifSun = 2
      end if
      !print*, "First cell check, ifSun = ", ifSun

      ! The can replace code block above
      isDaytimeInFirstCell = (prate(1,1) >= PHOTOLYSIS_THRESHOLD)
      ! TODO MRD: isDayTimeInFirstCell = T should be equivalent to 1 (see ifSun)

!     --------------------------------------------------------
!     If first cell is day time, and another cell is not, then
!     set idaynt to indicate day and night (idaynt=2)
!     --------------------------------------------------------
      idaynt = 1

      if (ifSun == 1) then

        ifSun1: do jloop = 2, ntLoopNcs(gasChemistryType)
          if (prate(jloop,1) .lt. PHOTOLYSIS_THRESHOLD) then
            idaynt = 2
            !print*, "ifSun = 1 and idaynt flipped"
            exit ifSun1
          end if
        end do ifSun1

!     --------------------------------------------------------
!     If first cell is night time, and another cell is not, then
!     set idaynt to indicate day and night (idaynt=2)
!     --------------------------------------------------------

      else if (ifSun == 2) then

        ifSun2: do jloop = 2, ntLoopNcs(gasChemistryType)
          if (prate(jloop,1) .gt. PHOTOLYSIS_THRESHOLD) then
            idaynt = 2
            !print*, "ifSun = 2 and idaynt flipped"
            exit ifSun2
          end if
        end do ifSun2

      end if

      ! The can replace the two code blocks above
      dayAndNight = .false. ! unless ...
      do jloop = 2, ntLoopNcs(gasChemistryType)
        isDaytime = (prate(jloop,1) >= PHOTOLYSIS_THRESHOLD)
        if (isDaytime .neqv. isDaytimeInFirstCell) then
            dayAndNight = .true.
            exit
        end if
      end do
      ! TODO MRD: dayAndNight = .true. should be equivalent to 2 (see idaynt)


!     --------------------------------------------------
!     Reorder cells and blocks then solve chemical odes.
!     --------------------------------------------------
      loreord = 2
      if ((doReOrder == 1) .and. (numZones > 1)) then
        loreord = 1
      end if

      do iday = 1, idaynt
        do ireord = loreord, 2

          if (ireord == loreord) then
!           =====================
            call Deter_Block_Size  &
!           =====================
     &        (savedVars, iday, idaynt, numZones, numGridCellsInBlock, gasChemistryType, doCellChem,  &
     &         ntLoopNcs, prate, ifSun, nblockuse, ntloopuse, jlowvar,  &
     &         ktlpvar, jreorder)
          else
            nblockuse = nreblock
          end if


          do jloopn = 1, ntloopuse
            lreorder(jloopn) = jreorder(jloopn)
          end do


!         ================
          call Solve_Block  &
!         ================
     &      (savedVars, doQqjkInChem, doSurfEmissInChem, prQqjk, prSmv2, ifSun,  &
     &       numLats, numLongs, numVert, airDensityIndex, nitrogenSpNum, oxygenSpNum, ireord,  &
     &       numZones, smvUnitNumber, nblockuse, gasChemistryType, fractionDecrease, maxTimeStepNight,  &
     &       ncOutPeriod, modelTimeStep, doCellChem, jPhotRate, numKineticRxns, numActInActGases,  &
     &       jlowvar, ktlpvar, origSpcNumber, npphotrat, arate, prate, yemis,  &
     &       jreorder, lreorder, errmx2, cx, yda, qqkda, qqjda, &
     &       i1, i2, ju1, j2, k1, k2, num_qks, num_qjs, num_active, &
     &       commuWorld)


          if (ireord == 1) then
!           =======================
            call Reorder_Grid_Cells  &
!           =======================
     &        (numZones, numGridCellsInBlock, ntloopuse, errmx2, jreorder, csuma,  &
     &         nreblock, lreorder, jlowvar, ktlpvar, csumc)
          end if

        end do
      end do

#endif

      return

      end


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Deter_Block_Size
!
! DESCRIPTION
!   This routine determines block sizes for chemistry or stiffness.
!
! ARGUMENTS
!   iday      : tbd
!   idaynt    : tbd
!   numZones    : # of zones (numLongs * numLats * numVert)
!   numGridCellsInBlock    : intended # of grid-cells in a grid-block
!   gasChemistryType       : identifies gas chemistry type (1..NCSGAS)
!   doCellChem : do chemistry for a particular cell?
!   ntLoopNcs : tbd
!   prate     : photolysis rate constants (s^-1)
!   ifSun     : 1 => for daytime, 2 => for nighttime
!   nblockuse : # of original blocks
!   ntloopuse : tbd
!   jlowvar   : tbd
!   ktlpvar   : tbd
!   jreorder  : gives original grid-cell from re-ordered grid-cell
!
!-----------------------------------------------------------------------------

      subroutine Deter_Block_Size  &
     &  (savedVars, iday, idaynt, numZones, numGridCellsInBlock, gasChemistryType, doCellChem, ntLoopNcs,  &
     &   prate, ifSun, nblockuse, ntloopuse, jlowvar, ktlpvar, jreorder)

      use GmiPrintError_mod, only : GmiPrintError
      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

      implicit none

#     include "smv2chem_par.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: iday
      integer, intent(in)  :: idaynt
      integer, intent(in)  :: numZones
      integer, intent(in)  :: numGridCellsInBlock
      integer, intent(in)  :: gasChemistryType
      logical, intent(in)  :: doCellChem(numZones)
      integer, intent(in)  :: ntLoopNcs(ICS)
      real*8,  intent(in)  :: prate(numZones, IPHOT)

      integer, intent(out) :: ifSun
      integer, intent(out) :: nblockuse
      integer, intent(out) :: ntloopuse
      integer, intent(out) :: jlowvar (MXBLOCK)
      integer, intent(out) :: ktlpvar (MXBLOCK)
      integer, intent(inout) :: jreorder(numZones)
      type(t_Smv2Saved), intent(inOut) :: savedVars


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: iavblok, iavgsize
      integer :: iremain
      integer :: jloop
      integer :: jlooplo  ! low ntloop grid-cell - 1 in a grid-block
      integer :: kblk
      integer :: nblock1

      real*8, parameter :: PHOTOLYSIS_THRESHOLD = 1.0d-80

!     ----------------
!     Begin execution.
!     ----------------

#ifndef nonZeroInd_tracers
      if (idaynt == 1) then

        ntloopuse = 0

        do jloop = 1, ntLoopNcs(gasChemistryType)
          if (doCellChem(jloop)) then
            ntloopuse           = ntloopuse + 1
            jreorder(ntloopuse) = jloop
          end if
        end do

      else

        ntloopuse = 0

        if (iday == 1) then

          ifSun = 1

          do jloop = 1, ntLoopNcs(gasChemistryType)
            if ((prate(jloop,1) > PHOTOLYSIS_THRESHOLD) .and.  &
     &          doCellChem(jloop)) then
              ntloopuse           = ntloopuse + 1
              jreorder(ntloopuse) = jloop
            end if
          end do

        else

          ifSun = 2

          do jloop = 1, ntLoopNcs(gasChemistryType)
            if ((prate(jloop,1) < PHOTOLYSIS_THRESHOLD) .and.  &
     &          doCellChem(jloop)) then
              ntloopuse           = ntloopuse + 1
              jreorder(ntloopuse) = jloop
            end if
          end do

        end if

      end if


      nblockuse = 1 + ntloopuse / (numGridCellsInBlock    + 0.0001d0)
      iavblok   = 1 + ntloopuse / (nblockuse + 0.0001d0)
      iavgsize  = Min (iavblok, numGridCellsInBlock)
      jlooplo   = 0
      nblock1   = nblockuse - 1


      if (nblockuse > savedVars%nblockuse_max) then
        savedVars%nblockuse_max = nblockuse
!c      Write (6,*) 'nblockuse_max:  ', savedVars%nblockuse_max
      end if

      if (savedVars%nblockuse_max > MXBLOCK) then
        Write (6,*) 'nblockuse > MXBLOCK', nblockuse
        call GmiPrintError ('Deter_Block_Size', .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if


      do kblk = 1, nblock1
        jlowvar(kblk) = jlooplo
        ktlpvar(kblk) = iavgsize
        jlooplo       = jlooplo + iavgsize
      end do

      iremain = Max (ntloopuse - (nblock1 * iavgsize), 0)

      jlowvar(nblockuse) = jlooplo
      ktlpvar(nblockuse) = iremain

#endif

      return

      end subroutine Deter_Block_Size


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Solve_Block
!
! DESCRIPTION
!   This routine solves the chemical odes for each block.
!
! ARGUMENTS
!   doQqjkInChem   : if prQqjk is on, should qqj's & qqk's be determined
!                      inside the chemistry solver, or outside?
!   doSurfEmissInChem : do surface emissions inside the chemistry solver, or
!                      outside?
!   prQqjk   : should the periodic qqjk output file be written?
!   prSmv2   : should the SmvgearII     output file be written
!               (non-parallel mode only)?
!   ifSun     : 1 => for daytime, 2 => for nighttime
!   numLats      : # of latitudes
!   numLongs     : # of longitudes
!   numVert     : # of vertical layers
!   airDensityIndex     : array index for air density
!   nitrogenSpNum : identifies spc # of nitrogen gas
!   oxygenSpNum   : identifies spc # of oxygen   gas
!   ireord    : 1 => reorder grid-cells and blocks for chemistry
!               2 => solve chemistry
!   numZones    : # of zones (numLongs * numLats * numVert)
!   smvUnitNumber    : logical unit number to write to when prSmv2 is true
!   nblockuse : (ireord == loreord) => # of original  blocks
!               (ireord /= loreord) => # of reordered blocks
!   gasChemistryType       : identifies gas chemistry type (1.NCSGAS)
!   fractionDecrease   : fraction time step is decreased in Smvgear if convergence
!               test fails
!   maxTimeStepNight   : max time step for night all chem (s)
!   ncOutPeriod     : NetCDF output period
!   modelTimeStep       : model time step (s)
!   doCellChem     : do chemistry for a particular cell?
!   jPhotRate  : tbd
!   numKineticRxns    : # of kinetic rxns (non-photo)
!   numActInActGases    : # of active + inactive gases
!   jlowvar   : tbd
!   ktlpvar   : tbd
!   origSpcNumber   : original spc # of each new jnew spc
!   npphotrat : tbd
!   arate     : thermal    rate constants (units vary)
!   prate     : photolysis rate constants (s^-1)
!   yemis     : surface emissions (molec/cm^3/sec)
!   jreorder  : gives original grid-cell from re-ordered grid-cell
!   lreorder  : gives original grid-cell from re-ordered cell, except when
!               cell is a virtual boundary cell, then it gives original
!               edge cell from re-ordered virtual boundary cell
!   errmx2    : tbd
!   cx        : spc conc (molec/cm^3)
!
!-----------------------------------------------------------------------------

      subroutine Solve_Block  &
     &  (savedVars, doQqjkInChem, doSurfEmissInChem, prQqjk, prSmv2, ifSun,  &
     &   numLats, numLongs, numVert, airDensityIndex, nitrogenSpNum, oxygenSpNum, ireord,  &
     &   numZones, smvUnitNumber, nblockuse, gasChemistryType, fractionDecrease, maxTimeStepNight,  &
     &   ncOutPeriod, modelTimeStep, doCellChem, jPhotRate, numKineticRxns, numActInActGases,  &
     &   jlowvar, ktlpvar, origSpcNumber, npphotrat, arate, prate, yemis,  &
     &   jreorder, lreorder, errmx2, cx, yda, qqkda, qqjda, &
     &   i1, i2, ju1, j2, k1, k2, num_qks, num_qjs, num_active, &
     &   commuWorld)

      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

      implicit none

!     ==========================
#     include "mpif.h"
!     ==========================
#     include "smv2chem_par.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: commuWorld
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: num_qks, num_qjs, num_active
      real*8 , intent(inout) :: qqjda(i1:i2, ju1:j2, k1:k2, num_qjs)
      real*8 , intent(inout) :: qqkda(i1:i2, ju1:j2, k1:k2, num_qks)
      real*8 , intent(inout) :: yda  (i1:i2, ju1:j2, k1:k2, num_active)
      logical, intent(in)  :: doQqjkInChem
      logical, intent(in)  :: doSurfEmissInChem
      logical, intent(in)  :: prQqjk
      logical, intent(in)  :: prSmv2
      integer, intent(in)  :: ifSun
      integer, intent(in)  :: numLats, numLongs, numVert
      integer, intent(in)  :: airDensityIndex
      integer, intent(in)  :: nitrogenSpNum
      integer, intent(in)  :: oxygenSpNum
      integer, intent(in)  :: ireord
      integer, intent(in)  :: numZones
      integer, intent(in)  :: smvUnitNumber
      integer, intent(in)  :: nblockuse
      integer, intent(in)  :: gasChemistryType
      real*8,  intent(in)  :: fractionDecrease
      real*8,  intent(in)  :: maxTimeStepNight
      real*8,  intent(in)  :: ncOutPeriod
      real*8,  intent(in)  :: modelTimeStep
      logical, intent(in)  :: doCellChem(numZones)
      integer, intent(in)  :: jPhotRate (ICS)
      integer, intent(in)  :: numKineticRxns   (ICS)
      integer, intent(in)  :: numActInActGases   (ICS)
      integer, intent(in)  :: jlowvar  (MXBLOCK)
      integer, intent(in)  :: ktlpvar  (MXBLOCK)
      integer, intent(in)  :: origSpcNumber  (MXGSAER, ICS)
      integer, intent(in)  :: npphotrat(IPHOT,   ICS)
      real*8,  intent(in)  :: arate    (numZones,  ITHERM)
      real*8,  intent(in)  :: prate    (numZones,  IPHOT)
      real*8,  intent(in)  :: yemis    (numLats*numLongs, IGAS)

      integer, intent(inout) :: jreorder(numZones)
      integer, intent(inout) :: lreorder(numZones)
      real*8,  intent(inout) :: errmx2  (numZones)
      real*8,  intent(inout) :: cx      (numZones, IGAS)

      type(t_Smv2Saved), intent(inOut) :: savedVars

!     -----------------------
!     Parameter declarations.
!     -----------------------

      logical, parameter :: DOWRT_SBDIAG = .false.


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ierr
      integer :: j
      integer :: jgas
      integer :: jloop
      integer :: jlooplo  ! low ntloop grid-cell - 1 in a grid-block
      integer :: jnew
      integer :: kblk
      integer :: kloop
      integer :: ktloop   ! # of grid-cells in a grid-block
      integer :: nallr    ! # of active rxns
      integer :: nfdh1    ! nfdh2 + # of rxns with one   active reactant
      integer :: nfdh2    ! nfdh3 + # of rxns with two   active reactants
      integer :: nfdh3    !         # of rxns with three active reactants
      integer :: nfdl1    ! nfdh2 + 1
      integer :: nfdl2    ! nfdh3 + 1
      integer :: nfdrep   ! nfdh3 + # of rxns with two active reactants that
                          ! are not followed by a rxn with the same reactants
      integer :: nfdrep1  ! nfdrep + 1
      integer :: np
      integer :: proc_num

!     ------------------------------------------------------------------
!     irma,b,c : spc # of each reactant; locates reordered active spc #s
!     ------------------------------------------------------------------

      integer :: irma(NMTRATE)
      integer :: irmb(NMTRATE)
      integer :: irmc(NMTRATE)

      real*8  :: denair(MXBLOCK)  ! density of air (molec/cm^3)

!     -----------------------------------------------------------------------
!     cc2    : array holding values of decomposed matrix
!     pratk1 : tbd
!     cblk   : gas-phase concs (molec/cm^3)
!     cnew   : stores conc (y (estimated)) (molec/cm^3)
!     corig  : original gas-phase concs used to restart Smvgear if a failure
!              occurs (molec/cm^3)
!     gloss  : value of first derivatives on output from Subfun; right-side
!              of eqn on input to Backsub; error term (solution from Backsub)
!              on output from Backsub
!     smvdm  : amount added to each spc at each grid-cell, for mass balance
!              accounting (# cm^-3 for gas chemistry (?))
!     vdiag  : 1 / current diagonal term of the decomposed matrix
!     rrate  : rate constants
!              (s^-1, cm^3/molec*s, cm^6/molec*s^2, or cm^9/molec*s^3)
!     trate  : rxn rate (moles l^-1-h2o s^-1 or # cm^-3 s^-1 (?))
!     -----------------------------------------------------------------------

      real*8  :: cc2   (KBLOOP, 0:MXARRAY)  = 0.0d0

      real*8  :: pratk1(KBLOOP, IPHOT)      = 0.0d0

      real*8  :: cblk  (KBLOOP, MXGSAER)    = 0.0d0
      real*8  :: cnew  (KBLOOP, MXGSAER)    = 0.0d0
      real*8  :: corig (KBLOOP, MXGSAER)    = 0.0d0
      real*8  :: gloss (KBLOOP, MXGSAER)    = 0.0d0
      real*8  :: smvdm (KBLOOP, MXGSAER)    = 0.0d0
      real*8  :: vdiag (KBLOOP, MXGSAER)    = 0.0d0

      real*8  :: rrate (KBLOOP, NMTRATE)    = 0.0d0
      real*8  :: dummy (KBLOOP, NMTRATE)    = 0.0d0

!     ===============================
!$    integer :: Omp_Get_Max_Threads
!$    integer :: Omp_Get_Num_Threads
!c    integer :: Omp_Get_Thread_Num

!$    integer :: max_thrds, num_thrds
!c    integer :: thrd_num
!     ===============================

!     =========================
!c    real*8  :: wclk1, wclk2
!c    real*8  :: wtime
!c    real*8  :: ctime(MXBLOCK)
!c    real*8  :: sumtime(0:15)
!     =========================


!     ----------------
!     Begin execution.
!     ----------------

      if (DOWRT_SBDIAG) then
!       ====================================================================
        call Mpi_Comm_Rank (commuWorld, proc_num, ierr)
        proc_num = proc_num + 1  ! change proc_num range from 0->N to 1->N+1

!$omp   parallel
!$      num_thrds = Omp_Get_Num_Threads ( )
!$omp   end parallel

!$      max_thrds = Omp_Get_Max_Threads ( )
!$      Write (6,900) proc_num, max_thrds, num_thrds, nblockuse
900   format ('Proc #, Max Thrds, # Thrds, nblockuse:  ', i6, i4, i4, i8)
!       ====================================================================
      end if

!c    call f_hpmstart (1, "chem block loop")

!c    call Mpi_Comm_Rank (commuWorld, proc_num, ierr)
!c    ctime  (:) = 0.0d0
!c    sumtime(:) = 0.0d0

!$omp   parallel do &
!$omp&  default(shared) &
!$omp&  schedule(runtime) &
!$omp&  private(j, jgas, jnew, kblk, np) &
!$omp&  private(jloop,jlooplo,kloop) &
!$omp&  firstprivate(ktloop) &
!$omp&  firstprivate(nallr, nfdrep, nfdrep1) &
!$omp&  firstprivate(nfdh1, nfdh2, nfdh3, nfdl1, nfdl2) &
!$omp&  firstprivate(cblk, cc2, cnew, corig) &
!$omp&  firstprivate(denair, gloss) &
!$omp&  firstprivate(irma, irmb, irmc) &
!$omp&  firstprivate(pratk1, smvdm, vdiag) &
!$omp&  firstprivate(dummy,rrate)


!     ======================
      do kblk = 1, nblockuse
!     ======================

!c      thrd_num  = Omp_Get_Thread_Num  ( )
!c      num_thrds = Omp_Get_Num_Threads ( )
!c      wclk1     = Mpi_Wtime (ierr)

        jlooplo = jlowvar(kblk)
        ktloop  = ktlpvar(kblk)

!       ================
        if (ktloop /= 0) then
!       ================

!         -------------------------------------
!         Set (and rearrange) photofrequencies.
!         -------------------------------------

          do j = 1, jPhotRate(gasChemistryType)
            np = npphotrat(j,gasChemistryType)
            do kloop = 1, ktloop
              jloop           = lreorder(jlooplo+kloop)
              pratk1(kloop,j) = prate(jloop,np)
            end do
          end do

!         -------------------------------------------------------
!         Place large domain gas array (molec/cm^3) into smaller
!         block array.
!
!         The urate=arate loop is where the 2D diurnal averaging
!         factors for the thermal reactions are passed in.  urate
!         later is used for another purpose.
!         K: urate replaced by "dummy"
!         -------------------------------------------------------

          do j = 1, numKineticRxns(gasChemistryType)
            do kloop = 1, ktloop
              jloop            = lreorder(jlooplo+kloop)
              dummy(kloop,j)   = arate(jloop,j)
            end do
          end do

!         --------------------------------------------------------------
!         Calculate rates and solve chemistry.
!
!         ireord = 1 : call Calcrate to find stiffness of each grid-cell
!         ireord = 2 : set chemistry rates and solve equations
!         --------------------------------------------------------------

!         =============
          call Calcrate  &
!         =============
     &      (ifSun, airDensityIndex, nitrogenSpNum, oxygenSpNum, numZones, jlooplo,  &
     &       ktloop, gasChemistryType, jreorder, numKineticRxns, numActInActGases, cx, dummy,  &
     &       denair, pratk1, cblk, rrate, nallr, nfdh1, nfdh2,  &
     &       nfdh3, nfdl1, nfdl2, nfdrep, nfdrep1, irma, irmb, irmc,  &
     &       corig, smvdm, savedVars)

!         --------------------
!         Solve chemical odes.
!         --------------------

!         ============
          call Smvgear  &
!         ============
     &      (savedVars, doQqjkInChem, doSurfEmissInChem, prQqjk, prSmv2,  &
     &       ifSun, numLats, numLongs, numVert, ireord, numZones, jlooplo,  &
     &       ktloop, smvUnitNumber, nallr, gasChemistryType, nfdh2, nfdh3, nfdl1, nfdl2,  &
     &       nfdrep, nfdrep1, fractionDecrease, maxTimeStepNight, ncOutPeriod, modelTimeStep,  &
     &       doCellChem, irma, irmb, irmc, jreorder, jPhotRate,  &
     &       numActInActGases, origSpcNumber, denair, corig, pratk1, yemis, smvdm,  &
     &       nfdh1, errmx2, cc2, cnew, gloss, vdiag, rrate,  &
     &       yda, qqkda, qqjda, &
     &       i1, i2, ju1, j2, k1, k2, num_qks, num_qjs, num_active)

!         -----------------------------------------------------
!         Replace block concentrations (molec/cm^3) into domain
!         concentrations, but only for active species.
!         -----------------------------------------------------

          if (ireord == 2) then

            do jnew = 1, numActInActGases(gasChemistryType)
              jgas = origSpcNumber(jnew,gasChemistryType)

              if (jgas <= SK_NACT) then

                do kloop = 1, ktloop
                  jloop          = jreorder(jlooplo+kloop)
                  cx(jloop,jgas) = Max (cnew(kloop,jnew), SMAL2)
                end do

              end if

            end do

          end if

!       ======
        end if
!       ======

!c      wclk2       = Mpi_Wtime (ierr)
!c      ctime(kblk) = wclk2 - wclk1
!c      Write (6,910) proc_num, thrd_num, kblk, ctime(kblk)
!910    format ('Proc #, Thrd, Block, Time:  ', i6, i4, i8, e20.8)
!c      sumtime(thrd_num) = sumtime(thrd_num) + ctime(kblk)

!     ======
      end do
!     ======

!c    call f_hpmstop (1)

!c    Write (6,920)
!c   &  proc_num, (sumtime(thrd_num), thrd_num=0,num_thrds-1)
!920  format ('Proc #, sumtime:  ', i6, 16e20.8)

      return

      end subroutine Solve_Block


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Calcrate
!
! DESCRIPTION
!   This routine calculates kinetic reaction and photorates
!   (s^-1, cm^3/s, or cm^6/s^2) and pressure- and temperature- dependence
!   for gas-phase chemical reactions.
!
! ARGUMENTS
!   ifSun     : 1 => for daytime, 2 => for nighttime
!   airDensityIndex     : array index for air density
!   nitrogenSpNum : identifies spc # of nitrogen gas
!   oxygenSpNum   : identifies spc # of oxygen   gas
!   numZones    : # of zones (numLongs * numLats * numVert)
!   jlooplo   : low ntloop grid-cell - 1 in a grid-block
!   ktloop    : # of grid-cells in a grid-block
!   gasChemistryType       : identifies gas chemistry type (1.NCSGAS)
!   jreorder  : gives original grid-cell from re-ordered grid-cell
!   numKineticRxns    : # of kinetic rxns (non-photo)
!   numActInActGases    : # of active + inactive gases
!   cx        : spc conc (molec/cm^3)
!   urate     : term of Jacobian (J) = partial derivative
!   denair    : density of air  (molec/cm^3)
!   pratk1    : tbd
!   cblk      : gas-phase concs (molec/cm^3)
!   rrate     : rate constants
!               (s^-1, cm^3/molec*s, cm^6/molec*s^2, or cm^9/molec*s^3)
!   trate     : rxn rate (moles l^-1-h2o s^-1 or # cm^-3 s^-1 (?))
!   nallr     : # of active rxns
!   nfdh1     : nfdh2 + # of rxns with one   active reactant
!   nfdh2     : nfdh3 + # of rxns with two   active reactants
!   nfdh3     :         # of rxns with three active reactants
!   nfdl1     : nfdh2 + 1
!   nfdl2     : nfdh3 + 1
!   nfdrep    : nfdh3 + # of rxns with two active reactants that are not
!               followed by a rxn with the same reactants
!   nfdrep1   : nfdrep + 1
!   irma,b,c  : spc # of each reactant; locates reordered active spc #s
!   corig     : original gas-phase concs used to restart Smvgear if a failure
!               occurs (molec/cm^3)
!   smvdm     : amount added to each spc at each grid-cell, for mass balance
!               accounting (# cm^-3 for gas chemistry (?))
!
!-----------------------------------------------------------------------------

      subroutine Calcrate  &
     &  (ifSun, airDensityIndex, nitrogenSpNum, oxygenSpNum, numZones, jlooplo, ktloop,  &
     &   gasChemistryType, jreorder, numKineticRxns, numActInActGases, cx, dummy, denair, pratk1,  &
     &   cblk, rrate, nallr, nfdh1, nfdh2, nfdh3, nfdl1, nfdl2,  &
     &   nfdrep, nfdrep1, irma, irmb, irmc, corig, smvdm, savedVars)

      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

      implicit none

#     include "smv2chem_par.h"
!#     include "smv2chem2.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: ifSun
      integer, intent(in)  :: airDensityIndex
      integer, intent(in)  :: nitrogenSpNum
      integer, intent(in)  :: oxygenSpNum
      integer, intent(in)  :: numZones
      integer, intent(in)  :: jlooplo
      integer, intent(in)  :: ktloop
      integer, intent(in)  :: gasChemistryType
      integer, intent(in)  :: jreorder(numZones)
      integer, intent(in)  :: numKineticRxns  (ICS)
      integer, intent(in)  :: numActInActGases  (ICS)
      real*8,  intent(in)  :: cx      (numZones, IGAS)
      !K: Fix the sloppy use of dummy
      real*8,  intent(inout)  :: dummy (KBLOOP, NMTRATE)

      real*8,  intent(inout) :: denair(ktloop)
      real*8,  intent(inout) :: pratk1(KBLOOP, IPHOT)
      real*8,  intent(inout) :: cblk  (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: rrate (KBLOOP, NMTRATE)

      integer, intent(out) :: nallr
      integer, intent(out) :: nfdh1,  nfdh2, nfdh3
      integer, intent(out) :: nfdl1,  nfdl2
      integer, intent(out) :: nfdrep, nfdrep1
      integer, intent(out) :: irma (NMTRATE)
      integer, intent(out) :: irmb (NMTRATE)
      integer, intent(out) :: irmc (NMTRATE)
      real*8,  intent(out) :: corig(KBLOOP, MXGSAER)
      real*8,  intent(out) :: smvdm(KBLOOP, MXGSAER)

      type(t_Smv2Saved), intent(inOut) :: savedVars

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: i, j
      integer :: jloop
      integer :: jnew
      integer :: jold            ! = mappl(jold) for inactive spc
      integer :: kloop
      integer :: ncsp            ! gasChemistryType       => for daytime   gas chemistry
                                 ! gasChemistryType + ICS => for nighttime gas chemistry
      integer :: nfdl0           ! nfdh1 + 1
      integer :: nh, nk, nkn

      real*8  :: concn2(ktloop)  ! nitrogen conc  (molec/cm^3)
      real*8  :: conco2(ktloop)  ! oxygen   conc  (molec/cm^3)


!     ----------------
!     Begin execution.
!     ----------------

!     ----------------------------
!     Load kinetic reaction rates.
!     ----------------------------

      do nk = 1, numKineticRxns(gasChemistryType)
        do kloop = 1, ktloop
          rrate(kloop,nk) = dummy(kloop,nk)
        end do
      end do


!     -------------------------------------------------------------------
!     Place large domain gas array (molec/cm^3) into smaller block array.
!     -------------------------------------------------------------------

      do jold = 1, numActInActGases(gasChemistryType)

        jnew = savedVars%mappl(jold,gasChemistryType)

        do kloop = 1, ktloop

          jloop = jreorder(jlooplo+kloop)

          cblk (kloop,jold) = cx(jloop,jold)
          corig(kloop,jnew) = cx(jloop,jold)
          smvdm(kloop,jnew) = 0.0d0

!         ---------------------------------
!         Load third body number densities.
!         ---------------------------------

          if (nitrogenSpNum > 0) concn2(kloop) = cx(jloop,nitrogenSpNum)
          if (oxygenSpNum   > 0) conco2(kloop) = cx(jloop,oxygenSpNum)
          if (airDensityIndex     > 0) denair(kloop) = cx(jloop,airDensityIndex)

        end do

      end do


!     ------------------------------------------------------
!     Multiply rates by constant species concentrations
!     (either M, O2, N2, or any active or inactive species).
!     ------------------------------------------------------

      do i = 1, savedVars%nmair(gasChemistryType)
        nk = savedVars%nreacair(i,gasChemistryType)
        do kloop = 1, ktloop
           rrate(kloop,nk) = rrate(kloop,nk) * denair(kloop)
        end do
      end do

      do i = 1, savedVars%nmo2(gasChemistryType)
        nk = savedVars%nreaco2(i,gasChemistryType)
        do kloop = 1, ktloop
          rrate(kloop,nk) = rrate(kloop,nk) * conco2(kloop)
        end do
      end do

      do i = 1, savedVars%nmn2(gasChemistryType)
        nk = savedVars%nreacn2(i,gasChemistryType)
        do kloop = 1, ktloop
          rrate(kloop,nk) = rrate(kloop,nk) * concn2(kloop)
        end do
      end do


!     ---------------------------------------------------------------
!     Multiply rate coefficient by any other third body concentration
!     (e.g., H2O); multiply by other inactive concentrations later.
!     ---------------------------------------------------------------

      do i = 1, savedVars%nm3bod(gasChemistryType)
        nk   = savedVars%nreac3b (i,gasChemistryType)
        jold = savedVars%lgas3bod(i,gasChemistryType)
        do kloop = 1, ktloop
          rrate(kloop,nk) = rrate(kloop,nk) * cblk(kloop,jold)
        end do
      end do


!     -----------------------------------------------------------
!     Multiply rate coefficient by other inactive concentrations.
!     This loop must occur after equilibrium reactions.
!     -----------------------------------------------------------

      do i = 1, savedVars%nmoth(gasChemistryType)
        nk   = savedVars%nreacoth(i,gasChemistryType)
        jold = savedVars%lgasbino(i,gasChemistryType)
        do kloop = 1, ktloop
          rrate(kloop,nk) = rrate(kloop,nk) * cblk(kloop,jold)
       end do
      end do


!     --------------------
!     Reorder rrate array.
!     --------------------

      nfdh3   = savedVars%ithrr(gasChemistryType)
      nfdl2   = nfdh3  + 1
      nfdrep  = savedVars%inorep(gasChemistryType)
      nfdrep1 = nfdrep + 1
      nfdh2   = nfdh3  + savedVars%itwor(gasChemistryType)
      nfdl1   = nfdh2  + 1
      nfdh1   = nfdh2  + savedVars%ioner(gasChemistryType)
      nfdl0   = nfdh1  + 1
      nallr   = savedVars%nallrat(gasChemistryType)

      do nkn = 1, nallr
        nk = savedVars%noldfnew(nkn,gasChemistryType)
        irma(nkn) = savedVars%irm2(1,nk,gasChemistryType)
        irmb(nkn) = savedVars%irm2(2,nk,gasChemistryType)
        irmc(nkn) = savedVars%irm2(3,nk,gasChemistryType)
      end do


!     ---------------------------------
!     trate here used as a dummy array.
!     ---------------------------------

      do nk = 1, savedVars%ntrates(gasChemistryType)
        do kloop = 1, ktloop
          dummy(kloop,nk) = rrate(kloop,nk)
        end do
      end do

      do nkn = 1, nallr
        nk = savedVars%noldfnew(nkn,gasChemistryType)
        do kloop = 1, ktloop
          rrate(kloop,nkn) = dummy(kloop,nk)
        end do
      end do




!     ---------------------------------------------------------------
!     Reset ncsp; multiply photorate coef. (s^-1) by concentration
!     (molec/cm^3) when photodissociating species is inactive.  Thus,
!     rate does not need to be multiplied later by concentration.
!     ---------------------------------------------------------------

      ncsp = (ifSun - 1) * ICS + gasChemistryType

      do i = 1, savedVars%nolosp(ncsp)
        nk   = savedVars%nknlosp(i,gasChemistryType)
        j    = savedVars%jphotnk(nk,gasChemistryType)
        jold = savedVars%losinacp(i,gasChemistryType)
        do kloop = 1, ktloop
          pratk1(kloop,j)  = pratk1(kloop,j) * cblk(kloop,jold)
        end do
      end do


      return

      end subroutine Calcrate


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Reorder_Grid_Cells
!
! DESCRIPTION
!   This routine reorders the  grid-cells from least to most stiff, with
!   those with similar stiffness grouped together.
!
!   From least to most stiff:
!     smaller errmx2 (csuma) => less stiff
!     (csuma = errmx2)
!
!   Sort using heapsort routine (numerical recipes), an n(logb2)n process.
!   This reordering scheme is very fast, although complicated.  errmx2 from
!   Smvgear:  denotes stiffness (larger value => more stiff).
!
! ARGUMENTS
!   numZones    : # of zones (numLongs * numLats * numVert)
!   numGridCellsInBlock    : intended # of grid-cells in a grid-block
!   ntloopuse : tbd
!   errmx2    : tbd
!   jreorder  : gives original grid-cell from re-ordered grid-cell
!   csuma     : tbd
!   nreblock  : tbd
!   lreorder  : gives original grid-cell from re-ordered cell, except when
!               cell is a virtual boundary cell, then it gives original
!               edge cell from re-ordered virtual boundary cell
!   jlowvar   : tbd
!   ktlpvar   : tbd
!   csumc     : tbd
!
!-----------------------------------------------------------------------------

      subroutine Reorder_Grid_Cells  &
     &  (numZones, numGridCellsInBlock, ntloopuse, errmx2, jreorder, csuma,  &
     &   nreblock, lreorder, jlowvar, ktlpvar, csumc)

      implicit none

#     include "smv2chem_par.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: numZones
      integer, intent(in)  :: numGridCellsInBlock
      integer, intent(in)  :: ntloopuse
      real*8,  intent(in)  :: errmx2(numZones)

      integer, intent(inout) :: jreorder(numZones)

      integer, intent(out) :: nreblock
      integer, intent(out) :: lreorder(numZones)
      integer, intent(out) :: jlowvar (MXBLOCK)
      integer, intent(out) :: ktlpvar (MXBLOCK)
      real*8,  intent(out) :: csuma(numZones)
      real*8,  intent(out) :: csumc(numZones)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: iavblok, iavgsize
      integer :: ipar, jpar, jpar1
      integer :: iradd, iradd1
      integer :: iremain
      integer :: irval, lval
      integer :: jllast
      integer :: jloop
      integer :: jlooplo    ! low ntloop grid-cell - 1 in a grid-block
      integer :: jreord
      integer :: kblk
      integer :: nblock1
      integer :: nblockrow  ! # blocks for each reorder group
                            ! (stiffness, sunrise, sunset)
      integer :: ncellrow
      integer :: nnorise

      real*8  :: vallow


!     ----------------
!     Begin execution.
!     ----------------

      do jloop = 1, ntloopuse
        lreorder(jloop) = jreorder(jloop)
        csumc   (jloop) = errmx2  (jloop)
      end do


      jllast = ntloopuse

      do jloop = 1, ntloopuse
        jreorder(jllast) = lreorder(jloop)
        csuma   (jllast) = csumc   (jloop)
        jllast           = jllast - 1
      end do


      nnorise = ntloopuse

      iradd   = 0
      lval    = iradd + (nnorise * 0.5d0) + 1
      irval   = iradd + nnorise
      iradd1  = iradd + 1

!     ===================
      if (irval > iradd1) then
!     ===================

!       =============
        OUTERLOOP: do
!       =============

          if (lval > iradd1) then

            lval            = lval - 1
            vallow          = csuma   (lval)
            jreord          = jreorder(lval)

          else

            vallow          = csuma   (irval)
            jreord          = jreorder(irval)
            csuma   (irval) = csuma   (iradd1)
            jreorder(irval) = jreorder(iradd1)
            irval           = irval - 1

            if (irval == iradd1) then
              csuma   (iradd1) = vallow
              jreorder(iradd1) = jreord
!             ==============
              exit OUTERLOOP
!             ==============
            end if

          end if

          ipar = lval
          jpar = lval + lval - iradd

!         ===================================
          INNERLOOP: do while (jpar <= irval)
!         ===================================

            if (jpar < irval) then
              jpar1 = jpar + 1
              if (csuma(jpar) < csuma(jpar1)) jpar = jpar1
            end if

            if (vallow < csuma(jpar)) then
              csuma   (ipar) = csuma   (jpar)
              jreorder(ipar) = jreorder(jpar)
              ipar           = jpar
              jpar           = jpar + jpar - iradd
!             ===============
              cycle INNERLOOP
!             ===============
            else
!             ==============
              exit INNERLOOP
!             ==============
            end if

!         ================
          end do INNERLOOP
!         ================

          csuma   (ipar) = vallow
          jreorder(ipar) = jreord

!       ================
        end do OUTERLOOP
!       ================

!     ======
      end if
!     ======


!     ---------------------------------------------------------------
!     Determine how many blocks of cells are needed after reordering.
!     ---------------------------------------------------------------

      jlooplo   = 0
      nblockrow = 0
      nreblock  = 0

      ncellrow = nnorise

      if (ncellrow /= 0) then

        nblockrow = 1 + (ncellrow / (numGridCellsInBlock + 0.0001d0))

        iavblok   = 1 + (ncellrow / (nblockrow + 0.0001d0))
        iavgsize  = Min (iavblok, numGridCellsInBlock)
        nblock1   = nblockrow - 1

        do kblk = 1, nblock1
          nreblock          = nreblock + 1
          jlowvar(nreblock) = jlooplo
          ktlpvar(nreblock) = iavgsize
          jlooplo           = jlooplo + iavgsize
        end do

        nreblock = nreblock + 1
        iremain  = Max (ncellrow - (nblock1 * iavgsize), 0)

        jlowvar(nreblock) = jlooplo
        ktlpvar(nreblock) = iremain

      end if


      return

      end subroutine Reorder_Grid_Cells
