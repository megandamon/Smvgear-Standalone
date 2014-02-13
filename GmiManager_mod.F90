! MRD: Notes from Tom
! two levels of management
! top level : GMI specific reordering
! below: list of boxes - focused on managing gear updates its order and timestep


module GmiManager_mod

   use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

   implicit none
   private

#     include "smv2chem_par.h"

   public :: Manager_type
   public :: resetGear
   public :: startTimeInterval
   public :: calculateErrorTolerances
   public :: calcInitialTimeStepSize
   public :: updateCoefficients
   public :: calculateTimeStep
   public :: tightenErrorTolerance
   public :: calculateNewRmsError
   public :: testAccumulatedError
   public :: estimateTimeStepRatio
   public :: resetBeforeUpdate
   public :: calcNewAbsoluteErrorTolerance
   public :: scaleDerivatives
   public :: updateChold
   public :: predictConcAndDerivatives
   public :: resetCnewDerivatives
   public :: updateDerivatives


! MRD: add type bound procedures here
! can remove "_type"
! need a constructor
   type Manager_type

       ! private
       integer :: numErrTolDecreases
       integer :: numFailOldJacobian ! of times corrector failed to converge while the Jacobian was old
       integer :: jFail ! of times correcter failed to converge after old "Pderiv" was called
       integer :: numFailErrorTest
       integer :: numFailAfterPredict
       integer :: numCallsPredict ! total # of times predictor is called
       integer :: numSuccessTdt ! numSuccessTdt    : total # of successful time steps taken
       integer :: numCallsVelocity ! total # of times velocity is called
       integer :: num1stOEqnsSolve !# of first-order eqns to solve, = # of spc = order of
               ! original matrix; num1stOEqnsSolve has a different value for day and
               ! night and for gas- and aqueous-phase chemistry;
               !  # spc with prod or loss terms in Smvgear (?)
       real*8  :: abtoler1
       real*8  :: order, order_inv
       real*8  :: chemTimeInterval !total chem time interval; same as chemintv (s)
       real*8  :: maxTimeStep ! max time step at a given time (s)
       real*8  :: failureFraction ! = 1 originially, but is decreased if excessive failures
               ! occur in order to reduce absolute error tolerance
       real*8  :: timeremain ! remaining time in an chem interval (s)
       real*8  :: iabove
       real*8  :: initialError, initialError_inv
       real*8  :: errmax_ncs_inv
       real*8  :: xelaps ! elapsed time in chem interval (s)
       real*8  :: told !stores last value of xelaps in case current step fails
       real*8  :: reltol1, reltol2, reltol3
       real*8  :: rmsError
       integer :: idoub ! records # of steps since the last change in step size or
               ! order; it must be at least kstep = nqq+1 before doubling is
               ! allowed
       integer :: nslp ! last time step # during which "Pderiv" was called
       integer :: jrestar ! counts # of times Smvgear starts over at order 1 because of
               ! excessive failures

       integer :: nqqold ! value of nqq during last time step
       integer :: nqq ! order of integration method; varies between 1 and MAXORD
       integer :: kstep ! nqq + 1
       real*8  :: hratio ! relative change in delt*aset(1) each change in step or order
               ! when Abs(hratio-1) > MAX_REL_CHANGE, reset jeval = 1 to call Pderiv
       real*8 :: asn1 ! value of aset(nqq,1)
       real*8 :: enqq ! pertst^2*order for current order
       real*8  :: conp1, conp2, conp3
       integer :: nqqisc ! nqq * num1stOEqnsSolve
       real*8 :: rdelt ! factor (time step ratio) by which delt is increased or decreased
       real*8 :: rdelmax ! max factor by which delt can be increased in a single step;
               !                 as in Lsodes, set it to 1d4 initially to compensate for the
               !                 small initial delt, but then set it to 10 after successful
               !                 steps and to 2 after unsuccessful steps
       real*8 :: der2max
       real*8  :: drate ! parameter which is used to determine whether convergence
               !                 has occurred
       real*8  :: dcon
       real*8  :: dtlos (KBLOOP, MXGSAER)! an array of length num1stOEqnsSolve, used for the accumulated corrections;
               ! on a successful return; dtlos(kloop,i) contains the estimated
               ! one step local error in cnew
      real*8  :: chold (KBLOOP, MXGSAER) ! 1 / (reltol * cnew + abtol); multiply chold by local errors in different error tests
      real*8  :: rdeltdn ! time step ratio at one order lower  than current order
      real*8  :: rdeltup ! time step ratio at one order higher than current order
      integer :: ifsuccess ! identifies whether step is successful (=1) or not (=0)
    end type Manager_type

contains

      subroutine updateDerivatives(this, cnewDerivatives, ktloop, savedVars)
         implicit none
         type (Manager_type) :: this
         real*8, intent(inout) :: cnewDerivatives(KBLOOP, MXGSAER*7)
         integer, intent(in) :: ktloop
         type(t_Smv2Saved), intent(inOut) :: savedVars

         integer :: i, i1, j, jspc, kloop
         real*8  :: asnqqj

         i1 = 1
         do j = 2, this%kstep
           i1 = i1 + this%num1stOEqnsSolve
           asnqqj = savedVars%aset(this%nqq,j)
           do jspc = 1, this%num1stOEqnsSolve
             i = jspc + i1 - 1
             do kloop = 1, ktloop
               cnewDerivatives(kloop,i) =  cnewDerivatives(kloop,i) + (asnqqj * this%dtlos(kloop,jspc))
             end do
           end do
         end do
      end subroutine updateDerivatives

      subroutine resetCnewDerivatives(this, cnewDerivatives, ktloop)
         implicit none

         type (Manager_type) :: this
         real*8, intent(inout) :: cnewDerivatives(KBLOOP, MXGSAER*7)
         integer, intent(in) :: ktloop
         integer :: i,i1,j,jb,kloop

         i1 = this%nqqisc + 1
         j = 0
         do jb = 1, this%nqq
            i1 = i1 - this%num1stOEqnsSolve
            do i = i1, this%nqqisc
               j = i + this%num1stOEqnsSolve
               do kloop = 1, ktloop
                  cnewDerivatives(kloop,i) = cnewDerivatives(kloop,i) - cnewDerivatives(kloop,j)
             end do
           end do
         end do

      end subroutine resetCnewDerivatives


!-----------------------------------------------------------------------------
!
! ROUTINE
!   predictConcAndDerivatives
! DESCRIPTION
!     Compute the predicted concentration and derivatives by multiplying
!     previous values by the pascal triangle matrix.
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
      subroutine predictConcAndDerivatives(this, conc, explic, ktloop)

         implicit none

         type (Manager_type) :: this
         real*8, intent(inout) :: conc(KBLOOP, MXGSAER*7)
         real*8, intent(out) :: explic(KBLOOP, MXGSAER)
         integer, intent(in) :: ktloop

         integer :: i,i1,j,jb,jspc,kloop

         !if (prDiag) Write(*,*) "Computing predicted conc and derivatives using pascal triangle matrix"
         i1 = this%nqqisc + 1

         do jb = 1, this%nqq - 1
           i1 = i1 - this%num1stOEqnsSolve
           do i = i1,  this%nqqisc
             j = i + this%num1stOEqnsSolve
             do kloop = 1, ktloop
               conc(kloop,i)  = conc(kloop,i) + conc(kloop,j)
             end do
           end do
         end do

         do jspc = 1,  this%num1stOEqnsSolve
           j = jspc + this%num1stOEqnsSolve
           do kloop = 1, ktloop
             conc  (kloop,jspc) = conc(kloop,jspc) + conc(kloop,j)
             explic(kloop,jspc) = conc(kloop,j)
           end do

         end do

         do i = this%num1stOEqnsSolve + 1, this%nqqisc
           j = i + this%num1stOEqnsSolve
           do kloop = 1, ktloop
             conc(kloop,i) = conc(kloop,i) + conc(kloop,j)
           end do
         end do

      end subroutine predictConcAndDerivatives

!-----------------------------------------------------------------------------
!
! ROUTINE
!   updateChold
! DESCRIPTION
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine updateChold (this, ktloop, cnew, yabst)
      type (Manager_type) :: this
      integer, intent(in) :: ktloop
      real*8,  intent(in) :: cnew  (KBLOOP, MXGSAER)
      real*8, intent(in)  :: yabst (KBLOOP)

      integer :: kloop, jspc

      do kloop = 1, ktloop
       do jspc = 1, this%num1stOEqnsSolve

         this%chold(kloop,jspc) =  &
  &        this%reltol3 /  &
  &        (Max (cnew(kloop,jspc), 0.0d0) +  &
  &         (yabst(kloop) * this%reltol2))

       end do
     end do

   end subroutine updateChold

!-----------------------------------------------------------------------------
!
! ROUTINE
!   scaleDerivatives
! DESCRIPTION
! Created by: Megan Rose Damon
! cnewDerivatives   : an array of length num1stOEqnsSolve*(MAXORD+1) that carries the
!   derivatives of cnew, scaled by delt^j/factorial(j), where j is
!   the jth derivative; j varies from 1 to nqq; e.g., conc(jspc,2)
!   stores delt*y' (estimated)
!-----------------------------------------------------------------------------
   subroutine scaleDerivatives (this, ktloop, cnewDerivatives)
      type (Manager_type) :: this
      integer, intent(in)  :: ktloop
      real*8, intent(inout)  :: cnewDerivatives  (KBLOOP, MXGSAER*7)

      real*8  :: rdelta
      integer :: i1, j, i, kloop

      rdelta = 1.0d0
      i1     = 1

      do j = 2, this%kstep
         rdelta = rdelta * this%rdelt
         i1 = i1 + this%num1stOEqnsSolve
          do i = i1, i1 + (this%num1stOEqnsSolve-1)
            do kloop = 1, ktloop
              cnewDerivatives(kloop,i) = cnewDerivatives(kloop,i) * rdelta
            end do
          end do
      end do

   end subroutine scaleDerivatives

!-----------------------------------------------------------------------------
!
! ROUTINE
!   calcNewAbsoluteErrorTolerance
! DESCRIPTION
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine  calcNewAbsoluteErrorTolerance (this, cnew, concAboveAbtolCount, &
      &  ktloop, yabst, ncs, savedVars)
      type (Manager_type) :: this
      real*8,  intent(in) :: cnew  (KBLOOP, MXGSAER)
      integer, intent(inout) :: concAboveAbtolCount(KBLOOP, 5)
      integer, intent(in)  :: ktloop
      real*8, intent(out)  :: yabst (KBLOOP)
      integer, intent(in)  :: ncs
      type(t_Smv2Saved), intent(inOut) :: savedVars

      integer :: jspc, kloop
      integer :: k1, k2, k3, k4, k5, k
      real*8  :: cnw

      do k = 1, 5
          do kloop = 1, ktloop
            concAboveAbtolCount(kloop,k) = 0
          end do
      end do

      do jspc = 1, this%num1stOEqnsSolve
          do kloop = 1, ktloop
            cnw = cnew(kloop,jspc)
            do k = 1, 5
               if (cnw > savedVars%abtol(k,ncs)) then
                 concAboveAbtolCount(kloop,k) = concAboveAbtolCount(kloop,k) + 1
                 exit
               end if
            end do
          end do
        end do

        do kloop = 1, ktloop

          k1 = concAboveAbtolCount(kloop,1)
          k2 = concAboveAbtolCount(kloop,2) + k1
          k3 = concAboveAbtolCount(kloop,3) + k2
          k4 = concAboveAbtolCount(kloop,4) + k3
          k5 = concAboveAbtolCount(kloop,5) + k4

          if (k1 > this%iabove) then
            yabst(kloop) = savedVars%abtol(1,ncs) ! MRD: these yabst should be passed in
          else if (k2 > this%iabove) then    ! does the driver pass them in?
            yabst(kloop) = savedVars%abtol(2,ncs) ! or does the mechanism specify them
          else if (k3 > this%iabove) then    ! tabled for now
            yabst(kloop) = savedVars%abtol(3,ncs)
          else if (k4 > this%iabove) then
            yabst(kloop) = savedVars%abtol(4,ncs)
          else if (k5 > this%iabove) then
            yabst(kloop) = savedVars%abtol(5,ncs)
          else
            yabst(kloop) = savedVars%abtol(6,ncs)
          end if

        end do

        end subroutine calcNewAbsoluteErrorTolerance
!-----------------------------------------------------------------------------
!
! ROUTINE
!   resetBeforeUpdate
! DESCRIPTION
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine resetBeforeUpdate (this)
      type (Manager_type) :: this

      this%hratio    = 0.0d0
      this%asn1      = 1.0d0
      this%ifsuccess = 1
      this%rdelmax   = 1.0d4

   end subroutine resetBeforeUpdate

!-----------------------------------------------------------------------------
!
! ROUTINE
!   estimateTimeStepRatio
! DESCRIPTION
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine estimateTimeStepRatio (this, ktloop, dely, conc, savedVars)

      ! ----------------------
      ! Argument declarations.
      ! ----------------------
      type (Manager_type) :: this
      integer, intent(in) :: ktloop
      real*8, intent(inout)  :: dely  (KBLOOP)
      real*8, intent(in)  :: conc  (KBLOOP, MXGSAER*7)
      type(t_Smv2Saved), intent(inOut) :: savedVars

      integer :: kloop, kstepisc, jspc, i
      real*8  :: errymax, der1max
      real*8  :: rdeltsm   ! time step ratio at current order

      !     Estimate the time step ratio (rdeltsm) at the current order.
      !     der2max was calculated during the error tests earlier.
      rdeltsm = 1.0d0 / ((this%conp2 * this%der2max**savedVars%enqq2(this%nqq)) + 1.2d-6)

      !     Estimate the time step ratio (rdeltdn) at one order lower than
      !     the current order.  if nqq = 1, then we cannot test a lower order
      if (this%nqq > 1) then

         do kloop = 1, ktloop
            dely(kloop) = 0.0d0
         end do

         kstepisc = (this%kstep - 1) * this%num1stOEqnsSolve

         do kloop = 1, ktloop
            do jspc = 1, this%num1stOEqnsSolve
               i = jspc + kstepisc
               errymax     = conc(kloop,i) * this%chold(kloop,jspc)
               dely(kloop) = dely(kloop) + (errymax * errymax)
            end do
         end do

        der1max = 0.0d0

        do kloop = 1, ktloop
          if (dely(kloop) > der1max) then
            der1max = dely(kloop)
          end if
        end do

        this%rdeltdn = 1.0d0 / &
            ((this%conp1 * der1max**savedVars%enqq1(this%nqq)) + 1.3d-6)

      else
        this%rdeltdn = 0.0d0
      end if

      !     Find the largest of the predicted time step ratios of each order.
      this%rdelt = Max (this%rdeltup, rdeltsm, this%rdeltdn)


   end subroutine estimateTimeStepRatio

!-----------------------------------------------------------------------------
!
! ROUTINE
!   testAccumulatedError
! DESCRIPTION
!     Test the accumulated error from the convergence process.
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine testAccumulatedError (this, ktloop, dely)

      ! ----------------------
      ! Argument declarations.
      ! ----------------------
      type (Manager_type) :: this
      integer, intent(in) :: ktloop
      real*8, intent(inout)  :: dely  (KBLOOP)

      integer:: jspc, kloop
      real*8  :: errymax

      do kloop = 1, ktloop
         dely(kloop) = 0.0d0
      end do

      do kloop = 1, ktloop
         do jspc = 1, this%num1stOEqnsSolve
            errymax     = this%dtlos(kloop,jspc) * this%chold(kloop,jspc)
            dely(kloop) = dely(kloop) + errymax * errymax
          end do
      end do

      this%der2max = 0.0d0

      do kloop = 1, ktloop
         if (dely(kloop) > this%der2max) then
            this%der2max = dely(kloop)
         end if
      end do

   end subroutine testAccumulatedError

!-----------------------------------------------------------------------------
! ROUTINE
!   calculateNewRmsError
! DESCRIPTION
!     Set the previous rms error and calculate the new rms error.
!     If dcon < 1, then sufficient convergence has occurred.  Otherwise,
!     if the ratio of the current to previous rmsError is decreasing,
!     iterate more.  If it is not, then the convergence test failed.
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine calculateNewRmsError (this, ktloop, dely, l3, savedVars)
      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved
      ! ----------------------
      ! Argument declarations.
      ! ----------------------
      type (Manager_type) :: this
      integer, intent(in) :: ktloop
      real*8, intent(in)  :: dely  (KBLOOP)
      integer, intent(inout) :: l3
      type(t_Smv2Saved), intent(inOut) :: savedVars


      real*8  :: rmsErrorPrevious, rmsrat
      integer :: kloop

      rmsErrorPrevious = this%rmsError
      this%der2max = 0.0d0

      ! make this a one line using maxval
      do kloop = 1, ktloop
        if (dely(kloop) > this%der2max) then
          this%der2max = dely(kloop)
        end if
      end do

      this%rmsError = Sqrt (this%der2max * this%order_inv)
      l3 = l3 + 1

      if (l3 > 1) then
        rmsrat = this%rmsError / rmsErrorPrevious
        this%drate  = Max (0.2d0*this%drate, rmsrat)
      else
        rmsrat = 1.0d0
      end if

      this%dcon = this%rmsError * Min (savedVars%conpst(this%nqq), savedVars%conp15(this%nqq)*this%drate)

   end subroutine calculateNewRmsError

!-----------------------------------------------------------------------------
! ROUTINE
!   tightenErrorTolerance
! DESCRIPTION
!     tighten absoloute error tolerance and restart integration
!     at beginning of time interva
! Created by: Megan Rose Damon
!   pr_smv2  : should the SmvgearII     output file be written
!   lunsmv   : logical unit number to write to when pr_smv2 is true
!   ncs      : identifies gas chemistry type (1..NCSGAS)
!   delt      : current time step (s)
! WARNING: this routine may not have been adequately tested
!-----------------------------------------------------------------------------
   subroutine tightenErrorTolerance (this, pr_smv2, lunsmv, ncs, delt, savedVars)

      use GmiPrintError_mod, only : GmiPrintError

      ! ----------------------
      ! Argument declarations.
      ! ----------------------
      type (Manager_type) :: this

      logical, intent(in)  :: pr_smv2
      integer, intent(in)  :: lunsmv
      integer, intent(in)  :: ncs ! ncs is argument to Smvgear
      real*8, intent(inout) :: delt
      type(t_Smv2Saved), intent(inOut) :: savedVars

      if (pr_smv2) then
         Write (lunsmv,950) delt, this%timeremain, this%failureFraction, savedVars%errmax(ncs)
      end if

      950    format ('Smvgear:  delt      = ', 1pe9.3, /,  '          timremain = ', 1pe9.3, /,  &
         &          '          failureFraction      = ', 1pe9.3, /,  '          errmax    = ', 1pe9.3)

      this%numErrTolDecreases = this%numErrTolDecreases + 1
      this%failureFraction     = this%failureFraction * 0.01d0

      ! handle magic number
      if (this%numErrTolDecreases == 10) then
         if (pr_smv2) then
            Write (lunsmv,960)
         end if

         960      format ('Smvgear:  too many decreases of failureFraction.')

         call GmiPrintError ('Problem in Smvgear', .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

   end subroutine tightenErrorTolerance

!-----------------------------------------------------------------------------
! ROUTINE
!   calculateTimeStep
! DESCRIPTION
!     Limit size of rdelt, then recalculate new time step and update
!     hratio.  Use hratio to determine whether or not the predictor
!     should be updated. (edit by MRD on 2/27/2013)
! Created by: Megan Rose Damon
!     delt      : current time step (s)
! Unit testing ideas: can't take a step bigger than what we have remaining
! Make this a method on a class?
! two routines: update time step and determine Jacobian (something like this)
!-----------------------------------------------------------------------------
   subroutine calculateTimeStep (this, delt, jeval, maxRelChange)

      ! ----------------------
      ! Argument declarations.
      ! ----------------------
      type (Manager_type) :: this
      real*8, intent(inout) :: delt
      integer, intent(out) :: jeval
      real*8, intent(in) :: maxRelChange

      real*8  :: hmtim

      hmtim  = Min (this%maxTimeStep, this%timeremain)
      this%rdelt  = Min (this%rdelt, this%rdelmax, hmtim/delt)
      delt   = delt   * this%rdelt

      this%hratio = this%hratio * this%rdelt
      this%xelaps = this%xelaps + delt
      ! rename nslp
      if ((Abs (this%hratio-1.0d0) > maxRelChange) .or. (this%numSuccessTdt >= this%nslp)) then
        jeval = 1 ! MRD: could be a boolean; this is signifying to whether or not to update Jacobian
      end if

   end subroutine calculateTimeStep


!-----------------------------------------------------------------------------
!
! ROUTINE
!   updateCoefficients
!
! DESCRIPTION
!     Update coefficients of (for?) the order; note that pertst2 is the original
!     pertst^2.
! MRD: Gear can be 1st, 2nd, etc. order
! MRD: this could be where it changes it order
! MRD: track these down, figure out what they are
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine updateCoefficients (this, savedVars)

      ! ----------------------
      ! Argument declarations.
      ! ----------------------
      type (Manager_type) :: this
      type(t_Smv2Saved), intent(inOut) :: savedVars

      real*8 :: eup ! pertst^2*order for one order higher than current order
      real*8  :: edwn ! pertst^2*order for one order lower  than current order

      this%nqqold = this%nqq
      this%kstep  = this%nqq + 1
      this%hratio = this%hratio * savedVars%aset(this%nqq,1) / this%asn1
      this%asn1   = savedVars%aset(this%nqq,1)
      this%enqq   = savedVars%pertst2(this%nqq,1) * this%order
      eup    = savedVars%pertst2(this%nqq,2) * this%order
      edwn   = savedVars%pertst2(this%nqq,3) * this%order
      this%conp3  = 1.4d0 /  (eup**savedVars%enqq3(this%nqq)) !eup is zero
      this%conp2  = 1.2d0 / (this%enqq**savedVars%enqq2(this%nqq)) !enqq is zero
      this%conp1  = 1.3d0 / (edwn**savedVars%enqq1(this%nqq)) !edwn is zero
      this%nqqisc = this%nqq * this%num1stOEqnsSolve

   end subroutine updateCoefficients

!-----------------------------------------------------------------------------
!
! ROUTINE
!   calcInitialTimeStepSize
!
! DESCRIPTION
!     Calculate initial time step size (s).
!     Sqrt (dely / [initialError * order]) =
!       rmsnorm of error scaled to initialError *
!       cnew + abtol / reltol
! This is a guess, later it will adapt
! Created by: Megan Rose Damon
!   ktloop   : # of grid-cells in a grid-block
!   dely     : TBD
!   delt      : current time step (s)
!   ncs      : identifies gas chemistry type (1..NCSGAS)
!-----------------------------------------------------------------------------
   subroutine calcInitialTimeStepSize (this, ktloop, dely, delt, ncs, savedVars)

      ! ----------------------
      ! Argument declarations.
      ! ----------------------
      type (Manager_type) :: this
      integer, intent(in)  :: ktloop
      real*8, intent(in)  :: dely  (KBLOOP)
      real*8, intent(out)    :: delt
      integer, intent(in)  :: ncs ! ncs is argument to Smvgear
      type(t_Smv2Saved), intent(inOut) :: savedVars

      integer :: kloop
      real*8  :: rmstop
      real*8  :: delt1

      rmstop = 0.0d0
      do kloop = 1, ktloop
         if (dely(kloop) > rmstop) rmstop = dely(kloop)
      end do

      delt1 = Sqrt (this%initialError / (savedVars%abst2(ncs) + (rmstop * this%order_inv)))
      delt  = Max  (Min (delt1, this%timeremain, this%maxTimeStep), HMIN)

   end subroutine calcInitialTimeStepSize

!-----------------------------------------------------------------------------
!
! ROUTINE
!   calculateErrorTolerances
!
! DESCRIPTION
! Use lowest absolute error tolerance when reordering.
! It is conceviable that there could be different norms or error criteria
! Created by: Megan Rose Damon
!   ktloop   : # of grid-cells in a grid-block
!   jlooplo  : low ntloop grid-cell - 1 in a grid-block
!   itloop   : # of zones (ilong * ilat * ivert)
!   cnew     : stores conc (y (estimated)) (molec/cm^3)
!   gloss    : value of first derivatives on output from velocity; right-side
!              of eqn on input to Backsub; error term (solution from Backsub)
!              on output from Backsub
!   dely     : TBD
!   errmx2   : measure of stiffness/nearness to convergence of each block
!              sum ydot/y for all species
!-----------------------------------------------------------------------------
   subroutine calculateErrorTolerances (this, ktloop, jlooplo, itloop, cnew, gloss, dely, errmx2)

      ! ----------------------
      ! Argument declarations.
      ! ----------------------
      type (Manager_type) :: this
      integer, intent(in)  :: ktloop
      integer, intent(in)  :: jlooplo
      integer, intent(in)  :: itloop
      real*8, intent(in) :: cnew  (KBLOOP, MXGSAER)
      real*8, intent(in) :: gloss (KBLOOP, MXGSAER)
      real*8, intent(inout)  :: dely  (KBLOOP)
      real*8, intent(inout) :: errmx2(itloop)

      integer :: kloop
      integer :: jspc
      real*8  :: errymax

      ! abtoler1 = failureFraction * abtol(6,ncs) / Min (errmax, 1.0d-03)
      do kloop = 1, ktloop
         do jspc = 1, this%num1stOEqnsSolve
            errymax     = gloss(kloop,jspc) / (cnew(kloop,jspc) + this%abtoler1)
            dely(kloop) = dely(kloop) + (errymax * errymax)
         end do
       end do

       do kloop = 1, ktloop
         errmx2(jlooplo+kloop) = dely(kloop)
       end do

   end subroutine calculateErrorTolerances

!-----------------------------------------------------------------------------
!
! ROUTINE
!   startTimeInterval
!
! DESCRIPTION
!   This routine starts time interval or is called after total failure.
! Created by: Megan Rose Damon
!   ncs      : identifies gas chemistry type (1..NCSGAS)
!-----------------------------------------------------------------------------
   subroutine startTimeInterval (this, ncs, savedVars)

      ! ----------------------
      ! Argument declarations.
      ! ----------------------
      type (Manager_type) :: this
      integer, intent(in)  :: ncs ! ncs is argument to Smvgear
      type(t_Smv2Saved), intent(inOut) :: savedVars

      this%idoub     = 2
      this%nslp      = MBETWEEN
      this%jrestar   = 0
      this%xelaps    = 0.0d0
      this%told      = 0.0d0
      this%timeremain = this%chemTimeInterval

      this%reltol1   = this%failureFraction * this%initialError_inv

      this%reltol2   = this%failureFraction * this%errmax_ncs_inv

      this%reltol3   = this%errmax_ncs_inv

      ! MRD: abtol in common block?
      ! put abtol in another structure
      this%abtoler1  = savedVars%abtol(6,ncs) * this%reltol1


   end subroutine startTimeInterval
!-----------------------------------------------------------------------------
!
! ROUTINE
!   resetGear
!
! DESCRIPTION
!   This routine is called once at the beginning of Smvgear
!
! Created by: Megan Rose Damon
!   ncs      : identifies gas chemistry type (1..NCSGAS)
!-----------------------------------------------------------------------------
      subroutine resetGear (this, ncsp, ncs, ifsun, hmaxnit, savedVars)

         implicit none

         ! ----------------------
         ! Argument declarations.
         ! ----------------------
         type (Manager_type) :: this
         integer, intent(out) :: ncsp
         integer, intent(in)  :: ncs ! ncs is argument to Smvgear
         integer, intent(in)  :: ifsun ! ifsun is an argument to Smvgear
         real*8,  intent(in)  :: hmaxnit
         type(t_Smv2Saved), intent(inOut) :: savedVars

         this%numFailOldJacobian     = 0
         this%jFail     = 0
         this%numFailErrorTest     = 0
         this%numFailAfterPredict     = 0
         this%numCallsPredict   = 0
         this%numSuccessTdt    = 0
         this%numCallsVelocity   = 0
         this%numErrTolDecreases  = 0

         ! MAX_REL_CHANGE moved to parameter
         this%rmsError    = 1.0d0

         ! MRD: derived from common block
         this%num1stOEqnsSolve    = savedVars%ischang(ncs)
         this%order     = this%num1stOEqnsSolve
         this%order_inv = 1.0d0 / this%num1stOEqnsSolve

         ! MRD: derived from common bloc
         this%chemTimeInterval = savedVars%timeintv(ncs)

         ! MRD: ICS is from common block
         ncsp      = (ifsun - 1) * ICS + ncs

         this%maxTimeStep = hmaxnit
         if (ifsun == 1) this%maxTimeStep = savedVars%hmaxday(ncs)

         this%failureFraction   = 1.0d0
         this%iabove = this%order * 0.4d0
         this%initialError     = Min (savedVars%errmax(ncs), 1.0d-03)
         this%initialError_inv = 1.0d0 / this%initialError
         this%errmax_ncs_inv = 1.0d0 / savedVars%errmax(ncs)

      end subroutine resetGear


end module GmiManager_mod
