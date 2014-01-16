module GmiMechanism_mod

   implicit none
   private

#     include "smv2chem_par.h"

   public :: velocity
   public :: calculateTermOfJacobian
   public :: Mechanism_type


! MRD: add type bound procedures here
! can remove "_type"
! need a constructor
! there will be a generic mechanism and generic matrix module
! but right now we have a sparse matrix linear system (not general)
! this is a GMI-specific verison of sparse matrix
   type Mechanism_type

       integer :: numRxns1, numRxns2, numRxns3 ! number of rxns
                     ! with 1, 2, and 3 reactants
       integer :: numRxns3Drep !numRxns3 + # of rxns with
                     ! two active reactants that are not
                       ! followed by a rxn with the same reactant
       integer :: numActiveReactants ! CHECK WITH TOM
       integer :: speciesNumberA    (NMTRATE)
       integer :: speciesNumberB    (NMTRATE)
       integer :: speciesNumberC    (NMTRATE)
       integer :: numGridCellsInBlock
       real*8  :: rateConstants (KBLOOP, NMTRATE) ! rate coefficient:
                              ! rates with 1 reactant:   s^-1
                              ! rates with 2 reactants:  l-h2o mole^-1 s^-1 or
                              !                         cm^3 #-1 s^-1 (?)
                              ! rates with 3 reactants:  l^2-h2o m-2 s-1  or
                              !                         cm^6 #-2 s-1 (?)
    end type Mechanism_type

contains

!-----------------------------------------------------------------------------
!
! ROUTINE
!   velocity
!
! DESCRIPTION
!   This routine evaluates the first derivative of each ordinary
!   differential equation (ODE).  It evaluates derivatives in the special
!   form f = y'(est) = f(x,y,estimated), where f is the right hand side of
!   the differential equation.
!
!   Example =>
!
!     Species:         A,   B,   C
!     Concentrations: (A), (B), (C)
!
!     Reactions:    1) A          --> B      J
!                   2) A  + B     --> C      K1
!                   3) A  + B + C --> D      K2
!
!     First         d(A) / dt = -J(A) - K1(A)(B) - K2(A)(B)(C)
!     Derivatives:  d(B) / dt = +J(A) - K1(A)(B) - K2(A)(B)(C)
!                   d(C) / dt =       + K1(A)(B) - K2(A)(B)(C)
!                   d(D) / dt =                  + K2(A)(B)(C)
!
! ARGUMENTS
!   ischan   : # of first-order eqns to solve, = # of spc = order of original
!              matrix; ischan has a different value for day and night, and for
!              gas- and aqueous-phase chemistry;
!              # spc with prod or loss terms in Smvgear (?)
!   ncsp     : ncs       => for daytime   gas chemistry
!              ncs + ICS => for nighttime gas chemistry
!   concentrationsNew     : init (and final) spc conc (# cm^-3-air or moles l^-1-h2o (?))
!   gloss    : first derivative = sum of prod. minus loss rates for a spc
!   reactionRates    : rxn rate (moles l^-1-h2o s^-1 or # cm^-3 s^-1 (?))
!   nfdh1    : nfdh2 + # of rxns with one   active reactant
!
!-----------------------------------------------------------------------------
! MRD: Subfun, which will part of the mechanism, needs to know sparseMatrix things.
! MRD: Subfun is similiar to Flow, in SmvgearTDD (or y dot, or velocity)

      subroutine velocity (this, ischan, ncsp, concentrationsNew, &
                        gloss, reactionRates, nfdh1, savedVars)

      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

      implicit none

#     include "smv2chem_par.h"
!#     include "smv2chem2.h"

!     ----------------------
!     Argument declarations.
!     ----------------------
      type (Mechanism_type) :: this
      integer, intent(in)  :: ischan ! derived from common block
      integer, intent(in)  :: ncsp   ! derived from common block
      real*8,  intent(in)  :: concentrationsNew (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: gloss(KBLOOP, MXGSAER)
      real*8,  intent(inout) :: reactionRates(KBLOOP, NMTRATE*2) ! production and loss are seperate reactions
      integer, intent(out) :: nfdh1
      type(t_Smv2Saved), intent(inOut) :: savedVars


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer ::  ja, jb, jc
      integer ::  jspc
      integer ::  k
      integer ::  n, nc, nh
      integer ::  nh1, nh2, nh3, nh4, nh5
      integer ::  nk0, nk1, nk2
      integer ::  nk3, nk4, nkn
      integer ::  nl1, nl2, nl3, nl4, nl5
      integer ::  npl

!     -----------------------------------------------------------------------
!     concmult : product of concs in a rate; if two consecutive rxns have the
!                same spc reacting (e.g., A + B --> C and A + B --> D + E),
!                then use the same value for both rxns
!     -----------------------------------------------------------------------

      real*8  :: concmult
      real*8  :: fracn

!     ----------------
!     Begin execution.
!     ----------------

!     Write (6,*) 'velocity called.'

!     ----------------------
!     Set rates of reaction.
!     ----------------------
      nfdh1   = this%numRxns2 + savedVars%ioner(ncsp)

!     ---------------------------------------------------------
!     First derivatives for rates with three active loss terms.
!     ---------------------------------------------------------
      do nkn = 1, this%numRxns3

        ja = this%speciesNumberA(nkn)
        jb = this%speciesNumberB(nkn)
        jc = this%speciesNumberC(nkn)

        nh = nkn + this%numActiveReactants

        do k = 1, this%numGridCellsInBlock
          reactionRates(k,nkn) =  this%rateConstants(k,nkn) * &
               & concentrationsNew(k,ja) * concentrationsNew(k,jb) * &
               & concentrationsNew(k,jc)
          reactionRates(k,nh)  = -reactionRates(k,nkn) !MRD: species d goes up
        end do

      end do

!     -------------------------------------------------------
!     First derivatives for rates with two active loss terms.
!     -------------------------------------------------------
      do nkn = this%numRxns3 + 1, this%numRxns3Drep

        ja = this%speciesNumberA(nkn)
        jb = this%speciesNumberB(nkn)

        nh = nkn + this%numActiveReactants

        do k = 1, this%numGridCellsInBlock
          reactionRates(k,nkn) = this%rateConstants(k,nkn) * &
               & concentrationsNew(k,ja) * concentrationsNew(k,jb)
          reactionRates(k,nh)  = -reactionRates(k,nkn)
        end do

      end do

!     -----------------------------------------------------------
!     First derivatives for rates with two active loss terms and
!     where the subsequent reaction has the same reactants, but a
!     different rate.
!     -----------------------------------------------------------

      do nkn = this%numRxns3Drep + 1, this%numRxns2, 2

        ja  = this%speciesNumberA(nkn)
        jb  = this%speciesNumberB(nkn)

        nk2 = nkn + 1
        nh  = nkn + this%numActiveReactants
        nh2 = nk2 + this%numActiveReactants

        do k = 1, this%numGridCellsInBlock
          concmult  =  concentrationsNew(k,ja) * concentrationsNew(k,jb)
          reactionRates(k,nkn) =  this%rateConstants(k,nkn) * concmult
          reactionRates(k,nk2) =  this%rateConstants(k,nk2) * concmult
          reactionRates(k,nh)  = -reactionRates(k,nkn)
          reactionRates(k,nh2) = -reactionRates(k,nk2)
        end do

      end do

!     ------------------------------------------------------
!     First derivatives for rates with one active loss term.
!     ------------------------------------------------------
      do nkn = this%numRxns2 + 1, nfdh1

        ja = this%speciesNumberA(nkn)

        nh = nkn + this%numActiveReactants

        do k = 1, this%numGridCellsInBlock
          reactionRates(k,nkn) =  this%rateConstants(k,nkn) * concentrationsNew(k,ja)
          reactionRates(k,nh)  = -reactionRates(k,nkn)
        end do

      end do

!     --------------------------------
!     Initialize first derivative = 0.
!     --------------------------------
      gloss(1:this%numGridCellsInBlock,1:ischan) = 0.0d0

!     ---------------------------------------------------------------
!     Sum net (not reproduced) kinetic and photo gains and losses for
!     each species.
!
!     Sum 1,2,3,4, or 5 terms at a time to improve vectorization.
!     ---------------------------------------------------------------
!     ==========================================
      NPLLOOP: do npl = savedVars%npllo(ncsp), savedVars%nplhi(ncsp)
!     ==========================================
        ! MRD: these arrays should be part of sparseMatrix
        jspc = savedVars%jspnpl(npl)

        nl5 = savedVars%npl5(npl)
        nh5 = savedVars%nph5(npl)
        nl4 = savedVars%npl4(npl)
        nh4 = savedVars%nph4(npl)
        nl3 = savedVars%npl3(npl)
        nh3 = savedVars%nph3(npl)
        nl2 = savedVars%npl2(npl)
        nh2 = savedVars%nph2(npl)
        nl1 = savedVars%npl1(npl)
        nh1 = savedVars%nph1(npl)

!       -- Sum 5 terms at a time. --

        do nc = nl5, nh5

          ! MRD: this is an out dated (unnecessary?) optimization
          nk0 = savedVars%lossra(nc)
          nk1 = savedVars%lossrb(nc)
          nk2 = savedVars%lossrc(nc)
          nk3 = savedVars%lossrd(nc)
          nk4 = savedVars%lossre(nc)

          ! up to five species
          do k = 1, this%numGridCellsInBlock
            gloss(k,jspc) =  &
     &        gloss(k,jspc) -  &
     &        reactionRates(k,nk0)  - reactionRates(k,nk1) - reactionRates(k,nk2) -  &
     &        reactionRates(k,nk3)  - reactionRates(k,nk4)
          end do

        end do

!       -- Sum 4 terms at a time. --

        do nc = nl4, nh4

          ! MRD: this is an out dated (unnecessary?) optimization
          nk0 = savedVars%lossra(nc)
          nk1 = savedVars%lossrb(nc)
          nk2 = savedVars%lossrc(nc)
          nk3 = savedVars%lossrd(nc)

          do k = 1, this%numGridCellsInBlock
            gloss(k,jspc) =  &
     &        gloss(k,jspc) -  &
     &        reactionRates(k,nk0)  - reactionRates(k,nk1) - reactionRates(k,nk2) -  &
     &        reactionRates(k,nk3)
          end do

        end do

!       -- Sum 3 terms at a time. --

        do nc = nl3, nh3

          nk0 = savedVars%lossra(nc)
          nk1 = savedVars%lossrb(nc)
          nk2 = savedVars%lossrc(nc)

          do k = 1, this%numGridCellsInBlock
            gloss(k,jspc) =  &
     &        gloss(k,jspc) -  &
     &        reactionRates(k,nk0)  - reactionRates(k,nk1)  - reactionRates(k,nk2)
          end do

        end do

!       -- Sum 2 terms at a time. --

        do nc = nl2, nh2

          nk0 = savedVars%lossra(nc)
          nk1 = savedVars%lossrb(nc)

          do k = 1, this%numGridCellsInBlock
            gloss(k,jspc) =  &
     &        gloss(k,jspc) -  &
     &        reactionRates(k,nk0)  - reactionRates(k,nk1)
          end do

        end do

!       -- Sum 1 term at a time. --

        do nc = nl1, nh1

          nk0 = savedVars%lossra(nc)

          do k = 1, this%numGridCellsInBlock
            gloss(k,jspc) =  &
     &        gloss(k,jspc) -  &
     &        reactionRates(k,nk0)
          end do

        end do

!     ==============
      end do NPLLOOP
!     ==============


!     --------------------------------------------------------------
!     Sum production term for reactions where products fractionated.
!     --------------------------------------------------------------

      do n = savedVars%nfrlo(ncsp), savedVars%nfrhi(ncsp)

        jspc  = savedVars%jspcnfr(n)
        nkn   = savedVars%nknfr  (n)
        fracn = savedVars%fracnfr(n)

        do k = 1, this%numGridCellsInBlock
          gloss(k,jspc) = gloss(k,jspc) + (fracn * reactionRates(k,nkn))
        end do

      end do


      return

      end subroutine velocity




!-----------------------------------------------------------------------------
!
! ROUTINE
!   calculateJacobian
!
! DESCRIPTION
!    This is calculating the Jacobian of the flow (jacobian) or
!    Partial derivatives for rates with three active loss terms.
!
!
! ARGUMENTS
!   numRxns1           : numRxns2 + # of rxns with one   active reactant
!   numRxns2           : numRxns3 + # of rxns with two   active reactants
!   numRxns3           :         # of rxns with three active reactants
!   speciesNumberA,B,C : spc # of each reactant; locates reordered active spc #s
!   numGridCellsInBlock   : # of grid-cells in a grid-block
!   rateConstants      : rate constants
!   concentrationsNew  : stores conc (y (estimated)) (molec/cm^3)
!   reactionRates :
!
! NOTES:
    ! Should the Jacobian change within time step?
    ! Chemical time step is small, and predictor takes multiple steps.
    ! rateConstants (reaction rates stay the same throughout the simulation)
    ! if chem_tdt changes or if the order changes, the Jacobian will change
!-----------------------------------------------------------------------------

      subroutine calculateTermOfJacobian (this, concentrationsNew, reactionRates)

         implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------
         type (Mechanism_type) :: this
         real*8,  intent(in)  :: concentrationsNew  (KBLOOP, MXGSAER)
         real*8,  intent(out) :: reactionRates (KBLOOP, NMTRATE, 3)

!     ----------------------
!     Variable declarations.
!     ----------------------
        integer :: nkn, block
        integer :: ja, jb, jc

!    Partial derivatives for rates with three active loss terms.
         do nkn = 1, this%numRxns3
            ja = this%speciesNumberA(nkn)
            jb = this%speciesNumberB(nkn)
            jc = this%speciesNumberC(nkn)

            do block = 1, this%numGridCellsInBlock
               reactionRates(block,nkn,1) = this%rateConstants(block,nkn) * concentrationsNew(block,jb) * concentrationsNew(block,jc)
               reactionRates(block,nkn,2) = this%rateConstants(block,nkn) * concentrationsNew(block,ja) * concentrationsNew(block,jc)
               reactionRates(block,nkn,3) = this%rateConstants(block,nkn) * concentrationsNew(block,ja) * concentrationsNew(block,jb)
            end do

          end do

!     ---------------------------------------------------------
!     Partial derivatives for rates with two active loss terms.
!     ---------------------------------------------------------
         do nkn = this%numRxns3+1, this%numRxns2

            ja = this%speciesNumberA(nkn)
            jb = this%speciesNumberB(nkn)

            do block = 1, this%numGridCellsInBlock
               reactionRates(block,nkn,1) = this%rateConstants(block,nkn) * concentrationsNew(block,jb)
               reactionRates(block,nkn,2) = this%rateConstants(block,nkn) * concentrationsNew(block,ja)
            end do

         end do

!     --------------------------------------------------------
!     Partial derivatives for rates with one active loss term.
!     --------------------------------------------------------
         do nkn = this%numRxns2+1, this%numRxns1
            do block = 1, this%numGridCellsInBlock
               reactionRates(block,nkn,1) = this%rateConstants(block,nkn)
            end do
         end do

      end subroutine calculateTermOfJacobian

   end module GmiMechanism_mod

