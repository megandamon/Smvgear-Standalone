module GmiMechanism_mod


   use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

   implicit none
   private

#     include "smv2chem_par.h"

   public :: initializeMechanism
   public :: setBoundaryConditions
   public :: velocity
   public :: calculateTermOfJacobian
   public :: updatePhotoDissRates
   public :: Mechanism_type


! MRD: add type bound procedures here
! can remove "_type"
! need a constructor
! there will be a generic mechanism and generic matrix module
! but right now we have a sparse matrix linear system (not general)
! this is a GMI-specific verison of sparse matrix
   type Mechanism_type

       integer :: numRxns3Drep !numRxns3 + # of rxns with
                     ! two active reactants that are not
                       ! followed by a rxn with the same reactant
       !integer :: numActiveReactants ! CHECK WITH TOM
       !K: These arrays should be removed as they're in ChemTable now
       integer :: speciesNumberA    (NMTRATE)
       integer :: speciesNumberB    (NMTRATE)
       integer :: speciesNumberC    (NMTRATE)
       integer :: numRxns1, numRxns2, numRxns3 ! number of rxns
                     ! with 1, 2, and 3 reactants
       integer :: numGridCellsInBlock ! probably all the same (check padding)

       ! MRD
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
!   initializeMechansim
!
! DESCRIPTION
!-----------------------------------------------------------------------------
      subroutine initializeMechanism (this, ktloop, irma, &
                                      & irmb, irmc, nfdh2, nfdh3, nfdrep, rrate)
         !     ----------------------
         !     Argument declarations.
         !     ----------------------
         type (Mechanism_type) :: this
         integer, intent(in)  :: irma    (NMTRATE)
         integer, intent(in)  :: irmb    (NMTRATE)
         integer, intent(in)  :: irmc    (NMTRATE)
         integer, intent(in)  :: ktloop
         integer, intent(in)  :: nfdh2,  nfdh3
         integer, intent(in)  :: nfdrep
         real*8,  intent(in) :: rrate (KBLOOP, NMTRATE)

         this%numGridCellsInBlock = ktloop
         this%speciesNumberA = irma ! MRD: these probably shouldn't be in object
         this%speciesNumberB = irmb
         this%speciesNumberC = irmc
         this%numRxns2 = nfdh2
         this%numRxns3 = nfdh3
         this%numRxns3Drep = nfdrep
         this%rateConstants = rrate
      end subroutine initializeMechanism

!-----------------------------------------------------------------------------
!
! ROUTINE
!   setBoundaryConditions
!
! DESCRIPTION
!     Zero first derviatives in surface zones for species with fixed
!     concentration boundary conditions (LLNL addition, PSC, 5/14/99).
!
! MRD: this will probably be called within mechanism and not directly
! ARGUMENTS
!   jreorder : gives original grid-cell from re-ordered grid-cell
!   jlooplo  : low ntloop grid-cell - 1 in a grid-block
!   ilat     : # of latitudes
!   ilong    : # of longitudes
!   ntspec   : # of active + inactive gases
!   ncs      : identifies gas chemistry type (1..NCSGAS)
!   inewold  : original spc # of each new jnew spc
!   do_semiss_inchem : do surface emissions inside the chemistry solver, or
!                      outside?
!   gloss    : value of first derivatives on output from velocity; right-side
!              of eqn on input to Backsub; error term (solution from Backsub)
!              on output from Backsub
!   yemis    : surface emissions (units?)
!-----------------------------------------------------------------------------
      subroutine setBoundaryConditions (this, itloop, jreorder, jlooplo, ilat, ilong, &
                  & ntspec, ncs, inewold, do_semiss_inchem, gloss, yemis, ibcb)


         implicit none

         !     ----------------------
         !     Argument declarations.
         !     ----------------------
         type (Mechanism_type) :: this
         integer, intent(in)  :: itloop
         integer, intent(in)  :: jreorder(itloop)
         integer, intent(in)  :: jlooplo
         integer, intent(in)  :: ilat, ilong
         integer, intent(in)  :: ntspec  (ICS)
         integer, intent(in)  :: ncs
         integer, intent(in)  :: inewold (MXGSAER, ICS)
         logical, intent(in)  :: do_semiss_inchem
         real*8,  intent(inout) :: gloss (KBLOOP, MXGSAER)
         real*8,  intent(in)  :: yemis   (ilat*ilong, IGAS)
         integer :: ibcb(IGAS)

         integer :: kloop, jspc, jgas


         !Write (6,*) 'setBoundaryConditions called'

         do kloop = 1, this%numGridCellsInBlock
            if (jreorder(jlooplo+kloop) <= (ilat*ilong)) then
               do jspc = 1, ntspec(ncs)
                  jgas = inewold(jspc,1)
                  if (ibcb(jgas) == 1) then
                     gloss(kloop,jspc) = 0.0d0
                  else if (do_semiss_inchem) then
                     gloss(kloop,jspc) =  gloss(kloop,jspc) +  &
                     &          yemis(jreorder(jlooplo+kloop),jgas)
                  end if
               end do
            end if
         end do

      end subroutine setBoundaryConditions

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

      subroutine velocity (this, ischan, ncsp, cx, &
                        gloss, nfdh1, savedVars)

      use ChemTable_mod
      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved
      implicit none

#     include "smv2chem_par.h"

      type (Mechanism_type) :: this
      integer, intent(in)  :: ischan ! derived from common block
      integer, intent(in)  :: ncsp   ! derived from common block
      real*8,  intent(in)  :: cx (KBLOOP, MXGSAER)
      real*8,  intent(out) :: gloss(KBLOOP, MXGSAER)
      integer, intent(out) :: nfdh1
      type(t_Smv2Saved), intent(inOut) :: savedVars


      integer :: i,j,k,spcnum,nf,ktloop,nin,ina,inb,inc
      integer :: nout, outa,outb,outc,outd
      real*8 :: frac,Rx(KBLOOP)

      !Write(*,*) 'Velocity called'
      ktloop = this%numGridCellsInBlock

      nf = 0
      gloss = 0.0

      !K: Why does this need to be sent back out?
      nfdh1 = GenChem%nfdh1

      do i=1,Nrxn
         nin = GenChem%RxnNumIn(i)
         select case (nin)
            case(2)
               ina = GenChem%RxnIn(1,i)
               inb = GenChem%RxnIn(2,i)
               do k=1,ktloop
                  Rx(k) = this%rateConstants(k,i)*cx(k,ina)*cx(k,inb)
                  gloss(k,ina) = gloss(k,ina) - Rx(k)
                  gloss(k,inb) = gloss(k,inb) - Rx(k)
               end do
            case(1)
               ina = GenChem%RxnIn(1,i)
               do k=1,ktloop
                  Rx(k) = this%rateConstants(k,i)*cx(k,ina)
                  gloss(k,ina) = gloss(k,ina) - Rx(k)
               end do
            case(0)
               do k=1,ktloop
                  Rx(k) = this%rateConstants(k,i)
               end do
            case(3)
               ina = GenChem%RxnIn(1,i)
               inb = GenChem%RxnIn(2,i)
               inc = GenChem%RxnIn(3,i)
               do k=1,ktloop
                  Rx(k) = this%rateConstants(k,i)*cx(k,ina)*cx(k,inb)*cx(k,inc)
                  gloss(k,ina) = gloss(k,ina) - Rx(k)
                  gloss(k,inb) = gloss(k,inb) - Rx(k)
                  gloss(k,inc) = gloss(k,inc) - Rx(k)
               end do
         end select

         nout = GenChem%RxnNumOut(i)
         select case (nout)
            case(2)
               outa = GenChem%RxnOut(1,i)
               outb = GenChem%RxnOut(2,i)
               do k=1,ktloop
                  gloss(k,outa) = gloss(k,outa) + Rx(k)
                  gloss(k,outb) = gloss(k,outb) + Rx(k)
               end do
            case(1)
               outa = GenChem%RxnOut(1,i)
               do k=1,ktloop
                  gloss(k,outa) = gloss(k,outa) + Rx(k)
               end do

            case(3)
               outa = GenChem%RxnOut(1,i)
               outb = GenChem%RxnOut(2,i)
               outc = GenChem%RxnOut(3,i)
               do k=1,ktloop
                  gloss(k,outa) = gloss(k,outa) + Rx(k)
                  gloss(k,outb) = gloss(k,outb) + Rx(k)
                  gloss(k,outc) = gloss(k,outc) + Rx(k)
               end do

            case default !Do general case
               do j=1,nout
                  spcnum = GenChem%RxnOut(j,i)
                  do k=1,ktloop
                     gloss(k,spcnum) = gloss(k,spcnum) + Rx(k)
                  end do
               end do
         end select

         if (GenChem%isFracRxn(i)) then
            nf = nf+1
            do j=1,GenChem%FracRxnNumOut(nf)
               spcnum = GenChem%FracRxnOut(j,nf)
               frac = GenChem%FracRate(j,nf)
               do k=1,ktloop
                  gloss(k,spcnum) = gloss(k,spcnum) + frac*Rx(k)
               end do
            end do
         end if

      end do
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
               reactionRates(block,nkn,1) = this%rateConstants(block,nkn) &
                    & * concentrationsNew(block,jb) * concentrationsNew(block,jc)
               reactionRates(block,nkn,2) = this%rateConstants(block,nkn) &
                    & * concentrationsNew(block,ja) * concentrationsNew(block,jc)
               reactionRates(block,nkn,3) = this%rateConstants(block,nkn) &
                    & * concentrationsNew(block,ja) * concentrationsNew(block,jb)
            end do

          end do

!     ---------------------------------------------------------
!     Partial derivatives for rates with two active loss terms.
!     ---------------------------------------------------------
         do nkn = this%numRxns3+1, this%numRxns2

            ja = this%speciesNumberA(nkn)
            jb = this%speciesNumberB(nkn)

            do block = 1, this%numGridCellsInBlock
               reactionRates(block,nkn,1) = this%rateConstants(block,nkn) &
                    & * concentrationsNew(block,jb)
               reactionRates(block,nkn,2) = this%rateConstants(block,nkn) &
                    & * concentrationsNew(block,ja)
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

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Update
!
! DESCRIPTION
!   This routine updates photodissociation rates.
!
!   Photorates are included in first and partial derivative equations.
!
! ARGUMENTS
!   ktloop   : # of grid-cells in a grid-block
!   numActiveReactants    : # of active rxns
!   ncs      : identifies gas chemistry type (1..NCSGAS)
!   ncsp     : ncs       => for daytime   gas chemistry
!              ncs + ICS => for nighttime gas chemistry
!   jphotrat : tbd
!   ptratk1  : tbd
!   rrate    : rate constants
!
!-----------------------------------------------------------------------------
! MRD: Update is probably setPhotolysisCoeffs and computeRateCoeffs
! MRD: Is going into the mechanism
! MRD: Solver will not call it

      subroutine updatePhotoDissRates  &
     &  (this, ktloop, numActiveReactants, ncs, ncsp, &
     &  jphotrat, pratk1, savedVars)

      implicit none
#     include "smv2chem_par.h"


!     ----------------------
!     Argument declarations.
!     ----------------------
      type (Mechanism_type) :: this
      integer, intent(in)  :: ktloop
      integer, intent(in)  :: numActiveReactants
      integer, intent(in)  :: ncs
      integer, intent(in)  :: ncsp
      integer, intent(in)  :: jphotrat(ICS)
      real*8,  intent(in)  :: pratk1  (KBLOOP, IPHOT)
      type(t_Smv2Saved), intent(inOut) :: savedVars


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: i, j
      integer :: kloop
      integer :: nh, nk, nkn


!     ----------------
!     Begin execution.
!     ----------------

!c    Write (6,*) 'Update called.'


!     -------------------------------
!     Load photolysis rate constants.
!     -------------------------------

      do j = 1, jphotrat(ncs)

        nkn = savedVars%nknphotrt(j,ncs)

        do kloop = 1, ktloop
          this%rateConstants(kloop,nkn) = pratk1(kloop,j)
        end do

      end do

      return

      end subroutine updatePhotoDissRates

   end module GmiMechanism_mod



