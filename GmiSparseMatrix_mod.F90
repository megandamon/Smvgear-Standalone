module GmiSparseMatrix_mod

   implicit none
   private

#     include "smv2chem_par.h"
!#     include "smv2chem2.h"

   public :: calculatePredictor
   public :: SparseMatrix_type


   type SparseMatrix_type

       integer :: holder

    end type SparseMatrix_type

contains
!-----------------------------------------------------------------------------
!
! ROUTINE
!   calculatePredictor
!
! DESCRIPTION
!     Calculates the predictor matrix: (P) = I - h * b * J:
!       J = Jacobian matrix of partial derivates
!       I = identity matrix
!       h = time step
!       b = coefficient of method
!       R = h * b = -r1delt
!
! ARGUMENTS
!   nondiag    : # of final matrix positions, excluding diagonal
                             ! terms, filled after all matrix processes
!   nondiag1   : nondiag + 1
!   iarry      : iarray(ncsp)
!   numGridCellsInBlock : # of grid-cells in a grid-block
!   npdh  :
!   npdl  :
!   r1delt : = -aset(nqq,1) * time_step = -coefficient_of_method * dt
!   jacobian : term of Jacobian (J) = partial derivative
!   predictor : array of iarray units holding values of each matrix
!              position actually used;
!              cc2 = P = I - delt * aset(nqq,1) * partial_derivatives
!
! NOTES: MRD: solver should not know number of reactions
!-----------------------------------------------------------------------------
      ! MRD: I + bT
      subroutine calculatePredictor (nondiag, iarry, mechObject, cx, &
      & npdh, npdl, r1delt, Predict, savedVars)

      use ChemTable_mod
      use GmiMechanism_mod
      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------
         integer, intent(in) :: nondiag     ! # of final matrix positions, excluding diagonal
                             ! terms, filled after all matrix processes
         integer, intent(in) :: iarry
         type (Mechanism_type), intent(in) :: mechObject
         real*8, intent(in)  :: cx(KBLOOP,MXGSAER)
         integer, intent(in) :: npdh, npdl

         real*8,  intent(in)  :: r1delt
         real*8,  intent(inout) :: Predict  (KBLOOP, 0:MXARRAY)
         type(t_Smv2Saved), intent(inOut) :: savedVars


!     ----------------------
!     Variable declarations.
!     ----------------------
         integer :: iar, k, n, nkn, ial,ktloop, nin, ina, inb
         real*8  :: fracr1

         ktloop = mechObject%numGridCellsInBlock

         ! list of non-zero values
         ! MRD: derived type that combines the predictor below,
         ! with this information below (next two loops)
         ! could be called sparseMatrix (stay on the solver type)
         ! is predictor and Jacobian the same size?
         do iar = 1, nondiag
            do k = 1, ktloop
               Predict(k,iar) = 0.0d0
            end do
         end do

         do iar = nondiag + 1, iarry
            do k = 1, ktloop
               Predict(k,iar) = 1.0d0
            end do
         end do

         do n = npdl, npdh

            nkn    = savedVars%nkpdterm(n) ! nkpdterm in common block
            iar    = savedVars%ipospd  (n) ! ipospd in common block
            ial    = savedVars%iialpd  (n) ! iialpd in common block
            fracr1 = savedVars%fracpl  (n) * r1delt ! fracpl in a common block

            nin = GenChem%RxnNumIn(nkn)
            select case (nin)
               case(2)
                  if (ial .eq. 1) then
                     ina = GenChem%RxnIn(2,nkn)
                  else
                     ina = GenChem%RxnIn(1,nkn)
                  end if
                  do k=1, ktloop
                     Predict(k,iar) = Predict(k,iar) + fracr1* &
                     & mechObject%rateConstants(k,nkn)*cx(k,ina)
                  end do
               case(1) !Only one var to take deriv wrt
                  do k=1, ktloop
                     Predict(k,iar) = Predict(k,iar)+ (fracr1*mechObject%rateConstants(k,nkn))
                  end do
               case(3)
                  if (ial .eq. 1) then
                     ina = GenChem%RxnIn(2,nkn)
                     inb = GenChem%RxnIn(3,nkn)
                  else if (ial .eq. 2) then
                     ina = GenChem%RxnIn(1,nkn)
                     inb = GenChem%RxnIn(3,nkn)
                  else
                     ina = GenChem%RxnIn(1,nkn)
                     inb = GenChem%RxnIn(2,nkn)
                  endif
                  do k=1, ktloop
                     Predict(k,iar) = Predict(k,iar) + fracr1* &
                     & mechObject%rateConstants(k,nkn)*cx(k,ina)*cx(k,inb)
                  end do
            end select


         end do

      end subroutine calculatePredictor
end module GmiSparseMatrix_mod
