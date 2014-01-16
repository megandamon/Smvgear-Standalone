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
      subroutine calculatePredictor (nondiag, iarry, numGridCellsInBlock, &
      & npdh, npdl, r1delt, jacobian, predictor, savedVars)

      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------
         integer, intent(in) :: nondiag     ! # of final matrix positions, excluding diagonal
                             ! terms, filled after all matrix processes
         integer, intent(in) :: iarry
         integer, intent(in)  :: numGridCellsInBlock
         integer, intent(in) :: npdh, npdl

         real*8,  intent(in)  :: r1delt
         real*8,  intent(in) :: jacobian(KBLOOP, NMTRATE, 3)
         real*8,  intent(inout) :: predictor  (KBLOOP, 0:MXARRAY)
         type(t_Smv2Saved), intent(inOut) :: savedVars


!     ----------------------
!     Variable declarations.
!     ----------------------
         integer :: iar, k, n, nkn, ial
         real*8  :: fracr1


         ! list of non-zero values
         ! MRD: derived type that combines the predictor below,
         ! with this information below (next two loops)
         ! could be called sparseMatrix (stay on the solver type)
         ! is predictor and Jacobian the same size?
         do iar = 1, nondiag
            do k = 1, numGridCellsInBlock
               predictor(k,iar) = 0.0d0
            end do
         end do

         do iar = nondiag + 1, iarry
            do k = 1, numGridCellsInBlock
               predictor(k,iar) = 1.0d0
            end do
         end do

         do n = npdl, npdh

            nkn    = savedVars%nkpdterm(n) ! nkpdterm in common block
            iar    = savedVars%ipospd  (n) ! ipospd in common block
            ial    = savedVars%iialpd  (n) ! iialpd in common block
            fracr1 = savedVars%fracpl  (n) * r1delt ! fracpl in a common block

            do k = 1, numGridCellsInBlock
               ! MRD: iar and ial need to get attached to sparseMatrix type
               predictor(k,iar) = predictor(k,iar) + (fracr1 * jacobian(k,nkn,ial))
            end do

         end do

      end subroutine calculatePredictor
end module GmiSparseMatrix_mod
