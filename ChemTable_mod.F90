module ChemTable_mod

   use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

   implicit none

      integer, parameter :: Nrxn = 402 !Total number reactions, kinetic+photolysis
      integer, parameter :: RxnMaxIn = 3 !Maximum number of input species for any reaction
      integer, parameter :: RxnMaxOut = 6 !Max number of output species (not incl. fractionated) for any reaction
      integer, parameter :: Nfrac = 69 !Number of fractionated reactions
      integer, parameter :: FracRxnMax = 11 !Max number of output species that are fractionated
   type ChemTable

      integer :: gCT !gasChemistryType
      integer :: nfdh3,nfdl2,nfdh2,nfdl1,nfdh1,nfdl0
      !Integer arrays, maybe change to short ints?
      integer :: RxnNumIn(Nrxn),RxnNumOut(Nrxn) !Number of in/out species for each reaction

      integer :: RxnIn(RxnMaxIn,Nrxn) !Specific species in/out for each reaction
      integer :: RxnOut(RxnMaxOut,Nrxn)

      !Fractionated reaction data
      integer :: FracRxns(Nfrac) !Rxn # of the ith fractionated reaction
      logical :: isFracRxn(Nrxn) !True if ith reaction is fractionated
      integer :: FracRxnNumOut(Nfrac) !Number of fractionated outputs for each fractionated reaction
      integer :: FracRxnOut(FracRxnMax,Nfrac) !species number of each fractionated output
      Real*8     :: FracRate(FracRxnMax,Nfrac) !fraction of reaction rate contributing to each output species

   end type ChemTable

   type (ChemTable) :: GenChem

   contains

   Subroutine CalcTabl (Tabl,gCT,ncs,ncsp, savedVars)

#include "smv2chem_par.h"

      type (ChemTable), intent(out) :: Tabl
      integer, intent(in) :: gCT, ncs, ncsp !gasChemistryType, ncs, ncsp
      type(t_Smv2Saved), intent(inOut) :: savedVars

      integer :: nfdh3,nfdl2,nfdh2,nfdl1,nfdh1,nfdl0,nallr
      integer :: nkn, nk, i,j,k, jspc, Sumin,Sumtemp
      !Temp arrays to invert reaction mapping
      integer :: SpcRxn(NNACT,NMTRATE),SpcRxnNum(NNACT),NumIntemp(NMTRATE),RxIntemp(RxnMaxIn,Nrxn)
      integer :: NumFracOuts(NMTRATE),FracOutSpc(FracRxnMax,NMTRATE)
      Real*8 :: FracRate(FracRxnMax,NMTRATE)
      SpcRxnNum = 0
      SpcRxn = 0
      NumIntemp = 0

      Tabl%gCT = gCT

      nfdh3   = savedVars%ithrr(gCT)
      nfdl2   = nfdh3  + 1
      nfdh2   = nfdh3  + savedVars%itwor(gCT)
      nfdl1   = nfdh2  + 1
      nfdh1   = nfdh2  + savedVars%ioner(gCT)
      nfdl0   = nfdh1  + 1
      nallr   = savedVars%nallrat(gCT)

      Tabl%nfdh3 = nfdh3
      Tabl%nfdl2 = nfdl2
      Tabl%nfdh2 = nfdh2
      Tabl%nfdh1 = nfdh1
      Tabl%nfdl1 = nfdl1
      Tabl%nfdl0 = nfdl0

      Tabl%RxnNumIn(1:nfdh3) = 3
      Tabl%RxnNumIn(nfdl2:nfdh2) = 2
      Tabl%RxnNumIn(nfdl1:nfdh1) = 1
      Tabl%RxnNumIn(nfdh1+1:) = 0
      Tabl%RxnNumOut = 0

      do nkn = 1, nallr
        nk = savedVars%noldfnew(nkn,gCT)
        Tabl%RxnIn(1,nkn) = savedVars%irm2(1,nk,gCT)
        Tabl%RxnIn(2,nkn) = savedVars%irm2(2,nk,gCT)
        Tabl%RxnIn(3,nkn) = savedVars%irm2(3,nk,gCT)
      end do

      !Calculate reaction table from spc->reaction mapping
      do i=savedVars%npllo(ncsp),savedVars%nplhi(ncsp)
         jspc = savedVars%jspnpl(i)
         do k=savedVars%npl5(i),savedVars%nph1(i)
            SpcRxnNum(jspc) = SpcRxnNum(jspc)+1
            SpcRxn(jspc,SpcRxnNum(jspc)) = savedVars%lossra(k)
         end do
         do k=savedVars%npl5(i),savedVars%nph2(i)
            SpcRxnNum(jspc) = SpcRxnNum(jspc)+1
            SpcRxn(jspc,SpcRxnNum(jspc)) = savedVars%lossrb(k)
         end do
         do k=savedVars%npl5(i),savedVars%nph3(i)
            SpcRxnNum(jspc) = SpcRxnNum(jspc)+1
            SpcRxn(jspc,SpcRxnNum(jspc)) = savedVars%lossrc(k)
         end do
         do k=savedVars%npl5(i),savedVars%nph4(i)
            SpcRxnNum(jspc) = SpcRxnNum(jspc)+1
            SpcRxn(jspc,SpcRxnNum(jspc)) = savedVars%lossrd(k)
         end do
         do k=savedVars%npl5(i),savedVars%nph5(i)
            SpcRxnNum(jspc) = SpcRxnNum(jspc)+1
            SpcRxn(jspc,SpcRxnNum(jspc)) = savedVars%lossre(k)
         end do
      end do

      !Convert original notation : Rxn # < nallr => input, Rxn # > nallr => output
      where (SpcRxn > nallr)
         SpcRxn = -(SpcRxn-nallr)
      end where
      SpcRxn = -SpcRxn
      !Now, negative values are input/positive values are output.  Abs value is reaction number

      do i=1,NNACT
         do j=1,SpcRxnNum(i)
            k = SpcRxn(i,j)
            if (k>0) then
               Tabl%RxnNumOut(k) = Tabl%RxnNumOut(k)+1
               Tabl%RxnOut(Tabl%RxnNumOut(k),k) = i
            end if
            if (k<0) then
               k = -k
               NumIntemp(k) = NumIntemp(k)+1
               !Keep track seperately to identify species with species on both sides
               RxIntemp(NumIntemp(k),k) = i
            end if

         end do
      end do

      !Now account for reactions where species is both input/output, ie Tabl%RxnNumIn>NumIntemp
      !Note, assumes that only one species is on both sides
      !TODO: generalize this

      do i=1,Nrxn
         if (Tabl%RxnNumIn(i) /= NumIntemp(i)) then
            Tabl%RxnNumOut(i) = Tabl%RxnNumOut(i)+1
            !Now find out which species is new output
            Sumtemp = sum(RxIntemp(1:NumIntemp(i),i))
            Sumin = sum(Tabl%RxnIn(1:Tabl%RxnNumIn(i),i))
            Tabl%RxnOut(Tabl%RxnNumOut(i),i) = Sumin-Sumtemp
         end if
      end do

      Tabl%isFracRxn = .false.
      NumFracOuts = 0
      FracOutSpc = 0
      FracRate = 0.0

      !Collect fractionated data from original structure
      do i=savedVars%nfrlo(ncsp),savedVars%nfrhi(ncsp)
         jspc = savedVars%jspcnfr(i)
         nkn = savedVars%nknfr(i)
         Tabl%isFracRxn(nkn) = .true.
         NumFracOuts(nkn) = NumFracOuts(nkn)+1
         FracOutSpc(NumFracOuts(nkn),nkn) = jspc
         FracRate(NumFracOuts(nkn),nkn) = savedVars%fracnfr(i)
      end do

      nkn = 0
      !Put fractionated data into new structure
      Tabl%FracRxns = 0
      Tabl%FracRxnNumOut = 0
      Tabl%FracRxnOut = 0
      Tabl%FracRate = 0.0

      do i=1,Nrxn
         if (Tabl%isFracRxn(i)) then
            nkn = nkn+1
            Tabl%FracRxns(nkn) = i
            Tabl%FracRxnNumOut(nkn) = NumFracOuts(i)
            Tabl%FracRxnOut(:,nkn) = FracOutSpc(:,i)
            Tabl%FracRate(:,nkn) = FracRate(:,i)
         end if
      end do

      !Call PrintChem(Tabl)
   End Subroutine CalcTabl

   Subroutine PrintChem(Tabl,Reorder, savedVars)
      type (ChemTable), intent(in) :: Tabl
      logical, optional, intent(in) :: Reorder
      type(t_Smv2Saved), intent(inOut) :: savedVars

      integer :: i,j,k,nfrac
      logical :: reord
#include "smv2chem_par.h"

      if (present(Reorder)) then
         reord = Reorder
      else
         reord = .true.
      endif

!Note, these are the reordered species numbers.  Can return back to original by using
!origSpcNumber(gasChemType=1,jreord) = savedVars%inewold from common

      nfrac = 0
      Write(*,*) 'Reaction table'
      do i=1,Nrxn
         if (Tabl%isFracRxn(i)) then
            nfrac = nfrac+1
            Print '("Reaction ", i4, " [",i2, " -> ",i2,"/",i2,"F]")', i,&
                 & Tabl%RxnNumIn(i),Tabl%RxnNumOut(i),Tabl%FracRxnNumOut(nfrac)
         else
            Print '("Reaction ", i4, " [",i2, " -> ",i2,"]")', i,&
                 & Tabl%RxnNumIn(i),Tabl%RxnNumOut(i)
         end if
         Print '("     In  : ",$)'
         do j=1,Tabl%RxnNumIn(i)
            k = Tabl%RxnIn(j,i)
            if (reord) k = savedVars%inewold(k,Tabl%gCT)
            Print '(i3,$)', k
            if (j /= Tabl%RxnNumIn(i)) Print '(" + ",$)'
         end do
         Print '("")'

         Print '("     Out : ",$)'
         do j=1,Tabl%RxnNumOut(i)
            k = Tabl%RxnOut(j,i)
            if (reord) k = savedVars%inewold(k,Tabl%gCT)
            Print '(i3,$)',k
            if (j /= Tabl%RxnNumOut(i)) Print '(" + ",$)'
         end do
         Print '("")'

         if (Tabl%isFracRxn(i)) then
            Print '("   F-Out : ",$)'
            do j=1,Tabl%FracRxnNumOut(nfrac)
               k = Tabl%FracRxnOut(j,nfrac)
               if (reord) k = savedVars%inewold(k,Tabl%gCT)
               Print '(i3,$)',k
               if (j /= Tabl%FracRxnNumOut(nfrac)) Print '(" + ",$)'
            end do
            Print '("")'
         end if
      end do


   End Subroutine PrintChem
end module ChemTable_mod

