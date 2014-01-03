!-----------------------------------------------------------------------------
!
! ROUTINE
!   doSmv2Solver
!
! DESCRIPTION
!   This is the main control routine for the ordinary differential equation
!   solver, "Smvgear II" (Sparse Matrix Vectorized Gear-type code).
! ! ARGUMENTS
!   doQqjkInchem   : if prQqjk is on, should qqj's & qqk's be determined
!                      inside the chemistry solver, or outside?
!   doSurfEmissInChem : do surface emissions inside the chemistry solver, or
!                      outside?
!   pr_diag          : print some diagnostic output to screen?
!   prQqjk          : should the periodic qqjk output file be written?
!   prSmv2          : should the SmvgearII     output file be written
!                      (non-parallel mode only)?
!   loc_proc         : local processor #
!   numLat             : # of latitudes
!   numLong            : # of longitudes
!   numVert            : # of vertical layers
!   itloop           : # of zones (numLong * numLat * numVert)
!   prNcPeriod     : NetCDF output period
!   timeStep              : model time step (s)
!   doCellChem     : do chemistry for a particular cell?
!   thermalRateConstants            : thermal    rate constants (units vary)
!   photolysisRateConstants            : photolysis rate constants (s^-1)
!   surfaceEmissions            : surface emissions (units?)
!   speciesConst               : spc conc (molec/cm^3)
!-----------------------------------------------------------------------------

      program doSmv2Solver
!  (savedVars, doQqjkInchem, &
!  doSurfEmissInChem, pr_diag, prQqjk, prSmv2,  &
!  loc_proc, numLat, numLong, numVert, itloop, prNcPeriod,&
!  timeStep, doCellChem, thermalRateConstants, photolysisRateConstants, surfaceEmissions, speciesConst, &
!  yda, qqkda, qqjda, i1, i2, ju1, j2, k1, k2, &
!  numQks, numQjs, numActive, commuWorld)

      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved
      use timing_mod

      implicit none

#     include "smv2chem_par.h"


!     ----------------------
!     Argument declarations in GMI
!     Local variables here
!     ----------------------

      integer :: commuWorld
      integer :: i1, i2, ju1, j2, k1, k2
      integer :: numQjo, numQks, numQjs, numActive

      real*8, allocatable :: qjGmi(:, :, :, :)
      real*8, allocatable :: qkGmi(:, :, :, :)
      real*8, allocatable :: qqjda(:,:,:,:)
      real*8, allocatable:: qqkda(:,:,:,:)
      real*8, allocatable :: yda (:,:,:,:)

      logical, save :: first = .true.
      logical  :: doQqjkInchem
      logical  :: doSurfEmissInChem
      logical  :: prQqjk
      logical  :: prSmv2
      integer  :: localProc
      integer  :: numLat, numLong, numVert
      integer  :: itloop
      real*8  :: prNcPeriod
      real*8  :: timeStep
      logical, allocatable  :: doCellChem(:)
      real*8, allocatable :: thermalRateConstants(:, :)
      real*8, allocatable :: photolysisRateConstants(:, :)
      real*8, allocatable :: surfaceEmissions(:, :)

      real*8, allocatable :: speciesConst(:, :)
      type(t_Smv2Saved) :: savedVars

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer, allocatable :: jreorder(:)
      integer, allocatable :: lreorder(:)

      real*8, allocatable  :: errmx2  (:)


!     Standalone-only declarations
      character(len=100) :: smv2Chem1Entry
      character(len=100) :: smv2Chem1Exit
      character(len=100) :: smv2Chem2Entry
      character(len=100) :: smv2Chem2Exit
      character(len=100) :: physProcEntry
      character(len=100) :: physProcExit
      character(:), allocatable :: rankString
      character(:), allocatable :: zeroString
      logical  :: prDiag
      integer :: rankSize, i

      integer, parameter :: MAX_RANK_SIZE = 4


!     ----------------
!     Begin execution.
!     ----------------
      call timingInit


      call get_command_argument(2, length=rankSize)
      print*, "rankSize: ", rankSize
      if (rankSize .gt. 3) then
         print*, "rank must be less than 1000"
         stop
      end if

      allocate(character(rankSize)::rankString)
      allocate(character(MAX_RANK_SIZE-rankSize)::zeroString)
      call get_command_argument(2, value=rankString)

      prDiag = .true.
      if (prDiag) then
        Write (6,*) 'doSmv2Solver called by ', rankString
      end if

      if (rankSize .eq. 1) zeroString = "000"
      if (rankSize .eq. 2) zeroString = "00"
      if (rankSize .eq. 3) zeroString = "0"

      smv2Chem1Entry = 'smv2chem1_entry.proc' // zeroString // rankString
      smv2Chem1Exit = 'smv2chem1_exit.proc' // zeroString // rankString
      smv2Chem2Entry = 'smv2chem2_entry.proc' // zeroString // rankString
      smv2Chem2Exit = 'smv2chem2_exit.proc' // zeroString // rankString
      physProcEntry = 'physproc_entry.proc' // zeroString // rankString
      physProcExit = 'physproc_exit.proc' // zeroString // rankString






      call initializeSavedVars(savedVars)
      call readSmv1Vars(savedVars,smv2Chem1Entry)
      call readSmv2Vars(savedVars,smv2Chem2Entry)



      print*, "reading: ", trim(physProcEntry)
      open(file=trim(physProcEntry),unit=27,form="formatted")
      read(27,*)
      read(27,*) i1, i2, ju1, j2, k1, k2
      read(27,*)
      read(27,*) numQjo, numQks, numQjs, numActive
      allocate(qjGmi(i1:i2, ju1:j2, k1:k2, numQjo))
      allocate(qkGmi(i1:i2, ju1:j2, k1:k2, numQks))
      allocate(qqjda(i1:i2, ju1:j2, k1:k2, numQjs))
      allocate(qqkda(i1:i2, ju1:j2, k1:k2, numQks))
      allocate(yda  (i1:i2, ju1:j2, k1:k2, numActive))
      read(27,*)
      read(27,*) qjGmi(i1:i2, ju1:j2, k1:k2, 1:numQjo)
      read(27,*)
      read(27,*) qkGmi(i1:i2, ju1:j2, k1:k2, 1:numQks)
      read(27,*)
      read(27,*) qqjda(i1:i2, ju1:j2, k1:k2, 1:numQjs)
      read(27,*)
      read(27,*) qqkda(i1:i2, ju1:j2, k1:k2, 1:numQks)
      read(27,*)
      read(27,*) yda  (i1:i2, ju1:j2, k1:k2, 1:numActive)
      read(27,*)
      read(27,*) doQqjkInchem
      read(27,*)
      read(27,*) doSurfEmissInChem
      read(27,*)
      read(27,*) prDiag
      read(27,*)
      read(27,*) prQqjk
      read(27,*)
      read(27,*) prSmv2
      read(27,*)
      read(27,*) localProc
      read(27,*)
      read(27,*) numLat, numLong, numVert
      read(27,*)
      read(27,*) itloop
      read(27,*)
      read(27,*) prNcPeriod
      read(27,*)
      read(27,*) timeStep

      allocate(doCellChem(itloop))
      allocate(thermalRateConstants(itloop, ITHERM))
      allocate(photolysisRateConstants(itloop, IPHOT))
      allocate(surfaceEmissions(numLat*numLong, IGAS))
      allocate(speciesConst(itloop, IGAS))

      read(27,*)
      read(27,*) doCellChem(1:itloop)
      read(27,*)
      read(27,*) thermalRateConstants(1:itloop, 1:ITHERM)
      read(27,*)
      read(27,*) photolysisRateConstants(1:itloop, 1:IPHOT)
      read(27,*)
      read(27,*) surfaceEmissions(1:numLat*numLong, 1:IGAS)
      read(27,*)
      read(27,*) speciesConst(1:itloop, 1:IGAS)
      close(27)



      if (first) then
         first = .false.
         Allocate (savedVars%csuma(itloop))
         Allocate (savedVars%csumc(itloop))
         savedVars%csuma = 0.0d0; savedVars%csumc = 0.0d0
      end if

      allocate (jreorder(itloop))
      allocate (lreorder(itloop))
      allocate (errmx2(itloop))
      jreorder(:) = 0; lreorder(:) = 0
      errmx2  (:) = 0.0d0

      commuWorld = 0

      print*, "Calling physproc"
      call Physproc (savedVars, doQqjkInchem, doSurfEmissInChem, &
               prQqjk, prSmv2, numLat, numLong, numVert, savedVars%ifreord,&
               savedVars%imgas, savedVars%initrogen, savedVars%ioxygen,&
               itloop, savedVars%kuloop, savedVars%lunsmv, &
               savedVars%ncs, savedVars%fracdec, savedVars%hmaxnit, &
               prNcPeriod, timeStep, doCellChem, savedVars%jphotrat, &
               savedVars%nrates, savedVars%ntloopncs, savedVars%ntspec,&
               savedVars%inewold, savedVars%npphotrat, thermalRateConstants, photolysisRateConstants, &
               surfaceEmissions, jreorder, lreorder, savedVars%csuma,  &
               savedVars%csumc, errmx2, speciesConst, yda, qqkda, qqjda, &
               i1, i2, ju1, j2, k1, k2, numQks, numQjs, numActive, &
               commuWorld)

      call writeSmv2Chem1Exit(smv2Chem1Exit,savedVars)
      call writeSmv2Chem2Exit(smv2Chem2Exit,savedVars)

      print*, "Writing to: ", physProcExit
      open(file=trim(physProcExit),unit=28,status="replace",form="unformatted")
      write(28) "i1, i2, ju1, j2, k1, k2"
      write(28) i1, i2, ju1, j2, k1, k2
      write(28) "num_qjo, num_qks, num_qjs, num_active"
      write(28) numQjo, numQks, numQjs, numActive
      write(28) "qjgmi(i1:i2, ju1:j2, k1:k2, num_qjo)"
      write(28) qjGmi(i1:i2, ju1:j2, k1:k2, 1:numQjo)
      write(28) "qkgmi(i1:i2, ju1:j2, k1:k2, num_qks)"
      write(28) qkGmi(:,:,:,:)
      write(28) "qqjda(i1:i2, ju1:j2, k1:k2, num_qjs)"
      write(28) qqjda(i1:i2, ju1:j2, k1:k2, 1:numQjs)
      write(28) "qqkda(i1:i2, ju1:j2, k1:k2, num_qks)"
      write(28) qqkda(:,:,:,:)
      write(28) "yda  (i1:i2, ju1:j2, k1:k2, num_active)"
      write(28) yda  (i1:i2, ju1:j2, k1:k2, 1:numActive)
      write(28) "do_qqjk_inchem"
      write(28) doQqjkInchem
      write(28) "do_semiss_inchem"
      write(28) doSurfEmissInChem
      write(28) "pr_diag"
      write(28) prDiag
      write(28) "pr_qqjk"
      write(28) prQqjk
      write(28) "pr_smv2"
      write(28) prSmv2
      write(28) "loc_proc"
      write(28) localProc
      write(28) "ilat, ilong, ivert"
      write(28) numLat, numLong, numVert
      write(28) "itloop"
      write(28) itloop
      write(28) "pr_nc_period"
      write(28) prNcPeriod
      write(28) "tdt"
      write(28) timeStep
      write(28) "do_cell_chem(itloop)"
      write(28) doCellChem(1:itloop)
      write(28) "arate(itloop, ITHERM)"
      write(28) thermalRateConstants(1:itloop, 1:ITHERM)
      write(28) "prate(itloop, IPHOT)"
      write(28) photolysisRateConstants(1:itloop, 1:IPHOT)
      write(28) "yemis(ilat*ilong, IGAS)"
      write(28) surfaceEmissions(1:numLat*numLong, 1:IGAS)
      write(28) "cx(itloop, IGAS)"
      write(28) speciesConst(1:itloop, 1:IGAS)
      close(28)

      deallocate (savedVars%csuma)
      deallocate (savedVars%csumc)
      deallocate (jreorder)
      deallocate (lreorder)
      deallocate (errmx2)
      deallocate(qjGmi)
      deallocate(qkGmi)
      deallocate(qqjda)
      deallocate(qqkda)
      deallocate(yda)


      deallocate(doCellChem)
      deallocate(thermalRateConstants)
      deallocate(photolysisRateConstants)
      deallocate(surfaceEmissions)
      deallocate(speciesConst)



      deallocate(savedVars%enqq1)
      deallocate(savedVars%enqq2)
      deallocate(savedVars%enqq3)
      deallocate(savedVars%conp15)
      deallocate(savedVars%conpst)

      deallocate(savedVars%aset)

      print*, "Exiting doSmv2Solver"

      end program doSmv2Solver

subroutine readSmv2Vars(savedVars,smv2Chem2Entry)

      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

      IMPLICIT NONE
#     include "smv2chem_par.h"

      type(t_Smv2Saved), intent(inOut) :: savedVars
      character(len=100), intent(in) :: smv2Chem2Entry

      print*, "reading: ", trim(smv2Chem2Entry)

!read smv2chem2 on entry to physproc to be used for standAlone code
      open(file=trim(smv2Chem2Entry),unit=25,form="formatted")
      read(25,*)
      read(25,*) savedVars%ioner
      read(25,*)
      read(25,*) savedVars%nallrat
      read(25,*)
      read(25,*) savedVars%inorep
      read(25,*)
      read(25,*) savedVars%ithrr
      read(25,*)
      read(25,*) savedVars%itwor
      read(25,*)
      read(25,*) savedVars%nm3bod
      read(25,*)
      read(25,*) savedVars%nmair
      read(25,*)
      read(25,*) savedVars%nmn2
      read(25,*)
      read(25,*) savedVars%nmo2
      read(25,*)
      read(25,*) savedVars%nmoth

      read(25,*)
      read(25,*) savedVars%ntrates
      read(25,*)
      read(25,*) savedVars%mappl
      read(25,*)
      read(25,*) savedVars%lgasbino
      read(25,*)
      read(25,*) savedVars%nreacoth
      read(25,*)
      read(25,*) savedVars%lgas3bod
      read(25,*)
      read(25,*) savedVars%losinacp
      read(25,*)
      read(25,*) savedVars%nreac3b
      read(25,*)
      read(25,*) savedVars%nreacair
      read(25,*)
      read(25,*) savedVars%nreacn2
      read(25,*)
      read(25,*) savedVars%nreaco2
      read(25,*)
      read(25,*) savedVars%jphotnk
      read(25,*)
      read(25,*) savedVars%noldfnew
      read(25,*)
      read(25,*) savedVars%irm2
      read(25,*)
      read(25,*) savedVars%ischang
      read(25,*)
      read(25,*) savedVars%kzthi
      read(25,*)
      read(25,*) savedVars%kztlo
      read(25,*)
      read(25,*) savedVars%ikztot
      read(25,*)
      read(25,*) savedVars%kbh1
      read(25,*)
      read(25,*) savedVars%kbh2
      read(25,*)
      read(25,*) savedVars%kbh3
      read(25,*)
      read(25,*) savedVars%kbh4
      read(25,*)
      read(25,*) savedVars%kbh5
      read(25,*)
      read(25,*) savedVars%kbl1
      read(25,*)
      read(25,*) savedVars%kbl2
      read(25,*)
      read(25,*) savedVars%kbl3
      read(25,*)
      read(25,*) savedVars%kbl4
      read(25,*)
      read(25,*) savedVars%kbl5
      read(25,*)
      read(25,*) savedVars%mbh1
      read(25,*)
      read(25,*) savedVars%mbh2
      read(25,*)
      read(25,*) savedVars%mbh3
      read(25,*)
      read(25,*) savedVars%mbh4
      read(25,*)
      read(25,*) savedVars%mbh5
      read(25,*)
      read(25,*) savedVars%mbl1
      read(25,*)
      read(25,*) savedVars%mbl2
      read(25,*)
      read(25,*) savedVars%mbl3
      read(25,*)
      read(25,*) savedVars%mbl4
      read(25,*)
      read(25,*) savedVars%mbl5
      read(25,*)
      read(25,*) savedVars%kzeroa
      read(25,*)
      read(25,*) savedVars%kzerob
      read(25,*)
      read(25,*) savedVars%kzeroc
      read(25,*)
      read(25,*) savedVars%kzerod
      read(25,*)
      read(25,*) savedVars%kzeroe
      read(25,*)
      read(25,*) savedVars%mzeroa
      read(25,*)
      read(25,*) savedVars%mzerob
      read(25,*)
      read(25,*) savedVars%mzeroc
      read(25,*)
      read(25,*) savedVars%mzerod
      read(25,*)
      read(25,*) savedVars%mzeroe
      read(25,*)
      read(25,*) savedVars%imztot
      read(25,*)
      read(25,*) savedVars%ijval
      read(25,*)
      read(25,*) savedVars%jzeroa
      read(25,*)
      read(25,*) savedVars%idh1
      read(25,*)
      read(25,*) savedVars%idh2
      read(25,*)
      read(25,*) savedVars%idh3
      read(25,*)
      read(25,*) savedVars%idh4
      read(25,*)
      read(25,*) savedVars%idh5
      read(25,*)
      read(25,*) savedVars%idl1
      read(25,*)
      read(25,*) savedVars%idl2
      read(25,*)
      read(25,*) savedVars%idl3
      read(25,*)
      read(25,*) savedVars%idl4
      read(25,*)
      read(25,*) savedVars%idl5
      read(25,*)
      read(25,*) savedVars%ikdeca
      read(25,*)
      read(25,*) savedVars%ikdecb
      read(25,*)
      read(25,*) savedVars%ikdecc
      read(25,*)
      read(25,*) savedVars%ikdecd
      read(25,*)
      read(25,*) savedVars%ikdece
      read(25,*)
      read(25,*) savedVars%kjdeca
      read(25,*)
      read(25,*) savedVars%kjdecb
      read(25,*)
      read(25,*) savedVars%kjdecc
      read(25,*)
      read(25,*) savedVars%kjdecd
      read(25,*)
      read(25,*) savedVars%kjdece
      read(25,*)
      read(25,*) savedVars%ijthi
      read(25,*)
      read(25,*) savedVars%ijtlo
      read(25,*)
      read(25,*) savedVars%jarrdiag
      read(25,*)
      read(25,*) savedVars%jhiz1
      read(25,*)
      read(25,*) savedVars%jloz1
      read(25,*)
      read(25,*) savedVars%iarray
      read(25,*)
      read(25,*) savedVars%npdhi
      read(25,*)
      read(25,*) savedVars%npdlo
      read(25,*)
      read(25,*) savedVars%iialpd
      read(25,*)
      read(25,*) savedVars%ipospd
      read(25,*)
      read(25,*) savedVars%nkpdterm
      read(25,*)
      read(25,*) savedVars%nfrhi
      read(25,*)
      read(25,*) savedVars%nfrlo
      read(25,*)
      read(25,*) savedVars%nplhi
      read(25,*)
      read(25,*) savedVars%npllo
      read(25,*)
      read(25,*) savedVars%jspcnfr
      read(25,*)
      read(25,*) savedVars%jspnpl
      read(25,*)
      read(25,*) savedVars%nknfr
      read(25,*)
      read(25,*) savedVars%lossra
      read(25,*)
      read(25,*) savedVars%lossrb
      read(25,*)
      read(25,*) savedVars%lossrc
      read(25,*)
      read(25,*) savedVars%lossrd
      read(25,*)
      read(25,*) savedVars%lossre
      read(25,*)
      read(25,*) savedVars%nph1
      read(25,*)
      read(25,*) savedVars%nph2
      read(25,*)
      read(25,*) savedVars%nph3
      read(25,*)
      read(25,*) savedVars%nph4
      read(25,*)
      read(25,*) savedVars%nph5
      read(25,*)
      read(25,*) savedVars%npl1
      read(25,*)
      read(25,*) savedVars%npl2
      read(25,*)
      read(25,*) savedVars%npl3
      read(25,*)
      read(25,*) savedVars%npl4
      read(25,*)
      read(25,*) savedVars%npl5
      read(25,*)
      read(25,*) savedVars%nolosp
      read(25,*)
      read(25,*) savedVars%newfold
      read(25,*)
      read(25,*) savedVars%nknlosp
      read(25,*)
      read(25,*) savedVars%nknphotrt
      read(25,*)
      read(25,*) savedVars%abst2
      read(25,*)
      read(25,*) savedVars%errmax
      read(25,*)
      read(25,*) savedVars%hmaxday
      read(25,*)
      read(25,*) savedVars%timeintv
      read(25,*)
      read(25,*) savedVars%abtol
      read(25,*)



      allocate(savedVars%enqq1 (MORDER))
      allocate(savedVars%enqq2 (MORDER))
      allocate(savedVars%enqq3 (MORDER))
      allocate(savedVars%conp15(MORDER))
      allocate(savedVars%conpst(MORDER))
      allocate(savedVars%pertst2(MORDER, 3))

      read(25,*) savedVars%enqq1
      read(25,*)
      read(25,*) savedVars%enqq2
      read(25,*)
      read(25,*) savedVars%enqq3
      read(25,*)
      read(25,*) savedVars%conp15
      read(25,*)
      read(25,*) savedVars%conpst
      read(25,*)
      read(25,*) savedVars%pertst2

      allocate(savedVars%aset(10, 8))
      !allocate(savedVars%fracpl (MXCOUNT2))
      !allocate(savedVars%fracnfr(MXCOUNT4))

      read(25,*)
      read(25,*) savedVars%aset
      read(25,*)

      read(25,*) savedVars%fracpl(:)
      read(25,*)
      read(25,*) savedVars%fracnfr(:)

      close(25)
   end subroutine readSmv2Vars

 !read smv2chem1 on entry to physproc to be used for standAlone code
  subroutine readSmv1Vars(savedVars,smv2Chem1Entry)

      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

      IMPLICIT NONE

      type(t_Smv2Saved), intent(inOut) :: savedVars
      character(len=100), intent(in) :: smv2Chem1Entry

      print*, "reading: ", trim(smv2Chem1Entry)

      open(file=trim(smv2Chem1Entry),unit=23,form="formatted")
      read(23,*)
      read(23,*) savedVars%ifreord
      read(23,*)
      read(23,*) savedVars%ih2o
      read(23,*)
      read(23,*) savedVars%imgas
      read(23,*)
      read(23,*) savedVars%initrogen
      read(23,*)
      read(23,*) savedVars%ioxygen
      read(23,*)
      read(23,*) savedVars%kuloop
      read(23,*)
      read(23,*) savedVars%lunsmv
      read(23,*)
      read(23,*) savedVars%ncs
      read(23,*)
      read(23,*) savedVars%jphotrat
      read(23,*)
      read(23,*) savedVars%nrates
      read(23,*)
      read(23,*) savedVars%ntloopncs
      read(23,*)
      read(23,*) savedVars%ntspec
      read(23,*)
      read(23,*) savedVars%inewold
      read(23,*)
      read(23,*) savedVars%npphotrat
      read(23,*)
      read(23,*) savedVars%fracdec
      read(23,*)
      read(23,*) savedVars%hmaxnit
      close(23)

  end subroutine readSmv1Vars


  subroutine printSmv1Vars(savedVars)

    use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

      IMPLICIT NONE

      type(t_Smv2Saved), intent(inOut) :: savedVars
      print*,  "ifreord: ", savedVars%ifreord
      print*,  "ih2o: ", savedVars%ih2o
      print*,  "imgas: ", savedVars%imgas
      print*,  "initrogen: ", savedVars%initrogen
      print*,  "ioxygen: ", savedVars%ioxygen
      print*,  "kuloop:", savedVars%kuloop
      print*,  "lunsmv: ", savedVars%lunsmv
      print*,  "ncs: ", savedVars%ncs
      print*,  "jphotrat:", savedVars%jphotrat
      print*,  "nrates: ", savedVars%nrates
      print*,  "ntloopncs: ", savedVars%ntloopncs
      print*,  "ntspec: ", savedVars%ntspec
      print*,  "inewold: ", savedVars%inewold
      print*,  "npphotrat: ", savedVars%npphotrat
      print*,  "fracdec: ", savedVars%fracdec
      print*,  "hmaxit: ", savedVars%hmaxnit

end subroutine printSmv1Vars

  subroutine printSmv2Vars(savedVars)

    use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

      IMPLICIT NONE

      type(t_Smv2Saved), intent(inOut) :: savedVars
      print*, "nallrat: ", savedVars%nallrat
      print*, "inorep: ", savedVars%inorep
      print*, "ithrr: ", savedVars%ithrr
      print*, "itwor: ", savedVars%itwor
      print*, "nm3bod: ", savedVars%nm3bod
      print*, "nmair: ", savedVars%nmair
      print*, "nmn2: ", savedVars%nmn2
      print*, "nmo2: ", savedVars%nmo2
      print*, "nmoth: ", savedVars%nmoth
      print*, "ntrates: ", savedVars%ntrates
      print*, "mappl: ", savedVars%mappl
      print*, "lgasbino: ", savedVars%lgasbino
      print*, "nreacoth: ", savedVars%nreacoth
      print*, "lgas3bod: ", savedVars%lgas3bod
      print*, "losinacp: ", savedVars%losinacp
      print*, "nreac3b: ", savedVars%nreac3b
      print*, "nreacair: ", savedVars%nreacair
      print*, "nreacn2: ", savedVars%nreacn2
      print*, "nreaco2: ", savedVars%nreaco2
      print*, "jphotnk: ", savedVars%jphotnk
      print*, "noldfnew: ", savedVars%noldfnew
      print*, "irm2: ", savedVars%irm2
      print*, "ischang: ", savedVars%ischang
      print*, "kzthi: ", savedVars%kzthi
      print*, "kztlo: ", savedVars%kztlo
      print*, "ikztot: ", savedVars%ikztot
      print*, "kbh1: ", savedVars%kbh1
      print*, "kbh2: ", savedVars%kbh2
      print*, "kbh3: ", savedVars%kbh3
      print*, "kbh4: ", savedVars%kbh4
      print*, "kbh5: ", savedVars%kbh5
      print*, "kbl1: ", savedVars%kbl1
      print*, "kbl2: ", savedVars%kbl2
      print*, "kbl3: ", savedVars%kbl3
      print*, "kbl4: ", savedVars%kbl4
      print*, "kbl5: ", savedVars%kbl5
      print*, "mbh1: ", savedVars%mbh1
      print*, "mbh2: ", savedVars%mbh2
      print*, "mbh3: ", savedVars%mbh3
      print*, "mbh4: ", savedVars%mbh4
      print*, "mbh5: ", savedVars%mbh5
      print*, "mbl1: ", savedVars%mbl1
      print*, "mbl2: ", savedVars%mbl2
      print*, "mbl3: ", savedVars%mbl3
      print*, "mbl4: ", savedVars%mbl4
      print*, "mbl5: ", savedVars%mbl5
      print*, "kzeroa: ", savedVars%kzeroa
      print*, "kzerob: ", savedVars%kzerob
      print*, "kzeroc: ", savedVars%kzeroc
      print*, "kzerod: ", savedVars%kzerod
      print*, "kzeroe: ", savedVars%kzeroe
      print*, "mzeroa: ", savedVars%mzeroa
      print*, "mzerob: ", savedVars%mzerob
      print*, "mzeroc: ", savedVars%mzeroc
      print*, "mzerod: ", savedVars%mzerod
      print*, "mzeroe: ", savedVars%mzeroe
      print*, "imztot: ", savedVars%imztot
      print*, "ijval: ", savedVars%ijval
      print*, "jzeroa: ", savedVars%jzeroa
      print*, "idh1: ", savedVars%idh1
      print*, "idh2: ", savedVars%idh2
      print*, "idh3: ", savedVars%idh3
      print*, "idh4: ", savedVars%idh4
      print*, "idh5: ", savedVars%idh5
      print*, "idl1: ", savedVars%idl1
      print*, "idl2: ", savedVars%idl2
      print*, "idl3: ", savedVars%idl3
      print*, "idl4: ", savedVars%idl4
      print*, "idl5: ", savedVars%idl5
      print*, "ikdeca: ", savedVars%ikdeca
      print*, "ikdecb: ", savedVars%ikdecb
      print*, "ikdecc: ", savedVars%ikdecc
      print*, "ikdecd: ", savedVars%ikdecd
      print*, "ikdece: ", savedVars%ikdece
      print*, "kjdeca: ", savedVars%kjdeca
      print*, "kjdecb: ", savedVars%kjdecb
      print*, "kjdecc: ", savedVars%kjdecc
      print*, "kjdecd: ", savedVars%kjdecd
      print*, "kjdece: ", savedVars%kjdece
      print*, "ijthi: ", savedVars%ijthi
      print*, "ijtlo: ", savedVars%ijtlo
      print*, "jarrdiag: ", savedVars%jarrdiag
      print*, "jhiz1: ", savedVars%jhiz1
      print*, "jloz1: ", savedVars%jloz1
      print*, "iarray: ", savedVars%iarray
      print*, "npdhi: ", savedVars%npdhi
      print*, "npdlo: ", savedVars%npdlo
      print*, "iialpd: ", savedVars%iialpd
      print*, "ipospd: ", savedVars%ipospd
      print*, "nkpdterm: ", savedVars%nkpdterm
      print*, "nfrhi: ", savedVars%nfrhi
      print*, "nfrlo: ", savedVars%nfrlo
      print*, "nplhi: ", savedVars%nplhi
      print*, "npllo: ", savedVars%npllo
      print*, "jspcnfr: ", savedVars%jspcnfr
      print*, "jspnpl: ", savedVars%jspnpl
      print*, "nknfr: ", savedVars%nknfr
      print*, "lossra: ", savedVars%lossra
      print*, "lossrb: ", savedVars%lossrb
      print*, "lossrc: ", savedVars%lossrc
      print*, "lossrd: ", savedVars%lossrd
      print*, "lossre: ", savedVars%lossre
      print*, "nph1: ", savedVars%nph1
      print*, "nph2: ", savedVars%nph2
      print*, "nph3: ", savedVars%nph3
      print*, "nph4: ", savedVars%nph4
      print*, "nph5: ", savedVars%nph5
      print*, "npl1: ", savedVars%npl1
      print*, "npl2: ", savedVars%npl2
      print*, "npl3: ", savedVars%npl3
      print*, "npl4: ", savedVars%npl4
      print*, "npl5: ", savedVars%npl5
      print*, "nolosp: ", savedVars%nolosp
      print*, "newfold: ", savedVars%newfold
      print*, "nknlosp: ", savedVars%nknlosp
      print*, "nknphotrt: ", savedVars%nknphotrt
      print*, "abst2: ", savedVars%abst2
      print*, "errmax: ", savedVars%errmax
      print*, "hmaxday: ", savedVars%hmaxday
      print*, "timeintv: ", savedVars%timeintv
      print*, "abtol: ", savedVars%abtol
      print*, "nallrat: ", savedVars%nallrat
      print*, "inorep: ", savedVars%inorep
      print*, "ithrr: ", savedVars%ithrr
      print*, "itwor: ", savedVars%itwor
      print*, "nm3bod: ", savedVars%nm3bod
      print*, "nmair: ", savedVars%nmair
      print*, "nmn2: ", savedVars%nmn2
      print*, "nmo2: ", savedVars%nmo2
      print*, "nmoth: ", savedVars%nmoth
      print*, "ntrates: ", savedVars%ntrates
      print*, "mappl: ", savedVars%mappl
      print*, "lgasbino: ", savedVars%lgasbino
      print*, "nreacoth: ", savedVars%nreacoth
      print*, "lgas3bod: ", savedVars%lgas3bod
      print*, "losinacp: ", savedVars%losinacp
      print*, "nreac3b: ", savedVars%nreac3b
      print*, "nreacair: ", savedVars%nreacair
      print*, "nreacn2: ", savedVars%nreacn2
      print*, "nreaco2: ", savedVars%nreaco2
      print*, "jphotnk: ", savedVars%jphotnk
      print*, "noldfnew: ", savedVars%noldfnew
      print*, "irm2: ", savedVars%irm2
      print*, "ischang: ", savedVars%ischang
      print*, "kzthi: ", savedVars%kzthi
      print*, "kztlo: ", savedVars%kztlo
      print*, "ikztot: ", savedVars%ikztot
      print*, "kbh1: ", savedVars%kbh1
      print*, "kbh2: ", savedVars%kbh2
      print*, "kbh3: ", savedVars%kbh3
      print*, "kbh4: ", savedVars%kbh4
      print*, "kbh5: ", savedVars%kbh5
      print*, "kbl1: ", savedVars%kbl1
      print*, "kbl2: ", savedVars%kbl2
      print*, "kbl3: ", savedVars%kbl3
      print*, "kbl4: ", savedVars%kbl4
      print*, "kbl5: ", savedVars%kbl5
      print*, "mbh1: ", savedVars%mbh1
      print*, "mbh2: ", savedVars%mbh2
      print*, "mbh3: ", savedVars%mbh3
      print*, "mbh4: ", savedVars%mbh4
      print*, "mbh5: ", savedVars%mbh5
      print*, "mbl1: ", savedVars%mbl1
      print*, "mbl2: ", savedVars%mbl2
      print*, "mbl3: ", savedVars%mbl3
      print*, "mbl4: ", savedVars%mbl4
      print*, "mbl5: ", savedVars%mbl5
      print*, "kzeroa: ", savedVars%kzeroa
      print*, "kzerob: ", savedVars%kzerob
      print*, "kzeroc: ", savedVars%kzeroc
      print*, "kzerod: ", savedVars%kzerod
      print*, "kzeroe: ", savedVars%kzeroe
      print*, "mzeroa: ", savedVars%mzeroa
      print*, "mzerob: ", savedVars%mzerob
      print*, "mzeroc: ", savedVars%mzeroc
      print*, "mzerod: ", savedVars%mzerod
      print*, "mzeroe: ", savedVars%mzeroe
      print*, "imztot: ", savedVars%imztot
      print*, "ijval: ", savedVars%ijval
      print*, "jzeroa: ", savedVars%jzeroa
      print*, "idh1: ", savedVars%idh1
      print*, "idh2: ", savedVars%idh2
      print*, "idh3: ", savedVars%idh3
      print*, "idh4: ", savedVars%idh4
      print*, "idh5: ", savedVars%idh5
      print*, "idl1: ", savedVars%idl1
      print*, "idl2: ", savedVars%idl2
      print*, "idl3: ", savedVars%idl3
      print*, "idl4: ", savedVars%idl4
      print*, "idl5: ", savedVars%idl5
      print*, "ikdeca: ", savedVars%ikdeca
      print*, "ikdecb: ", savedVars%ikdecb
      print*, "ikdecc: ", savedVars%ikdecc
      print*, "ikdecd: ", savedVars%ikdecd
      print*, "ikdece: ", savedVars%ikdece
      print*, "kjdeca: ", savedVars%kjdeca
      print*, "kjdecb: ", savedVars%kjdecb
      print*, "kjdecc: ", savedVars%kjdecc
      print*, "kjdecd: ", savedVars%kjdecd
      print*, "kjdece: ", savedVars%kjdece
      print*, "ijthi: ", savedVars%ijthi
      print*, "ijtlo: ", savedVars%ijtlo
      print*, "jarrdiag: ", savedVars%jarrdiag
      print*, "jhiz1: ", savedVars%jhiz1
      print*, "jloz1: ", savedVars%jloz1
      print*, "iarray: ", savedVars%iarray
      print*, "npdhi: ", savedVars%npdhi
      print*, "npdlo: ", savedVars%npdlo
      print*, "iialpd: ", savedVars%iialpd
      print*, "ipospd: ", savedVars%ipospd
      print*, "nkpdterm: ", savedVars%nkpdterm
      print*, "nfrhi: ", savedVars%nfrhi
      print*, "nfrlo: ", savedVars%nfrlo
      print*, "nplhi: ", savedVars%nplhi
      print*, "npllo: ", savedVars%npllo
      print*, "jspcnfr: ", savedVars%jspcnfr
      print*, "jspnpl: ", savedVars%jspnpl
      print*, "nknfr: ", savedVars%nknfr
      print*, "lossra: ", savedVars%lossra
      print*, "lossrb: ", savedVars%lossrb
      print*, "lossrc: ", savedVars%lossrc
      print*, "lossrd: ", savedVars%lossrd
      print*, "lossre: ", savedVars%lossre
      print*, "nph1: ", savedVars%nph1
      print*, "nph2: ", savedVars%nph2
      print*, "nph3: ", savedVars%nph3
      print*, "nph4: ", savedVars%nph4
      print*, "nph5: ", savedVars%nph5
      print*, "npl1: ", savedVars%npl1
      print*, "npl2: ", savedVars%npl2
      print*, "npl3: ", savedVars%npl3
      print*, "npl4: ", savedVars%npl4
      print*, "npl5: ", savedVars%npl5
      print*, "nolosp: ", savedVars%nolosp
      print*, "newfold: ", savedVars%newfold
      print*, "nknlosp: ", savedVars%nknlosp
      print*, "nknphotrt: ", savedVars%nknphotrt
      print*, "abst2: ", savedVars%abst2
      print*, "errmax: ", savedVars%errmax
      print*, "hmaxday: ", savedVars%hmaxday
      print*, "timeintv: ", savedVars%timeintv
      print*, "abtol: ", savedVars%abtol


end subroutine printSmv2Vars

   subroutine initializeSavedVars(savedVars)

      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

      IMPLICIT NONE

      type(t_Smv2Saved), intent(inOut) :: savedVars

      integer, parameter :: initVal = 0

      ! Begin initialize integer arrays to a default value
      savedVars%inewold   = initVal
      savedVars%npphotrat = initVal
      savedVars%mappl     = initVal
      savedVars%lgasbino  = initVal
      savedVars%nreacoth  = initVal
      savedVars%lgas3bod  = initVal
      savedVars%losinacp  = initVal
      savedVars%nreac3b   = initVal
      savedVars%nreacair  = initVal
      savedVars%nreacn2   = initVal
      savedVars%nreaco2   = initVal
      savedVars%jphotnk   = initVal
      savedVars%noldfnew  = initVal
      savedVars%irm2      = initVal

      savedVars%ikztot    = initVal

      savedVars%kbh1      = initVal
      savedVars%kbh2      = initVal
      savedVars%kbh3      = initVal
      savedVars%kbh4      = initVal
      savedVars%kbh5      = initVal
      savedVars%kbl1      = initVal
      savedVars%kbl2      = initVal
      savedVars%kbl3      = initVal
      savedVars%kbl4      = initVal
      savedVars%kbl5      = initVal

      savedVars%mbh1      = initVal
      savedVars%mbh2      = initVal
      savedVars%mbh3      = initVal
      savedVars%mbh4      = initVal
      savedVars%mbh5      = initVal
      savedVars%mbl1      = initVal
      savedVars%mbl2      = initVal
      savedVars%mbl3      = initVal
      savedVars%mbl4      = initVal
      savedVars%mbl5      = initVal

      savedVars%idh1      = initVal
      savedVars%idh2      = initVal
      savedVars%idh3      = initVal
      savedVars%idh4      = initVal
      savedVars%idh5      = initVal
      savedVars%idl1      = initVal
      savedVars%idl2      = initVal
      savedVars%idl3      = initVal
      savedVars%idl4      = initVal
      savedVars%idl5      = initVal

      savedVars%nph1      = initVal
      savedVars%nph2      = initVal
      savedVars%nph3      = initVal
      savedVars%nph4      = initVal
      savedVars%nph5      = initVal
      savedVars%npl1      = initVal
      savedVars%npl2      = initVal
      savedVars%npl3      = initVal
      savedVars%npl4      = initVal
      savedVars%npl5      = initVal

      savedVars%kzeroa    = initVal
      savedVars%kzerob    = initVal
      savedVars%kzeroc    = initVal
      savedVars%kzerod    = initVal
      savedVars%kzeroe    = initVal

      savedVars%mzeroa    = initVal
      savedVars%mzerob    = initVal
      savedVars%mzeroc    = initVal
      savedVars%mzerod    = initVal
      savedVars%mzeroe    = initVal

      savedVars%imztot    = initVal

      savedVars%ijval     = initVal
      savedVars%jzeroa    = initVal

      savedVars%ikdeca = initVal
      savedVars%ikdecb = initVal
      savedVars%ikdecc = initVal
      savedVars%ikdecd = initVal
      savedVars%ikdece = initVal

      savedVars%kjdeca = initVal
      savedVars%kjdecb = initVal
      savedVars%kjdecc = initVal
      savedVars%kjdecd = initVal
      savedVars%kjdece = initVal

      savedVars%ijthi    = initVal
      savedVars%ijtlo    = initVal
      savedVars%jarrdiag = initVal
      savedVars%jhiz1    = initVal
      savedVars%jloz1    = initVal
      savedVars%iialpd   = initVal
      savedVars%ipospd   = initVal
      savedVars%nkpdterm = initVal
      savedVars%jspcnfr  = initVal

      savedVars%jspnpl = initVal
      savedVars%nknfr  = initVal

      savedVars%lossra = initVal
      savedVars%lossrb = initVal
      savedVars%lossrc = initVal
      savedVars%lossrd = initVal
      savedVars%lossre = initVal

      savedVars%newfold   = initVal
      savedVars%nknlosp   = initVal
      savedVars%nknphotrt = initVal

      savedVars%abtol   = 0.0d0
      savedVars%fracpl  = 0.0d0
      savedVars%fracnfr = 0.0d0
      ! End initialize integer arrays to a default value

      RETURN

      end subroutine initializeSavedVars


      subroutine writeSmv2Chem1Exit (fileName, savedVars)
         use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

         IMPLICIT NONE

         type(t_Smv2Saved), intent(inOut) :: savedVars
         character(len=100), intent(in) :: fileName

         integer :: fileNumber

         print*, "Writing to: ", fileName
         open(newunit=fileNumber, file=trim(fileName),status="replace",form="unformatted")

         write(fileNumber) "ifreord   "
         write(fileNumber) savedVars%ifreord
         write(fileNumber) "ih2o      "
         write(fileNumber) savedVars%ih2o
         write(fileNumber) "imgas     "
         write(fileNumber) savedVars%imgas
         write(fileNumber) "initrogen "
         write(fileNumber) savedVars%initrogen
         write(fileNumber) "ioxygen   "
         write(fileNumber) savedVars%ioxygen
         write(fileNumber) "kuloop    "
         write(fileNumber) savedVars%kuloop
         write(fileNumber) "lunsmv    "
         write(fileNumber) savedVars%lunsmv
         write(fileNumber) "ncs       "
         write(fileNumber) savedVars%ncs
         write(fileNumber) "jphotrat (ICS)"
         write(fileNumber) savedVars%jphotrat
         write(fileNumber) "nrates   (ICS)            "
         write(fileNumber) savedVars%nrates
         write(fileNumber) "ntloopncs(ICS)            "
         write(fileNumber) savedVars%ntloopncs
         write(fileNumber) "ntspec   (ICS)            "
         write(fileNumber) savedVars%ntspec
         write(fileNumber) "inewold  (MXGSAER, ICS)   "
         write(fileNumber) savedVars%inewold
         write(fileNumber) "npphotrat(IPHOT,   ICS)   "
         write(fileNumber) savedVars%npphotrat
         write(fileNumber) "fracdec   "
         write(fileNumber) savedVars%fracdec
         write(fileNumber) "hmaxnit   "
         write(fileNumber) savedVars%hmaxnit

         close(fileNumber)

      end subroutine writeSmv2Chem1Exit


      subroutine writeSmv2Chem2Exit (fileName, savedVars)
         use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

         IMPLICIT NONE

         type(t_Smv2Saved), intent(inOut) :: savedVars
         character(len=100), intent(in) :: fileName

         integer :: fileNumber

         print*, "Writing to: ", fileName
         open(newunit=fileNumber, file=trim(fileName),status="replace",form="unformatted")

         write(fileNumber) "ioner(ICP)"
         write(fileNumber) savedVars%ioner
         write(fileNumber) "nallrat(ICP)"
         write(fileNumber) savedVars%nallrat
         write(fileNumber) "inorep(ICS)"
         write(fileNumber) savedVars%inorep
         write(fileNumber) "ithrr(ICS)"
         write(fileNumber) savedVars%ithrr
         write(fileNumber) "itwor(ICS)"
         write(fileNumber) savedVars%itwor
         write(fileNumber) "nm3bod(ICS)"
         write(fileNumber) savedVars%nm3bod
         write(fileNumber) "nmair(ICS)"
         write(fileNumber) savedVars%nmair
         write(fileNumber) "nmn2(ICS)"
         write(fileNumber) savedVars%nmn2
         write(fileNumber) "nmo2(ICS)"
         write(fileNumber) savedVars%nmo2
         write(fileNumber) "nmoth(ICS)"
         write(fileNumber) savedVars%nmoth
         write(fileNumber) "ntrates(ICS)"
         write(fileNumber) savedVars%ntrates
         write(fileNumber) "mappl(MXGSAER,ICS)"
         write(fileNumber) savedVars%mappl
         write(fileNumber) "lgasbino(MAXGL2,ICS)"
         write(fileNumber) savedVars%lgasbino
         write(fileNumber) "nreacoth(MAXGL2,ICS)"
         write(fileNumber) savedVars%nreacoth
         write(fileNumber) "lgas3bod(MAXGL3,ICS)"
         write(fileNumber) savedVars%lgas3bod
         write(fileNumber) "losinacp(MAXGL3,ICS)"
         write(fileNumber) savedVars%losinacp
         write(fileNumber) "nreac3b(MAXGL3,ICS)"
         write(fileNumber) savedVars%nreac3b
         write(fileNumber) "nreacair(MAXGL3,ICS)"
         write(fileNumber) savedVars%nreacair
         write(fileNumber) "nreacn2(MAXGL3,ICS)"
         write(fileNumber) savedVars%nreacn2
         write(fileNumber) "nreaco2(MAXGL3,ICS)"
         write(fileNumber) savedVars%nreaco2
         write(fileNumber) "jphotnk(NMTRATE,ICS)"
         write(fileNumber) savedVars%jphotnk
         write(fileNumber) "noldfnew(NMTRATE,ICS)"
         write(fileNumber) savedVars%noldfnew
         write(fileNumber) "irm2(NMRPROD,NMTRATE,ICS)"
         write(fileNumber) savedVars%irm2
         write(fileNumber) "ischang(ICS)"
         write(fileNumber) savedVars%ischang
         write(fileNumber) "kzthi(ICP)"
         write(fileNumber) savedVars%kzthi
         write(fileNumber) "kztlo(ICP)"
         write(fileNumber) savedVars%kztlo
         write(fileNumber) "ikztot(MXCOUNT4)"
         write(fileNumber) savedVars%ikztot
         write(fileNumber) "kbh1(MXCOUNT4)"
         write(fileNumber) savedVars%kbh1
         write(fileNumber) "kbh2(MXCOUNT4)"
         write(fileNumber) savedVars%kbh2
         write(fileNumber) "kbh3(MXCOUNT4)"
         write(fileNumber) savedVars%kbh3
         write(fileNumber) "kbh4(MXCOUNT4)"
         write(fileNumber) savedVars%kbh4
         write(fileNumber) "kbh5(MXCOUNT4)"
         write(fileNumber) savedVars%kbh5
         write(fileNumber) "kbl1(MXCOUNT4)"
         write(fileNumber) savedVars%kbl1
         write(fileNumber) "kbl2(MXCOUNT4)"
         write(fileNumber) savedVars%kbl2
         write(fileNumber) "kbl3(MXCOUNT4)"
         write(fileNumber) savedVars%kbl3
         write(fileNumber) "kbl4(MXCOUNT4)"
         write(fileNumber) savedVars%kbl4
         write(fileNumber) "kbl5(MXCOUNT4)"
         write(fileNumber) savedVars%kbl5
         write(fileNumber) "mbh1(MXCOUNT4)"
         write(fileNumber) savedVars%mbh1
         write(fileNumber) "mbh2(MXCOUNT4)"
         write(fileNumber) savedVars%mbh2
         write(fileNumber) "mbh3(MXCOUNT4)"
         write(fileNumber) savedVars%mbh3
         write(fileNumber) "mbh4(MXCOUNT4)"
         write(fileNumber) savedVars%mbh4
         write(fileNumber) "mbh5(MXCOUNT4)"
         write(fileNumber) savedVars%mbh5
         write(fileNumber) "mbl1(MXCOUNT4)"
         write(fileNumber) savedVars%mbl1
         write(fileNumber) "mbl2(MXCOUNT4)"
         write(fileNumber) savedVars%mbl2
         write(fileNumber) "mbl3(MXCOUNT4)"
         write(fileNumber) savedVars%mbl3
         write(fileNumber) "mbl4(MXCOUNT4)"
         write(fileNumber) savedVars%mbl4
         write(fileNumber) "mbl5(MXCOUNT4)"
         write(fileNumber) savedVars%mbl5
         write(fileNumber) "kzeroa(MXCOUNT4)"
         write(fileNumber) savedVars%kzeroa
         write(fileNumber) "kzerob(MXCOUNT4)"
         write(fileNumber) savedVars%kzerob
         write(fileNumber) "kzeroc(MXCOUNT4)"
         write(fileNumber) savedVars%kzeroc
         write(fileNumber) "kzerod(MXCOUNT4)"
         write(fileNumber) savedVars%kzerod
         write(fileNumber) "kzeroe(MXCOUNT4)"
         write(fileNumber) savedVars%kzeroe
         write(fileNumber) "mzeroa(MXCOUNT4)"
         write(fileNumber) savedVars%mzeroa
         write(fileNumber) "mzerob(MXCOUNT4)"
         write(fileNumber) savedVars%mzerob
         write(fileNumber) "mzeroc(MXCOUNT4)"
         write(fileNumber) savedVars%mzeroc
         write(fileNumber) "mzerod(MXCOUNT4)"
         write(fileNumber) savedVars%mzerod
         write(fileNumber) "mzeroe(MXCOUNT4)"
         write(fileNumber) savedVars%mzeroe
         write(fileNumber) "imztot(MXGSAER,ICP)"
         write(fileNumber) savedVars%imztot
         write(fileNumber) "ijval(MXCOUNT3)"
         write(fileNumber) savedVars%ijval
         write(fileNumber) "jzeroa(MXCOUNT3)"
         write(fileNumber) savedVars%jzeroa
         write(fileNumber) "idh1(MXCOUNT3)"
         write(fileNumber) savedVars%idh1
         write(fileNumber) "idh2(MXCOUNT3)"
         write(fileNumber) savedVars%idh2
         write(fileNumber) "idh3(MXCOUNT3)"
         write(fileNumber) savedVars%idh3
         write(fileNumber) "idh4(MXCOUNT3)"
         write(fileNumber) savedVars%idh4
         write(fileNumber) "idh5(MXCOUNT3)"
         write(fileNumber) savedVars%idh5
         write(fileNumber) "idl1(MXCOUNT3)"
         write(fileNumber) savedVars%idl1
         write(fileNumber) "idl2(MXCOUNT3)"
         write(fileNumber) savedVars%idl2
         write(fileNumber) "idl3(MXCOUNT3)"
         write(fileNumber) savedVars%idl3
         write(fileNumber) "idl4(MXCOUNT3)"
         write(fileNumber) savedVars%idl4
         write(fileNumber) "idl5(MXCOUNT3)"
         write(fileNumber) savedVars%idl5
         write(fileNumber) "ikdeca(MXCOUNT3)"
         write(fileNumber) savedVars%ikdeca
         write(fileNumber) "ikdecb(MXCOUNT3)"
         write(fileNumber) savedVars%ikdecb
         write(fileNumber) "ikdecc(MXCOUNT3)"
         write(fileNumber) savedVars%ikdecc
         write(fileNumber) "ikdecd(MXCOUNT3)"
         write(fileNumber) savedVars%ikdecd
         write(fileNumber) "ikdece(MXCOUNT3)"
         write(fileNumber) savedVars%ikdece
         write(fileNumber) "kjdeca(MXCOUNT3)"
         write(fileNumber) savedVars%kjdeca
         write(fileNumber) "kjdecb(MXCOUNT3)"
         write(fileNumber) savedVars%kjdecb
         write(fileNumber) "kjdecc(MXCOUNT3)"
         write(fileNumber) savedVars%kjdecc
         write(fileNumber) "kjdecd(MXCOUNT3)"
         write(fileNumber) savedVars%kjdecd
         write(fileNumber) "kjdece(MXCOUNT3)"
         write(fileNumber) savedVars%kjdece
         write(fileNumber) "ijthi(MXGSAER,ICP)"
         write(fileNumber) savedVars%ijthi
         write(fileNumber) "ijtlo(MXGSAER,ICP)"
         write(fileNumber) savedVars%ijtlo
         write(fileNumber) "jarrdiag(MXGSAER,ICP)"
         write(fileNumber) savedVars%jarrdiag
         write(fileNumber) "jhiz1(MXGSAER,ICP)"
         write(fileNumber) savedVars%jhiz1
         write(fileNumber) "jloz1(MXGSAER,ICP)"
         write(fileNumber) savedVars%jloz1
         write(fileNumber) "iarray(ICP)"
         write(fileNumber) savedVars%iarray
         write(fileNumber) "npdhi(ICP)"
         write(fileNumber) savedVars%npdhi
         write(fileNumber) "npdlo(ICP)"
         write(fileNumber) savedVars%npdlo
         write(fileNumber) "iialpd(MXCOUNT2)"
         write(fileNumber) savedVars%iialpd
         write(fileNumber) "ipospd(MXCOUNT2)"
         write(fileNumber) savedVars%ipospd
         write(fileNumber) "nkpdterm(MXCOUNT2)"
         write(fileNumber) savedVars%nkpdterm
         write(fileNumber) "nfrhi(ICP)"
         write(fileNumber) savedVars%nfrhi
         write(fileNumber) "nfrlo(ICP)"
         write(fileNumber) savedVars%nfrlo
         write(fileNumber) "nplhi(ICP)"
         write(fileNumber) savedVars%nplhi
         write(fileNumber) "npllo(ICP)"
         write(fileNumber) savedVars%npllo
         write(fileNumber) "jspcnfr(MXCOUNT4)"
         write(fileNumber) savedVars%jspcnfr
         write(fileNumber) "jspnpl(MXCOUNT4)"
         write(fileNumber) savedVars%jspnpl
         write(fileNumber) "nknfr(MXCOUNT4)"
         write(fileNumber) savedVars%nknfr
         write(fileNumber) "lossra(MXCOUNT4)"
         write(fileNumber) savedVars%lossra
         write(fileNumber) "lossrb(MXCOUNT4)"
         write(fileNumber) savedVars%lossrb
         write(fileNumber) "lossrc(MXCOUNT4)"
         write(fileNumber) savedVars%lossrc
         write(fileNumber) "lossrd(MXCOUNT4)"
         write(fileNumber) savedVars%lossrd
         write(fileNumber) "lossre(MXCOUNT4)"
         write(fileNumber) savedVars%lossre
         write(fileNumber) "nph1(MXCOUNT4)"
         write(fileNumber) savedVars%nph1
         write(fileNumber) "nph2(MXCOUNT4)"
         write(fileNumber) savedVars%nph2
         write(fileNumber) "nph3(MXCOUNT4)"
         write(fileNumber) savedVars%nph3
         write(fileNumber) "nph4(MXCOUNT4)"
         write(fileNumber) savedVars%nph4
         write(fileNumber) "nph5(MXCOUNT4)"
         write(fileNumber) savedVars%nph5
         write(fileNumber) "npl1(MXCOUNT4)"
         write(fileNumber) savedVars%npl1
         write(fileNumber) "npl2(MXCOUNT4)"
         write(fileNumber) savedVars%npl2
         write(fileNumber) "npl3(MXCOUNT4)"
         write(fileNumber) savedVars%npl3
         write(fileNumber) "npl4(MXCOUNT4)"
         write(fileNumber) savedVars%npl4
         write(fileNumber) "npl5(MXCOUNT4)"
         write(fileNumber) savedVars%npl5
         write(fileNumber) "nolosp(ICP)"
         write(fileNumber) savedVars%nolosp
         write(fileNumber) "newfold(NMTRATE*2,ICS)"
         write(fileNumber) savedVars%newfold
         write(fileNumber) "nknlosp(MAXGL3,ICS)"
         write(fileNumber) savedVars%nknlosp
         write(fileNumber) "nknphotrt(IPHOT,ICS)"
         write(fileNumber) savedVars%nknphotrt
         write(fileNumber) "abst2(ICS)"
         write(fileNumber) savedVars%abst2
         write(fileNumber) "errmax(ICS)"
         write(fileNumber) savedVars%errmax
         write(fileNumber) "hmaxday(ICS)"
         write(fileNumber) savedVars%hmaxday
         write(fileNumber) "timeintv(ICS)"
         write(fileNumber) savedVars%timeintv
         write(fileNumber) "abtol(6,ICS)"
         write(fileNumber) savedVars%abtol
         write(fileNumber) "enqq1(MORDER)"
         write(fileNumber) savedVars%enqq1
         write(fileNumber) "enqq2(MORDER)"
         write(fileNumber) savedVars%enqq2
         write(fileNumber) "enqq3(MORDER)"
         write(fileNumber) savedVars%enqq3
         write(fileNumber) "conp15(MORDER)"
         write(fileNumber) savedVars%conp15
         write(fileNumber) "conpst(MORDER)"
         write(fileNumber) savedVars%conpst
         write(fileNumber) "pertst2(MORDER,3)"
         write(fileNumber) savedVars%pertst2
         write(fileNumber) "aset(10,8)"
         write(fileNumber) savedVars%aset
         write(fileNumber) "fracpl(MXCOUNT2)"
         write(fileNumber) savedVars%fracpl
         write(fileNumber) "fracnfr(MXCOUNT4)"
         write(fileNumber) savedVars%fracnfr

         close(fileNumber)

      end subroutine writeSmv2Chem2Exit


