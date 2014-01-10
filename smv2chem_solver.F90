
program doSmv2Solver

      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved
      use Physproc_mod
      use timing_mod

      implicit none

#     include "smv2chem_par.h"
#     include 'mpif.h'


      integer :: commuWorld
      logical, save :: first = .true.
      type(t_Smv2Saved) :: savedVars
      type(Physproc_type) :: physprocVars


      integer, allocatable :: jreorder(:)
      integer, allocatable :: lreorder(:)
      real*8, allocatable  :: errmx2  (:)

      character(len=100) :: smv2Chem1Entry
      character(len=100) :: smv2Chem1Exit
      character(len=100) :: smv2Chem2Entry
      character(len=100) :: smv2Chem2Exit
      character(len=100) :: physProcEntry
      character(len=100) :: physProcExit
      character(:), allocatable :: rankString
      character(:), allocatable :: zeroString

      integer :: rankSize
      integer, parameter :: MAX_RANK_SIZE = 4

      integer mpiError, rc, numTasks
      !     ----------------
      !     Setup area
      !     ----------------

      call timingInit

      ! Initialize the MPI library:
      call MPI_INIT(mpiError)
      if (mpiError .ne. MPI_SUCCESS) then
         print *,'Error starting MPI program. Terminating.'
         call MPI_ABORT(MPI_COMM_WORLD, rc, mpiError)
      end if

      ! Get the number of processors this job is using:
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, mpiError)
      print*, "numTasks: ", numTasks
      if (numTasks .gt. 1) then
         print*, "This driver only accepts numTasks = 1"
         call MPI_FINALIZE(mpiError)
         stop
      endif

      call get_command_argument(2, length=rankSize)
      if (rankSize .gt. 3) then
         print*, ""
         print*, "Usage: ./Do_Smv2_Solver.exe -r processID"
         print*, ""
         print*, "e.g. ./Do_Smv2_Solver.exe -r 79 for process # 79"
         print*, ""
         print*, "Rules:"
         print*, "   - Processor decomposition must be less than 1,000"
         print*, "   - Process IDs must be expressed in fewer than 4 characters"
         print*, ""
         stop
      end if

      allocate(character(rankSize)::rankString)
      allocate(character(MAX_RANK_SIZE-rankSize)::zeroString)
      call get_command_argument(2, value=rankString)

      physprocVars%prDiag = .true.
      if (physprocVars%prDiag) then
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


      !     ----------------
      !    Prepare input data
      !     ----------------
      call readPhysprocVars(physprocVars, physProcEntry)
      call initializeSavedVars(savedVars)
      call readSmv1Vars(savedVars,smv2Chem1Entry)
      call readSmv2Vars(savedVars,smv2Chem2Entry)

      if (first) then
         first = .false.
         Allocate (savedVars%csuma(physprocVars%itloop))
         Allocate (savedVars%csumc(physprocVars%itloop))
         savedVars%csuma = 0.0d0; savedVars%csumc = 0.0d0
      end if

      allocate (jreorder(physprocVars%itloop))
      allocate (lreorder(physprocVars%itloop))
      allocate (errmx2(physprocVars%itloop))
      jreorder(:) = 0; lreorder(:) = 0
      errmx2  (:) = 0.0d0

      commuWorld = MPI_COMM_WORLD

      print*, "Calling Physproc"
      call Physproc (savedVars, physprocVars%doQqjkInchem, physprocVars%doSurfEmissInChem, &
               physprocVars%prQqjk, physprocVars%prSmv2, physprocVars%numLat, physprocVars%numLong, physprocVars%numVert, savedVars%ifreord,&
               savedVars%imgas, savedVars%initrogen, savedVars%ioxygen,&
               physprocVars%itloop, savedVars%kuloop, savedVars%lunsmv, &
               savedVars%ncs, savedVars%fracdec, savedVars%hmaxnit, &
               physprocVars%prNcPeriod, physprocVars%timeStep, physprocVars%doCellChem, savedVars%jphotrat, &
               savedVars%nrates, savedVars%ntloopncs, savedVars%ntspec,&
               savedVars%inewold, savedVars%npphotrat, physprocVars%thermalRateConstants, physprocVars%photolysisRateConstants, &
               physprocVars%surfaceEmissions, jreorder, lreorder, savedVars%csuma,  &
               savedVars%csumc, errmx2, physprocVars%speciesConst, physprocVars%yda, physprocVars%qqkda, physprocVars%qqjda, &
               physprocVars%i1, physprocVars%i2, physprocVars%ju1, physprocVars%j2, physprocVars%k1, physprocVars%k2, physprocVars%numQks, &
               physprocVars%numQjs, physprocVars%numActive, &
               commuWorld)
      print*, "Returned from Physproc"

      !     ----------------
      !    Write output data and deallocate variables
      !     ----------------
      call writeSmv2Chem1Exit(smv2Chem1Exit,savedVars)
      call writeSmv2Chem2Exit(smv2Chem2Exit,savedVars)
      call writePhysprocVars(physprocVars, physProcExit)

      deallocate (savedVars%csuma)
      deallocate (savedVars%csumc)
      deallocate (jreorder)
      deallocate (lreorder)
      deallocate (errmx2)

      ! Tell the MPI library to release all resources it is using:
      call MPI_FINALIZE(mpiError)

      print*, "Exiting doSmv2Solver"

end program doSmv2Solver

subroutine readSmv2Vars(savedVars,smv2Chem2Entry)

      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

      IMPLICIT NONE
#     include "smv2chem_par.h"

      type(t_Smv2Saved), intent(inOut) :: savedVars
      character(len=100), intent(in) :: smv2Chem2Entry

      integer :: fileNumber

      print*, "reading: ", trim(smv2Chem2Entry)

!read smv2chem2 on entry to physproc to be used for standAlone code
      open(file=trim(smv2Chem2Entry),newunit=fileNumber,form="formatted")
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ioner
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nallrat
      read(fileNumber,*)
      read(fileNumber,*) savedVars%inorep
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ithrr
      read(fileNumber,*)
      read(fileNumber,*) savedVars%itwor
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nm3bod
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nmair
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nmn2
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nmo2
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nmoth

      read(fileNumber,*)
      read(fileNumber,*) savedVars%ntrates
      read(fileNumber,*)
      read(fileNumber,*) savedVars%mappl
      read(fileNumber,*)
      read(fileNumber,*) savedVars%lgasbino
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nreacoth
      read(fileNumber,*)
      read(fileNumber,*) savedVars%lgas3bod
      read(fileNumber,*)
      read(fileNumber,*) savedVars%losinacp
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nreac3b
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nreacair
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nreacn2
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nreaco2
      read(fileNumber,*)
      read(fileNumber,*) savedVars%jphotnk
      read(fileNumber,*)
      read(fileNumber,*) savedVars%noldfnew
      read(fileNumber,*)
      read(fileNumber,*) savedVars%irm2
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ischang
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kzthi
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kztlo
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ikztot
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kbh1
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kbh2
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kbh3
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kbh4
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kbh5
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kbl1
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kbl2
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kbl3
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kbl4
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kbl5
      read(fileNumber,*)
      read(fileNumber,*) savedVars%mbh1
      read(fileNumber,*)
      read(fileNumber,*) savedVars%mbh2
      read(fileNumber,*)
      read(fileNumber,*) savedVars%mbh3
      read(fileNumber,*)
      read(fileNumber,*) savedVars%mbh4
      read(fileNumber,*)
      read(fileNumber,*) savedVars%mbh5
      read(fileNumber,*)
      read(fileNumber,*) savedVars%mbl1
      read(fileNumber,*)
      read(fileNumber,*) savedVars%mbl2
      read(fileNumber,*)
      read(fileNumber,*) savedVars%mbl3
      read(fileNumber,*)
      read(fileNumber,*) savedVars%mbl4
      read(fileNumber,*)
      read(fileNumber,*) savedVars%mbl5
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kzeroa
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kzerob
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kzeroc
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kzerod
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kzeroe
      read(fileNumber,*)
      read(fileNumber,*) savedVars%mzeroa
      read(fileNumber,*)
      read(fileNumber,*) savedVars%mzerob
      read(fileNumber,*)
      read(fileNumber,*) savedVars%mzeroc
      read(fileNumber,*)
      read(fileNumber,*) savedVars%mzerod
      read(fileNumber,*)
      read(fileNumber,*) savedVars%mzeroe
      read(fileNumber,*)
      read(fileNumber,*) savedVars%imztot
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ijval
      read(fileNumber,*)
      read(fileNumber,*) savedVars%jzeroa
      read(fileNumber,*)
      read(fileNumber,*) savedVars%idh1
      read(fileNumber,*)
      read(fileNumber,*) savedVars%idh2
      read(fileNumber,*)
      read(fileNumber,*) savedVars%idh3
      read(fileNumber,*)
      read(fileNumber,*) savedVars%idh4
      read(fileNumber,*)
      read(fileNumber,*) savedVars%idh5
      read(fileNumber,*)
      read(fileNumber,*) savedVars%idl1
      read(fileNumber,*)
      read(fileNumber,*) savedVars%idl2
      read(fileNumber,*)
      read(fileNumber,*) savedVars%idl3
      read(fileNumber,*)
      read(fileNumber,*) savedVars%idl4
      read(fileNumber,*)
      read(fileNumber,*) savedVars%idl5
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ikdeca
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ikdecb
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ikdecc
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ikdecd
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ikdece
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kjdeca
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kjdecb
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kjdecc
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kjdecd
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kjdece
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ijthi
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ijtlo
      read(fileNumber,*)
      read(fileNumber,*) savedVars%jarrdiag
      read(fileNumber,*)
      read(fileNumber,*) savedVars%jhiz1
      read(fileNumber,*)
      read(fileNumber,*) savedVars%jloz1
      read(fileNumber,*)
      read(fileNumber,*) savedVars%iarray
      read(fileNumber,*)
      read(fileNumber,*) savedVars%npdhi
      read(fileNumber,*)
      read(fileNumber,*) savedVars%npdlo
      read(fileNumber,*)
      read(fileNumber,*) savedVars%iialpd
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ipospd
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nkpdterm
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nfrhi
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nfrlo
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nplhi
      read(fileNumber,*)
      read(fileNumber,*) savedVars%npllo
      read(fileNumber,*)
      read(fileNumber,*) savedVars%jspcnfr
      read(fileNumber,*)
      read(fileNumber,*) savedVars%jspnpl
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nknfr
      read(fileNumber,*)
      read(fileNumber,*) savedVars%lossra
      read(fileNumber,*)
      read(fileNumber,*) savedVars%lossrb
      read(fileNumber,*)
      read(fileNumber,*) savedVars%lossrc
      read(fileNumber,*)
      read(fileNumber,*) savedVars%lossrd
      read(fileNumber,*)
      read(fileNumber,*) savedVars%lossre
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nph1
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nph2
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nph3
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nph4
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nph5
      read(fileNumber,*)
      read(fileNumber,*) savedVars%npl1
      read(fileNumber,*)
      read(fileNumber,*) savedVars%npl2
      read(fileNumber,*)
      read(fileNumber,*) savedVars%npl3
      read(fileNumber,*)
      read(fileNumber,*) savedVars%npl4
      read(fileNumber,*)
      read(fileNumber,*) savedVars%npl5
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nolosp
      read(fileNumber,*)
      read(fileNumber,*) savedVars%newfold
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nknlosp
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nknphotrt
      read(fileNumber,*)
      read(fileNumber,*) savedVars%abst2
      read(fileNumber,*)
      read(fileNumber,*) savedVars%errmax
      read(fileNumber,*)
      read(fileNumber,*) savedVars%hmaxday
      read(fileNumber,*)
      read(fileNumber,*) savedVars%timeintv
      read(fileNumber,*)
      read(fileNumber,*) savedVars%abtol
      read(fileNumber,*)



      allocate(savedVars%enqq1 (MORDER))
      allocate(savedVars%enqq2 (MORDER))
      allocate(savedVars%enqq3 (MORDER))
      allocate(savedVars%conp15(MORDER))
      allocate(savedVars%conpst(MORDER))
      allocate(savedVars%pertst2(MORDER, 3))

      read(fileNumber,*) savedVars%enqq1
      read(fileNumber,*)
      read(fileNumber,*) savedVars%enqq2
      read(fileNumber,*)
      read(fileNumber,*) savedVars%enqq3
      read(fileNumber,*)
      read(fileNumber,*) savedVars%conp15
      read(fileNumber,*)
      read(fileNumber,*) savedVars%conpst
      read(fileNumber,*)
      read(fileNumber,*) savedVars%pertst2

      allocate(savedVars%aset(10, 8))
      !allocate(savedVars%fracpl (MXCOUNT2))
      !allocate(savedVars%fracnfr(MXCOUNT4))

      read(fileNumber,*)
      read(fileNumber,*) savedVars%aset
      read(fileNumber,*)

      read(fileNumber,*) savedVars%fracpl(:)
      read(fileNumber,*)
      read(fileNumber,*) savedVars%fracnfr(:)

      close(fileNumber)
   end subroutine readSmv2Vars

 !read smv2chem1 on entry to physproc to be used for standAlone code
  subroutine readSmv1Vars(savedVars,smv2Chem1Entry)

      use GmiSolver_SavedVariables_mod, only : t_Smv2Saved

      IMPLICIT NONE

      type(t_Smv2Saved), intent(inOut) :: savedVars
      character(len=100), intent(in) :: smv2Chem1Entry
      integer :: fileNumber

      print*, "reading: ", trim(smv2Chem1Entry)

      open(file=trim(smv2Chem1Entry),newunit=fileNumber,form="formatted")
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ifreord
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ih2o
      read(fileNumber,*)
      read(fileNumber,*) savedVars%imgas
      read(fileNumber,*)
      read(fileNumber,*) savedVars%initrogen
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ioxygen
      read(fileNumber,*)
      read(fileNumber,*) savedVars%kuloop
      read(fileNumber,*)
      read(fileNumber,*) savedVars%lunsmv
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ncs
      read(fileNumber,*)
      read(fileNumber,*) savedVars%jphotrat
      read(fileNumber,*)
      read(fileNumber,*) savedVars%nrates
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ntloopncs
      read(fileNumber,*)
      read(fileNumber,*) savedVars%ntspec
      read(fileNumber,*)
      read(fileNumber,*) savedVars%inewold
      read(fileNumber,*)
      read(fileNumber,*) savedVars%npphotrat
      read(fileNumber,*)
      read(fileNumber,*) savedVars%fracdec
      read(fileNumber,*)
      read(fileNumber,*) savedVars%hmaxnit
      close(fileNumber)

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

         deallocate(savedVars%enqq1)
         deallocate(savedVars%enqq2)
         deallocate(savedVars%enqq3)
         deallocate(savedVars%conp15)
         deallocate(savedVars%conpst)

         deallocate(savedVars%aset)

      end subroutine writeSmv2Chem2Exit


