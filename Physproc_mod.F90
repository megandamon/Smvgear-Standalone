module Physproc_mod

   implicit none

    public :: readPhysprocVars

   type Physproc_type

      real*8, allocatable :: qjGmi(:, :, :, :)
      real*8, allocatable :: qkGmi(:, :, :, :)
      real*8, allocatable :: qqjda(:,:,:,:)
      real*8, allocatable:: qqkda(:,:,:,:)
      real*8, allocatable :: yda (:,:,:,:)
      real*8, allocatable :: thermalRateConstants(:, :)
      real*8, allocatable :: photolysisRateConstants(:, :)
      real*8, allocatable :: surfaceEmissions(:, :)
      real*8, allocatable :: speciesConst(:, :)
      real*8  :: prNcPeriod
      real*8  :: timeStep

      logical, allocatable  :: doCellChem(:)

      logical  :: doQqjkInchem
      logical  :: doSurfEmissInChem
      logical  :: prQqjk
      logical  :: prSmv2
      logical  :: prDiag

      integer  :: localProc
      integer  :: numLat, numLong, numVert
      integer  :: itloop
      integer :: numQjo, numQks, numQjs, numActive
      integer :: i1, i2, ju1, j2, k1, k2


   end type Physproc_type

   contains

   subroutine readPhysprocVars (this, physProcEntry)
#     include "smv2chem_par.h"

      type(Physproc_type), intent(inout) :: this
      character(len=100), intent(in) :: physProcEntry
      integer :: i1, i2, ju1, j2, k1, k2
      integer :: fileNumber

      print*, "reading: ", trim(physProcEntry)
      open(file=trim(physProcEntry),newunit=fileNumber,form="formatted")
      read(fileNumber,*)
      read(fileNumber,*) i1, i2, ju1, j2, k1, k2
      read(fileNumber,*)
      read(fileNumber,*) this%numQjo, this%numQks, this%numQjs, this%numActive

      allocate(this%qjGmi(i1:i2, ju1:j2, k1:k2, this%numQjo))
      allocate(this%qkGmi(i1:i2, ju1:j2, k1:k2, this%numQks))
      allocate(this%qqjda(i1:i2, ju1:j2, k1:k2, this%numQjs))
      allocate(this%qqkda(i1:i2, ju1:j2, k1:k2, this%numQks))
      allocate(this%yda  (i1:i2, ju1:j2, k1:k2, this%numActive))
      read(fileNumber,*)
      read(fileNumber,*) this%qjGmi(i1:i2, ju1:j2, k1:k2, 1:this%numQjo)
      read(fileNumber,*)
      read(fileNumber,*) this%qkGmi(i1:i2, ju1:j2, k1:k2, 1:this%numQks)
      read(fileNumber,*)
      read(fileNumber,*) this%qqjda(i1:i2, ju1:j2, k1:k2, 1:this%numQjs)
      read(fileNumber,*)
      read(fileNumber,*) this%qqkda(i1:i2, ju1:j2, k1:k2, 1:this%numQks)
      read(fileNumber,*)
      read(fileNumber,*) this%yda  (i1:i2, ju1:j2, k1:k2, 1:this%numActive)
      read(fileNumber,*)
      read(fileNumber,*) this%doQqjkInchem
      read(fileNumber,*)
      read(fileNumber,*) this%doSurfEmissInChem
      read(fileNumber,*)
      read(fileNumber,*) this%prDiag
      read(fileNumber,*)
      read(fileNumber,*) this%prQqjk
      read(fileNumber,*)
      read(fileNumber,*) this%prSmv2
      read(fileNumber,*)
      read(fileNumber,*) this%localProc
      read(fileNumber,*)
      read(fileNumber,*) this%numLat, this%numLong, this%numVert
      read(fileNumber,*)
      read(fileNumber,*) this%itloop
      read(fileNumber,*)
      read(fileNumber,*) this%prNcPeriod
      read(fileNumber,*)
      read(fileNumber,*) this%timeStep

      allocate(this%doCellChem(this%itloop))
      allocate(this%thermalRateConstants(this%itloop, ITHERM))
      allocate(this%photolysisRateConstants(this%itloop, IPHOT))
      allocate(this%surfaceEmissions(this%numLat*this%numLong, IGAS))
      allocate(this%speciesConst(this%itloop, IGAS))

      read(fileNumber,*)
      read(fileNumber,*) this%doCellChem(1:this%itloop)
      read(fileNumber,*)
      read(fileNumber,*) this%thermalRateConstants(:,:)
      read(fileNumber,*)
      read(fileNumber,*) this%photolysisRateConstants(1:this%itloop, 1:IPHOT)
      read(fileNumber,*)
      read(fileNumber,*) this%surfaceEmissions(1:this%numLat*this%numLong, 1:IGAS)
      read(fileNumber,*)
      read(fileNumber,*) this%speciesConst(1:this%itloop, 1:IGAS)
      close(fileNumber)

      this%i1 = i1
      this%i2 = i2
      this%ju1 = ju1
      this%j2 = j2
      this%k1 = k1
      this%k2 = k2

   end subroutine     readPhysprocVars

   subroutine     writePhysprocVars (this, physProcExit)

#     include "smv2chem_par.h"

      type(Physproc_type), intent(inout) :: this
      character(len=100), intent(in) :: physProcExit
      integer :: i1, i2, ju1, j2, k1, k2

      integer :: fileNumber

      i1 = this%i1
      i2 = this%i2
      ju1 = this%ju1
      j2 = this%j2
      k1 = this%k1
      k2 = this%k2

      print*, "Writing to: ", physProcExit
      open(file=trim(physProcExit),newunit=fileNumber,status="replace",form="unformatted")

      write(fileNumber) "i1, i2, ju1, j2, k1, k2"
      write(fileNumber) i1, i2, ju1, j2, k1, k2
      write(fileNumber) "num_qjo, num_qks, num_qjs, num_active"
      write(fileNumber) this%numQjo, this%numQks, this%numQjs, this%numActive
      write(fileNumber) "qjgmi(i1:i2, ju1:j2, k1:k2, num_qjo)"
      write(fileNumber) this%qjGmi(i1:i2, ju1:j2, k1:k2, 1:this%numQjo)
      write(fileNumber) "qkgmi(i1:i2, ju1:j2, k1:k2, num_qks)"
      write(fileNumber) this%qkGmi(:,:,:,:)
      write(fileNumber) "qqjda(i1:i2, ju1:j2, k1:k2, num_qjs)"
      write(fileNumber) this%qqjda(i1:i2, ju1:j2, k1:k2, 1:this%numQjs)
      write(fileNumber) "qqkda(i1:i2, ju1:j2, k1:k2, num_qks)"
      write(fileNumber) this%qqkda(:,:,:,:)
      write(fileNumber) "yda  (i1:i2, ju1:j2, k1:k2, num_active)"
      write(fileNumber) this%yda  (i1:i2, ju1:j2, k1:k2, 1:this%numActive)
      write(fileNumber) "do_qqjk_inchem"
      write(fileNumber) this%doQqjkInchem
      write(fileNumber) "do_semiss_inchem"
      write(fileNumber) this%doSurfEmissInChem
      write(fileNumber) "pr_diag"
      write(fileNumber) this%prDiag
      write(fileNumber) "pr_qqjk"
      write(fileNumber) this%prQqjk
      write(fileNumber) "pr_smv2"
      write(fileNumber) this%prSmv2
      write(fileNumber) "loc_proc"
      write(fileNumber) this%localProc
      write(fileNumber) "ilat, ilong, ivert"
      write(fileNumber) this%numLat, this%numLong, this%numVert
      write(fileNumber) "itloop"
      write(fileNumber) this%itloop
      write(fileNumber) "pr_nc_period"
      write(fileNumber) this%prNcPeriod
      write(fileNumber) "tdt"
      write(fileNumber) this%timeStep
      write(fileNumber) "do_cell_chem(itloop)"
      write(fileNumber) this%doCellChem(1:this%itloop)
      write(fileNumber) "arate(itloop, ITHERM)"
      write(fileNumber) this%thermalRateConstants(:,:)
      write(fileNumber) "prate(itloop, IPHOT)"
      write(fileNumber) this%photolysisRateConstants(1:this%itloop, 1:IPHOT)
      write(fileNumber) "yemis(ilat*ilong, IGAS)"
      write(fileNumber) this%surfaceEmissions(1:this%numLat*this%numLong, 1:IGAS)
      write(fileNumber) "cx(itloop, IGAS)"
      write(fileNumber) this%speciesConst(1:this%itloop, 1:IGAS)
      close(fileNumber)

      deallocate(this%qjGmi)
      deallocate(this%qkGmi)
      deallocate(this%qqjda)
      deallocate(this%qqkda)
      deallocate(this%yda)


      deallocate(this%doCellChem)
      deallocate(this%thermalRateConstants)
      deallocate(this%photolysisRateConstants)
      deallocate(this%surfaceEmissions)
      deallocate(this%speciesConst)

   end subroutine writePhysprocVars
end module Physproc_mod
