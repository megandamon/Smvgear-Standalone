!=======================================================================
!
! $Id: setkin_smv2par.h,v 1.8 2011-08-09 22:52:11 mrdamon Exp $
!
! FILE
!   setkin_smv2par.h
!   12 JUN 02 - PSC
!
! DESCRIPTION
!   This include file sets symbolic constants for the photochemical
!   mechanism.
!
!  Chemistry input file:    10/2006
!  Reaction dictionary:     GMI_Combo_rxns_124species_SO2_JPL06.db
!  Setkin files generated:  Thu Nov 19 21:21:30 2009
!
!========1=========2=========3=========4=========5=========6=========7==

      integer &
     &  SK_IGAS &
     & ,SK_IPHOT &
     & ,SK_ITHERM &
     & ,SK_NACT

      parameter (SK_IGAS   = 121)
      parameter (SK_IPHOT  =  81)
      parameter (SK_ITHERM = 321)
      parameter (SK_NACT   = 117)

!                                  --^--

