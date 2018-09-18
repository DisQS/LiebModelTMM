!********************************************************************
!
! TMSE2DCOE - Transfer matrix method for the Anderson
! model with diagonal disorder in two dimensions
!
!********************************************************************

!********************************************************************
!
! $Header: /home/cvs/phsht/tmseXd/src/eigen.f90,v 1.1.1.1 2012/09/07 10:00:09 phsht Exp $
!
!********************************************************************

!**************************************************************************
!
! $Log: eigen.f90,v $
! Revision 1.1.1.1  2012/09/07 10:00:09  phsht
! tmseXd
!
! Revision 1.1  2011/05/06 08:13:09  phsht
! 1st installement
!
! Revision 1.3  2010/10/26 14:28:07  phrkaj
! Replaced !! with !, fixed the OpenOutputGamma of negative energy in inout.f90
!
! Revision 1.2  2010/10/26 09:43:39  phrkaj
! Deleted naive debugging statements, got rid of ILevelflag and IConvflag, deleted old logs
!
! Revision 1.1.1.1  2010/10/22 12:23:38  phsht
! ParaTMM
!
!**************************************************************************

SUBROUTINE EigenSpectrum(isize, matU, evals, IErr)

  USE MyNumbers
  
  INTEGER(KIND=IKIND) isize, IErr
  !REAL(KIND=RKIND) matA(isize,isize), matB(isize,isize)

  COMPLEX(KIND=RKIND) matU(isize,isize), matUU(isize,isize), &
                      invmatU(isize,isize), evals(isize)

  INTEGER(KIND=IKIND) LWORK, LRWORK, LIWORK
  COMPLEX(KIND=RKIND), DIMENSION(:), ALLOCATABLE :: WORK
  REAL(KIND=RKIND),    DIMENSION(:), ALLOCATABLE :: RWORK
  REAL(KIND=IKIND),    DIMENSION(:), ALLOCATABLE :: IWORK
  !EXTERNAL ZGEEV

  !PRINT*,"DBG: EigenSpectrum()"

  IErr=0
  
  ! ------------------------------------------------
  ! find optimum size of arrays
  ! ------------------------------------------------
  LWORK=1
  ALLOCATE(WORK(LWORK), STAT = IErr)
  ALLOCATE(RWORK(2*isize), STAT = IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"EigenSpectrum: error in ALLOCATE() for work arrays (query stage)"
     RETURN
  ENDIF

  LWORK=-1
  CALL ZGEEV('N','N', isize, matU, isize, evals, 0,1, 0,1, WORK, LWORK, RWORK, IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"EigenSpectrum: error in ZGEEV determining work arrays"
     RETURN
  ENDIF

  LWORK = INT(WORK(1))
  !PRINT*,"DBG:   LWORK=", LWORK

  ! ------------------------------------------------
  ! ALLOCATE necessary memory
  ! ------------------------------------------------
  DEALLOCATE(WORK)
  ALLOCATE(WORK(LWORK), STAT = IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"EigenSpectrum: error in ALLOCATE() for work arrays (final stage)"
     RETURN
  ENDIF

  ! ------------------------------------------------
  ! do the actual call to get the spectrum
  ! ------------------------------------------------
  CALL ZGEEV('N','N', isize, matU, isize, evals, 0,1, 0,1, WORK, LWORK, RWORK, IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"EigenSpectrum: error ", IErr, " in ZGEEV"
     RETURN
  ENDIF
  DEALLOCATE(WORK,RWORK)

  !PRINT*,"DBG: evals=", evals

  RETURN

END SUBROUTINE EigenSpectrum
