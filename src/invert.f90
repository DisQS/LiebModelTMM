!*****************************************************************************
!
! THERMO CALCULATES THE KINETIC TRANSPORT PROPERTIES 
! OF 2&3 DIMENSIONAL STRIPS USING A RECURSIVE METHOD.
!
! A. MacKinnon, R.A. Roemer, C. Villagonzalo
!
! $Id: invert.f90,v 1.1.1.1 2012/09/07 10:00:09 phsht Exp $
!***************************************************************************

!***************************************************************************
! $Log: invert.f90,v $
! Revision 1.1.1.1  2012/09/07 10:00:09  phsht
! tmseXd
!
! Revision 1.1  2011/05/06 08:13:09  phsht
! 1st installement
!
! Revision 1.2  2010/10/26 09:43:39  phrkaj
! Deleted naive debugging statements, got rid of ILevelflag and IConvflag, deleted old logs
!
! Revision 1.1.1.1  2010/10/22 12:23:38  phsht
! ParaTMM
!
!***************************************************************************

SUBROUTINE INVERT(M,A,B,INFO)  
!
!   Invert an M*M Complex Matrix
!   A: the Matrix (Destroyed)
!   B: the Inverse
   USE MyNUMBERS
   IMPLICIT NONE

   INTEGER :: M, LWORK, INFO, I
   COMPLEX(KIND=RKIND), DIMENSION(1:M,1:M) :: A
   COMPLEX(KIND=RKIND), DIMENSION(1:M,1:M) :: B

   INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
   COMPLEX(KIND=RKIND), DIMENSION(:), ALLOCATABLE :: WORK

   PRINT*,"DGB: INVERT()"

   B = CZERO 
   DO I=1,M
      B(I,I) = CONE
   END DO
   INFO=0

   ALLOCATE(IPIV(M),STAT=INFO)
   IF( INFO.NE.0 ) THEN
      PRINT*,"Invert(): ERR in ALLOCATE(IPIV(M)) statement, M=", M
      STOP
   ENDIF

   CALL ZGETRF(M,M,A,M,IPIV,INFO)
   LWORK = M*M

   ALLOCATE(WORK(LWORK),STAT=INFO)   
   IF( INFO.NE.0 ) THEN
      PRINT*,"Invert(): ERR in ALLOCATE(WORK(LWORK)) statement, LWORK=", LWORK
      RETURN
   ENDIF
   
   CALL ZGETRI(M,A,M,IPIV,WORK,LWORK,INFO)
   IF ( INFO.NE.0 ) THEN
      PRINT *,'Inversion Error: IFAIL=',INFO
      RETURN
   END IF
   DEALLOCATE(IPIV,WORK)
   B = A  
   RETURN
END subroutine INVERT


