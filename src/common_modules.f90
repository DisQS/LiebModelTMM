!     Last change:  RAR  18 Sep 1999   10:54 am
!********************************************************************
!
! TMSE2D - Transfer matrix method for the Anderson
! model with diagonal disorder in two dimensions
!
!********************************************************************
       
!********************************************************************
!
! $Header: /home/cvs/phsht/tmseXd/src/common_modules.f90,v 1.1 2016/10/02 19:13:05 phsht Exp $
!
!********************************************************************

!**************************************************************************
!$Log: common_modules.f90,v $
!Revision 1.1  2016/10/02 19:13:05  phsht
!now gfortran version
!
!Revision 1.1.1.1  2012/09/07 10:00:09  phsht
!tmseXd
!
!Revision 1.1  2011/05/06 08:13:09  phsht
!1st installement
!
!Revision 1.4  2010/10/26 14:28:07  phrkaj
!Replaced !! with !, fixed the OpenOutputGamma of negative energy in inout.f90
!
!Revision 1.3  2010/10/25 15:41:33  phsht
!small changes to remove a "malloc/glibc" error
!
!Revision 1.2  2010/10/24 14:23:59  phrkaj
!Added DZERO parameter
!
!Revision 1.1.1.1  2010/10/22 12:23:38  phsht
!ParaTMM
!
!**************************************************************************

MODULE MyNumbers     
  IMPLICIT NONE

  INTEGER, PARAMETER :: IKIND = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(15,307)
  INTEGER, PARAMETER :: CKIND = SELECTED_REAL_KIND(15,307)

  REAL(KIND=RKIND) :: PI, TWOPI, ONEPLS, ONEMNS

  REAL(KIND=RKIND), PARAMETER :: ZERO = 0.0, ONE = 1.0 ,TWO = 2.0, THREE = 3.0, FOUR = 4.0, t=1.0
  COMPLEX(KIND=RKIND), PARAMETER :: CZERO = (0.0d0,0.0d0), CONE = (1.0d0,0.0d0), &
       CIMAGONE= (0.0d0,1.0d0)            

  REAL (KIND=RKIND), PARAMETER :: HALF = 0.5D0, QUARTER = 0.25D0, EIGHTH = 0.125D0

  REAL(KIND=RKIND) :: TINY= 1.0D-9
  
CONTAINS
  SUBROUTINE INIT_NUMBERS
    PI     = 4.0D0* ATAN(1.0D0)
    TWOPI  = 8.0D0* ATAN(1.0D0)
    ONEMNS = SQRT(EPSILON(ONEMNS))
    ONEPLS = ONE + ONEMNS
    ONEMNS = ONE - ONEMNS
  END SUBROUTINE INIT_NUMBERS

  FUNCTION ARG(X,Y)
    
    REAL(KIND=RKIND) ARG, X, Y
    
    IF( X > 0. ) THEN 
       ARG= ATAN(Y/X)
    ELSE IF ( (X == 0.) .and. (Y > 0. )) THEN 
       ARG = PI/2.0D0
    ELSE IF ( (X == 0.) .and. (Y < 0. )) THEN 
       ARG = -PI/2.0D0
    ELSE IF ( (X < 0. ) .and. (Y >= 0.)) THEN 
       ARG = PI + ATAN(Y/X)
    ELSE IF ( (X < 0. ) .and. (Y < 0. )) THEN 
       ARG = -PI + ATAN(Y/X)
    ENDIF
    
    RETURN
  END FUNCTION ARG

END MODULE MyNumbers



