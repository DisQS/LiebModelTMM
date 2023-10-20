!********************************************************************
!
! TMSEXD - Transfer matrix method for the Anderson
! model with diagonal disorder in X dimensions
!
!********************************************************************
       
!********************************************************************
!
! $Header: /home/cvs/phsht/tmseXd/src/application_modules.f90,v 1.2 2016/10/03 17:33:32 phsht Exp $
!
!********************************************************************

!**************************************************************************
!$Log: application_modules.f90,v $
!Revision 1.2  2016/10/03 17:33:32  phsht
!1st attempt at correct 1D version
!
!Revision 1.1  2016/10/02 19:13:05  phsht
!now gfortran version
!
!Revision 1.2  2012/09/07 10:44:10  phsht
!removed MagFlux and hence all CMPLX,COMPLEX definitions;
!included a final column in AVG output to indicate convergence (or not)
!
!Revision 1.1  2012/09/07 10:03:27  phsht
!initially forgotten modules file
!
!Revision 1.2  2011/05/31 13:53:43  ccspam
!*** empty log message ***
!
!Revision 1.1  2011/05/06 08:13:09  phsht
!1st installement
!
!Revision 1.1  2010/11/11 11:16:25  phrkaj
!Renamed files to get rid of any mention of SB.
!
!
!**************************************************************************

!--------------------------------------------------------------------
MODULE CConstants
  CHARACTER*18, PARAMETER :: RStr= "$Revision: 1.2 $ "
  CHARACTER*30, PARAMETER :: DStr= "$Date: 2016/10/03 17:33:32 $ "
  CHARACTER*16, PARAMETER :: AStr= "$Author: phsht $ "
END MODULE CConstants

! MAXGamma needs to be equal to MAXWidth, as we need to find ALL
! Lyapunov exponents, so do not change!
MODULE IConstants
  INTEGER, PARAMETER :: MAXWidth= 1000, MAXGamma= MAXWidth, MAXIter=2147483646
  INTEGER, PARAMETER :: MAXKeepFlag= 3, MAXWriteFlag= 4, MAXFluxFlag= 3, MAXRNGFlag=30
  INTEGER, PARAMETER :: MAXSortFlag=1, MAXBCFlag=2
  INTEGER, PARAMETER :: MINDimenFlag=1, MAXDimenFlag=3
  INTEGER, PARAMETER :: MAXFiles= 5, MINIter=3
END MODULE IConstants

!--------------------------------------------------------------------
MODULE IPara
  INTEGER :: ISeed, NOfIter, NOfOrtho, NOfPrint, NOfGamma
  INTEGER :: Width0, Width1, dWidth, IKeepFlag, IWriteFlag, ISortFlag
  INTEGER :: IFluxFlag, IBCFlag, IRNGFlag
  INTEGER :: IDimenFlag
END MODULE IPara

!--------------------------------------------------------------------
MODULE DPara
  USE MyNumbers
  REAL(KIND=RKIND) :: DiagDis0,DiagDis1,dDiagDis
  REAL(KIND=RKIND) :: tW0,tW1,dtW
  REAL(KIND=RKIND) :: Energy0,Energy1,dEnergy
  REAL(KIND=RKIND) :: Kappa, Epsilon
  
  !REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE :: PsiA, PsiB
END MODULE DPara

!--------------------------------------------------------------------
!      Input- and Outputchannels
MODULE IChannels
  INTEGER, PARAMETER :: IChInp= 40, IChOut= 41, IChOutGam= 42, &
       IChOutPsi= 43, IChOutRHO=44,IChOutWrite= 39, IChOutAvgRHO= 45, &
       IChOutAvgRHO1= 46, IChOutAvgRHOL= 47, &
       ICHtmp= 48, IChOutHAV= 49
END MODULE IChannels






