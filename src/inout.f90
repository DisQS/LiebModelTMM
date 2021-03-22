! ********************************************************************
!       
! TMSEXD - Transfer matrix method for the Anderson
! model with diagonal disorder in X dimensions
!
! ********************************************************************
      
       
! ********************************************************************
!       
! $Header: /home/cvs/phsht/tmseXd/src/inout.f90,v 1.4 2016/10/03 20:02:31 phsht Exp $
!
! ********************************************************************

! **************************************************************************
! $Log: inout.f90,v $
! Revision 1.4  2016/10/03 20:02:31  phsht
! corrected the output filename handling
!
! Revision 1.3  2016/10/03 19:32:32  phsht
! new filenames for .raw files
!
! Revision 1.2  2012/09/07 10:44:10  phsht
! removed MagFlux and hence all CMPLX,COMPLEX definitions;
! included a final column in AVG output to indicate convergence (or not)
!
! Revision 1.1.1.1  2012/09/07 10:00:09  phsht
! tmseXd
!
! Revision 1.2  2011/05/31 13:53:43  ccspam
! *** empty log message ***
!
! Revision 1.1  2011/05/06 08:13:09  phsht
! 1st installement
!
! Revision 1.11  2010/11/15 15:30:00  phrkaj
! Added IDimenFlag to openoutputavg so that .raw file shows IDimenFlag
!
! Revision 1.10  2010/11/11 11:16:25  phrkaj
! Renamed files to get rid of any mention of SB.
!
! Revision 1.9  2010/11/03 10:06:18  phrkaj
! Uncommented wave function stuff, reverted back to original .raw file formatting
!
! Revision 1.8  2010/10/27 16:43:51  phrkaj
! Modified the output of the .raw file so that only the data is printed
!
! Revision 1.7  2010/10/26 14:28:07  phrkaj
! Replaced !! with !, fixed the OpenOutputGamma of negative energy in inout.f90
!
! Revision 1.6  2010/10/26 12:47:59  phrkaj
! Got rid of more SB, Conv and Level stuff
!
! Revision 1.5  2010/10/26 09:43:39  phrkaj
! Deleted naive debugging statements, got rid of ILevelflag and IConvflag, deleted old logs
!
! Revision 1.4  2010/10/25 15:41:33  phsht
! small changes to remove a "malloc/glibc" error
!
! Revision 1.3  2010/10/25 14:59:09  phrkaj
! *** empty log message ***
!
! Revision 1.2  2010/10/25 11:09:24  phrkaj
! Deleted SB stuff
!
! Revision 1.1.1.1  2010/10/22 12:23:38  phsht
! ParaTMM
!
! **************************************************************************

! --------------------------------------------------------------------
!	Input:
!
!	IErr	error code

SUBROUTINE Input( IErr )

  USE MyNumbers
  
  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara
  
  USE IChannels
  
  INTEGER IErr, ILine
  REAL(KIND=RKIND) ROfIter
  
  !	PRINT*,"DBG: Input()"
  
  IErr = 0
  ILine= 0
  
  OPEN(UNIT= IChInp, ERR= 120, FILE= "tmseLMxD.inp",&
       STATUS= 'OLD')
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) ISeed
  !PRINT*,"ISeed        = ",ISeed
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) ROfIter
  NOfIter= NINT(MIN(ROfIter,DBLE(MAXIter)))
  ! PRINT*,"ROfIter      = ", ROfIter
  ! PRINT*,"NOfIter      = ", NOfIter
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) NOfOrtho
  !PRINT*,"NOfOrtho     = ", NOfOrtho
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) NOfPrint
  !PRINT*,"NOfPrint     = ", NOfPrint
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) NOfGamma
  !PRINT*,"NOfGamma     = ", NOfGamma
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IDimenFlag
  PRINT*,"IDimenFlag   = ", IDimenFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IBCFlag
  !PRINT*,"IBCFlag      = ", IBCFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IRNGFlag
  !PRINT*,"IRNGFlag      = ", IRNGFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IKeepFlag
  !PRINT*,"IKeepFlag    = ", IKeepFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IWriteFlag
  !PRINT*,"IWriteFlag   = ", IWriteFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) ISortFlag
  !PRINT*,"ISortFlag    = ", ISortFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IFluxFlag
  !PRINT*,"IFluxFlag    = ", IFluxFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) Width0
  !PRINT*,"Width0       = ",Width0
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) Width1
  !PRINT*,"Width1       = ", Width1
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) dWidth
  !PRINT*,"dWidth       = ", dWidth
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) DiagDis0
  !PRINT*,"DiagDis0     = ", DiagDis0
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) DiagDis1
  !PRINT*,"DiagDis1     = ", DiagDis1
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) dDiagDis
  !PRINT*,"dDiagDis     = ", dDiagDis
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) Energy0
  !PRINT*,"Energy0      = ", Energy0
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) Energy1
  !PRINT*,"Energy1      = ", Energy1
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) dEnergy
  !PRINT*,"dEnergy      = ", dEnergy

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) Kappa
  !PRINT*,"Kappa        = ", Kappa
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) Epsilon
  !PRINT*,"Epsilon      = ", Epsilon

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) SmallDisOnRim
  !PRINT*,"SmallDisOnRim      = ", SmallDisOnRim
  
10 FORMAT(16X,I15.1)
  ! 10	FORMAT("IMAXIteration= ",I15.1)
15 FORMAT(16X,F18.9)
  ! 15     FORMAT("IMAXIteration= ",F18.9)
  
  CLOSE(IChInp, ERR=130)
  
  ! check the parameters for validity
  
  IF( ROfIter.GT.MAXIter ) THEN
     PRINT*,"Input(): NOfIter=",ROfIter,&
          " > MAXIter=", MAXIter
     PRINT*,"Input(): NOfIter set to", NOfIter
  ENDIF
  
  IF( NOfIter.LE.0 ) THEN
     PRINT*,"Input(): NOfIter <= 0"
     IErr= 1
  ENDIF
  
  IF( NOfGamma.GT.MAXGamma ) THEN
     PRINT*,"Input(): NOfGamma > MAXGamma (=",MAXGamma,")"
     IErr= 1
  ENDIF
  
  IF( Width0.LE.0 ) THEN
     PRINT*,"Input(): Width0 <= 0"
     IErr= 1
  ENDIF
  
  IF( Width0.GT.MAXWidth ) THEN
     PRINT*,"Input(): Width0 > MAXWidth (=",MAXWidth,")"
     IErr= 1
  ENDIF
  
  IF( Width1.GT.MAXWidth ) THEN
     PRINT*,"Input(): Width1 > MAXWidth (=",MAXWidth,")"
     IErr= 1
  ENDIF
  
  IF( (Width0.GT.Width1) .AND. (dWidth.GT.0) ) THEN
     PRINT*,"Input(): Width0 > Width1 and dWidth>0"
     IErr= 1
  ENDIF

  SELECT CASE(IDimenFlag)
  CASE(21,22,23,24,31,32,33,34)
     CONTINUE
  !CASE(32)
    ! PRINT*,"Input(): IDimenFlag=",IDimenFlag," outside allowed range!"
     !IErr= 1
  CASE DEFAULT
     PRINT*,"Input(): IDimenFlag=",IDimenFlag," outside allowed range!"
     IErr= 1
  END SELECT

  IF( IBCFlag.GT.MAXBCFlag ) THEN
     PRINT*,"Input(): IBCFlag > MAXBCFlag (=",&
          MAXBCFlag,")"
     IErr= 1
  ENDIF
  
  IF( IRNGFlag.GT.MAXRNGFlag ) THEN
     PRINT*,"Input(): IRNGFlag > MAXRNGFlag (=",&
          MAXRNGFlag,")"
     IErr= 1
  ENDIF
  
  IF( IKeepFlag.GT.MAXKeepFlag ) THEN
     PRINT*,"Input(): IKeepFlag > MAXKeepFlag (=",&
          MAXKeepFlag,")"
     IErr= 1
  ENDIF
  
  IF( IWriteFlag.GT.MAXWriteFlag ) THEN
     PRINT*,"Input(): IWriteFlag > MAXWriteFlag (=",&
          MAXWriteFlag,")"
     IErr= 1
  ENDIF
  
  IF( ISortFlag.GT.MAXSortFlag ) THEN
     PRINT*,"Input(): ISortFlag > MAXSortFlag (=",&
          MAXSortFlag,")"
     IErr= 1
  ENDIF
  
  IF( IFluxFlag.GT.MAXFluxFlag ) THEN
     PRINT*,"Input(): IFluxFlag > MAXFluxFlag (=",&
          MAXFluxFlag,")"
     IErr= 1
  ENDIF
  
  IF( Epsilon.LE.0.0D0) THEN
     PRINT*,"Input(): nonpositive Epsilon"
     IErr= 1
  ENDIF

  IF( SmallDisOnRim.LE.0.0D0) THEN
     PRINT*,"Input(): nonpositive SmallDisOnRim"
     IErr= 1
  ENDIF

  RETURN
  
  !	error in OPEN detected
120 PRINT*,"Input(): ERR in OPEN"
  GOTO 1000
  
  !	error in CLOSE detected
130 PRINT*,"Input(): ERR in CLOSE"
  GOTO 1000
  
  !	error in READ detected
20 PRINT*,"Input(): ERR in READ at line", ILine
  GOTO 1000
  
  !	EOF in READ occured prematurely
30 PRINT*,"Input(): EOF in READ at line", ILine
  
  ! dump the input help
  
1000 &
  PRINT*,"Input parameters:          ; explanation:"
  PRINT*,"--------------------------------------------------------------------"
  PRINT*,"ISeed         = 123456     ; (1) seed"
  
  PRINT*,"NOfIter       = 1          ; (2) # steps for SINGLE config."
  PRINT*,"NOfOrtho      = 10         ; (3) # steps until reorthonormalization"
  PRINT*,"NOfPrint      = 10         ; (4) # steps until printout"
  PRINT*,"NOfGamma      = 1          ; (5) # smallest Lyapunov exponents"

  PRINT*,"IDimenFlag    = 2          ; (8) 2/3 = 2/3 dimensions"
  PRINT*,"IBCFlag       = 0          ; (9) 0/1/2 = hard wall/periodic/antiperiodic BC"
  PRINT*,"IRNGFlag      = 0          ; (10) 0/1/2 = uniform/rescaled uniform/Gaussian disorder"
  PRINT*,"IKeepFlag     = 0          ; (11) 0/1 = yes/no data overwrite"
  PRINT*,"IWriteFlag    = 1          ; (12) 0/1/2/3/4 = no/log/category/wave fcn/Gamma/RHO output"
  PRINT*,"ISortFlag     = 0          ; (13) 0/1 = no/yes ReSort()"
  PRINT*,"IFluxFlag     = 1          ; (14) 0/1 = DiagDis/Energy loop"
  
  PRINT*,"Width0        = 0          ; (15) minimal width"
  PRINT*,"Width1        = 0          ; (16) maximal width"
  PRINT*,"dWidth        = 0          ; (17) width increment"
  
  PRINT*,"DiagDis0      = 0.         ; (18) minimal diagonal disorder"
  PRINT*,"DiagDis1      = 5.         ; (19) maximal  diagonal disorder"
  PRINT*,"dDiagDis      = 0.5        ; (20) increment of diagonal disorder"
  
  PRINT*,"Energy0       = 0.         ; (21) minimal energy"
  PRINT*,"Energy1       = 5.         ; (22) maximal energy"
  PRINT*,"dEnergy       = 0.5        ; (23) increment of energy"
  
  PRINT*,"Kappa         = 1.0        ; (24) inter-layer hopping"
  PRINT*,"Epsilon       = 5.0E-2     ; (26) accuracy goal of iteration"
  PRINT*,"SmallDisOnRim = 0.01       ; (27) small disorder on rim sites of Lieb models to avoid divergence"
  
  IErr= 1
  RETURN
  
END SUBROUTINE Input

!	--------------------------------------------------------------------
!	FileNameAvg
!
!	IErr	error code

CHARACTER*50 FUNCTION FileNameAvg(IWidth, fluxStr, fluxVal)
  USE MyNumbers
  
  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara

  INTEGER IWidth
  REAL(RKIND) fluxVal
  CHARACTER*1 fluxStr

  IF( fluxVal.GE.-1.0D-10 ) THEN

     WRITE(FileNameAvg, 100) &
          "L",IDimenFlag,"_",&
          IWidth,"_",&
          fluxStr,NINT(100.0D0*ABS(fluxVal)),&
          ".raw"
100  FORMAT(A1,I2.2,A1,I4.4,A1,A1,I4.4,A4)

  ELSE

     WRITE(FileNameAvg, 200)                         &
          "L",IDimenFlag,"_",&
          IWidth,"_",&
          fluxStr,"-",NINT(100.0D0*ABS(fluxVal)),&
          ".raw"
200  FORMAT(A1,I2.2,A1,I4.4,A1,A1,A1,I4.4,A4)

  ENDIF

END FUNCTION FileNameAvg

!	--------------------------------------------------------------------
!	CheckOutputAvg:
!
!	IErr	error code

SUBROUTINE CheckOutputAvg( filename, IErr )

  USE MyNumbers
  
  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara
  
  USE IChannels
  
  INTEGER IErr
  CHARACTER*50 filename
  
  !	PRINT*,"DBG: CheckOutputAvg()"
  
  IErr= 0
  
  !   WRITE out the input parameter

  OPEN(UNIT= IChOut, ERR= 10, STATUS= 'NEW',&
       FILE=TRIM(filename))
  
  IErr= 0
  
20 CLOSE(UNIT= IChOut, ERR= 100)
  
  RETURN
  
10 WRITE(*,15) filename
15 FORMAT(" CheckOutputAvg(): ", A10,&
        " exists -- skipped!")
  
  IErr= 2
  GOTO 20
  
  !  ERR in CLOSE detected
100 &
  PRINT*,"CheckOutputAvg(): ERR in CLOSE()"
  IErr= 1
  RETURN
  
END SUBROUTINE CheckOutputAvg


!	--------------------------------------------------------------------
!	OpenOutputAvg:
!
!	IErr	error code

SUBROUTINE OpenOutputAvg( filename, IWidth, IErr )
  
  USE MyNumbers
  
  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara
  
  USE IChannels
  
  INTEGER IWidth, IErr
  CHARACTER*50 filename
    
  INTEGER ICh, IChList(MAXFiles)

  PRINT*,"DBG: OpenOutputAvg()"
  
  IErr= 0
  
  OPEN(UNIT= IChOut, ERR= 10, STATUS= 'UNKNOWN',&
       FILE=TRIM(filename))

  IChList(1)= IChOut
  IChList(2)= IChOutPsi
  
  DO 1000 ICh= 1, 1
     
     WRITE(IChList(ICh),90,ERR=20) RStr,DStr,AStr
90   FORMAT("(* ",3A," *)")
     
     WRITE(IChList(ICh),100,ERR=20) ISeed
100  FORMAT("ISeed        = ", I15.1, "; (* + Config *)")
     
     WRITE(IChList(ICh),120,ERR=20) NOfIter
120  FORMAT("NOfIter      = ", I15.1, ";")
     
     WRITE(IChList(ICh),140,ERR=20) NOfOrtho
140  FORMAT("NOfOrtho     = ", I15.1, ";")
     
     WRITE(IChList(ICh),145,ERR=20) NOfPrint
145  FORMAT("NOfPrint     = ", I15.1, ";")
     
     WRITE(IChList(ICh),150,ERR=20) NOfGamma
150  FORMAT("NOfGamma     = ", I15.1, ";")

     WRITE(IChList(ICh),152,ERR=20) IDimenFlag
152  FORMAT("IDimenFlag   = ", I15.1, ";")
     
     WRITE(IChList(ICh),155,ERR=20) IBCFlag
155  FORMAT("IBCFlag      = ", I15.1, ";")
     
     WRITE(IChList(ICh),160,ERR=20) IRNGFlag
160  FORMAT("IRNGFlag     = ", I15.1, ";")
     
     WRITE(IChList(ICh),170,ERR=20) IKeepFlag
170  FORMAT("IKeepFlag    = ", I15.1, ";")
     
     WRITE(IChList(ICh),180,ERR=20) IWriteFlag
180  FORMAT("IWriteFlag   = ", I15.1, ";")
     
     WRITE(IChList(ICh),190,ERR=20) ISortFlag
190  FORMAT("ISortFlag    = ", I15.1, ";")
     
     WRITE(IChList(ICh),195,ERR=20) IFluxFlag
195  FORMAT("IFluxFlag    = ", I15.1, ";")
     
     WRITE(IChList(ICh),200,ERR=20) Width0
200  FORMAT("Width0       = ", I15.1, ";")
     
     WRITE(IChList(ICh),210,ERR=20) Width1
210  FORMAT("Width1       = ", I15.1, ";")
     
     WRITE(IChList(ICh),220,ERR=20) dWidth
220  FORMAT("dWidth       = ", I15.1, ";")
     
     WRITE(IChList(ICh),290,ERR=20) DiagDis0
290  FORMAT("DiagDis0     = ", G18.9, ";")
     
     WRITE(IChList(ICh),293,ERR=20) DiagDis1
293  FORMAT("DiagDis1     = ", G18.9, ";")
     
     WRITE(IChList(ICh),296,ERR=20) dDiagDis
296  FORMAT("dDiagDis     = ", G18.9, ";")
     
     WRITE(IChList(ICh),300,ERR=20) Energy0
300  FORMAT("energy0      = ", G18.9, ";")
     
     WRITE(IChList(ICh),302,ERR=20) Energy1
302  FORMAT("energy1      = ", G18.9, ";")
     
     WRITE(IChList(ICh),304,ERR=20) dEnergy
304  FORMAT("denergy      = ", G18.9, ";")
     
     WRITE(IChList(ICh),310,ERR=20) Epsilon
310  FORMAT("epsilon      = ", G18.9, ";")

     WRITE(IChList(ICh),315,ERR=20) SmallDisOnRim
315  FORMAT("SmallDisOnRim= ", G18.9, ";")
     
     WRITE(IChList(ICh),400,ERR=20) IWidth
400  FORMAT("Width        = ", I15.1, ";")
     
     WRITE(IChList(ICh),500,ERR=20)
500  FORMAT("data= {")
     
1000 CONTINUE

  RETURN
     
  !	error in OPEN detected
10 PRINT*,"OpenOutputAvg(): ERR in OPEN()"
  IErr= 1
  RETURN
  
  !	error in WRITE detected
20 PRINT*,"OpenOutputAvg(): ERR in WRITE()"
  IErr= 1
  RETURN
  
END SUBROUTINE OpenOutputAvg


! --------------------------------------------------------------------
! ReOpenOutputAvg:
!
! IErr	error code

SUBROUTINE ReOpenOutputAvg( filename, IErr )
  
  USE MyNumbers
  
  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara
  
  USE IChannels
  
  INTEGER IErr
  CHARACTER*50 filename
  
  INTEGER ICh, IChList(MAXFiles)

  ! PRINT*,"DBG: ReOpenOutputAvg()"

  IErr= 0
  
  OPEN(UNIT= IChOut, ERR= 10, STATUS= 'OLD',&
       FILE=filename, POSITION='APPEND')
  
  RETURN
  
  ! error in OPEN detected
10 PRINT*,"ReOpenOutputAvg(): ERR in OPEN()"
  PRINT*, "file ", filename, " does not exist --- REOPEN not possible!"
  IErr= 1
  RETURN
  
  ! error in WRITE detected
20 PRINT*,"ReOpenOutputAvg(): ERR in WRITE()"
  IErr= 1
  RETURN
  
END SUBROUTINE ReOpenOutputAvg


! --------------------------------------------------------------------
! WriteOutputAvg:
!
! IErr	error code

SUBROUTINE WriteOutputAvg(&
     IWidth, 	&
     IConvergence, &
     DiagDis, 	&
     Energy,	&
     gam, var, NOfL,  &
     psi,	&
     IErr )
        
  USE MyNumbers
  
  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara
  
  USE IChannels
  
  INTEGER IWidth, IConvergence, NOfL, IErr
  REAL(KIND=RKIND) DiagDis, Energy,&
       psi(IWidth,IWidth)
  
  REAL(KIND=RKIND) gam(NOfL), var(NOfL)
  
  INTEGER iState, jSite, iL
  
  !PRINT*,"DBG: WriteOutputAvg()"
  
  IErr= 0
  
  ! average Lyapunov exponent Gamma
  
  DO 100 iL= 1,NOfL
     
     WRITE(IChOut,410,ERR=10) &
          iL, &
          DiagDis, &
          Energy, &
          gam(iL), var(iL), IConvergence
     
410  FORMAT("{ ", I7.1, &
          ", ", F15.6, &
          ", ", F15.6, &
          ", ", F25.16, ", ", F25.16, &
          ", ", I2.1, &
          " },")

!410  FORMAT(I7.1,F15.6,F15.6,F25.16,F25.16)

100 &
  CONTINUE
     
  !	wave functions
     
     IF( IWriteFlag.GE.MAXWriteFlag+1 ) THEN
        
        DO 1000 iState= 1, IWidth
           
           WRITE(IChOutPsi,510,ERR=10) 
510        FORMAT("{ ")
           
           DO 2000 jSite= 1, IWidth
              WRITE(IChOutPsi,550,ERR=10) iState,jSite,&
                   psi(jSite,IWidth+1-iState)**2
550           FORMAT( I4.1, " ", I4.1, " ", 1(G25.15) )
2000       CONTINUE
              
           WRITE(IChOutPsi,570,ERR=10) 
570        FORMAT("}, ")
              
1000    CONTINUE
              
      ENDIF
           
      RETURN
           
      !	ERR in Write detected
10    PRINT*,"WriteOutputAvg(): ERR in WRITE()"
      IErr= 1
      RETURN
           
END SUBROUTINE WriteOutputAvg

!--------------------------------------------------------------------
! WriteOutputPsi: Write Lyapunov exponent data if converged!
!
! IErr	error code
!---------------------------------------------------------------------

SUBROUTINE WriteOutputPsi(&
     Ilayer, IWidth, 	&
     psi,	&
     IErr)
        
  USE MyNumbers
  
  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara
  
  USE IChannels
  
  INTEGER Ilayer, IErr, IWidth
  COMPLEX(KIND=CKIND) psi(IWidth,IWidth)
  
  INTEGER jState, iSite

  ! wavefunctions

    IF( IWriteFlag.GE.MAXWriteFlag+1 ) THEN
        
        DO jState= 1, 1
           
           DO iSite= 1, IWidth
              WRITE(IChOutPsi,550,ERR=10) Ilayer, iSite,&
                   psi(IWidth+1-jState,iSite)
550           FORMAT( I6.1, " ", I4.1, " ", 1(F25.15), " ", 1(F25.15))
           ENDDO
              
        ENDDO
              
      ENDIF
           

      RETURN
           
      !	ERR in Write detected
10    PRINT*,"WriteOutputAvg(): ERR in WRITE()"
      IErr= 1
      RETURN
           
END SUBROUTINE WriteOutputPsi

! --------------------------------------------------------------------
! CloseOutputAvg:
!
! IErr	error code

SUBROUTINE CloseOutputAvg( IErr )

  USE MyNumbers
  
  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara
  
  USE IChannels
  
  INTEGER IErr
  
  INTEGER ICh, IChList(MAXFiles)
  
  !PRINT*,"DBG: CloseOutputAvg()"
  
  IErr= 0
  
  IChList(1)= IChOut
  IChList(2)= IChOutPsi
  ICHList(3)= IChOutAvgRHO1
  ICHList(4)= IChOutAvgRHOL
  ICHList(5)= IChOutHAV
  
  DO ICh= 1, MAXFiles
     
     SELECT CASE(IChList(ICh))
     CASE(IChOut)
        WRITE(IChList(ICh),100,ERR=20)
100     FORMAT('};')
     CASE DEFAULT
        CONTINUE
     END SELECT
     
     CLOSE(UNIT= IChList(ICh), ERR= 10)
     
  ENDDO
  RETURN
     
  ! ERR in CLOSE detected
10 PRINT*,"CloseOutputAvg(): ERR in CLOSE()"
  IErr= 1
  RETURN
  
  ! error in WRITE detected
20 PRINT*,"CloseOutputAvg(): ERR in WRITE()"
  IErr= 1
  RETURN
  
END SUBROUTINE CloseOutputAvg


!	--------------------------------------------------------------------
!	OpenOutputGamma:
!
!	IErr	error code

SUBROUTINE OpenOutputGamma( IWidth, 				  &
     DiagDis, Energy,  &
     IErr )

  USE MyNumbers

  USE CConstants
  USE IConstants

  USE IPara
  USE DPara

  USE IChannels

  INTEGER IWidth, IErr
  REAL(KIND=RKIND) DiagDis, Energy

  INTEGER ICh

  CHARACTER*28 CNameP
  CHARACTER*29 CNameM

  ! PRINT*,"DBG: OpenOutputGamma()"

  IErr= 0
  ICh = IChOutGam

  !	the filename is different for this gamma logging

  IF( Energy.GE.-1.0D-10 ) THEN

     WRITE(CNameP, 100) &
          "L",IDimenFlag,"_", &
          IWidth,".", &
          NINT(100.0D0*ABS(DiagDis)),".", &
          NINT(100.0D0*ABS(Energy))
100  FORMAT(A1,I2.2,A1,I4.4,A1,I4.4,A1,I4.4)

     OPEN(UNIT= ICh, ERR= 10, STATUS= 'UNKNOWN', &
          FILE=CNameP)
  ELSE

     WRITE(CNameM, 200) &
          "L",IDimenFlag,"_", &
          IWidth,".", &
          NINT(100.0D0*ABS(DiagDis)),".-", &
          NINT(100.0D0*ABS(Energy))
200  FORMAT(A1,I2.2,A1,I4.4,A1,I4.4,A2,I4.4) 

     OPEN(UNIT= ICh, ERR= 10, STATUS= 'UNKNOWN', &
          FILE=CNameM)
  ENDIF

  RETURN

  !	error in OPEN detected
10 PRINT*,"OpenOutputGamma(): ERR in OPEN()"
  IErr= 1
  RETURN

END SUBROUTINE OpenOutputGamma


!	--------------------------------------------------------------------
!	WriteOutputGamma:
!
!	IErr	error code

SUBROUTINE WriteOutputGamma( index, gam, var, NOfL, IErr )
        
   USE MyNumbers

	USE CConstants
	USE IConstants

	USE IPara
	USE DPara

	USE IChannels

	INTEGER index, NOfL, IErr
	REAL(KIND=RKIND) gam(NOfL), var(NOfL)

	INTEGER ICh, iL

! 	PRINT*,"DBG: WriteOutputGamma()"

	IErr= 0
	ICh = IChOutGam

!	Lyapunov exponent Gamma

	DO 100 iL=1,NOfL
		WRITE(ICh,410,ERR=10) index, iL, gam(iL), var(iL)
410	FORMAT(" ",I7.1, " ", I4.1, " ", G25.16, " ", G25.16)
100 &
	CONTINUE

	RETURN

!	ERR in Write detected
10	PRINT*,"WriteOutputGamma(): ERR in WRITE()"
	IErr= 1
	RETURN

END SUBROUTINE WriteOutputGamma


!	--------------------------------------------------------------------
!	CloseOutputGamma:
!
!	IErr	error code

SUBROUTINE CloseOutputGamma( IErr )

	USE CConstants
	USE IConstants

	USE IPara
	USE DPara

	USE IChannels

	INTEGER IErr

	INTEGER ICh

!	PRINT*,"DBG: CloseOutputGamma()"

	IErr= 0
	ICh = IChOutGam

	CLOSE(UNIT= ICh, ERR= 10)

	RETURN

!	ERR in CLOSE detected
10	PRINT*,"CloseOutputGamma(): ERR in CLOSE()"
	IErr= 1
	RETURN

END SUBROUTINE CloseOutputGamma


!	--------------------------------------------------------------------
!	WriteOutputAvgRHO:
!
!	IErr	error code

SUBROUTINE WriteOutputAvgRHO(  		&
     IChannelOut, IWidth, IChannelMax,	&
     DiagDis, Energy,			&
     rhoMat, NOfL, NSamples,		&
     IErr )
  
   USE MyNumbers
   
   USE CConstants
   USE IConstants

   USE IPara
   USE DPara
	
   USE IChannels

   INTEGER IWidth, IChannelOut, IChannelMax, NOfL, NSamples, IErr
   REAL(KIND=RKIND) DiagDis, Energy

   REAL(KIND=RKIND) rhoMat(IWidth,0:IChannelMax) 

   INTEGER iState, jChannel, iL

! PRINT*,"DBG: WriteOutputAvgRHO()"

	IErr= 0

!	average Lyapunov exponent Gamma

   DO jChannel=0,IChannelMax
       DO iState= 1,IWidth

		   WRITE(IChannelOut,410,ERR=10) &
                jChannel, iState, &
                DiagDis, Energy,	&
                rhoMat(iState,jChannel)/REAL(NSamples)

410	   FORMAT(" ", I7.1, " ", I7.1, &
                " ", F15.6, " ", F15.6, &
                " ", G25.16)
      ENDDO
   ENDDO

	RETURN

!	ERR in Write detected
10	PRINT*,"WriteOutputAvgRHO(): ERR in WRITE()"
	IErr= 1
	RETURN

END SUBROUTINE WriteOutputAvgRHO

!	--------------------------------------------------------------------
!	WriteOutputAvgRHOX:
!
!	IErr	error code

SUBROUTINE WriteOutputAvgRHOX(  	&
     IChannelOut, IWidth, IChannelMax,	&
     DiagDis, Energy,			&
     rhoMat, NOfL, NSamples,		&
     IErr )
  
   USE MyNumbers
   
   USE CConstants
   USE IConstants

   USE IPara
   USE DPara
	
	USE IChannels

	INTEGER IWidth, IChannelOut, IChannelMax, NOfL, NSamples, IErr
   REAL(KIND=RKIND) DiagDis, Energy

	REAL(KIND=RKIND) rhoMat(0:IChannelMax,0:IChannelMax) 

   INTEGER iChannel, jChannel, iL

! PRINT*,"DBG: WriteOutputAvgRHO()"

	IErr= 0

!	average Lyapunov exponent Gamma

   DO iChannel=0,IChannelMax
      DO jChannel= 0,IChannelMax
         
         WRITE(IChannelOut,410,ERR=10) &
              iChannel, jChannel, &
              DiagDis, Energy,	&
              rhoMat(iChannel,jChannel)/REAL(NSamples)
         
410      FORMAT(" ", I7.1, " ", I7.1, &
              " ", F15.6, " ", F15.6, &
              " ", G25.16)
      ENDDO
   ENDDO
   
   RETURN
   
   !	ERR in Write detected
10 PRINT*,"WriteOutputAvgRHOX(): ERR in WRITE()"
   IErr= 1
   RETURN

END SUBROUTINE WriteOutputAvgRHOX

! --------------------------------------------------------------------
! SaveCurrentParameters:
!

SUBROUTINE SaveCurrentParameters( &
  IWidth, flux, IErr)
  
  USE MyNumbers
  
  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara
  
  USE IChannels
  
  INTEGER IWidth, IErr
  REAL(KIND=RKIND) flux

  CHARACTER*8 CName
  
  ! PRINT*,"DBG: ReOpenOutputAvg()"

  IErr= 0
  
  ! WRITE out the input parameter
  
  OPEN(UNIT= IChtmp, ERR= 10, STATUS= 'UNKNOWN',&
       FILE="tmseLMxD.tmp")

  WRITE(IChtmp,100,ERR=20) IWidth, flux
100 FORMAT(I15.1,F18.9)

  CLOSE(UNIT= IChtmp, ERR=40)

  RETURN

  ! ERR in Open detected
10 PRINT*,"SaveCurrentParameters(): ERR in OPEN()"
  IErr= 1
  RETURN
  
  ! error in WRITE detected
20 PRINT*,"SaveCurrentParameters(): ERR in WRITE()"
  IErr= 1
  RETURN
  
  ! ERR in CLOSE
40 PRINT*,"SaveCurrentParameters(): ERR in CLOSE()"
  IErr= 1
  RETURN

END SUBROUTINE SaveCurrentParameters

! --------------------------------------------------------------------
! ReloadCurrentParameters:
!

SUBROUTINE ReloadCurrentParameters( &
  IWidth, flux, IErr)
  
  USE MyNumbers
  
  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara
  
  USE IChannels
  
  INTEGER IWidth, IErr
  REAL(KIND=RKIND) flux

  CHARACTER*8 CName
  
  ! PRINT*,"DBG: ReOpenOutputAvg()"

  IErr= 0
  
  ! WRITE out the input parameter
        
  OPEN(UNIT= IChtmp, ERR= 10, STATUS= 'OLD',&
       FILE="tmseLMxD.tmp")

  READ(IChtmp,100,ERR=20,END=30) IWidth, flux
100 FORMAT(I15.1,F18.9)

  CLOSE(UNIT= IChtmp, ERR=40)

  RETURN

  ! ERR in Open detected
10 PRINT*,"ReloadCurrentParameters(): ERR in OPEN()"
  IErr= 1
  RETURN
  
  ! error in READ detected
20 PRINT*,"ReloadCurrentParameters(): ERR in READ()"
  IErr= 1
  RETURN
  
  ! EOF in READ occured prematurely
30 PRINT*,"ReloadCurrentParameters(): EOF in READ()"
  IErr= 1
  RETURN
  
  ! ERR in CLOSE
40 PRINT*,"ReloadCurrentParameters(): ERR in CLOSE()"
  IErr= 1
  RETURN

END SUBROUTINE ReloadCurrentParameters

! --------------------------------------------------------------------
! DeleteCurrentParameters:
!

SUBROUTINE DeleteCurrentParameters( &
  IWidth, IErr)
  
  USE MyNumbers
  
  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara
  
  USE IChannels
  
  INTEGER IWidth, IErr

  ! PRINT*,"DBG: ReOpenOutputAvg()"

  IErr= 0
  
  ! DELETE the parameter file
   
  OPEN(UNIT=IChtmp)
  CLOSE(UNIT=IChtmp,STATUS='DELETE',ERR=10)

  RETURN

  ! ERR in Delete detected
10 PRINT*,"DeleteCurrentParameters(): ERR in DELETE()"
  RETURN
  
END SUBROUTINE DeleteCurrentParameters

!	--------------------------------------------------------------------
!	OpenOutputRHO:
!
!	IErr	error code

SUBROUTINE OpenOutputRHO( IWidth, &
           DiagDis, Energy,  &
           IErr )

  USE MyNumbers
  
  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara
  
  USE IChannels
  
  INTEGER IWidth, IErr
  REAL(KIND=RKIND) DiagDis, Energy
  
  INTEGER ICh
  
  CHARACTER*28 CNameP
  CHARACTER*29 CNameM
  
  ! PRINT*,"DBG: OpenOutputRHO()"
  
  IErr= 0
  ICh = IChOutRHO
  
  !	the filename is different for this RHO logging
  
  IF( Energy.GE.-1.0D-10 ) THEN
     
     WRITE(CNameP, 100) &
          IWidth,".",	&
          NINT(100.0D0*ABS(DiagDis)),".", &
          NINT(100.0D0*ABS(Energy)),".rho"
100  FORMAT(I4.4,A1,I4.4,A1,I4.4,A4)
     
     OPEN(UNIT= ICh, ERR= 10, STATUS= 'UNKNOWN', &
          FILE=CNameP)
  ELSE
     
     WRITE(CNameM, 200) &
          IWidth,".",	&
          NINT(100.0D0*ABS(DiagDis)),".-", &
          NINT(100.0D0*ABS(Energy)),".rho"
200  FORMAT(I4.4,A1,I4.4,A2,I4.4,A4) 
     
     OPEN(UNIT= ICh, ERR= 10, STATUS= 'UNKNOWN', &
          FILE=CNameM)
  ENDIF
  
  RETURN
  
  !	error in OPEN detected
10 PRINT*,"OpenOutputRHO(): ERR in OPEN()"
  IErr= 1
  RETURN
  
END SUBROUTINE OpenOutputRHO


!	--------------------------------------------------------------------
!	WriteOutputRHO:
!
!	IErr	error code

SUBROUTINE WriteOutputRHO(  &
     IChannelOut, IWidth, IChannelMax,	&
     DiagDis, Energy,	&
     rhoMat, NOfL, NSamples,	&
     IErr )
  
  USE MyNumbers
  
  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara
  
  USE IChannels
  
  INTEGER IWidth, IChannelOut, IChannelMax, NOfL, NSamples, IErr
  REAL(KIND=RKIND) DiagDis, Energy
  
  REAL(KIND=RKIND) rhoMat(IWidth,0:IChannelMax) 
  
  INTEGER iState, jChannel, iL
  
  ! PRINT*,"DBG: WriteOutputRHO()"
  
  IErr= 0
  
  !	average Lyapunov exponent Gamma
  
  DO jChannel=0,IChannelMax
     DO iState= 1,IWidth
        
        WRITE(IChannelOut,410,ERR=10) &
             NSamples,jChannel, iState, &
             DiagDis, Energy,	&
             rhoMat(iState,jChannel)/REAL(NSamples)
        
410     FORMAT(" ", I7.1, " ", I7.1, " ", I7.1, &
             " ", F15.6, " ", F15.6, &
             " ", G25.16)
     ENDDO
  ENDDO
  
  RETURN
  
  !	ERR in Write detected
10 PRINT*,"WriteOutputAvgRHO(): ERR in WRITE()"
  IErr= 1
  RETURN
  
END SUBROUTINE WriteOutputRHO


!	--------------------------------------------------------------------
!	CloseOutputRHO:
!
!	IErr	error code

SUBROUTINE CloseOutputRHO( IErr )
  
  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara
  
  USE IChannels
  
  INTEGER IErr
  
  INTEGER ICh
  
  !	PRINT*,"DBG: CloseOutputRHO()"
  
  IErr= 0
  ICh = IChOutRHO
  
  CLOSE(UNIT= ICh, ERR= 10)
  
  RETURN
  
  !	ERR in CLOSE detected
10 PRINT*,"CloseOutputRHO(): ERR in CLOSE()"
  IErr= 1
  RETURN
  
END SUBROUTINE CloseOutputRHO

! --------------------------------------------------------------------
! WriteOutputHAV:
!
! IErr	error code

SUBROUTINE WriteOutputHAV(&
     IWidth, 	&
     DiagDis, 	&
     Energy,	&
     gam, NOfL,  &
     sumhavg,	&
     IErr )
  
  USE MyNumbers
  
  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara
  
  USE IChannels
  
  INTEGER IWidth, NOfL, IErr
  REAL(KIND=RKIND) DiagDis, Energy
  REAL(KIND=RKIND) gam(NOfL), sumhavg

  REAL(KIND=RKIND) sumgam
  
  INTEGER iState, jSite, iL
  
  ! PRINT*,"DBG: WriteOutputHAV()"
  
  IErr= 0
  
  ! sum of Lyapunov exponents

  sumgam  = 0.0D0
  DO iL= 1,NOfL
     sumgam   = sumgam + gam(iL)
  ENDDO

  WRITE(IChOutHAV,410,ERR=10) &
     NOfL, &
     DiagDis, &
     Energy, &
     sumgam, sumhavg
     
410 FORMAT(I7.1, &
         ", ", F15.6, &
         ", ", F15.6, &
         ", ", F25.16, ", ", F25.16, ", ", F25.16)
  
  RETURN
           
  !	ERR in Write detected
10 PRINT*,"WriteOutputHAV(): ERR in WRITE()"
  IErr= 1
  RETURN
         
END SUBROUTINE WriteOutputHAV

!	--------------------------------------------------------------------
!	OutputEVals:
!
!	IErr	error code

SUBROUTINE OutputEVals( ILayer, IWidth, IElements, Disorder,Energy, EVals, Csur, IErr )
  
  USE MyNumbers
  
  !USE CConstants
  !USE IConstants
  
  USE IPara
  !USE DPara
  
  USE IChannels
  
  INTEGER ILayer, IWidth, IElements, iState, IErr
  REAL(KIND=RKIND) Disorder, Energy
  CHARACTER*4 Csur
  COMPLEX(KIND=CKIND) EVals(IWidth)
  !EXTERNAL ARG

  CHARACTER*38 CNameP
  CHARACTER*39 CNameM
  CHARACTER*1 Cpre
  
  !PRINT*,"DBG: OutputEVals()"
  
  IErr= 0
  
  ! WRITE out the input parameter

  !IF(IStripeFlag.LT.0) THEN
  !   WRITE(Cpre, '(A1)') "m"
  !ELSE
     WRITE(Cpre, '(A1)') "-"
  !ENDIF
   
  IF( Energy.GE.-1.0D-10 ) THEN
     
     WRITE(CNameP, '(I4.4,A1,I6.6,A1,I6.6,A1,I6.6,A1,I5.5,A4)') &
          IWidth,".", &
          ISeed, ".", &
          NINT(10000.0D0*ABS(Disorder)),".",		&
          NINT(10000.0D0*ABS(Energy)),".",ILayer/(NOfPrint*NOfOrtho),Csur
     
     !PRINT*,CNameP

     OPEN(UNIT= IChOutPsi, ERR= 10, STATUS= 'UNKNOWN', &
          FILE=CNameP)
  ELSE
     
     WRITE(CNameM, '(I4.4,A1,I6.6,A1,I6.6,A2,I6.6,A1,I5.5,A4) ') &
          IWidth,".", &
          ISeed, ".", &
          NINT(10000.0D0*ABS(Disorder)),".-",		 &
          NINT(10000.0D0*ABS(Energy)),".",ILayer/(NOfPrint*NOfOrtho),Csur
     
     !PRINT*,CNameM

     OPEN(UNIT= IChOutPsi, ERR= 10, STATUS= 'UNKNOWN', &
          FILE=CNameM)
  ENDIF

  DO iState= 1, IElements
     
     WRITE(IChOutPsi,550,ERR=30) iState,&
          DBLE(EVals(iState)), AIMAG(EVals(iState)), &
          !-ATAN(AIMAG(EVals(iState))/DBLE(EVals(iState)))
          ARG(AIMAG(EVals(iState)),DBLE(EVals(iState)))
550  FORMAT( I4.1, " ", 3(G25.15) )
  ENDDO
     
  CLOSE(UNIT=IChOutPsi,ERR=40)
     
  RETURN
     
  !	error in OPEN detected
10 PRINT*,"OutputEVals(): ERR in OPEN()"
  IErr= 1
  RETURN
  
  !	error in WRITE detected
20 PRINT*,"OutputEVals(): ERR in WRITE()"
  IErr= 1
  RETURN

  !	ERR in Write detected
30 PRINT*,"OutputEVals(): ERR in WRITE()"
  IErr= 1
  RETURN

  ! ERR in CLOSE detected
40 PRINT*,"OutputEVals(): ERR in CLOSE()"
  IErr= 1
  RETURN

END SUBROUTINE OutputEVals

! --------------------------------------------------------------------
! WriteDataC

SUBROUTINE WriteDataC( &
     prefix, surname, data, size, step, IErr)

  USE MyNumbers

  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara

  USE IChannels

  CHARACTER*40 surname
  CHARACTER*6 prefix,postfix
  INTEGER(KIND=IKIND) size, step, IErr
  COMPLEX(KIND=CKIND) data(size)

  CHARACTER*50 filename
  INTEGER index

  !PRINT*,"WriteDataC()"

  WRITE(filename,'(A6,A40,A4)') prefix,surname,".raw"

  IF(IWriteFlag.EQ.MAXWriteFlag) PRINT*,"WriteDataC() for ", filename

  IErr=0

  OPEN(UNIT= IChOutWrite, ERR= 10, STATUS= 'UNKNOWN',&
          FILE=filename)

  DO index=1,size,step
     !PRINT*,"DBG: index=", index
	 IF (ABS(data(index)).GE. TINY) THEN
     WRITE(IChOutWrite,100) DBLE(data(index)),AIMAG(data(index)), &
          !-ATAN(AIMAG(EVals(iState))/DBLE(EVals(iState)))
          !ARG(AIMAG(data(index)),DBLE(data(index)))
          ARG(DBLE(data(index)),AIMAG(data(index)))
	 ELSE
     WRITE(IChOutWrite,100) DBLE(data(index)),AIMAG(data(index)), &
          !-ATAN(AIMAG(EVals(iState))/DBLE(EVals(iState)))
          0.0D0
	 ENDIF
	 
100  FORMAT(3(G25.15))
  ENDDO
     
  CLOSE(IChoutWrite,ERR=30)
  
  RETURN

! error in OPEN detected
10 PRINT*,"WriteDataC(): ERR in OPEN()"
  PRINT*, "file ", filename, " does not exist --- REOPEN not possible!"
  IErr= 1
  RETURN
  
  ! error in WRITE detected
20 PRINT*,"WriteDataC(): ERR in WRITE() for file=", filename
  IErr= 1
  RETURN

  ! ERR in CLOSE detected
30 PRINT*,"WriteDataC(): ERR in CLOSE() for file=", filename
  IErr= 1
  RETURN

END SUBROUTINE WriteDataC

! --------------------------------------------------------------------
! WriteDataR

SUBROUTINE WriteDataR( &
     prefix, surname, data, size, step, IErr)

  USE MyNumbers

  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara

  USE IChannels

  CHARACTER*40 surname
  CHARACTER*6 prefix,postfix
  INTEGER(KIND=IKIND) size, step, IErr
  REAL(KIND=RKIND) data(size)

  CHARACTER*50 filename
  INTEGER index

  !PRINT*,"WriteDataR()"

  WRITE(filename,'(A6,A40,A4)') prefix,surname,".raw"

  IF(IWriteFlag.EQ.MAXWriteFlag) PRINT*,"WriteDataR() for ", filename

  IErr=0

  OPEN(UNIT= IChOutWrite, ERR= 10, STATUS= 'UNKNOWN',&
          FILE=filename)!, POSITION='APPEND')

  DO index=1,size,step
     !PRINT*,"DBG: index=", index
     WRITE(IChOutWrite,100) data(index)
100  FORMAT(G25.15)
  ENDDO
     
  CLOSE(IChoutWrite,ERR=30)
  
  RETURN

! error in OPEN detected
10 PRINT*,"WriteDataR(): ERR in OPEN()"
  PRINT*, "file ", filename, " does not exist --- REOPEN not possible!"
  IErr= 1
  RETURN
  
  ! error in WRITE detected
20 PRINT*,"WriteDataR(): ERR in WRITE() for file=", filename
  IErr= 1
  RETURN

  ! ERR in CLOSE detected
30 PRINT*,"WriteDataR(): ERR in CLOSE() for file=", filename
  IErr= 1
  RETURN

END SUBROUTINE WriteDataR

! --------------------------------------------------------------------
! DBGWritePsi

SUBROUTINE DBGWritePsi( &
     PsiA,PsiB,M, iter1,iter2, ilayer, IErr)

  USE MyNumbers

  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara

  USE IChannels

  CHARACTER*20 filename
  INTEGER(KIND=IKIND) M, iter1,iter2, ilayer, IErr
  REAL(KIND=RKIND) PsiA(M,M), PsiB(M,M)

  INTEGER index, jndex

  PRINT*,"DBGWritePsiAB()"

  ! PsiA
  
  WRITE(filename,'(A5,I4.4,A4)') "PsiA-",ilayer,".raw"

  OPEN(UNIT= IChOutPsi, STATUS= 'UNKNOWN',&
          FILE=filename)!, POSITION='APPEND')

  DO jndex=1,M
     DO index=1,M
        WRITE(IChOutPsi,100) ilayer, jndex, index, PsiA(index,jndex)
100  FORMAT("A",I4.1,I4.1,I4.1,G25.15)
     END DO
  END DO
     
  CLOSE(IChoutPsi)
  
  ! PsiB
  
  WRITE(filename,'(A5,I4.4,A4)') "PsiB-",ilayer,".raw"

  OPEN(UNIT= IChOutPsi, STATUS= 'UNKNOWN',&
          FILE=filename)!, POSITION='APPEND')

  DO jndex=1,M
     DO index=1,M
        WRITE(IChOutPsi,200) iter1,iter2,ilayer, jndex, index, PsiB(index,jndex)
200  FORMAT("B",I4.1,I4.1,I4.1,I4.1,I4.1,G25.15)
     END DO
  END DO
     
  CLOSE(IChoutPsi)
  
  RETURN
END SUBROUTINE DBGWritePsi
