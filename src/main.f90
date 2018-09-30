!********************************************************************
!
! TMSEXD - Transfer matrix method for the Anderson
! model with diagonal disorder in X dimensions
!
!********************************************************************

!********************************************************************
!
! $Header: /home/cvs/phsht/tmseXd/src/main.f90,v 1.4 2016/10/03 19:32:32 phsht Exp $
!
!********************************************************************

!**************************************************************************
!
! $Log: main.f90,v $
! Revision 1.4  2016/10/03 19:32:32  phsht
! new filenames for .raw files
!
! Revision 1.3  2016/10/03 17:33:32  phsht
! 1st attempt at correct 1D version
!
! Revision 1.2  2012/09/07 10:44:10  phsht
! removed MagFlux and hence all CMPLX,COMPLEX definitions;
! included a final column in AVG output to indicate convergence (or not)
!
! Revision 1.1.1.1  2012/09/07 10:00:09  phsht
! tmseXd
!
! Revision 1.1  2011/05/06 08:13:09  phsht
! 1st installement
!
! Revision 1.14  2010/11/17 15:56:11  phrkaj
! Neatened the code a bit
!
! Revision 1.13  2010/11/16 12:07:28  phrkaj
! Displays time taken whatever the writeflag is
!
! Revision 1.12  2010/11/11 11:16:25  phrkaj
! Renamed files to get rid of any mention of SB.
!
! Revision 1.11  2010/11/10 10:26:38  phrkaj
! So far working in 1D
!
! Revision 1.10  2010/11/04 15:14:02  phrkaj
! Added 3D stuff. Pre-debug
!
! Revision 1.9  2010/10/27 16:43:51  phrkaj
! Modified the output of the .raw file so that only the data is printed
!
! Revision 1.8  2010/10/26 09:43:39  phrkaj
! Deleted naive debugging statements, got rid of ILevelflag and IConvflag, deleted old logs
!
! Revision 1.7  2010/10/25 15:41:33  phsht
! small changes to remove a "malloc/glibc" error
!
! Revision 1.6  2010/10/25 14:59:04  phrkaj
! *** empty log message ***
!
! Revision 1.5  2010/10/25 14:31:07  phrkaj
! Debugging
!
! Revision 1.4  2010/10/24 14:17:29  phrkaj
! Got rid of some more SB stuff, changed a couple of things which were coming up with errors (probably because I'm using gfortran at home?)
!
! Revision 1.3  2010/10/22 13:10:42  phrkaj
! deleted some more old SB lines
!
! Revision 1.2  2010/10/22 12:44:02  phrkaj
! deleted old SB parts, not fully done yet
!
! Revision 1.1.1.1  2010/10/22 12:23:38  phsht
! ParaTMM
!**************************************************************************

!********************************************************************
!
!      Comments:
!
!      IBCFlag         0        hardwall or
!                      1        periodic or
!                      2        antiperiodic boundary conditions
!
!      IKeepFlag       0        overwrite old .raw files
!                      1        keep old .raw files and skip to next
!                               configuration
!
!      IWriteFlag      0        no protocolling
!                      1        write out final wave function for each
!                               strip width and flux (disorder/U)
!                      2        additionally, write out each nGamma
!                               value after reorthonormalisation
!
!      ISortFlag       0        NO sorting of eigenvalues/vectors after
!                               each reorthonormalization
!                      1        YES, sorting is done.
!
!      IFluxFlag       0        DiagDis loop
!                      1        Energy loop
!
!      Notation:
!
!      see B. Kramer and M. Schreiber, "Transfer-Matrix Methods and Finite-
!	Size Scaling for Disordered Systems" and K. Frahm et al., EPL 31,
!      169 (1995).
!
!      Also, Martin Hennecke, "Anderson Uebergang im Magnetfeld", PTB
!      Bericht, PTB-PG-6, Braunschweig, Dec 94.
!
!********************************************************************

PROGRAM TMSEXd

  !--------------------------------------------------------------------
  ! parameter and global variable definitions
  !--------------------------------------------------------------------
  
  USE MyNumbers
  
  USE CConstants
  USE IConstants
  USE IChannels
  
  USE IPara
  USE DPara
  
  USE RNG
  
  !--------------------------------------------------------------------
  ! local variable definitions
  !--------------------------------------------------------------------
  
  IMPLICIT NONE
  
  INTEGER IWidth, IWidthSquared, IWidthRL, Width0_X, &
       IChannelMax, index,jndex,kndex,lndex, jjndex,kkndex, andex,bndex, &
       Iter1,Iter2, iStep, NOfG,iG, Ilayer, TMM_CONVERGED
  
  REAL(KIND=RKIND) flux, flux0,flux1,dflux, &
       fluxRL,flux0_X, fluxVal
  REAL(KIND=RKIND) DiagDis, Energy, sum
  
  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE ::             &
       nGamma, gamma, gamma2, acc_variance
  
  REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE :: PsiA, PsiB

  CHARACTER*1 fluxStr
  CHARACTER*50 filenameAVG, filename
  CHARACTER*3 EnStr
  EXTERNAL FileNameAvg
  
  ! timing variables
  REAL(4), DIMENSION(1:3,1:1) :: TIME, TIMESUM
  REAL(4) STARTTIME(2), ENDTIME(2) 
  REAL(4) T, T2, ETIME
  
  EXTERNAL ETIME
  
  REAL(KIND=RKIND) dummy, dummyB
  INTEGER :: IErr = 0
  
  !-------------------------------------------------------------------
  ! constants
  !-------------------------------------------------------------------

  CALL Init_Numbers
  
  !--------------------------------------------------------------------
  ! protocal feature startup
  !--------------------------------------------------------------------
  
  PRINT*,"TMSEXd ASYMP (double precision)"
  !      PRINT*,"$Revision: 1.4 $ $Date: 2016/10/03 19:32:32 $ $Author: phsht $"
  PRINT*,RStr,DStr,AStr
  !      PRINT*,"$Revision: 1.4 $ $Date: 2016/10/03 19:32:32 $ $Author: phsht $"
  PRINT*,"--------------------------------------------------------------"
  
  !--------------------------------------------------------------------
  ! DBG mode
  !--------------------------------------------------------------------
  
  !MS$IF DEFINED (DBG)
  !PRINT*,"DBG mode"
  !MS$ENDIF
  
  !--------------------------------------------------------------------
  ! input handling
  !--------------------------------------------------------------------
  
  CALL Input( IErr )
  !PRINT*, "DBG: IErr=", IErr
  IF( IErr.NE.0 ) THEN
     PRINT*,"main: error in Input()"
     STOP
  ENDIF
  
  IF (ISortFlag.EQ.1) THEN
     PRINT*,"ReSort() is done after each renormalization !"
     PRINT*,"-------------------------------------------------------"
  ENDIF
  
!!$  !--------------------------------------------------------------------
!!$  ! making sure that D=1 is handled correctly
!!$  !--------------------------------------------------------------------
!!$
!!$  IF (IDimenFlag.EQ.1) THEN
!!$     IF(Width0.NE.1 .OR. Width1.NE.1) THEN
!!$        PRINT*,"main: ERR, Width0=",Width0, &
!!$             " & Width1=", Width1, &
!!$             " not permissible when IDimenFlag=", IDimenFlag
!!$        IErr=IErr+1
!!$     ENDIF
!!$     IF(IBCFlag.GT.0) THEN
!!$        PRINT*,"main: ERR, IBCGlag=",IBCFlag, &
!!$             " not permissible when IDimenFlag=", IDimenFlag
!!$        IErr=IErr+1
!!$     ENDIF
!!$     IF(IERR.GT.0) STOP 
!!$  ENDIF
  
  !--------------------------------------------------------------------
  ! select the quantity to be varied
  !--------------------------------------------------------------------
  
  IF (IFluxFlag.EQ.0) THEN
     flux0     = DiagDis0
     flux1     = DiagDis1
     dflux     = dDiagDis
     Energy    = Energy0
  ELSE
     flux0     = Energy0
     flux1     = Energy1
     dflux     = dEnergy
     DiagDis   = DiagDis0
  ENDIF
  
  !--------------------------------------------------------------------
  ! reload last saved parameters
  !--------------------------------------------------------------------
  
  SELECT CASE(IKeepFlag)
  CASE(3)
     CALL ReloadCurrentParameters( IWidthRL, fluxRL, IErr)
     SELECT CASE(IErr)
     CASE(0)
        PRINT*,"tmseXdXX: running with reloaded parameters !"
        Width0_X = IWidthRL
        flux0_X  = fluxRL
     CASE(1)
        PRINT*,"tmseXdXX: reloading of parameters not possible --- full rerun mode !"
        Width0_X = Width0
        flux0_X  = flux0
        PRINT*,"tmseXdXX: IKeepFlag=1 choosen."
        IKeepFlag= 1
     END SELECT
  CASE(2)
     CALL ReloadCurrentParameters( IWidthRL, fluxRL, IErr)
     SELECT CASE(IErr)
     CASE(0)
        PRINT*,"tmseXdXX: running with reloaded parameters !"
        Width0_X = IWidthRL
        flux0_X  = fluxRL
     CASE(1)
        PRINT*,"tmseXdXX: reloading of parameters not possible --- full rerun mode !"
        Width0_X = Width0
        flux0_X  = flux0
        PRINT*,"tmseXdXX: IKeepFlag=0 choosen."
        IKeepFlag= 0
     END SELECT
  CASE DEFAULT
     Width0_X = Width0
     flux0_X  = flux0
  END SELECT
  
  !--------------------------------------------------------------------
  ! main parameter sweep
  !--------------------------------------------------------------------
  
  width_loop: &
  DO IWidth= Width0_X,Width1,dWidth

     ! --------------------------------------------------     
     ! get time at start of the process
     ! --------------------------------------------------
     T = ETIME(STARTTIME)
  
     !--------------------------------------------------------------
     ! the # Lyapunov exponents is maximally .EQ. to IWidth
     !--------------------------------------------------------------
     
     NOfG= MIN( NOfGamma, IWidth )
     
     !--------------------------------------------------------------
     ! Calculate IWidthSquared if doing 3D TMM
     !--------------------------------------------------------------

     SELECT CASE(IDimenFlag)
     CASE(21,22)
        CONTINUE
     CASE(31,32)
        IWidthSquared = IWidth*IWidth
     CASE DEFAULT
        PRINT*,"tmseLMxD: IDimenFlag=", IDimenFlag, " not yet implemented --- aborting!"
        STOP
     END SELECT
     
     !--------------------------------------------------------------
     ! open the AVG files
     !--------------------------------------------------------------

     IF (IFluxFlag.EQ.0) THEN
        fluxStr = "E"
        fluxVal = Energy0
     ELSE
        fluxStr = "D"
        fluxVal = DiagDis0
     ENDIF
     filename= FileNameAvg(IWidth,fluxStr,fluxVal)
     !PRINT*,TRIM(filename)
     
     SELECT CASE(IKeepFlag)
     CASE(1)
        CALL CheckOutputAvg( filename, IErr )
        ! PRINT*, "DBG: IErr=", IErr
        IF( IErr.EQ.1 ) THEN
           PRINT*,"main: error in CheckOutputAvg()"
           STOP
        ELSE IF (IErr.EQ.2 ) THEN
           CYCLE width_loop
        ENDIF
        CALL OpenOutputAvg( filename, IWidth, IErr )
        ! PRINT*, "DBG: IErr=", IErr
        IF( IErr.EQ.1 ) THEN
           PRINT*,"main: error in OpenOutputAvg()"
           STOP
        ENDIF
     CASE(2,3)
        CALL ReOpenOutputAvg( filename, IErr )
        ! PRINT*, "DBG: IErr=", IErr
        IF( IErr.EQ.1 ) THEN
           PRINT*,"main: error in ReOpenOutputAvg()"
           STOP
        ENDIF
     CASE DEFAULT
        CALL OpenOutputAvg( filename, IWidth, IErr )
        ! PRINT*, "DBG: IErr=", IErr
        IF( IErr.EQ.1 ) THEN
           PRINT*,"main: error in OpenOutputAvg()"
           STOP
        ENDIF
     END SELECT
     
     !--------------------------------------------------------------
     ! allocate the memory for the TMM, gamma and RND vectors
     !--------------------------------------------------------------

     SELECT CASE(IDimenFlag)
     CASE(21,22)
        ALLOCATE(PsiA(IWidth,IWidth), STAT = IErr)
        ALLOCATE(PsiB(IWidth,IWidth), STAT = IErr)
        
        ALLOCATE(nGamma(IWidth), STAT = IErr)
     
        ALLOCATE(gamma(IWidth), STAT = IErr)
        ALLOCATE(gamma2(IWidth), STAT = IErr)
        ALLOCATE(acc_variance(IWidth), STAT = IErr)
     CASE(31,32)
        ALLOCATE(PsiA(IWidthSquared,IWidthSquared), STAT = IErr)
        ALLOCATE(PsiB(IWidthSquared,IWidthSquared), STAT = IErr)
        
        ALLOCATE(nGamma(IWidthSquared), STAT = IErr)
     
        ALLOCATE(gamma(IWidthSquared), STAT = IErr)
        ALLOCATE(gamma2(IWidthSquared), STAT = IErr)
        ALLOCATE(acc_variance(IWidthSquared), STAT = IErr)
     CASE DEFAULT
        PRINT*,"tmseLMxD: ERR, IDimenFlag=", IDimenFLag, " is not implemented --- aborting!"
        STOP
     END SELECT
        
!!$     SELECT CASE(IDimenFlag)
!!$     !3D case
!!$     CASE(3)
!!$        ALLOCATE(PsiA(IWidthSquared,IWidthSquared), STAT = IErr)
!!$        ALLOCATE(PsiB(IWidthSquared,IWidthSquared), STAT = IErr)
!!$        
!!$        ALLOCATE(nGamma(IWidthSquared), STAT = IErr)
!!$     
!!$        ALLOCATE(gamma(IWidthSquared), STAT = IErr)
!!$        ALLOCATE(gamma2(IWidthSquared), STAT = IErr)
!!$        ALLOCATE(acc_variance(IWidthSquared), STAT = IErr)
!!$
!!$     !2D/1D case
!!$     CASE DEFAULT
!!$        ALLOCATE(PsiA(IWidth,IWidth), STAT = IErr)
!!$        ALLOCATE(PsiB(IWidth,IWidth), STAT = IErr)
!!$        
!!$        ALLOCATE(nGamma(IWidth), STAT = IErr)
!!$     
!!$        ALLOCATE(gamma(IWidth), STAT = IErr)
!!$        ALLOCATE(gamma2(IWidth), STAT = IErr)
!!$        ALLOCATE(acc_variance(IWidth), STAT = IErr)
!!$     END SELECT
     
     !PRINT*, "DBG: IErr=", IErr
     IF( IErr.NE.0 ) THEN
        PRINT*,"main: error in ALLOCATE()"
        STOP
     ENDIF
     
     !--------------------------------------------------------------
     ! flux loop
     !--------------------------------------------------------------
     
flux_loop: &
     DO flux= flux0_X,flux1,dflux

        !--------------------------------------------------------------
        ! set values for the physical quantities
        !--------------------------------------------------------------

        IF (IFluxFlag.EQ.0) THEN
           DiagDis   = flux
        ELSE
           Energy    = flux
        ENDIF
        
        !--------------------------------------------------------------
        ! save the current parameters
        !--------------------------------------------------------------
        
        ! Segmentation fault here
        CALL SaveCurrentParameters(IWidth, flux, IErr)
        IF(IErr.NE.0) THEN
           PRINT*,"main: ERR in SaveCurrentParameters()"
        ENDIF
        
        !--------------------------------------------------------------
        ! open the nGamma file
        !--------------------------------------------------------------
        
        IF(IWriteFlag.GE.3) THEN
           CALL OpenOutputGamma( IWidth, DiagDis,Energy, IErr )
           !PRINT*, "DBG: IErr=", IErr
           IF( IErr.NE.0 ) THEN
              PRINT*,"main: error in OpenOutputGamma()"
              STOP
           ENDIF
        ENDIF

        !--------------------------------------------------------------
        ! protocoll feature
        !--------------------------------------------------------------
        
2500    WRITE(*,2510) IWidth, DiagDis,Energy,Kappa
2510    FORMAT("START @ IW= ",I4.1,    &
             ", DD= ", G10.3,          &
             ", En= ", G10.3,          &
             ", Ka= ", G10.3)
        
        !--------------------------------------------------------------
        ! initialize the wave vectors and the gamma sums
        !--------------------------------------------------------------
        
        ! reset the wave vectors
        PsiA        = ZERO
        PsiB        = ZERO
        !reset the gamma sum approximants for a SINGLE
        !configuration
        gamma       = ZERO
        gamma2      = ZERO
        acc_variance= ZERO
        
        SELECT CASE(IDimenFlag)
        CASE(21,22)
           DO index=1,IWidth
              PsiA(index,index)= ONE
           ENDDO
        CASE(31,32)
           DO index=1,IWidthSquared
              PsiA(index,index)= ONE
           ENDDO        
        CASE DEFAULT
           PRINT*,"tmseLMxD: ERR, IDimenFlag=", IDimenFLag, " is not implemented --- aborting!"
           STOP           
        END SELECT
        
        !PRINT *, 'DBG: PsiA=', PsiA
        !PRINT *, 'DBG: PsiB=', PsiB
        
        !--------------------------------------------------------------
        ! initialize the random number generator. NOTE that this is
        ! done with the same ISeed for each value of IWidth and
        ! Disorder in order to make these runs INDIVIDUALLY reprodu-
        ! cible
        !--------------------------------------------------------------
        
        CALL SRANDOM(ISeed)
        
        !--------------------------------------------------------------
        ! select right size of step in NOrtho loop
        !--------------------------------------------------------------
        SELECT CASE(IDimenFlag)
        CASE(21,31)
           iStep=2
        CASE (22,32)
           iStep=3
        CASE DEFAULT
           PRINT*,"tmseLMxD: ERR, IDimenFlag=", IDimenFLag, " is not implemented --- aborting!"
           STOP
        END SELECT

        IF(MOD(NOfOrtho,iStep).NE.0) THEN
           PRINT*,"tmseLMxD: WRNG, iStep=", iStep, " is not integer divisor of NOrtho=", NOfOrtho
           NOfOrtho= MAX(1,(NOfOrtho/iStep)*iStep)
           PRINT*,"tmseLMxD: NOrtho=", NOfOrtho, " taken as new value."
        END IF

        !--------------------------------------------------------------
        ! iteration loop
        !--------------------------------------------------------------

tmm_loop: &
        DO Iter1= 1, MAX( (NOfIter)/(NOfOrtho), 1)
        
northo_loop: &
           DO Iter2= 1, NOfOrtho, iStep

              ! Ilayer is the current horizontal position
              Ilayer= (Iter1-1)*NOfOrtho+Iter2

              ! do the TM multiplication
              SELECT CASE(IDimenFlag)
              CASE(21)
                 CALL TMMultLieb2DAtoB( PsiA, PsiB, Ilayer, &
                      Energy, DiagDis, IWidth)              
                 CALL TMMultLieb2DBtoA( PsiB, PsiA, Ilayer+1, &
                      Energy, DiagDis, IWidth)
              CASE (22)
                 CALL TMMultLieb2DAtoB1( PsiA, PsiB, Ilayer, &
                      Energy, DiagDis, IWidth)              
                 CALL TMMultLieb2DB1toB2( PsiB, PsiA, Ilayer+1, &
                      Energy, DiagDis, IWidth)
                 CALL TMMultLieb2DB2toA( PsiA, PsiB, Ilayer+2, &
                      Energy, DiagDis, IWidth)
                 CALL Swap( PsiA, PsiB, IWidth)
              CASE(31)
                 PRINT*,"DBG: WRNG! --- TMMultLieb3DAtoB/BtoA() not yet implemented, using old TMMult3D()"
                 CALL TMMultLieb3DAtoB( PsiA, PsiB, Ilayer, &
                      Energy, DiagDis, IWidth)
                 CALL TMMultLieb3DBtoA( PsiB, PsiA, Ilayer+1, &
                      Energy, DiagDis, IWidth)
              CASE(32)
                 PRINT*,"DBG: WRNG! --- TMMultLieb3DAtoB/BtoA() not yet implemented, using old TMMult3D()"
                 CALL TMMultLieb3DAtoB( PsiA, PsiB, Ilayer, &
                      Energy, DiagDis, IWidth)
                 CALL TMMultLieb3DBtoA( PsiB, PsiA, Ilayer+1, &
                      Energy, DiagDis, IWidth)
!!$                 CALL TMMultLieb3DAtoB1( PsiA, PsiB, Ilayer, &
!!$                      Energy, DiagDis, IWidth)              
!!$                 CALL TMMultLieb3DB1toB2( PsiB, PsiA, Ilayer+1, &
!!$                      Energy, DiagDis, IWidth)
!!$                 CALL TMMultLieb3DB2toA( PsiA, PsiB, Ilayer+2, &
!!$                      Energy, DiagDis, IWidth)
!!$                 CALL Swap( PsiA, PsiB, IWidth)
              CASE DEFAULT
                 PRINT*,"tmseLMxD: ERR, IDimenFlag=", IDimenFLag, " is not implemented --- aborting!"
                 STOP
              END SELECT
              
              !CALL Swap( PsiA, PsiB, IWidth)
              
           ENDDO northo_loop
            
           !-------------------------------------------------------
           ! renormalize via Gram-Schmidt
           !-------------------------------------------------------
           SELECT CASE(IDimenFlag)
           CASE(21,22)
              CALL ReNorm(PsiA,PsiB,gamma,gamma2,IWidth)
           CASE(31,32)
              CALL ReNorm(PsiA,PsiB,gamma,gamma2,IWidthSquared)
           CASE DEFAULT
              PRINT*,"tmseLMxD: IDimenFlag=", IDimenFlag, " not yet implemented --- aborting!"
              STOP
           END SELECT
           
!!$              CALL WriteOutputPsi( Iter1, IWidth, &
!!$                    PsiA, IErr)
!!$                 !PRINT*, "DBG: IErr=", IErr
!!$                 IF( IErr.NE.0 ) THEN
!!$                     PRINT*,"main: error in WriteOutputAvg()"
!!$                  STOP
!!$                 ENDIF 
           
           IF(IWriteFlag.GE.MAXWriteFLAG+1) THEN   
              PRINT*,"DBG: Orth,iL,PsiA", iLayer, PsiA
              PRINT*,"DBG: Orth,iL,PsiB", iLayer, PsiB; !PAUSE
           ENDIF

!!$           !-------------------------------------------------------
!!$           ! sort the eigenvalues by LARGEST first AND also
!!$           ! resort the eigenvectors accordingly
!!$           !-------------------------------------------------------           
!!$           
!!$           IF (ISortFlag.EQ.1) THEN
!!$           
!!$              SELECT CASE(IDimenFlag)  
!!$              !3D case
!!$              CASE(3)
!!$                 CALL ReSort(PsiA,PsiB,gamma,gamma2,IWidthSquared)
!!$              !2D/1D case
!!$              CASE DEFAULT
!!$                 CALL ReSort(PsiA,PsiB,gamma,gamma2,IWidth)
!!$              END SELECT           
!!$              
!!$           ENDIF
           
           !--------------------------------------------------------
           ! "Iter1" counts the number of actual renormalizations
           ! of the transfer matrix.
           !--------------------------------------------------------
           
           !-----------------------------------------------------------
           ! do the gamma computations
           !-----------------------------------------------------------
           
           SELECT CASE(IDimenFlag)
           CASE(21,22)
              DO iG=1, IWidth
                
                 nGamma(IWidth+1-iG)= gamma(iG)/REAL(NOfOrtho*Iter1)
                
                 acc_variance(IWidth+1-iG)=             &
                      SQRT( ABS(                        &
                      (gamma2(iG)/REAL(Iter1) -         &
                      (gamma(iG)/REAL(Iter1))**2 )      &
                      / REAL( MAX(Iter1-1,1) )          &
                      )) / ABS( gamma(iG)/REAL(Iter1) )
              ENDDO           
           CASE(31,32)              
              DO iG=1, IWidthSquared
                
                 nGamma(IWidthSquared+1-iG)= gamma(iG)/REAL(NOfOrtho*Iter1)
                
                 acc_variance(IWidthSquared+1-iG)=             &
                      SQRT( ABS(                        &
                      (gamma2(iG)/REAL(Iter1) -         &
                      (gamma(iG)/REAL(Iter1))**2 )      &
                      / REAL( MAX(Iter1-1,1) )          &
                      )) / ABS( gamma(iG)/REAL(Iter1) )
              ENDDO           
           CASE DEFAULT
              PRINT*,"tmseLMxD: IDimenFlag=", IDimenFlag, " not yet implemented --- aborting!"
              STOP
           END SELECT
           
           !-----------------------------------------------------------
           ! write the nGamma data to file
           !-----------------------------------------------------------
           
           IF(IWriteFlag.GE.3) THEN
              CALL WriteOutputGamma( Iter1, & 
                   nGamma, acc_variance, & 
                   NOfG, IErr )
              !PRINT*, "DBG: IErr=", IErr
              IF( IErr.NE.0 ) THEN
                 PRINT*,"main: error in WriteOutputGamma()"
                 STOP
              ENDIF
           ENDIF
            
           !--------------------------------------------------------
           ! write the nGamma data to stdout
           !--------------------------------------------------------
           
           IF(IWriteFlag.GE.1 .AND. MOD(Iter1,NOfPrint).EQ.0 ) THEN
              WRITE(*,4210) Iter1, nGamma(1), acc_variance(1)
4210          FORMAT(I7.1, G15.7, G15.7)
           ENDIF
           !PAUSE
           
           ! HP F77 compiler
           !CALL fortranfflush()
           
           !-----------------------------------------------------------
           ! check accuracy and dump the result
           !-----------------------------------------------------------
           
           IF( Iter1.GE.IWidth  .AND. &
                Iter1.GE.MINIter ) THEN
              
              IF( acc_variance(1).LE.Epsilon .AND. &
                   acc_variance(1).GE.TINY) THEN
                 IF(IWriteFlag.GE.4) PRINT*,"main: TMM_convergence in NOfIter-loop!"
                 TMM_CONVERGED=0
                 GOTO 4000
              ENDIF
              
           ENDIF
           
        ENDDO tmm_loop
        
        !-------------------------------------------------------------
        ! continue through here if convergence for a single configura-
        ! tion is NOT achieved, reset Iter1
        !-------------------------------------------------------------
        
        PRINT*,"NO convergence in NOfIter-loop:"
        TMM_CONVERGED=1
        Iter1= Iter1-1
        
        !-------------------------------------------------------------
        ! jump to this label if convergence for a single configura-
        ! tion is achieved.
        !-------------------------------------------------------------

4000    CONTINUE

        !--------------------------------------------------------------
        ! write the AVG data
        !--------------------------------------------------------------
        
        CALL WriteOutputAvg( IWidth, TMM_CONVERGED, &
             DiagDis,Energy, &
             nGamma, acc_variance, &
             NOfG, PsiA, IErr )
        !PRINT*, "DBG: IErr=", IErr
        IF( IErr.NE.0 ) THEN
           PRINT*,"main: error in WriteOutputAvg()"
           STOP
        ENDIF
        
        !-------------------------------------------------------------
        ! jump to this label if parabolic 
        !-------------------------------------------------------------
        
5000    WRITE(*,5010) Iter1, &
             DiagDis,Energy, Kappa
        WRITE(*,5012) nGamma(1), acc_variance(1)

5010    FORMAT("END @ ", I7.1, &
             ",", G15.7, ",", G15.7,",", G15.7)
5012    FORMAT("     ", &
             G15.7, ",", G15.7)
        
        !PAUSE

        !--------------------------------------------------------------
        ! close the nGamma file
        !--------------------------------------------------------------
        
        IF(IWriteFlag.GE.3) THEN
           CALL CloseOutputGamma( IErr )
           !PRINT*, "DBG: IErr=", IErr
           IF( IErr.NE.0 ) THEN
              PRINT*,"main: error in CloseOutputGamma()"
              STOP
           ENDIF
        ENDIF
        
        !--------------------------------------------------------------
        ! end of flux loop
        !--------------------------------------------------------------
        
     ENDDO flux_loop
     
     !-----------------------------------------------------------------
     ! close the AVG files
     !-----------------------------------------------------------------
     
     CALL CloseOutputAvg( IErr )
     !PRINT*, "DBG: IErr=", IErr
     IF( IErr.NE.0 ) THEN
        PRINT*,"main: error in CloseOutputAvg()"
        STOP
     ENDIF
     
     !--------------------------------------------------------------
     ! delete the current parameters
     !--------------------------------------------------------------
     
     CALL DeleteCurrentParameters(IWidth, IErr)
     IF(IErr.NE.0) THEN
        PRINT*,"main: ERR in DeleteCurrentParameters()"
     ENDIF
     
     !--------------------------------------------------------------
     ! DEallocate the memory for the TMM, gamma and RND vectors
     !--------------------------------------------------------------
     
     DEALLOCATE(PsiA,PsiB,nGamma,gamma,gamma2,acc_variance)
     
     ! --------------------------------------------------
     ! get time at the end of the process
     ! --------------------------------------------------
     
     T2 = ETIME(ENDTIME) 
     
     TIME(1,1) = T2 -T
     TIME(2,1) = ENDTIME(1)-STARTTIME(1)
     TIME(3,1) = ENDTIME(2)-STARTTIME(2)
     
     !IF(IWriteFlag.GE.2) THEN
        T  = TIME(1,1)
        T2 = TIME(2,1)
        
        WRITE(*,'(A39,6F12.4)') &
             "SINGLE PROC --> TIME(USR,SYS,DIFF): ", &
             & TIME(1,1), TIME(2,1), TIME(3,1)
        
     !ENDIF
     
     !-----------------------------------------------------------------
     ! end of width loop
     !-----------------------------------------------------------------
     
  ENDDO width_loop
  
  STOP "TMSEXd $Revision: 1.4 $"
  
END PROGRAM TMSEXd
