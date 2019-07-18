! --------------------------------------------------------------------
! TMMultLieb3DAtoB:
!
! 3D version of TMMult2D. Extra boundary conditions

SUBROUTINE TMMultLieb3DAtoB(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

  USE MyNumbers
  USE IPara
  USE RNG
  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis,&! diagonal disorder
       En                    ! energy

  REAL(KIND=RKIND) PSI_A(M*M,M*M),PSI_B(M*M,M*M),OnsitePotVec(2*M,2*M)
  
  INTEGER jState, ISeedDummy, xSiteS,ySiteS, xSiteL,ySiteL, indexK
  REAL(KIND=RKIND) OnsitePot, OnsiteRight, OnsiteLeft, OnsiteUp, OnsiteDown
  REAL(KIND=RKIND) new, PsiLeft, PsiRight, PsiUp, PsiDown, stub

  INTEGER, PARAMETER :: LiebSpacer=2

  INTEGER Co2InL31
  EXTERNAL Co2InL31

  !PRINT*,"DBG: TMMultLieb3DAtoB()"
     
  ! create the new onsite potential
  DO xSiteS=1,2*M
     DO ySiteS=1,2*M
        SELECT CASE(IRNGFlag)
        CASE(0)
           OnsitePotVec(xSiteS,ySiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
        CASE(1)
           OnsitePotVec(xSiteS,ySiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
        CASE(2)
           OnsitePotVec(xSiteS,ySiteS)= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
        END SELECT

        stub= OnsitePotVec(xSiteS,ySiteS)
        IF( ABS(stub).LT.TINY) OnsitePotVec(xSiteS,ySiteS)= SIGN(TINY,stub)
     END DO
  END DO
  
  !PRINT*,"iS,pL,RndVec", xSite,pLevel,RndVec((pLevel-1)*M+xSite)

  !to the TMM

  DO xSiteL=1,M
     DO ySiteL=1,M
        
        xSiteS= (xSiteL-1)*LiebSpacer + 1
        ySiteS= (ySiteL-1)*LiebSpacer + 1

        indexK=(ySiteS/2)*M+(xSiteS+1)/2

        OnsitePot=OnsitePotVec(xSiteS,ySiteS)

       ! PRINT*,"DBG1: xSiteS, ySiteS, indexK", xSiteS, ySiteS, indexK
        
        DO jState=1,M*M
           
           !PsiLeft
           IF (xSiteL.LE.1) THEN !mod(indexK,M).EQ.1)
              SELECT CASE(IBCFlag)
              CASE(-1,0) ! hard wall BC
                 OnsiteLeft=ZERO      
                 PsiLeft=ZERO
              CASE(1) ! periodic BC
                 stub= OnsitePotVec(LiebSpacer*M,ySiteS)    
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteLeft=1.0D0 /stub
                 PsiLeft=PSI_A(Co2InL31(M,M,ySiteL),jState) /stub
              CASE(2) ! antiperiodic BC
                 stub= OnsitePotVec(LiebSpacer*M,ySiteS)    
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteLeft= 1.0D0 /stub
                 PsiLeft= -PSI_A(Co2InL31(M,M,ySiteL),jState) /stub
              CASE DEFAULT
                 PRINT*,"TMMultLieb2DAtoB(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              stub= OnsitePotVec(xSiteS-1,ySiteS)
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteLeft=1.0D0/stub
              PsiLeft=PSI_A(Co2InL31(M,xSiteL-1,ySiteL),jState) /stub
           END IF
           
           !PsiRight
           IF (xSiteL.GE.M) THEN !(mod(indexK,M).EQ.0) 
              SELECT CASE(IBCFLag)
              CASE(-1) ! hard wall BC with stubs
                 stub= OnsitePotVec(xSiteS+1,ySiteS)  
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteRight=1.0D0 /stub
                 PsiRight= ZERO
              CASE(0) ! hard wall BC
                 OnsiteRight= ZERO
                 PsiRight= ZERO
              CASE(1) ! periodic BC
                 stub= OnsitePotVec(xSiteS+1,ySiteS)  
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteRight= 1.0D0 /stub
                 PsiRight= PSI_A(Co2InL31(M,1,ySiteL),jState) /stub
              CASE(2) ! antiperiodic BC
                 stub= OnsitePotVec(xSiteS+1,ySiteS)  
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteRight= 1.0D0 /stub
                 PsiRight= -PSI_A(Co2InL31(M,1,ySiteL),jState) /stub
              CASE DEFAULT
                 PRINT*,"TMMultLieb2DAtoB(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              stub= OnsitePotVec(xSiteS+1,ySiteS)
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteRight=1.0D0 /stub
              PsiRight= PSI_A(Co2InL31(M,xSiteL+1,ySiteL),jState) /stub
           END IF

           !PsiDown
           IF (ySiteS.GE.M) THEN !(indexK>M*(M-1))
              SELECT CASE(IBCFlag)
              CASE(-1) ! hard wall BC with stubs
                 stub= OnsitePotVec(xSiteS,ySiteS+1)    
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteDown= 1.0D0 /stub
                 PsiDown= ZERO
              CASE(0) ! hard wall BC
                 OnsiteDown= ZERO
                 PsiDown= ZERO
              CASE(1) ! periodic BC
                 stub= OnsitePotVec(xSiteS,ySiteS+1)    
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteDown= 1.0D0 /stub
                 !   PsiDown= 1.0D0/OnsitePotVec(xSiteS,ySiteS+1)*PSI_A(jState,mod(Co2InL31(M,xSiteL,ySiteL),M))
                 PsiDown= PSI_A(Co2InL31(M,xSiteL,1),jState) /stub
              CASE(2) ! antiperiodic BC
                 stub= OnsitePotVec(xSiteS,ySiteS+1)    
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteDown= 1.0D0 /stub
                 !   PsiDown= 1.0D0/OnsitePotVec(xSiteS,ySiteS+1)*PSI_A(jState,mod(indexK,M))
                 PsiDown=-PSI_A(Co2InL31(M,xSiteL,1),jState) /stub
              CASE DEFAULT
                 PRINT*,"TMMultLieb2DAtoB(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              stub= OnsitePotVec(xSiteS,ySiteS+1)
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteDown= 1.0D0 /stub
              PsiDown= PSI_A(Co2InL31(M,xSiteL,ySiteL+1),jState) /stub
           END IF
           
           !PsiUp
           IF (ySiteL.LE.1) THEN !(indexK<=M)
              SELECT CASE(IBCFlag)
              CASE(-1,0) ! hard wall BC
                 PsiUp= ZERO            
                 OnsiteUp= ZERO
              CASE(1) ! periodic BC
                 stub= OnsitePotVec(xSiteS,LiebSpacer*M)
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteUp= 1.0D0 /stub
                 PsiUp= PSI_A(Co2InL31(M,xSiteL,M),jState) /stub
              CASE(2) ! antiperiodic BC
                 stub= OnsitePotVec(xSiteS,LiebSpacer*M)
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteUp= 1.0D0 /stub
                 PsiUp= -PSI_A(Co2InL31(M,xSiteL,M),jState) /stub
              CASE DEFAULT
                 PRINT*,"TMMultLieb2DAtoB(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              stub= OnsitePotVec(xSiteS,ySiteS-1)
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteUp= 1.0D0 /stub
              PsiUp= PSI_A(Co2InL31(M,xSiteL,ySiteL-1),jState) /stub
           END IF
           
           !PRINT*,"DBG2: jState,xSiteS, ySiteS, indexK", jState, xSiteS, ySiteS, indexK
           new= ( OnsitePot - OnsiteLeft - OnsiteRight - OnsiteUp - OnsiteDown ) * PSI_A(Co2InL31(M,xSiteL,ySiteL),jState)&
                - Kappa * ( PsiLeft + PsiRight + PsiUp + PsiDown  ) &
                - PSI_B(Co2InL31(M,xSiteL,ySiteL),jState) 
           
           !PRINT*,"xSiteS,ySiteS,En, OP, PL, PR, PA,PB, PN"
           !PRINT*, xSiteS, jState, En, OnsitePot, PsiLeft, PsiRight,
           !        PSI_A(xSiteL,jState), PSI_B(xSite,jState),
           !        new
           
           PSI_B(Co2InL31(M,xSiteL,ySiteL),jState)= new
           
        ENDDO !jState
        
     ENDDO ! xSiteL
  ENDDO !ySiteL
  
  RETURN

END SUBROUTINE TMMultLieb3DAtoB

! --------------------------------------------------------------------
! convert i,j coordinates to an index
! used in lieb31
FUNCTION Co2InL31(M, xSiteS, ySiteS)
  INTEGER Co2InL31, M, xSiteS, ySiteS

  !old: indexK=(jSite/2)*M+(iSite+1)/2
  Co2InL31= (ySiteS-1)*M + xSiteS

  
  RETURN
END FUNCTION Co2InL31
  
! --------------------------------------------------------------------
! TMMultLieb3DBtoA:
!
! 3D version of TMMult2D. Extra boundary conditions

SUBROUTINE TMMultLieb3DBtoA(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

  USE MyNumbers
  USE IPara
  USE RNG
  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis,&! diagonal disorder
       En                    ! energy
  
  REAL(KIND=CKIND) PSI_A(M*M,M*M), PSI_B(M*M,M*M)
  
  INTEGER iSite, jState, ISeedDummy
  REAL(KIND=RKIND) OnsitePot
  REAL(KIND=CKIND) new
  
  !PRINT*,"DBG: TMMultLieb3DBtoA()"
  
  DO iSite=1,M*M
     
     ! create the new onsite potential
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
     CASE(1)
        OnsitePot= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
     CASE(2)
        OnsitePot= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
     
     !PRINT*,"iS,pL,RndVec", iSite,pLevel,RndVec((pLevel-1)*M+iSite)
     
     DO jState=1,M*M
        
        !PRINT*,"jState, iSite", jState, iSite,
        
        new= ( OnsitePot * PSI_A(iSite,jState) &
             - PSI_B(iSite,jState) )
        
        !PRINT*,"i,jSite,En, OP, PL, PR, PA,PB, PN"
        !PRINT*, iSite, jState, En, OnsitePot, PsiLeft, PsiRight,
        !        PSI_A(iSite,jState), PSI_B(iSite,jState),
        !        new
        
        PSI_B(iSite,jState)= new
        
     ENDDO ! jState
  ENDDO ! iSite
  
  RETURN
END SUBROUTINE TMMultLieb3DBtoA

