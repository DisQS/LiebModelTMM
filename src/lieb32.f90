! --------------------------------------------------------------------
! TMMultLieb3DAtoB:
!
! 3D version of TMMult2D. Extra boundary conditions

SUBROUTINE TMMultLieb3DAtoB5(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

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

  REAL(KIND=RKIND) PSI_A(M*M,M*M),PSI_B(M*M,M*M),OnsitePotVec(3*M,3*M)

  INTEGER jState, ISeedDummy,xSiteS,ySiteS, xSiteL,ySiteL
  REAL(KIND=RKIND) OnsitePot, OnsiteRight, OnsiteLeft, OnsiteUp, OnsiteDown
  REAL(KIND=RKIND) NEW, PsiLeft, PsiRight, PsiUp, PsiDown, stub

  INTEGER, PARAMETER :: LiebSpacer=3 

  INTEGER Co2InL32
  EXTERNAL Co2InL32

  !PRINT*,"DBG: TMMultLieb3DAtoB()"

  ! create the new onsite potential
  DO xSiteS=1,LiebSpacer*M
     DO ySiteS=1,LiebSpacer*M

        !indexS= (xSiteS-1)*LiebSpacer*M + ySiteS

        SELECT CASE(IRNGFlag)
        CASE(0)
           OnsitePotVec(xSiteS,ySiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
        CASE(1)
           OnsitePotVec(xSiteS,ySiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
        CASE(2)
           OnsitePotVec(xSiteS,ySiteS)= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
        END SELECT
     END DO
  END DO

  !PRINT*,"iS,pL,RndVec", xSite,pLevel,RndVec((pLevel-1)*M+xSite)

  !new TMM
  DO xSiteL=1,M
     DO ySiteL=1,M

        xSiteS= (xSiteL-1)*LiebSpacer + 1
        ySiteS= (ySiteL-1)*LiebSpacer + 1
        
        !indexL= (ySiteL-1)*M + xSiteL
        !indexS= (ySiteS-1)*M + xSiteS
        
!!$        PRINT*,"iSL,jSL, iSS, jSS, C2I", &
!!$             xSiteL,ySiteL, xSiteS,ySiteS, Co2InL32(M,xSiteL,ySiteL)

        OnsitePot=OnsitePotVec(xSiteS,ySiteS)

        DO jState=1,M*M
           
           !PsiLeft
           IF (xSiteL.LE.1) THEN
              SELECT CASE(IBCFLag)
              CASE(-1,0) ! hard wall BC
                 OnsiteLeft= ZERO
                 PsiLeft= ZERO
              CASE(1) ! periodic BC
                 stub= OnsitePotVec(LiebSpacer*M,ySiteS)*OnSitePotVec(LiebSpacer*M-1,ySiteS)-1.0D0
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteLeft= OnsitePotVec(LiebSpacer*M-1,ySiteS) /stub
                 PsiLeft= Psi_A(Co2InL32(M,M,ySiteL),jState) /stub
!              CASE(2) ! antiperiodic BC
              CASE DEFAULT
                 PRINT*,"TMMultLieb3DAtoB5(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              stub= OnsitePotVec(xSiteS-1,ySiteS)*OnSitePotVec(xSiteS-2,ySiteS)-1.0D0
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteLeft= OnsitePotVec(xSiteS-2,ySiteS) /stub
              PsiLeft= Psi_A(Co2InL32(M,xSiteL-1,ySiteL),jState) /stub
           END IF

           !PsiRight
           IF (xSiteL.GE.M) THEN
              SELECT CASE(IBCFLag)
              CASE(-1) ! hard wall BC + stubs
                 stub= OnsitePotVec(xSiteS+1,ySiteS)*OnSitePotVec(xSiteS+2,ySiteS)-1.0D0
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteRight= OnsitePotVec(xSiteS+2,ySiteS)/stub
                 PsiRight= ZERO
              CASE(0) ! hard wall
                 OnsiteRight= ZERO
                 PsiRight= ZERO
              CASE(1) ! periodic BC
                 stub= OnsitePotVec(xSiteS+1,ySiteS)*OnSitePotVec(xSiteS+2,ySiteS)-1.0D0
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteRight= OnsitePotVec(xSiteS+2,ySiteS) /stub
                 PsiRight= Psi_A(Co2InL32(M,1,ySiteL),jState) /stub
!              CASE(2) ! antiperiodic BC
              CASE DEFAULT
                 PRINT*,"TMMultLieb3DAtoB5(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              stub= OnsitePotVec(xSiteS+1,ySiteS)*OnSitePotVec(xSiteS+2,ySiteS)-1.0D0
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteRight= OnsitePotVec(xSiteS+2,ySiteS) /stub
              PsiRight= Psi_A(Co2InL32(M,xSiteL+1,ySiteL),jState) /stub
           END IF

           !PsiDown
           IF (ySiteL.GE.M) THEN
              SELECT CASE(IBCFlag)
              CASE(-1) ! hard wall BC + stubs
                 stub= OnsitePotVec(xSiteS,ySiteS+1)*OnSitePotVec(xSiteS,ySiteS+2)-1.0D0
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteDown= OnsitePotVec(xSiteS,ySiteS+2) /stub
                 PsiDown= ZERO
              CASE(0) ! hard wall
                 OnsiteUp=ZERO
                 PsiDown= ZERO
              CASE(1) ! periodic BC
                 stub= OnsitePotVec(xSiteS,ySiteS+1)*OnSitePotVec(xSiteS,ySiteS+2)-1.0D0
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteDown= OnsitePotVec(xSiteS,ySiteS+2) /stub
                 PsiDown= Psi_A(Co2InL32(M,xSiteL,1),jState) /stub
!              CASE(2) ! antiperiodic BC
              CASE DEFAULT
                 PRINT*,"TMMultLieb3DAtoB5(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              stub= OnsitePotVec(xSiteS,ySiteS+1)*OnSitePotVec(xSiteS,ySiteS+2)-1.0D0
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteDown= OnsitePotVec(xSiteS,ySiteS+2) /stub
              PsiDown= Psi_A(Co2InL32(M,xSiteL,ySiteL+1),jState) /stub
           END IF

           !PsiUp
           IF (ySiteL.LE.1) THEN
              SELECT CASE(IBCFlag)
              CASE(-1,0) ! hard wall BC
                 OnsiteUp= ZERO
                 PsiUp= ZERO
              CASE(1) ! periodic BC
                 stub= (OnsitePotVec(xSiteS,LiebSpacer*M)*OnSitePotVec(xSiteS,LiebSpacer*M-1)-1.0D0)
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteUp= OnsitePotVec(xSiteS,LiebSpacer*M-1) /stub
                 PsiUp=  Psi_A(Co2InL32(M,xSiteL,M),jState) /stub
              CASE(2) ! antiperiodic BC
              CASE DEFAULT
                 PRINT*,"TMMultLieb3DAtoB5(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              stub= OnsitePotVec(xSiteS,ySiteS-1)*OnSitePotVec(xSiteS,ySiteS-2)-1.0D0
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteUp= OnsitePotVec(xSiteS,ySiteS-2) /stub
              PsiUp= Psi_A(Co2InL32(M,xSiteL,ySiteL-1),jState) /stub
           END IF

           NEW= ( OnsitePot - OnsiteLeft - OnsiteRight - OnsiteUp - OnsiteDown ) * &
                Psi_A(Co2InL32(M,xSiteL,ySiteL),jState) &
                - Kappa * ( PsiLeft + PsiRight + PsiUp + PsiDown ) &
                - PSI_B(Co2InL32(M,xSiteL,ySiteL),jState)
           
           PSI_B(Co2InL32(M,xSiteL,ySiteL),jState)= NEW
        END DO !jState
        
     END DO !xSiteL
  END DO !ySiteL
  RETURN

END SUBROUTINE TMMultLieb3DAtoB5

! --------------------------------------------------------------------
! convert i,j coordinates to an index
! used in lieb32/33
FUNCTION Co2InL32(M, xSiteL, ySiteL)
  INTEGER Co2InL32, M, xSiteL, ySiteL
  
  Co2InL32= (ySiteL-1)*M + xSiteL
  
  RETURN
END FUNCTION Co2InL32

! --------------------------------------------------------------------
! TMMultLieb3DBtoA:
!
! 3D version of TMMult2D. Extra boundary conditions

SUBROUTINE TMMultLieb3DB5toB6(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

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
  REAL(KIND=CKIND) NEW
  
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
!!$        
!!$        NEW= ( OnsitePot * PSI_A(jState,iSite) &
!!$             - PSI_B(jState,iSite) )
        NEW=  OnsitePot * PSI_A(iSite,jState) &
             - PSI_B(iSite,jState) 
        
        !PRINT*,"i,jSite,En, OP, PL, PR, PA,PB, PN"
        !PRINT*, iSite, jState, En, OnsitePot, PsiLeft, PsiRight,
        !        PSI_A(iSite,jState), PSI_B(iSite,jState),
        !        new
        
!!$        PSI_B(jState,iSite)= NEW
        PSI_B(iSite,jState)= NEW
        
     ENDDO ! jState
  ENDDO ! iSite
  
  RETURN
END SUBROUTINE TMMultLieb3DB5toB6

SUBROUTINE TMMultLieb3DB6toA(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

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
  
  CALL TMMultLieb3DB5toB6(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

  RETURN
END SUBROUTINE TMMultLieb3DB6toA
