! --------------------------------------------------------------------
! TMMultLieb3DAtoB:
!
! 3D version of TMMult2D. Extra boundary conditions

SUBROUTINE TMMultLieb3D_AtoD1(PSI_A,PSI_B, Ilayer, En, DiagDis, RimDis, CubeConstPoten, LiebConstPoten, M )

  USE MyNumbers
  USE IPara
  USE RNG
!  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis,&! diagonal disorder
       RimDis, CubeConstPoten, LiebConstPoten, &
       En                    ! energy

  REAL(KIND=RKIND) PSI_A(M*M,M*M),PSI_B(M*M,M*M),OnsitePotVec(4*M,4*M)
  
  INTEGER xSiteL,ySiteL, xSiteS,ySiteS, jState, ISeedDummy
  REAL(KIND=RKIND) OnsitePot, OnsiteRight, OnsiteLeft, OnsiteUp, OnsiteDown
  REAL(KIND=RKIND) new, PsiLeft, PsiRight, PsiUp, PsiDown,stub, Kappa

  INTEGER, PARAMETER :: LiebSpacer=4

  INTEGER C2IL3
  EXTERNAL C2IL3

  !PRINT*,"DBG: TMMultLieb3DAtoB()"
     
  ! create the new onsite potential
!!$  DO xSiteS=1,LiebSpacer*M
!!$     DO ySiteS=1,LiebSpacer*M
!!$        SELECT CASE(IRNGFlag)
!!$        CASE(0)
!!$           OnsitePotVec(xSiteS,ySiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
!!$        CASE(1)
!!$           OnsitePotVec(xSiteS,ySiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
!!$        CASE(2)
!!$           OnsitePotVec(xSiteS,ySiteS)= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
!!$        END SELECT
!!$     END DO
!!$  END DO

  ! create the new onsite potential
  !  IRNGFlag=(xy)
  !  xy=0x      presenting conditions that all positions are disorder
  !  xy=x0,     presenting conditions that only central positions are disorder 
   
  DO xSiteS=1,LiebSpacer*M
     Do ySiteS=1,LiebSpacer*M
        SELECT CASE(IRNGFlag)
        CASE(01)
           OnsitePotVec(xSiteS,ySiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
        CASE(02)
           OnsitePotVec(xSiteS,ySiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
        CASE(03)
           OnsitePotVec(xSiteS,ySiteS)= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
        CASE(10)
           IF(Mod(xSiteS,LiebSpacer)==1 .AND. Mod(ySiteS,LiebSpacer)==1) THEN
              OnsitePotVec(xSiteS,ySiteS)= -En + CubeConstPoten + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
           ELSE
              OnsitePotVec(xSiteS,ySiteS)= -En + LiebConstPoten + RimDis*(DRANDOM(ISeedDummy)-0.5D0)
           END IF
        CASE(20)
           IF(Mod(xSiteS,LiebSpacer)==1 .AND. Mod(ySiteS,LiebSpacer)==1) THEN
              OnsitePotVec(xSiteS,ySiteS)= -En + CubeConstPoten + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
           ELSE
              OnsitePotVec(xSiteS,ySiteS)= -En + LiebConstPoten + RimDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
           END IF
        CASE(30)
           IF(Mod(xSiteS,LiebSpacer)==1 .AND. Mod(ySiteS,LiebSpacer)==1) THEN
              OnsitePotVec(xSiteS,ySiteS)= -En + CubeConstPoten + GRANDOM(ISeedDummy,0.0D0,DiagDis)
           ELSE
              OnsitePotVec(xSiteS,ySiteS)= -En + LiebConstPoten + GRANDOM(ISeedDummy,0.0D0,RimDis)
           END IF
        END SELECT
     END DO
  END DO
  
  !to the TMM
  DO xSiteL=1,M
     DO ySiteL=1,M

        xSiteS= (xSiteL-1)*LiebSpacer + 1
        ySiteS= (ySiteL-1)*LiebSpacer + 1

        OnsitePot=OnsitePotVec(xSiteS,ySiteS)
        
        DO jState=1,M*M
           
           !PsiLeft
           IF (xSiteL.LE.1) THEN
              SELECT CASE(IBCFlag)
              CASE(-1,0) ! hard wall BC
                 PsiLeft= ZERO           
                 OnsiteLeft= ZERO
              CASE(1) ! periodic BC
                 stub= ( OnsitePotVec(LiebSpacer*M-1,ySiteS)*&
                      OnsitePotVec(LiebSpacer*M-2,ySiteS)*OnsitePotVec(LiebSpacer*M,ySiteS) &
                      - OnsitePotVec(LiebSpacer*M,ySiteS)-OnsitePotVec(LiebSpacer*M-2,ySiteS) )
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteLeft=( OnsitePotVec(LiebSpacer*M-1,ySiteS)*OnsitePotVec(LiebSpacer*M-2,ySiteS)-1.0D0) /stub
                 PsiLeft=PSI_A(C2IL3(M,M,ySiteL),jState) /stub
              !CASE(2) ! antiperiodic BC
              CASE DEFAULT
                 PRINT*,"TMMultLieb3D_AtoD1(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              stub= ( OnsitePotVec(xSiteS-2,ySiteS)*OnsitePotVec(xSiteS-3,ySiteS)* &
                   OnsitePotVec(xSiteS-1,ySiteS) &
                   - OnsitePotVec(xSiteS-1,ySiteS)-OnsitePotVec(xSiteS-3,ySiteS) )
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteLeft=( OnsitePotVec(xSiteS-2,ySiteS)*OnsitePotVec(xSiteS-3,ySiteS)-1.0D0)/stub 
              PsiLeft=PSI_A(C2IL3(M,xSiteL-1,ySiteL),jState)/stub
           END IF
           
           !PsiRight
           IF (xSiteL.GE.M) THEN
              SELECT CASE(IBCFlag)
              CASE(-1) ! hard wall BC + stubs
                 stub=( OnsitePotVec(xSiteS+2,ySiteS)*OnsitePotVec(xSiteS+3,ySiteS)* &
                      OnsitePotVec(xSiteS+1,ySiteS) &
                      - OnsitePotVec(xSiteS+1,ySiteS)-OnsitePotVec(xSiteS+3,ySiteS) )
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteRight=( OnsitePotVec(xSiteS+2,ySiteS)*OnsitePotVec(xSiteS+3,ySiteS)-1.0D0)/stub  
                 PsiRight= ZERO
              CASE(0) ! hard wall
                 OnsiteRight= ZERO
                 PsiRight= ZERO
              CASE(1) ! periodic BC
                 stub=( OnsitePotVec(xSiteS+2,ySiteS)*OnsitePotVec(xSiteS+3,ySiteS)* &
                      OnsitePotVec(xSiteS+1,ySiteS) &
                      - OnsitePotVec(xSiteS+1,ySiteS)-OnsitePotVec(xSiteS+3,ySiteS) )
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteRight=( OnsitePotVec(xSiteS+2,ySiteS)*OnsitePotVec(xSiteS+3,ySiteS)-1.0D0)/ stub
                 PsiRight=PSI_A(C2IL3(M,1,ySiteL),jState)/ stub
              !CASE(2) ! antiperiodic BC 
              CASE DEFAULT
                 PRINT*,"TMMultLieb3D_AtoD1(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              stub=( OnsitePotVec(xSiteS+2,ySiteS)*OnsitePotVec(xSiteS+3,ySiteS)* &
                   OnsitePotVec(xSiteS+1,ySiteS) &
                   - OnsitePotVec(xSiteS+1,ySiteS)-OnsitePotVec(xSiteS+3,ySiteS) )
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteRight=( OnsitePotVec(xSiteS+2,ySiteS)*OnsitePotVec(xSiteS+3,ySiteS)-1.0D0)/ stub
              PsiRight=PSI_A(C2IL3(M,xSiteL+1,ySiteL),jState)/ stub
           END IF
           
           !PsiDown
           IF (ySiteL.LE.1) THEN
              SELECT CASE(IBCFLag)
              CASE(-1,0) ! hard wall BC
                 OnsiteDown=ZERO      
                 PsiDown=ZERO
              CASE(1) ! periodic BC
                 stub=( OnsitePotVec(xSiteS,LiebSpacer*M-1)*&
                      OnsitePotVec(xSiteS,LiebSpacer*M-2)*OnsitePotVec(xSiteS,LiebSpacer*M) &
                      -OnsitePotVec(xSiteS,LiebSpacer*M)-OnsitePotVec(xSiteS,LiebSpacer*M-2))
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteDown=(OnsitePotVec(xSiteS,LiebSpacer*M-1)*OnsitePotVec(xSiteS,LiebSpacer*M-2)-1.0D0)/stub
                 PsiDown=PSI_A(C2IL3(M,xSiteL,M),jState) /stub      
              !CASE(2) ! antiperiodic BC
              CASE DEFAULT
                 PRINT*,"TMMultLieb3D_AtoD1(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT                 
           ELSE
              stub=( OnsitePotVec(xSiteS,ySiteS-2)*OnsitePotVec(xSiteS,ySiteS-3)* &
                   OnsitePotVec(xSiteS,ySiteS-1) &
                   -OnsitePotVec(xSiteS,ySiteS-1)-OnsitePotVec(xSiteS,ySiteS-3))
              IF( ABS(stub).LT.TINY)stub= SIGN(TINY,stub)
              OnsiteDown=(OnsitePotVec(xSiteS,ySiteS-2)*OnsitePotVec(xSiteS,ySiteS-3)-1.0D0)/stub
              PsiDown=PSI_A(C2IL3(M,xSiteL,ySiteL-1),jState) /stub                   
           END IF
           
           !PsiUp
           IF (ySiteL.GE.M) THEN
              SELECT CASE(IBCFlag)
              CASE(-1) ! hard wall BC with stubs
                 stub= ( OnsitePotVec(xSiteS,ySiteS+2)*OnsitePotVec(xSiteS,ySiteS+3)* &
                      OnsitePotVec(xSiteS,ySiteS+1) &
                      -OnsitePotVec(xSiteS,ySiteS+1)-OnsitePotVec(xSiteS,ySiteS+3))
                 IF( ABS(stub).LT.TINY)stub= SIGN(TINY,stub)
                 OnsiteUp=(OnsitePotVec(xSiteS,ySiteS+2)*OnsitePotVec(xSiteS,ySiteS+3)-1.0D0) /stub    
                 PsiUp=ZERO
              CASE(0) ! hard wall BC
                 OnsiteUp=ZERO
                 PsiUp=ZERO
              CASE(1) ! periodic BC
                 stub= ( OnsitePotVec(xSiteS,ySiteS+2)*OnsitePotVec(xSiteS,ySiteS+3)* &
                      OnsitePotVec(xSiteS,ySiteS+1) &
                      -OnsitePotVec(xSiteS,ySiteS+1)-OnsitePotVec(xSiteS,ySiteS+3))
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteUp=(OnsitePotVec(xSiteS,ySiteS+2)*OnsitePotVec(xSiteS,ySiteS+3)-1.0D0 )/stub
                 PsiUp=PSI_A(C2IL3(M,xSiteL,1),jState)/stub
              !CASE(2) ! antiperiodic BC
              CASE DEFAULT
                 PRINT*,"TMMultLieb3D_AtoD1(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT                 
           ELSE
              stub= ( OnsitePotVec(xSiteS,ySiteS+2)*OnsitePotVec(xSiteS,ySiteS+3)* &
                   OnsitePotVec(xSiteS,ySiteS+1) &
                   -OnsitePotVec(xSiteS,ySiteS+1)-OnsitePotVec(xSiteS,ySiteS+3))
              IF( ABS(stub).LT.TINY)stub= SIGN(TINY,stub)
              OnsiteUp=(OnsitePotVec(xSiteS,ySiteS+2)*OnsitePotVec(xSiteS,ySiteS+3)-1.0D0 )/stub
              PsiUp=PSI_A(C2IL3(M,xSiteL,ySiteL+1),jState)/stub
           END IF

           new= ( OnsitePot - OnsiteLeft - OnsiteRight - OnsiteDown - OnsiteUp ) * &
                PSI_A(C2IL3(M,xSiteL,ySiteL),jState)&
                - Kappa * ( PsiLeft + PsiRight + PsiDown + PsiUp ) &
                - PSI_B(C2IL3(M,xSiteL,ySiteL),jState) 
           
           PSI_B(C2IL3(M,xSiteL,ySiteL),jState)= new
           
        ENDDO !jState
        
     ENDDO ! iSite
  ENDDO!jSite
  
  RETURN

END SUBROUTINE TMMultLieb3D_AtoD1

!!$! --------------------------------------------------------------------
!!$! convert i,j coordinates to an index
!!$! used in lieb32/33
!!$FUNCTION Co2InL33(M, xSiteL, ySiteL)
!!$  INTEGER Co2InL33, M, xSiteL, ySiteL
!!$
!!$  INTEGER Co2InL31
!!$  EXTERNAL Co2InL31
!!$
!!$  Co2InL33= Co2InL31(M,xSiteL,ySiteL)
!!$  
!!$  RETURN
!!$END FUNCTION Co2InL33
!!$
! --------------------------------------------------------------------
! TMMultLieb3DBtoA:
!
! 3D version of TMMult2D. Extra boundary conditions

SUBROUTINE TMMultLieb3D_D1toD2(PSI_A,PSI_B, Ilayer, En, DiagDis, RimDis, CubeConstPoten, LiebConstPoten, M )

  USE MyNumbers
  USE IPara
  USE RNG
!  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis,&! diagonal disorder
       RimDis, CubeConstPoten, LiebConstPoten, &
       En                    ! energy
  
  REAL(KIND=CKIND) PSI_A(M*M,M*M), PSI_B(M*M,M*M)
  
  INTEGER iSite, jState, ISeedDummy
  REAL(KIND=RKIND) OnsitePot
  REAL(KIND=CKIND) new
  
  !PRINT*,"DBG: TMMultLieb3DBtoA()"
  
  DO iSite=1,M*M
     
     ! create the new onsite potential
!!$     SELECT CASE(IRNGFlag)
!!$     CASE(0)
!!$        OnsitePot= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
!!$     CASE(1)
!!$        OnsitePot= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
!!$     CASE(2)
!!$        OnsitePot= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
!!$     END SELECT

     SELECT CASE(IRNGFlag)
     CASE(01)
        OnsitePot= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
     CASE(02)
        OnsitePot= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
     CASE(03)
        OnsitePot= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
     CASE(10)
        OnsitePot= -En + LiebConstPoten + RimDis*(DRANDOM(ISeedDummy)-0.5D0)
     CASE(20)
        OnsitePot= -En + LiebConstPoten + RimDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
     CASE(30)
        OnsitePot= -En + LiebConstPoten + GRANDOM(ISeedDummy,0.0D0,RimDis)
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
END SUBROUTINE TMMultLieb3D_D1toD2


SUBROUTINE TMMultLieb3D_D2toD3(PSI_A,PSI_B, Ilayer, En, DiagDis, RimDis, CubeConstPoten, LiebConstPoten, M )

  USE MyNumbers
  USE IPara
  USE RNG
!  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis,&! diagonal disorder
       RimDis, CubeConstPoten, LiebConstPoten, &
       En                    ! energy
  
  REAL(KIND=CKIND) PSI_A(M*M,M*M), PSI_B(M*M,M*M)
  
  CALL TMMultLieb3D_D1toD2(PSI_A,PSI_B, Ilayer, En, DiagDis, RimDis, CubeConstPoten, LiebConstPoten, M )

  RETURN
END SUBROUTINE TMMultLieb3D_D2toD3

SUBROUTINE TMMultLieb3D_D3toA(PSI_A,PSI_B, Ilayer, En, DiagDis, RimDis, CubeConstPoten, LiebConstPoten, M )

  USE MyNumbers
  USE IPara
  USE RNG
!  USE DPara
  
  ! wave functions:
  !       
  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
  
  IMPLICIT NONE
  
  INTEGER Ilayer,           &! current # TM multiplications
       M                     ! strip width
  
  REAL(KIND=RKIND)  DiagDis,&! diagonal disorder
       RimDis, CubeConstPoten, LiebConstPoten, &
       En                    ! energy
  
  REAL(KIND=CKIND) PSI_A(M*M,M*M), PSI_B(M*M,M*M)
  
  CALL TMMultLieb3D_D1toD2(PSI_A,PSI_B, Ilayer, En, DiagDis, RimDis, CubeConstPoten, LiebConstPoten, M )

  RETURN
END SUBROUTINE TMMultLieb3D_D3toA
