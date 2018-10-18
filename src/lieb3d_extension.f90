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

  INTEGER iSite,jSite,indexK,jstate, ISeedDummy
  REAL(KIND=RKIND) OnsitePot, OnsiteRight, OnsiteLeft, OnsiteUp, OnsiteDown
  REAL(KIND=RKIND) new, PsiLeft, PsiRight, PsiUp, PsiDown, stub

  !PRINT*,"DBG: TMMultLieb3DAtoB()"

  ! create the new onsite potential
  DO iSite=1,3*M
     DO jSite=1,3*M
        SELECT CASE(IRNGFlag)
        CASE(0)
           OnsitePotVec(iSite,jSite)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
        CASE(1)
           OnsitePotVec(iSite,jSite)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
        CASE(2)
           OnsitePotVec(iSite,jSite)= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
        END SELECT

        IF( ABS(OnsitePotVec(iSite,jSite)).LT.TINY) THEN
           OnsitePotVec(iSite,jSite)= SIGN(TINY,OnsitePotVec(iSite,jSite))
        END IF
     END DO
  END DO


  !PRINT*,"iS,pL,RndVec", iSite,pLevel,RndVec((pLevel-1)*M+iSite)

  !to the TMM
  DO jSite=1,3*M,3
     DO iSite=1,3*M,3

        indexK=(iSite/3+1)+(jSite/3)*M
        OnsitePot=OnsitePotVec(iSite,jSite)
        ! PRINT*,"DBG1: iSite, jSite, indexK", iSite, jSite, indexK

        DO jstate=1,M*M

           !PsiLeft
           IF (indexK<=M) THEN

              IF (IBCFlag.EQ.0) THEN
                 PsiLeft= ZERO            ! hard wall BC
                 OnsiteLeft= ZERO
              ELSE IF (IBCFlag.EQ.1) THEN
                 stub= (OnsitePotVec(iSite,3*M-1)*OnsitePotVec(iSite,3*M)-1.0D0)
                 IF( ABS(stub).LT.TINY) THEN
                    PRINT*,"DBG: iSite, jSite, jState, stub(OL,PBC)", &
                         iSite, jSite, JState, stub
                    stub= SIGN(TINY,stub)
                 ENDIF
                 OnsiteLeft=OnsitePotVec(iSite,3*M-1) &
                      /stub !(OnsitePotVec(iSite,3*M-1)*OnsitePotVec(iSite,3*M)-1.0D0)  ! periodic BC
                 PsiLeft=Psi_A(jstate,M*M-mod(indexK,M)) &
                      /stub !(OnsitePotVec(iSite,3*M-1)*OnsitePotVec(iSite,3*M)-1.0D0)
              ELSE IF (IBCFlag.EQ.2) THEN
                 stub= (OnsitePotVec(iSite,3*M-1)*OnsitePotVec(iSite,3*M)-1.0D0)
                 IF( ABS(stub).LT.TINY) THEN
                    PRINT*,"DBG: iSite, jSite, jState, stub(OL,APBC)", &
                         iSite, jSite, jState, stub
                    stub= SIGN(TINY,stub)
                 ENDIF
                 OnsiteLeft=OnsitePotVec(iSite,3*M-1) &
                      /stub !(OnsitePotVec(iSite,3*M-1)*OnsitePotVec(iSite,3*M)-1.0D0)   ! antiperiodic BC
                 PsiLeft=-Psi_A(jstate,M*M-mod(indexK,M)) &
                      /stub !(OnsitePotVec(iSite,3*M-1)*OnsitePotVec(iSite,3*M)-1.0D0)    
              ENDIF
           ELSE
              stub= (OnsitePotVec(iSite,jSite-1)*OnsitePotVec(iSite,jSite-2)-1.0D0)
              IF( ABS(stub).LT.TINY) THEN
                 PRINT*,"DBG: iSite, jSite, jState, stub(OL)", &
                      iSite, jSite, jState, stub
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteLeft=OnsitePotVec(iSite,jSite-2) &
                   /stub !(OnsitePotVec(iSite,jSite-1)*OnsitePotVec(iSite,jSite-2)-1.0D0)
              PsiLeft=Psi_A(jstate,indexK-M) &
                   /stub !(OnsitePotVec(iSite,jSite-1)*OnsitePotVec(iSite,jSite-2)-1.0D0)
           END IF

           !PsiRight
           IF (indexK>M*(M-1)) THEN

              IF (IBCFlag.EQ.0) THEN
                 stub= (OnsitePotVec(iSite,jSite+1)*OnsitePotVec(iSite,jSite+2)-1.0D0)
                 IF( ABS(stub).LT.TINY) THEN
                    PRINT*,"DBG: iSite, jSite, jState, stub(OR,HW)", &
                         iSite, jSite, jState, stub
                    stub= SIGN(TINY,stub)
                 ENDIF
                 OnsiteRight=OnsitePotVec(iSite,jSite+2) &
                      /stub !(OnsitePotVec(iSite,jSite+1)*OnsitePotVec(iSite,jSite+2)-1.0D0)   ! hard wall BC
                 PsiRight= ZERO        
              ELSE IF (IBCFlag.EQ.1) THEN
                 stub= (OnsitePotVec(iSite,jSite+1)*OnsitePotVec(iSite,jSite+2)-1.0D0)
                 IF( ABS(stub).LT.TINY) THEN
                    PRINT*,"DBG: iSite, jSite, jState, stub(OR,PBC)", &
                         iSite, jSite, jState, stub
                    stub= SIGN(TINY,stub)
                 ENDIF
                 OnsiteRight=OnsitePotVec(iSite,jSite+2) &
                      /stub !(OnsitePotVec(iSite,jSite+1)*OnsitePotVec(iSite,jSite+2)-1.0D0)  ! periodic BC
                 PsiRight=Psi_A(jstate,mod(indexK,M)) &
                      /stub !(OnsitePotVec(iSite,jSite+1)*OnsitePotVec(iSite,jSite+2)-1.0D0)
              ELSE IF (IBCFlag.EQ.2) THEN
                 stub= (OnsitePotVec(iSite,jSite+1)*OnsitePotVec(iSite,jSite+2)-1.0D0)
                 IF( ABS(stub).LT.TINY) THEN
                    PRINT*,"DBG: iSite, jSite, jState, stub(OR,APBS)", &
                         iSite, jSite, jState, stub
                    stub= SIGN(TINY,stub)
                 ENDIF
                 OnsiteRight=OnsitePotVec(iSite,jSite+2) &
                      /stub !(OnsitePotVec(iSite,jSite+1)*OnsitePotVec(iSite,jSite+2)-1.0D0)   ! antiperiodic BC
                 PsiRight=-Psi_A(jstate,mod(indexK,M)) &
                      /stub !(OnsitePotVec(iSite,jSite+1)*OnsitePotVec(iSite,jSite+2)-1.0D0)
              ENDIF
           ELSE
              stub= (OnsitePotVec(iSite,jSite+1)*OnsitePotVec(iSite,jSite+2)-1.0D0)
              IF( ABS(stub).LT.TINY) THEN
                 PRINT*,"DBG: iSite, jSite, jState, stub(OR)", &
                      iSite, jSite, jState, stub
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteRight=OnsitePotVec(iSite,jSite+2) &
                   /stub !(OnsitePotVec(iSite,jSite+1)*OnsitePotVec(iSite,jSite+2)-1.0D0)
              PsiRight=Psi_A(jstate,indexK+M) &
                   /stub !(OnsitePotVec(iSite,jSite+1)*OnsitePotVec(iSite,jSite+2)-1.0D0)
           END IF

           !PsiUp
           IF (mod(indexK,M).EQ.1) THEN

              IF (IBCFlag.EQ.0) THEN
                 OnsiteUp=ZERO      ! hard wall BC
                 PsiUp=ZERO
              ELSE IF (IBCFlag.EQ.1) THEN
                 stub= (OnsitePotVec(3*M,jSite)*OnsitePotVec(3*M-1,jSite)-1.0D0)
                 IF( ABS(stub).LT.TINY) THEN
                    PRINT*,"DBG: iSite, jSite, jState, stub(OU,HW)", &
                         iSite, jSite, jState, stub
                    stub= SIGN(TINY,stub)
                 ENDIF
                 OnsiteUp=OnsitePotVec(3*M-1,jSite) &
                      /stub !(OnsitePotVec(3*M,jSite)*OnsitePotVec(3*M-1,jSite)-1.0D0)   ! periodic BC
                 PsiUp=Psi_A(jstate,indexK+M-1) &
                      /stub !(OnsitePotVec(3*M,jSite)*OnsitePotVec(3*M-1,jSite)-1.0D0)

              ELSE IF (IBCFlag.EQ.2) THEN
                 stub= (OnsitePotVec(3*M,jSite)*OnsitePotVec(3*M-1,jSite)-1.0D0)
                 IF( ABS(stub).LT.TINY) THEN
                    PRINT*,"DBG: iSite, jSite, jState, stub(OU,APBC)", &
                         iSite, jSite, jState, stub
                    stub= SIGN(TINY,stub)
                 ENDIF
                 OnsiteUp=OnsitePotVec(3*M-1,jSite) &
                      /stub !(OnsitePotVec(3*M,jSite)*OnsitePotVec(3*M-1,jSite)-1.0D0)   ! antiperiodic BC
                 PsiUp=-Psi_A(jstate,indexK+M-1) &
                      /stub !(OnsitePotVec(3*M,jSite)*OnsitePotVec(3*M-1,jSite)-1.0D0)
              ENDIF
           ELSE
              stub= (OnsitePotVec(iSite-1,jSite)*OnsitePotVec(iSite-2,jSite)-1.0D0)
              IF( ABS(stub).LT.TINY) THEN
                 PRINT*,"DBG: iSite, jSite, jState, stub(OU)", &
                      iSite, jSite, jState, stub
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteUp=OnsitePotVec(iSite-2,jSite) &
                   /stub !(OnsitePotVec(iSite-1,jSite)*OnsitePotVec(iSite-2,jSite)-1.0D0)
              PsiUp=Psi_A(jstate,indexK-1) &
                   /stub !(OnsitePotVec(iSite-1,jSite)*OnsitePotVec(iSite-2,jSite)-1.0D0)
           END IF

           !PsiDown
           IF (mod(indexK,M).EQ.0) THEN

              IF (IBCFlag.EQ.0) THEN
                 stub= (OnsitePotVec(iSite+1,jSite)*OnsitePotVec(iSite+2,jSite)-1.0D0)
                 IF( ABS(stub).LT.TINY) THEN
                    PRINT*,"DBG: iSite, jSite, jState, stub(OD,HW)", &
                         iSite, jSite, jState, stub
                    stub= SIGN(TINY,stub)
                 ENDIF
                 OnsiteDown=OnsitePotVec(iSite+2,jSite) &
                      /stub !(OnsitePotVec(iSite+1,jSite)*OnsitePotVec(iSite+2,jSite)-1.0D0)   ! hard wall BC
                 PsiDown=ZERO                 
              ELSE IF (IBCFlag.EQ.1) THEN
                 stub= (OnsitePotVec(iSite+1,jSite)*OnsitePotVec(iSite+2,jSite)-1.0D0)
                 IF( ABS(stub).LT.TINY) THEN
                    PRINT*,"DBG: iSite, jSite, jState, stub(OD,PBC)", &
                         iSite, jSite, jState, stub
                    stub= SIGN(TINY,stub)
                 ENDIF
                 OnsiteDown=OnsitePotVec(iSite+2,jSite) &
                      /stub !(OnsitePotVec(iSite+1,jSite)*OnsitePotVec(iSite+2,jSite)-1.0D0)  ! periodic BC
                 PsiDown=Psi_A(jstate,((indexK/M-1)*M+1)) &
                      /stub !(OnsitePotVec(iSite+1,jSite)*OnsitePotVec(iSite+2,jSite)-1.0D0)
              ELSE IF (IBCFlag.EQ.2) THEN
                 stub= (OnsitePotVec(iSite+1,jSite)*OnsitePotVec(iSite+2,jSite)-1.0D0)
                 IF( ABS(stub).LT.TINY) THEN
                    PRINT*,"DBG: iSite, jSite, jState, stub(OU,APBC)", &
                         iSite, jSite, jState, stub
                    stub= SIGN(TINY,stub)
                 ENDIF
                 OnsiteDown=OnsitePotVec(iSite+2,jSite) &
                      /stub !(OnsitePotVec(iSite+1,jSite)*OnsitePotVec(iSite+2,jSite)-1.0D0)   ! antiperiodic BC
                 PsiDown=-Psi_A(jstate,((indexK/M-1)*M+1)) &
                      /stub !(OnsitePotVec(iSite+1,jSite)*OnsitePotVec(iSite+2,jSite)-1.0D0)
              ENDIF
           ELSE
              stub= (OnsitePotVec(iSite+1,jSite)*OnsitePotVec(iSite+2,jSite)-1.0D0)
              IF( ABS(stub).LT.TINY) THEN
                 PRINT*,"DBG: iSite, jSite, jState, stub(OU)", &
                      iSite, jSite, jState, stub
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteDown=OnsitePotVec(iSite+2,jSite) &
                   /stub !(OnsitePotVec(iSite+1,jSite)*OnsitePotVec(iSite+2,jSite)-1.0D0)
              PsiDown=Psi_A(jstate,indexK+1) &
                   /stub !(OnsitePotVec(iSite+1,jSite)*OnsitePotVec(iSite+2,jSite)-1.0D0)
           END IF

           !PRINT*,"DBG2: jState,iSite, jSite, indexK", jState, iSite, jSite, indexK
           new= ( OnsitePot + OnsiteLeft + OnsiteRight + OnsiteUp +  OnsiteDown ) * Psi_A(jstate,indexK)&
                - Kappa * ( PsiLeft + PsiRight + PsiUp + PsiDown  ) &
                - PSI_B(jstate,indexK) 

           !PRINT*,"iSite,jSite,En, OP, PL, PR, PA,PB, PN"
           !PRINT*, iSite, jState, En, OnsitePot, PsiLeft, PsiRight,
           !        PSI_A(iSite,jState), PSI_B(iSite,jState),
           !        new

           PSI_B(jstate,indexK)= new

        ENDDO !jState

     ENDDO ! iSite
  ENDDO!jSite

  !PRINT*,"PSIA(1,1),(2,1),(1,2),(2,2)",&
  !      PSI_A(1,1),PSI_A(2,1),PSI_A(1,2),PSI_A(2,2)

  RETURN

END SUBROUTINE TMMultLieb3DAtoB5

  
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
        
        new= ( OnsitePot * PSI_A(jState,iSite) &
             - PSI_B(jState,iSite) )
        
        !PRINT*,"i,jSite,En, OP, PL, PR, PA,PB, PN"
        !PRINT*, iSite, jState, En, OnsitePot, PsiLeft, PsiRight,
        !        PSI_A(iSite,jState), PSI_B(iSite,jState),
        !        new
        
        PSI_B(jState,iSite)= new
        
     ENDDO ! jState
  ENDDO ! iSite
  
  !PRINT*,"PSIA(1,1),(1,2),(1,3),(1,4)",&
        !PSI_A(1,1),PSI_A(1,2),PSI_A(1,3),PSI_A(1,4)

  !PRINT*,"PSIB(1,1),(1,2),(1,3),(1,4)",&
        !PSI_B(1,1),PSI_B(1,2),PSI_B(1,3),PSI_B(1,4)
  
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
  
  INTEGER iSite, jState, ISeedDummy
  REAL(KIND=RKIND) OnsitePot
  REAL(KIND=CKIND) new

  CALL TMMultLieb3DB5toB6(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

  RETURN
END SUBROUTINE TMMultLieb3DB6toA
