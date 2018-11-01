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
  
  INTEGER iSite,jSite,indexK,jstate, ISeedDummy
  REAL(KIND=RKIND) OnsitePot, OnsiteRight, OnsiteLeft, OnsiteUp, OnsiteDown
  REAL(KIND=RKIND) new, PsiLeft, PsiRight, PsiUp, PsiDown

  !PRINT*,"DBG: TMMultLieb3DAtoB()"
     
  ! create the new onsite potential
  DO iSite=1,2*M
     DO jSite=1,2*M
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
  DO jSite=1,2*M,2
     DO iSite=1,2*M,2
        
        indexK=(jSite/2)*M+(iSite+1)/2
        OnsitePot=OnsitePotVec(iSite,jSite)
        PRINT*,"DBG1: iSite, jSite, indexK", iSite, jSite, indexK
        
        DO jstate=1,M*M
           
           !PsiLeft
           IF (indexK<=M) THEN
              
              IF (IBCFlag.EQ.0) THEN
                 PsiLeft= ZERO            ! hard wall BC
                 OnsiteLeft= ZERO
              ELSE IF (IBCFlag.EQ.1) THEN
                 OnsiteLeft=1.0D0/OnsitePotVec(iSite,2*M)    ! periodic BC
                 PsiLeft=1.0D0/OnsitePotVec(iSite,2*M)*Psi_A(jstate,indexK+M**2-M) 
              ELSE IF (IBCFlag.EQ.2) THEN
                 OnsiteLeft=1.0D0/OnsitePotVec(iSite,2*M)   ! antiperiodic BC
                 PsiLeft=-1.0D0/OnsitePotVec(iSite,2*M)*Psi_A(jstate,indexK+M**2-M)    
              ENDIF
           ELSE
              OnsiteLeft=1.0D0/OnsitePotVec(iSite,jSite-1)
              PsiLeft=1.0D0/OnsitePotVec(iSite,jSite-1)*Psi_A(jstate,indexK-M)
           END IF
           
           !PsiRight
           IF (indexK>M*(M-1)) THEN
              
              IF (IBCFlag.EQ.0) THEN
                 OnsiteRight= 1.0D0/OnsitePotVec(iSite,jSite+1)    ! hard wall BC
                 PsiRight= ZERO        
              ELSE IF (IBCFlag.EQ.1) THEN
                 OnsiteRight= 1.0D0/OnsitePotVec(iSite,jSite+1)    ! periodic BC
                 PsiRight= 1.0D0/OnsitePotVec(iSite,jSite+1)*Psi_A(jstate,mod(indexK,M))
              ELSE IF (IBCFlag.EQ.2) THEN
                 OnsiteRight= 1.0D0/OnsitePotVec(iSite,jSite+1)    ! antiperiodic BC
                 PsiRight= -1.0D0/OnsitePotVec(iSite,jSite+1)*Psi_A(jstate,mod(indexK,M))  
              ENDIF
           ELSE
              OnsiteRight=1.0D0/OnsitePotVec(iSite,jSite+1)
              PsiRight=1.0D0/OnsitePotVec(iSite,jSite+1)*Psi_A(jstate,indexK+M)
           END IF
           
           !PsiUp
           IF (mod(indexK,M).EQ.1) THEN
              
              IF (IBCFlag.EQ.0) THEN
                 OnsiteUp=ZERO      ! hard wall BC
                 PsiUp=ZERO
              ELSE IF (IBCFlag.EQ.1) THEN
                 OnsiteUp=1.0D0/OnsitePotVec(2*M,jSite)     ! periodic BC
                 PsiUp=1.0D0/OnsitePotVec(2*M,jSite)*Psi_A(jstate,indexK+M-1)
              ELSE IF (IBCFlag.EQ.2) THEN
                 OnsiteUp=1.0D0/OnsitePotVec(2*M,jSite)     ! antiperiodic BC
                 PsiUp=-1.0D0/OnsitePotVec(2*M,jSite)*Psi_A(jstate,indexK+M-1)
                 
              ENDIF
           ELSE
              OnsiteUp=1.0D0/OnsitePotVec(iSite-1,jSite)
              PsiUp=1.0D0/OnsitePotVec(iSite-1,jSite)*Psi_A(jstate,indexK-1)
           END IF
           
           !PsiDown
           IF (mod(indexK,M).EQ.0) THEN
              
              IF (IBCFlag.EQ.0) THEN
                 OnsiteDown=1.0D0/OnsitePotVec(iSite+1,jSite)  ! hard wall BC
                 PsiDown=ZERO                 
              ELSE IF (IBCFlag.EQ.1) THEN
                 OnsiteDown=1.0D0/OnsitePotVec(iSite+1,jSite)   ! periodic BC
                 PsiDown=1.0D0/OnsitePotVec(iSite+1,jSite)*Psi_A(jstate,indexK-M+1)
              ELSE IF (IBCFlag.EQ.2) THEN
                 OnsiteDown=1.0D0/OnsitePotVec(iSite+1,jSite)   ! antiperiodic BC
                 PsiDown=-1.0D0/OnsitePotVec(iSite+1,jSite)*Psi_A(jstate,indexK-M+1)                 
              ENDIF
           ELSE
              OnsiteDown=1.0D0/OnsitePotVec(iSite+1,jSite)
              PsiDown=1.0D0/OnsitePotVec(iSite+1,jSite)*Psi_A(jstate,indexK+1)
           END IF

           !PRINT*,"DBG2: jState,iSite, jSite, indexK", jState, iSite, jSite, indexK
           new= ( OnsitePot - OnsiteLeft - OnsiteRight - OnsiteUp - OnsiteDown ) * Psi_A(jstate,indexK)&
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

END SUBROUTINE TMMultLieb3DAtoB

  
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
END SUBROUTINE TMMultLieb3DBtoA

