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

  REAL(KIND=RKIND) PSI_A(M*M,M*M),PSI_B(M*M,M*M),OnsitePotVec(3*M*3*M)

  INTEGER jState, ISeedDummy,iSiteS,iSiteL
  REAL(KIND=RKIND) OnsitePot, OnsiteRight, OnsiteLeft, OnsiteUp, OnsiteDown
  REAL(KIND=RKIND) NEW, PsiLeft, PsiRight, PsiUp, PsiDown, stub

  !PRINT*,"DBG: TMMultLieb3DAtoB()"

  ! create the new onsite potential
  DO iSiteS=1,3*M*3*M
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePotVec(iSiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
     CASE(1)
        OnsitePotVec(iSiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
     CASE(2)
        OnsitePotVec(iSiteS)= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
  END DO

  !PRINT*,"iS,pL,RndVec", iSite,pLevel,RndVec((pLevel-1)*M+iSite)

  !new TMM
  DO iSiteL=1,M*M
     IF (MOD(iSiteL,M)==0) THEN
        iSiteS=(iSiteL/M-1)*6*M+3*iSiteL-2
     ELSE
        iSiteS=(iSiteL/M)*9*M+MOD(iSiteL,M)*3-2
     ENDIF
     
     OnsitePot=OnsitePotVec(iSiteS)

     DO jState=1,M*M
        !PsiLeft
        IF (iSiteL.LE.M) THEN
           IF (IBCFlag.EQ.0) THEN       ! hard wall BC
              OnsiteLeft= ZERO
              PsiLeft= ZERO               
           ELSE IF (IBCFlag.EQ.1) THEN  ! periodic BC
              stub= OnsitePotVec(9*M*M-3*M+iSiteS)*OnsitePotVec(9*M*M-6*M+iSiteS)-1.0D0
              IF( ABS(stub).LT.TINY) THEN
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteLeft= OnsitePotVec(9*M*M-6*M+iSiteS)/stub
              PsiLeft= Psi_A(jState,M*M-M+MOD(iSiteL,M))/stub
           ELSE IF (IBCFlag.EQ.2) THEN  ! antiperiodic BC
              stub= OnsitePotVec(9*M*M-3*M+iSiteS)*OnsitePotVec(9*M*M-6*M+iSiteS)-1.0D0
              IF( ABS(stub).LT.TINY) THEN            
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteLeft= OnsitePotVec(9*M*M-6*M+iSiteS)/stub
              PsiLeft= -Psi_A(jState,M*M-M+MOD(iSiteL,M))/stub
           ENDIF
        ELSE
           stub= OnsitePotVec(iSiteS-3*M)*OnSitePotVec(iSiteS-6*M)-1.0D0
           IF( ABS(stub).LT.TINY) THEN           
              stub= SIGN(TINY,stub)
           ENDIF
           OnsiteLeft= OnsitePotVec(iSiteS-6*M)/stub
           PsiLeft= Psi_A(jState,iSiteL-M)/stub
        END IF
 
        !PsiRight
        IF (iSiteL.GT.M*(M-1)) THEN
           IF (IBCFlag.EQ.0) THEN        ! hard wall BC
              stub= (OnsitePotVec(iSiteS+3*M)*OnSitePotVec(iSiteS+6*M)-1.0D0) 
              IF( ABS(stub).LT.TINY) THEN
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteRight= OnsitePotVec(iSiteS+6*M)/stub 
              PsiRight= ZERO        
           ELSE IF (IBCFlag.EQ.1) THEN   ! periodic BC
              stub= (OnsitePotVec(iSiteS+3*M)*OnSitePotVec(iSiteS+6*M)-1.0D0) 
              IF( ABS(stub).LT.TINY) THEN
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteRight= OnsitePotVec(iSiteS+6*M)/stub 
              PsiRight= Psi_A(jState,MOD(iSiteL,M))/stub
           ELSE IF (IBCFlag.EQ.2) THEN   ! antiperiodic BC             
              stub= (OnsitePotVec(iSiteS+3*M)*OnSitePotVec(iSiteS+6*M)-1.0D0) 
              IF( ABS(stub).LT.TINY) THEN
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteRight= OnsitePotVec(iSiteS+6*M)/stub 
              PsiRight= -Psi_A(jState,MOD(iSiteL,M))/stub
           ENDIF
        ELSE
           stub= (OnsitePotVec(iSiteS+3*M)*OnSitePotVec(iSiteS+6*M)-1.0D0)
           IF( ABS(stub).LT.TINY) THEN             
              stub= SIGN(TINY,stub)
           ENDIF
           OnsiteRight= OnsitePotVec(iSiteS+6*M)/stub
           PsiRight= Psi_A(jState,iSiteL+M)/stub
        END IF
        
        !PsiUp
        IF (MOD(iSiteL,M).EQ.1) THEN
           IF (IBCFlag.EQ.0) THEN        ! hard wall BC
              OnsiteUp=ZERO      
              PsiUp=ZERO
           ELSE IF (IBCFlag.EQ.1) THEN   ! periodic BC
              stub= (OnsitePotVec(iSiteS+3*M-1)*OnSitePotVec(iSiteS+3*M-2)-1.0D0)
              IF( ABS(stub).LT.TINY) THEN
                 stub= SIGN(TINY,stub)
              ENDIF 
              OnsiteUp= OnsitePotVec(iSiteS+3*M-2)/stub
              PsiUp= Psi_A(jState,iSiteL+M-1)/stub
           ELSE IF (IBCFlag.EQ.2) THEN   ! antiperiodic BC
              stub= (OnsitePotVec(iSiteS+3*M-1)*OnSitePotVec(iSiteS+3*M-2)-1.0D0)
              IF( ABS(stub).LT.TINY) THEN
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteUp= OnsitePotVec(iSiteS+3*M-2)/stub
              PsiUp= -Psi_A(jState,iSiteL+M-1)/stub
           ENDIF
        ELSE
           stub= (OnsitePotVec(iSiteS-1)*OnSitePotVec(iSiteS-2)-1.0D0)
           IF( ABS(stub).LT.TINY) THEN
              stub= SIGN(TINY,stub)
           ENDIF
           OnsiteUp= OnsitePotVec(iSiteS-2)/stub
           PsiUp= Psi_A(jState,iSiteL-1)/stub
        END IF
           
        !PsiDown
        IF (MOD(iSiteL,M).EQ.0) THEN
           IF (IBCFlag.EQ.0) THEN       ! hard wall BC
              stub= (OnsitePotVec(iSiteS+1)*OnSitePotVec(iSiteS+2)-1.0D0)
              IF( ABS(stub).LT.TINY) THEN
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteDown= OnsitePotVec(iSiteS+2)/stub
              PsiDown= ZERO                 
           ELSE IF (IBCFlag.EQ.1) THEN   ! periodic BC
              stub= (OnsitePotVec(iSiteS+1)*OnSitePotVec(iSiteS+2)-1.0D0)
              IF( ABS(stub).LT.TINY) THEN
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteDown= OnsitePotVec(iSiteS+2)/stub
              PsiDown= Psi_A(jState,iSiteL-M+1)/stub
           ELSE IF (IBCFlag.EQ.2) THEN   ! antiperiodic BC
              stub= (OnsitePotVec(iSiteS+1)*OnSitePotVec(iSiteS+2)-1.0D0)
              IF( ABS(stub).LT.TINY) THEN
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteDown= OnsitePotVec(iSiteS+2)/stub
              PsiDown= -Psi_A(jState,iSiteL-M+1)/stub
           ENDIF
        ELSE
           stub= (OnsitePotVec(iSiteS+1)*OnSitePotVec(iSiteS+2)-1.0D0)
           IF( ABS(stub).LT.TINY) THEN
              stub= SIGN(TINY,stub)
           ENDIF
           OnsiteDown= OnsitePotVec(iSiteS+2)/(OnsitePotVec(iSiteS+1)*OnSitePotVec(iSiteS+2)-1.0D0)
           PsiDown=  Psi_A(jState,iSiteL+1)/(OnsitePotVec(iSiteS+1)*OnSitePotVec(iSiteS+2)-1.0D0)
        END IF
                
        NEW= ( OnsitePot - OnsiteLeft - OnsiteRight - OnsiteUp - OnsiteDown ) * Psi_A(jState,iSiteL)&
                - Kappa * ( PsiLeft + PsiRight + PsiUp + PsiDown  ) &
                - PSI_B(jState,iSiteL) 
        PSI_B(jState,iSiteL)= NEW
     ENDDO !jState
  ENDDO !iSiteL
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
        
        NEW= ( OnsitePot * PSI_A(jState,iSite) &
             - PSI_B(jState,iSite) )
        
        !PRINT*,"i,jSite,En, OP, PL, PR, PA,PB, PN"
        !PRINT*, iSite, jState, En, OnsitePot, PsiLeft, PsiRight,
        !        PSI_A(iSite,jState), PSI_B(iSite,jState),
        !        new
        
        PSI_B(jState,iSite)= NEW
        
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
  REAL(KIND=CKIND) NEW

  CALL TMMultLieb3DB5toB6(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

  RETURN
END SUBROUTINE TMMultLieb3DB6toA
