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
  
  INTEGER iSite,jSite,indexK,jState, ISeedDummy
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
       ! PRINT*,"DBG1: iSite, jSite, indexK", iSite, jSite, indexK
        
        DO jState=1,M*M
           
           !PsiLeft
           IF (indexK<=M) THEN
              SELECT CASE(IBCFlag)
              CASE(-1,0) ! hard wall BC
                 PsiLeft= ZERO            
                 OnsiteLeft= ZERO
              CASE(1) ! periodic BC
                 OnsiteLeft=1.0D0/OnsitePotVec(iSite,2*M)    
                 PsiLeft=1.0D0/OnsitePotVec(iSite,2*M)*PSI_A(indexK+M**2-M,jState)
              CASE(2) ! antiperiodic BC
                 OnsiteLeft=1.0D0/OnsitePotVec(iSite,2*M)  
                 PsiLeft=-1.0D0/OnsitePotVec(iSite,2*M)*PSI_A(indexK+M**2-M,jState)    
              CASE DEFAULT
                 PRINT*,"TMMultLieb2DAtoB(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              OnsiteLeft=1.0D0/OnsitePotVec(iSite,jSite-1)
              PsiLeft=1.0D0/OnsitePotVec(iSite,jSite-1)*PSI_A(indexK-M,jState)
           END IF
           
           !PsiRight
           IF (indexK>M*(M-1)) THEN
              SELECT CASE(IBCFlag)
              CASE(-1) ! hard wall BC with stubs
                 OnsiteRight= 1.0D0/OnsitePotVec(iSite,jSite+1)    
                 PsiRight= ZERO
              CASE(0) ! hard wall BC
                 OnsiteRight= ZERO
                 PsiRight= ZERO
              CASE(1) ! periodic BC
                 OnsiteRight= 1.0D0/OnsitePotVec(iSite,jSite+1)    
                 !   PsiRight= 1.0D0/OnsitePotVec(iSite,jSite+1)*PSI_A(jState,mod(indexK,M))
                 PsiRight= 1.0D0/OnsitePotVec(iSite,jSite+1)*PSI_A(indexK-M*(M-1),jState)
              CASE(2) ! antiperiodic BC
                 OnsiteRight= 1.0D0/OnsitePotVec(iSite,jSite+1)   
                 PsiRight= -1.0D0/OnsitePotVec(iSite,jSite+1)*PSI_A(indexK-M*(M-1),jState)  
              CASE DEFAULT
                 PRINT*,"TMMultLieb2DAtoB(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              OnsiteRight=1.0D0/OnsitePotVec(iSite,jSite+1)
              PsiRight=1.0D0/OnsitePotVec(iSite,jSite+1)*PSI_A(indexK+M,jState)
           END IF
           
           !PsiUp
           IF (mod(indexK,M).EQ.1) THEN
              SELECT CASE(IBCFlag)
              CASE(-1,0) ! hard wall BC
                 OnsiteUp=ZERO      
                 PsiUp=ZERO
              CASE(1) ! periodic BC
                 OnsiteUp=1.0D0/OnsitePotVec(2*M,jSite)    
                 PsiUp=1.0D0/OnsitePotVec(2*M,jSite)*PSI_A(indexK+M-1,jState)
              CASE(2) ! antiperiodic BC
                 OnsiteUp=1.0D0/OnsitePotVec(2*M,jSite)    
                 PsiUp=-1.0D0/OnsitePotVec(2*M,jSite)*PSI_A(indexK+M-1,jState)
              CASE DEFAULT
                 PRINT*,"TMMultLieb2DAtoB(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              OnsiteUp=1.0D0/OnsitePotVec(iSite-1,jSite)
              PsiUp=1.0D0/OnsitePotVec(iSite-1,jSite)*PSI_A(indexK-1,jState)
           END IF
           
           !PsiDown
           IF (mod(indexK,M).EQ.0) THEN
              SELECT CASE(IBCFLag)
              CASE(-1) ! hard wall BC with stubs
                 OnsiteDown=1.0D0/OnsitePotVec(iSite+1,jSite)  
                 PsiDown=ZERO
              CASE(0) ! hard wall BC
                 OnsiteDown=ZERO
                 PsiDown=ZERO
              CASE(1) ! periodic BC
                 OnsiteDown=1.0D0/OnsitePotVec(iSite+1,jSite)  
                 PsiDown=1.0D0/OnsitePotVec(iSite+1,jSite)*PSI_A(indexK-M+1,jState)
              CASE(2) ! antiperiodic BC
                 OnsiteDown=1.0D0/OnsitePotVec(iSite+1,jSite)  
                 PsiDown=-1.0D0/OnsitePotVec(iSite+1,jSite)*PSI_A(indexK-M+1,jState)                 
              CASE DEFAULT
                 PRINT*,"TMMultLieb2DAtoB(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              OnsiteDown=1.0D0/OnsitePotVec(iSite+1,jSite)
              PsiDown=1.0D0/OnsitePotVec(iSite+1,jSite)*PSI_A(indexK+1,jState)
           END IF

           !PRINT*,"DBG2: jState,iSite, jSite, indexK", jState, iSite, jSite, indexK
           new= ( OnsitePot - OnsiteLeft - OnsiteRight - OnsiteUp - OnsiteDown ) * PSI_A(indexK,jState)&
                - Kappa * ( PsiLeft + PsiRight + PsiUp + PsiDown  ) &
                - PSI_B(indexK,jState) 
           
           !PRINT*,"iSite,jSite,En, OP, PL, PR, PA,PB, PN"
           !PRINT*, iSite, jState, En, OnsitePot, PsiLeft, PsiRight,
           !        PSI_A(iSite,jState), PSI_B(iSite,jState),
           !        new
           
           PSI_B(indexK,jState)= new
           
        ENDDO !jState
        
     ENDDO ! iSite
  ENDDO!jSite
  
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

