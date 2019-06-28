! --------------------------------------------------------------------
! Multiplication Of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_B,PSI_A) so that the structure of the transfer matrix 
! can be exploited

SUBROUTINE TMMultLieb2D_AtoB1(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

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
  
  REAL(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  
  INTEGER iSite, jSite,jState, indexK,ISeedDummy
  REAL(KIND=RKIND) OnsitePot, OnsiteDown, OnsiteUp, OnsitePotVec(4*M,1)
  REAL(KIND=CKIND) new , PsiUp, PsiDown,stub

  ! create the new onsite potential
   jSite=1

  DO iSite=1,4*M   !label different row
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePotVec(iSite,jSite)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
     CASE(1)
        OnsitePotVec(iSite,jSite)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
     CASE(2)
        OnsitePotVec(iSite,jSite)= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
  END DO

  ! to the TMM
  DO iSite=1,4*M,4  ! lable the the row of A atom.
     
     indexK=(iSite/4+1)      !label the A atom
     OnsitePot=OnsitePotVec(iSite,jSite)
     
     DO jState=1,M
        ! PSI UP
        IF (indexK==1) THEN
           
           IF (IBCFlag.EQ.0) THEN ! hard wall BC
              OnsiteUp=ZERO      
              PsiUp=ZERO              
           ELSE IF (IBCFlag.EQ.1) THEN ! periodic BC
              CONTINUE              
           ELSE IF (IBCFlag.EQ.2) THEN ! antiperiodic BC
              CONTINUE                      
           ENDIF
        ELSE
           stub= (OnsitePotVec(iSite-1,jSite)*OnsitePotVec(iSite-2,jSite)*OnsitePotVec(iSite-3,jSite) &
           -OnsitePotVec(iSite-1,jSite)-OnsitePotVec(iSite-3,jSite))
           IF( ABS(stub).LT.TINY) THEN
!!$              PRINT*,"DBG: iSite, jSite, jState, stub(OU)", iSite, jSite, jState, stub
              stub= SIGN(TINY,stub)
           ENDIF
           OnsiteUp=(OnsitePotVec(iSite-2,jSite)*OnsitePotVec(iSite-3,jSite)-1.0D0)/stub
           PsiUp=PSI_A(indexK-1,jState)/stub 
        END IF
        
        !PSI DOWN
        IF (indexK==M) THEN
           
           IF (IBCFlag.EQ.0) THEN ! hard wall BC
              stub= (OnsitePotVec(iSite+1,jSite)*OnsitePotVec(iSite+2,jSite)*OnsitePotVec(iSite+3,jSite) &
                   -OnsitePotVec(iSite+1,jSite)-OnsitePotVec(iSite+3,jSite))
              IF( ABS(stub).LT.TINY) THEN
!!$                 PRINT*,"DBG: iSite, jSite, jState, stub(OU)", iSite, jSite, jState, stub
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteDown=(OnsitePotVec(iSite+2,jSite)*OnsitePotVec(iSite+3,jSite)-1.0D0)/stub 
              PsiDown=ZERO                     
           ELSE IF (IBCFlag.EQ.1) THEN ! periodic BC
              CONTINUE
           ELSE IF (IBCFlag.EQ.2) THEN ! antiperiodic BC
              CONTINUE          
           ENDIF
        ELSE
           stub= (OnsitePotVec(iSite+1,jSite)*OnsitePotVec(iSite+2,jSite)*OnsitePotVec(iSite+3,jSite) &
                -OnsitePotVec(iSite+1,jSite)-OnsitePotVec(iSite+3,jSite))
           IF( ABS(stub).LT.TINY) THEN
!!$              PRINT*,"DBG: iSite, jSite, jState, stub(OU)", iSite, jSite, jState, stub
              stub= SIGN(TINY,stub)
           ENDIF
           OnsiteDown=(OnsitePotVec(iSite+2,jSite)*OnsitePotVec(iSite+3,jSite)-1.0D0)/stub  
           PsiDown=PSI_A(indexK+1,jState) /stub 
        END IF
        
        new =(( OnsitePot-OnsiteUp-OnsiteDown )*PSI_A(indexK,jState) &
             - Kappa * ( PsiUp + PsiDown ) &
             - PSI_B(indexK,jState) )
         PSI_B(indexK,jState)= new
        
      ENDDO ! jState
   ENDDO ! iSite
  
   RETURN
    
 END SUBROUTINE TMMultLieb2D_AtoB1

 SUBROUTINE TMMultLieb2D_B1toB2(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

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
  
  REAL(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  
  INTEGER iSite, jState, ISeedDummy
  REAL(KIND=RKIND) OnsitePot
  REAL(KIND=CKIND) new, PsiLeft, PsiRight
  
  DO iSite=1,M
     
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
     
     DO jState=1,M
        
        !PRINT*,"jState, iSite", jState, iSite,
        
        new= ( OnsitePot * PSI_A(iSite,jState) &
             - PSI_B(iSite,jState) )
        
        PSI_B(iSite,jState)= new
        
     ENDDO ! jState
  ENDDO ! iSite
  
  RETURN
END SUBROUTINE TMMultLieb2D_B1toB2


SUBROUTINE TMMultLieb2D_B2toB3(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

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
  
  REAL(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  
  INTEGER iSite, jState, ISeedDummy
  REAL(KIND=RKIND) OnsitePot
  REAL(KIND=CKIND) new, PsiLeft, PsiRight
  
  !PRINT*,"DBG: TMMultLieb2DB1toB2()"

  CALL TMMultLieb2D_B1toB2(PSI_A,PSI_B, Ilayer, En, DiagDis, M )
  RETURN
END SUBROUTINE TMMultLieb2D_B2toB3


SUBROUTINE TMMultLieb2D_B3toA(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

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
  
  REAL(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  
  INTEGER iSite, jState, ISeedDummy
  REAL(KIND=RKIND) OnsitePot
  REAL(KIND=CKIND) new, PsiLeft, PsiRight
  
  !PRINT*,"DBG: TMMultLieb2DB1toB2()"

  CALL TMMultLieb2D_B1toB2(PSI_A,PSI_B, Ilayer, En, DiagDis, M )
  RETURN
END SUBROUTINE TMMultLieb2D_B3toA
