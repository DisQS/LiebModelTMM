! --------------------------------------------------------------------
! Multiplication Of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_B,PSI_A) so that the structure of the transfer matrix 
! can be exploited

SUBROUTINE TMMultLieb2DAtoB1(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

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
  REAL(KIND=RKIND) OnsitePot, OnsiteDown, OnsiteUp, OnsitePotVec(3*M,1)
  REAL(KIND=CKIND) new , PsiUp, PsiDown,stub
  
  !PRINT*,"DBG: TMMultLieb2DAtoB1()"

  ! create the new onsite potential
   
  ! DO jSite=1  !label different column
  jSite=1

  DO iSite=1,3*M   !label different row
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePotVec(iSite,jSite)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
     CASE(1)
        OnsitePotVec(iSite,jSite)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
     CASE(2)
        OnsitePotVec(iSite,jSite)= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
     
!!$        IF( ABS(OnsitePotVec(iSite,jSite)).LT.TINY) THEN
!!$           OnsitePotVec(iSite,jSite)= SIGN(TINY,OnsitePotVec(iSite,jSite))
!!$        END IF
  END DO
!  END DO

  
!!$  DO iSiteS=1,3*M   
!!$     SELECT CASE(IRNGFlag)
!!$     CASE(0)
!!$        OnsitePotVec(iSiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
!!$     CASE(1)
!!$        OnsitePotVec(iSiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
!!$     CASE(2)
!!$        OnsitePotVec(iSiteS)= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
!!$     END SELECT
!!$
!!$     IF( ABS(OnsitePotVec(iSiteS)).LT.TINY) THEN
!!$        OnsitePotVec(iSiteS)= SIGN(TINY,OnsitePotVec(iSiteS))
!!$     END IF
!!$  END DO

  
  ! to the TMM


  ! DO jSite=1   !label the first colum
  !  jSite=1
  DO iSite=1,3*M,3  ! lable the the row of A atom.
     
     indexK=(iSite/3+1)      !label the A atom
     OnsitePot=OnsitePotVec(iSite,jSite)
     
     DO jState=1,M
        ! PSI UP
        IF (indexK==1) THEN
           
           IF (IBCFlag.EQ.0) THEN ! hard wall BC
              OnsiteUp=ZERO      
              PsiUp=ZERO
              
           ELSE IF (IBCFlag.EQ.1) THEN ! periodic BC
              stub= (OnsitePotVec(3*M,jSite)*OnsitePotVec(3*M-1,jSite)-1.0D0)
              IF( ABS(stub).LT.TINY) THEN
                 PRINT*,"DBG: iSite, jSite, jState, stub(OU)", iSite, jSite, jState, stub
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteUp=OnsitePotVec(3*M-1,jSite) &
                   /stub !(OnsitePotVec(iSite-1,jSite)*OnsitePotVec(iSite-2,jSite)-1.0D0)
              PsiUp=PSI_A(M,jState) &
                   /stub !(OnsitePotVec(iSite-1,jSite)*OnsitePotVec(iSite-2,jSite)-1.0D0)
              
           ELSE IF (IBCFlag.EQ.2) THEN ! antiperiodic BC
              stub= (OnsitePotVec(3*M,jSite)*OnsitePotVec(3*M-1,jSite)-1.0D0)
              IF( ABS(stub).LT.TINY) THEN
                 PRINT*,"DBG: iSite, jSite, jState, stub(OU)", iSite, jSite, jState, stub
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteUp=OnsitePotVec(3*M-1,jSite) &
                   /stub !(OnsitePotVec(iSite-1,jSite)*OnsitePotVec(iSite-2,jSite)-1.0D0)
              PsiUp=-PSI_A(M,jState) &
                   /stub !(OnsitePotVec(iSite-1,jSite)*OnsitePotVec(iSite-2,jSite)-1.0D0)                       
           ENDIF
        ELSE
           stub= (OnsitePotVec(iSite-1,jSite)*OnsitePotVec(iSite-2,jSite)-1.0D0)
           IF( ABS(stub).LT.TINY) THEN
              PRINT*,"DBG: iSite, jSite, jState, stub(OU)", iSite, jSite, jState, stub
              stub= SIGN(TINY,stub)
           ENDIF
           OnsiteUp=OnsitePotVec(iSite-2,jSite) &
                /stub !(OnsitePotVec(iSite-1,jSite)*OnsitePotVec(iSite-2,jSite)-1.0D0)
           PsiUp=PSI_A(indexK-1,jState) &
                /stub !(OnsitePotVec(iSite-1,jSite)*OnsitePotVec(iSite-2,jSite)-1.0D0)
        END IF
        
        !PSI DOWN
        IF (indexK==M) THEN
           
           IF (IBCFlag.EQ.0) THEN ! hard wall BC
              stub= (OnsitePotVec(iSite+1,jSite)*OnsitePotVec(iSite+2,jSite)-1.0D0)
              IF( ABS(stub).LT.TINY) THEN
                 PRINT*,"DBG: iSite, jSite, jState, stub(OU)", iSite, jSite, jState, stub
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteDown=OnsitePotVec(iSite+2,jSite) &
                   /stub !(OnsitePotVec(iSite-1,jSite)*OnsitePotVec(iSite-2,jSite)-1.0D0)
              PsiDown=ZERO        
              
           ELSE IF (IBCFlag.EQ.1) THEN ! periodic BC
              stub= (OnsitePotVec(iSite+1,jSite)*OnsitePotVec(iSite+2,jSite)-1.0D0)
              IF( ABS(stub).LT.TINY) THEN
                 PRINT*,"DBG: iSite, jSite, jState, stub(OU)", iSite, jSite, jState, stub
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteDown=OnsitePotVec(iSite+2,jSite) &
                   /stub !(OnsitePotVec(iSite-1,jSite)*OnsitePotVec(iSite-2,jSite)-1.0D0)
              PsiDown=PSI_A(1,jState) &
                   /stub !(OnsitePotVec(iSite-1,jSite)*OnsitePotVec(iSite-2,jSite)-1.0D0)
           ELSE IF (IBCFlag.EQ.2) THEN ! antiperiodic BC
              stub= (OnsitePotVec(iSite+1,jSite)*OnsitePotVec(iSite+2,jSite)-1.0D0)
              IF( ABS(stub).LT.TINY) THEN
                 PRINT*,"DBG: iSite, jSite, jState, stub(OU)", iSite, jSite, jState, stub
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteDown=OnsitePotVec(iSite+2,jSite) &
                   /stub !(OnsitePotVec(iSite-1,jSite)*OnsitePotVec(iSite-2,jSite)-1.0D0)
              PsiDown=-PSI_A(1,jState) &
                   /stub !(OnsitePotVec(iSite-1,jSite)*OnsitePotVec(iSite-2,jSite)-1.0D0)                  
           ENDIF
        ELSE
           stub= (OnsitePotVec(iSite+1,jSite)*OnsitePotVec(iSite+2,jSite)-1.0D0)
           IF( ABS(stub).LT.TINY) THEN
              PRINT*,"DBG: iSite, jSite, jState, stub(OU)", iSite, jSite, jState, stub
              stub= SIGN(TINY,stub)
           ENDIF
           OnsiteDown=OnsitePotVec(iSite+2,jSite) &
                /stub !(OnsitePotVec(iSite-1,jSite)*OnsitePotVec(iSite-2,jSite)-1.0D0)
           PsiDown=PSI_A(indexK+1,jState) &
                /stub !(OnsitePotVec(iSite-1,jSite)*OnsitePotVec(iSite-2,jSite)-1.0D0)
        END IF
        
        new =(( OnsitePot-OnsiteUp-OnsiteDown )*PSI_A(indexK,jState) &
             - Kappa * ( PsiUp + PsiDown ) &
             - PSI_B(indexK,jState) )
        
!!$        PRINT*,"i,j,En, OP, PL, PR, PA,PB, PN"
!!$        PRINT*, iSiteL, jState, En, OnsitePot, PsiLeft, PsiRight, &
!!$                PSI_A(iSiteL,jState), PSI_B(iSiteL,jState), &
!!$                new
        
        PSI_B(indexK,jState)= new
        
     ENDDO ! jState
  ENDDO ! iSite
  
  RETURN
    
END SUBROUTINE TMMultLieb2DAtoB1

! --------------------------------------------------------------------
! TMMultLieb2DBtoA:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_B,PSI_A) so that the structure of the transfer matrix 
! can be exploited

SUBROUTINE TMMultLieb2DB2toA(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

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
  
  !PRINT*,"DBG: TMMultLieb2DB2toA()"
  
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
        
        !PRINT*,"i,j,En, OP, PL, PR, PA,PB, PN"
        !PRINT*, iSite, jState, En, OnsitePot, PsiLeft, PsiRight,
        !        PSI_A(iSite,jState), PSI_B(iSite,jState),
        !        new
        
        PSI_B(iSite,jState)= new
        
     ENDDO ! jState
  ENDDO ! iSite
  
  RETURN
END SUBROUTINE TMMultLieb2DB2toA

! --------------------------------------------------------------------
! TMMultLieb2DBtoA:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_B,PSI_A) so that the structure of the transfer matrix 
! can be exploited

SUBROUTINE TMMultLieb2DB1toB2(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

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

  CALL TMMultLieb2DB2toA(PSI_A,PSI_B, Ilayer, En, DiagDis, M )
    
!!$  DO iSite=1,M
!!$     
!!$     ! create the new onsite potential
!!$     SELECT CASE(IRNGFlag)
!!$     CASE(0)
!!$        OnsitePot= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
!!$     CASE(1)
!!$        OnsitePot= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
!!$     CASE(2)
!!$        OnsitePot= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
!!$     END SELECT
!!$     
!!$     !PRINT*,"iS,pL,RndVec", iSite,pLevel,RndVec((pLevel-1)*M+iSite)
!!$     
!!$     DO jState=1,M
!!$        
!!$        !PRINT*,"jState, iSite", jState, iSite,
!!$        
!!$        new= ( OnsitePot * PSI_A(jState,iSite) &
!!$             - PSI_B(jState,iSite) )
!!$        
!!$        !PRINT*,"i,j,En, OP, PL, PR, PA,PB, PN"
!!$        !PRINT*, iSite, jState, En, OnsitePot, PsiLeft, PsiRight,
!!$        !        PSI_A(iSite,jState), PSI_B(iSite,jState),
!!$        !        new
!!$        
!!$        PSI_B(jState,iSite)= new
!!$        
!!$     ENDDO ! jState
!!$  ENDDO ! iSite
!!$  
!!$  !PRINT*,"PSIA(1,1),(1,2),(1,3),(1,4)",&
!!$        !PSI_A(1,1),PSI_A(1,2),PSI_A(1,3),PSI_A(1,4)
!!$
!!$  !PRINT*,"PSIB(1,1),(1,2),(1,3),(1,4)",&
!!$        !PSI_B(1,1),PSI_B(1,2),PSI_B(1,3),PSI_B(1,4)
  
  RETURN
END SUBROUTINE TMMultLieb2DB1toB2
