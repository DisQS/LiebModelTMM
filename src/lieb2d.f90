! --------------------------------------------------------------------
! TMMultLieb2DAtoB:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_B,PSI_A) so that the structure of the transfer matrix 
! can be exploited

SUBROUTINE TMMultLieb2DAtoB(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

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
  REAL(KIND=RKIND) OnsitePot, OnsiteRight, OnsiteLeft, OnsitePotVec(2*M)
  REAL(KIND=CKIND) new, PsiLeft, PsiRight
  
  !PRINT*,"DBG: TMMultLieb2DAtoB()"

  ! create the new onsite potential
  DO iSite=1,M   
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePotVec(iSite)= En - DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
     CASE(1)
        OnsitePotVec(iSite)= En - DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
     CASE(2)
        OnsitePotVec(iSite)= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
  END DO
     
  ! to the TMM
  DO iSite=1,M
     
     !PRINT*,"iS,pL,RndVec", iSite,pLevel,RndVec((pLevel-1)*M+iSite)

     OnsitePot= &
          OnsitePotVec(iSite) !+ 1.D0/(OnsitePotVec(iSite-1) + 1.D0/(OnsitePotVec(iSite+1)
     
     DO jState=1,M
        
        !PRINT*,"jState, iSite", jState, iSite,
        
        SELECT CASE(iSite)
        CASE(1)
           IF (IBCFlag.EQ.0) THEN
              PsiLeft= ZERO            ! hard wall BC
              OnsiteLeft= ZERO
           ELSE IF (IBCFlag.EQ.1) THEN
              PsiLeft= PSI_A(jState,M)  ! periodic BC
              OnsiteLeft= OnsitePotVec(2*M)
           ELSE IF (IBCFlag.EQ.2) THEN
              PsiLeft= -PSI_A(jState,M) ! antiperiodic BC
              OnsiteLeft= OnsitePotVec(2*M)
           ENDIF
        CASE DEFAULT
           PsiLeft= PSI_A(jState,iSite-1)/OnsitePotVec(2*iSite-1 -1)
           OnsiteLeft= OnsiteRight + 1.D0/OnsitePotVec(2*(iSite-1))
        END SELECT
        
        SELECT CASE(iSite)
        CASE(1)
           IF (IBCFlag.EQ.0) THEN
              PsiRight= ZERO            ! hard wall BC
              OnsiteRight= ZERO
           ELSE IF (IBCFlag.EQ.1) THEN
              PsiRight= PSI_A(jState,1)  ! periodic BC
              OnsiteRight= OnsitePotVec(2*M)
           ELSE IF (IBCFlag.EQ.2) THEN
              PsiRight= -PSI_A(jState,1) ! antiperiodic BC
              OnsiteRight= OnsitePotVec(2*M)
           ENDIF
        CASE DEFAULT
           PsiRight= PSI_A(jState,iSite+1)/OnsitePotVec(2*iSite-1 +1)
           OnsiteRight= OnsiteRight + 1.D0/OnsitePotVec(2*(iSite+1))
        END SELECT
        
        new= ( OnsitePot * PSI_A(jState,iSite) &
             + Kappa * ( PsiLeft + PsiRight ) &
             - PSI_B(jState,iSite) )
        
        !PRINT*,"i,j,En, OP, PL, PR, PA,PB, PN"
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
END SUBROUTINE TMMultLieb2DAtoB

! --------------------------------------------------------------------
! TMMultLieb2DBtoA:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_B,PSI_A) so that the structure of the transfer matrix 
! can be exploited

SUBROUTINE TMMultLieb2DBtoA(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

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
  
  !PRINT*,"DBG: TMMultLieb2DBtoA()"
  
  DO iSite=1,M
     
     ! create the new onsite potential
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot= En - DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
     CASE(1)
        OnsitePot= En - DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
     CASE(2)
        OnsitePot= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
     
     !PRINT*,"iS,pL,RndVec", iSite,pLevel,RndVec((pLevel-1)*M+iSite)
     
     DO jState=1,M
        
        !PRINT*,"jState, iSite", jState, iSite,
        
        new= ( OnsitePot * PSI_A(jState,iSite) &
             - PSI_B(jState,iSite) )
        
        !PRINT*,"i,j,En, OP, PL, PR, PA,PB, PN"
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
END SUBROUTINE TMMultLieb2DBtoA

