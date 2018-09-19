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
  
  REAL(KIND=RKIND) PSI_A(M*M,M*M), PSI_B(M*M,M*M)
  
  INTEGER iSite, jState, ISeedDummy
  REAL(KIND=RKIND) OnsitePot
  REAL(KIND=RKIND) new, PsiLeft, PsiRight, PsiUp, PsiDown
  
  !PRINT*,"DBG: TMMultLieb3DAtoB()"
  
  DO iSite=1,M*M
     
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
     
     DO jState=1,M*M
        
        !PRINT*,"jState, iSite", jState, iSite,
         
        ! PsiLeft
        !IF (iSite.EQ.1) THEN
        IF (MOD(iSite,M).EQ.1) THEN
           
           IF (IBCFlag.EQ.0) THEN
              PsiLeft= CZERO            ! hard wall BC
           ELSE IF (IBCFlag.EQ.1) THEN
              PsiLeft= PSI_A(jState,iSite+M-1)  ! periodic BC
           ELSE IF (IBCFlag.EQ.2) THEN
              PsiLeft= -PSI_A(jState,iSite+M-1) ! antiperiodic BC
           ENDIF
           
        ELSE
           !PRINT*,"DBG: iSite=",iSite
           PsiLeft= PSI_A(jState,iSite-1)
        ENDIF

        ! PsiRight        
        !IF (iSite.EQ.M) THEN
        IF (MOD(iSite,M).EQ.0) THEN
           
           IF (IBCFlag.EQ.0) THEN
              PsiRight= CZERO            ! hard wall BC
           ELSE IF (IBCFlag.EQ.1) THEN
              PsiRight= PSI_A(jState,iSite-M+1)  ! periodic BC
           ELSE IF (IBCFlag.EQ.2) THEN
              PsiRight= -PSI_A(jState,iSite-M+1) ! antiperiodic BC
           ENDIF
           
        ELSE
           PsiRight= PSI_A(jState,iSite+1)
        ENDIF

        ! PsiUp
        IF (iSite.GT.(M-1)*M) THEN
           
           IF (IBCFlag.EQ.0) THEN
              PsiUp= CZERO                ! hard wall BC
           ELSE IF (IBCFlag.EQ.1) THEN
              PsiUp= PSI_A(jState,iSite-(M-1)*M)  ! periodic BC
           ELSE IF (IBCFlag.EQ.2) THEN
              PsiUp= -PSI_A(jState,iSite-(M-1)*M) ! antiperiodic BC
           ENDIF
           
        ELSE
           PsiUp= PSI_A(jState,iSite+M)
        ENDIF

        ! PsiDown        
        IF (iSite.LT.(M+1)) THEN
           
           IF (IBCFlag.EQ.0) THEN
              PsiDown= CZERO            ! hard wall BC
           ELSE IF (IBCFlag.EQ.1) THEN
              PsiDown= PSI_A(jState,iSite+(M-1)*M)  ! periodic BC
           ELSE IF (IBCFlag.EQ.2) THEN
              PsiDown= -PSI_A(jState,iSite+(M-1)*M) ! antiperiodic BC
           ENDIF
           
        ELSE
           PsiDown= PSI_A(jState,iSite-M)
        ENDIF        
        
        new= ( OnsitePot * PSI_A(jState,iSite) &
             - Kappa * ( PsiLeft + PsiRight + PsiUp + PsiDown  ) &
             - PSI_B(jState,iSite) )
        
        !PRINT*,"i,j,En, OP, PL, PR, PA,PB, PN"
        !PRINT*, iSite, jState, En, OnsitePot, PsiLeft, PsiRight,
        !        PSI_A(iSite,jState), PSI_B(iSite,jState),
        !        new
        
        PSI_B(jState,iSite)= new
        
     ENDDO ! jState
  ENDDO ! iSite
  
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
  
  REAL(KIND=RKIND) PSI_A(M*M,M*M), PSI_B(M*M,M*M)
  
  INTEGER iSite, jState, ISeedDummy
  REAL(KIND=RKIND) OnsitePot
  REAL(KIND=RKIND) new, PsiLeft, PsiRight, PsiUp, PsiDown
  
  !PRINT*,"DBG: TMMultLieb3DBtoA()"
  
  DO iSite=1,M*M
     
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
     
     DO jState=1,M*M
        
        !PRINT*,"jState, iSite", jState, iSite,
         
        ! PsiLeft
        !IF (iSite.EQ.1) THEN
        IF (MOD(iSite,M).EQ.1) THEN
           
           IF (IBCFlag.EQ.0) THEN
              PsiLeft= CZERO            ! hard wall BC
           ELSE IF (IBCFlag.EQ.1) THEN
              PsiLeft= PSI_A(jState,iSite+M-1)  ! periodic BC
           ELSE IF (IBCFlag.EQ.2) THEN
              PsiLeft= -PSI_A(jState,iSite+M-1) ! antiperiodic BC
           ENDIF
           
        ELSE
           !PRINT*,"DBG: iSite=",iSite
           PsiLeft= PSI_A(jState,iSite-1)
        ENDIF

        ! PsiRight        
        !IF (iSite.EQ.M) THEN
        IF (MOD(iSite,M).EQ.0) THEN
           
           IF (IBCFlag.EQ.0) THEN
              PsiRight= CZERO            ! hard wall BC
           ELSE IF (IBCFlag.EQ.1) THEN
              PsiRight= PSI_A(jState,iSite-M+1)  ! periodic BC
           ELSE IF (IBCFlag.EQ.2) THEN
              PsiRight= -PSI_A(jState,iSite-M+1) ! antiperiodic BC
           ENDIF
           
        ELSE
           PsiRight= PSI_A(jState,iSite+1)
        ENDIF

        ! PsiUp
        IF (iSite.GT.(M-1)*M) THEN
           
           IF (IBCFlag.EQ.0) THEN
              PsiUp= CZERO                ! hard wall BC
           ELSE IF (IBCFlag.EQ.1) THEN
              PsiUp= PSI_A(jState,iSite-(M-1)*M)  ! periodic BC
           ELSE IF (IBCFlag.EQ.2) THEN
              PsiUp= -PSI_A(jState,iSite-(M-1)*M) ! antiperiodic BC
           ENDIF
           
        ELSE
           PsiUp= PSI_A(jState,iSite+M)
        ENDIF

        ! PsiDown        
        IF (iSite.LT.(M+1)) THEN
           
           IF (IBCFlag.EQ.0) THEN
              PsiDown= CZERO            ! hard wall BC
           ELSE IF (IBCFlag.EQ.1) THEN
              PsiDown= PSI_A(jState,iSite+(M-1)*M)  ! periodic BC
           ELSE IF (IBCFlag.EQ.2) THEN
              PsiDown= -PSI_A(jState,iSite+(M-1)*M) ! antiperiodic BC
           ENDIF
           
        ELSE
           PsiDown= PSI_A(jState,iSite-M)
        ENDIF        
        
        new= ( OnsitePot * PSI_A(jState,iSite) &
             - Kappa * ( PsiLeft + PsiRight + PsiUp + PsiDown  ) &
             - PSI_B(jState,iSite) )
        
        !PRINT*,"i,j,En, OP, PL, PR, PA,PB, PN"
        !PRINT*, iSite, jState, En, OnsitePot, PsiLeft, PsiRight,
        !        PSI_A(iSite,jState), PSI_B(iSite,jState),
        !        new
        
        PSI_B(jState,iSite)= new
        
     ENDDO ! jState
  ENDDO ! iSite
  
  !PRINT*,"PSIA(1,1),(2,1),(1,2),(2,2)",&
  !      PSI_A(1,1),PSI_A(2,1),PSI_A(1,2),PSI_A(2,2)
  
  RETURN
END SUBROUTINE TMMultLieb3DBtoA

  
