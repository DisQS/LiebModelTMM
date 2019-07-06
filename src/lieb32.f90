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

  INTEGER jState, ISeedDummy,iSiteS,jSiteS, iSiteL,jSiteL
  REAL(KIND=RKIND) OnsitePot, OnsiteRight, OnsiteLeft, OnsiteUp, OnsiteDown
  REAL(KIND=RKIND) NEW, PsiLeft, PsiRight, PsiUp, PsiDown, stub

  INTEGER Coord2IndexL
  EXTERNAL Coord2IndexL

  !PRINT*,"DBG: TMMultLieb3DAtoB()"

  ! create the new onsite potential
  DO iSiteS=1,3*M
     DO jSiteS=1,3*M

        !indexS= (iSiteS-1)*3*M + jSiteS

        SELECT CASE(IRNGFlag)
        CASE(0)
           OnsitePotVec(iSiteS,jSiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
        CASE(1)
           OnsitePotVec(iSiteS,jSiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
        CASE(2)
           OnsitePotVec(iSiteS,jSiteS)= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
        END SELECT
     END DO
  END DO

  !PRINT*,"iS,pL,RndVec", iSite,pLevel,RndVec((pLevel-1)*M+iSite)

  !new TMM
  DO iSiteL=1,M
     DO jSiteL=1,M

        iSiteS= (iSiteL-1)*3 + 1
        jSiteS= (jSiteL-1)*3 + 1
        
        !indexL= (jSiteL-1)*M + iSiteL
        !indexS= (jSiteS-1)*M + iSiteS
        
!!$        PRINT*,"iSL,jSL, iSS, jSS, C2I", &
!!$             iSiteL,jSiteL, iSiteS,jSiteS, Coord2IndexL(M,iSiteL,jSiteL)

        OnsitePot=OnsitePotVec(iSiteS,jSiteS)

        DO jState=1,M*M
           
           !PsiLeft
           IF (iSiteL.LE.1) THEN
              SELECT CASE(IBCFLag)
              CASE(-1,0) ! hard wall BC
                 OnsiteLeft= ZERO
                 PsiLeft= ZERO
              CASE(1) ! periodic BC
                 stub= OnsitePotVec(3*M,jSiteS)*OnSitePotVec(3*M-1,jSiteS)-1.0D0
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteLeft= OnsitePotVec(3*M-1,jSiteS) /stub
                 PsiLeft= Psi_A(Coord2IndexL(M,M,jSiteL),jState) /stub
!              CASE(2) ! antiperiodic BC
              CASE DEFAULT
                 PRINT*,"TMMultLieb3DAtoB5(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              stub= OnsitePotVec(iSiteS-1,jSiteS)*OnSitePotVec(iSiteS-2,jSiteS)-1.0D0
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteLeft= OnsitePotVec(iSiteS-2,jSiteS) /stub
              PsiLeft= Psi_A(Coord2IndexL(M,iSiteL-1,jSiteL),jState) /stub
           END IF

           !PsiRight
           IF (iSiteL.GE.M) THEN
              SELECT CASE(IBCFLag)
              CASE(-1) ! hard wall BC + stubs
                 stub= OnsitePotVec(iSiteS+1,jSiteS)*OnSitePotVec(iSiteS+2,jSiteS)-1.0D0
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteRight= OnsitePotVec(iSiteS+2,jSiteS)/stub
                 PsiRight= ZERO
              CASE(0) ! hard wall
                 OnsiteRight= ZERO
                 PsiRight= ZERO
              CASE(1) ! periodic BC
                 stub= OnsitePotVec(iSiteS+1,jSiteS)*OnSitePotVec(iSiteS+2,jSiteS)-1.0D0
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteRight= OnsitePotVec(iSiteS+2,jSiteS) /stub
                 PsiRight= Psi_A(Coord2IndexL(M,1,jSiteL),jState) /stub
!              CASE(2) ! antiperiodic BC
              CASE DEFAULT
                 PRINT*,"TMMultLieb3DAtoB5(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              stub= OnsitePotVec(iSiteS+1,jSiteS)*OnSitePotVec(iSiteS+2,jSiteS)-1.0D0
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteRight= OnsitePotVec(iSiteS+2,jSiteS) /stub
              PsiRight= Psi_A(Coord2IndexL(M,iSiteL+1,jSiteL),jState) /stub
           END IF

           !PsiDown
           IF (jSiteL.GE.M) THEN
              SELECT CASE(IBCFlag)
              CASE(-1) ! hard wall BC + stubs
                 stub= OnsitePotVec(iSiteS,jSiteS+1)*OnSitePotVec(iSiteS,jSiteS+2)-1.0D0
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteDown= OnsitePotVec(iSiteS,jSiteS+2) /stub
                 PsiDown= ZERO
              CASE(0) ! hard wall
                 OnsiteUp=ZERO
                 PsiDown= ZERO
              CASE(1) ! periodic BC
                 stub= OnsitePotVec(iSiteS,jSiteS+1)*OnSitePotVec(iSiteS,jSiteS+2)-1.0D0
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteDown= OnsitePotVec(iSiteS,jSiteS+2) /stub
                 PsiDown= Psi_A(Coord2IndexL(M,iSiteL,1),jState) /stub
!              CASE(2) ! antiperiodic BC
              CASE DEFAULT
                 PRINT*,"TMMultLieb3DAtoB5(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              stub= OnsitePotVec(iSiteS,jSiteS+1)*OnSitePotVec(iSiteS,jSiteS+2)-1.0D0
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteDown= OnsitePotVec(iSiteS,jSiteS+2) /stub
              PsiDown= Psi_A(Coord2IndexL(M,iSiteL,jSiteL+1),jState) /stub
           END IF

           !PsiUp
           IF (jSiteL.LE.1) THEN
              SELECT CASE(IBCFlag)
              CASE(-1,0) ! hard wall BC
                 OnsiteUp= ZERO
                 PsiUp= ZERO
              CASE(1) ! periodic BC
                 stub= (OnsitePotVec(iSiteS,3*M)*OnSitePotVec(iSiteS,3*M-1)-1.0D0)
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 OnsiteUp= OnsitePotVec(iSiteS,3*M-1) /stub
                 PsiUp=  Psi_A(Coord2IndexL(M,iSiteL,M),jState) /stub
              CASE(2) ! antiperiodic BC
              CASE DEFAULT
                 PRINT*,"TMMultLieb3DAtoB5(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              stub= OnsitePotVec(iSiteS,jSiteS-1)*OnSitePotVec(iSiteS,jSiteS-2)-1.0D0
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteUp= OnsitePotVec(iSiteS,jSiteS-2) /stub
              PsiUp= Psi_A(Coord2IndexL(M,iSiteL,jSiteL-1),jState) /stub
           END IF

           NEW= ( OnsitePot - OnsiteLeft - OnsiteRight - OnsiteUp - OnsiteDown ) * &
                Psi_A(Coord2IndexL(M,iSiteL,jSiteL),jState) &
                - Kappa * ( PsiLeft + PsiRight + PsiUp + PsiDown ) &
                - PSI_B(Coord2IndexL(M,iSiteL,jSiteL),jState)
           
           PSI_B(Coord2IndexL(M,iSiteL,jSiteL),jState)= NEW
        END DO !jState
        
     END DO !iSiteL
  END DO !jSiteL
  RETURN

END SUBROUTINE TMMultLieb3DAtoB5

!!$! --------------------------------------------------------------------
!!$! convert i,j coordinates to an index
!!$FUNCTION Coord2IndexL(isize, iSite, jSite)
!!$  INTEGER Coord2IndexL, isize, iSite, jSite
!!$  
!!$  Coord2IndexL= (jSite-1)*isize + iSite
!!$  
!!$  RETURN
!!$END FUNCTION Coord2IndexL

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
