! --------------------------------------------------------------------
! TMMultLieb2DAtoB:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_B,PSI_A) so that the structure of the transfer matrix 
! can be exploited

SUBROUTINE TMMultLieb2DAtoB(PSI_A,PSI_B, Ilayer, En, DiagDis, RimDis, CubeConPot, LiebConstPoten, M )

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
       RimDis, CubeConPot, LiebConstPoten, &
       En                    ! energy
  
  REAL(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  
  INTEGER xSiteL,xSiteS, jState, ISeedDummy
  REAL(KIND=RKIND) OnsitePot, OnsiteRight, OnsiteLeft, OnsitePotVec(2*M)
  REAL(KIND=CKIND) new , PsiLeft, PsiRight, stub, Kappa
  
  INTEGER, PARAMETER :: LiebSpacer=2

  !PRINT*,"DBG: TMMultLieb2DAtoB()"

  ! create the new onsite potential
  !  IRNGFlag=(xy)
  !  xy=0x      presenting conditions that all positions are disorder
  !  xy=x0,     presenting conditions that only central positions are disorder   
  
  DO xSiteS=1,LiebSpacer*M   
     SELECT CASE(IRNGFlag)
     CASE(01)
        OnsitePotVec(xSiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
     CASE(02)
        OnsitePotVec(xSiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
     CASE(03)
        OnsitePotVec(xSiteS)= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
     CASE(10)
        IF(Mod(xSiteS,LiebSpacer)==1) THEN
           OnsitePotVec(xSiteS)= -En + CubeConPot + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
        ELSE
           OnsitePotVec(xSiteS)= -En + LiebConPot + RimDis*(DRANDOM(ISeedDummy)-0.5D0)
        END IF
     CASE(20)
        IF(Mod(xSiteS,LiebSpacer)==1) THEN
           OnsitePotVec(xSiteS)= -En + CubeConPot + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
        ELSE
           OnsitePotVec(xSiteS)= -En + LiebConPot + RimDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
        END IF
     CASE(30)
        IF(Mod(xSiteS,LiebSpacer)==1) THEN
           OnsitePotVec(xSiteS)= -En + CubeConPot + GRANDOM(ISeedDummy,0.0D0,DiagDis)
        ELSE
           OnsitePotVec(xSiteS)= -En + LiebConPot + GRANDOM(ISeedDummy,0.0D0,RimDis)
        END IF
     END SELECT
  END DO
     
  ! to the TMM
  DO xSiteL=1,M
     
     xSiteS= (xSiteL-1)*LiebSpacer + 1
     
     OnsitePot= OnsitePotVec(xSiteS) 
     
     DO jState=1,M
        
        !PsiLeft
        IF (xSiteL.EQ.1) THEN
           SELECT CASE(IBCFlag)
           CASE(-1,0) ! hard wall BC
              PsiLeft= ZERO            
              OnsiteLeft= ZERO
           CASE(1) ! periodic BC
              stub= OnsitePotVec(LiebSpacer*M)
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteLeft= 1.D0 /stub
              PsiLeft= PSI_A(M,jState) /stub
           CASE(2) ! antiperiodic BC
              stub= OnsitePotVec(LiebSpacer*M)
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteLeft= 1.D0 /stub
              PsiLeft= -PSI_A(M,jState) /stub
           CASE DEFAULT
              PRINT*,"TMMultLieb2DAtoB(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
           END SELECT
        ELSE
           stub= OnsitePotVec(xSiteS -1)
           IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
           OnsiteLeft= 1.D0 /stub
           PsiLeft= PSI_A(xSiteL-1,jState) /stub
        END IF

        !PsiRight
        IF (xSiteL.EQ.M) THEN
           SELECT CASE(IBCFlag)
           CASE(-1) ! hard wall BC + STUBS
              stub= OnsitePotVec(xSiteS +1)
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteRight= 1.D0 /stub
              PsiRight= ZERO            
           CASE(0) ! hard wall BC
              PsiRight= ZERO            
              OnsiteRight= ZERO
           CASE(1) ! periodic BC
              stub= OnsitePotVec(xSiteS +1)
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteRight= 1.D0 /stub
              PsiRight= PSI_A(1,jState) /stub
           CASE(2) ! antiperiodic BC
              stub= OnsitePotVec(xSiteS +1)
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteRight= 1.D0 /stub
              PsiRight= -PSI_A(1,jState) /stub
           CASE DEFAULT
              PRINT*,"TMMultLieb2DAtoB(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
           END SELECT
        ELSE
!!$           PsiRight= PSI_A(xSiteL+1,jState)/OnsitePotVec(xSiteS +1)
!!$           OnsiteRight= 1.D0/OnsitePotVec(xSiteS +1)
           stub=OnsitePotVec(xSiteS +1)
           IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
           PsiRight= PSI_A(xSiteL+1,jState)/stub
           OnsiteRight= 1.D0/stub
        END IF
        
        new =(( OnsitePot - OnsiteLeft - OnsiteRight ) * PSI_A(xSiteL,jState) &
             - Kappa * ( PsiLeft + PsiRight ) &
             - PSI_B(xSiteL,jState) )
        
        !PRINT*,"i,j,En, OP, PL, PR, PA,PB, PN"
        !PRINT*, xSite, jState, En, OnsitePot, PsiLeft, PsiRight,
        !        PSI_A(xSite,jState), PSI_B(xSite,jState),
        !        new
        
        PSI_B(xSiteL,jState)= new
        
     ENDDO ! jState
  ENDDO ! xSite
  
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

SUBROUTINE TMMultLieb2DBtoA(PSI_A,PSI_B, Ilayer, En, DiagDis, RimDis, CubeConPot, LiebConPot, M )

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
       RimDis, CubeConPot, LiebConPot, &
       En                    ! energy
  
  REAL(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  
  INTEGER xSite, jState, ISeedDummy
  REAL(KIND=RKIND) OnsitePot
  REAL(KIND=CKIND) new
  
  !PRINT*,"DBG: TMMultLieb2DBtoA()"
  
  DO xSite=1,M
     
     ! create the new onsite potential
     SELECT CASE(IRNGFlag)
     CASE(01)
        OnsitePot= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
     CASE(02)
        OnsitePot= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
     CASE(03)
        OnsitePot= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
     CASE(10)
        OnsitePot= -En + LiebConPot + RimDis*(DRANDOM(ISeedDummy)-0.5D0)
     CASE(20)
        OnsitePot= -En + LiebConPot + RimDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
     CASE(30)
        OnsitePot= -En + LiebConPot + GRANDOM(ISeedDummy,0.0D0,RimDis)
     END SELECT
     
     !PRINT*,"iS,pL,RndVec", xSite,pLevel,RndVec((pLevel-1)*M+xSite)
     
     DO jState=1,M
        
        !PRINT*,"jState, xSite", jState, xSite,
        
        new= ( OnsitePot * PSI_A(xSite,jState) &
             - PSI_B(xSite,jState) )
        
        !PRINT*,"i,j,En, OP, PL, PR, PA,PB, PN"
        !PRINT*, xSite, jState, En, OnsitePot, PsiLeft, PsiRight,
        !        PSI_A(xSite,jState), PSI_B(xSite,jState),
        !        new
        
        PSI_B(xSite,jState)= new
        
     ENDDO ! jState
  ENDDO ! xSite
  
  !PRINT*,"PSIA(1,1),(1,2),(1,3),(1,4)",&
        !PSI_A(1,1),PSI_A(1,2),PSI_A(1,3),PSI_A(1,4)

  !PRINT*,"PSIB(1,1),(1,2),(1,3),(1,4)",&
        !PSI_B(1,1),PSI_B(1,2),PSI_B(1,3),PSI_B(1,4)
  
  RETURN
END SUBROUTINE TMMultLieb2DBtoA

