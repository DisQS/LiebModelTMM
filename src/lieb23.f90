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
  
  INTEGER xSiteS,xSiteL, jSite,jState, indexK,ISeedDummy
  REAL(KIND=RKIND) OnsitePot, OnsiteLeft, OnsiteRight, OnsitePotVec(4*M,1)
  REAL(KIND=CKIND) new , PsiLeft, PsiRight,stub

  INTEGER, PARAMETER :: LiebSpacer=4

  ! create the new onsite potential
   jSite=1

  DO xSiteS=1,LiebSpacer*M   !label different row
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePotVec(xSiteS,jSite)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
     CASE(1)
        OnsitePotVec(xSiteS,jSite)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
     CASE(2)
        OnsitePotVec(xSiteS,jSite)= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
  END DO

  ! to the TMM
  DO xSiteL=1,M  ! lable the the row of A atom.
     
     !indexK=(xSiteS/4+1)      !label the A atom

     xSiteS= (xSiteL-1)*LiebSpacer + 1

     OnsitePot=OnsitePotVec(xSiteS,jSite)
     
     DO jState=1,M

        !PsiLeft
        IF (xSiteL.LE.1) THEN
           SELECT CASE(IBCFlag)
           CASE(-1,0) ! hard wall BC
              OnsiteLeft=ZERO      
              PsiLeft=ZERO
           CASE(1) ! periodic BC
              stub= (OnsitePotVec(LiebSpacer*M,jSite)*OnsitePotVec(LiebSpacer*M-1,jSite)*OnsitePotVec(LiebSpacer*M-2,jSite) &
                   -OnsitePotVec(LiebSpacer*M,jSite)-OnsitePotVec(LiebSpacer*M-2,jSite))
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteLeft= (OnsitePotVec(LiebSpacer*M-1,jSite)*OnsitePotVec(LiebSpacer*M-2,jSite)-1.0D0) /stub
              PsiLeft= PSI_A(M,jState) /stub
           !CASE(2)
           CASE DEFAULT
              PRINT*,"TMMultLieb2DAtoB1(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
           END SELECT
        ELSE
           stub= (OnsitePotVec(xSiteS-1,jSite)*OnsitePotVec(xSiteS-2,jSite)*OnsitePotVec(xSiteS-3,jSite) &
                -OnsitePotVec(xSiteS-1,jSite)-OnsitePotVec(xSiteS-3,jSite))
           IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
           OnsiteLeft= (OnsitePotVec(xSiteS-2,jSite)*OnsitePotVec(xSiteS-3,jSite)-1.0D0) /stub
           PsiLeft= PSI_A(xSiteL-1,jState) /stub 
        END IF
        
        !PsiRight
        IF (xSiteL.GE.M) THEN
           SELECT CASE(IBCFlag)
           CASE(-1) ! hard wall BC + stubs
              stub= (OnsitePotVec(xSiteS+1,jSite)*OnsitePotVec(xSiteS+2,jSite)*OnsitePotVec(xSiteS+3,jSite) &
                   -OnsitePotVec(xSiteS+1,jSite)-OnsitePotVec(xSiteS+3,jSite))
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteRight= (OnsitePotVec(xSiteS+2,jSite)*OnsitePotVec(xSiteS+3,jSite)-1.0D0) /stub 
              PsiRight= ZERO
           CASE(0) ! hard wall BC 
              OnsiteRight=ZERO
              PsiRight=ZERO
           CASE(1) ! periodic BC
              stub= (OnsitePotVec(xSiteS+1,jSite)*OnsitePotVec(xSiteS+2,jSite)*OnsitePotVec(xSiteS+3,jSite) &
                   -OnsitePotVec(xSiteS+1,jSite)-OnsitePotVec(xSiteS+3,jSite))
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              OnsiteRight=(OnsitePotVec(xSiteS+2,jSite)*OnsitePotVec(xSiteS+3,jSite)-1.0D0) /stub  
              PsiRight=PSI_A(1,jState) /stub 
           !CASE(2) ! antiperiodic BC
           CASE DEFAULT
              PRINT*,"TMMultLieb2DAtoB1(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
           END SELECT
        ELSE
           stub= (OnsitePotVec(xSiteS+1,jSite)*OnsitePotVec(xSiteS+2,jSite)*OnsitePotVec(xSiteS+3,jSite) &
                -OnsitePotVec(xSiteS+1,jSite)-OnsitePotVec(xSiteS+3,jSite))
           IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
           OnsiteRight=(OnsitePotVec(xSiteS+2,jSite)*OnsitePotVec(xSiteS+3,jSite)-1.0D0) /stub  
           PsiRight=PSI_A(xSiteL+1,jState) /stub 
        END IF
        
        new =(( OnsitePot - OnsiteLeft - OnsiteRight )*PSI_A(xSiteL,jState) &
             - Kappa * ( PsiLeft + PsiRight ) &
             - PSI_B(xSiteL,jState) )

        PSI_B(xSiteL,jState)= new
        
     ENDDO ! jState
  ENDDO ! xSiteS
  
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
  
  INTEGER xSite, jState, ISeedDummy
  REAL(KIND=RKIND) OnsitePot
  REAL(KIND=CKIND) new
  
  DO xSite=1,M
     
     ! create the new onsite potential
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePot= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
     CASE(1)
        OnsitePot= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
     CASE(2)
        OnsitePot= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT
     
     !PRINT*,"iS,pL,RndVec", xSite,pLevel,RndVec((pLevel-1)*M+xSite)
     
     DO jState=1,M
        
        !PRINT*,"jState, xSite", jState, xSite,
        
        new= ( OnsitePot * PSI_A(xSite,jState) &
             - PSI_B(xSite,jState) )
        
        PSI_B(xSite,jState)= new
        
     ENDDO ! jState
  ENDDO ! xSite
  
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
  
  !PRINT*,"DBG: TMMultLieb2DB1toB2()"

  CALL TMMultLieb2D_B1toB2(PSI_A,PSI_B, Ilayer, En, DiagDis, M )
  RETURN
END SUBROUTINE TMMultLieb2D_B3toA
