! --------------------------------------------------------------------
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
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
  
  INTEGER iSiteL,iSiteS, jState, ISeedDummy
  REAL(KIND=RKIND) OnsitePot, OnsiteRight, OnsiteLeft, OnsitePotVec(3*M)
  REAL(KIND=CKIND) new , PsiLeft, PsiRight
  
  PRINT*,"DBG: TMMultLieb2DAtoB1()"

  ! create the new onsite potential
  DO iSiteS=1,3*M   
     SELECT CASE(IRNGFlag)
     CASE(0)
        OnsitePotVec(iSiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
     CASE(1)
        OnsitePotVec(iSiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
     CASE(2)
        OnsitePotVec(iSiteS)= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
     END SELECT

     IF( ABS(OnsitePotVec(iSiteS)).LT.TINY) THEN
        OnsitePotVec(iSiteS)= SIGN(TINY,OnsitePotVec(iSiteS))
     END IF
  END DO
     
  ! to the TMM
  DO iSiteL=1,M
     
     !PRINT*,"iS,pL,RndVec", iSite,pLevel,RndVec((pLevel-1)*M+iSite)
     iSiteS= 3*iSiteL-2
     
     OnsitePot= &
          OnsitePotVec(iSiteS) !+ 1.D0/(OnsitePotVec(iSite-1) + 1.D0/(OnsitePotVec(iSite+1)    
     
     DO jState=1,M
        
        !PRINT*,"jState, iSite", jState, iSite,

        IF (iSiteL.EQ.1) THEN

           IF (IBCFlag.EQ.0) THEN
              PsiLeft= ZERO            ! hard wall BC
              OnsiteLeft= ZERO
           ELSE IF (IBCFlag.EQ.1) THEN
              PsiLeft= 1.0D0*PSI_A(jState,M)    &
                     /( OnsitePotVec(3*M -1)*OnsitePotVec(3*M)-1.0D0 )  ! periodic BC
              OnsiteLeft= 1.0d0*OnsitePotVec(3*M -1)      &
                        /( OnsitePotVec(3*M -1)*OnsitePotVec(3*M) - 1.0D0 )
           ELSE IF (IBCFlag.EQ.2) THEN
              PsiLeft= 1.0D0*PSI_A(jState,M)    &
                     /( OnsitePotVec(3*M -1)*OnsitePotVec(3*M)-1.0D0 ) ! antiperiodic BC
              OnsiteLeft= 1.0d0*OnsitePotVec(3*M -1)      &
                        /( OnsitePotVec(3*M -1)*OnsitePotVec(3*M) - 1.0D0 )
           ENDIF
        ELSE
           PsiLeft=1.0D0*PSI_A(jState,iSiteL-1)    &
                   /( OnsitePotVec(iSiteS -2)*OnsitePotVec(iSiteS -1) - 1.0D0 )
           OnsiteLeft= 1.0d0*OnsitePotVec(isiteS -2)      &
                      /( OnsitePotVec(iSiteS -2)*OnsitePotVec(iSiteS -1) - 1.0D0 )
        END IF

        IF (iSiteL.EQ.M) THEN

           IF (IBCFlag.EQ.0) THEN
              PsiRight= ZERO            ! hard wall BC
              OnsiteRight= 1.0d0*OnsitePotVec(iSiteS +2)     &
                       /( OnsitePotVec(iSiteS +2)*OnsitePotVec(iSiteS +1) - 1.0D0 )
           ELSE IF (IBCFlag.EQ.1) THEN
              PsiRight= 1.0D0*PSI_A(jState,1)    &
                    /( OnsitePotVec(iSiteS +1)*OnsitePotVec(iSiteS +2) - 1.0D0 )    ! periodic BC
              OnsiteRight= 1.0d0*OnsitePotVec(iSiteS +2)     &
                       /( OnsitePotVec(iSiteS +2)*OnsitePotVec(iSiteS +1) - 1.0D0 )
           ELSE IF (IBCFlag.EQ.2) THEN
              PsiRight= - 1.0D0*PSI_A(jState,1)    &
                    /( OnsitePotVec(iSiteS +1)*OnsitePotVec(iSiteS +2)-1.0D0 )! antiperiodic BC
              OnsiteRight=1.0d0*OnsitePotVec(iSiteS +2)     &
                       /( OnsitePotVec(iSiteS +2)*OnsitePotVec(iSiteS +1) - 1.0D0 )
           ENDIF
        ELSE
           PsiRight= 1.0D0*PSI_A(jState,iSiteL+1)    &
                    /( OnsitePotVec(iSiteS +2)*OnsitePotVec(iSiteS +1)-1.0D0 )           
           OnsiteRight= 1.0d0*OnsitePotVec(iSiteS +2)     &
                       /( OnsitePotVec(iSiteS +2)*OnsitePotVec(iSiteS +1) - 1.0D0 )
        END IF
        
        new =(( OnsitePot+OnsiteLeft+OnsiteRight )/PSI_A(jState,iSiteL) &
             - Kappa * ( PsiLeft + PsiRight ) &
             - PSI_B(jState,iSiteL) )
        
        !PRINT*,"i,j,En, OP, PL, PR, PA,PB, PN"
        !PRINT*, iSite, jState, En, OnsitePot, PsiLeft, PsiRight,
        !        PSI_A(iSite,jState), PSI_B(iSite,jState),
        !        new
        
        PSI_B(jState,iSiteL)= new
        
     ENDDO ! jState
  ENDDO ! iSite
  
  !PRINT*,"PSIA(1,1),(1,2),(1,3),(1,4)",&
        !PSI_A(1,1),PSI_A(1,2),PSI_A(1,3),PSI_A(1,4)

  !PRINT*,"PSIB(1,1),(1,2),(1,3),(1,4)",&
        !PSI_B(1,1),PSI_B(1,2),PSI_B(1,3),PSI_B(1,4)
  
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
  
  PRINT*,"DBG: TMMultLieb2DB2toA()"
  
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
  
  PRINT*,"DBG: TMMultLieb2DB1toB2()"

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
