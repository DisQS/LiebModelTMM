! --------------------------------------------------------------------
! TMMultLieb3DAtoB:
!
! 3D version of TMMult2D. Extra boundary conditions

SUBROUTINE TMMultLieb3D_AtoD1(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

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

  REAL(KIND=RKIND) PSI_A(M*M,M*M),PSI_B(M*M,M*M),OnsitePotVec(4*M,4*M)
  
  INTEGER iSiteL,jSiteL,iSiteS,jSiteS,indexK,jState, ISeedDummy
  REAL(KIND=RKIND) OnsitePot, OnsiteRight, OnsiteLeft, OnsiteUp, OnsiteDown
  REAL(KIND=RKIND) new, PsiLeft, PsiRight, PsiUp, PsiDown,stub

  INTEGER Coord2IndexL
  EXTERNAL Coord2IndexL

  !PRINT*,"DBG: TMMultLieb3DAtoB()"
     
  ! create the new onsite potential
  DO iSiteS=1,4*M
     DO jSiteS=1,4*M
        SELECT CASE(IRNGFlag)
        CASE(0)
           OnsitePotVec(iSiteS,jSiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
        CASE(1)
           OnsitePotVec(iSiteS,jSiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
        CASE(2)
           OnsitePotVec(iSiteS,jSiteS)= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
        END SELECT

        IF( ABS(OnsitePotVec(iSiteS,jSiteS)).LT.TINY) THEN
           OnsitePotVec(iSiteS,jSiteS)= SIGN(TINY,OnsitePotVec(iSiteS,jSiteS))
        END IF
     END DO
  END DO
  
  !PRINT*,"iS,pL,RndVec", iSite,pLevel,RndVec((pLevel-1)*M+iSite)

  !to the TMM
  DO iSiteL=1,M
     DO jSiteL=1,M

        iSiteS= (iSiteL-1)*4 + 1
        jSiteS= (jSiteL-1)*4 + 1

        OnsitePot=OnsitePotVec(iSiteS,jSiteS)
!        PRINT*,"DBG1: iSite, jSite, indexK", iSite, jSite, indexK
        
        DO jState=1,M*M
           
           !PsiLeft
           IF (iSiteL.EQ.1) THEN
              
              IF (IBCFlag.EQ.0) THEN
                 PsiLeft= ZERO            ! hard wall BC
                 OnsiteLeft= ZERO
                 
              ELSE IF (IBCFlag.EQ.1) THEN
!!$                 CONTINUE                 ! periodic BC
                 stub= ( OnsitePotVec(4*M-1,jSiteS)*OnsitePotVec(4*M-2,jSiteS)*OnsitePotVec(4*M,jSiteS) &
                      - OnsitePotVec(4*M,jSiteS)-OnsitePotVec(4*M-2,jSiteS) )
                 IF( ABS(stub).LT.TINY) THEN           
                    stub= SIGN(TINY,stub)
                 ENDIF
                 OnsiteLeft=( OnsitePotVec(4*M-1,jSiteS)*OnsitePotVec(4*M-2,jSiteS)-1.0D0)/stub                  
                 PsiLeft=PSI_A(Coord2IndexL(M,M,jSiteL),jState)/stub
                 
              ELSE IF (IBCFlag.EQ.2) THEN
                 CONTINUE                 ! antiperiodic BC    
              ENDIF
           ELSE
              stub= ( OnsitePotVec(iSiteS-2,jSiteS)*OnsitePotVec(iSiteS-3,jSiteS)*OnsitePotVec(iSiteS-1,jSiteS) &
                   - OnsitePotVec(iSiteS-1,jSiteS)-OnsitePotVec(iSiteS-3,jSiteS) )
              IF( ABS(stub).LT.TINY) THEN           
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteLeft=( OnsitePotVec(iSiteS-2,jSiteS)*OnsitePotVec(iSiteS-3,jSiteS)-1.0D0)/stub                  
              PsiLeft=PSI_A(Coord2IndexL(M,iSiteL-1,jSiteL),jState)/stub
           END IF
           
           !PsiRight
           IF (iSiteL.EQ.M) THEN
              
              IF (IBCFlag.EQ.0) THEN
                 stub=( OnsitePotVec(iSiteS+2,jSiteS)*OnsitePotVec(iSiteS+3,jSiteS)*OnsitePotVec(iSiteS+1,jSiteS) &
                      - OnsitePotVec(iSiteS+1,jSiteS)-OnsitePotVec(iSiteS+3,jSiteS) )
                 IF( ABS(stub).LT.TINY) THEN           
                    stub= SIGN(TINY,stub)
                 ENDIF
                 OnsiteRight=( OnsitePotVec(iSiteS+2,jSiteS)*OnsitePotVec(iSiteS+3,jSiteS)-1.0D0)/stub  ! hard wall BC                     
                 PsiRight= ZERO
                 
              ELSE IF (IBCFlag.EQ.1) THEN
!!$                 CONTINUE                                          ! periodic BC
                 stub=( OnsitePotVec(iSiteS+2,jSiteS)*OnsitePotVec(iSiteS+3,jSiteS)*OnsitePotVec(iSiteS+1,jSiteS) &
                      - OnsitePotVec(iSiteS+1,jSiteS)-OnsitePotVec(iSiteS+3,jSiteS) )
                 IF( ABS(stub).LT.TINY) THEN           
                    stub= SIGN(TINY,stub)
                 ENDIF
                 OnsiteRight=( OnsitePotVec(iSiteS+2,jSiteS)*OnsitePotVec(iSiteS+3,jSiteS)-1.0D0)/ stub
                 PsiRight=PSI_A(Coord2IndexL(M,1,jSiteL),jState)/ stub
                 
              ELSE IF (IBCFlag.EQ.2) THEN
                 CONTINUE                                          ! antiperiodic BC 
              ENDIF
           ELSE
              stub=( OnsitePotVec(iSiteS+2,jSiteS)*OnsitePotVec(iSiteS+3,jSiteS)*OnsitePotVec(iSiteS+1,jSiteS) &
                   - OnsitePotVec(iSiteS+1,jSiteS)-OnsitePotVec(iSiteS+3,jSiteS) )
              IF( ABS(stub).LT.TINY) THEN           
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteRight=( OnsitePotVec(iSiteS+2,jSiteS)*OnsitePotVec(iSiteS+3,jSiteS)-1.0D0)/ stub
              PsiRight=PSI_A(Coord2IndexL(M,iSiteL+1,jSiteL),jState)/ stub
           END IF
           
           !PsiUp
           IF (jSiteL.EQ.1) THEN
              
              IF (IBCFlag.EQ.0) THEN
                 OnsiteUp=ZERO      ! hard wall BC
                 PsiUp=ZERO
                 
              ELSE IF (IBCFlag.EQ.1) THEN
!!$                 CONTINUE           ! periodic BC
                 stub=( OnsitePotVec(iSiteS,4*M-1)*OnsitePotVec(iSiteS,4*M-2)*OnsitePotVec(iSiteS,4*M) &
                      -OnsitePotVec(iSiteS,4*M)-OnsitePotVec(iSiteS,4*M-2))
                 IF( ABS(stub).LT.TINY) THEN           
                    stub= SIGN(TINY,stub)
                 ENDIF
                 OnsiteUp=(OnsitePotVec(iSiteS,4*M-1)*OnsitePotVec(iSiteS,4*M-2)-1.0D0)/stub
                 PsiUp=PSI_A(Coord2IndexL(M,iSiteL,M),jState) /stub      
                 
              ELSE IF (IBCFlag.EQ.2) THEN
                 CONTINUE           ! antiperiodic BC
              ENDIF
           ELSE
              stub=( OnsitePotVec(iSiteS,jSiteS-2)*OnsitePotVec(iSiteS,jSiteS-3)*OnsitePotVec(iSiteS,jSiteS-1) &
                   -OnsitePotVec(iSiteS,jSiteS-1)-OnsitePotVec(iSiteS,jSiteS-3))
              IF( ABS(stub).LT.TINY) THEN           
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteUp=(OnsitePotVec(iSiteS,jSiteS-2)*OnsitePotVec(iSiteS,jSiteS-3)-1.0D0)/stub
              PsiUp=PSI_A(Coord2IndexL(M,iSiteL,jSiteL-1),jState) /stub                   
           END IF
           
           !PsiDown
           IF (jSiteL.EQ.M) THEN
              
              IF (IBCFlag.EQ.0) THEN                                  
                 stub= ( OnsitePotVec(iSiteS,jSiteS+2)*OnsitePotVec(iSiteS,jSiteS+3)*OnsitePotVec(iSiteS,jSiteS+1) &
                      -OnsitePotVec(iSiteS,jSiteS+1)-OnsitePotVec(iSiteS,jSiteS+3))
                 IF( ABS(stub).LT.TINY) THEN           
                    stub= SIGN(TINY,stub)
                 ENDIF
                 OnsiteDown=(OnsitePotVec(iSiteS,jSiteS+2)*OnsitePotVec(iSiteS,jSiteS+3)-1.0D0) /stub    ! hard wall BC
                     
                 PsiDown=ZERO
                 
              ELSE IF (IBCFlag.EQ.1) THEN
!!$                 CONTINUE                                         ! periodic BC
                 stub= ( OnsitePotVec(iSiteS,jSiteS+2)*OnsitePotVec(iSiteS,jSiteS+3)*OnsitePotVec(iSiteS,jSiteS+1) &
                      -OnsitePotVec(iSiteS,jSiteS+1)-OnsitePotVec(iSiteS,jSiteS+3))
                 IF( ABS(stub).LT.TINY) THEN           
                    stub= SIGN(TINY,stub)
                 ENDIF
                 OnsiteDown=(OnsitePotVec(iSiteS,jSiteS+2)*OnsitePotVec(iSiteS,jSiteS+3)-1.0D0 )/stub
                 PsiDown=PSI_A(Coord2IndexL(M,iSiteL,1),jState)/stub
                 
              ELSE IF (IBCFlag.EQ.2) THEN
                 CONTINUE                                         ! antiperiodic BC               
              ENDIF
           ELSE
              stub= ( OnsitePotVec(iSiteS,jSiteS+2)*OnsitePotVec(iSiteS,jSiteS+3)*OnsitePotVec(iSiteS,jSiteS+1) &
                   -OnsitePotVec(iSiteS,jSiteS+1)-OnsitePotVec(iSiteS,jSiteS+3))
              IF( ABS(stub).LT.TINY) THEN           
                 stub= SIGN(TINY,stub)
              ENDIF
              OnsiteDown=(OnsitePotVec(iSiteS,jSiteS+2)*OnsitePotVec(iSiteS,jSiteS+3)-1.0D0 )/stub
              PsiDown=PSI_A(Coord2IndexL(M,iSiteL,jSiteL+1),jState)/stub
           END IF

           new= ( OnsitePot - OnsiteLeft - OnsiteRight - OnsiteUp - OnsiteDown ) * PSI_A(Coord2IndexL(M,iSiteL,jSiteL),jState)&
                - Kappa * ( PsiLeft + PsiRight + PsiUp + PsiDown  ) &
                - PSI_B(Coord2IndexL(M,iSiteL,jSiteL),jState) 
           
           PSI_B(Coord2IndexL(M,iSiteL,jSiteL),jState)= new
           
        ENDDO !jState
        
     ENDDO ! iSite
  ENDDO!jSite
  
  RETURN

END SUBROUTINE TMMultLieb3D_AtoD1

! --------------------------------------------------------------------
! convert i,j coordinates to an index
!!$FUNCTION Coord2IndexL(isize, iSite, jSite)
!!$  INTEGER Coord2IndexL, isize, iSite, jSite
!!$  
!!$  Coord2IndexL= (jSite-1)*isize + iSite
!!$  
!!$  RETURN
!!$END FUNCTION Coord2IndexL
!!$  
! --------------------------------------------------------------------
! TMMultLieb3DBtoA:
!
! 3D version of TMMult2D. Extra boundary conditions

SUBROUTINE TMMultLieb3D_D1toD2(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

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
END SUBROUTINE TMMultLieb3D_D1toD2


SUBROUTINE TMMultLieb3D_D2toD3(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

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

  CALL TMMultLieb3D_D1toD2(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

  RETURN
END SUBROUTINE TMMultLieb3D_D2toD3

SUBROUTINE TMMultLieb3D_D3toA(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

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

  CALL TMMultLieb3D_D1toD2(PSI_A,PSI_B, Ilayer, En, DiagDis, M )

  RETURN
END SUBROUTINE TMMultLieb3D_D3toA
