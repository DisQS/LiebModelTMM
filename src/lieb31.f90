! --------------------------------------------------------------------
! TMMultLieb3DAtoB:
!
! 3D version of TMMult2D. Extra boundary conditions

SUBROUTINE TMMultLieb3DAtoB(PSI_A,PSI_B, Ilayer, En, DiagDis, M, tz0, tHop, tW )

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
       En,                  &! energy
       tHop,tW

  INTEGER, PARAMETER :: LiebSpacer=2
  
  REAL(KIND=RKIND) PSI_A(M*M,M*M),PSI_B(M*M,M*M),OnsitePotVec(2*M,2*M)
  REAL(KIND=RKIND) txy((LiebSpacer*M)**2,(LiebSpacer*M)**2)
  REAL(KIND=RKIND) tz1(M*M),tz0(M*M)
  
  INTEGER jState, ISeedDummy, xSiteS,ySiteS, xSiteL,ySiteL, indexK,  &
       pos, X0, XL, XR, YU, YD, numerator, ri, rf
  REAL(KIND=RKIND) OnsitePot, OnsiteRight, OnsiteLeft, OnsiteUp, OnsiteDown
  REAL(KIND=RKIND) new, PsiLeft, PsiRight, PsiUp, PsiDown, stub

!!$  INTEGER, PARAMETER :: LiebSpacer=2

  INTEGER C2IL3
  EXTERNAL C2IL3
  
  txy=0.0
  
!!$  Print*,"debug---AtoD----the hopping matrix tz0 in lieb31.f90 passed by main.f90"
!!$
!!$  indexK=0 
!!$  Do pos=1,M**2
!!$
!!$     Print*,tz0(pos)
!!$     indexK=indexK+1
!!$     If(mod(pos,M)==0) Print*,"-------------------------"
!!$  End Do
!!$
!!$  Pause

  !PRINT*,"DBG: TMMultLieb3DAtoB()"
     
  ! create the new onsite potential
!!$  DO xSiteS=1,LiebSpacer*M
!!$     DO ySiteS=1,LiebSpacer*M
!!$        SELECT CASE(IRNGFlag)
!!$        CASE(0)
!!$           OnsitePotVec(xSiteS,ySiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
!!$        CASE(1)
!!$           OnsitePotVec(xSiteS,ySiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
!!$        CASE(2)
!!$           OnsitePotVec(xSiteS,ySiteS)= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
!!$        END SELECT
!!$     END DO
!!$  END DO

  ! create the new onsite potential
  !  IRNGFlag=(xy)
  !  xy=0x      presenting conditions that all positions are disorder
  !  xy=x0,     presenting conditions that only central positions are disorder 
   
  DO xSiteS=1,LiebSpacer*M
     Do ySiteS=1,LiebSpacer*M
        SELECT CASE(IRNGFlag)
        CASE(01)
           OnsitePotVec(xSiteS,ySiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
        CASE(02)
           OnsitePotVec(xSiteS,ySiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
        CASE(03)
           OnsitePotVec(xSiteS,ySiteS)= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
        CASE(10)
           IF(Mod(xSiteS,LiebSpacer)==1 .AND. Mod(ySiteS,LiebSpacer)==1) THEN
              OnsitePotVec(xSiteS,ySiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
           ELSE
              OnsitePotVec(xSiteS,ySiteS)= -En + 0.0D0
           END IF
        CASE(20)
           IF(Mod(xSiteS,LiebSpacer)==1 .AND. Mod(ySiteS,LiebSpacer)==1) THEN
              OnsitePotVec(xSiteS,ySiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
           ELSE
              OnsitePotVec(xSiteS,ySiteS)= -En + 0.0D0
           END IF
        CASE(30)
           IF(Mod(xSiteS,LiebSpacer)==1 .AND. Mod(ySiteS,LiebSpacer)==1) THEN
              OnsitePotVec(xSiteS,ySiteS)= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
           ELSE
              OnsitePotVec(xSiteS,ySiteS)= -En + 0.0D0
           END IF
        END SELECT
     END DO
  END DO


  ! create the hopping disorder term

  DO xSiteS=1,LiebSpacer*M
     Do ySiteS=1,LiebSpacer*M

        X0= (ySiteS-1)*LiebSpacer*M + xSiteS

        If(xSiteS==1)Then
           XL= ySiteS*LiebSpacer*M
        Else
           XL= (ySiteS-1)*LiebSpacer*M + xSiteS -1
        End If

        If(xSiteS==LiebSpacer*M)Then
           XR= (ySiteS-1)*LiebSpacer*M +1
        Else
           XR= (ySiteS-1)*LiebSpacer*M + xSiteS +1
        End If

        If(ySiteS==1)Then
           YU= (LiebSpacer*M-1)*LiebSpacer*M + xSiteS
        Else
           YU= (ySiteS-2)*LiebSpacer*M + xSiteS
        End If

        If(ySiteS==LiebSpacer*M)Then
           YD= xSiteS
        Else
           YD= ySiteS*LiebSpacer*M + xSiteS
        End If

        txy(X0,XL)= tHop*(1.0 + tW*(DRANDOM(ISeedDummy)-0.5D0))
        txy(XL,X0)= txy(X0,XL)
        txy(X0,XR)= tHop*(1.0 + tW*(DRANDOM(ISeedDummy)-0.5D0))
        txy(XR,X0)= txy(X0,XR)

        txy(X0,YU)= tHop*(1.0 + tW*(DRANDOM(ISeedDummy)-0.5D0))
        txy(YU,X0)= txy(X0,YU)
        txy(X0,YD)= tHop*(1.0 + tW*(DRANDOM(ISeedDummy)-0.5D0))
        txy(YD,X0)= txy(X0,YD)

     END DO
  END DO

!!$  Print*,"debug----------- in lieb31.f90-----------"
!!$  Print*,"tHop=",tHop
!!$
!!$  Pause

  
!!$  Print*,"debug-------the hopping matrix txy generated in lieb31.f90-----------"
!!$
!!$  Do xSiteS=1,LiebSpacer*M
!!$     Do ySiteS=1,LiebSpacer*M
!!$
!!$        Print*,txy(xSiteS,ySiteS)
!!$        
!!$     End Do
!!$     Print*,"--------------------------"
!!$  End Do
!!$
!!$  Pause

  
  !to the TMM
  DO xSiteL=1,M
     DO ySiteL=1,M
        
        xSiteS= (xSiteL-1)*LiebSpacer + 1
        ySiteS= (ySiteL-1)*LiebSpacer + 1

!!$        pos=(ySiteS-1)*LiebSpacer*M + xSiteS
        pos=(ySiteL-1)*M + xSiteL
        
        tz1(pos)= tz0(pos) ! tz1 stores the past(n-1_n plane) hopping terms
        tz0(pos)= tHop*(1.0 + tW*(DRANDOM(ISeedDummy)-0.5D0)) ! tz0 stores the future(n_n+1 plane) hopping terms

        OnsitePot=OnsitePotVec(xSiteS,ySiteS)/tz0(pos)

        DO jState=1,M*M
           
           !PsiLeft
           IF (xSiteL.LE.1) THEN !mod(indexK,M).EQ.1)
              SELECT CASE(IBCFlag)
              CASE(-1,0) ! hard wall BC
                 OnsiteLeft=ZERO      
                 PsiLeft=ZERO
              CASE(1) ! periodic BC
                 stub= OnsitePotVec(LiebSpacer*M,ySiteS)*tz0(pos)    
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 Call tmap(LiebSpacer*M,1,ySiteS,ySiteS,ri,rf,LiebSpacer,M)
                 numerator=txy(ri,rf)
                 OnsiteLeft=numerator**2/stub
                 Call tmap(LiebSpacer*M,LiebSpacer*M-1,ySiteS,ySiteS,ri,rf,LiebSpacer,M)
                 numerator=numerator*txy(ri,rf)
                 
                 PsiLeft=numerator*PSI_A(C2IL3(M,M,ySiteL),jState)/stub
              CASE(2) ! antiperiodic BC
                 stub= OnsitePotVec(LiebSpacer*M,ySiteS)*tz0(pos)    
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 Call tmap(LiebSpacer*M,xSiteS,ySiteS,ySiteS,ri,rf,LiebSpacer,M)
                 numerator=txy(ri,rf)
                 OnsiteLeft=numerator**2/stub                 
                 Call tmap(LiebSpacer*M,LiebSpacer*M-1,ySiteS,ySiteS,ri,rf,LiebSpacer,M)
                 numerator=numerator*txy(ri,rf)

                 PsiLeft= -numerator*PSI_A(C2IL3(M,M,ySiteL),jState) /stub
              CASE DEFAULT
                 PRINT*,"TMMultLieb3DAtoB(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
!!$              stub= OnsitePotVec(xSiteS-1,ySiteS)
              
              stub= OnsitePotVec(xSiteS-1,ySiteS)*tz0(pos)
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
!!$              OnsiteLeft=1.0D0/stub
              Call tmap(xSiteS-1,xSiteS,ySiteS,ySiteS,ri,rf,LiebSpacer,M)
              numerator=txy(ri,rf)
              OnsiteLeft=numerator**2/stub
              Call tmap(xSiteS-2,xSiteS-1,ySiteS,ySiteS,ri,rf,LiebSpacer,M)
              numerator=numerator*txy(ri,rf)              

              PsiLeft= numerator*PSI_A(C2IL3(M,xSiteL-1,ySiteL),jState)/stub
           END IF
           
           !PsiRight
           IF (xSiteL.GE.M) THEN !(mod(indexK,M).EQ.0) 
              SELECT CASE(IBCFLag)
              CASE(-1) ! hard wall BC with stubs
                 stub= OnsitePotVec(xSiteS+1,ySiteS)*tz0(pos)  
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 Call tmap(xSiteS+1,xSiteS,ySiteS,ySiteS,ri,rf,LiebSpacer,M)
                 OnsiteRight=txy(ri,rf)**2/stub
               
                 PsiRight= ZERO
              CASE(0) ! hard wall BC
                 OnsiteRight= ZERO
                 PsiRight= ZERO
              CASE(1) ! periodic BC
                 stub= OnsitePotVec(xSiteS+1,ySiteS)*tz0(pos)   
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 Call tmap(xSiteS+1,xSiteS,ySiteS,ySiteS,ri,rf,LiebSpacer,M)
                 numerator=txy(ri,rf)
                 OnsiteRight=numerator**2/stub

                 Call tmap(1,LiebSpacer*M,ySiteS,ySiteS,ri,rf,LiebSpacer,M)
                 numerator=numerator*txy(ri,rf)
                 
                 PsiRight= numerator*PSI_A(C2IL3(M,1,ySiteL),jState)/stub
              CASE(2) ! antiperiodic BC
                 stub= OnsitePotVec(xSiteS+1,ySiteS)*tz0(pos)  
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 Call tmap(xSiteS+1,xSiteS,ySiteS,ySiteS,ri,rf,LiebSpacer,M)                 
                 numerator=txy(ri,rf)
                 OnsiteRight=numerator**2/stub
                
                 Call tmap(1,LiebSpacer*M,ySiteS,ySiteS,ri,rf,LiebSpacer,M)
                 numerator=numerator*txy(ri,rf)
                 
                 PsiRight= -numerator*PSI_A(C2IL3(M,1,ySiteL),jState) /stub
              CASE DEFAULT
                 PRINT*,"TMMultLieb3DAtoB(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              
              stub= OnsitePotVec(xSiteS+1,ySiteS)*tz0(pos)
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              Call tmap(xSiteS+1,xSiteS,ySiteS,ySiteS,ri,rf,LiebSpacer,M)
              numerator=txy(ri,rf)
              OnsiteRight=numerator**2/stub

              Call tmap(xSiteS+2,xSiteS+1,ySiteS,ySiteS,ri,rf,LiebSpacer,M)
              numerator=numerator*txy(ri,rf)
              
              PsiRight= numerator*PSI_A(C2IL3(M,xSiteL+1,ySiteL),jState) /stub
           END IF

           !PsiDown
           IF (ySiteL.LE.1) THEN !(indexK<=M)
              SELECT CASE(IBCFlag)
              CASE(-1,0) ! hard wall BC
                 PsiDown= ZERO            
                 OnsiteDown= ZERO
              CASE(1) ! periodic BC
                 stub= OnsitePotVec(xSiteS,LiebSpacer*M)*tz0(pos)
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 Call tmap(xSiteS,xSiteS,LiebSpacer*M,1,ri,rf,LiebSpacer,M)
                 numerator=txy(ri,rf)
                 OnsiteLeft=numerator**2/stub
                 
                 Call tmap(xSiteS,xSiteS,LiebSpacer*M-1,LiebSpacer*M,ri,rf,LiebSpacer,M)
                 numerator=numerator*txy(ri,rf)
                 
                 PsiDown= numerator*PSI_A(C2IL3(M,xSiteL,M),jState) /stub
              CASE(2) ! antiperiodic BC
                 stub= OnsitePotVec(xSiteS,LiebSpacer*M)*tz0(pos)
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 Call tmap(xSiteS,xSiteS,LiebSpacer*M,1,ri,rf,LiebSpacer,M)
                 numerator=txy(ri,rf)
                 OnsiteDown= numerator**2/stub
                 
                 Call tmap(xSiteS,xSiteS,LiebSpacer*M-1,LiebSpacer*M,ri,rf,LiebSpacer,M)
                 numerator=numerator*txy(ri,rf)
                 
                 PsiDown= -numerator*PSI_A(C2IL3(M,xSiteL,M),jState) /stub
              CASE DEFAULT
                 PRINT*,"TMMultLieb3DAtoB(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              
              stub= OnsitePotVec(xSiteS,ySiteS-1)*tz0(pos)
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              Call tmap(xSiteS,xSiteS,ySiteS-1,ySiteS,ri,rf,LiebSpacer,M)
              numerator=txy(ri,rf)
              OnsiteDown= numerator**2/stub

              Call tmap(xSiteS,xSiteS,ySiteS-2,ySiteS-1,ri,rf,LiebSpacer,M)
              numerator=numerator*txy(ri,rf)
              
              PsiDown= numerator*PSI_A(C2IL3(M,xSiteL,ySiteL-1),jState)/stub
           END IF
           
           !PsiUp
           IF (ySiteL.GE.M) THEN !(indexK>M*(M-1))
              SELECT CASE(IBCFlag)
              CASE(-1) ! hard wall BC with stubs
                 stub= OnsitePotVec(xSiteS,ySiteS+1)*tz0(pos)    
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 Call tmap(xSiteS,xSiteS,ySiteS+1,ySiteS,ri,rf,LiebSpacer,M)
                 OnsiteUp=txy(ri,rf)**2/stub
                 
                 PsiUp= ZERO
              CASE(0) ! hard wall BC
                 OnsiteUp= ZERO
                 PsiUp= ZERO
              CASE(1) ! periodic BC
                 stub= OnsitePotVec(xSiteS,ySiteS+1)*tz0(pos)    
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 Call tmap(xSiteS,xSiteS,ySiteS+1,ySiteS,ri,rf,LiebSpacer,M)
                 numerator=txy(ri,rf)
                 OnsiteUp=numerator**2/stub                 
                 
                 Call tmap(xSiteS,xSiteS,1,LiebSpacer*M,ri,rf,LiebSpacer,M)
                 numerator=numerator*txy(ri,rf)
                 !   PsiUp= 1.0D0/OnsitePotVec(xSiteS,ySiteS+1)*PSI_A(jState,mod(C2IL3(M,xSiteL,ySiteL),M))
                 PsiUp=numerator*PSI_A(C2IL3(M,xSiteL,1),jState) /stub
              CASE(2) ! antiperiodic BC
                 stub= OnsitePotVec(xSiteS,ySiteS+1)*tz0(pos)    
                 IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
                 Call tmap(xSiteS,xSiteS,ySiteS+1,ySiteS,ri,rf,LiebSpacer,M)
                 numerator=txy(ri,rf)
                 OnsiteUp=numerator**2/stub                 
              
                 Call tmap(xSiteS,xSiteS,1,ySiteS+1,ri,rf,LiebSpacer,M)
                 numerator=numerator*txy(ri,rf)

                 !   PsiUp= 1.0D0/OnsitePotVec(xSiteS,ySiteS+1)*PSI_A(jState,mod(indexK,M))
                 PsiUp=-numerator*PSI_A(C2IL3(M,xSiteL,1),jState) /stub
              CASE DEFAULT
                 PRINT*,"TMMultLieb3DAtoB(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
              END SELECT
           ELSE
              
              stub= OnsitePotVec(xSiteS,ySiteS+1)*tz0(pos)
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              Call tmap(xSiteS,xSiteS,ySiteS+1,ySiteS,ri,rf,LiebSpacer,M)
              numerator=txy(ri,rf)
              OnsiteUp=numerator**2/stub

              Call tmap(xSiteS,xSiteS,ySiteS+2,ySiteS+1,ri,rf,LiebSpacer,M)
              numerator=numerator*txy(ri,rf)

              PsiUp= numerator*PSI_A(C2IL3(M,xSiteL,ySiteL+1),jState) /stub
           END IF
           
           !PRINT*,"DBG2: jState,xSiteS, ySiteS, indexK", jState, xSiteS, ySiteS, indexK
           new= ( OnsitePot - OnsiteLeft - OnsiteRight - OnsiteDown - OnsiteUp ) * PSI_A(C2IL3(M,xSiteL,ySiteL),jState)&
                - ( PsiLeft + PsiRight + PsiDown + PsiUp ) &
                - tz1(pos)*PSI_B(C2IL3(M,xSiteL,ySiteL),jState)/tz0(pos) 
           
           !PRINT*,"xSiteS,ySiteS,En, OP, PL, PR, PA,PB, PN"
           !PRINT*, xSiteS, jState, En, OnsitePot, PsiLeft, PsiRight,
           !        PSI_A(xSiteL,jState), PSI_B(xSite,jState),
           !        new
           
           PSI_B(C2IL3(M,xSiteL,ySiteL),jState)= new
           
        ENDDO !jState
        
     ENDDO ! xSiteL
  ENDDO !ySiteL


!!$  Print*,"debug---AtoD----the hopping matrix tz0 generated by lieb31.f90"
!!$
!!$  indexK=0 
!!$  Do pos=1,M**2
!!$
!!$     Print*,tz0(pos)
!!$     indexK=indexK+1
!!$     If(mod(pos,M)==0) Print*,"-------------------------"
!!$  End Do
!!$
!!$  Pause
  
  RETURN

END SUBROUTINE TMMultLieb3DAtoB

! --------------------------------------------------------------------
! convert i,j coordinates to an index
! used in lieb31/32/33
FUNCTION C2IL3(M, xSiteL, ySiteL)
  INTEGER C2IL3, M, xSiteL, ySiteL

  !old: indexK=(jSite/2)*M+(iSite+1)/2
  C2IL3= (ySiteL-1)*M + xSiteL
  
  RETURN
END FUNCTION C2IL3


!------------------------------------------------------------------
! Convert from two number labele to one number labele 

Subroutine tmap(xSiteSi,xSiteSf,ySiteSi,ySiteSf,ri,rf,LiebSpacer,M)
  
  Integer xSiteSi,xSiteSf,ySiteSi,ySiteSf,LiebSpacer,M
  Integer ri,rf    
  
  ri= (ySiteSi-1)*LiebSpacer*M + xSiteSi
  rf= (ySiteSf-1)*LiebSpacer*M + xSiteSf
  
End Subroutine tmap


! --------------------------------------------------------------------
! TMMultLieb3DBtoA:
!
! 3D version of TMMult2D. Extra boundary conditions

SUBROUTINE TMMultLieb3DBtoA(PSI_A,PSI_B, Ilayer, En, DiagDis, M ,tz0, tHop, tW)

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
       En,                  &! energy
       tHop, tW
  
  REAL(KIND=CKIND) PSI_A(M*M,M*M), PSI_B(M*M,M*M)
  REAL(KIND=RKIND) tz1(M*M),tz0(M*M)
  
  INTEGER xSiteL, ySiteL, xSiteS, ySiteS, ISeedDummy,jState, pos,indexK
  REAL(KIND=RKIND) OnsitePot
  REAL(KIND=CKIND) new

  INTEGER, PARAMETER :: LiebSpacer=2

  
  !PRINT*,"DBG: TMMultLieb3DBtoA()"
!!$
!!$  DO xSiteS=1,LiebSpacer*M
!!$     DO ySiteS=1,LiebSpacer*M
!!$
!!$        X0= (ySiteS-1)*LiebSpacer*M + xSiteS
!!$        tz0(X0)=tHop*(1.0 + tW*(DRANDOM(ISeedDummy)-0.5D0))
!!$
!!$     END DO
!!$  END DO

!!$  Print*,"debug---DtoA----the hopping matrix tz0 in lieb31.f90 passed by AtoD"
!!$
!!$  indexK=0 
!!$  Do pos=1,M**2
!!$
!!$     Print*,tz0(pos)
!!$     indexK=indexK+1
!!$     If(mod(pos,M)==0) Print*,"-------------------------"
!!$  End Do
!!$
!!$  Pause
  
  
  DO xSiteL=1,M
     DO ySiteL=1,M

        xSiteS= (xSiteL-1)*LiebSpacer + 1
        ySiteS= (ySiteL-1)*LiebSpacer + 1

!!$        pos=(ySiteS-1)*LiebSpacer*M + xSiteS
        pos=(ySiteL-1)*M + xSiteL

        tz1(pos)= tz0(pos) ! tz1 stores the past(n-1_n plane) hopping terms
        tz0(pos)=tHop*(1.0 + tW*(DRANDOM(ISeedDummy)-0.5D0))
     
     ! create the new onsite potential
!!$     SELECT CASE(IRNGFlag)
!!$     CASE(0)
!!$        OnsitePot= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
!!$     CASE(1)
!!$        OnsitePot= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
!!$     CASE(2)
!!$        OnsitePot= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
!!$     END SELECT

        SELECT CASE(IRNGFlag)
        CASE(01)
           OnsitePot= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
        CASE(02)
           OnsitePot= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
        CASE(03)
           OnsitePot= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
        CASE(10)
           OnsitePot= -En + 0.0D0
        CASE(20)
           OnsitePot= -En + 0.0D0
        CASE(30)
           OnsitePot= -En + 0.0D0
        END SELECT

        OnsitePot=OnsitePot/tz0(pos)
     
     !PRINT*,"iS,pL,RndVec", iSite,pLevel,RndVec((pLevel-1)*M+iSite)
     
        DO jState=1,M*M

           !PRINT*,"jState, iSite", jState, iSite,

           new= ( OnsitePot * PSI_A(pos,jState) &
                - tz1(pos)*PSI_B(pos,jState)/tz0(pos) )

           !PRINT*,"i,jSite,En, OP, PL, PR, PA,PB, PN"
           !PRINT*, iSite, jState, En, OnsitePot, PsiLeft, PsiRight,
           !        PSI_A(iSite,jState), PSI_B(iSite,jState),
           !        new

           PSI_B(pos,jState)= new

        ENDDO ! jState
     ENDDO ! jSiteL
  ENDDO !iSiteL

!!$  Print*,"debug---DtoA----the hopping matrix tz0 generated in lieb31.f90"
!!$
!!$  indexK=0 
!!$  Do pos=1,M**2
!!$
!!$     Print*,tz0(pos)
!!$     indexK=indexK+1
!!$     If(mod(pos,M)==0) Print*,"-------------------------"
!!$  End Do
!!$
!!$  Pause
  
  RETURN
END SUBROUTINE TMMultLieb3DBtoA

