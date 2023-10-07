! --------------------------------------------------------------------
! TMMultLieb2DAtoB:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_B,PSI_A) so that the structure of the transfer matrix 
! can be exploited

SUBROUTINE TMMultLieb2DAtoB(PSI_A,PSI_B, Ilayer, En, DiagDis, M, tz0, tHop, tW )

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
  
  REAL(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  
  INTEGER xSiteL,xSiteS, jState, ISeedDummy
  REAL(KIND=RKIND) OnsitePot, OnsiteRight, OnsiteLeft, OnsitePotVec(2*M)
  REAL(KIND=CKIND) new , PsiLeft, PsiRight, stub
  
  INTEGER, PARAMETER :: LiebSpacer=2
  REAL(KIND=RKIND) PSI_A(M,M),PSI_B(M,M),OnsitePotVec(2*M)
  REAL(KIND=RKIND) txy(LiebSpacer*M,LiebSpacer*M)
  REAL(KIND=RKIND) tz1(M),tz0(M)
  
  INTEGER jState, ISeedDummy, xSiteS,ySiteS, xSiteL,ySiteL, indexK, pos, X0, XL, XR,  numerator, ri, rf
  REAL(KIND=RKIND) OnsitePot, OnsiteRight, OnsiteLeft, OnsiteUp, OnsiteDown
  REAL(KIND=RKIND) new, PsiLeft, PsiRight, PsiUp, PsiDown, stub

!!$  INTEGER, PARAMETER :: LiebSpacer=2

  
  txy=0.0
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
           OnsitePotVec(xSiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
        ELSE
           OnsitePotVec(xSiteS)= -En + 0.0D0
        END IF
     CASE(20)
        IF(Mod(xSiteS,LiebSpacer)==1) THEN
           OnsitePotVec(xSiteS)= -En + DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
        ELSE
           OnsitePotVec(xSiteS)= -En + 0.0D0
        END IF
     CASE(30)
        IF(Mod(xSiteS,LiebSpacer)==1) THEN
           OnsitePotVec(xSiteS)= -En + GRANDOM(ISeedDummy,0.0D0,DiagDis)
        ELSE
           OnsitePotVec(xSiteS)= -En + 0.0D0
        END IF
     END SELECT
  END DO
   ! create the hopping disorder term

  DO xSiteS=1,LiebSpacer*M
     

     X0=  xSiteS

     If(xSiteS==1)Then
        XL= 1
     Else
        XL=  xSiteS -1
     End If

     If(xSiteS==LiebSpacer*M)Then
        XR= 1
     Else
        XR=  xSiteS +1
     End If

     


     txy(X0,XL)= tHop*(1.0 + tW*(DRANDOM(ISeedDummy)-0.5D0))
     txy(XL,X0)= txy(X0,XL)
     txy(X0,XR)= tHop*(1.0 + tW*(DRANDOM(ISeedDummy)-0.5D0))
     txy(XR,X0)= txy(X0,XR)
     
     
  END DO    
  ! to the TMM
  DO xSiteL=1,M
     
     xSiteS= (xSiteL-1)*LiebSpacer + 1
     pos= xSiteL
     tz1(pos)= tz0(pos) ! tz1 stores the past(n-1_n plane) hopping terms
     tz0(pos)= tHop*(1.0 + tW*(DRANDOM(ISeedDummy)-0.5D0)) ! tz0 stores the future(n_n+1 plane) hopping terms
    
     
     OnsitePot= OnsitePotVec(xSiteS) /tz0(pos) 
     
     DO jState=1,M
        
        !PsiLeft
        IF (xSiteL.EQ.1) THEN
           SELECT CASE(IBCFlag)
           CASE(-1,0) ! hard wall BC
              PsiLeft= ZERO            
              OnsiteLeft= ZERO
           CASE(1) ! periodic BC
              stub= OnsitePotVec(LiebSpacer*M)*tz0(pos)  
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)

              Call tap(LiebSpacer*M,1,ri,rf,LiebSpacer,M)
              numerator=txy(ri,rf)
              OnsiteLeft=numerator**2/stub
              Call tap(LiebSpacer*M,LiebSpacer*M-1,ri,rf,LiebSpacer,M)
              numerator=numerator*txy(ri,rf)
              PsiLeft=numerator*PSI_A(M,jState)/stub
           CASE(2) ! antiperiodic BC
              stub= OnsitePotVec(LiebSpacer*M)*tz0(pos)    
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)

              Call tap(LiebSpacer*M,xSiteS,ri,rf,LiebSpacer,M)
              numerator=txy(ri,rf)
              OnsiteLeft=numerator**2/stub                 
              Call tap(LiebSpacer*M,LiebSpacer*M-1,ri,rf,LiebSpacer,M)
              numerator=numerator*txy(ri,rf)

              PsiLeft= -numerator*PSI_A(M,jState) /stub
           CASE DEFAULT
              PRINT*,"TMMultLieb2DAtoB(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
           END SELECT
        ELSE
           stub= OnsitePotVec(xSiteS -1)*tz0(pos)
           IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)

           Call tap(xSiteS-1,xSiteS,ri,rf,LiebSpacer,M)
           numerator=txy(ri,rf)
           OnsiteLeft=numerator**2/stub
           Call tap(xSiteS-2,xSiteS-1,ri,rf,LiebSpacer,M)
           numerator=numerator*txy(ri,rf)              

           PsiLeft= numerator*PSI_A(xSiteL-1,jState)/stub
        END IF

        !PsiRight
        IF (xSiteL.EQ.M) THEN
           SELECT CASE(IBCFlag)
           CASE(-1) ! hard wall BC + STUBS
              stub= OnsitePotVec(xSiteS +1)*tz0(pos)
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
              Call tap(xSiteS+1,xSiteS,ri,rf,LiebSpacer,M)
              OnsiteRight=txy(ri,rf)**2/stub
              PsiRight= ZERO            
           CASE(0) ! hard wall BC
              PsiRight= ZERO            
              OnsiteRight= ZERO
           CASE(1) ! periodic BC
              stub= OnsitePotVec(xSiteS+1)*tz0(pos)*tz0(pos)   
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)

              Call tap(xSiteS+1,xSiteS,ri,rf,LiebSpacer,M)
              numerator=txy(ri,rf)
              OnsiteRight=numerator**2/stub

              Call tap(1,LiebSpacer*M,ri,rf,LiebSpacer,M)
              numerator=numerator*txy(ri,rf)
              
              PsiRight= numerator*PSI_A(1,jState)/stub
           CASE(2) ! antiperiodic BC
              stub= OnsitePotVec(xSiteS +1)*tz0(pos)
              IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)

              Call tap(xSiteS+1,xSiteS,ri,rf,LiebSpacer,M)                 
              numerator=txy(ri,rf)
              OnsiteRight=numerator**2/stub
             
              Call tap(1,LiebSpacer*M,ri,rf,LiebSpacer,M)
              numerator=numerator*txy(ri,rf)
              
              PsiRight= -numerator*PSI_A(1,jState) /stub
           CASE DEFAULT
              PRINT*,"TMMultLieb2DAtoB(): IBCFlag=", IBCFlag, " not implemented --- WRNG!"
           END SELECT
        ELSE
!!$           PsiRight= PSI_A(xSiteL+1,jState)/OnsitePotVec(xSiteS +1)
!!$           OnsiteRight= 1.D0/OnsitePotVec(xSiteS +1)
           stub= OnsitePotVec(xSiteS+1)*tz0(pos)
           IF( ABS(stub).LT.TINY) stub= SIGN(TINY,stub)
           Call tap(xSiteS+1,xSiteS,ri,rf,LiebSpacer,M)
           numerator=txy(ri,rf)
           OnsiteRight=numerator**2/stub

           Call tap(xSiteS+2,xSiteS+1,ri,rf,LiebSpacer,M)
           numerator=numerator*txy(ri,rf)
           
           PsiRight= numerator*PSI_A(xSiteL+1,jState) /stub
        END IF
        
        new =(( OnsitePot - OnsiteLeft - OnsiteRight ) * PSI_A(xSiteL,jState) &
             - Kappa * ( PsiLeft + PsiRight ) &
             -  tz1(pos)*PSI_B(xSiteL,jState)/tz0(pos))
        
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
! Convert from two number labele to one number labele 

Subroutine tap(xSiteSi,xSiteSf,ri,rf,LiebSpacer,M)
  
  Integer xSiteSi,xSiteSf,LiebSpacer,M
  Integer ri,rf    
  
  ri= xSiteSi
  rf= xSiteSf
  
End Subroutine tap
! --------------------------------------------------------------------
! TMMultLieb2DBtoA:
!
! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
! giving (PSI_B,PSI_A) so that the structure of the transfer matrix 
! can be exploited

SUBROUTINE TMMultLieb2DBtoA(PSI_A,PSI_B, Ilayer, En, DiagDis,M, tz0, tHop, tW )

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
       tHop, tW
  REAL(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
  REAL(KIND=RKIND) tz1(M),tz0(M)
  

  INTEGER xSiteL, ySiteL, xSiteS, ySiteS, ISeedDummy,jState, pos,indexK
  INTEGER xSite, jState, ISeedDummy
  REAL(KIND=RKIND) OnsitePot
  REAL(KIND=CKIND) new
  
  !PRINT*,"DBG: TMMultLieb2DBtoA()"
  INTEGER, PARAMETER :: LiebSpacer=2
  DO xSiteL=1,M
     xSiteS= (xSiteL-1)*LiebSpacer + 1
     pos= xSiteL
     tz1(pos)= tz0(pos) ! tz1 stores the past(n-1_n plane) hopping terms
     tz0(pos)= tHop*(1.0 + tW*(DRANDOM(ISeedDummy)-0.5D0)) ! tz0 stores the future(n_n+1 plane) hopping terms
     
     ! create the new onsite potential
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
     !PRINT*,"iS,pL,RndVec", xSite,pLevel,RndVec((pLevel-1)*M+xSite)
     
     DO jState=1,M
        
        !PRINT*,"jState, xSite", jState, xSite,
        
        new= ( OnsitePot * PSI_A(pos,jState) &

             - tz1(pos)*PSI_B(pos,jState)/tz0(pos))
        
        !PRINT*,"i,j,En, OP, PL, PR, PA,PB, PN"
        !PRINT*, xSite, jState, En, OnsitePot, PsiLeft, PsiRight,
        !        PSI_A(xSite,jState), PSI_B(xSite,jState),
        !        new
        
        PSI_B(pos,jState)= new
        
     ENDDO ! jState
  ENDDO ! xSite
  
  !PRINT*,"PSIA(1,1),(1,2),(1,3),(1,4)",&
        !PSI_A(1,1),PSI_A(1,2),PSI_A(1,3),PSI_A(1,4)

  !PRINT*,"PSIB(1,1),(1,2),(1,3),(1,4)",&
        !PSI_B(1,1),PSI_B(1,2),PSI_B(1,3),PSI_B(1,4)
  
  RETURN
END SUBROUTINE TMMultLieb2DBtoA

