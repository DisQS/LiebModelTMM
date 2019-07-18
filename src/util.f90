! ********************************************************************
!       
! TMSE2D - Transfer matrix method for the Anderson
! model with diagonal disorder in two dimensions
!
! ********************************************************************
       
! ********************************************************************
!     
! $Header: /home/cvs/phsht/tmseXd/src/util.f90,v 1.2 2012/09/07 10:44:10 phsht Exp $
!
! ********************************************************************
! $Log: util.f90,v $
! Revision 1.2  2012/09/07 10:44:10  phsht
! removed MagFlux and hence all CMPLX,COMPLEX definitions;
! included a final column in AVG output to indicate convergence (or not)
!
! Revision 1.1.1.1  2012/09/07 10:00:09  phsht
! tmseXd
!
! Revision 1.1  2011/05/06 08:13:09  phsht
! 1st installement
!
! Revision 1.8  2010/11/17 16:17:09  phsht
! small typo corrected
!
! Revision 1.7  2010/11/10 10:26:38  phrkaj
! So far working in 1D
!
! Revision 1.6  2010/11/04 15:14:02  phrkaj
! Added 3D stuff. Pre-debug
!
! Revision 1.5  2010/10/26 14:28:07  phrkaj
! Replaced !! with !, fixed the OpenOutputGamma of negative energy in inout.f90
!
! Revision 1.4  2010/10/26 12:47:59  phrkaj
! Got rid of more SB, Conv and Level stuff
!
! Revision 1.3  2010/10/25 15:41:33  phsht
! small changes to remove a "malloc/glibc" error
!
! Revision 1.2  2010/10/25 11:11:01  phrkaj
! Deleted old logs
!
! Revision 1.1.1.1  2010/10/22 12:23:38  phsht
! ParaTMM
!
! ********************************************************************

!!$! --------------------------------------------------------------------
!!$! TMMult2D:
!!$!
!!$! Multiplication of the transfer matrix onto the vector (PSI_A,PSI_B), 
!!$! giving (PSI_B,PSI_A) so that the structure of the transfer matrix 
!!$! can be exploited
!!$
!!$SUBROUTINE TMMult2D(PSI_A,PSI_B, Ilayer, En, DiagDis, M )
!!$
!!$  USE MyNumbers
!!$  USE IPara
!!$  USE RNG
!!$  USE DPara
!!$  
!!$  ! wave functions:
!!$  !       
!!$  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
!!$  
!!$  IMPLICIT NONE
!!$  
!!$  INTEGER Ilayer,           &! current # TM multiplications
!!$       M                     ! strip width
!!$  
!!$  REAL(KIND=RKIND)  DiagDis,&! diagonal disorder
!!$       En                    ! energy
!!$  
!!$  REAL(KIND=CKIND) PSI_A(M,M), PSI_B(M,M)
!!$  
!!$  INTEGER iSite, jState, ISeedDummy
!!$  REAL(KIND=RKIND) OnsitePot
!!$  REAL(KIND=CKIND) new, PsiLeft, PsiRight
!!$  
!!$  !PRINT*,"DBG: TMMult2D()"
!!$  
!!$  DO iSite=1,M
!!$     
!!$     ! create the new onsite potential
!!$     SELECT CASE(IRNGFlag)
!!$     CASE(0)
!!$        OnsitePot= En - DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
!!$     CASE(1)
!!$        OnsitePot= En - DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
!!$     CASE(2)
!!$        OnsitePot= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
!!$     END SELECT
!!$     
!!$     !PRINT*,"iS,pL,RndVec", iSite,pLevel,RndVec((pLevel-1)*M+iSite)
!!$     
!!$     DO jState=1,M
!!$        
!!$        !PRINT*,"jState, iSite", jState, iSite,
!!$        
!!$        IF (iSite.EQ.1) THEN
!!$           
!!$           IF (IBCFlag.EQ.0) THEN
!!$              PsiLeft= CZERO            ! hard wall BC
!!$           ELSE IF (IBCFlag.EQ.1) THEN
!!$              PsiLeft= PSI_A(jState,M)  ! periodic BC
!!$           ELSE IF (IBCFlag.EQ.2) THEN
!!$              PsiLeft= -PSI_A(jState,M) ! antiperiodic BC
!!$           ENDIF
!!$           
!!$        ELSE
!!$           PsiLeft= PSI_A(jState,iSite-1)
!!$        ENDIF
!!$        
!!$        IF (iSite.EQ.M) THEN
!!$           
!!$           IF (IBCFlag.EQ.0) THEN
!!$              PsiRight= CZERO            ! hard wall BC
!!$           ELSE IF (IBCFlag.EQ.1) THEN
!!$              PsiRight= PSI_A(jState,1)  ! periodic BC
!!$           ELSE IF (IBCFlag.EQ.2) THEN
!!$              PsiRight= -PSI_A(jState,1) ! antiperiodic BC
!!$           ENDIF
!!$           
!!$        ELSE
!!$           PsiRight= PSI_A(jState,iSite+1)
!!$        ENDIF
!!$        
!!$        new= ( OnsitePot * PSI_A(jState,iSite) &
!!$             - Kappa * ( PsiLeft + PsiRight ) &
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
!!$  
!!$  RETURN
!!$END SUBROUTINE TMMult2D
!!$
!!$! --------------------------------------------------------------------
!!$! TMMult3D:
!!$!
!!$! 3D version of TMMult2D. Extra boundary conditions
!!$
!!$SUBROUTINE TMMult3D(PSI_A,PSI_B, Ilayer, En, DiagDis, M )
!!$
!!$  USE MyNumbers
!!$  USE IPara
!!$  USE RNG
!!$  USE DPara
!!$  
!!$  ! wave functions:
!!$  !       
!!$  ! (PSI_A, PSI_B) on input, (PSI_B,PSI_A) on output
!!$  
!!$  IMPLICIT NONE
!!$  
!!$  INTEGER Ilayer,           &! current # TM multiplications
!!$       M                     ! strip width
!!$  
!!$  REAL(KIND=RKIND)  DiagDis,&! diagonal disorder
!!$       En                    ! energy
!!$  
!!$  REAL(KIND=RKIND) PSI_A(M*M,M*M), PSI_B(M*M,M*M)
!!$  
!!$  INTEGER iSite, jState, ISeedDummy
!!$  REAL(KIND=RKIND) OnsitePot
!!$  REAL(KIND=RKIND) new, PsiLeft, PsiRight, PsiUp, PsiDown
!!$  
!!$  !PRINT*,"DBG: TMMult3D()"
!!$  
!!$  DO iSite=1,M*M
!!$     
!!$     ! create the new onsite potential
!!$     SELECT CASE(IRNGFlag)
!!$     CASE(0)
!!$        OnsitePot= En - DiagDis*(DRANDOM(ISeedDummy)-0.5D0)
!!$     CASE(1)
!!$        OnsitePot= En - DiagDis*(DRANDOM(ISeedDummy)-0.5D0)*SQRT(12.0D0)
!!$     CASE(2)
!!$        OnsitePot= En - GRANDOM(ISeedDummy,0.0D0,DiagDis)
!!$     END SELECT
!!$     
!!$     !PRINT*,"iS,pL,RndVec", iSite,pLevel,RndVec((pLevel-1)*M+iSite)
!!$     
!!$     DO jState=1,M*M
!!$        
!!$        !PRINT*,"jState, iSite", jState, iSite,
!!$         
!!$        ! PsiLeft
!!$        !IF (iSite.EQ.1) THEN
!!$        IF (MOD(iSite,M).EQ.1) THEN
!!$           
!!$           IF (IBCFlag.EQ.0) THEN
!!$              PsiLeft= CZERO            ! hard wall BC
!!$           ELSE IF (IBCFlag.EQ.1) THEN
!!$              PsiLeft= PSI_A(jState,iSite+M-1)  ! periodic BC
!!$           ELSE IF (IBCFlag.EQ.2) THEN
!!$              PsiLeft= -PSI_A(jState,iSite+M-1) ! antiperiodic BC
!!$           ENDIF
!!$           
!!$        ELSE
!!$           !PRINT*,"DBG: iSite=",iSite
!!$           PsiLeft= PSI_A(jState,iSite-1)
!!$        ENDIF
!!$
!!$        ! PsiRight        
!!$        !IF (iSite.EQ.M) THEN
!!$        IF (MOD(iSite,M).EQ.0) THEN
!!$           
!!$           IF (IBCFlag.EQ.0) THEN
!!$              PsiRight= CZERO            ! hard wall BC
!!$           ELSE IF (IBCFlag.EQ.1) THEN
!!$              PsiRight= PSI_A(jState,iSite-M+1)  ! periodic BC
!!$           ELSE IF (IBCFlag.EQ.2) THEN
!!$              PsiRight= -PSI_A(jState,iSite-M+1) ! antiperiodic BC
!!$           ENDIF
!!$           
!!$        ELSE
!!$           PsiRight= PSI_A(jState,iSite+1)
!!$        ENDIF
!!$
!!$        ! PsiUp
!!$        IF (iSite.GT.(M-1)*M) THEN
!!$           
!!$           IF (IBCFlag.EQ.0) THEN
!!$              PsiUp= CZERO                ! hard wall BC
!!$           ELSE IF (IBCFlag.EQ.1) THEN
!!$              PsiUp= PSI_A(jState,iSite-(M-1)*M)  ! periodic BC
!!$           ELSE IF (IBCFlag.EQ.2) THEN
!!$              PsiUp= -PSI_A(jState,iSite-(M-1)*M) ! antiperiodic BC
!!$           ENDIF
!!$           
!!$        ELSE
!!$           PsiUp= PSI_A(jState,iSite+M)
!!$        ENDIF
!!$
!!$        ! PsiDown        
!!$        IF (iSite.LT.(M+1)) THEN
!!$           
!!$           IF (IBCFlag.EQ.0) THEN
!!$              PsiDown= CZERO            ! hard wall BC
!!$           ELSE IF (IBCFlag.EQ.1) THEN
!!$              PsiDown= PSI_A(jState,iSite+(M-1)*M)  ! periodic BC
!!$           ELSE IF (IBCFlag.EQ.2) THEN
!!$              PsiDown= -PSI_A(jState,iSite+(M-1)*M) ! antiperiodic BC
!!$           ENDIF
!!$           
!!$        ELSE
!!$           PsiDown= PSI_A(jState,iSite-M)
!!$        ENDIF        
!!$        
!!$        new= ( OnsitePot * PSI_A(jState,iSite) &
!!$             - Kappa * ( PsiLeft + PsiRight + PsiUp + PsiDown  ) &
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
!!$  !PRINT*,"PSIA(1,1),(2,1),(1,2),(2,2)",&
!!$  !      PSI_A(1,1),PSI_A(2,1),PSI_A(1,2),PSI_A(2,2)
!!$  
!!$  RETURN
!!$END SUBROUTINE TMMult3D

!	--------------------------------------------------------------------
!	ReNorm:
!
!	Gram-Schmidt orthonormalization of the wave vectors (PSI_A,PSI_B),
!	see e.g. Horn/Johnson, "Matrix Analysis", pp. 15. 
!
!	PLUS reading off the Lyapunov exponents
!
!	(PSI_A,PSI_B) is the incoming vector of eigenvectors,
!	GAMMA, GAMMA2 are Lyapunov exponent and its square,
!	M is the width of the strip,
!
!	IVec/JVec label different vectors, 
!	KIndex is the component index of a SINGLE vector,
!
!	dummy, sum and norm are just different names for the same
!	double precision workspace

SUBROUTINE ReNorm(PSI_A,PSI_B,GAMMA,GAMMA2,M)
  
  USE MyNumbers
  USE IConstants
  USE IPara
  
  INTEGER M
  
  REAL(KIND=RKIND) PSI_A(M,M), PSI_B(M,M)
  REAL(KIND=RKIND) GAMMA(M), GAMMA2(M)
  
  INTEGER IVec,JVec,KIndex
  
  REAL(KIND=RKIND) sum
  REAL(KIND=RKIND) dummy,norm
  EQUIVALENCE (dummy,norm)

  !make the local variables static
  !SAVE
  
  !PRINT*,"DBG: ReNorm()"
  
  DO 100 IVec=1,M
     
     DO 200 JVec=1,IVec-1
        
        sum= ZERO
        
        DO 300 KIndex=1,M
           
!!$           sum= sum + (PSI_A(JVec,KIndex))*PSI_A(IVec,KIndex) &
!!$                + (PSI_B(JVec,KIndex))*PSI_B(IVec,KIndex)
           sum= sum + (PSI_A(KIndex,JVec))*PSI_A(KIndex,IVec) &
                + (PSI_B(KIndex,JVec))*PSI_B(KIndex,IVec)
300     ENDDO
        
        DO 400 KIndex=1,M
           
!!$           PSI_A(IVec,KIndex)= PSI_A(IVec,KIndex) - &
!!$                sum * PSI_A(JVec,KIndex)
!!$           PSI_B(IVec,KIndex)= PSI_B(IVec,KIndex) - &
!!$                sum * PSI_B(JVec,KIndex)
           PSI_A(KIndex,IVec)= PSI_A(KIndex,IVec) - &
                sum * PSI_A(KIndex,JVec)
           PSI_B(KIndex,IVec)= PSI_B(KIndex,IVec) - &
                sum * PSI_B(KIndex,JVec)
           
400     ENDDO
        
200  ENDDO
     
     ! calculation of norm
     
     norm= REAL(0.D0,RKIND)
     DO 500 KIndex=1,M                      
!!$        norm= norm + (PSI_A(IVec,KIndex)) * PSI_A(IVec,KIndex) &
!!$             + (PSI_B(IVec,KIndex)) * PSI_B(IVec,KIndex)
        norm= norm + (PSI_A(KIndex,IVec)) * PSI_A(KIndex,IVec) &
             + (PSI_B(KIndex,IVec)) * PSI_B(KIndex,IVec)
500  ENDDO
     dummy= 1.D0/SQRT(norm)
     DO 600 KIndex=1,M
!!$        PSI_A(IVec,KIndex)= dummy * PSI_A(IVec,KIndex)
!!$        PSI_B(IVec,KIndex)= dummy * PSI_B(IVec,KIndex)
        PSI_A(KIndex,IVec)= dummy * PSI_A(KIndex,IVec)
        PSI_B(KIndex,IVec)= dummy * PSI_B(KIndex,IVec)
600  ENDDO

  !Print*,"PSIA(1,1),(1,2),(1,3),(1,4)",&
        !PSI_A(1,1),PSI_A(1,2),PSI_A(1,3),PSI_A(1,4)
     
     !	----------------------------------------------------------------
     !	gammadummy is ordered s.t. the vector with the smallest overall 
     !	norm is last and the vector with the largest overall norm first. 
     !	We sum this so that the smallest value is used for GAMMA(M)
     !	and the largest for GAMMA(1). Same for GAMMA2(). Therefore, the
     !	inverse of the smallest Lyapunov exponent is obtained by taking
     !	LAMBDA(1) = N/GAMMA(M)
     dummy       = LOG(dummy)
     GAMMA(IVec) = GAMMA(IVec) - dummy
     GAMMA2(IVec)= GAMMA2(IVec) + dummy*dummy
          
     ! ----------------------------------------------------------------
     ! check orthogonality if desired

     IF(IWriteFlag.GE.MAXWriteFlag) THEN

        DO JVec=1,IVec-1
           sum= REAL(0.D0,RKIND)
           DO KIndex=1,M
              sum= sum + (PSI_A(KIndex,JVec))*PSI_A(KIndex,IVec) &
                   + (PSI_B(KIndex,JVec))*PSI_B(KIndex,IVec)
           ENDDO
           PRINT*,"Renorm: <",JVec,"|",IVec,">=",sum
        ENDDO
        
     ENDIF
     
100 ENDDO
  
  RETURN
  
END SUBROUTINE ReNorm

!	--------------------------------------------------------------------
!	Swap:
!
!	(PSI_A,PSI_B)= (old,new) is the incoming vector, this is swapped
!	into (PSI_A,PSI_B)= (new,old)

SUBROUTINE Swap( PSI_A, PSI_B, M)

  USE MyNumbers
  
  INTEGER M
  REAL(KIND=RKIND) PSI_A(M,M), PSI_B(M,M), PSI_dummy(M,M)
  
  INTEGER jState, index
  REAL(KIND=RKIND) dummy
  
  !	PRINT*,"DBG: Swap()"
  
!!$  DO jState=1,M
!!$     DO index=1,M
!!$        
!!$        dummy              = PSI_B(index,jState)
!!$        PSI_B(index,jState)= PSI_A(index,jState)
!!$        PSI_A(index,jState)= dummy
!!$        
!!$     ENDDO
!!$  ENDDO

  PSI_dummy= PSI_A
  PSI_A    = PSI_B
  PSI_B    = PSI_dummy
  
  RETURN
  
END SUBROUTINE Swap


!	--------------------------------------------------------------------
!	ReSort:
!
!	sort the Lyapunov eigenvalues s.t. the largest comes first. RESORT()
!	is ShellSort taken from NumRec, SHELL().

SUBROUTINE ReSort( PSI_A, PSI_B, array0, array1, N )

  USE MyNumbers
  
  INTEGER N
  REAL(KIND=RKIND) PSI_A(N,N),PSI_B(N,N)
  REAL(KIND=RKIND) array0(N), array1(N)
  
  REAL(KIND=RKIND) ALN2I, LocalTINY
  PARAMETER (ALN2I=1.4426950D0, LocalTINY=1.D-5)
  
  INTEGER NN,M,L,K,J,I,LOGNB2, index
  REAL(KIND=RKIND) dummyA, dummyB
  
  !	PRINT*,"DBG: ReSort()"
  
  !	PRINT*,"array0(1),array0(N)",array0(1),array0(N)
  
  LOGNB2=INT(LOG(REAL(N))*ALN2I+LocalTINY)
  M=N
  DO 12 NN=1,LOGNB2
     M=M/2
     K=N-M
     DO 11 J=1,K
        I=J
3       CONTINUE
        L=I+M
        IF(array0(L).GT.array0(I)) THEN
           
           dummyA   = array0(I)
           array0(I)= array0(L)
           array0(L)= dummyA
           
           dummyB   = array1(I)
           array1(I)= array1(L)
           array1(L)= dummyB
           
           DO 100 index=1,N
              dummyA        = PSI_A(index,I)
              dummyB        = PSI_B(index,I)
              
              PSI_A(index,I)= PSI_A(index,L)
              PSI_B(index,I)= PSI_B(index,L)
              
              PSI_A(index,L)= dummyA
              PSI_B(index,L)= dummyB
100        ENDDO
        
           I=I-M
           IF(I.GE.1) GOTO 3
        ENDIF
11   ENDDO
12 ENDDO
  
  !	PRINT*,"array0(1),array0(N)",array0(1),array0(N)
  RETURN

END SUBROUTINE ReSort

