MODULE TRICUBSPLINEUTILS
  ! WARNING: this whole module is a misnomer. All it does is local tricubic interpolation (*not* global splines!)
  ! but still gives a smooth interpolation with smooth derivatives
  ! utilities for tricubic splines in 3D

  IMPLICIT NONE
  DOUBLE PRECISION :: BINVMAT(64,64)  
  LOGICAL :: SPLINETEST  
  TYPE SPDATA3D
     ! this is a spline data object
     ! storing a 3D table of data and the corresponding spline coefficients
     ! spline coefficients and data matrix
     DOUBLE PRECISION, POINTER :: SPCOEFF(:,:,:,:), DATAMAT(:,:,:)
     ! axis coordinates
     DOUBLE PRECISION, POINTER :: X1(:),X2(:),X3(:)
     ! number of coords in each dimension
     INTEGER :: N1, N2, N3
     ! have arrays been allocated?
     LOGICAL :: ARRAYSET = .FALSE.
  END TYPE SPDATA3D  

CONTAINS
  SUBROUTINE EVALUATESPLINE(SP,X1,X2,X3,F,DF1,DF2,DF3)
    USE GENUTILS, ONLY : INTERP1
    ! using cubic spline interpolation, get the function  and its derivatives
    ! SP is the spline data object to use
    IMPLICIT NONE
    TYPE(SPDATA3D), POINTER :: SP
    DOUBLE PRECISION ,INTENT(IN) :: X1,X2,X3
    DOUBLE PRECISION, INTENT(OUT) :: F
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: DF1,DF2,DF3
    INTEGER :: I, J, K, IND1, IND2, IND3
    DOUBLE PRECISION :: COEFF, Z1, Z2, Z3,DX,DY,DZ

    IF (.NOT.SP%ARRAYSET) THEN
       PRINT*, 'ERROR IN EVALUATESPLINE: spline coefficients not yet set up!'
       stop 1
    endif

    ! get appropriate array indices
    CALL INTERP1(SP%X1,SP%N1,X1,IND1)
    CALL INTERP1(SP%X2,SP%N2,X2,IND2)
    CALL INTERP1(SP%X3,SP%N3,X3,IND3)

    IF (IND1.LE.0.OR.IND1.GE.SP%N1.OR.IND2.LE.0.OR.IND2.GE.SP%N2 &
         & .OR.IND3.LE.0.OR.IND3.GE.SP%N3) THEN
       PRINT*, 'ERROR IN EVALUATESPLINE: value out of bounds', IND1, IND2, IND3
       PRINT*, X1,SP%X1(1),SP%X1(SP%N1), X1.EQ.SP%X1(1)
       PRINT*, X2,SP%X2(1),SP%X2(SP%N2),  X2.EQ.SP%X2(1)
       PRINT*, X3,SP%X3(1),SP%X3(SP%N3),  X3.EQ.SP%X3(1)
       STOP 1
       IND1 = MIN(MAX(1,IND1),SP%N1)
       IND2 = MIN(MAX(1,IND2),SP%N2)
       IND3 = MIN(MAX(1,IND3),SP%N3)
    ENDIF

    DX = SP%X1(IND1+1)-SP%X1(IND1)
    DY = SP%X2(IND2+1)-SP%X2(IND2)
    DZ = SP%X3(IND3+1)-SP%X3(IND3)

    Z1 = (X1 - SP%X1(IND1))/DX
    Z2 = (X2 - SP%X2(IND2))/DY
    Z3 = (X3 - SP%X3(IND3))/DZ

    F = 0D0
    IF (PRESENT(DF1)) DF1 = 0D0
    IF (PRESENT(DF2)) DF2 = 0D0
    IF (PRESENT(DF3)) DF3 = 0D0

    DO I = 0,3
       DO J = 0,3
          DO K = 0,3
             COEFF = SP%SPCOEFF(1+I+4*J+16*K,IND1,IND2,IND3)
             
             F = F + COEFF*(Z1**I)*(Z2**J)*(Z3**K)

             IF (PRESENT(DF1).AND.I.GT.0) THEN
                DF1 = DF1 + I*COEFF*(Z1**(I-1))*(Z2**J)*(Z3**K)
             ENDIF
             IF (PRESENT(DF2).AND.J.GT.0) THEN
                DF2 = DF2 + J*COEFF*(Z1**I)*(Z2**(J-1))*(Z3**K)
             ENDIF
             IF (PRESENT(DF3).AND.K.GT.0) THEN
                DF3 = DF3 + K*COEFF*(Z1**I)*(Z2**J)*(Z3**(K-1))
             ENDIF

          ENDDO
       ENDDO
    ENDDO   

    IF (PRESENT(DF1)) DF1 = DF1/DX
    IF (PRESENT(DF2)) DF2 = DF2/DY
    IF (PRESENT(DF3)) DF3 = DF3/DZ

    IF ((.NOT.F.GE.0.AND..NOT.F.LE.0).OR.SPLINETEST) THEN
       PRINT*, 'BAD SPLINE VALUE:', F
       PRINT*, IND1, IND2, IND3
       PRINT*, SP%X1(IND1), SP%X2(IND2), SP%X3(IND3)
       PRINT*, 'SPCOEFF:'
       PRINT*, SP%SPCOEFF(:,IND1,IND2,IND3)
       PRINT*, 'MAT VALUES:'
       PRINT*, SP%DATAMAT(IND1,IND2,IND3)
       PRINT*, SP%DATAMAT(IND1+1,IND2,IND3)
       PRINT*, SP%DATAMAT(IND1,IND2+1,IND3)
       PRINT*, SP%DATAMAT(IND1,IND2,IND3+1)
       PRINT*, SP%DATAMAT(IND1+1,IND2+1,IND3)
       PRINT*, SP%DATAMAT(IND1+1,IND2,IND3+1)
       PRINT*, SP%DATAMAT(IND1,IND2+1,IND3+1)
       PRINT*, SP%DATAMAT(IND1+1,IND2+1,IND3+1)
       STOP 1
    ENDIF
  END SUBROUTINE EVALUATESPLINE

  SUBROUTINE CREATESPLINE(SP)
    ! get all the spline coefficients for interpolating the distance matrix
    ! N1,N2,N3 are the number of points in the cos(theta),alpha,beta dimensions
    ! X1,X2,X3 are the grid points in the three dimensions
    ! DISTMAT is the matrix of minimal distances
    IMPLICIT NONE
    TYPE(SPDATA3D), POINTER :: SP
    INTEGER :: I, J, K, IM, IP, JM, JP, KM, KP, C,N1,N2,N3
    DOUBLE PRECISION :: DI, DJ, DK
    DOUBLE PRECISION, DIMENSION(SP%N1,SP%N2,SP%N3) :: DX, DY, DZ, DXY, DXZ, DYZ, DXYZ
    INTEGER :: ILIST(8), JLIST(8), KLIST(8)
    DOUBLE PRECISION :: BVEC(64), COEFF(64), TEST,KRANGE

    N1 = SP%N1; N2 = SP%N2; N3 = SP%N3

    IF (MINVAL((/N1,N2,N3/)).LT.2) THEN
       PRINT*, 'ERROR IN SETUPSLINE: need at least 2 points in each dimension'
       STOP 1
    ENDIF

    IF (.NOT.SP%ARRAYSET) THEN
       PRINT*, 'ARRAYS NOT YET ALLOCATED. CANNOT CREATE SPLINE.'
       STOP 1
    ENDIF

    ! first get all the partial derivative estimates at each grid point
    ! using centered differencing
    ! 0 second derivatives at the edges in each individual dimension    
    PRINT*, 'GETTING DERIVATIVES...'
    KRANGE = SP%X3(N3)-SP%X3(1)
    DO I = 1,N1
       IM = MAX(I-1,1)
       IP = MIN(I+1,N1)
       DI = SP%X1(IP)-SP%X1(IM)

       DO J = 1,N2
          JM = MAX(J-1,1)
          JP = MIN(J+1,N2)
          DJ = SP%X2(JP)-SP%X2(JM)

          DO K = 1,N3
             ! want periodic boundaries in phi
             IF (K.EQ.1) THEN
                KM = N3-1; KP = K+1
                DK = SP%X3(KP)-SP%X3(KM)+KRANGE
             ELSEIF(K.EQ.N3) THEN
                KP=2; KM = K-1
                DK = SP%X3(KP)-SP%X3(KM)+KRANGE              
             ELSE
                KM = K-1; KP = K+1
                DK = SP%X3(KP)-SP%X3(KM)
             ENDIF
                          
             
             ! first derivatives
             DX(I,J,K) = (SP%DATAMAT(IP,J,K)-SP%DATAMAT(IM,J,K))/DI             
             DY(I,J,K) = (SP%DATAMAT(I,JP,K)-SP%DATAMAT(I,JM,K))/DJ
             DZ(I,J,K) = (SP%DATAMAT(I,J,KP)-SP%DATAMAT(I,J,KM))/DK             

             ! second derivatives
             DXY(I,J,K) = (SP%DATAMAT(IP,JP,K)-SP%DATAMAT(IP,JM,K)&
                  & -SP%DATAMAT(IM,JP,K) + SP%DATAMAT(IM,JM,K))/DI/DJ
             DXZ(I,J,K) = (SP%DATAMAT(IP,J,KP)-SP%DATAMAT(IP,J,KM)&
                  & -SP%DATAMAT(IM,J,KP) + SP%DATAMAT(IM,J,KM))/DI/DK
             DYZ(I,J,K) = (SP%DATAMAT(I,JP,KP)-SP%DATAMAT(I,JP,KM)&
                  & -SP%DATAMAT(I,JM,KP) + SP%DATAMAT(I,JM,KM))/DJ/DK

             ! third derivative
             DXYZ(I,J,K) = (SP%DATAMAT(IP,JP,KP) - SP%DATAMAT(IP,JP,KM) - SP%DATAMAT(IP,JM,KP) & 
                  & - SP%DATAMAT(IM,JP,KP) + SP%DATAMAT(IP,JM,KM) + SP%DATAMAT(IM,JP,KM) &
                  & + SP%DATAMAT(IM,JM,KP) - SP%DATAMAT(IM,JM,KM))/DI/DJ/DK
          ENDDO
       ENDDO
    ENDDO

    PRINT*, 'GETTING SPLINE COEFFICIENTS...'
    ! for each spline cell, get the coefficients 
    DO I = 1,N1-1       
       ILIST= (/I,I+1,I,I+1,I,I+1,I,I+1/) 
       DI = SP%X1(I+1)-SP%X1(I)
       DO J = 1,N2-1
          JLIST = (/J,J,J+1,J+1,J,J,J+1,J+1/)
          DJ = SP%X2(J+1)-SP%X2(J)
          DO K = 1,N3-1
             KLIST = (/K,K,K,K,K+1,K+1,K+1,K+1/)
             DK = SP%X3(K+1)-SP%X3(K)

             ! set up the vector of 64 derivatives
             ! NOTE: these are scaled by box size since algorithm works with unit cubes
             DO C = 1,8
                BVEC(C) = SP%DATAMAT(ILIST(C),JLIST(C),KLIST(C))
                BVEC(8+C) = DX(ILIST(C),JLIST(C),KLIST(C))*DI
                BVEC(16+C) = DY(ILIST(C),JLIST(C),KLIST(C))*DJ
                BVEC(24+C) = DZ(ILIST(C),JLIST(C),KLIST(C))*DK
                BVEC(32+C) = DXY(ILIST(C),JLIST(C),KLIST(C))*DI*DJ
                BVEC(40+C) = DXZ(ILIST(C),JLIST(C),KLIST(C))*DI*DK
                BVEC(48+C) = DYZ(ILIST(C),JLIST(C),KLIST(C))*DJ*DK
                BVEC(56+C) = DXYZ(ILIST(C),JLIST(C),KLIST(C))*DI*DJ*DK
             ENDDO
             TEST = SUM(BVEC)
             IF (.NOT.TEST.LE.0.AND..NOT.TEST.GE.0) THEN
                PRINT*, 'BAD B VECTOR IN CUBIC INTERPOLATION:', TEST, I, J, K
                DO C = 1,64
                   PRINT*, C, BVEC(C)
                ENDDO
                STOP 1
             ENDIF

             ! multiply by the stored inverse B matrix to get the coefficients
             CALL DGEMV('N',64,64,1D0,BINVMAT,64,BVEC,1,0D0,COEFF,1)
             !if (splinetest) print*, 'testx1:', i,j,k,coeff
             SP%SPCOEFF(:,I,J,K) = COEFF
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE CREATESPLINE

  SUBROUTINE SETUPSPLINE(SP,N1,N2,N3)
    ! set up the various arrays for the cubic spline stuff
    IMPLICIT NONE
    TYPE(SPDATA3D), POINTER :: SP
    INTEGER, INTENT(IN) :: N1,N2,N3

    IF (SP%ARRAYSET) CALL CLEANUPSPLINE(SP)
    
    SP%N1 = N1; SP%N2 = N2; SP%N3 = N3

    ! allocate the arrays
    ALLOCATE(SP%SPCOEFF(64,N1-1,N2-1,N3-1),SP%DATAMAT(N1,N2,N3))
    ALLOCATE(SP%X1(N1),SP%X2(N2),SP%X3(N3))
    SP%ARRAYSET = .TRUE.

  END SUBROUTINE SETUPSPLINE

  SUBROUTINE CLEANUPSPLINE(SP)
    ! deallocate various spline arrays
    IMPLICIT NONE
    TYPE(SPDATA3D), POINTER :: SP

    IF (.NOT.SP%ARRAYSET) RETURN
    DEALLOCATE(SP%SPCOEFF,SP%DATAMAT,SP%X1,SP%X2,SP%X3)
    NULLIFY(SP%SPCOEFF,SP%DATAMAT,SP%X1,SP%X2,SP%X3)
    SP%ARRAYSET = .FALSE.
  END SUBROUTINE CLEANUPSPLINE

  SUBROUTINE SETBINVMAT
    ! set the inverse B matrix elements
    ! obtained using splinecoeff.pl from coeff.h downloaded from 
    ! http://gyre.cds.caltech.edu/pub/software/tricubic//source 
    ! (src/libtricubic/coeff.h once unpacked)

    BINVMAT(1,:) = (/ & 
         & 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(2,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(3,:) = (/ & 
         & -3,  3,  0,  0,  0,  0,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(4,:) = (/ & 
         & 2, -2,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(5,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(6,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(7,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(8,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(9,:) = (/ & 
         & -3,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(10,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(11,:) = (/ & 
         & 9, -9, -9,  9,  0,  0,  0,  0,  6,  3, -6, -3,  0,  0,  0,  0, & 
         &  6, -6,  3, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  4,  2,  2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(12,:) = (/ & 
         & -6,  6,  6, -6,  0,  0,  0,  0, -3, -3,  3,  3,  0,  0,  0,  0, & 
         & -4,  4, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -2, -2, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(13,:) = (/ & 
         & 2,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(14,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(15,:) = (/ & 
         & -6,  6,  6, -6,  0,  0,  0,  0, -4, -2,  4,  2,  0,  0,  0,  0, & 
         & -3,  3, -3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -2, -1, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(16,:) = (/ & 
         & 4, -4, -4,  4,  0,  0,  0,  0,  2,  2, -2, -2,  0,  0,  0,  0, & 
         &  2, -2,  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(17,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(18,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(19,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(20,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(21,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(22,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(23,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -3,  3,  0,  0,  0,  0,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(24,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  2, -2,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(25,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(26,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -1,  0,  0,  0,  0,  0 /)
    BINVMAT(27,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  9, -9, -9,  9,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  6,  3, -6, -3,  0,  0,  0,  0, & 
         &  6, -6,  3, -3,  0,  0,  0,  0,  4,  2,  2,  1,  0,  0,  0,  0 /)
    BINVMAT(28,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0, -6,  6,  6, -6,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0, -3, -3,  3,  3,  0,  0,  0,  0, & 
         & -4,  4, -2,  2,  0,  0,  0,  0, -2, -2, -1, -1,  0,  0,  0,  0 /)
    BINVMAT(29,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(30,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0 /)
    BINVMAT(31,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0, -6,  6,  6, -6,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0, -4, -2,  4,  2,  0,  0,  0,  0, & 
         & -3,  3, -3,  3,  0,  0,  0,  0, -2, -1, -2, -1,  0,  0,  0,  0 /)
    BINVMAT(32,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  4, -4, -4,  4,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  2,  2, -2, -2,  0,  0,  0,  0, & 
         &  2, -2,  2, -2,  0,  0,  0,  0,  1,  1,  1,  1,  0,  0,  0,  0 /)
    BINVMAT(33,:) = (/ & 
         & -3,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0, -2,  0,  0,  0, -1,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(34,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  0,  0,  3,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0, -2,  0,  0,  0, -1,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(35,:) = (/ & 
         & 9, -9,  0,  0, -9,  9,  0,  0,  6,  3,  0,  0, -6, -3,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  3, -3,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  4,  2,  0,  0,  2,  1,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(36,:) = (/ & 
         & -6,  6,  0,  0,  6, -6,  0,  0, -3, -3,  0,  0,  3,  3,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0, -4,  4,  0,  0, -2,  2,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0, -2, -2,  0,  0, -1, -1,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(37,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -3,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -2,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(38,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -3,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0, -2,  0,  0,  0, -1,  0,  0,  0 /)
    BINVMAT(39,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  9, -9,  0,  0, -9,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  6,  3,  0,  0, -6, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  6, -6,  0,  0,  3, -3,  0,  0,  4,  2,  0,  0,  2,  1,  0,  0 /)
    BINVMAT(40,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -6,  6,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -3, -3,  0,  0,  3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -4,  4,  0,  0, -2,  2,  0,  0, -2, -2,  0,  0, -1, -1,  0,  0 /)
    BINVMAT(41,:) = (/ & 
         & 9,  0, -9,  0, -9,  0,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  6,  0,  3,  0, -6,  0, -3,  0,  6,  0, -6,  0,  3,  0, -3,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  4,  0,  2,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(42,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  9,  0, -9,  0, -9,  0,  9,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  6,  0,  3,  0, -6,  0, -3,  0,  6,  0, -6,  0,  3,  0, -3,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  2,  0,  2,  0,  1,  0 /)
    BINVMAT(43,:) = (/ & 
         & -27, 27, 27, -27, 27, -27, -27, 27, -18, -9, 18,  9, 18,  9, -18, -9, & 
         & -18, 18, -9,  9, 18, -18,  9, -9, -18, 18, 18, -18, -9,  9,  9, -9, & 
         & -12, -6, -6, -3, 12,  6,  6,  3, -12, -6, 12,  6, -6, -3,  6,  3, & 
         & -12, 12, -6,  6, -6,  6, -3,  3, -8, -4, -4, -2, -4, -2, -2, -1 /)
    BINVMAT(44,:) = (/ & 
         & 18, -18, -18, 18, -18, 18, 18, -18,  9,  9, -9, -9, -9, -9,  9,  9, & 
         & 12, -12,  6, -6, -12, 12, -6,  6, 12, -12, -12, 12,  6, -6, -6,  6, & 
         &  6,  6,  3,  3, -6, -6, -3, -3,  6,  6, -6, -6,  3,  3, -3, -3, & 
         &  8, -8,  4, -4,  4, -4,  2, -2,  4,  4,  2,  2,  2,  2,  1,  1 /)
    BINVMAT(45,:) = (/ & 
         & -6,  0,  6,  0,  6,  0, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -3,  0, -3,  0,  3,  0,  3,  0, -4,  0,  4,  0, -2,  0,  2,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -2,  0, -2,  0, -1,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(46,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0, -6,  0,  6,  0,  6,  0, -6,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -3,  0, -3,  0,  3,  0,  3,  0, -4,  0,  4,  0, -2,  0,  2,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -2,  0, -1,  0, -1,  0 /)
    BINVMAT(47,:) = (/ & 
         & 18, -18, -18, 18, -18, 18, 18, -18, 12,  6, -12, -6, -12, -6, 12,  6, & 
         &  9, -9,  9, -9, -9,  9, -9,  9, 12, -12, -12, 12,  6, -6, -6,  6, & 
         &  6,  3,  6,  3, -6, -3, -6, -3,  8,  4, -8, -4,  4,  2, -4, -2, & 
         &  6, -6,  6, -6,  3, -3,  3, -3,  4,  2,  4,  2,  2,  1,  2,  1 /)
    BINVMAT(48,:) = (/ & 
         & -12, 12, 12, -12, 12, -12, -12, 12, -6, -6,  6,  6,  6,  6, -6, -6, & 
         & -6,  6, -6,  6,  6, -6,  6, -6, -8,  8,  8, -8, -4,  4,  4, -4, & 
         & -3, -3, -3, -3,  3,  3,  3,  3, -4, -4,  4,  4, -2, -2,  2,  2, & 
         & -4,  4, -4,  4, -2,  2, -2,  2, -2, -2, -2, -2, -1, -1, -1, -1 /)
    BINVMAT(49,:) = (/ & 
         & 2,  0,  0,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(50,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0, -2,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(51,:) = (/ & 
         & -6,  6,  0,  0,  6, -6,  0,  0, -4, -2,  0,  0,  4,  2,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0, -3,  3,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0, -2, -1,  0,  0, -2, -1,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(52,:) = (/ & 
         & 4, -4,  0,  0, -4,  4,  0,  0,  2,  2,  0,  0, -2, -2,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  2, -2,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  1,  1,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(53,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  2,  0,  0,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(54,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  2,  0,  0,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0 /)
    BINVMAT(55,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -6,  6,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -4, -2,  0,  0,  4,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -3,  3,  0,  0, -3,  3,  0,  0, -2, -1,  0,  0, -2, -1,  0,  0 /)
    BINVMAT(56,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  4, -4,  0,  0, -4,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  2,  2,  0,  0, -2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  2, -2,  0,  0,  2, -2,  0,  0,  1,  1,  0,  0,  1,  1,  0,  0 /)
    BINVMAT(57,:) = (/ & 
         & -6,  0,  6,  0,  6,  0, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -4,  0, -2,  0,  4,  0,  2,  0, -3,  0,  3,  0, -3,  0,  3,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -2,  0, -1,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(58,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0, -6,  0,  6,  0,  6,  0, -6,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         & -4,  0, -2,  0,  4,  0,  2,  0, -3,  0,  3,  0, -3,  0,  3,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -1,  0, -2,  0, -1,  0 /)
    BINVMAT(59,:) = (/ & 
         & 18, -18, -18, 18, -18, 18, 18, -18, 12,  6, -12, -6, -12, -6, 12,  6, & 
         & 12, -12,  6, -6, -12, 12, -6,  6,  9, -9, -9,  9,  9, -9, -9,  9, & 
         &  8,  4,  4,  2, -8, -4, -4, -2,  6,  3, -6, -3,  6,  3, -6, -3, & 
         &  6, -6,  3, -3,  6, -6,  3, -3,  4,  2,  2,  1,  4,  2,  2,  1 /)
    BINVMAT(60,:) = (/ & 
         & -12, 12, 12, -12, 12, -12, -12, 12, -6, -6,  6,  6,  6,  6, -6, -6, & 
         & -8,  8, -4,  4,  8, -8,  4, -4, -6,  6,  6, -6, -6,  6,  6, -6, & 
         & -4, -4, -2, -2,  4,  4,  2,  2, -3, -3,  3,  3, -3, -3,  3,  3, & 
         & -4,  4, -2,  2, -4,  4, -2,  2, -2, -2, -1, -1, -2, -2, -1, -1 /)
    BINVMAT(61,:) = (/ & 
         & 4,  0, -4,  0, -4,  0,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  2,  0,  2,  0, -2,  0, -2,  0,  2,  0, -2,  0,  2,  0, -2,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  1,  0,  1,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
    BINVMAT(62,:) = (/ & 
         & 0,  0,  0,  0,  0,  0,  0,  0,  4,  0, -4,  0, -4,  0,  4,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & 
         &  2,  0,  2,  0, -2,  0, -2,  0,  2,  0, -2,  0,  2,  0, -2,  0, & 
         &  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0,  1,  0 /)
    BINVMAT(63,:) = (/ & 
         & -12, 12, 12, -12, 12, -12, -12, 12, -8, -4,  8,  4,  8,  4, -8, -4, & 
         & -6,  6, -6,  6,  6, -6,  6, -6, -6,  6,  6, -6, -6,  6,  6, -6, & 
         & -4, -2, -4, -2,  4,  2,  4,  2, -4, -2,  4,  2, -4, -2,  4,  2, & 
         & -3,  3, -3,  3, -3,  3, -3,  3, -2, -1, -2, -1, -2, -1, -2, -1 /)
    BINVMAT(64,:) = (/ & 
         & 8, -8, -8,  8, -8,  8,  8, -8,  4,  4, -4, -4, -4, -4,  4,  4, & 
         &  4, -4,  4, -4, -4,  4, -4,  4,  4, -4, -4,  4,  4, -4, -4,  4, & 
         &  2,  2,  2,  2, -2, -2, -2, -2,  2,  2, -2, -2,  2,  2, -2, -2, & 
         &  2, -2,  2, -2,  2, -2,  2, -2,  1,  1,  1,  1,  1,  1,  1,  1 /)

  END SUBROUTINE SETBINVMAT
END MODULE TRICUBSPLINEUTILS
