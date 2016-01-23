MODULE OPTIMIZEUTILS
  ! subroutines for optimizing chain configuration using quasi-newton method
CONTAINS
  SUBROUTINE BFGS(CHAINP,TOL,SUCCESS)
    ! find a local minimum for the configuration of the chain
    ! where the mean square force is less than TOL
    ! returns success=true if local minimum successfully found
    ! NOTE: this works with the single-vector representation of the chain
    ! the BEADS, and BEADQ arrays are only set for output and at the end
    ! of the optimization
    ! Follows algorithm from Nocedal & Wright, 1999, algorithm 8.1 )
    ! except no line search: simply insist the energy not rise beyond MAXEJUMP in each step

    USE CHAINUTILS, ONLY : CHAIN, OUTPUTCHAIN
    USE KEYS, ONLY : MAXOPTSTEP, DGUESS,MAXSTEPSIZE,MAXEJUMP,MAXDCR,STEPDCR,&
         & OPTFLEXTAILSONLY, USEFLEXTAILS, FIXHELIXPARAM, USEFIXDIAM, VERBOSE, &
         & OPTPRINTFREQ, DUMPCURRENT, DUMPSTEPS, DUMPFILE, DUMPNBEND, MAXOPTATTEMPT
    USE ENERGYUTILS, ONLY : HELENERGYGRAD
    USE HELIXUTILS, ONLY : FROMHELIXREP, REGULARIZEHELCOORDS, OUTPUTREPLIC, CHECKGIMBAL, ROTATECHAIN
    USE GENUTILS, ONLY : REPLACESUBSTR, RANDOMAXIS
    USE QUATUTILS, ONLY : QUATERNION, ROTQUAT, PI
    USE MT19937, ONLY : GRND

    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(IN) :: TOL
    LOGICAL, INTENT(OUT) :: SUCCESS
    DOUBLE PRECISION,DIMENSION(CHAINP%NCRD,CHAINP%NCRD) :: HMAT, SYHYS,SYH,HYS,SS
    DOUBLE PRECISION, DIMENSION(CHAINP%NCRD) :: NRDIR,HY,YH,YHYS, PREVGRAD,PREVVEC, S, Y
    DOUBLE PRECISION :: RMS,NG,PREVE, YS,YY, SLENGTH, STP, STARTSTP, RO,GAMMA,GRADPROJ
    INTEGER :: K, I, J, D, B, IND, ATTEMPT
    CHARACTER*158 ::FILENAME
    CHARACTER*6 :: TMPSTR
    DOUBLE PRECISION :: STARTVEC(CHAINP%NCRD), ROTAX(3), ANG
    LOGICAL :: BAD, GIMBAL,BEADHIT
    TYPE(QUATERNION) :: QROT

    STARTVEC = CHAINP%VEC
    SUCCESS = .FALSE.

    DO ATTEMPT = 1,MAXOPTATTEMPT
       PRINT*, 'Optimization attempt: ', ATTEMPT
       CALL HELENERGYGRAD(CHAINP)          

       ! fix certain helix parameters
       DO I = 1,6
          IF (FIXHELIXPARAM(I)) CHAINP%GRAD(I) = 0D0
       ENDDO

       IF (OPTFLEXTAILSONLY.AND.USEFLEXTAILS) THEN
          ! optimizing only the flexible tails
          CHAINP%GRAD(1:CHAINP%NCRDDNA) = 0D0
       END IF

       ! set the original estimate for HMAT (hessian matrix)
       HMAT = 0D0
       DO I = 1,CHAINP%NCRD
          HMAT(I,I) = DGUESS
       ENDDO

       STP = 1D0; STARTSTP = 1D0
       DO K = 0,MAXOPTSTEP      
          ! check for convergence
          NG = DOT_PRODUCT(CHAINP%GRAD,CHAINP%GRAD)
          RMS = NG/CHAINP%NCRD ! rms square force

          IF (MOD(K,OPTPRINTFREQ).EQ.0) THEN
             print '(A,I10,2F20.10,G30.10,F20.10)', 'K,STP,ENERGY,RMS, H:', K,STARTSTP, CHAINP%ENERGY,RMS, chainp%vec(1)
          ENDIF

          IF (RMS.LT.TOL) THEN
             ! successfully optimized
             CALL HELENERGYGRAD(CHAINP)          
             SUCCESS = .TRUE.
             EXIT
          ENDIF

          !dump current configuration to file if desired
          IF (DUMPCURRENT.AND.(MOD(K,DUMPSTEPS) == 0)) THEN  
             CALL FROMHELIXREP(CHAINP)

             WRITE(TMPSTR,'(I6)') K
             FILENAME = DUMPFILE
             CALL REPLACESUBSTR(FILENAME,'#',TRIM(ADJUSTL(TMPSTR)))

             IF (DUMPNBEND.LE.2) THEN
                CALL OUTPUTCHAIN(CHAINP,FILENAME)
             ELSE
                CALL OUTPUTREPLIC(CHAINP,DUMPNBEND,FILENAME)
             ENDIF
          ENDIF

          ! get the search direction
          DO I = 1,CHAINP%NCRD
             NRDIR(I) = -DOT_PRODUCT(HMAT(I,1:CHAINP%NCRD),CHAINP%GRAD)          
          ENDDO

          ! if search step has become too small, try again
          ! IF (DOT_PRODUCT(NRDIR,NRDIR).LT.EPSILON(1D0)) THEN
          !    IF (VERBOSE) PRINT*, 'Search direction has become too small. &
          !         & Reset hessian estimate and optimize again.', &
          !         & DOT_PRODUCT(NRDIR,NRDIR), DOT_PRODUCT(CHAINP%GRAD,CHAINP%GRAD)
          !    EXIT
          ! ENDIF

          ! if diameter is fixed, then R is a dependent variable of the other coords
          IF (USEFIXDIAM) NRDIR(3) = 0D0

          ! fix certain helix parameters
          DO I = 1,6
             IF (FIXHELIXPARAM(I)) NRDIR(I) = 0D0
          ENDDO

          IF (OPTFLEXTAILSONLY.AND.USEFLEXTAILS) THEN
             ! Only move the flexible tails
             NRDIR(1:CHAINP%NCRDDNA) = 0D0
          ENDIF

          ! projection of step direction onto gradient
          SLENGTH = SQRT(DOT_PRODUCT(NRDIR,NRDIR))
          GRADPROJ = DOT_PRODUCT(NRDIR,CHAINP%GRAD)/&
               & (SLENGTH*SQRT(DOT_PRODUCT(CHAINP%GRAD,CHAINP%GRAD)))       

          IF (GRADPROJ.GT.0) THEN
             IF (VERBOSE) print*, 'NRDIR has insufficient projection to negative gradient.Reversing step.'
             NRDIR = -NRDIR
          ENDIF

          ! make sure you're not taking more than the maximum step
          STP = MIN(MAXSTEPSIZE/SLENGTH,STARTSTP)       
          ! take the step
          PREVVEC = CHAINP%VEC
          PREVE = CHAINP%ENERGY       
          PREVGRAD = CHAINP%GRAD   

          ! decrease step size until energy no longer rising
          DO D = 1,MAXDCR
             BAD = .TRUE.
             DO WHILE (BAD)
                ! update the step vector
                S = STP*NRDIR
                ! do not allow negative radius
                BAD = (.NOT.FIXHELIXPARAM(3).AND.CHAINP%VEC(3)+S(3).lt.0)
                IF (BAD) THEN
                   IF (VERBOSE) PRINT*, 'decreasing step to avoid negative radius', &
                        & STP, CHAINP%VEC(1), CHAINP%VEC(3)
                   STP = STP*STEPDCR
                   IF (STP.LT.EPSILON(1D0)) THEN
                      PRINT*, 'Failed to take sufficiently small steps to avoid negative  radius'
                      print*, chainp%vec(1:6)
                      SUCCESS = .FALSE.
                      CHAINP%VEC = STARTVEC
                      CALL HELENERGYGRAD(CHAINP)
                      RETURN
                   ENDIF
                ENDIF
             ENDDO

             ! update chain config
             CHAINP%VEC = CHAINP%VEC + S
             !print*, 'testxa:', sqrt(dot_product(s,s)), sqrt(dot_product(nrdir,nrdir)), stp

             !IF (CHAINP%VEC(1).EQ.0D0.AND.CHAINP%VEC(2).EQ.0D0) THEN
             ! Fiber is essentially a single nucleosome
             ! check for gimbal lock if any linker segments lie along negative z axis
             CALL CHECKGIMBAL(CHAINP,GIMBAL,BEADHIT)
             IF (GIMBAL.OR.BEADHIT) THEN
                IF (GIMBAL) THEN
                   PRINT*, 'Chain hit gimbal lock.'
                ELSE
                   PRINT*, 'Two beads hit.'
                ENDIF
                CHAINP%VEC = CHAINP%VEC-S
                IF (CHAINP%VEC(1).EQ.0D0.AND.CHAINP%VEC(2).EQ.0D0) THEN
                   PRINT*, 'Rotate coordinate system and retry.'
                   ! rotate entire chain
                   CALL RANDOMAXIS((/1D0,0D0,0D0/),2D0,ROTAX)
                   ANG = GRND()*2*PI
                   QROT = ROTQUAT(ANG,ROTAX)
                   CALL ROTATECHAIN(CHAINP,QROT)
                ELSE
                   PRINT*, 'Reset hessian matrix estimate and retry'
                ENDIF
                SUCCESS = .FALSE.
                BAD = .TRUE.
                EXIT
             ENDIF
             !ENDIF


             ! get new energy and gradient
             CALL HELENERGYGRAD(CHAINP)

             ! fix certain helix parameters
             DO I = 1,6
                IF (FIXHELIXPARAM(I)) CHAINP%GRAD(I) = 0D0
             ENDDO

             IF (OPTFLEXTAILSONLY.AND.USEFLEXTAILS) THEN
                ! optimizing only the flexible tails
                CHAINP%GRAD(1:CHAINP%NCRDDNA) = 0D0
             END IF

             IF (CHAINP%ENERGY-PREVE.GT.MAXEJUMP) THEN
                IF (VERBOSE) print*, 'Energy increased, decreasing step size.', D, CHAINP%ENERGY, PREVE
                CHAINP%VEC = PREVVEC; 
                STP = STP*STEPDCR
             ELSE
                EXIT
             ENDIF
          ENDDO
          IF (BAD) EXIT

          IF (D.GT.MAXDCR) THEN
             ! too many decreases: reset hermitian estimate 
             CHAINP%GRAD = PREVGRAD
             CHAINP%VEC = PREVVEC          
             print*, 'LBFGS: failed to find lower energy -- too many decreases! Try again from this starting point.'
             SUCCESS=.FALSE.
             EXIT
          ENDIF

          ! update the initial guess for the stepsize
          IF (D.GT.1) THEN
             ! had to decrease step more than once: decrease initial stepsize
             STARTSTP = STARTSTP*STEPDCR
          ELSE IF (D.LE.1.AND.STARTSTP.LT.1D0) THEN
             ! did not have to decrease step: raise initial stepsize
             STARTSTP = STARTSTP/STEPDCR          
          ENDIF

          ! fix certain helix parameters
          DO I = 1,6
             IF (FIXHELIXPARAM(I)) CHAINP%GRAD(I) = 0D0
          ENDDO
          IF (OPTFLEXTAILSONLY.AND.USEFLEXTAILS) THEN
             CHAINP%GRAD(1:CHAINP%NCRDDNA) = 0D0
          ENDIF

          ! update the change in gradient
          Y = CHAINP%GRAD - PREVGRAD 

          ! update the hessian matrix estimate H_k+1
          YS = DOT_PRODUCT(S,Y)  
          IF (ABS(YS).LT.EPSILON(1D0)) YS = 1D0
          RO = 1D0/YS

          IF (K.EQ.0) THEN
             YY = DOT_PRODUCT(Y,Y)
             IF (YY.LT.EPSILON(1D0)) YY = 1D0
             ! fix the H_0 estimate (Nocedal 8.20)
             GAMMA = YS/YY
             DO I = 1,CHAINP%NCRD
                HMAT(I,I) = GAMMA
             ENDDO
          ENDIF

          DO I = 1,CHAINP%NCRD
             HY(I)  = DOT_PRODUCT(HMAT(I,:),Y)
             DO J = 1,CHAINP%NCRD
                HYS(I,J) = HY(I)*S(J)
             ENDDO
          ENDDO

          DO J = 1,CHAINP%NCRD
             YHYS(J) = DOT_PRODUCT(Y,HYS(:,J))
             YH(J) = DOT_PRODUCT(Y,HMAT(:,J))
             DO I = 1,CHAINP%NCRD
                SYHYS(I,J) = S(I)*YHYS(J)
                SYH(I,J) = S(I)*YH(J)
                SS(I,J) = S(I)*S(J)
             ENDDO
          ENDDO

          HMAT = HMAT + RO**2*SYHYS - RO*(SYH+HYS) + RO*SS

          RMS = DOT_PRODUCT(CHAINP%GRAD,CHAINP%GRAD)/CHAINP%NCRD
          IF (.NOT.RMS.GE.0) THEN
             PRINT*, 'ERROR IN BFGS: bad rms', RMS
             stop 1
          ENDIF
       ENDDO

       IF (SUCCESS) EXIT

       IF (K.GT.MAXOPTSTEP) THEN
          PRINT*, 'BFGS FAILED TO MINIMIZE IN THE MAXIMUM NUMBER OF STEPS.',K
          CALL FROMHELIXREP(CHAINP)
          SUCCESS = .FALSE.
          CALL HELENERGYGRAD(CHAINP)
          ! CALL TESTCURRENTGRADIENT(CHAINP)
          ! STOP 1
          RETURN          
       ENDIF
    ENDDO

    ! regularize the helix coordinates (no negative heights, etc)
    !CALL REGULARIZEHELCOORDS(CHAINP%vec(1:6))
    ! get final energy and gradient
    CALL HELENERGYGRAD(CHAINP)
    ! fix certain helix parameters
    DO I = 1,6
       IF (FIXHELIXPARAM(I)) CHAINP%GRAD(I) = 0D0
    ENDDO

    IF (OPTFLEXTAILSONLY.AND.USEFLEXTAILS) THEN
       ! optimizing only the flexible tails
       CHAINP%GRAD(1:CHAINP%NCRDDNA) = 0D0
    END IF

    PRINT '(a,2G20.10)', 'Energy, RMS force after optimization:', CHAINP%ENERGY, &
         & SQRT(DOT_PRODUCT(CHAINP%GRAD,CHAINP%GRAD)/CHAINP%NCRD)
    ! convert to output coordinates
    CALL FROMHELIXREP(CHAINP)

    ! if (.not.success) then
    !    CALL TESTCURRENTGRADIENT(CHAINP)
    !    STOP 1
    ! endif
  END SUBROUTINE BFGS

  SUBROUTINE TESTCURRENTGRADIENT(CHAINP)
    use keys, only : outfile
    USE CHAINUTILS, ONLY : CHAIN, OUTPUTCHAIN
    USE ENERGYUTILS, ONLY : HELENERGYGRAD, HELENERGYTEST
    USE INPUTSTRUCTS, ONLY : READOUTFILE
    USE HELIXUTILS, ONLY : DUMPHELIXCRD
    ! test the current gradient of the chain numerically
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION :: GRAD(CHAINP%NCRD), ENERGY
    DOUBLE PRECISION, PARAMETER :: TINY = 1D-8
    INTEGER :: I
    DOUBLE PRECISION :: STP,PREVVEC(CHAINP%NCRD)
    DOUBLE PRECISION :: RVAL, DRVAL(6), DEDR


    ! print*, 'current vector:'
    !  do i = 1,chainp%ncrd
    !     print*, i, chainp%vec(i)
    !  enddo
    !  stop 1

    CALL OUTPUTCHAIN(CHAINP,'badconfig.out')
    CALL readoutfile(chainp,'badconfig.out',2)

    PRINT*, 'TESTING GRADIENT NUMERICALLY:'

    HELENERGYTEST = .TRUE.
    CALL HELENERGYGRAD(CHAINP)
    HELENERGYTEST = .FALSE.

    PRINT*, 'ENERGY:', CHAINP%ENERGY
    print*, 'energyparts', chainp%energyparts, chainp%ncrd   
    ENERGY = CHAINP%ENERGY
    GRAD(1:CHAINP%NCRD) = CHAINP%GRAD(1:CHAINP%NCRD)
    CALL DUMPHELIXCRD(CHAINP,'badcrd.out')

    PREVVEC = CHAINP%VEC(:)       
    DO I = 1,chainp%ncrd
       CHAINP%VEC(I) = CHAINP%VEC(I) + TINY          
       CALL HELENERGYGRAD(CHAINP)          

       PRINT*, I, (CHAINP%ENERGY-ENERGY)/TINY, grad(i), energy, chainp%energy
       CHAINP%VEC(I) = CHAINP%VEC(I) - TINY              
    ENDDO


    call helenergygrad(chainp)
    PRINT*, 'NOW TAKE INCREASINGLY SMALLER STEPS ALONG -GRADIENT:'
    print*, 'start energy:', chainp%energy,chainp%energyparts
    STP = 1D-3; PREVVEC = CHAINP%VEC; ENERGY = CHAINP%ENERGY
    GRAD(1:CHAINP%NCRD) = CHAINP%GRAD(1:CHAINP%NCRD)
    DO I = 1,10
       CHAINP%VEC = PREVVEC
       CHAINP%VEC = CHAINP%VEC - STP*GRAD
       CALL HELENERGYGRAD(CHAINP)
       PRINT*, I, STP, ENERGY, CHAINP%ENERGY,chainp%energyparts
       STP = STP*0.1D0
    ENDDO

  END SUBROUTINE TESTCURRENTGRADIENT

END MODULE OPTIMIZEUTILS
