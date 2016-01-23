MODULE DRIVERS
  ! driver subroutines for running various calculations

  IMPLICIT NONE

CONTAINS 
  SUBROUTINE BASINHOPDRIVER(CHAINP)    
    ! run basinhopping calculation starting with the given chain configuration
    USE KEYS, ONLY : READDATAFILE, DATAFILE, BHTOL, NBASINHOP,NCONFIGSAVE,&
         & DODBMINIMIZE, DODBSORT, OPTOL, NCONFIGDUMP, DUMPFILE, REPLICFILE,&
         & NBENDREPLIC, DBCONFIGMAXE, PRINTEPARTS, PRINTNUCDIST, PRINTFIBERDIAM,&
         & NEWDATAFILE, REPLICATESTRUCT, ENERGYRECALC,HELCRDSET,&
         & HELH, HELT, HELR, HELA, HELB, HELG
    USE CHAINUTILS, ONLY : CHAIN
    USE DATABASEUTILS, ONLY : DATACHAIN, SETUPDATABASE, CLEANUPDATABASE, &
         & READDATABASE, DBDUMPCONFIGS, DBDUMPINFO, DBMINIMIZE, DBSORT, &
         & DUMPDATABASE
    USE BASINHOPUTILS, ONLY : BASINHOP
    USE ENERGYUTILS, ONLY : HELENERGYGRAD
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP    
    TYPE(DATACHAIN), TARGET :: DB
    TYPE(DATACHAIN), POINTER :: DBPT
    LOGICAL :: FILEEXISTS
    INTEGER :: STARTHOP, C

    DBPT => DB

    CALL SETUPDATABASE(DBPT,NCONFIGSAVE,DBCONFIGMAXE)

    STARTHOP = 0
    IF (READDATAFILE) THEN       
       ! read in the database file
       INQUIRE(FILE=DATAFILE,EXIST=FILEEXISTS)
       IF (FILEEXISTS) THEN
          PRINT*, 'Reading in database from ', TRIM(ADJUSTL(DATAFILE))
          CALL READDATABASE(DBPT,DATAFILE,STARTHOP,HELCRDSET,(/HELH,HELT,HELR,HELA,HELB,HELG/))
          PRINT*, 'Number of hops already done, number of desired hops ', STARTHOP, NBASINHOP

          IF (ENERGYRECALC) THEN
             ! recalculate energy of each configuration
             PRINT*, 'Recalculating energies for all configurations in database...'
             DBPT%MINE = HUGE(1D0)
             DO C = 1,DBPT%NCHAIN
                CHAINP=>DBPT%CP(C)          
                CALL HELENERGYGRAD(CHAINP)
                IF (CHAINP%ENERGY.LT.DBPT%MINE) DBPT%MINE = CHAINP%ENERGY
             ENDDO
          ENDIF
       ENDIF
    ENDIF

    IF (STARTHOP.LT.NBASINHOP) THEN
       PRINT*, 'Running basinhopping calculation...'
       CALL BASINHOP(CHAINP,DBPT,NBASINHOP,BHTOL, STARTHOP+1)
    ENDIF

    IF (DODBMINIMIZE) THEN
       ! minimize each structure
       CALL DBMINIMIZE(DBPT,OPTOL)
    ENDIF

    IF (DODBSORT) THEN
       ! resort the database
       CALL DBSORT(DBPT)
    ENDIF

    ! dump out configurations to file
    IF (REPLICATESTRUCT) THEN
       CALL DBDUMPCONFIGS(DBPT,NCONFIGDUMP,DUMPFILE,REPLICFILE,NBENDREPLIC,DBCONFIGMAXE)
    ELSE
       CALL DBDUMPCONFIGS(DBPT,NCONFIGDUMP,DUMPFILE,MAXE=DBCONFIGMAXE)
    ENDIF

    ! print out info
    CALL DBDUMPINFO(DBPT,PRINTEPARTS,PRINTNUCDIST,PRINTFIBERDIAM)

    IF (DODBMINIMIZE.OR.DODBSORT) THEN
       CALL DUMPDATABASE(DBPT,NEWDATAFILE,NBASINHOP)
    ENDIF

    ! deallocate database arrays
    CALL CLEANUPDATABASE(DBPT)

  END SUBROUTINE BASINHOPDRIVER

  SUBROUTINE DBPARSEDRIVER
    ! --------------------------
    ! parse a previously generated database, output configurations, and print data
    ! -------------------------
    USE KEYS, ONLY : NCONFIGDUMP,DUMPFILE,REPLICFILE,NBENDREPLIC,REPLICATESTRUCT,&
         & PRINTEPARTS,PRINTNUCDIST,PRINTFIBERDIAM, NCONFIGSAVE,DATAFILE,&
         & DBCONFIGMAXE,OPTOL,DODBMINIMIZE,HELCRDSET,HELH,HELT,HELR,HELA,HELB,HELG,&
         & NEWDATAFILE,DODBSORT,ENERGYRECALC
    USE DATABASEUTILS, ONLY : DBDUMPCONFIGS, DBDUMPINFO, SETUPDATABASE,&
         & CLEANUPDATABASE, DATACHAIN,READDATABASE,DBMINIMIZE,DBSORT, DUMPDATABASE
    USE CHAINUTILS, ONLY : CHAIN
    USE ENERGYUTILS, ONLY : HELENERGYGRAD
    IMPLICIT NONE
    TYPE(DATACHAIN), TARGET :: DB
    TYPE(DATACHAIN), POINTER :: DBPT
    DOUBLE PRECISION :: HELCRD(6)
    INTEGER :: I,c,EXTRAINT
    TYPE(CHAIN), POINTER :: CHAINP

    DBPT=>DB

    IF (ENERGYRECALC.OR.PRINTNUCDIST.OR.DODBMINIMIZE.OR.PRINTFIBERDIAM) CALL SETUPDRIVER

    ! set up the database and allocate all arrays
    CALL SETUPDATABASE(DBPT,NCONFIGSAVE,DBCONFIGMAXE)

    ! read in database
    CALL READDATABASE(DBPT,DATAFILE,EXTRAINT)

    HELCRD = (/HELH,HELT,HELR,HELA,HELB,HELG/)

    ! fix any present helix coords
    DO I = 1,6       
       IF (HELCRDSET(I)) THEN
          DO C = 1,DBPT%NCHAIN
             DBPT%CP(C)%VEC(I) = HELCRD(I)
          ENDDO
       ENDIF
    ENDDO

    IF (ENERGYRECALC) THEN
       ! recalculate energy of each configuration
       PRINT*, 'Recalculating energies for all configurations in database...'
       DBPT%MINE = HUGE(1D0)
       DO C = 1,DBPT%NCHAIN
          CHAINP=>DBPT%CP(C)          
          CALL HELENERGYGRAD(CHAINP)
          IF (CHAINP%ENERGY.LT.DBPT%MINE) DBPT%MINE = CHAINP%ENERGY
       ENDDO
    ENDIF

    IF (DODBMINIMIZE) THEN
       ! minimize each structure
       CALL DBMINIMIZE(DBPT,OPTOL)
    ENDIF

    IF (DODBSORT) THEN
       ! resort the database
       CALL DBSORT(DBPT)
    ENDIF

    IF (DODBMINIMIZE.OR.DODBSORT.OR.ENERGYRECALC) THEN
       CALL DUMPDATABASE(DBPT,NEWDATAFILE,EXTRAINT)
    ENDIF

    ! dump out configurations to file
    IF (REPLICATESTRUCT) THEN
       CALL DBDUMPCONFIGS(DBPT,NCONFIGDUMP,DUMPFILE,REPLICFILE,NBENDREPLIC,DBCONFIGMAXE)
    ELSE
       CALL DBDUMPCONFIGS(DBPT,NCONFIGDUMP,DUMPFILE, MAXE=DBCONFIGMAXE)
    ENDIF

    ! print out info
    CALL DBDUMPINFO(DBPT,PRINTEPARTS,PRINTNUCDIST,PRINTFIBERDIAM)

    ! dump out the database
    CALL DUMPDATABASE(DBPT,NEWDATAFILE,EXTRAINT)

    ! deallocate database arrays
    CALL CLEANUPDATABASE(DBPT)
    
     IF (ENERGYRECALC.OR.PRINTNUCDIST) CALL CLEANUPDRIVER

  END SUBROUTINE DBPARSEDRIVER

  SUBROUTINE OPTIMIZEDRIVER(CHAINP)
    ! -------------
    ! Find a local minimum for the chain structure (using BFGS optimization)
    ! -------------
    USE KEYS, ONLY : PRINTEPARTS, PRINTNUCDIST,PRINTFIBERDIAM, OPTOL, INTERMOPTOL, DISCORAMP,&
         & NOPTRUN,DISCORANGE,OPTFLEXTAILSONLY, OUTFILE
    USE CHAINUTILS, ONLY : CHAIN, OUTPUTCHAIN
    USE ENERGYUTILS, ONLY : HELENERGYGRAD
    USE OPTIMIZEUTILS, ONLY : BFGS
    USE HELIXUTILS, ONLY : GETFIBERDIAM, GETNUCDISTS
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION :: DISTS(CHAINP%NINTERACT)
    CHARACTER*5 :: NUMSTR
    CHARACTER*20 :: FMTSTR
    LOGICAL :: SUCCESS, TMPSWITCH
    INTEGER :: I
    DOUBLE PRECISION :: DISCOSTRENGTH

    print '(A,6G20.10)', 'Starting helix params (h, theta, R, alpha, beta, gamma):', chainp%vec(1:6)
    ! print out the starting energy
    CALL HELENERGYGRAD(CHAINP)
    PRINT '(A,G25.15)', ' Initial chain energy: ', CHAINP%ENERGY
    IF (PRINTEPARTS) THEN
       PRINT '(A,6G12.5)', ' Initial energy parts :(bend, twist, stretch, nuc-nuc sterics, &
            & seg-seg sterics, nuc-seg sterics) ', CHAINP%ENERGYPARTS
    ENDIF

    IF (DISCORAMP) THEN
       PRINT*, 'Ramping up disco potential in successive optimizations'
       DO I = 1,NOPTRUN 
          SUCCESS = .FALSE.
          IF (NOPTRUN.EQ.1) THEN
             DISCOSTRENGTH = 1D0
          ELSE
             DISCOSTRENGTH = 10d0**(-DISCORANGE + (DISCORANGE*(I-1))/(NOPTRUN-1))
          ENDIF          
          IF (I.LT.NOPTRUN) THEN
             print*, 'BFGS run number ', I, ' with disco strength ', DISCOSTRENGTH, ', tolerance ', INTERMOPTOL
             CALL BFGS(CHAINP,INTERMOPTOL,SUCCESS)               
          ELSE
             print*, 'BFGS run number ', I, ' with disco strength ', DISCOSTRENGTH, ', tolerance ', OPTOL
             CALL BFGS(CHAINP,OPTOL,SUCCESS)
          ENDIF 
       ENDDO
    ELSE
       print '(A,G25.15)', 'Running BFGS optimization on the structure. Tolerance: ', OPTOL
       CALL BFGS(CHAINP,OPTOL,SUCCESS)
    ENDIF
    IF (SUCCESS) print*, 'Successfully optimized!'
    PRINT '(A,6G12.5)', 'Final helix parameters:', CHAINP%VEC(1:6)
    print '(A,G25.15)', 'Final energy: ', CHAINP%ENERGY
    print*, 'Outputing final structure into file: ', OUTFILE
    CALL OUTPUTCHAIN(CHAINP,OUTFILE)

    IF (PRINTEPARTS.and..NOT.OPTFLEXTAILSONLY) THEN
       CALL HELENERGYGRAD(CHAINP)
       PRINT '(A,6G12.5)', ' Energy parts (bend, twist, stretch, &
            & nuc-nuc sterics, seg-seg sterics, nuc-seg sterics):', CHAINP%ENERGYPARTS
    ENDIF

    IF (PRINTFIBERDIAM) THEN
       CALL GETFIBERDIAM(CHAINP)
       PRINT '(A,G25.15)', ' Overall fiber diameter : ', chainp%diam
    ENDIF

    ! if optimized flexible tails only, print the overall energy as well
    TMPSWITCH = .FALSE.
    IF (OPTFLEXTAILSONLY) THEN
       TMPSWITCH = .TRUE.
       OPTFLEXTAILSONLY = .FALSE.
       CALL HELENERGYGRAD(CHAINP)
       PRINT '(A,2G20.10)', 'Final full energy: ', CHAINP%ENERGY     
       IF (PRINTEPARTS) THEN
          PRINT '(A,6G12.5)', ' Energy parts (bend, twist, stretch, &
               & nuc-nuc sterics, seg-seg sterics, nuc-seg sterics):', CHAINP%ENERGYPARTS
       ENDIF
    ENDIF
    
    IF (TMPSWITCH) OPTFLEXTAILSONLY = .NOT.OPTFLEXTAILSONLY

    ! print distances to other nucleosome shells
    IF (PRINTNUCDIST.AND.CHAINP%NINTERACT.GT.0) THEN
       CALL GETNUCDISTS(CHAINP,DISTS)
       PRINT*, 'Nucleosome distances:'
       WRITE(NUMSTR,'(I5)') CHAINP%NINTERACT
       FMTSTR = '('//TRIM(ADJUSTL(NUMSTR))//'G15.5)'
       PRINT FMTSTR, DISTS
    ENDIF
    
  END SUBROUTINE OPTIMIZEDRIVER

  SUBROUTINE GETSTRUCTDRIVER(CHAINP)
    ! ----------------
    ! dump out starting structure and various info about it
    ! ------------------
    ! chainp should already be set up and initialized
    USE CHAINUTILS, ONLY : CHAIN, OUTPUTCHAIN
    USE HELIXUTILS, ONLY : GETFIBERDIAM, GETNUCDISTS
    USE ENERGYUTILS, ONLY : HELENERGYGRAD
    USE KEYS, ONLY : OUTFILE, PRINTFIBERDIAM, PRINTNUCDIST,PRINTEPARTS
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION :: DISTS(CHAINP%NINTERACT)
    CHARACTER*5 :: NUMSTR
    CHARACTER*20 :: FMTSTR

    PRINT*, 'Getting the starting structure and dumping into file: ', OUTFILE
    CALL OUTPUTCHAIN(CHAINP,OUTFILE)
    PRINT '(A,6G12.5)', ' Regular helix coordinates (h,theta,R,alpha,beta,gamma):  ', CHAINP%VEC(1:6)
    IF (PRINTFIBERDIAM) THEN
       CALL GETFIBERDIAM(CHAINP)
       PRINT '(A,G25.15)', ' Overall fiber diameter : ', chainp%diam
    ENDIF

    ! print out the energy
    CALL HELENERGYGRAD(CHAINP)
    PRINT '(A,G25.15)', ' Chain energy: ', CHAINP%ENERGY
    IF (PRINTEPARTS) THEN        
       PRINT '(A,6G20.10)', ' Energy parts (bend, twist, stretch, &
            & nuc-nuc sterics, seg-seg sterics, nuc-seg sterics):', CHAINP%ENERGYPARTS
    ENDIF

    ! print distances to other nucleosome shells
    IF (PRINTNUCDIST.AND.CHAINP%NINTERACT.GT.0) THEN
       CALL GETNUCDISTS(CHAINP,DISTS)
       PRINT*, 'Nucleosome distances:'
       WRITE(NUMSTR,'(I5)') CHAINP%NINTERACT
       FMTSTR = '('//TRIM(ADJUSTL(NUMSTR))//'G15.5)'
       PRINT FMTSTR, DISTS
    ENDIF

  END SUBROUTINE GETSTRUCTDRIVER

  SUBROUTINE SETUPDRIVER(CHAINP)
    ! Do all the setup necessary to allocate and the various arrays used in calculations
    ! as well as initializing the chain configuration and getting the cylinder steric tables
    ! if a pointer to a chain object is not supplied, then only do the cylinder tables

    USE KEYS, ONLY : USEFIXDIAM,RESTARTDUMP, DUMPFILE, RESTARTFROMFILE, &
         & RESTARTFILE,RESTARTBEADS, PRINTNUCDIST, PRINTFIBERDIAM,BENDSTERICS, &
         & SEGSTERICS, BENDSEGSTERICS, NSEGPERLINK, BENDSEP, STARTFROMHELIX, &
         & HELH, HELT, HELR, HELA, HELB, HELG, HELCRDSET, STERRADIUS, STERHEIGHT, DNARAD
    USE CHAINUTILS, ONLY : CHAIN, CREATECHAIN,OUTPUTCHAIN
    USE HELIXUTILS, ONLY : MAKEHELIX, GETHELIXREP, FROMHELIXREP, MAKESTRAIGHTLINKERHELIX
    USE INPUTSTRUCTS, ONLY : READOUTFILE
    USE CYLINDERUTILS, ONLY : SETUPCYLARRAYS

    IMPLICIT NONE
    TYPE(CHAIN), POINTER, OPTIONAL :: CHAINP
    INTEGER :: I
    LOGICAL :: CHAINDONE, DUMPFILEEXIST
    DOUBLE PRECISION :: HELCRD(6)

    ! set up the cylinder arrays
    CALL SETUPCYLARRAYS(BENDSTERICS.OR.PRINTNUCDIST.OR.USEFIXDIAM.OR.PRINTFIBERDIAM,&
         & SEGSTERICS,BENDSEGSTERICS,STERRADIUS,STERHEIGHT,DNARAD,BENDSEP/NSEGPERLINK)  

    IF (.NOT.PRESENT(CHAINP)) RETURN

    ! Set up the chain arrays
    CALL CREATECHAIN(CHAINP)

    HELCRD = (/HELH,HELT,HELR,HELA,HELB,HELG/)

    ! initialize the chain configuration
    CHAINDONE = .FALSE.
    IF (RESTARTDUMP) THEN
       ! read in from DUMPFILE if such a file exists
       INQUIRE(FILE=DUMPFILE,EXIST=DUMPFILEEXIST)
       IF (DUMPFILEEXIST) THEN
          PRINT*, 'Will restart from dump file ', trim(adjustl(DUMPFILE))
          CALL READOUTFILE(CHAINP,DUMPFILE,RESTARTBEADS)
          CALL GETHELIXREP(CHAINP)
          CHAINDONE = .TRUE.
       ENDIF
    ENDIF
    IF (.NOT.CHAINDONE) THEN
       IF (RESTARTFROMFILE) THEN
          ! read in initial chain configuration from a file
          CALL READOUTFILE(CHAINP,RESTARTFILE,RESTARTBEADS)
          CALL GETHELIXREP(CHAINP)

       ELSEIF (STARTFROMHELIX) THEN
          ! start from a particular helix configuration
          CALL MAKEHELIX(CHAINP,HELH,HELT,HELR,(/HELA,HELB,HELG/))
          CALL GETHELIXREP(CHAINP,.TRUE.)
       ELSE
          ! start with straight-linker configuration
          CALL MAKESTRAIGHTLINKERHELIX(CHAINP,HELCRDSET,HELCRD)
       ENDIF
    ENDIF
    ! overwrite any desired helix coordinates    
    DO I = 1,6
       IF (HELCRDSET(I)) CHAINP%VEC(I) = HELCRD(I)
    ENDDO

    PRINT*, 'Chain has been set up with helix coordinates:', CHAINP%VEC(1:6)
  END SUBROUTINE SETUPDRIVER

  SUBROUTINE CLEANUPDRIVER(CHAINP)
    ! -------------
    ! clean up cylinder arrays and chain
    ! -------------

    USE CYLINDERUTILS, ONLY : CLEANUPCYLARRAYS
    USE CHAINUTILS, ONLY : CHAIN, CLEANUPCHAIN

    IMPLICIT NONE
    TYPE(CHAIN), POINTER, OPTIONAL :: CHAINP

    CALL CLEANUPCYLARRAYS
    IF (PRESENT(CHAINP)) THEN
       CALL CLEANUPCHAIN(CHAINP)
    ENDIF    

  END SUBROUTINE CLEANUPDRIVER

END MODULE DRIVERS
