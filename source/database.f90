MODULE DATABASEUTILS
  USE CHAINUTILS

  IMPLICIT NONE

  TYPE DATACHAIN
     !   object corresponding to a database of chain structures

     ! list of chain objects
     TYPE(CHAIN), POINTER :: CP(:)
     ! maximum number of chains kept in database
     INTEGER :: SIZEDB
     ! current number of chains in database
     INTEGER :: NCHAIN
     ! current minimum energy
     DOUBLE PRECISION :: MINE
     ! maximum energy allowed in this database
     DOUBLE PRECISION :: MAXESAVE

     ! various arrays are allocated
     LOGICAL :: ARRAYSET = .FALSE.

  END TYPE DATACHAIN

CONTAINS
  SUBROUTINE DBSORT(DBPT)
    USE KEYS, ONLY : SAMEECUT, SAMEDISTCUT
    ! resort the configurations in a database
    TYPE(DATACHAIN), POINTER :: DBPT
    TYPE(DATACHAIN), TARGET :: DB2
    TYPE(DATACHAIN), POINTER :: DB2PT
    TYPE(CHAIN), POINTER :: CHAINP2,CHAINP
    INTEGER :: C,RANK
    
    print*, 'Sorting database...'

    DB2PT=>DB2
    CALL SETUPDATABASE(DB2PT,DBPT%SIZEDB)
    DO C = 1,DBPT%NCHAIN
       CHAINP=>DBPT%CP(C)
       CALL ORDEREDINSERT(DB2PT,CHAINP,RANK,.TRUE.,SAMEECUT,SAMEDISTCUT)
       IF (RANK.EQ.0) PRINT*, 'Removing redundant structure ', C
    ENDDO

    DBPT%NCHAIN = DB2PT%NCHAIN
    DO C = 1,DBPT%NCHAIN
       CHAINP=>DBPT%CP(C)
       CHAINP2=>DB2PT%CP(C)
       CALL COPYCHAIN(CHAINP2,CHAINP,.TRUE.)
    ENDDO
    CALL CLEANUPDATABASE(DB2PT)

  END SUBROUTINE DBSORT

  SUBROUTINE DBMINIMIZE(DBPT,TOL)    
    USE OPTIMIZEUTILS, ONLY : BFGS
    ! minimize each config in the database
    TYPE(DATACHAIN), POINTER :: DBPT
    DOUBLE PRECISION, INTENT(IN) :: TOL
    LOGICAL :: SUCCESS
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER :: C

    print '(A,G10.5)', 'Minimizing all configurations in database to tolerance ', TOL

    DO C = 1,DBPT%NCHAIN
       CHAINP=>DBPT%CP(C)
       print*, 'Minimizing configuration ', C
       CALL BFGS(CHAINP,TOL,SUCCESS)
       IF (.NOT.SUCCESS) PRINT*, 'Failed to optimize configuration ', C       
    ENDDO
  END SUBROUTINE DBMINIMIZE

  SUBROUTINE DBDUMPCONFIGS(DBPT,NDUMP,DUMPFILE,REPLICFILE,NBENDREPLIC,MAXE)
    USE HELIXUTILS, ONLY : OUTPUTREPLIC
    ! dump out the first NDUMP configurations in a database
    ! the stored configurations are dumped into DUMPFILE
    ! the replicated configurations with NBENDREPLIC are dumped into REPLICFILE
    ! for both DUMPFILE and REPLICFILE, the last instance of # is replaced with the config index in the database
    TYPE(DATACHAIN), POINTER :: DBPT
    INTEGER, INTENT(IN) :: NDUMP
    CHARACTER(LEN=*), INTENT(IN) :: DUMPFILE
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: REPLICFILE
    INTEGER, INTENT(IN), OPTIONAL :: NBENDREPLIC
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: MAXE
    CHARACTER*10 :: NUMSTR
    CHARACTER*100 :: FILENAME,REPLICFILENAME
    INTEGER :: C, NBEND
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION :: CONFIGMAXE

    IF (PRESENT(MAXE)) THEN
       CONFIGMAXE = MAXE
    ELSE
       CONFIGMAXE = HUGE(1D0)
    ENDIF

    print*, 'Dumping out ', MIN(NDUMP,DBPT%NCHAIN), ' configs from database.'
    print*, 'Dump file is: ', dumpfile    

    IF (PRESENT(NBENDREPLIC)) THEN
       NBEND = NBENDREPLIC
    ELSE
       NBEND = 20
    ENDIF

    IF (PRESENT(REPLICFILE)) THEN
       PRINT*, 'Replicated structures dumped into ', trim(adjustl(replicfile)), ' with ', NBEND, ' bends.'
    ENDIF

    DO C = 1,MIN(NDUMP,DBPT%NCHAIN)
       CHAINP=>DBPT%CP(C)
       IF (CHAINP%ENERGY.GT.CONFIGMAXE) CYCLE       
       WRITE(NUMSTR,'(I10)') C
       FILENAME = DUMPFILE
       CALL REPLACESUBSTR(FILENAME,'#',TRIM(ADJUSTL(NUMSTR)))
       CALL OUTPUTCHAIN(CHAINP,FILENAME)

       IF (PRESENT(REPLICFILE)) THEN
          REPLICFILENAME = REPLICFILE
          CALL REPLACESUBSTR(REPLICFILENAME,'#',TRIM(ADJUSTL(NUMSTR)))
          CALL OUTPUTREPLIC(CHAINP,NBEND,REPLICFILENAME)
       ENDIF       
    ENDDO

  END SUBROUTINE DBDUMPCONFIGS

  SUBROUTINE DBDUMPINFO(DBPT, PRINTEPARTS,PRINTNUCDIST,PRINTFIBERDIAM)
    ! dump out information about the configurations in the database
    USE HELIXUTILS, ONLY : GETNUCDISTS, GETFIBERDIAM
    TYPE(DATACHAIN), POINTER :: DBPT
    LOGICAL, INTENT(IN), OPTIONAL :: PRINTNUCDIST, PRINTFIBERDIAM, PRINTEPARTS
    LOGICAL :: PNUCDIST, PFIBERDIAM,PEPARTS
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER :: C
    CHARACTER*5 :: NUMSTR
    CHARACTER*100 :: FMTSTR
    DOUBLE PRECISION :: NUCDISTS(DBPT%NCHAIN,DBPT%CP(1)%NINTERACT)
    DOUBLE PRECISION :: ESAVE


    PNUCDIST = .FALSE.; PFIBERDIAM = .FALSE.; PEPARTS = .FALSE.
    IF (PRESENT(PRINTEPARTS)) THEN
       PEPARTS = PRINTEPARTS
    ENDIF
    IF (PRESENT(PRINTNUCDIST)) THEN
       PNUCDIST = PRINTNUCDIST
    ENDIF
    IF (PRESENT(PRINTFIBERDIAM)) THEN
       PFIBERDIAM = PRINTFIBERDIAM
    ENDIF

    DO C = 1,DBPT%NCHAIN
       CHAINP=>DBPT%CP(C)       

       ! get the nucleosome distances for each configuration
       IF (PNUCDIST) THEN
          CALL GETNUCDISTS(CHAINP,NUCDISTS(C,:))
       ENDIF

       ! get the fiber diameter for each config
       IF (PFIBERDIAM) THEN
          CALL GETFIBERDIAM(CHAINP)
       ENDIF
    ENDDO


    ! Print out energies of each configuration
    PRINT*, 'CONFIG INFORMATION:'
    PRINT*, 'CHAIN,    ENERGY'
    PRINT*, '------------------------------------------'
    DO C = 1,DBPT%NCHAIN      
       print '(I10,G25.8)', C, DBPT%CP(C)%ENERGY
    ENDDO
    ! print out helical coordinates for each config
    IF (PFIBERDIAM) THEN
       PRINT*, '------------------------------------------'
       PRINT*, 'HEIGHT, ANGLE, RADIUS, ALPHA, BETA, GAMMA, DIAMETER'
       PRINT*, '------------------------------------------'
       DO C = 1,DBPT%NCHAIN
          print '(I10,7G15.5)', C, DBPT%CP(C)%VEC(1:6),DBPT%CP(C)%DIAM
       ENDDO
    ELSE
       PRINT*, '------------------------------------------'
       PRINT*, 'HEIGHT, ANGLE, RADIUS, ALPHA,BETA,GAMMA'
       PRINT*, '------------------------------------------'
       DO C = 1,DBPT%NCHAIN
          print '(I10,6G15.5)', C, DBPT%CP(C)%VEC(1:6)
       ENDDO
    ENDIF
    ! print out energy parts
    PRINT*, '------------------------------------------'
    IF (PEPARTS) THEN
       PRINT*, 'BENDING, TWISTING, STRETCHING, NUC-NUC STERICS, DNA-DNA STERICS, NUC-DNA STERICS'
       PRINT*, '------------------------------------------'
       DO C = 1,DBPT%NCHAIN
          print '(I10,6G15.5)', C, DBPT%CP(C)%ENERGYPARTS
       ENDDO
       PRINT*, '------------------------------------------'
    ENDIF
    IF (PNUCDIST) THEN
       PRINT*, 'NUCLEOSOME DISTANCES'
       PRINT*, '-------------------------------------------' 
       DO C = 1,DBPT%NCHAIN
          WRITE(NUMSTR,'(I5)') DBPT%CP(C)%NINTERACT
          FMTSTR = '(I10,'//TRIM(ADJUSTL(NUMSTR))//'G15.4)'
          print FMTSTR, C, NUCDISTS(C,:)
       ENDDO
       PRINT*, '------------------------------------------'
    ENDIF

  END SUBROUTINE DBDUMPINFO

  SUBROUTINE READDATABASE(DBPT,FILENAME,EXTRAINT,HELCRDFIX,HELCRD)
    ! read back in a database of configurations as spat out in DUMPDATABASE
    ! DBPT and the associated chain objects should be all set up
    ! only the actual configurations are altered
    ! WARNING: assumes that aside from the configurations, the parameters for the chains in the database were the same as those for the chains currently in DBPT
    ! returns an extra integer from the file as well if EXTRAINT provided
    ! also, assumes the datafile actually exists
    ! if HELCRDFIX is provided, some helix coordinates are fixed to the
    ! values given in HELCRD

    USE HELIXUTILS, ONLY : GETHELIXREP
    IMPLICIT NONE
    TYPE(DATACHAIN), POINTER :: DBPT
    CHARACTER (LEN=*), INTENT(IN) :: FILENAME
    INTEGER, INTENT(OUT), OPTIONAL :: EXTRAINT
    LOGICAL, INTENT(IN), OPTIONAL :: HELCRDFIX(6)
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: HELCRD(6)
    INTEGER :: C, B, TMPNBEAD, TMPNBEND, SIZEDBIN, DUMI, I
    TYPE(CHAIN), POINTER :: CHAINP
    TYPE(BEND), POINTER :: BP

    OPEN(UNIT=55,FILE=FILENAME,STATUS='OLD')
    IF (PRESENT(EXTRAINT)) THEN
       READ(55,'(2I10)',ERR=100,END=100) SIZEDBIN, EXTRAINT
    ELSE
       READ(55,'(2I10)',ERR=100,END=100) SIZEDBIN, DUMI
    ENDIF

    IF (SIZEDBIN.GT.DBPT%SIZEDB) THEN
       PRINT*, 'WARNING IN READDATABASE: more configurations are in database file &
            & than can fit in the current database. Will ignore the extra configs.', DBPT%SIZEDB, SIZEDBIN
       DBPT%NCHAIN = DBPT%SIZEDB
    ELSE
       DBPT%NCHAIN = SIZEDBIN
    ENDIF

    DO C = 1,DBPT%NCHAIN       
       CHAINP=>DBPT%CP(C)

       ! read in the number of beads, bends, and energy
       CHAINP%ENERGYPARTS = 0D0
       READ(55,'(2I10,7G25.8)',ERR=100,END=100) TMPNBEAD,TMPNBEND,CHAINP%ENERGY,CHAINP%ENERGYPARTS

       IF (CHAINP%ENERGY.LT.DBPT%MINE) DBPT%MINE = CHAINP%ENERGY

       IF (TMPNBEAD.NE.CHAINP%NBEAD.OR.TMPNBEND.NE.CHAINP%NBEND) THEN 
          PRINT*, 'ERROR IN READDATABASE: chain has wrong number of beads or bends.'
          PRINT*, 'Current number of beads, bends: ', CHAINP%NBEAD, CHAINP%NBEND
          PRINT*, 'File number of beads, bends:', TMPNBEAD, TMPNBEND
          STOP 1
       ENDIF

       ! read in the bead positions and orientations
       DO B = 1,CHAINP%NBEAD
          READ(55,'(7F25.8)',ERR=100,END=100) &
               & CHAINP%BEADS(B,:),CHAINP%BEADQ(B)%W, CHAINP%BEADQ(B)%X, CHAINP%BEADQ(B)%Y, CHAINP%BEADQ(B)%Z
       ENDDO

       ! read in bend positions and orientations
       DO B = 1,CHAINP%NBEND
          BP=>CHAINP%BENDS(B)
          READ(55,'(7F25.8)',ERR=100,END=100), &
               & BP%CENT, BP%ORIENT%W, BP%ORIENT%X, BP%ORIENT%Y, BP%ORIENT%Z
       ENDDO

       ! convert to the single-vector regular helix coordinates used in minimization and energy function
       CALL GETHELIXREP(CHAINP)

       IF (PRESENT(HELCRDFIX)) THEN
          IF (PRESENT(HELCRD)) THEN
             DO I = 1,6
                IF (HELCRDFIX(I)) THEN
                   CHAINP%VEC(I) = HELCRD(I)
                ENDIF
             ENDDO
          ELSE
             PRINT*, 'ERROR IN READDATABASE: cannot provided HELCRDFIX without HELCRD.'
             STOP 1          
          ENDIF
       ENDIF
    ENDDO
    CLOSE(55)    

    RETURN

100 CLOSE(55)
    PRINT*, 'ERROR IN READDATABASE: problem reading database file. Reached a premature end of file or an error'
    STOP 1
  END SUBROUTINE READDATABASE

  SUBROUTINE DUMPDATABASE(DBPT,FILENAME,EXTRAINT)
    ! dump out the configurations of a database into a single data file
    ! optionally, include an extra integer in the first intro line
    ! format is:
    ! number of configs, extrainteger (-1 if not surprised)    
    ! for each config: 
    ! number of beads, number of bends, energy
    ! then a line for the position and quaternion of each bead
    ! then a line for the center and orientation of each bend
    IMPLICIT NONE
    TYPE(DATACHAIN), POINTER :: DBPT
    CHARACTER (LEN=*), INTENT(IN) :: FILENAME
    INTEGER, INTENT(IN), OPTIONAL :: EXTRAINT
    INTEGER :: C, B
    TYPE(CHAIN), POINTER :: CHAINP
    TYPE(BEND), POINTER :: BP

    OPEN(UNIT=55,FILE=FILENAME,STATUS='UNKNOWN')
    ! write out the total number of structures and the extra integer
    IF (PRESENT(EXTRAINT)) THEN
       WRITE(55,'(2I10)') DBPT%NCHAIN, EXTRAINT
    ELSE
       WRITE(55,'(2I10)') DBPT%NCHAIN, -1
    ENDIF

    DO C = 1,DBPT%NCHAIN
       CHAINP=>DBPT%CP(C)
       WRITE(55,'(2I10,7G25.8)') CHAINP%NBEAD, CHAINP%NBEND,CHAINP%ENERGY, CHAINP%ENERGYPARTS

       ! write out bead data as position and quaternion
       DO B = 1,CHAINP%NBEAD
          WRITE(55,'(7F25.8)') CHAINP%BEADS(B,:),CHAINP%BEADQ(B)%W, CHAINP%BEADQ(B)%X, CHAINP%BEADQ(B)%Y, CHAINP%BEADQ(B)%Z
       ENDDO

       ! write out center and orientation for each bend
       DO B = 1,CHAINP%NBEND
          BP=>CHAINP%BENDS(B)
          WRITE(55,'(7F25.8)'), BP%CENT, BP%ORIENT%W, BP%ORIENT%X, BP%ORIENT%Y, BP%ORIENT%Z
       ENDDO

    ENDDO

    CLOSE(55)
  END SUBROUTINE DUMPDATABASE

  SUBROUTINE ORDEREDINSERT(DBPT,CHAINP,RANK,CHECKSAME,ECUT,RMSCUT)
    ! insert a new chain (CHAINP) into the database (DBPT), maintaining
    ! the order from lowest to highest energy
    ! if CHECKSAME is set, then also check whether CHAINP
    ! is the same as other configurations in the database
    ! preventing any duplicates in the resulting database
    ! DBPT%NCHAIN is updated accordingly on exit
    ! RANK contains the position of the new chain in the database    
    ! if a chain is the same as something in the database but of lower energy then
    ! it replaces its duplicate in the database; if it is the same as something in the database
    ! but of higher energy then it is not inserted
    ! if RANK == 0, then the chain was not inserted because it was 
    ! a higher-energy duplicate, or because there is no space for it in the database, 
    ! or because the chain energy is above the maximum allowed

    ! notation note: downstream goes to higher energies

    IMPLICIT NONE
    TYPE(DATACHAIN), POINTER :: DBPT
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(OUT) :: RANK
    LOGICAL, INTENT(IN) :: CHECKSAME
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: ECUT, RMSCUT
    LOGICAL :: SAMELOW, SAMEHIGH
    INTEGER :: RANK2, STARTIND, C
    TYPE(CHAIN), POINTER :: FROMCHAIN, TOCHAIN

    IF (CHECKSAME.AND.(.NOT.PRESENT(ECUT).OR..NOT.PRESENT(RMSCUT))) THEN
       PRINT*, 'ERROR IN ORDEREDINSERT: cannot check for duplicates if ecut and rmscut are not supplied'
       STOP 1
    ENDIF

    IF (CHAINP%ENERGY.GT.DBPT%MAXESAVE) THEN
       ! energy is too high; do not insert into database
       RANK = 0
       RETURN
    ENDIF

    SAMELOW = .FALSE.; SAMEHIGH = .FALSE.

    ! get the index right below where the chain would be inserted
    RANK = 0
    DO RANK = DBPT%NCHAIN,1,-1
       IF (CHAINP%ENERGY > DBPT%CP(RANK)%ENERGY) EXIT
    ENDDO   

    ! check if chain is the same as the lower energy config
    IF (CHECKSAME.AND.RANK.GT.0) THEN
       ! just in case, keep going down in energy and chains can no longer be the same
       ! due to energy cutoff
       DO RANK2 = RANK,1,-1          
          FROMCHAIN=>DBPT%CP(RANK2)
          CALL SAMECONFIG(FROMCHAIN,CHAINP,ECUT,RMSCUT,SAMELOW)
          IF (SAMELOW.OR.DBPT%CP(RANK2)%ENERGY.LT.CHAINP%ENERGY-ECUT) EXIT
       ENDDO
    ENDIF    

    ! If chain is same as something with lower energy, do not add to database
    IF (SAMELOW) THEN
       RANK = 0
    ELSE
       IF (CHECKSAME.AND.RANK.LT.DBPT%SIZEDB) THEN
          ! Check if chain is the same as the higher energy config
          DO RANK2 = RANK+1,DBPT%NCHAIN
             FROMCHAIN=>DBPT%CP(RANK2)
             CALL SAMECONFIG(FROMCHAIN,CHAINP,ECUT,RMSCUT,SAMEHIGH)
             IF (SAMEHIGH.OR.DBPT%CP(RANK2)%ENERGY.GT.CHAINP%ENERGY+ECUT) EXIT
          ENDDO
       ENDIF

       IF (SAMEHIGH) THEN
          ! chain is same as something of higher energy

          ! insert chain into database and remove the higher energy duplicate
          DO C = RANK2-1,RANK+1,-1
             FROMCHAIN=>DBPT%CP(C); TOCHAIN => DBPT%CP(C+1)
             CALL COPYCHAIN(FROMCHAIN,TOCHAIN,.FALSE.)
          ENDDO
          TOCHAIN=>DBPT%CP(RANK+1)
          CALL COPYCHAIN(CHAINP,TOCHAIN,.FALSE.)
          RANK = RANK + 1
       ELSE                    
          ! insert chain into database          
          STARTIND = MIN(DBPT%NCHAIN,DBPT%SIZEDB-1)
          DO C = STARTIND,RANK+1,-1                
             FROMCHAIN=>DBPT%CP(C); TOCHAIN => DBPT%CP(C+1)
             CALL COPYCHAIN(FROMCHAIN,TOCHAIN,.FALSE.)
          ENDDO
          IF (RANK.LT.DBPT%SIZEDB) THEN
             TOCHAIN=>DBPT%CP(RANK+1)
             CALL COPYCHAIN(CHAINP,TOCHAIN,.FALSE.)
             IF (DBPT%NCHAIN.LT.DBPT%SIZEDB) DBPT%NCHAIN = DBPT%NCHAIN + 1
             RANK = RANK + 1
          ELSE
             RANK = 0
          ENDIF
       ENDIF
    ENDIF
    
    IF (RANK.EQ.1) DBPT%MINE = CHAINP%ENERGY
  END SUBROUTINE ORDEREDINSERT

  SUBROUTINE SETUPDATABASE(DBPT,MAXNCHAIN,MAXESAVE)
    ! Allocate everything including the chain objects
    ! does not initialize chain conformations
    ! set it up to be able to hold a maximum of MAXNCHAIN chain objects
    ! MAXESAVE sets the corresponding field in the database, which is 
    ! used to limit the maximal energy of any configurations saved in this database
    IMPLICIT NONE
    TYPE(DATACHAIN), POINTER :: DBPT
    INTEGER, INTENT(IN) :: MAXNCHAIN
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: MAXESAVE
    INTEGER :: C
    TYPE(CHAIN), POINTER :: CHAINP

    DBPT%SIZEDB = MAXNCHAIN

    ALLOCATE(DBPT%CP(MAXNCHAIN))
    DO C = 1,MAXNCHAIN
       CHAINP=>DBPT%CP(C)       
       CALL CREATECHAIN(CHAINP)
    ENDDO

    DBPT%NCHAIN = 0
    DBPT%MINE = HUGE(1D0)
    IF (PRESENT(MAXESAVE)) THEN
       DBPT%MAXESAVE = MAXESAVE
    ELSE
       DBPT%MAXESAVE = HUGE(1D0)
    ENDIF

    DBPT%ARRAYSET = .TRUE.
  END SUBROUTINE SETUPDATABASE

  SUBROUTINE CLEANUPDATABASE(DBPT)
    TYPE(DATACHAIN), POINTER :: DBPT
    INTEGER :: C
    TYPE(CHAIN), POINTER :: CHAINP

    DO C = 1,DBPT%SIZEDB
       CHAINP=>DBPT%CP(C)
       CALL CLEANUPCHAIN(CHAINP)
    ENDDO

    DEALLOCATE(DBPT%CP)
    DBPT%ARRAYSET = .FALSE.
  END SUBROUTINE CLEANUPDATABASE
END MODULE DATABASEUTILS
