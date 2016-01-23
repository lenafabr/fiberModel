SUBROUTINE READKEY
  ! this subroutine reads in keywords from a parameter file
  ! it sets the various global variables defined in KEYS module
  ! name of the parameter file is param.* where * is a keyword argument
  ! if no keyword argument is supplied, the default is just a file called param
  ! The EXTRAPARAMFILES keyword will allow extra parameter files to be 
  ! read in as well

  USE KEYS
  USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO
  USE GENUTILS
  USE QUATUTILS, ONLY : PI

  IMPLICIT NONE

  ! ---- stuff for inputing the parameter file in free format --------
  CHARACTER*100 :: ARG ! command line argument
  INTEGER :: NUMARG ! number of command line arguments
  INTEGER :: NITEMS ! number of items on the line in the parameter file
  INTEGER :: PF ! input file unit
  LOGICAL :: FILEEND=.FALSE. ! done reading file?
  CHARACTER*100 :: WORD ! keyword
  ! -------------- for reading multiple parameter files --------  
  INTEGER, PARAMETER :: MAXNFILES = 10
  CHARACTER*100 :: PARAMFILES(MAXNFILES)
  INTEGER :: NPARAMFILES, NPARAMREAD
  ! ------ for initializing stuff
  INTEGER :: TIMEVAL(8), SEED
  DOUBLE PRECISION :: ROTMAT(3,3)
  ! ---------------- temporary variables ---------------
  INTEGER :: DUMI, I, TMPI, DUMI1, DUMI2, DUMI3
  LOGICAL :: FILLINBEADSET=.FALSE., DATAFILENEWSET = .FALSE., NUCSEGCYLFILESET=.false.
  LOGICAL :: SPTWISTSET = .FALSE., LDUM, NSEGSET = .FALSE.
  CHARACTER*100 :: DUMSTR
  ! ------------------------
  ! set variable defaults
  ! ------------------------
  ACTION = 'NONE'
  RNGSEED = 0
  BPLEN = 0.34D0
  VERBOSE = .FALSE.

  ! geometry parameters
  NTAIL1SEG = 0
  NTAIL2SEG = 0
  BENDSEP = 17D0
  TWIST = 2*PI/(10.46*BPLEN)
  HELH = 1D0
  HELT = 1D0
  HELR = 12D0
  HELA = 0D0
  HELB = 0.0D0
  HELG = 0D0
  HELCRDSET = .FALSE.
  BENDTANM = (/0D0,0D0,1D0/)
  BENDTANP = (/0D0,0D0,1D0/)
  XAXM = (/1D0,0D0,0D0/)
  XAXP = (/1D0,0D0,0D0/)
  POSM = (/-1D0,0D0,0D0/)
  POSP = (/1D0,0D0,0D0/)
  SAMEECUT = 0.5D0
  SAMEDISTCUT = 1D0

  FILLINBEADS = 0
  FIXDIAMETER = -1D0
  USEFIXDIAM = .FALSE.
  FIXHELIXPARAM = .FALSE.
  SPIRALBENDS = .FALSE.
  SPHAND = -1
  SPTWIST = TWIST
  SPTWIST0 = 0D0
  SPRADIUS = 4.18D0
  SPHEIGHT = 2.39D0
  SPLEN = 146*BPLEN

  ! energetic parameters  
  LP = 50D0
  LTW = 110D0
  LSTRETCH = 268.2927D0
  TENSION = 0D0

  MAXBENDINTER = 20
  CYLINDERFILE = 'cylinder?.bin'
  CYLINDERFROMFILE = .TRUE.
  CYLINDERPTS = 50
  SEGCYLFROMFILE = .TRUE.
  NUCSEGCYLFROMFILE = .TRUE.
  SEGCYLFILE = 'cylinder?.bin'
  NUCSEGCYLFILE = 'cylinder?.bin'
  DNARAD = 1D0
  BENDSTERICS = .TRUE.
  SEGSTERICS = .TRUE.
  BENDSEGSTERICS = .TRUE.
  STERRADIUS = 5.2D0
  STERHEIGHT = 5.5D0
  ESTERIC = 1D3

  ! DiSCO potential
  DISCOPOTENTIAL = .FALSE.
  USEFLEXTAILS = .FALSE.
  DISCORAMP = .FALSE.
  DISCORANGE = 10  
  NCHARGE = 0
  NTAIL = 0
  TAILNBEAD = 0
  TAILSTARTSET = .FALSE.
  TAILPOS = 0D0
  TAILBEADINFO = 0D0

  ! Langowski potential
  LANGOWSKIPOT = .FALSE.
  SIG0 = 10.3D0
  EPS0 = 0.25D0
  LJX = -0.506D0
  LJXP = 0.383D0

  ! optimization
  OPTPRINTFREQ = 100
  DUMPCURRENT = .FALSE.
  DUMPSTEPS = 100
  DUMPNBEND = 2
  NOPTRUN = 1
  INTERMOPTOL = 1D-4
  OPTOL = 1D-4
  DGUESS = 1D-2
  STEPDCR = 0.1D0
  MAXDCR = 10
  MAXEJUMP = 1D-9
  MAXOPTSTEP = 10000
  MAXSTEPSIZE = 1D0
  OPTFLEXTAILSONLY = .FALSE.
  MAXOPTATTEMPT = 10

  ! basin hopping
  SEGHOPEVERY = -1
  FLIPNUCMC = .FALSE.
  FLIPNUCPROB = 0.5D0
  MAXHOPTRY = 10
  SEGHOPEVERY = 0
  BHFACC = 0.5D0
  BHFSAME = 0.2D0
  BHCHECKRANGE = 50
  NBASINHOP = 100
  SAMEECUT = 0.5D0
  SAMEDISTCUT = 1D0
  BHTOL = 1D-4
  STARTRRANGE = 1D0
  STARTTRANGE = 2*PI
  STARTARANGE = 2D0
  BHTEMP = 1D0
  READDATAFILE = .false.

  ! parsing databases
  DBCONFIGMAXE = HUGE(1D0)
  NCONFIGSAVE = 1000
  NCONFIGDUMP = 1
  DATAFILE = '*.data.out'
  NEWDATAFILE = '*.datanew.out'
  DODBMINIMIZE = .FALSE.
  DODBSORT = .FALSE.
  ENERGYRECALC = .FALSE.

  ! input/output
  DUMPFILE = '*.dump.out'
  OUTFILE = '*.out'
  PRINTNUCDIST = .FALSE.
  PRINTEPARTS = .FALSE.
  PRINTFIBERDIAM = .FALSE.
  REPLICATESTRUCT = .FALSE.
  NBENDREPLIC = 20
  REPLICFILE = '*.replic.out'
  RESTARTFROMFILE = .FALSE.
  RESTARTDUMP = .FALSE.
  VERBOSE = .FALSE.
  RANDSTART = .FALSE.
  RSRANGE = 1D-4
  RESTARTBEADS = 2

  ! -------------------------
  ! Read in all parameter files, starting with the ones specified on command line
  ! --------------------------

  PF = 55 ! i/o unit number to be used for parameter files

  ! get input parameter files from command line
  NPARAMFILES = 0
  NUMARG = COMMAND_ARGUMENT_COUNT()  
  IF (NUMARG==0) THEN
     NPARAMFILES = 1
     PARAMFILES(1) = 'param'
     ARG = ''
  ELSE
     DO I = 1,NUMARG
        CALL GETARG(I, ARG)
        NPARAMFILES = NPARAMFILES + 1
        WRITE(DUMSTR,'(A)') 'param.' //TRIM(ADJUSTL(ARG))
        PARAMFILES(NPARAMFILES) = DUMSTR
     ENDDO
     ! reset arg to its original value
     IF (NUMARG.GT.1) CALL GETARG(1,ARG)
  ENDIF

  NPARAMREAD = 0 ! keep track of how many files have been read
  DO WHILE (NPARAMREAD.LT.NPARAMFILES)
     NPARAMREAD = NPARAMREAD + 1

     PRINT*, 'Reading parameter file: ', PARAMFILES(NPARAMREAD)
     INQUIRE(FILE=PARAMFILES(NPARAMREAD),EXIST=LDUM)
     IF (.NOT.LDUM) THEN
        PRINT*, 'ERROR in READKEY: Parameter file ', TRIM(ADJUSTL(PARAMFILES(NPARAMREAD))), ' does not exist.'
        STOP 1
     ENDIF
     OPEN(UNIT=PF, FILE=PARAMFILES(NPARAMREAD), STATUS='OLD')

     ! read in the keywords one line at a time
     DO 
        CALL READLINE(PF,FILEEND,NITEMS)
        IF (FILEEND.and.nitems.eq.0) EXIT

        ! skip empty lines
        IF (NITEMS.EQ.0) CYCLE

        ! Read in the keyword for this line
        CALL READA(WORD,CASESET=1)

        ! Skip any empty lines or any comment lines
        IF (WORD(1:1).EQ.'#') CYCLE

        SELECT CASE(WORD) ! pick which keyword
        CASE('ACTION')
           CALL READA(ACTION, CASESET=1)           
        CASE('BENDTANM')
           CALL READF(BENDTANM(1))
           CALL READF(BENDTANM(2))
           CALL READF(BENDTANM(3))
        CASE('BENDTANP')
           CALL READF(BENDTANP(1))
           CALL READF(BENDTANP(2))
           CALL READF(BENDTANP(3))             
        CASE('BHCHECKRANGE')
           CALL READI(BHCHECKRANGE)
        CASE('BHFACC')
           CALL READF(BHFACC)
        CASE('BHFSAME')
           CALL READF(BHFSAME)
        CASE('BHTEMP')
           CALL READF(BHTEMP)
        CASE('BHTOL')
           CALL READF(BHTOL)
        CASE('BPLEN')
           CALL READF(BPLEN)     
        CASE('CHECKRESTART')
           CHECKRESTART = .TRUE.
        CASE('CYLINDERPTS')           
           CALL READI(DUMI1); CYLINDERPTS(1) = DUMI1
           CALL READI(DUMI2); CYLINDERPTS(2) = DUMI2
           CALL READI(DUMI3); CYLINDERPTS(3) = DUMI3
        CASE('DATAFILE')
           CALL READA(DATAFILE)
           IF (NITEMS.GT.2) CALL READA(NEWDATAFILE)
        CASE('DBCONFIGMAXE')
           CALL READF(DBCONFIGMAXE)
        CASE('DBMINIMIZE')
           DODBMINIMIZE = .TRUE.
        CASE('DBSORT')
           DODBSORT = .TRUE.
        CASE('DISCOPOTENTIAL')
           DISCOPOTENTIAL = .TRUE.
           IF (NITEMS.GT.1) CALL READF(LDEBYEINV)
           IF (NITEMS.GT.2) CALL READF(EPSEV)
           IF (NITEMS.GT.3) CALL READF(EPSEVT)
           IF (NITEMS.GT.4) CALL READF(SIGCC)
           IF (NITEMS.GT.5) CALL READF(SIGCL)
           IF (NITEMS.GT.6) CALL READF(SIGTC)
           IF (NITEMS.GT.7) CALL READF(SIGTL)
           IF (NITEMS.GT.8) CALL READF(SIGTT)
        CASE ('DISCOCHARGE') ! charge position and magnitude
           NCHARGE = NCHARGE + 1
           IF (NCHARGE.GT.MAXNCHARGE) THEN
              PRINT*, 'ERROR: TOO MANY DISCO CHARGES', NCHARGE
              STOP 1
           ENDIF
           CALL READF(CHARGEPOS(NCHARGE,1))
           CALL READF(CHARGEPOS(NCHARGE,2))
           CALL READF(CHARGEPOS(NCHARGE,3))
           CALL READF(CHARGEMAG(NCHARGE))
        CASE('DISCORAMP')
           DISCORAMP = .TRUE.
           IF (NITEMS.GT.1) CALL READF(DISCORANGE)
           IF (NITEMS.GT.2) CALL READI(NOPTRUN)
           IF (NITEMS.GT.3) CALL READF(INTERMOPTOL)
        CASE('DNARAD')
           CALL READF(DNARAD)
        CASE('DUMPCURRENT')
           DUMPCURRENT = .TRUE.
           IF (NITEMS.GT.1) CALL READI(DUMPSTEPS)
           IF (NITEMS.GT.2) CALL READI(DUMPNBEND)
        CASE('DUMPFILE')
           CALL READA(DUMPFILE)   
        CASE('ENERGYRECALC')
           ENERGYRECALC = .TRUE.
        CASE('ESTERIC') 
           CALL READF(ESTERIC)
        CASE('EXTRAPARAMFILES')
           DO DUMI = 1,MAXNFILES-1
              IF (NITEMS.GT.DUMI) THEN
                 NPARAMFILES = NPARAMFILES + 1
                 CALL READA(PARAMFILES(NPARAMFILES))
                 CALL REPLACESUBSTR(PARAMFILES(NPARAMFILES),'*',TRIM(ADJUSTL(ARG)))
              ELSE
                 EXIT
              ENDIF
           ENDDO
        CASE('FILLBEAD')
           FILLINBEADS = FILLINBEADS + 1
           IF (FILLINBEADS.GT.MAXFILLBEADS) THEN
              print*, 'ERROR: too many fill beads. Maximum is ', MAXFILLBEADS, ' . Ignoring this one.'
           ELSE
              CALL READF(FILLBEADLIST(FILLINBEADS,1))
              CALL READF(FILLBEADLIST(FILLINBEADS,2))
              CALL READF(FILLBEADLIST(FILLINBEADS,3))
              CALL READF(FILLBEADLIST(FILLINBEADS,4))
              CALL READF(FILLBEADLIST(FILLINBEADS,5))
              CALL READF(FILLBEADLIST(FILLINBEADS,6))
           ENDIF
           FILLINBEADSET = .TRUE.
        CASE('FIXDIAMETER')
           USEFIXDIAM = .TRUE.
           CALL READF(FIXDIAMETER)
        CASE('FIXHELIXPARAM')
           DO I = 1,NITEMS-1
              CALL READI(TMPI)
              IF (TMPI.LT.1.OR.TMPI.GT.6) THEN
                 PRINT*, 'ERROR: CANNOT FIX HELIX PARAMETER ', TMPI
              ELSE
                 FIXHELIXPARAM(TMPI) = .TRUE.
              ENDIF
           ENDDO
        CASE('FLIPNUCMC')
           FLIPNUCMC = .TRUE.
           IF (NITEMS.GT.1) CALL READF(FLIPNUCPROB)
        CASE('LANGOWSKI')
           LANGOWSKIPOT = .TRUE.
           IF (NITEMS.GT.1) THEN
              CALL READF(SIG0)
              CALL READF(EPS0)
              CALL READF(LJX)
              CALL READF(LJXP)
           ENDIF
        CASE('LINKLEN')
           CALL READF(BENDSEP)
        CASE('LP')
           CALL READF(LP)           
        CASE('LSTRETCH')
           CALL READF(LSTRETCH)
        CASE('LTW') 
           CALL READF(LTW) 
        CASE('MARKPOINT')
           IF (NMARKPOINT.GT.MAXMARKPOINT) THEN
              PRINT*, 'Too many mark points specified. Maximum number of points is ', &
                   & MAXMARKPOINT, '. Ignoring all the rest.' 
           ELSE
              NMARKPOINT = NMARKPOINT + 1
              CALL READA(MARKPOINTNAME(NMARKPOINT),1)
              DO I = 1,3
                 CALL READF(MARKPOINTS(NMARKPOINT,I))
              ENDDO
           ENDIF              
        CASE('MAXBENDINTER')
           CALL READI(MAXBENDINTER)
        CASE('MAXDCR')
           CALL READI(MAXDCR)
        CASE('MAXEJUMP')
           CALL READF(MAXEJUMP)
        CASE('MAXHOPTRY')
           CALL READI(MAXHOPTRY)
        CASE('MAXOPTATTEMPT')
           CALL READI(MAXOPTATTEMPT)
        CASE('MAXOPTSTEP')
           CALL READI(MAXOPTSTEP)
        CASE('MAXSTEPSIZE')
           CALL READF(MAXSTEPSIZE)
        CASE('MCSTARTRANGE')
           CALL READF(STARTRRANGE)
           CALL READF(STARTTRANGE)
           CALL READF(STARTARANGE)
        CASE('NBASINHOP')
           CALL READI(NBASINHOP)
        CASE('NCONFIGSAVE')
           CALL READI(NCONFIGSAVE)              
        CASE('NCONFIGDUMP')
           CALL READI(NCONFIGDUMP)
        CASE('NONUCSTERICS')
           BENDSTERICS = .FALSE.  
        CASE('NONUCSEGSTERICS')
           BENDSEGSTERICS = .FALSE.
        CASE('NOSEGSTERICS')
           SEGSTERICS = .FALSE.
        CASE('NSEGPERLINK')
           NSEGSET = .TRUE.
           CALL READI(NSEGPERLINK)
        CASE('NTAILSEG')
           CALL READI(NTAIL1SEG)
           IF (NITEMS.GT.2) THEN
              CALL READI(NTAIL2SEG)
           ELSE
              NTAIL2SEG = NTAIL1SEG
           ENDIF
        CASE('NUCCYLFILE')
           CALL READA(CYLINDERFILE)           
           IF (NITEMS.GT.2) THEN
              CALL READO(CYLINDERFROMFILE)
           ENDIF
        CASE('NUCSEGCYLFILE')
           CALL READA(NUCSEGCYLFILE)
           IF (NITEMS.GT.2) THEN
              CALL READO(NUCSEGCYLFROMFILE)
           ENDIF
        CASE('OPTFLEXTAILSONLY')
           OPTFLEXTAILSONLY = .TRUE.
        CASE('OPTPRINTFREQ')
           CALL READI(OPTPRINTFREQ)
        CASE('OPTOL')
           CALL READF(OPTOL)
        CASE('OUTFILE')
           CALL READA(OUTFILE)
           print*, 'outfile is:', outfile              
        CASE('POSM') 
           CALL READF(POSM(1))
           CALL READF(POSM(2))
           CALL READF(POSM(3))
        CASE('POSP') 
           CALL READF(POSP(1))
           CALL READF(POSP(2))
           CALL READF(POSP(3))
        CASE('PRINTFIBERDIAM')
           PRINTFIBERDIAM = .TRUE.
        CASE('PRINTNUCDIST')
           PRINTNUCDIST = .TRUE.
        CASE('PRINTENERGYPARTS')
           PRINTEPARTS = .TRUE.
        CASE('RANDSTART')
           RANDSTART = .TRUE.
           IF (NITEMS.GT.1) CALL READF(RSRANGE)
        CASE('READDATAFILE')
           READDATAFILE = .TRUE.
        CASE('REPLICATESTRUCT')
           REPLICATESTRUCT = .TRUE.
           CALL READI(NBENDREPLIC)
           IF (NITEMS.GT.2) CALL READA(REPLICFILE)
        CASE('RESTART')
           RESTARTFROMFILE=.TRUE.
           CALL READA(RESTARTFILE)
           IF (NITEMS.GT.2) THEN
              CALL READI(RESTARTBEADS)
           ENDIF
        CASE('RESTARTDUMP')              
           RESTARTDUMP=.TRUE.
           RESTARTFROMFILE = .TRUE.
        CASE('RNGSEED')
           CALL READI(RNGSEED)
        CASE('SAMECONF')
           CALL READF(SAMEECUT)
           CALL READF(SAMEDISTCUT)
        CASE('SEGCYLFILE')
           CALL READA(SEGCYLFILE)
           IF (NITEMS.GT.2) THEN
              CALL READO(SEGCYLFROMFILE)
           ENDIF
        CASE('SEGHOPEVERY')
           CALL READI(SEGHOPEVERY)           
        CASE('SPIRALBENDS')              
           SPIRALBENDS = .TRUE.
           IF (NITEMS.GT.1) CALL READF(SPLEN)
           IF (NITEMS.GT.2) THEN
              CALL READI(FILLINBEADS)
              FILLINBEADSET = .TRUE.
           ENDIF
           IF (NITEMS.GT.3) CALL READF(SPRADIUS)
           IF (NITEMS.GT.4) CALL READF(SPHEIGHT)
           IF (NITEMS.GT.5) THEN
              CALL READF(SPTWIST)
              SPTWISTSET = .TRUE.
           ENDIF
           IF (NITEMS.GT.6) CALL READF(SPTWIST0)              
        CASE('SPIRALHAND')
           CALL READI(SPHAND)
        CASE('STARTHELIX')
           STARTFROMHELIX = .TRUE.
           CALL READF(HELH)
           CALL READF(HELT)
           CALL READF(HELR)
           CALL READF(HELA)
           CALL READF(HELB)
           CALL READF(HELG)
        CASE('STARTHELH')
           HELCRDSET(1) = .TRUE.
           CALL READF(HELH)
        CASE('STARTHELT')
           HELCRDSET(2) = .TRUE.
           CALL READF(HELT)
        CASE('STARTHELR')
           HELCRDSET(3) = .TRUE.
           CALL READF(HELR)
        CASE('STARTHELA')
           HELCRDSET(4) = .TRUE.
           CALL READF(HELA)
        CASE('STARTHELB')
           HELCRDSET(5) = .TRUE.
           CALL READF(HELB)
        CASE('STARTHELG')
           HELCRDSET(6) = .TRUE.
           CALL READF(HELG)
        CASE('STEPDCR')
           CALL READF(STEPDCR)
        CASE('STERICSHAPE')
           CALL READF(STERRADIUS)
           IF (NITEMS.GT.2) CALL READF(STERHEIGHT)       
        CASE('TAIL')
           USEFLEXTAILS = .TRUE.
           NTAIL = NTAIL + 1
           CALL READF(TAILPOS(NTAIL,1))
           CALL READF(TAILPOS(NTAIL,2))
           CALL READF(TAILPOS(NTAIL,3))
        CASE('TAILBEAD')
           TAILNBEAD(NTAIL) = TAILNBEAD(NTAIL)+1
           DO DUMI = 1,5
              CALL READF(TAILBEADINFO(NTAIL,TAILNBEAD(NTAIL),DUMI))
           ENDDO
           IF (NITEMS.GT.6) THEN
              TAILSTARTSET = .TRUE.
              DO DUMI = 6,8
                 CALL READF(TAILBEADINFO(NTAIL,TAILNBEAD(NTAIL),DUMI))      
              ENDDO
           ENDIF
        CASE('TENSION')
           CALL READF(TENSION)
        CASE('TWIST')
           CALL READF(TWIST)
        CASE('VERBOSE')
           VERBOSE = .TRUE.
        CASE('XAXM')
           CALL READF(XAXM(1))
           CALL READF(XAXM(2))
           CALL READF(XAXM(3))
        CASE('XAXP')
           CALL READF(XAXP(1))
           CALL READF(XAXP(2))
           CALL READF(XAXP(3))           
        CASE DEFAULT
           print*, 'ERROR: unidentified keyword ', TRIM(WORD), " Will ignore."
        END SELECT
     ENDDO
     CLOSE(PF)
  ENDDO

  ! ----- set some more defaults -----
  IF (.NOT.SPTWISTSET) THEN
     SPTWIST = TWIST
  ENDIF
  IF (SPIRALBENDS.AND..NOT.FILLINBEADSET) THEN
     FILLINBEADS = SPLEN/BPLEN+1
  ENDIF
  IF (.NOT.NSEGSET) NSEGPERLINK = FLOOR(BENDSEP/BPLEN)

  ! -----------------
  ! check validity of some values and adjust as necessary
  ! -----------------  

  IF (.NOT.ABS(SPHAND).EQ.1) THEN
     PRINT*, 'ERROR IN READKEY: sphand must be 1 or -1. Current value: ', SPHAND
  ENDIF

  ! normalize the axes:
  CALL NORMALIZE(BENDTANM); CALL NORMALIZE(BENDTANP)
  CALL NORMALIZE(XAXM); CALL NORMALIZE(XAXP)  
  ! orthogonalize the axes
  IF (ABS(DOT_PRODUCT(XAXM,BENDTANM)).GT.10*EPSILON(1D0)) THEN
     !print*, 'WARNING: the chain x-axis andtangent &
     !     & are not perpendicular at the (-) bend edge. Overlap:', &
     !     & DOT_PRODUCT(XAXM,BENDTANM)
     !print*, 'Fixing this now by modifying XAXM.'
     XAXM = XAXM - DOT_PRODUCT(XAXM,BENDTANM)*BENDTANM
     CALL NORMALIZE(XAXM)
  ENDIF
  IF (ABS(DOT_PRODUCT(XAXP,BENDTANP)).GT.10*EPSILON(1D0)) THEN
     !print*, 'WARNING: the chain x-axis and chain-tangent &
     !     & are not perpendicular at (+) bend edge. Overlap: ', &
     !     & DOT_PRODUCT(XAXP,BENDTANP)
     !print*, 'Fixing this now by modifying XAXP.'
     XAXP = XAXP - DOT_PRODUCT(XAXP,BENDTANP)*BENDTANP
     CALL NORMALIZE(XAXP)
  ENDIF

  ! ----- place single filler bead if none specified --------
  IF (FILLINBEADS.EQ.0) THEN
     FILLINBEADS = 1
     ! place single filler bead half way between bend edges
     FILLBEADLIST(1,1:3) = (POSP + POSM)/2
     ! give it an arbitrary x-axis branch
     CALL COORDS2ROTMAT(0D0,POSP-FILLBEADLIST(1,1:3),ROTMAT)
     FILLBEADLIST(1,4:6) = FILLBEADLIST(1,1:3) + ROTMAT(:,1)
  ENDIF

  ! ----------- fix file names -----------
  CALL REPLACESUBSTR(OUTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(DUMPFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(RESTARTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(REPLICFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(DATAFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(NEWDATAFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(CYLINDERFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(NUCSEGCYLFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(SEGCYLFILE,'*',TRIM(ADJUSTL(ARG)))
  ! ---------------------------

  ! check if restart file exists
  IF (RESTARTFROMFILE.AND.CHECKRESTART) THEN
     INQUIRE(FILE=RESTARTFILE,EXIST=LDUM)
     IF (.NOT.LDUM) THEN
        print*, "File ", TRIM(RESTARTFILE), ' does not exist. Will not restart'
        RESTARTFROMFILE = .FALSE.
     ENDIF
  ENDIF

  ! Initiate random number generator 
  IF (RNGSEED.EQ.0) THEN
     ! use the current time of day in milliseconds
     CALL DATE_AND_TIME(VALUES=TIMEVAL)
     SEED = TIMEVAL(5)*3600*1000 + TIMEVAL(6)*60*1000 + TIMEVAL(7)*1000 + TIMEVAL(8)
  ELSEIF (RNGSEED.EQ.-1) THEN
     ! use the last 5 characters in the command-line argument
     SEED = STRING2NUM(TRIM(ADJUSTL(ARG)))    
  ELSEIF (RNGSEED.EQ.-2) THEN
     ! use the last 4 characters in the command-line argument 
     ! and additionally the millisecond time 
     CALL DATE_AND_TIME(VALUES=TIMEVAL)
     SEED = STRING2NUM(TRIM(ADJUSTL(ARG)),TIMEVAL(8))
  ELSE
     ! use this seed directly
     SEED = RNGSEED
  ENDIF

  print*, 'Initiating Mersenne twister random number generator with seed:', SEED
  CALL SGRND(SEED)

  print*, '------------Parameter values : -------------------'
  print*, 'ACTION: ', TRIM(ADJUSTL(ACTION))
  print*, 'Output file: ', TRIM(OUTFILE)  
  print*, 'Linker length, number of segments:', BENDSEP, NSEGPERLINK
  print*, 'Nucleosome geometry (+ end):'
  print '(A30,3G20.10)', 'DNA attachment position:', POSP
  print '(A30,3G20.10)', 'DNA tangent:', BENDTANP
  print '(A30,3G20.10)', 'DNA x-axis:', XAXP
  print*, 'Nucleosome geometry (- end):'
  print '(A30,3G20.10)', 'DNA attachment position:', POSM
  print '(A30,3G20.10)', 'DNA tangent:', BENDTANM
  print '(A30,3G20.10)', 'DNA x-axis:', XAXM
  PRINT*, 'DNA twist density:', TWIST
  IF (SPIRALBENDS) THEN
     PRINT'(A,4G20.10)', 'Spiral bound DNA. Length bound, radius, height, handedness:', SPLEN, SPRADIUS, SPHEIGHT, SPHAND
  ENDIF
  IF (ANY(FIXHELIXPARAM)) THEN 
     PRINT*, 'Fixed helix parameters:', FIXHELIXPARAM
  ENDIF

  IF (RESTARTFROMFILE) THEN
     PRINT*, 'Restarting structure from file: ', RESTARTFILE
     IF (RESTARTBEADS.EQ.0) THEN
        PRINT*, 'Only nucleosome positions are read from restart file.'
     ELSEIF(RESTARTBEADS.EQ.1) THEN
        PRINT*, 'Linker segments will be placed along piecewise linear path defined in restart file.'
     ELSE
        PRINT*, 'All linker segment and nucleosome positions will be read directly from restart file.'
     ENDIF
  ENDIF

  IF (RANDSTART) THEN
     PRINT*, 'Randomize starting structure. Range: ', RSRANGE
  ENDIF
  PRINT '(A,3G20.10)',' Bend, twist, stretch moduli:', LP, LTW, LSTRETCH
  PRINT*, 'External tension:', TENSION
  IF (.NOT.BENDSTERICS) PRINT*, 'Not including nucleosome-nucleosome sterics'
  IF (.NOT.SEGSTERICS) PRINT*, 'Not including linker-linker sterics'
  IF (.NOT.BENDSEGSTERICS) PRINT*, 'Not including nucleosome-linker sterics'
  print*, 'Nucleosome steric radius and height:', STERRADIUS, STERHEIGHT
  print*, 'Maximum repeat number for interactions', MAXBENDINTER
  PRINT*, 'Output file, dump file: ', TRIM(OUTFILE), ' ', TRIM(DUMPFILE)
  IF (REPLICATESTRUCT) PRINT*, 'Replicate structure to ', NBENDREPLIC, ' and output to file ', REPLICFILE

  IF (ACTION.EQ.'OPTIMIZE') THEN
     PRINT*, 'Optimization tolerance: ', OPTOL
     IF (DUMPCURRENT) PRINT*, 'Will dump out config every ', DUMPSTEPS, ' steps, with ', DUMPNBEND, ' nucleosomes.'
  ENDIF

  IF (ACTION.EQ.'BASINHOP'.OR.ACTION.EQ.'DATABASEPARSE') THEN
     PRINT*, 'Input, output data files: ', TRIM(DATAFILE),' ', TRIM(NEWDATAFILE)     
     PRINT*, 'Total number of configurations saved in database: ', NCONFIGSAVE  
     print*, 'Number of configs dumped from database: ', NCONFIGDUMP
     IF (DBCONFIGMAXE.LT.HUGE(1D0)/2) THEN
        PRINT*, 'Energy cutoff for storing configs in database: ', DBCONFIGMAXE
     ENDIF
  ENDIF

  IF (ACTION.EQ.'BASINHOP') THEN
     PRINT*, 'Number of basin hops to execute: ', NBASINHOP
     IF (READDATAFILE) PRINT*, 'Read in database file before starting'
     print*, 'Initial ranges for Monte Carlo steps: ', STARTRRANGE, STARTTRANGE, STARTARANGE
     IF (FLIPNUCMC) PRINT*, 'Flip nucleosomes on a MC hop with probability: ', FLIPNUCPROB
     PRINT*, 'Desired fraction accepted, fraction same: ', BHFACC, BHFSAME
     PRINT*, 'Adjust hop sizes every ', BHCHECKRANGE, ' steps.'
     print*, 'After each hop, minimize to tolerance: ', BHTOL
     PRINT*, 'Energy and distance cutoffs for determining identical structures: ', SAMEECUT, SAMEDISTCUT
     PRINT*, 'Do a hop that moves only linker segments every ', SEGHOPEVERY, ' steps.'
  END IF

  IF (ACTION.EQ.'DATABASEPARSE') THEN
     IF (ENERGYRECALC) PRINT*, 'Will recalculate database energies'
  ENDIF

  DO I = 1,NMARKPOINT
     PRINT '(A,1X,A,1X,3F10.5)', 'Marked point on nucleosome: ', MARKPOINTNAME(I), MARKPOINTS(I,:)
  ENDDO
  print*, '----------------------------------------------------'


END SUBROUTINE READKEY
