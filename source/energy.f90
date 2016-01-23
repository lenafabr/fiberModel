MODULE ENERGYUTILS
  ! utilities for calculating the energy associated with a chain configuration
  ! These use the single-vector regular helix representation stored in CHAINP%VEC
  USE QUATUTILS
  USE CHAINUTILS, ONLY : CHAIN, BEND

  IMPLICIT NONE
  LOGICAL :: HELENERGYTEST = .FALSE.

CONTAINS
  SUBROUTINE HELENERGYGRAD(CHAINP)
    ! get the energy and gradient for a chain (stored in chainp%energy, chainp%grad)
    ! CHAINP%ENERGY and CHAINP%GRAD are set on exit
    ! CHAINP%ENERGYPARTS will contain on exit:
    ! (bend, twist, stretch, nucleosome sterics, segment sterics) if not using any internucleosomal potential
    ! (bend, twist, stretch, langowski energy, segment sterics) if using langowski internucleosmal potential
    ! (bend, twist, stretch, sterics, electrostatics) if using DiSCO potential
    USE KEYs, ONLY : OPTFLEXTAILSONLY, LANGOWSKIPOT, DISCOPOTENTIAL, BPLEN,USEFIXDIAM

    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION :: PBOT(3),MATBOT(3,3),DPBOT(3,6),DMATBOT(3,3,6)
    DOUBLE PRECISION :: PTOP(3),MATTOP(3,3),DPTOP(3,6),DMATTOP(3,3,6)
    DOUBLE PRECISION :: MAT1(3,3),MAT2(3,3),DMAT1(3,3,14),DMAT2(3,3,14)
    DOUBLE PRECISION :: ZVEC(3), CBETA, ALPHAGAMMA,DCBETA(14),DAG(14), PHI
    DOUBLE PRECISION :: DIST, DDIST(10)
    DOUBLE PRECISION :: LPLS, LTLS, LSTLS
    DOUBLE PRECISION :: BENDE, TWISTE, STRETCHE, ELECTE, LJE, TAILE, NUCSTERE, SEGSTERE, NUCSEGSTERE
    DOUBLE PRECISION :: GRADNUCSTER(CHAINP%NCRDDNA), GRADSEGSTER(CHAINP%NCRDDNA), GRADNUCSEGSTER(CHAINP%NCRDDNA)
    DOUBLE PRECISION :: GRADBND(CHAINP%NCRDDNA),GRADST(CHAINP%NCRDDNA),GRADTW(CHAINP%NCRDDNA)
    DOUBLE PRECISION :: ELECTGRAD(CHAINP%NCRD), LJGRAD(CHAINP%NCRD), TAILGRAD(CHAINP%NCRD)
    DOUBLE PRECISION :: RVAL, DRVAL(6)
    INTEGER :: I, J, CRDI, IND
    DOUBLE PRECISION :: TMP, DEDR

    IF (OPTFLEXTAILSONLY.AND.DISCOPOTENTIAL) THEN
       ! only get the energy and gradient components for the flexible tails
       PRINT*, 'ERROR IN ENERGYUTILS: flexible tail energies only not set up yet'
       STOP 1
    ENDIF

    IF (USEFIXDIAM) THEN
       ! treat helix radius R as a dependent variable of the other helix coords
       CALL SETFIBERDIAM(CHAINP%VEC(1:6),CHAINP%DIAM, CHAINP%BENDS(1)%RADIUS, RVAL,DRVAL)
       CHAINP%VEC(3) = RVAL       
    ENDIF

    ! Scale the persistence lengths
    LPLS = CHAINP%LP/CHAINP%LS; LTLS = CHAINP%LTW/CHAINP%LS/2;
    LSTLS = CHAINP%LSTRETCH/CHAINP%LS/2

    ! Initialize energies and grads
    BENDE = 0D0; TWISTE = 0D0; STRETCHE = 0D0; 
    GRADBND = 0D0; GRADTW = 0D0; GRADST = 0D0; 
    ELECTE = 0D0; ELECTGRAD = 0D0
    LJE = 0D0; LJGRAD = 0D0; TAILE = 0D0; TAILGRAD = 0D0
    NUCSTERE = 0D0; SEGSTERE = 0D0; NUCSEGSTERE = 0D0
    GRADNUCSTER = 0D0; GRADSEGSTER = 0D0; GRADNUCSEGSTER = 0D0
    CHAINP%GRAD = 0D0

    ! position and orientation of linker start and end
    DMATTOP = 0D0; DMATBOT = 0D0
    CALL GETTOPEDGE(CHAINP,PTOP,MATTOP,DPTOP,DMATTOP)
    CALL GETBOTTOMEDGE(CHAINP,PBOT,MATBOT,DPBOT,DMATBOT) 

    DMAT1 = 0D0
    MAT1 = MATTOP; DMAT1(:,:,1:6) = DMATTOP

    ! ---------------------------------
    ! top edge of 1st bend 
    ! --------------------------------
    ! get rotation matrix for 1st segment
    ZVEC = CHAINP%VEC(8:10)-PTOP
    DMAT2 = 0D0
    CALL COORDS2ROTMAT(CHAINP%VEC(7),ZVEC,MAT2,DMAT2(:,:,7),DMAT2(:,:,8:10))
    ! derivatives wrt helix coordinates
    DO J = 1,3
       CALL DGEMM('N','N',3,6,3,-1D0,DMAT2(:,J,8:10),3,DPTOP,3,0D0,DMAT2(:,J,1:6),3)
    ENDDO

    ! get relative rotation matrix
    CALL ROTMAT2RELEULER(MAT1,MAT2,10,DMAT1(:,:,1:10),DMAT2(:,:,1:10),CBETA,ALPHAGAMMA,DCBETA(1:10),DAG(1:10))

    ! bend energy
    BENDE = LPLS*(1D0-CBETA)    
    GRADBND(1:10) = -LPLS*DCBETA(1:10)
    ! twist energy
    PHI = ALPHAGAMMA - CHAINP%TWIST*(CHAINP%LS/2+BPLEN/2)
    PHI = ANGLE0(PHI)
    !print*, 'testxa:', phi, alphagamma

    TWISTE =  LTLS*PHI**2     
    GRADTW(1:10) = 2*LTLS*PHI*DAG(1:10)
    ! stretch energy
    DIST = SQRT(DOT_PRODUCT(ZVEC,ZVEC))
    DDIST = 0D0
    DDIST(8:10) = ZVEC/DIST
    CALL DGEMV('T',3,6,-1D0,dPTOP,3,DDIST(8:10),1,0D0,DDIST(1:6),1)
    STRETCHE = LSTLS*(DIST-CHAINP%LS)**2
    GRADST(1:10) = 2*LSTLS*(DIST-CHAINP%LS)*DDIST    

    ! ----- energy at each joint along the linker -------------
    CRDI = 1
    DO I = 1,CHAINP%NSEG-1
       MAT1 = MAT2
       DMAT1 = 0D0
       DMAT1(:,:,1:14-CRDI+1) = DMAT2(:,:,CRDI:14)

       IF (I.EQ.1) THEN
          ! index at which coordinates affecting 2nd segment start
          ! in the temporary derivative vectors
          CRDI = 8
       ELSE
          CRDI = 5
       ENDIF

       IND = 7+4*(I-1)

       ! rotation matrix for next segment    
       DMAT2 = 0D0       
       IF (I.EQ.CHAINP%NSEG-1) THEN
          ZVEC = PBOT-CHAINP%VEC(IND+1:IND+3)
          CALL COORDS2ROTMAT(CHAINP%VEC(ind+4),ZVEC,MAT2,DMAT2(:,:,CRDI+3),DMAT2(:,:,CRDI:CRDI+2))
          DMAT2(:,:,CRDI:CRDI+2) = -DMAT2(:,:,CRDI:CRDI+2)
          ! derivatives wrt helix coordinates
          DO J = 1,3
             CALL DGEMM('N','N',3,6,3,-1D0,DMAT2(:,J,CRDI:CRDI+2),3,DPBOT,3,0D0,DMAT2(:,J,CRDI+4:CRDI+9),3)
          ENDDO

       ELSE
          ZVEC = CHAINP%VEC(IND+5:IND+7)-CHAINP%VEC(IND+1:IND+3)          

          CALL COORDS2ROTMAT(CHAINP%VEC(ind+4),ZVEC,MAT2,DMAT2(:,:,CRDI+3),DMAT2(:,:,CRDI+4:CRDI+6))          

          DMAT2(:,:,CRDI:CRDI+2) = -DMAT2(:,:,CRDI+4:CRDI+6)

       ENDIF

       ! get relative rotation matrix
       IF (I.EQ.1.OR.I.EQ.CHAINP%NSEG-1) THEN         
          CALL ROTMAT2RELEULER(MAT1,MAT2,14,DMAT1,DMAT2,CBETA,ALPHAGAMMA,DCBETA,DAG)
       ELSE                    
          CALL ROTMAT2RELEULER(MAT1,MAT2,11,DMAT1(:,:,1:11),DMAT2(:,:,1:11),&
               & CBETA,ALPHAGAMMA,DCBETA(1:11),DAG(1:11))
       ENDIF

       ! bend energy      
       BENDE = BENDE + LPLS*(1D0-CBETA)    

       ! twist energy
       PHI = ALPHAGAMMA - CHAINP%TWIST*CHAINP%LS
       PHI = ANGLE0(PHI)
       TWISTE =  TWISTE + LTLS*PHI**2           
      
       IF (I.EQ.1) THEN
          GRADBND(1:14) = GRADBND(1:14)-LPLS*DCBETA          
          GRADTW(1:14) = GRADTW(1:14) + 2*LTLS*PHI*DAG
       ELSEIF(I.EQ.CHAINP%NSEG-1) THEN
          GRADBND(1:6) = GRADBND(1:6) -LPLS*DCBETA(9:14)
          GRADBND(IND-3:IND+4) = GRADBND(IND-3:IND+4) -LPLS*DCBETA(1:8)
          GRADTW(1:6) = GRADTW(1:6) + 2*LTLS*PHI*DAG(9:14)
          GRADTW(IND-3:IND+4) = GRADTW(IND-3:IND+4) + 2*LTLS*PHI*DAG(1:8)
       ELSE
          GRADBND(IND-3:IND+7) = GRADBND(IND-3:IND+7)-LPLS*DCBETA(1:11)
          GRADTW(IND-3:IND+7) = GRADTW(IND-3:IND+7) + 2*LTLS*PHI*DAG(1:11)
       ENDIF       

       ! stretch energy
       DIST = SQRT(DOT_PRODUCT(ZVEC,ZVEC))
       TMP = DIST-CHAINP%LS
       STRETCHE = STRETCHE + LSTLS*TMP**2       

!       print*, 'testxb:', i,phi,alphagamma/dist
        
!      print*, 'testxa:', i, ind, cbeta, phi,dist

       IF (I.EQ.CHAINP%NSEG-1) THEN
          DDIST = 0D0; DDIST(1:3) = -ZVEC/DIST
          CALL DGEMV('T',3,6,-1D0,dPBOT,3,DDIST(1:3),1,0D0,DDIST(4:9),1)
          GRADST(ind+1:ind+3) = GRADST(ind+1:ind+3)+2*LSTLS*TMP*DDIST(1:3)
          GRADST(1:6) = GRADST(1:6) + 2*LSTLS*TMP*DDIST(4:9)
       ELSE
          GRADST(ind+5:ind+7) = GRADST(ind+5:ind+7)+2*LSTLS*TMP*ZVEC/DIST
          GRADST(IND+1:IND+3) = GRADST(IND+1:IND+3)-2*LSTLS*TMP*ZVEC/DIST
       ENDIF
    ENDDO

    ! ---------------------------------
    ! bottom edge of 2nd bend 
    ! --------------------------------
    MAT1 = MAT2
    DMAT1(:,:,1:10) = DMAT2(:,:,5:14)

    MAT2 = MATBOT
    DMAT2 = 0D0; DMAT2(:,:,5:10) = DMATBOT
    ! get relative rotation matrix
    CALL ROTMAT2RELEULER(MAT1,MAT2,10,DMAT1(:,:,1:10),DMAT2(:,:,1:10),CBETA,ALPHAGAMMA,DCBETA(1:10),DAG(1:10))        

    ! bend energy
    BENDE = BENDE+LPLS*(1D0-CBETA)    
    GRADBND(1:6) = GRADBND(1:6)-LPLS*DCBETA(5:10)
    GRADBND(4*CHAINP%NSEG:4*CHAINP%NSEG+3) = GRADBND(4*CHAINP%NSEG:4*CHAINP%NSEG+3)-LPLS*DCBETA(1:4)
    ! twist energy
    PHI = ALPHAGAMMA - CHAINP%TWIST*(CHAINP%LS/2+BPLEN/2)
    PHI = ANGLE0(PHI)
    TWISTE =  TWISTE + LTLS*PHI**2     
    GRADTW(1:6) = GRADTW(1:6) + 2*LTLS*PHI*DAG(5:10)
    GRADTW(4*CHAINP%NSEG:4*CHAINP%NSEG+3) = GRADTW(4*CHAINP%NSEG:4*CHAINP%NSEG+3) + 2*LTLS*PHI*DAG(1:4)           

    !print*, 'testxc:', i, phi

    ! Steric energies
    IF (CHAINP%BENDSTERICS.OR.CHAINP%BENDSEGSTERICS.OR.CHAINP%SEGSTERICS) THEN
       CALL GETSTERICENERGY(CHAINP,PTOP,PBOT,DPTOP,DPBOT,NUCSTERE,SEGSTERE,NUCSEGSTERE,&
            & GRADNUCSTER,GRADSEGSTER,GRADNUCSEGSTER)
    ENDIF
    

    IF (LANGOWSKIPOT) THEN
       PRINT*, 'ERROR IN HELENERGYGRAD: langowski potential not yet set up'
       stop 1
    ELSEIF (DISCOPOTENTIAL) THEN
       PRINT*, 'ERROR IN HELENERGYGRAD: disco potential not yet set up'
       STOP 1
    ENDIF
   
    ! --------------------------------
    CHAINP%ENERGY = BENDE + TWISTE + STRETCHE+NUCSTERE+SEGSTERE+NUCSEGSTERE
    CHAINP%GRAD(1:CHAINP%NCRDDNA) = GRADBND + GRADTW + GRADST+GRADNUCSTER+GRADSEGSTER+GRADNUCSEGSTER

    IF (USEFIXDIAM) THEN
       DEDR = CHAINP%GRAD(3)
       CHAINP%GRAD(1:6) = CHAINP%GRAD(1:6) + DEDR*DRVAL
       CHAINP%GRAD(3) = 0D0
    ENDIF

    CHAINP%ENERGYPARTS = (/BENDE,TWISTE,STRETCHE,NUCSTERE,SEGSTERE,NUCSEGSTERE/)

    ! add on external tension
    IF (CHAINP%TENSION.NE.0) THEN
       CHAINP%ENERGY = CHAINP%ENERGY - CHAINP%TENSION*CHAINP%VEC(1)
       CHAINP%GRAD(1) = CHAINP%GRAD(1) - CHAINP%TENSION
    ENDIF

  END SUBROUTINE HELENERGYGRAD

  
  SUBROUTINE GETSTERICENERGY(CHAINP,PTOP,PBOT,DPTOP,DPBOT,&
       & NUCNUCE,SEGSEGE,NUCSEGE,NUCNUCGRAD,SEGSEGGRAD,NUCSEGGRAD)
    ! get the steric energies for the chain
    ! NUCNUCE has steric energies between nucleosomes
    ! SEGSEGE has steric energies between DNA segments
    ! NUCSEGE has steric energies btwn DNA and nucleosomes    
    ! GRAD arrays have the derivatives

    USE CYLINDERUTILS, ONLY : SEGCYL, NUCCYL, NUCSEGCYL, SEGCYLSET, NUCCYLSET, NUCSEGCYLSET
    USE TRICUBSPLINEUTILS, ONLY : EVALUATESPLINE
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(IN) :: PTOP(3), PBOT(3), DPTOP(3,6), DPBOT(3,6)
    DOUBLE PRECISION, INTENT(OUT) :: NUCNUCE, SEGSEGE, NUCSEGE
    DOUBLE PRECISION, INTENT(OUT) :: NUCNUCGRAD(CHAINP%NCRDDNA), SEGSEGGRAD(CHAINP%NCRDDNA), &
         & NUCSEGGRAD(CHAINP%NCRDDNA)
    DOUBLE PRECISION :: OVERLAP, PBN(CHAINP%NSEG+1,3),DPBN(CHAINP%NSEG+1,3,6)
    DOUBLE PRECISION :: IJ(3), KJ(3), LK(3), KJN(3), CA, SA, CB, SB, SNTA, CNTA
    DOUBLE PRECISION :: H, R, THETA, CNT, SNT, NTA
    DOUBLE PRECISION :: CT1, CT2, PHI,DIST,MINDIST
    DOUBLE PRECISION, DIMENSION(CHAINP%NCRDDNA) :: DCT1, DCT2, DPHI, DDIST, DMIN
    DOUBLE PRECISION :: DCT1DIJ(3),DCT1DKJ(3), DCT2DJK(3),DCT2DLK(3)
    DOUBLE PRECISION :: DPDIJ(3),DPDJK(3),DPDLK(3)
    DOUBLE PRECISION :: NUC0Z(3),NUC0CENT(3),NUCZ(3),NUCCENT(3)
    DOUBLE PRECISION :: DNUC0Z(3,6),DNUC0CENT(3,6),DNUCZ(3,6),DNUCCENT(3,6), DKJ(3,6)
    DOUBLE PRECISION :: QVROT(4), DQROT(3,4), DPNQ(3,4),DPNP(3,3)
    DOUBLE PRECISION :: DMDC1,DMDC2,DMDP,DIFF(3)
    INTEGER :: MINSEG, S1, S2, I, NNUC, NSEG
    DOUBLE PRECISION :: NNMAXSEP, NSMAXSEP, SSMAXSEP

    NSEG = CHAINP%NSEG

    NUCNUCE = 0D0; SEGSEGE = 0D0; NUCSEGE = 0D0; 
    NUCNUCGRAD = 0D0; SEGSEGGRAD = 0D0; NUCSEGGRAD = 0D0;
    
    ! The first MINSEG segments attached to a nucleosome do not interact with that nucleosome
    MINSEG = INT(2*CHAINP%DNARAD/CHAINP%LS) + 2

    IF (chainp%NSEG.EQ.1) THEN
       PRINT*, 'ERROR IN getstericenergy: CURRENTLY NOT SET UP TO DEAL WITH LINKERS THAT HAVE ONLY 1 SEGMENT'
       STOP 1
    ENDIF

    IF ((CHAINP%SEGSTERICS.AND..NOT.SEGCYLSET).OR.&
         & (CHAINP%BENDSTERICS.AND..NOT.NUCCYLSET).OR.&
         & (CHAINP%BENDSEGSTERICS.AND..NOT.NUCSEGCYLSET)) THEN
       PRINT*, 'ERROR IN SEGMENTOVERLAP: cylindrical steric data not set up', NUCCYLSET, SEGCYLSET
       PRINT*, CHAINP%BENDSTERICS, CHAINP%SEGSTERICS, CHAINP%BENDSEGSTERICS
       STOP 1
    ENDIF

    ! maximal separation for nuc-nuc, nuc-seg, seg-seg at which cylinders can still overlap
    ! plus a little bit of an offset to avoid numerical issues 
    ! due to the interpolation of the cylinder tables
    NNMAXSEP = SQRT(CHAINP%BENDS(1)%RADIUS**2 + (CHAINP%BENDS(1)%HEIGHT/2)**2)*2*1.1D0
    SSMAXSEP = SQRT(CHAINP%DNARAD**2 + (CHAINP%LS/2)**2)*2*1.1D0
    NSMAXSEP = (NNMAXSEP+SSMAXSEP)/2


    IF (CHAINP%BENDSTERICS.OR.CHAINP%BENDSEGSTERICS) THEN
       ! get position and z axis of 0th nucleosome and the derivatives
       H = CHAINP%VEC(1); THETA = CHAINP%VEC(2); R = CHAINP%VEC(3); 
       NUC0CENT = (/R,0D0,0D0/)
       DNUC0CENT = 0D0; DNUC0CENT(1,3) = 1D0
       CA = COS(CHAINP%VEC(4)); SA = SIN(CHAINP%VEC(4)); CB = COS(CHAINP%VEC(5)); SB = SIN(CHAINP%VEC(5))
       NUC0Z = (/SB*SA,-SB*CA,CB/)             
       DNUC0Z = 0D0
       DNUC0Z(:,4) = (/SB*CA,SB*SA,0D0/)
       DNUC0Z(:,5) = (/CB*SA,-CB*CA,-SB/)
    ENDIF

    DO NNUC = 0,CHAINP%NINTERACT
       IF (CHAINP%BENDSTERICS.OR.CHAINP%BENDSEGSTERICS) THEN          
          IF (NNUC.EQ.0) THEN
             NUCCENT = NUC0CENT
             DNUCCENT = DNUC0CENT
             NUCZ = NUC0Z
             DNUCZ = DNUC0Z
          ELSE
             ! get position and z axis of Nth nucleosome and the derivatives
             H = CHAINP%VEC(1); THETA = CHAINP%VEC(2); R = CHAINP%VEC(3); 
             CNT = COS(NNUC*THETA); SNT = SIN(NNUC*THETA)
             NUCCENT = (/R*CNT,R*SNT,H*NNUC/)

             DNUCCENT = 0D0
             DNUCCENT(:,1) = (/0D0,0D0,dble(NNUC)/)
             DNUCCENT(:,2) = (/-R*NNUC*SNT, R*NNUC*CNT, 0D0/)
             DNUCCENT(:,3) = (/CNT, SNT, 0D0/)

             NTA = NNUC*CHAINP%VEC(2) + CHAINP%VEC(4); SNTA = SIN(NTA); CNTA = COS(NTA)

             NUCZ = (/SB*SNTA,-SB*CNTA,CB/)             
             DNUCZ = 0D0
             DNUCZ(:,2) = (/SB*CNTA*NNUC,SB*SNTA*NNUC,0D0/)
             DNUCZ(:,4) = (/SB*CNTA,SB*SNTA,0D0/)
             DNUCZ(:,5) = (/CB*SNTA,-CB*CNTA,-SB/)    


             IF (CHAINP%BENDSTERICS) THEN
                ! Interaction btwn 0th nucleosome and n-th nucleosome
                KJ = NUC0CENT - NUCCENT
                DKJ = DNUC0CENT-DNUCCENT

                DIST = SQRT(DOT_PRODUCT(KJ,KJ)); KJN = KJ/DIST                

                DDIST = 0D0; DCT1 = 0D0; DCT2 = 0D0; DPHI = 0D0
                IF (DIST.LT.NNMAXSEP) THEN
                   CALL GETANGLE(NUC0Z,KJ,CT1,DCT1DIJ,DCT1DKJ)
                   CALL GETANGLE(-KJ,NUCZ,CT2,DCT2DJK,DCT2DLK)
                   CALL GETDIHEDRAL(NUC0Z,-KJ,NUCZ,PHI,DPDIJ,DPDJK,DPDLK)

                   CALL DGEMV('T',3,6,1D0,DKJ,3,KJN,1,0D0,DDIST(1:6),1)   
                   CALL DGEMV('T',3,6,1D0,DNUC0Z,3,DCT1DIJ,1,0D0, DCT1(1:6),1)
                   CALL DGEMV('T',3,6,1D0,DKJ,3,DCT1DKJ,1,1D0, DCT1(1:6),1)
                   CALL DGEMV('T',3,6,-1D0,DKJ,3,DCT2DJK,1,0D0,DCT2(1:6),1) 
                   CALL DGEMV('T',3,6,1D0,DNUCZ,3,DCT2DLK,1,1D0,DCT2(1:6),1) 
                   CALL DGEMV('T',3,6,1D0,DNUC0Z,3,DPDIJ,1,1D0,DPHI(1:6),1)
                   CALL DGEMV('T',3,6,-1D0,DKJ,3,DPDJK,1,1D0,DPHI(1:6),1)
                   CALL DGEMV('T',3,6,1D0,DNUCZ,3,DPDLK,1,1D0,DPHI(1:6),1)

                   CALL EVALUATESPLINE(NUCCYL,CT1,CT2,PHI,MINDIST,DMDC1,DMDC2,DMDP)
                   IF (.NOT.MINDIST.GE.0) THEN
                      PRINT*, 'ERROR IN GETSTERICENERGY: bad mindist between nucleosomes', MINDIST, CT1, CT2, PHI, nnuc
                      STOP 1
                   ENDIF
                   
                   IF (DIST.LT.MINDIST) THEN
                      DMIN = DMDC1*DCT1 + DMDC2*DCT2 + DMDP*DPHI 
                      OVERLAP =  MINDIST-DIST
                      NUCNUCE = NUCNUCE + CHAINP%ESTERIC*OVERLAP**2
                      NUCNUCGRAD = NUCNUCGRAD +  2* CHAINP%ESTERIC*OVERLAP*(DMIN-DDIST)
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
       ENDIF


       IF (CHAINP%SEGSTERICS.OR.CHAINP%BENDSEGSTERICS) THEN
          PBN = 0D0; DPBN = 0D0; 
          ! get positions, orientations, and derivatives for the segments in the n-th linker
          ! rotation vector for rotation around helix axis
          CALL ROTQV(NNUC*CHAINP%VEC(2),(/0D0,0D0,1D0/),QVROT,DQROT)

          ! first bead in linker (top edge of 1st bead)
          ! derivatives: 1-6 wrt helix coordinates
          CALL QVPTMULT(QVROT,PTOP,PBN(1,:),DPNQ,DPNP)
          CALL DGEMM('N','N',3,6,3,1D0,DPNP,3,DPTOP,3,0D0,DPBN(1,:,:),3) 
          ! add on the derivative wrt theta
          CALL DGEMV('N',3,4,dble(nnuc),DPNQ,3,DQROT,1,1D0,dPBN(1,:,2),1) 

          ! internal beads 
          ! derivatives: 1 wrt helix height, 2 wrt helix angle 3-5 wrt bead coords
          DO I = 2,NSEG
             ! rotate the bead around the fiber axis
             CALL QVPTMULT(QVROT,CHAINP%VEC(4+4*(I-1):6+4*(I-1)),PBN(I,:),DPNQ,DPBN(I,:,3:5))

             ! get the derivatives wrt theta
             CALL DGEMV('N',3,4,1D0,DPNQ,3,DQROT,1,0D0,dPBN(I,:,2),1)  
             dPBN(I,:,2) = dPBN(I,:,2)*NNUC             
          ENDDO

          ! last bead in linker (bottom edge of 2nd bead)
          I = NSEG+1
          CALL QVPTMULT(QVROT,PBOT,PBN(I,:),DPNQ,DPNP)
          CALL DGEMM('N','N',3,6,3,1D0,DPNP,3,DPBOT,3,0D0,DPBN(I,:,:),3) 
          ! add on the derivative wrt theta
          CALL DGEMV('N',3,4,dble(NNUC),DPNQ,3,DQROT,1,1D0,dPBN(I,:,2),1) 

          ! shift beads upward and get derivative wrt height
          DO I = 1,NSEG+1
             PBN(I,3) = PBN(I,3) + NNUC*CHAINP%VEC(1)       
             dPBN(I,3,1) = dPBN(I,3,1) + NNUC  
          ENDDO
       ENDIF

       
       ! now go through and get the overlap for each segment after the 0th nucleosome (with other linkers and  nucleosomes)
       DO S1 = 1,NSEG      

          IF (S1.EQ.1) THEN  
             IJ = (PTOP - CHAINP%VEC(4+4*S1:6+4*S1) )/2          
          ELSE IF (S1.EQ.NSEG) THEN
             IJ = (CHAINP%VEC(4+4*(S1-1):6+4*(S1-1)) - PBOT )/2         
          ELSE
             IJ = (CHAINP%VEC(4+4*(S1-1):6+4*(S1-1)) - CHAINP%VEC(4+4*S1:6+4*S1) )/2         
          ENDIF


          ! S1 segment after 0th nucleosome overlap with NNUC nucleosome
          IF (CHAINP%BENDSEGSTERICS.AND..NOT.&
               & ((NNUC.EQ.0.AND.S1.LE.MINSEG).OR.(NNUC.EQ.1.AND.S1.GE.NSEG-MINSEG+1))) THEN

             DDIST = 0D0; DCT1 = 0D0; DCT2 = 0D0; DPHI = 0D0     
             IF (S1.EQ.1) THEN               
                KJ = NUCCENT - ( CHAINP%VEC(4+4*S1:6+4*S1) + PTOP)/2
             ELSE IF (S1.EQ.NSEG) THEN             
                KJ = NUCCENT - (PBOT + CHAINP%VEC(4+4*(S1-1):6+4*(S1-1)))/2
             ELSE             
                KJ = NUCCENT - (CHAINP%VEC(4+4*S1:6+4*S1) + CHAINP%VEC(4+4*(S1-1):6+4*(S1-1)))/2
             ENDIF

             DIST = SQRT(DOT_PRODUCT(KJ,KJ)); KJN = KJ/DIST*0.5D0             

             IF (DIST.LE.NSMAXSEP) THEN
                DCT1 = 0D0; DCT2 = 0D0; DPHI = 0D0; DDIST = 0D0

                ! get angles and dihedrals
                CALL GETANGLE(IJ,KJ,CT1,DCT1DIJ,DCT1DKJ)
                CALL GETANGLE(-KJ,NUCZ,CT2,DCT2DJK,DCT2DLK)
                CALL GETDIHEDRAL(IJ,-KJ,NUCZ,PHI,DPDIJ,DPDJK,DPDLK)
                ! derivatives wrt endpoints of linker segment
                IF (S1.EQ.1) THEN
                   CALL DGEMV('T',3,6,-1D0,DPTOP,3,KJN,1,1D0,DDIST(1:6),1)   
                   CALL DGEMV('T',3,6,0.5D0,DPTOP,3,DCT1DIJ-DCT1DKJ,1,1D0, DCT1(1:6),1)
                   CALL DGEMV('T',3,6,0.5D0,DPTOP,3,DCT2DJK,1,1D0,DCT2(1:6),1) 
                   CALL DGEMV('T',3,6,0.5D0,DPTOP,3,DPDIJ+DPDJK,1,1D0,Dphi(1:6),1)
                ELSE
                   DDIST(4+4*(S1-1):6+4*(S1-1)) = DDIST(4+4*(S1-1):6+4*(S1-1)) - KJN
                   DCT1(4+4*(S1-1):6+4*(S1-1)) = DCT1(4+4*(S1-1):6+4*(S1-1)) + 0.5D0*(DCT1DIJ-DCT1DKJ)
                   DCT2(4+4*(S1-1):6+4*(S1-1)) = DCT2(4+4*(S1-1):6+4*(S1-1)) + 0.5D0*DCT2DJK
                   DPHI(4+4*(S1-1):6+4*(S1-1)) = DPHI(4+4*(S1-1):6+4*(S1-1)) + 0.5D0*(DPDIJ+DPDJK)
                ENDIF

                IF (S1.EQ.NSEG) THEN
                   CALL DGEMV('T',3,6,-1D0,DPBOT,3,KJN,1,1D0,DDIST(1:6),1)   
                   CALL DGEMV('T',3,6,0.5D0,DPBOT,3,-DCT1DIJ-DCT1DKJ,1,1D0, DCT1(1:6),1)
                   CALL DGEMV('T',3,6,0.5D0,DPBOT,3,DCT2DJK,1,1D0,DCT2(1:6),1) 
                   CALL DGEMV('T',3,6,0.5D0,DPBOT,3,-DPDIJ+DPDJK,1,1D0,Dphi(1:6),1)
                ELSE
                   DDIST(4+4*S1:6+4*S1) = DDIST(4+4*S1:6+4*S1) - KJN             
                   DCT1(4+4*S1:6+4*S1) = DCT1(4+4*S1:6+4*S1) + 0.5D0*(-DCT1DIJ - DCT1DKJ)             
                   DCT2(4+4*S1:6+4*S1) = DCT2(4+4*S1:6+4*S1) + 0.5D0*DCT2DJK             
                   DPHI(4+4*S1:6+4*S1) = DPHI(4+4*S1:6+4*S1) + 0.5D0*(-DPDIJ +DPDJK)             
                ENDIF

                ! derivatives wrt Nth nucleosome position and orientation
                CALL DGEMV('T',3,6,2D0,DNUCCENT,3,KJN,1,1D0,DDIST(1:6),1)   
                CALL DGEMV('T',3,6,1D0,DNUCCENT,3,DCT1DKJ,1,1D0, DCT1(1:6),1)
                CALL DGEMV('T',3,6,-1D0,DNUCCENT,3,DCT2DJK,1,1D0,DCT2(1:6),1) 
                CALL DGEMV('T',3,6,-1D0,DNUCCENT,3,DPDJK,1,1D0,Dphi(1:6),1)

                CALL DGEMV('T',3,6,1D0,DNUCZ,3,DCT2DLK,1,1D0,DCT2(1:6),1) 
                CALL DGEMV('T',3,6,1D0,DNUCZ,3,DPDLK,1,1D0,Dphi(1:6),1)

                CALL EVALUATESPLINE(NUCSEGCYL,CT1,CT2,PHI,MINDIST,DMDC1,DMDC2,DMDP)
                IF (.NOT.MINDIST.GE.0) THEN
                   PRINT*, 'ERROR IN calculating overlap between segment &
                        & after 0th nucleosome and Nth nucleosome', MINDIST, CT1, CT2, PHI, S1
                   STOP 1
                ENDIF

                DMIN = DMDC1*DCT1 + DMDC2*DCT2 + DMDP*DPHI          
                IF (DIST.LT.MINDIST) THEN
                   OVERLAP =  MINDIST-DIST
                   NUCSEGE = NUCSEGE + CHAINP%ESTERIC*OVERLAP**2
                   NUCSEGGRAD = NUCSEGGRAD +  2* CHAINP%ESTERIC*OVERLAP*(DMIN-DDIST)
                ENDIF
             ENDIF
          ENDIF

          IF (NNUC.EQ.0) CYCLE

          IF (CHAINP%BENDSEGSTERICS) THEN
             DDIST = 0D0; DCT1 = 0D0; DCT2 = 0D0; DPHI = 0D0     

             ! get interactions between linkers after n-th nucleosome and 1st nucleosome
             KJ = (PBN(S1+1,:) + PBN(S1,:))/2 - NUC0CENT
             LK = (PBN(S1+1,:) - PBN(S1,:))/2

             ! get the actual separation distance                
             DIST = SQRT(DOT_PRODUCT(KJ,KJ)); KJN = KJ/DIST*0.5D0             

             IF (DIST.LE.NSMAXSEP) THEN
                ! get angles and dihedrals
                CALL GETANGLE(NUC0Z,KJ,CT1,DCT1DIJ,DCT1DKJ)
                CALL GETANGLE(-KJ,LK,CT2,DCT2DJK,DCT2DLK)
                CALL GETDIHEDRAL(NUC0Z,-KJ,LK,PHI,DPDIJ,DPDJK,DPDLK)

                ! derivatives based on nucleosome position and orientation
                CALL DGEMV('T',3,6,-2D0,DNUC0CENT,3,KJN,1,1D0,DDIST(1:6),1)   
                CALL DGEMV('T',3,6,1D0,DNUC0CENT,3,-DCT1DKJ,1,1D0, DCT1(1:6),1)
                CALL DGEMV('T',3,6,1D0,DNUC0CENT,3,DCT2DJK,1,1D0,DCT2(1:6),1)
                CALL DGEMV('T',3,6,1D0,DNUC0CENT,3,DPDJK,1,1D0,Dphi(1:6),1)

                CALL DGEMV('T',3,6,1d0,DNUC0Z,3,DCT1DIJ,1,1D0, DCT1(1:6),1)
                CALL DGEMV('T',3,6,1D0,DNUC0Z,3,DPDIJ,1,1D0,Dphi(1:6),1)          

                ! derivatives wrt endpoints of 2nd segment
                ! derivatives based on bottom bead
                IF (S1.EQ.1) THEN
                   CALL DGEMV('T',3,6,1D0,DPBN(S1,:,:),3,KJN,1,1D0,DDIST(1:6),1)   
                   CALL DGEMV('T',3,6,0.5D0,DPBN(S1,:,:),3,DCT1DKJ,1,1D0, DCT1(1:6),1)
                   CALL DGEMV('T',3,6,0.5D0,DPBN(S1,:,:),3,-DCT2DLK-DCT2DJK,1,1D0,DCT2(1:6),1) 
                   CALL DGEMV('T',3,6,0.5D0,DPBN(S1,:,:),3,-DPDJK-DPDLK,1,1D0,Dphi(1:6),1) 
                ELSE
                   CALL DGEMV('T',3,3,1D0,DPBN(S1,:,3:5),3,KJN,1,1D0,DDIST(4+4*(S1-1):6+4*(S1-1)),1) 
                   CALL DGEMV('T',3,2,1D0,DPBN(S1,:,1:2),3,KJN,1,1D0, DDIST(1:2),1)

                   CALL DGEMV('T',3,3,0.5D0,DPBN(S1,:,3:5),3,DCT1DKJ,1,1D0,DCT1(4+4*(S1-1):6+4*(S1-1)),1) 
                   CALL DGEMV('T',3,2,0.5D0,DPBN(S1,:,1:2),3,DCT1DKJ,1,1D0, DCT1(1:2),1)

                   CALL DGEMV('T',3,3,0.5D0,DPBN(S1,:,3:5),3,-DCT2DLK-DCT2DJK,1,1D0,DCT2(4+4*(S1-1):6+4*(S1-1)),1)
                   CALL DGEMV('T',3,2,0.5D0,DPBN(S1,:,1:2),3,-DCT2DLK-DCT2DJK,1,1D0,DCT2(1:2),1) 

                   CALL DGEMV('T',3,3,0.5D0,DPBN(S1,:,3:5),3,-DPDJK-DPDLK,1,1D0,DPHI(4+4*(S1-1):6+4*(S1-1)),1) 
                   CALL DGEMV('T',3,2,0.5D0,DPBN(S1,:,1:2),3,-DPDJK-DPDLK,1,1D0,Dphi(1:2),1) 
                ENDIF

                ! derivatives based on upper bead
                IF (S1.EQ.NSEG) THEN             
                   CALL DGEMV('T',3,6,1D0,DPBN(S1+1,:,:),3,KJN,1,1D0, DDIST(1:6),1)
                   CALL DGEMV('T',3,6,0.5D0,DPBN(S1+1,:,:),3,DCT1DKJ,1,1D0, DCT1(1:6),1) 
                   CALL DGEMV('T',3,6,0.5D0,DPBN(S1+1,:,:),3,DCT2DLK - DCT2DJK,1,1D0, DCT2(1:6),1) 
                   CALL DGEMV('T',3,6,0.5D0,DPBN(S1+1,:,:),3,-DPDJK+DPDLK,1,1D0, Dphi(1:6),1)
                ELSE
                   CALL DGEMV('T',3,3,1D0,DPBN(S1+1,:,3:5),3,KJN,1,1D0, DDIST(4+4*S1:6+4*S1),1)
                   CALL DGEMV('T',3,2,1D0,DPBN(S1+1,:,1:2),3,KJN,1,1D0, DDIST(1:2),1)

                   CALL DGEMV('T',3,3,0.5D0,DPBN(S1+1,:,3:5),3,DCT1DKJ,1,1D0, DCT1(4+4*S1:6+4*S1),1)
                   CALL DGEMV('T',3,2,0.5D0,DPBN(S1+1,:,1:2),3,DCT1DKJ,1,1D0, DCT1(1:2),1) 

                   CALL DGEMV('T',3,3,0.5D0,DPBN(S1+1,:,3:5),3,DCT2DLK - DCT2DJK,1,1D0, DCT2(4+4*S1:6+4*S1),1)
                   CALL DGEMV('T',3,2,0.5D0,DPBN(S1+1,:,1:2),3,DCT2DLK - DCT2DJK,1,1D0, DCT2(1:2),1) 

                   CALL DGEMV('T',3,3,0.5D0,DPBN(S1+1,:,3:5),3,-DPDJK+DPDLK,1,1D0, DPHI(4+4*S1:6+4*S1),1)
                   CALL DGEMV('T',3,2,0.5D0,DPBN(S1+1,:,1:2),3,-DPDJK+DPDLK,1,1D0, Dphi(1:2),1) 
                ENDIF

                ! interpolate the minimal separation distance and its derivatives
                CT1 = SIGN(MIN(ABS(CT1),1D0),CT1)
                CT2 = SIGN(MIN(ABS(CT2),1D0),CT2)

                IF (.NOT.((CT1.GT.0.OR.CT1.LE.0).AND.(CT2.GT.0.OR.CT2.LE.0)&
                     & .AND.(PHI.GT.0.OR.PHI.LE.0))) THEN
                   PRINT*, 'ERROR IN GETSTERICENERGY (overlap between 0th nucleosome &
                        & and linker following n-th nucleosome):  bad CT1, CT2, PHI', CT1, CT2, PHI                   
                   STOP 1
                ENDIF

                !note: in nucsegcyl, first angle is always for the segment and 2nd for the nucleosome
                CALL EVALUATESPLINE(nucSEGCYL,CT2,CT1,PHI,MINDIST,DMDC2,DMDC1,DMDP)
                IF (.NOT.MINDIST.GE.0) THEN
                   PRINT*, 'ERROR IN GETSTERICENERGY (overlap between 0th nucleosome &
                        & and linker following n-th nucleosome): bad mindist', MINDIST, CT1, CT2, PHI, S1
                   STOP 1
                ENDIF

                DMIN = DMDC1*DCT1 + DMDC2*DCT2 + DMDP*DPHI                                      

                IF (DIST.LT.MINDIST) THEN
                   OVERLAP = MINDIST-DIST
                   NUCSEGE = NUCSEGE + CHAINP%ESTERIC*OVERLAP**2
                   NUCSEGGRAD = NUCSEGGRAD +  2* CHAINP%ESTERIC*OVERLAP*(DMIN-DDIST)
                ENDIF
             ENDIF
          ENDIF


          ! interactions between segments on first linker and those on linker following n-th nucleosome
          IF (CHAINP%SEGSTERICS) THEN
             DO S2 = 1,NSEG
                ! get the segments and the angles between them
                ! as well as derivatives wrt helix and bead coordinates
                DDIST = 0D0; DCT1 = 0D0; DCT2 = 0D0; DPHI = 0D0          

                IF (S1.EQ.1) THEN               
                   KJ = (PBN(S2+1,:) + PBN(S2,:) - CHAINP%VEC(4+4*S1:6+4*S1) &
                        & - PTOP)/2
                ELSE IF (S1.EQ.NSEG) THEN             
                   KJ = (PBN(S2+1,:) + PBN(S2,:) - PBOT &
                        & - CHAINP%VEC(4+4*(S1-1):6+4*(S1-1)))/2
                ELSE             
                   KJ = (PBN(S2+1,:) + PBN(S2,:) - CHAINP%VEC(4+4*S1:6+4*S1) &
                        & - CHAINP%VEC(4+4*(S1-1):6+4*(S1-1)))/2
                ENDIF
                LK = (PBN(S2+1,:) - PBN(S2,:))/2

                ! get the actual separation distance                
                DIST = SQRT(DOT_PRODUCT(KJ,KJ)); KJN = KJ/DIST*0.5D0
                IF (DIST.GT.SSMAXSEP) CYCLE

                ! get angles and dihedrals
                CALL GETANGLE(IJ,KJ,CT1,DCT1DIJ,DCT1DKJ)
                CALL GETANGLE(-KJ,LK,CT2,DCT2DJK,DCT2DLK)
                CALL GETDIHEDRAL(IJ,-KJ,LK,PHI,DPDIJ,DPDJK,DPDLK)

                ! Now for the derivatives
                ! derivatives wrt endpoints of 1st segment
                IF (S1.EQ.1) THEN
                   CALL DGEMV('T',3,6,-1D0,DPTOP,3,KJN,1,1D0,DDIST(1:6),1)   
                   CALL DGEMV('T',3,6,0.5D0,DPTOP,3,DCT1DIJ-DCT1DKJ,1,1D0, DCT1(1:6),1)
                   CALL DGEMV('T',3,6,0.5D0,DPTOP,3,DCT2DJK,1,1D0,DCT2(1:6),1) 
                   CALL DGEMV('T',3,6,0.5D0,DPTOP,3,DPDIJ+DPDJK,1,1D0,Dphi(1:6),1)
                ELSE
                   DDIST(4+4*(S1-1):6+4*(S1-1)) = DDIST(4+4*(S1-1):6+4*(S1-1)) - KJN
                   DCT1(4+4*(S1-1):6+4*(S1-1)) = DCT1(4+4*(S1-1):6+4*(S1-1)) + 0.5D0*(DCT1DIJ-DCT1DKJ)
                   DCT2(4+4*(S1-1):6+4*(S1-1)) = DCT2(4+4*(S1-1):6+4*(S1-1)) + 0.5D0*DCT2DJK
                   DPHI(4+4*(S1-1):6+4*(S1-1)) = DPHI(4+4*(S1-1):6+4*(S1-1)) + 0.5D0*(DPDIJ+DPDJK)
                ENDIF

                IF (S1.EQ.NSEG) THEN
                   CALL DGEMV('T',3,6,-1D0,DPBOT,3,KJN,1,1D0,DDIST(1:6),1)   
                   CALL DGEMV('T',3,6,0.5D0,DPBOT,3,-DCT1DIJ-DCT1DKJ,1,1D0, DCT1(1:6),1)
                   CALL DGEMV('T',3,6,0.5D0,DPBOT,3,DCT2DJK,1,1D0,DCT2(1:6),1) 
                   CALL DGEMV('T',3,6,0.5D0,DPBOT,3,-DPDIJ+DPDJK,1,1D0,Dphi(1:6),1)
                ELSE
                   DDIST(4+4*S1:6+4*S1) = DDIST(4+4*S1:6+4*S1) - KJN             
                   DCT1(4+4*S1:6+4*S1) = DCT1(4+4*S1:6+4*S1) + 0.5D0*(-DCT1DIJ - DCT1DKJ)             
                   DCT2(4+4*S1:6+4*S1) = DCT2(4+4*S1:6+4*S1) + 0.5D0*DCT2DJK             
                   DPHI(4+4*S1:6+4*S1) = DPHI(4+4*S1:6+4*S1) + 0.5D0*(-DPDIJ +DPDJK)             
                ENDIF

                ! derivatives wrt endpoints of 2nd segment
                ! derivatives based on bottom bead
                IF (S2.EQ.1) THEN
                   CALL DGEMV('T',3,6,1D0,DPBN(S2,:,:),3,KJN,1,1D0,DDIST(1:6),1)   
                   CALL DGEMV('T',3,6,0.5D0,DPBN(S2,:,:),3,DCT1DKJ,1,1D0, DCT1(1:6),1)
                   CALL DGEMV('T',3,6,0.5D0,DPBN(S2,:,:),3,-DCT2DLK-DCT2DJK,1,1D0,DCT2(1:6),1) 
                   CALL DGEMV('T',3,6,0.5D0,DPBN(S2,:,:),3,-DPDJK-DPDLK,1,1D0,Dphi(1:6),1) 
                ELSE
                   CALL DGEMV('T',3,3,1D0,DPBN(S2,:,3:5),3,KJN,1,1D0,DDIST(4+4*(S2-1):6+4*(S2-1)),1) 
                   CALL DGEMV('T',3,2,1D0,DPBN(S2,:,1:2),3,KJN,1,1D0, DDIST(1:2),1)

                   CALL DGEMV('T',3,3,0.5D0,DPBN(S2,:,3:5),3,DCT1DKJ,1,1D0,DCT1(4+4*(S2-1):6+4*(S2-1)),1) 
                   CALL DGEMV('T',3,2,0.5D0,DPBN(S2,:,1:2),3,DCT1DKJ,1,1D0, DCT1(1:2),1)

                   CALL DGEMV('T',3,3,0.5D0,DPBN(S2,:,3:5),3,-DCT2DLK-DCT2DJK,1,1D0,DCT2(4+4*(S2-1):6+4*(S2-1)),1)
                   CALL DGEMV('T',3,2,0.5D0,DPBN(S2,:,1:2),3,-DCT2DLK-DCT2DJK,1,1D0,DCT2(1:2),1) 

                   CALL DGEMV('T',3,3,0.5D0,DPBN(S2,:,3:5),3,-DPDJK-DPDLK,1,1D0,DPHI(4+4*(S2-1):6+4*(S2-1)),1) 
                   CALL DGEMV('T',3,2,0.5D0,DPBN(S2,:,1:2),3,-DPDJK-DPDLK,1,1D0,Dphi(1:2),1) 
                ENDIF

                ! derivatives based on upper bead
                IF (S2.EQ.NSEG) THEN             
                   CALL DGEMV('T',3,6,1D0,DPBN(S2+1,:,:),3,KJN,1,1D0, DDIST(1:6),1)
                   CALL DGEMV('T',3,6,0.5D0,DPBN(S2+1,:,:),3,DCT1DKJ,1,1D0, DCT1(1:6),1) 
                   CALL DGEMV('T',3,6,0.5D0,DPBN(S2+1,:,:),3,DCT2DLK - DCT2DJK,1,1D0, DCT2(1:6),1) 
                   CALL DGEMV('T',3,6,0.5D0,DPBN(S2+1,:,:),3,-DPDJK+DPDLK,1,1D0, Dphi(1:6),1)
                ELSE
                   CALL DGEMV('T',3,3,1D0,DPBN(S2+1,:,3:5),3,KJN,1,1D0, DDIST(4+4*S2:6+4*S2),1)
                   CALL DGEMV('T',3,2,1D0,DPBN(S2+1,:,1:2),3,KJN,1,1D0, DDIST(1:2),1)

                   CALL DGEMV('T',3,3,0.5D0,DPBN(S2+1,:,3:5),3,DCT1DKJ,1,1D0, DCT1(4+4*S2:6+4*S2),1)
                   CALL DGEMV('T',3,2,0.5D0,DPBN(S2+1,:,1:2),3,DCT1DKJ,1,1D0, DCT1(1:2),1) 

                   CALL DGEMV('T',3,3,0.5D0,DPBN(S2+1,:,3:5),3,DCT2DLK - DCT2DJK,1,1D0, DCT2(4+4*S2:6+4*S2),1)
                   CALL DGEMV('T',3,2,0.5D0,DPBN(S2+1,:,1:2),3,DCT2DLK - DCT2DJK,1,1D0, DCT2(1:2),1) 

                   CALL DGEMV('T',3,3,0.5D0,DPBN(S2+1,:,3:5),3,-DPDJK+DPDLK,1,1D0, DPHI(4+4*S2:6+4*S2),1)
                   CALL DGEMV('T',3,2,0.5D0,DPBN(S2+1,:,1:2),3,-DPDJK+DPDLK,1,1D0, Dphi(1:2),1) 
                ENDIF

                ! interpolate the minimal separation distance and its derivatives
                CT1 = SIGN(MIN(ABS(CT1),1D0),CT1)
                CT2 = SIGN(MIN(ABS(CT2),1D0),CT2)

                IF (.NOT.((CT1.GT.0.OR.CT1.LE.0).AND.(CT2.GT.0.OR.CT2.LE.0)&
                     & .AND.(PHI.GT.0.OR.PHI.LE.0))) THEN
                   PRINT*, 'ERROR IN GETSTERICENERGY (interaction between 2 linker segments): &
                        & bad CT1, CT2, PHI', CT1, CT2, PHI
                   STOP 1
                ENDIF

                CALL EVALUATESPLINE(SEGCYL,CT1,CT2,PHI,MINDIST,DMDC1,DMDC2,DMDP)
                IF (.NOT.MINDIST.GE.0) THEN
                   PRINT*, 'ERROR IN GETSTERICENERGY (interaction between 2 linker segments): &
                        & bad mindist', MINDIST, CT1, CT2, PHI, S1, S2
                   STOP 1
                ENDIF

                DMIN = DMDC1*DCT1 + DMDC2*DCT2 + DMDP*DPHI          

                IF (DIST.LT.MINDIST) THEN
                   OVERLAP = MINDIST-DIST   
                   SEGSEGE = SEGSEGE + CHAINP%ESTERIC*OVERLAP**2
                   SEGSEGGRAD = SEGSEGGRAD +  2* CHAINP%ESTERIC*OVERLAP*(DMIN-DDIST)
                ENDIF
             ENDDO
          ENDIF
       ENDDO
    ENDDO
   
  END SUBROUTINE GETSTERICENERGY


  SUBROUTINE GETTOPEDGE(CHAINP,PTOP,MATTOP,dPTOP,dMATTOP)
    ! get position and orientation of the top edge of the 1st bend
    ! also get derivatives wrt the 6 coordinates defining the regular helix
    ! orientation is given as a rotation matrix
    USE KEYS, ONLY : BPLEN
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(OUT) :: PTOP(3),MATTOP(3,3)
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: dPTOP(3,6),dMATTOP(3,3,6)
    INTEGER :: I
    DOUBLE PRECISION :: CT,ST,RAD,THETA
    DOUBLE PRECISION :: MAT(3,3),DMAT(3,3,3),MP(3,3)

    THETA = CHAINP%VEC(2)
    RAD = CHAINP%VEC(3)


    ! nucleosome orientation wrt canonical axes
    IF (PRESENT(DPTOP).OR.PRESENT(DMATTOP)) THEN
       CALL EUL2ROTMAT(CHAINP%VEC(4:6),MAT,DMAT)
    ELSE
       CALL EUL2ROTMAT(CHAINP%VEC(4:6),MAT)
    ENDIF

    ! TOP end orientation relative to nucleosome orientation
    CALL QUAT2ROTMAT(CHAINP%BENDS(1)%TP,MP)

    ! overall top end orientation
    CALL DGEMM('N','N',3,3,3,1D0,MAT,3,MP,3,0D0,MATTOP,3)

    IF (PRESENT(DMATTOP).OR.PRESENT(DPTOP)) THEN
       DMATTOP = 0D0
       ! derivatives of bottom end orientation
       DO I = 4,6
          CALL DGEMM('N','N',3,3,3,1D0,DMAT(:,:,I-3),3,MP,3,0D0,DMATTOP(:,:,I),3)
       ENDDO
    ENDIF

    ! position of minus end 
    PTOP = (/RAD,0D0,0D0/) ! bend center
    CALL DGEMV('N',3,3,1D0,MAT,3,CHAINP%BENDS(1)%POSP,1,1D0,PTOP,1)    

    ! derivatives of minus end position
    IF (PRESENT(DPTOP)) THEN
       DPTOP = 0D0
       DO I = 4,6
          CALL DGEMV('N',3,3,1D0,DMAT(:,:,I-3),3,CHAINP%BENDS(2)%POSP,1,0D0,DPTOP(:,I),1)
       ENDDO
       DPTOP(:,3) = (/1D0,0D0,0D0/)
    ENDIF

    ! shift by 1/2 bp
    PTOP = PTOP + BPLEN/2*MATTOP(:,3)
    IF (PRESENT(DPTOP)) THEN
       DPTOP = DPTOP + BPLEN/2*DMATTOP(:,3,:)
    ENDIF
  END SUBROUTINE GETTOPEDGE

  SUBROUTINE GETBOTTOMEDGE(CHAINP,PBOT,MATBOT,dPBOT,dMATBOT)
    ! get position and orientation of the bottom edge of the 2nd bend
    ! also get derivatives wrt the 6 coordinates defining the regular helix
    ! orientation is given as a rotation matrix
    USE KEYS, ONLY : BPLEN
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(OUT) :: PBOT(3),MATBOT(3,3)
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: dPBOT(3,6),dMATBOT(3,3,6)
    INTEGER :: I
    DOUBLE PRECISION :: CT,ST,RAD,THETA, EUL(3)
    DOUBLE PRECISION :: MAT(3,3),DMAT(3,3,3),MM(3,3)

    THETA = CHAINP%VEC(2)
    RAD = CHAINP%VEC(3)
    
    ! nucleosome orientation wrt canonical axes
    EUL = (/CHAINP%VEC(2)+CHAINP%VEC(4),CHAINP%VEC(5),CHAINP%VEC(6)/)
    IF (PRESENT(DPBOT).OR.PRESENT(DMATBOT)) THEN
       CALL EUL2ROTMAT(EUL,MAT,DMAT)
    ELSE
       CALL EUL2ROTMAT(EUL,MAT)
    ENDIF
    
    ! bottom end orientation relative to nucleosome orientation
    CALL QUAT2ROTMAT(CHAINP%BENDS(2)%TM,MM)

    ! overall bottom end orientation
    CALL DGEMM('N','N',3,3,3,1D0,MAT,3,MM,3,0D0,MATBOT,3)

    IF (PRESENT(DMATBOT).OR.PRESENT(DPBOT)) THEN
       DMATBOT = 0D0
       ! derivatives of bottom end orientation
       DO I = 4,6
          CALL DGEMM('N','N',3,3,3,1D0,DMAT(:,:,I-3),3,MM,3,0D0,DMATBOT(:,:,I),3)
       ENDDO
       DMATBOT(:,:,2) = DMATBOT(:,:,4)
    ENDIF

    ! position of minus end 
    CT = COS(THETA); ST = SIN(THETA)
    PBOT = (/RAD*CT,RAD*ST,CHAINP%VEC(1)/) ! bend center
    CALL DGEMV('N',3,3,1D0,MAT,3,CHAINP%BENDS(2)%POSM,1,1D0,PBOT,1)
    

    ! derivatives of minus end position
    IF (PRESENT(DPBOT)) THEN
       DPBOT = 0D0
       DO I = 4,6
          CALL DGEMV('N',3,3,1D0,DMAT(:,:,I-3),3,CHAINP%BENDS(2)%POSM,1,0D0,DPBOT(:,I),1)
       ENDDO
       DPBOT(:,1) = (/0D0,0D0,1D0/)
       DPBOT(:,2) = DPBOT(:,4) + (/-RAD*ST,RAD*CT,0D0/)
       DPBOT(:,3) = (/CT,ST,0D0/)
    ENDIF 
   
    ! shift by 1/2 bp
    PBOT = PBOT - BPLEN/2*MATBOT(:,3)
    IF (PRESENT(DPBOT)) THEN
       DPBOT = DPBOT - BPLEN/2*DMATBOT(:,3,:)
    ENDIF
  END SUBROUTINE GETBOTTOMEDGE

  SUBROUTINE SETFIBERDIAM(CRDS,DIAM,NUCRAD,RVAL,DRVAL)
    ! express the radius as a function of the other helix coordinates so as to
    ! fix the overall diameter of the helix (using cylinder steric tabulations)
    ! RVAL returns the radius, and DRVAL returns its derivatives wrt
    ! the other helix coords
    ! DIAM is the desired overall diameter out to the furthest point of the nucleosome
    ! NUCRAD is the steric radius of the nucleosomes

    USE CYLINDERUTILS, ONLY : NUCCYL, NUCCYLSET
    USE TRICUBSPLINEUTILS, ONLY : EVALUATESPLINE
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: CRDS(6),DIAM, NUCRAD
    DOUBLE PRECISION, INTENT(OUT) :: RVAL, DRVAL(6)
    DOUBLE PRECISION :: R, T, H, A, B
    DOUBLE PRECISION :: CA, SA, CB, SB, CNT, SNT
    DOUBLE PRECISION :: X, DXDT, DXDA,DXDB,Y
    DOUBLE PRECISION :: CT1, CT2, PHI,MINDIST,DMDC1,DMDC2,DMDP,TMP

    R = CRDS(3); T = CRDS(2); H = CRDS(1)
    A = CRDS(4); B = CRDS(5)
    CA = COS(A); SA = SIN(A); CB = COS(B); SB = SIN(B)


    ! place a reference cylinder so it just touches the fiber on the x-axis
    ! find the closest approach distance for the first nucleosome

    IF (.NOT.NUCCYLSET) THEN
       PRINT*, 'Cannot set fiber diameter if nucleosome cylinder tables have not been preset'
       STOP 1
    ENDIF

    CT2 = 0D0
    CT1 = SB*SA
    PHI = ATAN2(-SB*CA,CB)

    CALL EVALUATESPLINE(NUCCYL,CT1,CT2,PHI,MINDIST,DMDC1,DMDC2,DMDP)
    RVAL = DIAM/2 - MINDIST+NUCRAD    

    DRVAL = 0D0
    TMP = CB*CB+SB*SB*CA*CA
    DRVAL(5) = -(CB*SA*DMDC1 - CA / TMP*DMDP)
    DRVAL(4) = -(SB*CA*DMDC1 + CB*SB*SA/TMP*DMDP       )

    IF (.NOT.RVAL.GT.0) THEN
       PRINT*, 'ERROR IN SETFIBERDIAM: bad R', RVAL   
       print*, CRDS
       STOP 1
    ENDIF
  END SUBROUTINE SETFIBERDIAM

  SUBROUTINE ROTMAT2RELEULER(MAT1,MAT2,NCRD,DMAT1,DMAT2,CBETA,ALPHAGAMMA,DCBETA,DAG)
    ! given 2 rotation matrices and their derivatives
    ! get the relative cos(beta) and alpha+gamma euler angles to go from 1st coord system
    ! to 2nd
    ! also get derivatives of these wrt the original coordinates    
    ! NCRD is the number of coordinates for derivatives
    ! WARNING: assumes the same coordinates are involved in MAT1 and MAT2 derivatives
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NCRD
    DOUBLE PRECISION,INTENT(IN) :: MAT1(3,3),MAT2(3,3),DMAT1(3,3,NCRD),DMAT2(3,3,NCRD)
    DOUBLE PRECISION, INTENT(OUT) :: CBETA,ALPHAGAMMA,DCBETA(NCRD),DAG(NCRD)
    DOUBLE PRECISION :: Z1,Z2,X3,Y3,NZ,X1,X2
    DOUBLE PRECISION :: DZ1(NCRD),DZ2(NCRD),DX3(NCRD),DY3(NCRD),DX1(NCRD),DX2(NCRD)

    DCBETA = 0D0; DAG = 0D0

    ! get the relative euler angles
    CBETA = DOT_PRODUCT(MAT1(:,3),MAT2(:,3))
    Z1 = DOT_PRODUCT(MAT1(:,1),MAT2(:,3))
    Z2 = DOT_PRODUCT(MAT1(:,2),MAT2(:,3))

    NZ = Z1*Z1+Z2*Z2
    IF (NZ.LT.NZTINY) THEN       
       X1 = DOT_PRODUCT(MAT1(:,1),MAT2(:,1))
       X2 = DOT_PRODUCT(MAT1(:,2),MAT2(:,1))
       ALPHAGAMMA = ATAN2(X2,X1)
    ELSE       
       X3 = DOT_PRODUCT(MAT1(:,3),MAT2(:,1))
       Y3 = DOT_PRODUCT(MAT1(:,3),MAT2(:,2))      
       ALPHAGAMMA = ATAN2(Z1*Y3-Z2*X3,-Z2*Y3-Z1*X3)              
    ENDIF

    ! derivatives of euler angles
    CALL DGEMV('T',3,NCRD,1D0,DMAT1(:,3,:),3,MAT2(:,3),1,0D0,DCBETA,1)
    CALL DGEMV('T',3,NCRD,1D0,DMAT2(:,3,:),3,MAT1(:,3),1,1D0,DCBETA,1)       

    IF (NZ.LT.NZTINY) THEN            
       CALL DGEMV('T',3,NCRD,1D0,DMAT1(:,1,:),3,MAT2(:,1),1,0D0,DX1,1)
       CALL DGEMV('T',3,NCRD,1D0,DMAT2(:,1,:),3,MAT1(:,1),1,1D0,DX1,1)  
       CALL DGEMV('T',3,NCRD,1D0,DMAT1(:,2,:),3,MAT2(:,1),1,0D0,DX2,1)
       CALL DGEMV('T',3,NCRD,1D0,DMAT2(:,1,:),3,MAT1(:,2),1,1D0,DX2,1)  
       DAG = (X1*DX2 - X2*DX1)/(X1*X1+X2*X2)
    ELSE
       CALL DGEMV('T',3,NCRD,1D0,DMAT1(:,1,:),3,MAT2(:,3),1,0D0,DZ1,1)
       CALL DGEMV('T',3,NCRD,1D0,DMAT2(:,3,:),3,MAT1(:,1),1,1D0,DZ1,1)  
       CALL DGEMV('T',3,NCRD,1D0,DMAT1(:,2,:),3,MAT2(:,3),1,0D0,DZ2,1)
       CALL DGEMV('T',3,NCRD,1D0,DMAT2(:,3,:),3,MAT1(:,2),1,1D0,DZ2,1)  
       CALL DGEMV('T',3,NCRD,1D0,DMAT1(:,3,:),3,MAT2(:,1),1,0D0,DX3,1)
       CALL DGEMV('T',3,NCRD,1D0,DMAT2(:,1,:),3,MAT1(:,3),1,1D0,DX3,1)  
       CALL DGEMV('T',3,NCRD,1D0,DMAT1(:,3,:),3,MAT2(:,2),1,0D0,DY3,1)
       CALL DGEMV('T',3,NCRD,1D0,DMAT2(:,2,:),3,MAT1(:,3),1,1D0,DY3,1)  
       DAG = (-Z2*DZ1+Z1*DZ2+Y3*DX3-X3*DY3)/NZ
    ENDIF

  END SUBROUTINE ROTMAT2RELEULER

END MODULE ENERGYUTILS
