      PROGRAM cluster
C
C     version 3.2
C     Apr 23, 1999
C
C     Written by Michael Schaefer, (c) 1998, 1999
C                ULP Strasbourg
C                schaefer@brel.u-strasbg.fr
C
      IMPLICIT NONE
C     MXNX: Max No of atoms in each structure to be clustered;
C     MXNC: Max No of conformations in the trajectory.
      INTEGER      MXNX,MXNC
      PARAMETER    (MXNX=2000,MXNC=200000)
      INTEGER      NX,NC,IClust(MXNC),NCount(MXNC)
      REAL         X(3,MXNX),WX(MXNX),RMSIC(MXNC)
      LOGICAL      IsRest(MXNC)
C
C     MXNC0 is the maximum No of structures that can be used
C     to determine the cluster centers, and the maximum No of
C     cluster centers that can be stored:
      INTEGER      MXNBIN,MXNRMS,MXNC0
      PARAMETER    (MXNBIN=200,MXNRMS=126800000,MXNC0=13000)
      INTEGER      NRMSD(MXNBIN)
      INTEGER      N0,N0Inp,NC0
      REAL         X0(3,MXNX,MXNC0),del,RMSCut
      REAL         XClust(3,MXNX,MXNC0)
      PARAMETER    (del=0.1)
      INTEGER      NClust,LastN,I2RMSD(0:MXNC0),JRMSD(MXNRMS)
      INTEGER      ICl0(MXNC0),IC02IC(MXNC0)
C
      INTEGER      MXNCL
      PARAMETER    (MXNCL=MXNC)
      INTEGER      ICC(MXNCL),NMemb(MXNCL)
      REAL         RMSDAV(MXNCL)
C
      INTEGER      MODE,IER
      REAL         U(3,3),T(3),RMS
C
      CHARACTER*80 fnam,fbin
      CHARACTER*3  key
      REAL         R,RMSMin
      INTEGER      IC,JC,IX,I,J,J1,JMin,IPass,IRest,MaxIt,N10
      DATA         N0/10/,N0Inp/10/,MODE/0/,MaxIt/5/,RMSCut/2.0/
C-----------------------------------------------------------------------
      WRITE(6,*)
      WRITE(6,'(1X,A)') '  XYZ file >'
      READ(5,'(A)') fnam
      WRITE(6,'(1X,A)') fnam
      OPEN(71,FILE=fnam,STATUS='OLD',ERR=999)
      WRITE(6,'(1X,A)') '  temporary binary XYZ file >'
      READ(5,'(A)') fbin
      WRITE(6,'(1X,A)') fbin
      OPEN(72,FILE=fbin,STATUS='UNKNOWN',FORM='UNFORMATTED',ERR=998)
      WRITE(6,'(1X,A)') '  No of coordinate sets >'
      READ(5,*) NC
      WRITE(6,'(1X,I8)') NC
      IF (NC.GT.MXNC) GOTO 997
      WRITE(6,'(1X,A)') '  No of atoms per coordinate set >'
      READ(5,*) NX
      WRITE(6,'(1X,I8)') NX
      IF (NX.GT.MXNX) GOTO 996
      DO IX=1,NX
         WX(IX) = 1.0/REAL(NX)
      ENDDO
      WRITE(6,'(1X,A)')
     &     '  use every Nth structure for clustering >'
      WRITE(6,'(1X,A)')
     &     '  (0 = use maximum number of structures) >'
      READ(5,*) N0Inp
      WRITE(6,'(1X,I8)') N0Inp
      N0 = 1+INT(REAL(NC-1)/REAL(MXNC0))
      IF (N0Inp.GT.N0) N0=N0Inp
      WRITE(6,'(1X,A)')
     &     '  maximum number of clustering iterations >'
      READ(5,*) MaxIt
      WRITE(6,'(1X,I8)') MaxIt
      WRITE(6,'(1X,A)')
     &     '  cluster RMSD cutoff [A] >'
      READ(5,*) RMSCut
      WRITE(6,'(1X,F8.3)') RMSCut
      WRITE(6,*)
      WRITE(6,'(1X,A,F7.3,A)')
     &     '  CLUSTER: RMSD CUTOFF FOR CLUSTERS ',RMSCut,' A'
      WRITE(6,*)
C
      WRITE(6,*) '  reading coordinates, writing to binary file...'
      DO IC=1,NC
         READ(71,*) key,I
C         IF ((key.NE.'XYZ').OR.(I.NE.IC)) GOTO 998
         IF ((key.NE.'XYZ')) GOTO 998
         DO IX=1,NX
            READ(71,*) (X(I,IX),I=1,3)
         ENDDO
         WRITE(72) X
      ENDDO
      CLOSE(71)
      REWIND(72)
C
      DO I=1,MXNBIN
         NRMSD(I) = 0
      ENDDO
      DO IC=1,NC
         IClust(IC) = 0
      ENDDO
      DO I=1,MXNCL
         ICC(I) = 0
      ENDDO
      WRITE(6,*)
C
      IPass = 1
      NClust = 0
100   CONTINUE
      WRITE(6,*) '  clustering pass No ',IPass
      IF (IPass.GT.MaxIt) THEN
         WRITE(6,*) '  exit clustering because MaxIt ',MaxIt,' exceeded'
         GOTO 200
      ENDIF
      I = 0
      NC0 = 0
      DO IC=1,NC
         READ(72) X
         IF (IClust(IC).EQ.0) THEN
            I=I+1
            IF (MOD(I,N0).EQ.0) THEN
               NC0 = NC0+1
               IsRest(IC) = .FALSE.
               IC02IC(NC0) = IC
               DO IX=1,NX
                  X0(1,IX,NC0) = X(1,IX)
                  X0(2,IX,NC0) = X(2,IX)
                  X0(3,IX,NC0) = X(3,IX)
               ENDDO
            ELSE
C           this conformation is not used in the clustering call, and
C           it is not yet assigned to any of the clusters:
               IsRest(IC) = .TRUE.
            ENDIF
         ELSE
            IsRest(IC) = .FALSE.
         ENDIF
      ENDDO
      REWIND(72)
      WRITE(6,'(A,I3,A,I8,A)')'         using every ',N0
     &     ,' th structure (total of ',NC0,')...'
      LastN = NClust
      DO IC=1,NC0
         ICl0(IC) = 0
      ENDDO
      CALL Clustn(NC0,NX,MXNX,X0,WX,ICl0,RMSCut,
     &     I2RMSD,JRMSD,MXNRMS,
     &     NRMSD,MXNBIN,del,IPass,
     &     NCount,NClust,ICC,MXNC0)
      WRITE(6,*) '       # of clusters ',NCLust
C     stop if there are no new clusters determined:
      IF (NClust.EQ.LastN) THEN
         WRITE(6,*) '  exit clustering because No of clusters unchanged'
         GOTO 200
      ENDIF
C     store new cluster centers:
      DO J=LastN+1,NClust
         JC = ICC(J)
         DO IX=1,NX
            XClust(1,IX,J) = X0(1,IX,JC)
            XClust(2,IX,J) = X0(2,IX,JC)
            XClust(3,IX,J) = X0(3,IX,JC)
         ENDDO
      ENDDO
C     translate cluster index based on IC0 index back to IC:
      DO IC=1,NC0
         IClust(IC02IC(IC)) = ICl0(IC)
      ENDDO
C     translate conformation index of the new cluster centers:
      DO I=LastN+1,NClust
         ICC(I) = IC02IC(ICC(I))
      ENDDO
C     check which structures can be assigned to the new clusters:
      N10 = MAX(NC/10,1)
      J1 = LastN+1
      DO IC=1,NC
         READ(72) X
         IF (IsRest(IC)) THEN
            DO J=J1,NClust
               CALL U3BEST(WX,X,XClust(1,1,J),NX,MODE,RMS,U,T,IER)
               RMS = SQRT(RMS)
               IF (RMS.LT.RMSCut) THEN
                  IClust(IC) = J
                  GOTO 150
               ENDIF
            ENDDO
         ENDIF
150      CONTINUE
C       IF (MOD(IC,N10).EQ.0) WRITE(6,*) '  done with IC ',IC
      ENDDO
      REWIND(72)
      IF (NClust.GE.MXNC0) THEN
         WRITE(6,*)
     &        '  exit clustering because Max No of clusters reached'
         GOTO 200
      ENDIF
C     check for termination criterion:
      IRest = 0
      DO IC=1,NC
         IF (IClust(IC).EQ.0) IRest=IRest+1
      ENDDO
      WRITE(6,*) '       # unassigned ',IRest
      IF (N0.EQ.1) THEN
         WRITE(6,'(2A)') '  exit clustering, all unassigned structures '
     &        ,'used in last clustering call'
         GOTO 200
      ENDIF
      IF (IRest.GT.0) THEN
         IPass = IPass+1
         N0 = 1+INT(REAL(IRest-1)/REAL(MXNC0))
         I = N0Inp/(2**(IPass-1))
         IF (I.GT.N0) N0=I
         GOTO 100
      ENDIF
200   CONTINUE
C
      WRITE(6,*)
      WRITE(6,*) '  total No of cluster centers ',NClust
      WRITE(6,*)
      WRITE(6,*) '  histogram of RMSDs (intervals of 0.1A):'
      DO I=1,200
         R = 0.05001+REAL(I-1)*0.1
         WRITE(6,'(1X,F6.2,1X,I16)') R,NRMSD(I)
      ENDDO
      WRITE(6,*)
C
C     assign all structures to the nearest cluster center:
      DO IC=1,NC
         IClust(IC) = 0
         RMSIC(IC) = -1.0
      ENDDO
      DO I=1,MXNCL
         NMemb(I) = 0
         RMSDAV(I) = 0.0
      ENDDO
      WRITE(6,*)
C
C     first the cluster centers:
      DO J=1,NClust
         ICLust(ICC(J)) = J
         RMSIC(ICC(J)) = 0.0
      ENDDO
C     then all other structures:
      N10 = MAX(NC/10,1)
      DO IC=1,NC
         READ(72) X
         IF (IClust(IC).EQ.0) THEN
            DO J=1,NClust
               CALL U3BEST(WX,X,XClust(1,1,J),NX,MODE,RMS,U,T,IER)
               IF ((J.EQ.1).OR.(RMS.LT.RMSMin)) THEN
                  RMSMin = RMS
                  JMin = J
               ENDIF
            ENDDO
            RMSMin = SQRT(RMSMin)
            RMSIC(IC) = RMSMin
            IF (RMSMin.LT.RMSCut) THEN
               IClust(IC) = JMin
               NMemb(JMin) = NMemb(JMin)+1
               RMSDAV(JMin) = RMSDAV(JMin)+RMSMin
            ELSE
               IClust(IC) = -JMin
            ENDIF
         ENDIF
C       IF (MOD(IC,N10).EQ.0) WRITE(6,*) '  done with IC ',IC
      ENDDO
      CLOSE(72)
C
      WRITE(6,'(1X,A)') '  cluster assignment file >'
      READ(5,'(A)') fnam
      WRITE(6,'(1X,A)') fnam
      OPEN(81,FILE=fnam,STATUS='UNKNOWN',ERR=999)
      WRITE(81,'(A)')
     &     '# Nstruct  Cluster  Rmsic'
      DO IC=1,NC
         WRITE(81,'(I6,I9,1X,F7.3)') IC,IClust(IC),RMSIC(IC)
      ENDDO
      CLOSE(81)
      WRITE(6,*)
C
      WRITE(6,'(1X,A)') '  cluster center file >'
      READ(5,'(A)') fnam
      WRITE(6,'(1X,A)') fnam
      OPEN(81,FILE=fnam,STATUS='UNKNOWN',ERR=999)
      WRITE(81,'(A)')
     &     '#     Nclust  struct   aveRMSD   Nmembrs(excl. center)'
      DO J=1,NClust
         WRITE(81,'(1X,I9,1X,I9,1X,F9.3,1X,I9)')
     &        J,ICC(J),RMSDAV(J)/REAL(NMemb(J)),NMemb(J)
      ENDDO
      CLOSE(81)
      WRITE(6,*)
C
      STOP
995   WRITE(6,'(1X,A)') '  cluster: ERROR opening binary file!'
      STOP
996   WRITE(6,'(1X,A,I8)') '  cluster: ERROR MXNX exceeded! '
      STOP
997   WRITE(6,'(1X,A,I8)') '  cluster: ERROR MXNC exceeded! '
      STOP
998   WRITE(6,'(1X,A,I8)') '  cluster: ERROR reading input file! IC ',IC
      STOP
999   WRITE(6,'(1X,A)') '  cluster: ERROR opening input file!'
      STOP
C
      END
C
C
C=======================================================================
      SUBROUTINE Clustn(NC,NX,MXNX,X,WX,IClust,RMSCut,
     &     I2RMSD,JRMSD,MXNRMS,
     &     NRMSD,MXNBIN,del,IPass,
     &     NCount,NClust,ICC,MXNCLU)
      IMPLICIT NONE
      INTEGER   NC,NX,MXNX,IClust(*)
      REAL      X(3,MXNX,*),WX(*),RMSCut,del
      INTEGER   I2RMSD(0:*),JRMSD(*),MXNRMS
      INTEGER   NRMSD(*),MXNBIN,IPass
      INTEGER   NCount(*),NClust,ICC(*),MXNCLU
C
      INTEGER   MODE,IER
      REAL      U(3,3),T(3),RMS
C
      INTEGER   I2,NFAIL,IC,JC,IND,J,ICMax,NMax,N10
      DATA      MODE/0/
C-----------------------------------------------------------------------
      N10 = MAX(NC/10,1)
      NFAIL = 0
      I2 = 0
      I2RMSD(0) = 0
      DO IC=1,NC
         DO JC=IC+1,NC
            CALL U3BEST(WX,X(1,1,IC),X(1,1,JC),NX,MODE,RMS,U,T,IER)
            RMS = SQRT(RMS)
            IF (IPass.EQ.1) THEN
               IND = 1+INT(RMS/del)
               IF (IND.LE.MXNBIN) THEN
                  NRMSD(IND) = NRMSD(IND)+1
               ELSE
                  NFAIL = NFAIL+1
               ENDIF
            ENDIF
            IF (RMS.LT.RMSCut) THEN
               IF (I2.GE.MXNRMS) THEN
                  WRITE(6,*) '  clustn: ERROR, MXNRMS exceeded!'
                  STOP
               ENDIF
               I2 = I2+1
               JRMSD(I2) = JC
            ENDIF
         ENDDO
         I2RMSD(IC) = I2
C       IF (MOD(IC,N10).EQ.0) WRITE(6,*) '  done with IC ',IC
      ENDDO
      IF (NFAIL.GT.0) WRITE(6,*) '  clustn: WARNING, ',NFAIL,
     &     ' RMSDs outside array boundary!'
C
100   CONTINUE
      DO IC=1,NC
         NCount(IC) = 0
      ENDDO
      DO IC=1,NC
         IF (IClust(IC).EQ.0) THEN
            DO J=I2RMSD(IC-1)+1,I2RMSD(IC)
               JC = JRMSD(J)
               IF (IClust(JC).EQ.0) THEN
                  NCount(IC) = NCount(IC)+1
                  NCount(JC) = NCount(JC)+1
               ENDIF
            ENDDO
         ENDIF
      ENDDO
C     search for structure with highest RMSD<RMSCut count:
      ICMax = 0
      NMax = 0
      DO IC=1,NC
         IF (NCount(IC).GT.NMax) THEN
            NMax = NCount(IC)
            ICMax = IC
         ENDIF
      ENDDO
C     if there is not at least 2 structures left with RMSD<1, stop here:
      IF (ICMax.LE.0) GOTO 200
C     if maximum No of clusters reached, stop here:
      IF (NClust.GE.MXNCLU) GOTO 200
      NClust = NClust+1
      ICC(NClust) = ICMax
C     take members of cluster out of the game (set IClust to NClust):
C     first the cluster center itself:
      IClust(ICMax) = NClust
      DO IC=1,ICMax-1
         IF (IClust(IC).EQ.0) THEN
            DO J=I2RMSD(IC-1)+1,I2RMSD(IC)
               IF (JRMSD(J).EQ.ICMax) THEN
                  IClust(IC) = NClust
                  GOTO 150
               ENDIF
            ENDDO
         ENDIF
150      CONTINUE
      ENDDO
      DO J=I2RMSD(ICMax-1)+1,I2RMSD(ICMax)
         JC = JRMSD(J)
         IF (IClust(JC).EQ.0) THEN
            IClust(JC) = NCLust
         ENDIF
      ENDDO
      GOTO 100      
200   CONTINUE
C
      RETURN
      END
C
C
