C ***********************************************************************************************************   
C ** This code was prepared for calculating stress-strain relationship based on Mohr-Coulomb with associated*
C ** flow rule for ABAQUS/Explicit                                                                          *      
C ** In this code, the material parameters are as follows:                                                  *
C ** E (Elastic modulus)                                                                                    *
C ** v (Nu)                                                                                                 *
C ** C (Cohesion)                                                                                           *
C ** Max fi (to calculate FoS, this term is friction angle of soil)                                         *
C ** min fi (this term was eliminated for FoS calculation)                                                  *
C ** Sai (dilation angle)                                                                                   *
C ** Plastic shear strain (this term relates to softening/hardening, please ignore this term)               *      
C ** In order to monitor the output, please make sure that SDV output is available from Step menu           *
C ** SDV 5 is the magnitude of shear strain in model                                                        *
C ** SDV 6 is the magnitude of friction angle                                                               *
C ** In order to use PYTHON file, FORTRAN VUMAT, and ABAQUS ".cae" file, please follow the instruction      *
C ***********************************************************************************************************
C ***********************************************************************************************************
C ** This code was written by                                                                               *
C ** Yousef Javanmardi (Y.javanmardi@ucl.ac.uk)                                                             *                                                                                      *
C ** Morteza Naeij (morteza.naeij@ut.ac.ir)                                                                 *
C ** and is a part of supplementary data for manuscript with title:                                         *
C ** A novel procedure for finite element method to calculate the minimum factor of safety                  *
C **                                                                                                        *
C ***********************************************************************************************************  
C ** Please make sure that VUMAT file (.FOR) introduces to the ABAQUS properly from                         *
C ** job\General\user subroutine                                                                            *
C ***********************************************************************************************************
C ***********************************************************************************************************       
      SUBROUTINE VUMAT(
C     Read only -
     1 NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL,
     2 STEPTIME, TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH,
     3 PROPS, DENSITY, STRAININC, RELSPININC,
     4 TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD,
     5 STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD,
     6 TEMPNEW, STRETCHNEW, DEFGRADNEW, FIELDNEW,
C     Write only -
     7 STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW)
C
C      INCLUDE 'VABA_PARAM.INC'
C
      DIMENSION PROPS(NPROPS), DENSITY(NBLOCK), COORDMP(NBLOCK),
     1 CHARLENGTH(NBLOCK), STRAININC(NBLOCK, NDIR+NSHR),
     2 RELSPININC(NBLOCK, NSHR), TEMPOLD(NBLOCK),
     3 STRETCHOLD(NBLOCK, NDIR+NSHR),DEFGRADOLD(NBLOCK,NDIR+NSHR+NSHR),
     4 FIELDOLD(NBLOCK, NFIELDV), STRESSOLD(NBLOCK, NDIR+NSHR),
     5 STATEOLD(NBLOCK, NSTATEV), ENERINTERNOLD(NBLOCK),
     6 ENERINELASOLD(NBLOCK), TEMPNEW(NBLOCK),
     7 STRETCHNEW(NBLOCK, NDIR+NSHR),DEFGRADNEW(NBLOCK,NDIR+NSHR+NSHR),
     8 FIELDNEW(NBLOCK, NFIELDV), STRESSNEW(NBLOCK,NDIR+NSHR),
     9 STATENEW(NBLOCK, NSTATEV), ENERINTERNNEW(NBLOCK),
     1 ENERINELASNEW(NBLOCK)
C
      CHARACTER*80 CMNAME
      PARAMETER( ONE = 1.D0, TWO = 2.D0, PI=3.1415926535D0)
      DIMENSION DELAS(4,4),STRESSNEWINT(4),STRESSP(3),
     1 STRESSPR(3), LOCA(3),DE(3,3), RG_RSI(3), RF_RSI(3),STRESSPB(3),
     2 STRAININCINT(4),STRESSINT(4), XLAM1(3),EPSTEN(13),RGAMA_REPS(3),
     3 RF_REPSP(3)
      INTEGER BLOCK_COUNT
      DOUBLE PRECISION F1,F2,F,GAMA_F
C     MODEL CONSTANTS
      E=PROPS(1)
      V=PROPS(2)
      C=PROPS(3)
      PHI_P=PROPS(4)*PI/180.D0
      PHI_RES=PROPS(5)*PI/180.D0
      PSI_P=PROPS(6)*PI/180.D0
      GAMA_F=PROPS(7)
C     ELASTIC CONSTITUTIVE MATRIX      
      F1=E*(1-V)/((ONE+V)*(ONE-TWO*V))
      F2=E*V/((ONE+V)*(ONE-TWO*V))
      DELAS(1,1:4)=(/F1, F2, F2, 0.D0/)
      DELAS(2,1:4)=(/F2, F1, F2, 0.D0/)
      DELAS(3,1:4)=(/F2, F2, F1, 0.D0/)      
      DELAS(4,1:4)=(/0.D0, 0.D0, 0.D0, F1-F2/)
C     ELASTIC MATRIX FOR PRINCIPAL STRESSES
      DE(1,1:3)=(/F1, F2, F2/)
      DE(2,1:3)=(/F2, F1, F2/)
      DE(3,1:3)=(/F2, F2, F1/)     
C     BEGIN MAIN LOOP      
      DO BLOCK_COUNT=1, NBLOCK
C     GETTING STATE VARIABLE
      EPSTEN(1)=STATEOLD(BLOCK_COUNT,1)
      EPSTEN(2)=STATEOLD(BLOCK_COUNT,2)
      EPSTEN(3)=STATEOLD(BLOCK_COUNT,3)
      EPSTEN(4)=STATEOLD(BLOCK_COUNT,4)
      EPSTEN(5)=STATEOLD(BLOCK_COUNT,5)
      EPSTEN(6)=STATEOLD(BLOCK_COUNT,6)
      EPSTEN(9)=STATEOLD(BLOCK_COUNT,9)
      EPSTEN(10)=STATEOLD(BLOCK_COUNT,10)
      EPSTEN(11)=STATEOLD(BLOCK_COUNT,11)
C     CALCULATION OF GAMA_P MOBILIZED
      GAMA_P=EPSTEN(5)+(2.D0*SQRT(EPSTEN(1)**2+EPSTEN(2)**2-EPSTEN(1)*
     1 EPSTEN(2)+3.D0*EPSTEN(4)**2)/SQRT(3.D0))
C     STRESS AND STRAIN INCREMENT SUBSTITUTION
      STRESSINT(1:4)=STRESSOLD(BLOCK_COUNT,1:4)
      STRAININCINT(1:4)=STRAININC(BLOCK_COUNT,1:4)
C     ELASTIC ESTIMATION OF NEW STRESSES      
      STRESSNEWINT=STRESSINT+MATMUL(DELAS,STRAININCINT)
C     PRINCIPAL VALUES OF ELASTIC ESTIMATION STRESSES
      STRESSMEAN=0.5D0*(STRESSNEWINT(1)+STRESSNEWINT(2))
      OFFSET=SQRT(0.25D0*(STRESSNEWINT(1)-STRESSNEWINT(2))**2+
     1 STRESSNEWINT(4)**2)
      STRESSP(1)=STRESSMEAN-OFFSET
      STRESSP(2)=STRESSMEAN+OFFSET
      STRESSP(3)=STRESSNEWINT(3)
C     VALUE OF ROTATION OF PRINCIPAL AXIS      
      TWOTHETA=ATAN2(STRESSNEWINT(4),0.5D0*(STRESSNEWINT(1)-
     1 STRESSNEWINT(2)))
C     SORTING PRINCIPLAL STRESSES AND FINDIND THEIR LOCATION
      LOCA(1)=MINLOC(STRESSP,1)
      STRESSPR(1)=MIN(STRESSP(1),STRESSP(2),STRESSP(3))
      LOCA(3)=MAXLOC(STRESSP,1)
      STRESSPR(3)=MAX(STRESSP(1),STRESSP(2),STRESSP(3))
      LOCA(2)=3-MOD(LOCA(1)+LOCA(3),3)
      STRESSPR(2)=STRESSP(LOCA(2))
C     UPDATING PHI AND PSI
 
        PHI=PHI_P
        PSI=PSI_P
      
      SPHI=SIN(PHI)
      SPSI=SIN(PSI) 
C     ROUND F/ ROUND FI
      RF_RPHI=0
C     Round PHI/ ROUND GAMA

        RPHI_RGAMA=0
 
                  
C     CONTROLLING THE OCCURANCE OF YIELDING
      F=0.5D0*(STRESSPR(3)-STRESSPR(1))+0.5D0*(STRESSPR(3)+STRESSPR(1))
     1 *SPHI-C*COS(PHI)
C     ELASTIC BEHAVIOR
      IF ((F .LE. -1.D-6) .OR. (TOTALTIME .EQ. 0.0D0)) THEN
        STRESSNEW(BLOCK_COUNT,1:4)=STRESSNEWINT
      ELSE
C     PLASTIC BEHAVIOR

C     AND 3 FOR TRIAXIAL EXTENSION
C     ROUND GAMA/ ROUND_EPSILON

            RGAMA_REPS=0.D0

        RF_REPSP=RF_RPHI*RPHI_RGAMA*RGAMA_REPS
C     DETERMINATION OF MODE, 1 FOR TRIAXIAL COMPRESSION, 2 FOR REGULAR STRESS  
        IAREA=2
        IF (STRESSPR(2) .EQ. STRESSPR(3)) THEN
            IAREA=1
            LOCA(1:3)=(/1,2,3/)
        ELSEIF (STRESSPR(1) .EQ. STRESSPR(2)) THEN
            IAREA=3
            LOCA(1:3)=(/1,3,2/)
        ENDIF
C       NORMAL VECTOR TO YIELD SURFACE AND PLASTIC POTENTIAL SURFACE
        IF (IAREA .EQ. 1) THEN
            RF_RSI(1:3)=(/0.5D0*(-1.D0+SPHI),0.25D0*(1.D0+SPHI),
     1                  0.25D0*(1.D0+SPHI)/)
            RG_RSI(1:3)=(/0.5D0*(-1.D0+SPSI),0.25D0*(1.D0+SPSI),
     1                  0.25D0*(1.D0+SPSI)/)
        ELSEIF (IAREA .EQ. 2) THEN
            RF_RSI(1:3)=(/0.5D0*(-1.D0+SPHI),0.0 D0,0.5D0*(1.D0+SPHI)/)
            RG_RSI(1:3)=(/0.5D0*(-1.D0+SPSI),0.0 D0,0.5D0*(1.D0+SPSI)/)
        ELSEIF (IAREA .EQ. 3) THEN
            RF_RSI(1:3)=(/0.25D0*(-1.D0+SPHI),0.25D0*(-1.D0+SPHI),
     1                  0.5D0*(1.D0+SPHI)/)
            RG_RSI(1:3)=(/0.25D0*(-1.D0+SPSI),0.25D0*(-1.D0+SPSI),
     1                  0.5D0*(1.D0+SPSI)/)
        ENDIF  
        H_M1=-1.D0*(RF_REPSP(1)*RG_RSI(1)+ RF_REPSP(2)*RG_RSI(2)+  
     1         RF_REPSP(3)*RG_RSI(3)) 
        XLAM1=MATMUL(DE,RG_RSI)
       H_M=H_M1+XLAM1(1)*RF_RSI(1)+XLAM1(2)*RF_RSI(2)+XLAM1(3)*RF_RSI(3)
       SOFTTERM=F*XLAM1(3)/H_M
        STRESSPR=STRESSPR-F*XLAM1/H_M
        DO I_PRE=1,3
            STRESSPB(LOCA(I_PRE))=STRESSPR(I_PRE)
        ENDDO
        STRESSMEANB=0.5D0 * (STRESSPB(1)+STRESSPB(2))
        OFFSETB=0.5D0 * (STRESSPB(2)-STRESSPB(1))
        STRESSNEW(BLOCK_COUNT,1)=STRESSMEANB+OFFSETB*COS(TWOTHETA)
        STRESSNEW(BLOCK_COUNT,2)=STRESSMEANB-OFFSETB*COS(TWOTHETA)
        STRESSNEW(BLOCK_COUNT,3)=STRESSPB(3)
        STRESSNEW(BLOCK_COUNT,4)=OFFSETB*SIN(TWOTHETA)
C     UPDATING STATE VARIABLE
      STATENEW(BLOCK_COUNT,1)=STRAININCINT(1)
      STATENEW(BLOCK_COUNT,2)=STRAININCINT(2)
      STATENEW(BLOCK_COUNT,3)=STRAININCINT(3)
      STATENEW(BLOCK_COUNT,4)=STRAININCINT(4)
      STATENEW(BLOCK_COUNT,8)=H_M1
      STATENEW(BLOCK_COUNT,9)= EPSTEN(9)+STRAININCINT(1)
      STATENEW(BLOCK_COUNT,10)= EPSTEN(10)+STRAININCINT(2)
      STATENEW(BLOCK_COUNT,11)= EPSTEN(11)+STRAININCINT(4)
      ENDIF
C     update the Gamma_P
      STATENEW(BLOCK_COUNT,5)=GAMA_P
      STATENEW(BLOCK_COUNT,6)=phi*180.D0/PI
      STATENEW(BLOCK_COUNT,7)=PSI*180.D0/PI
      STATENEW(BLOCK_COUNT,12)=f
      STATENEW(BLOCK_COUNT,13)=xlam1(1)      
      ENDDO
      RETURN
      END SUBROUTINE VUMAT
        
        
            
        
        
      
             