C********************************************************************************
C*      ESTE PROGRAMA RESOLVE ESPALHAMENTO de duas particulas de spin zero
C*      COM POTENCIAL DO TIPO WOODS-SAXON 
C*
C*	Ele é uma adaptação da rotina principal do programa global.f 
C*	desenvolvido Por Luiz Carlos Chamon.
C*							Renato Negrão Fev 2010
C********************************************************************************

C___________________Inicio do Programa Principal______________________________

      IMPLICIT REAL*8(A-H,O-Z)
      dimension ff(0:800),ffc(0:800),vf(40002),vc(40002),VREAL(40002),
     *VIMAG(40002),ro1(0:800),ro2(0:800),roc1(0:800),roc2(0:800),
     *rcoef(9),qcoef(9),ENUCLEON(13),SIGN(13),SIGP(13)
      COMMON/um/ ZR(40002),ZI(40002),SR(0:2001),SI(0:2001),HR(40002),
     *HI(40002)
      COMMON/dois/ DR,FK,ECM,ETA,NPR,L
      COMMON/tres/ FMP,FMA,ZP,ZA,LMAX
      COMMON/quatro/ Vef(40002)

C
C     Definição das Constantes
C
      hbar=197.329d0
      pi=3.141592653d0
      alpha=1.d0/137.0359895d0	
C
C     Entrada de dados
C
      OPEN(UNIT=5,FILE="input.in",STATUS="unknown")
      READ(5,*) ELAB,FMP,ZP,FMA,ZA
      READ(5,*) RMAX,LMAX,DR,R0C
      READ(5,*) V0,RR0,AR,W0,RI0,AI
      READ(5,*) imprime,lef,londa,dang,nexp
C
C     Verifica se as Particulas são identicas ou não
C
      if(zp.eq.za.and.fmp.eq.fma) then
      IOPART=1
      LDELTA=2
      LMAX=(LMAX/2)
      LMAX=LMAX*2
      else
      IOPART=0
      LDELTA=1
      endif
      if(IOPART.eq.1) write(6,1011)
 1011 format(/' Identical Particles',/)
C
C     Saida
C
      OPEN(UNIT=6,FILE="output.out",STATUS="unknown")
      WRITE(6,1000) ELAB,FMP,ZP,FMA,ZA
 1000 FORMAT(1X,'ELAB=',F7.2,' MP=',F5.1,' ZP=',F4.1,' MA=',F5.1,
     *' ZA=',F4.1)
      FMI=FMP*FMA/(FMP+FMA)
      ECM=ELAB*FMA/(FMA+FMP)
      FK=DSQRT(ECM*FMI/20.913067D0)
      WRITE(6,1001) FK,RMAX,DR,R0C
 1001 FORMAT(1X,'K=',F5.2,' RMAX=',F5.2,' DR=',F5.3,' RC=',F5.2)
      ZZE2=hbar*alpha*ZP*ZA
      ETA=ZZE2/DSQRT(4.d0*ECM*20.913067d0/FMI)
C
C      Estimativa de um bom valor para DR
C
      EDR=1.d0/(5.d0*FK)
      IF(DR.gt.EDR) WRITE(6,1030) EDR
 1030 FORMAT(/,' Warning: DR seems too large.',/,
     *' Try a DR value smaller than ',f6.4,/)

      NPR=IDINT(RMAX/DR)+2
      IF(NPR.GT.40002) THEN
      WRITE(6,1014)
 1014 FORMAT(' Rmax/dr is too large')
      CLOSE(UNIT=5)
      CLOSE(UNIT=6)
      STOP
      ENDIF
      IF(imprime.eq.1) WRITE(6,1015)
 1015 FORMAT(/,1x,' R     Vc               V              W',/)
      DO 230 IR=1,NPR
      R=DFLOAT(IR)*DR

C
C     Potencial Coulombiano - Esfera de com distribuicao de carga uniforme 
C
      RC=R0C*(FMP**(1.d0/3.d0)+FMA**(1.d0/3.d0))
      IF(R.gt.RC) THEN
      VC(IR)=ZP*ZA/R*hbar/137.035d0
      ELSE
      VC(IR)=ZP*ZA*hbar/137.035d0*(3.d0*rc**2-R**2)/2.d0/rc**3
      ENDIF
C
C     Definicao do Potencial Nuclear
C
      RI=RI0*(FMP**(1.d0/3.d0)+FMA**(1.d0/3.d0))
      XI=(R-RI)/AI
      VIMAG(IR)=-W0/(1.d0+DEXP(XI))

      RR=RR0*(FMP**(1.d0/3.d0)+FMA**(1.d0/3.d0))
      XR=(R-RR)/AR
      VREAL(IR)=-V0/(1.d0+DEXP(XR))

C
C      CALCULO DE HI
C
      HI(IR)=-VIMAG(IR)/ECM
      IF(imprime.eq.1) THEN
      WRITE(6,1016) R,VC(IR),VREAL(IR),VIMAG(IR)
 1016 FORMAT(1x,f5.2,3(1x,e14.7))
      ENDIF
 230  CONTINUE
C
C      Verificando se Rmax esta razoavel
C
      IF(VC(NPR).ne.0..and.dabs(VREAL(NPR)/VC(NPR)).gt.0.001d0) 
     *WRITE(6,1018)
 1018 FORMAT(/,' Warning: Rmax seems to be too small',/)
C
C      RESOLVENDO PARA TODAS AS ONDAS PARCIAIS
C
      DO 280 L=0,LMAX,LDELTA
      IF(l.eq.lef) WRITE(6,1019) lef
 1019 FORMAT(/,1x,'Effective potential for l=',i3,//,
     *1x,'  R        Vef',/)
C
C      Calculo de HR para cada valor de L
C
      DO 290 IR=1,NPR
      R=DFLOAT(IR)*DR
      RO=R*FK
      HR(IR)=1.d0-DFLOAT(L*(L+1))/(RO**2)-(vc(IR)+vreal(IR))/ECM
      IF(l.eq.lef) THEN
      WRITE(6,1020) R,-ECM*(HR(IR)-1.d0)
 1020 FORMAT(1x,f5.2,1x,e14.7)
      ENDIF
 290  CONTINUE
      ZR(1)=1.d0
      ZI(1)=0.d0
      CALL INTEG(2)
      IF(LONDA.eq.l) THEN
      WRITE(6,1021) LONDA
 1021 FORMAT(/,1x,'Partial wave function for l=',i3,//,
     *1x,'   R      |U_L(R)|',/)
      DO IR=1,NPR
      R=DFLOAT(IR)*DR
      FMOD=dsqrt(ZR(IR)**2+ZI(IR)**2)
      WRITE(6,1022) r,fmod
 1022 FORMAT(1x,f5.2,1x,e14.7)
      ENDDO
      ENDIF
 280  CONTINUE
C
C      CALCULO DAS SECOES DE CHOQUE de Espalhamento Elastico e Fusao
C
      CALL CROSS(ETA,ECM,ZZE2,FK,RMAX,LMAX,SR,SI,IOPART,DANG,NEXP)
      IF(ETA.eq.0.) THEN
      CLOSE(UNIT=5)
      CLOSE(UNIT=6)
      STOP
      ENDIF
      SECWKB=0.d0
      SECEF=0.d0
      TLANT=1.d0
      LCR=-1
      DO L=0,LMAX,LDELTA
      IF(LCR.eq.-1) THEN
      DO 500 IR=1,NPR
      R=DFLOAT(IR)*DR
      RO=R*FK
      VEF(IR)=ECM*DFLOAT(L*(L+1))/(RO**2)+VC(IR)+VREAL(IR)
 500  CONTINUE
      CALL FUSAO(TLANT,SWKB,SEF,LCR)
      SECWKB=SECWKB+SWKB
      SECEF=SECEF+SEF
      ENDIF
      ENDDO
      IF(IOPART.eq.1) THEN
      LCR=LCR-2
      SECWKB=2.d0*SECWKB
      SECEF=2.d0*SECEF
      ELSE
      LCR=LCR-1
      ENDIF
      WRITE(6,5000) LCR,SECWKB,SECEF
 5000 FORMAT(/,1x,' critical L value = ',i3,
     *' (negative value = you should increase Rmax or Lmax)',
     */,1x,' BPM cross section = ',e14.7,' mb',/,1x,
     *' Effective Fusion cross section = ',e14.7,' mb')
      CLOSE(UNIT=5)
      CLOSE(UNIT=6)
      STOP
      END

C___________________Fim do Programa Principal_________________________________

C*****************************************************************************
C*    SUBROUTINE INTEG
C*
C*****************************************************************************
      SUBROUTINE INTEG(IOPT)
C
C      RESOLVE A EQUACAO DIFERENCIAL
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION FC(2001),FCP(2001),GC(2001),GCP(2001),SIGMA(2001)
      COMMON/um/ ZR(40002),ZI(40002),SR(0:2001),SI(0:2001),HR(40002),
     *HI(40002)
      COMMON/dois/ DR,FK,ECM,ETA,NPR,L
      PI=4.0d0*datan(1.0d0)
      DRO2=(DR*FK)**2/12.d0
C
C     IOPT=0        NAO NORMALIZA, SO INTEGRA
C     IOPT=OUTRO    INTEGRA E NORMALIZA
C
      X=2.d0*ZR(1)-10.d0*DRO2*(HR(1)*ZR(1)-HI(1)*ZI(1))
      Y=2.d0*ZI(1)-10.d0*DRO2*(HR(1)*ZI(1)+HI(1)*ZR(1))
      IF(L.EQ.1) THEN
      X=X-DRO2*HR(1)*ZR(1)
      Y=Y-DRO2*HR(1)*ZI(1)
      ENDIF
      Z=1.d0+DRO2*HR(2)
      T=-DRO2*HI(2)
      ZR(2)=(X*Z-Y*T)/(Z**2+T**2)
      ZI(2)=(X*T+Y*Z)/(Z**2+T**2)
      DO 20 IR=3,NPR
      X=2.d0*ZR(IR-1)-ZR(IR-2)-10.d0*DRO2*(HR(IR-1)*ZR(IR-1)-
     *HI(IR-1)*ZI(IR-1))-DRO2*(HR(IR-2)*ZR(IR-2)-HI(IR-2)*ZI(IR-2))
      Y=2.d0*ZI(IR-1)-ZI(IR-2)-10.d0*DRO2*(HR(IR-1)*ZI(IR-1)+
     *HI(IR-1)*ZR(IR-1))-DRO2*(HR(IR-2)*ZI(IR-2)+HI(IR-2)*ZR(IR-2))
      Z=1.d0+DRO2*HR(IR)
      T=-DRO2*HI(IR)
      U=Z**2+T**2
      ZR(IR)=(Z*X-T*Y)/U
      ZI(IR)=(X*T+Y*Z)/U
      FMOD=DSQRT(ZR(IR)**2+ZI(IR)**2)
      IF(FMOD.GT.1.d20) THEN
      DO 22 IRR=1,IR
      ZR(IRR)=ZR(IRR)*1.d-20
      ZI(IRR)=ZI(IRR)*1.d-20
 22   CONTINUE
      ENDIF
 20   CONTINUE
C
C      OBTENDO AS DERIVADAS DAS FUNCOES DE ONDA
C
      IF(IOPT.EQ.0) RETURN
      R=DFLOAT(NPR-2)*DR
      RO=R*FK
      DRO=DR*FK
      ULR=(ZR(NPR-4)-8.d0*ZR(NPR-3)+8.d0*ZR(NPR-1)-ZR(NPR))/
     *(12.d0*DRO)
      ULI=(ZI(NPR-4)-8.d0*ZI(NPR-3)+8.d0*ZI(NPR-1)-ZI(NPR))/
     *(12.d0*DRO)
C
C      OBTENDO AS DERIVADAS DE FC E GC
C
      RK1=RO-2.d0*DRO
      RK2=RO-DRO
      RK3=RO
      RK4=RO+DRO
      RK5=RO+2.d0*DRO
      LM=L+1
      IF(LM.EQ.1) LM=LM+1
      CALL COU(RK1,RK2,ETA,LM,FC,FCP,GC,GCP,SIGMA)
      FC1=FC(L+1)
      FC2=FCP(L+1)
      GC1=GC(L+1)
      GC2=GCP(L+1)
      CALL COU(RK4,RK5,ETA,LM,FC,FCP,GC,GCP,SIGMA)
      FC4=FC(L+1)
      FC5=FCP(L+1)
      GC4=GC(L+1)
      GC5=GCP(L+1)
      CALL COU(RK3,RK4,ETA,LM,FC,FCP,GC,GCP,SIGMA)
      FL=(FC1-8.d0*FC2+8.d0*FC4-FC5)/(12.d0*DRO)
      GL=(GC1-8.d0*GC2+8.d0*GC4-GC5)/(12.d0*DRO)
C
C      OBTENDO A MATRIZ S
C
      FMOD=ZR(NPR-2)**2+ZI(NPR-2)**2
      IF(FMOD.EQ.0.d0) THEN
      WRITE(6,1000) L
 1000 FORMAT(//,' ZEROU A FUNCAO DE ONDA EM R=RMAXIMO PARA L=',I3)
      CLOSE(UNIT=5)
      CLOSE(UNIT=6)
      STOP
      ENDIF
      X=(ULR*ZR(NPR-2)+ULI*ZI(NPR-2))/FMOD
      Y=(ULI*ZR(NPR-2)-ULR*ZI(NPR-2))/FMOD
      A=FL-X*FC(L+1)+Y*GC(L+1)
      B=GL-Y*FC(L+1)-X*GC(L+1)
      C=X*FC(L+1)+Y*GC(L+1)-FL
      D=Y*FC(L+1)-X*GC(L+1)+GL
      SR(L)=(A*C+B*D)/(C**2+D**2)
      SI(L)=(B*C-D*A)/(C**2+D**2)
C
C      NORMALIZANDO A FUNCAO DE ONDA
C
      UBR=(FC(L+1)*(1.d0+SR(L))+GC(L+1)*SI(L))*DCOS(SIGMA(L+1))/2.d0+
     *(GC(L+1)*(SR(L)-1.d0)-FC(L+1)*SI(L))*DSIN(SIGMA(L+1))/2.d0
      UBI=(GC(L+1)*(1.d0-SR(L))+FC(L+1)*SI(L))*DCOS(SIGMA(L+1))/2.d0+
     *(FC(L+1)*(SR(L)+1.d0)+GC(L+1)*SI(L))*DSIN(SIGMA(L+1))/2.d0
      FMOD2=UBR**2+UBI**2
      RNORM=DSQRT(FMOD2/FMOD)
      FASE=0.d0
      IF(ZR(NPR-2).NE.0.d0) THEN
      FASE=DATAN(ZI(NPR-2)/ZR(NPR-2))
      IF(ZR(NPR-2).LT.0.d0) FASE=FASE+PI
      ELSE
      IF(ZI(NPR-2).GT.0.d0) FASE=PI/2.d0
      IF(ZI(NPR-2).LT.0.d0) FASE=3.d0*PI/2.d0
      ENDIF
      IF(FASE.LT.0.d0) FASE=FASE+2.d0*PI
      FASE2=0.d0
      IF(UBR.NE.0.d0) THEN
      FASE2=DATAN(UBI/UBR)
      IF(UBR.LT.0.d0) FASE2=FASE2+PI
      ELSE
      IF(UBI.GT.0.d0) FASE2=PI/2.d0
      IF(UBI.LT.0.d0) FASE2=3.d0*PI/2.d0
      ENDIF
      IF(FASE2.LT.0.d0) FASE2=FASE2+2.d0*PI
      FI=FASE2-FASE
      COR=RNORM*DCOS(FI)
      COI=RNORM*DSIN(FI)
      DO 30 IR=1,NPR
      SALR=ZR(IR)
      SALI=ZI(IR)
      ZR(IR)=SALR*COR-SALI*COI
 30   ZI(IR)=SALI*COR+SALR*COI
      RETURN
      END

C*****************************************************************************
C*    SUBROUTINE CROSS: CALCULA SECOES DE CHOQUE
C*
C*****************************************************************************
      SUBROUTINE CROSS(ETA,ECM,ZZE2,FK,RMAX,LMAX,SR,SI,IOPT,dang,nexp)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 UM,CI,FCOU,F,SUM,SL(0:2001),FCOUMENOS,FMENOS,SUMENOS
      DIMENSION PL(2001),SR(0:2001),SI(0:2001),SIGMA(2001),FC(2001),
     *FCP(2001),GC(2001),GCP(2001),PMENOS(2001)
C
C      IOPT= 0 - PARTICULAS NAO IDENTICAS
C            1 - PARTICULAS IDENTICAS
C
      IF(IOPT.NE.0.AND.IOPT.NE.1) THEN
      CLOSE(UNIT=5)
      CLOSE(UNIT=6)
      STOP
      ENDIF
      IF(IOPT.eq.1) LDELTA=2
      IF(IOPT.eq.0) LDELTA=1
      CI=DCMPLX(0.d0,1.d0)
      UM=DCMPLX(1.d0,0.d0)
      PI=4.d0*datan(1.0d0)
      WRITE(6,100)
 100  FORMAT(/,1X,'  L     !S!          FI',5x,
     *'-  S(L)=!S!*EXP(2*I*FI)',/)
      DO L=0,LMAX,LDELTA
      SL(L)=DCMPLX(SR(L),SI(L))
      FA=0.d0
      IF(SR(L).NE.0.d0) THEN
      FA=DATAN(SI(L)/SR(L))
      IF(SR(L).LT.0.d0) FA=FA+PI
      ELSE
      IF(SI(L).GT.0.d0) FA=PI/2.d0
      IF(SI(L).LT.0.d0) FA=3.d0*PI/2.d0
      ENDIF
      IF(FA.LT.0.d0) FA=FA+2.d0*PI
      FA=FA*90.d0/PI
      FMOD=DSQRT(SR(L)**2+SI(L)**2)
      WRITE(6,111) L,FMOD,FA
 600  CONTINUE
 111  FORMAT(1X,I4,2x,E10.4,2x,E10.4)
      ENDDO
C
C      Verificando se Lmax esta razoavel - Elemento da matriz de espalhamento
C      Deve tender a 1 para eventos perifericos com grande momento angular
C
      X=dabs(fa-180.d0)
      y=dabs(fa)
      IF(y.le.x) x=y
      If(dabs(1.d0-fmod).gt.0.001d0.or.x.gt.0.1d0) THEN
      WRITE(6,9000)
 9000 FORMAT(' Lmax is too small')
      CLOSE(UNIT=5)
      CLOSE(UNIT=6)
      STOP
      ENDIF
C
C      CALCULO DOS DESLOCAMENTOS DE FASE DE COULOMB, POLINOMIOS
C      DE LEGENDRE E FCOULOMB
C
      LL=LMAX+1
      R=RMAX*FK
      RP=R+.1d0
      CALL COU(R,RP,ETA,LL,FC,FCP,GC,GCP,SIGMA)
      WRITE(6,200)
 200  FORMAT(/,1X,'ELASTIC SCATTERING ANGULAR DISTRIBUTION',//,
     *' TETAcm  SIG (mb/sr) SIG/SIG_Ruth SIG_Ruth',/)
      nang=idint(180./dang)-1

      OPEN(UNIT=7,FILE="xmgrace.plot",STATUS="unknown")
      WRITE(7,201)
 201  FORMAT(/,1X,'@g0 type LOGY',/)	
      DO 310 J=1,nang
      T=dfloat(j)*dang*PI/180.d0
      CALL PLEGO(T,PL)
      TMENOS=PI-T
      CALL PLEGO(TMENOS,PMENOS)
      FCOU=-(ETA/(2.d0*FK*DSIN(T/2.d0)**2))*CDEXP(2.d0*CI*SIGMA(1)-
     *CI*ETA*DLOG(DSIN(T/2.d0)**2))
      FCOUMENOS=-(ETA/(2.d0*FK*DSIN(TMENOS/2.d0)**2))*
     *CDEXP(2.d0*CI*SIGMA(1)-CI*ETA*DLOG(DSIN(TMENOS/2.d0)**2))
      SUM=(0.d0,0.d0)
      SUMENOS=(0.d0,0.d0)
      DO 320 L=0,LMAX,LDELTA
      SUM=SUM+(2.d0*DFLOAT(L)+1.d0)*CDEXP(2.d0*CI*SIGMA(L+1))*
     *(1.d0-SL(L))*PL(L+1)
      SUMENOS=SUMENOS+(2.d0*DFLOAT(L)+1.d0)*CDEXP(2.d0*CI*SIGMA(L+1))*
     *(1.d0-SL(L))*PMENOS(L+1)
 320  CONTINUE
      F=FCOU+SUM*CI/(2.d0*FK)
      FMENOS=FCOUMENOS+SUMENOS*CI/(2.d0*FK)
      IF(ETA.ne.0.d0.and.IOPT.EQ.0) ST2=(CDABS(F)/CDABS(FCOU))**2
      IF(ETA.ne.0.d0.and.IOPT.EQ.1) ST2=
     *(CDABS(F+FMENOS)/CDABS(FCOU+FCOUMENOS))**2
      IF(IOPT.eq.0) ST=10.d0*CDABS(F)**2
      IF(IOPT.EQ.1) ST=10.d0*CDABS(F+FMENOS)**2
      SRUTH=0.d0
      IF(IOPT.eq.0.and.ETA.ne.0.d0) 
     *SRUTH=10.d0*(ZZE2/(4.d0*ECM))**2/dsin(T/2.d0)**4
      IF(IOPT.eq.1.and.ETA.ne.0.d0) SRUTH=10.d0*(ETA/(2.d0*FK))**2*
     *(1.d0/dsin(T/2.d0)**4+1.d0/dcos(T/2.d0)**4+
     *8.d0*dcos(2.d0*ETA*dlog(dsin(T/2.d0)/dcos(T/2.d0)))/dsin(T)**2)
      TETA=dfloat(j)*DANG
      WRITE(7,9500) TETA,ST2
 310  WRITE(6,9500) TETA,ST,ST2,SRUTH
      WRITE(7,202)
 202  FORMAT(/,1X,'end',/)
      CLOSE(UNIT=7)
C
C      CALCULO DO CHI-QUADRADO
C
      CHI=0.d0
      IF(NEXP.EQ.0) GO TO 500
      NEXP=0
      OPEN(UNIT=4,FILE="Global.ENT",STATUS="unknown")
 400  READ(4,*,END=410) T,SEXP,ERRO
      NEXP=NEXP+1
      T=T*PI/180.d0
      TMENOS=PI-T
      CALL PLEGO(T,PL)
      CALL PLEGO(TMENOS,PMENOS)
      FCOU=-(ETA/(2.d0*FK*DSIN(T/2.d0)**2))*CDEXP(2.d0*CI*SIGMA(1)-
     *CI*ETA*DLOG(DSIN(T/2.d0)**2))
      SUM=(0.d0,0.d0)
      FCOUMENOS=-(ETA/(2.d0*FK*DSIN(TMENOS/2.d0)**2))*
     *CDEXP(2.d0*CI*SIGMA(1)-CI*ETA*DLOG(DSIN(TMENOS/2.d0)**2))
      SUMENOS=(0.d0,0.d0)
      DO 420 L=0,LMAX,LDELTA
      SUM=SUM+(2.d0*DFLOAT(L)+1.)*CDEXP(2.d0*CI*SIGMA(L+1))*
     *(1.d0-SL(L))*PL(L+1)
      SUMENOS=SUMENOS+(2.d0*DFLOAT(L)+1.d0)*CDEXP(2.d0*CI*SIGMA(L+1))*
     *(1.d0-SL(L))*PMENOS(L+1)
 420  CONTINUE
      F=FCOU+SUM*CI/(2.d0*FK)
      FMENOS=FCOUMENOS+SUMENOS*CI/(2.d0*FK)
      IF(ETA.ne.0.d0.and.IOPT.EQ.0) ST=(CDABS(F)/CDABS(FCOU))**2
      IF(ETA.ne.0.d0.and.IOPT.EQ.1) ST=
     *(CDABS(F+FMENOS)/CDABS(FCOU+FCOUMENOS))**2
      IF(ETA.eq.0.d0.and.IOPT.EQ.0) ST=10.d0*CDABS(F)**2
      IF(ETA.eq.0.d0.and.IOPT.EQ.1) ST=10.d0*CDABS(F+FMENOS)**2
      CHI=CHI+((ST-SEXP)/ERRO)**2
      TETA=T*180.d0/PI
      ERROU=100.d0*Dabs((ST-SEXP)/ERRO)
      IF(nexp.eq.1) WRITE(6,9003)
 9003 FORMAT(/)
      WRITE(6,9002) TETA,100.d0*(ST-SEXP)/SEXP,ERROU
 9002 FORMAT(1X,'TETA=',F6.2,'   (Steo-Sexp)/Sexp= ',e9.2,' %',
     *' (Steo-Sexp)/DS= ',E8.2,' %')
      go to 400
 410  CONTINUE
      CHIN=CHI/dFLOAT(NEXP)
      WRITE(6,9022) NEXP,CHI,CHIN
 9022 FORMAT(/,1X,'Nexp=',I3,' Total CHI2=',E9.3,
     *' CHI2/NEXP=',E9.3,/)
 500  CONTINUE
C
C      CALCULO DAS SECOES DE CHOQUE DE REACAO
C
      IF(eta.eq.0.d0) THEN
      WRITE(6,9501)
      ELSE
      WRITE(6,9601)
      ENDIF
 9501 FORMAT(/,' Partial Reaction and Elastic Cross Sections',/,
     *'   L    SIGR(L)    SIGE(L)',/)
 9601 FORMAT(/,' Partial Reaction Cross Sections',/,
     *'   L    SIGR(L)',/)
      REACAO=0.d0
      ELASTIC=0.d0
      FPART=1.d0
      IF(IOPT.EQ.1) FPART=2.d0
      DO 330 L=0,LMAX,LDELTA
      SIGRL=(10.d0*PI/FK**2)*(2.d0*DFLOAT(L)+1.d0)*
     *(1.d0-CDABS(SL(L))**2)
      SIGRL=FPART*SIGRL
      REACAO=REACAO+SIGRL
      SIGEL=(10.d0*PI/FK**2)*(2.d0*DFLOAT(L)+1.d0)*
     *CDABS(1.d0-SL(L))**2
      SIGEL=FPART*SIGEL
      ELASTIC=ELASTIC+SIGEL
      IF(ETA.eq.0.d0) THEN
      WRITE(6,1000) L,SIGRL,SIGEL
      ELSE
      WRITE(6,1100) L,SIGRL
      ENDIF
 330  CONTINUE
      WRITE(6,501) REACAO
      IF(eta.eq.0.d0) THEN
      T=0.d0
      CALL PLEGO(T,PL)
      TMENOS=PI-T
      CALL PLEGO(TMENOS,PMENOS)
      SUM=(0.d0,0.d0)
      SUMENOS=(0.d0,0.d0)
      DO 620 L=0,LMAX,LDELTA
      SUM=SUM+(2.d0*DFLOAT(L)+1.d0)*(1.d0-SL(L))*PL(L+1)
      SUMENOS=SUMENOS+(2.d0*DFLOAT(L)+1.d0)*(1.d0-SL(L))
     **PMENOS(L+1)
 620  CONTINUE
      F=SUM*CI/(2.d0*FK)
      IF(IOPT.eq.1) F=F+SUMENOS*CI/(2.d0*FK)
      TOTAL=(4.d0*10.d0*pi/FK)*dimag(F)
      WRITE(6,502) ELASTIC,TOTAL
      ENDIF
 9500 FORMAT(1X,F6.2,3(2x,E10.4))
 1000 FORMAT(1X,I4,2x,E10.4,2x,e10.4)
 1100 FORMAT(1X,I4,2x,E10.4)
 501  FORMAT(/,1X,'Total Reaction Cross Section = ',E10.4,' MB')
 502  FORMAT(1X,'Total Elastic Cross Section = ',E10.4,' MB',/,
     *1x,'Total Cross Section = ',e10.4,' MB')
      RETURN
      END

C*****************************************************************************
C*	SUBROUTINE FUSAO
C*
C*****************************************************************************
      SUBROUTINE FUSAO(TLANT,SWKB,SEF,LCR)
C
C      Calcula secao de choque de fusao para cada onda parcial 
C      com o modelo de penetracao de barreiras
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/dois/ DR,FK,ECM,ETA,NPR,L
      COMMON/tres/ FMP,FMA,ZP,ZA,LMAX
      Common/quatro/ VEF(40002)
      PI=4.d0*datan(1.0d0)      
      FMI=FMP*FMA/(FMP+FMA)
C
C     Verificando se Ecm e' menor que Vef(NPR)
C
      IF(ECM.le.VEF(NPR)) THEN
      IF(L.eq.0) THEN
      WRITE(6,1001)
 1001 FORMAT(1x,' Numerical problems in the fusion calculations')
      CLOSE(UNIT=5)
      CLOSE(UNIT=6)
      STOP
      ENDIF
      IF(TLANT.gt.0.01d0) THEN
      WRITE(6,1001)
      CLOSE(UNIT=5)
      CLOSE(UNIT=6)
      STOP
      ENDIF
      TLWKB=0.d0
      TLEF=0.d0
      SWKB=0.d0
      SEF=0.d0
      RETURN
      ENDIF
C
C     Calculo do raio da barreira
C
      IR=NPR
 20   IR=IR-1
      If(VEF(IR).gt.VEF(IR+1).and.IR.ne.5) go to 20
      IRL=IR+1
      R=dfloat(IRL)*DR
C
C     Se nao encontrou a barreira, imprime Lcritico e para
C
      IF(r.le.1.0d0) THEN
      LCR=L
      RETURN
      ENDIF
C
C      Parametros das Barreiras
C
      RBL=dfloat(IRL)*DR
      VBL=vef(IRL)
C
C     Segunda derivada de VEF
C
      DV=(-VEF(IRL-2)+16.d0*VEF(IRL-1)-30.d0*VBL+16.d0*
     *VEF(IRL+1)-VEF(IRL+2))/(12.d0*DR**2)
      HWL=dsqrt(dabs((197.3d0**2/(FMI*931.5d0))*DV))
 10   CONTINUE
C
C     Calculo dos demais HW
C
      IF(ECM.ge.(VBL-0.2d0*HWL)) THEN
      HWWKB=HWL
      go to 29
      ENDIF
C
C     Calculo dos pontos de retorno
C
      IMAX=NPR
      DO 25 IR=NPR,IRL,-1
      IF((VEF(IR)-ECM).lt.0.d0) IMAX=IR-1
 25   CONTINUE
 35   IR=IR-1
      IF((VEF(IR)-ECM).gt.0.d0.and.IR.ge.2) go to 35 
      IMIN=IR+1
      SOMA_JWKB = 0.d0
      DO 45 IR=IMIN,IMAX
      SOMA_JWKB=SOMA_JWKB+dsqrt(0.1914d0*FMI*(VEF(IR)-ECM))
 45   CONTINUE
      HWWKB=2.d0*pi*(VBL-ECM)/(SOMA_JWKB*DR)
 29   CONTINUE
C
C     Parametros de barreira efetivos
C
      IF(fmi.gt.8.d0) THEN
      HWEF=HWWKB*(1.d0+0.1d0*(FMI-8.d0))
      ELSE
      HWEF=HWWKB
      ENDIF
 30   CONTINUE
C
C     Calculo dos TL
C
      X=2.d0*pi*(VBL-ECM)/HWWKB
      IF(X.ge.-49.d0.and.X.le.49.d0) TLWKB=1.d0/(1.d0+dexp(x))
      IF(X.lt.-49.d0) TLWKB=1.d0
      IF(X.gt.49.d0) TLWKB=0.d0
      x=2.d0*pi*(VBL-ECM)/HWEF
      IF(X.ge.-49.d0.and.X.le.49.d0) TLEF=1.d0/(1.d0+dexp(x))
      IF(X.lt.-49.d0) TLEF=1.d0
      IF(X.gt.49.d0) TLEF=0.d0
      TLANT=TLWKB
C
C     Calculo das secoes de choque
C
      SWKB=(2.d0*dfloat(L)+1)*TLWKB*(pi/FK**2)*10.d0
      SEF=(2.d0*dfloat(L)+1)*TLEF*(pi/FK**2)*10.d0
      IF(L.eq.0) Write(6,1000)
 1000 FORMAT(/,' The calculation of BPM cross sections is a possible',
     */,' source of numerical errors. Check the convergence by ',/,
     *' changing parameters, such as: Rmax, Lmax, dr, etc.'
     *//,1x,'  L   RB     VB      HWwkb        HWef        TLwkb',
     *'        TLef',/)
      WRITE(6,1002) L,RBL,VBL,HWWKB,HWEF,TLWKB,TLEF
 1002 FORMAT(1x,i4,1x,f5.2,1x,f6.2,4(1x,e12.5))
      RETURN
      END

C*****************************************************************************
C*    SUBROUTINE COU: CALCULA AS FUNCOES DE ONDA COULOMBIANAS REGULAR(F) 
C*    E IRREGULAR(G) E TAMBEM AS DEFASAGENS COULOMBIANAS (S)
C*****************************************************************************
      SUBROUTINE COU(R,RP,E,L,F,FP,G,GP,S)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(2001),FP(2001),G(2001),GP(2001),S(2001)
      H=0.1d0
      TE=2.d0*E
      TF=E**2
      LL=L 
      IF(LL-100)20,35,35
   20 ELP=100.d0
      J=100
      GO TO 45
   35 ELP=LL
      J=LL
   45 A=DATAN (E/ELP)
      B=DSQRT (TF+ELP**2)
      Y=A*(ELP-0.5d0)+E*(DLOG(B)-1.d0)-DSIN (A)/(12.d0*B)
     1+DSIN (3.d0*A)/(360.d0*B**3)-DSIN (5.d0*A)/(1260.d0*B**5)
     2+DSIN (7.d0*A)/(1680.d0*B**7)-DSIN (9.d0*A)/(1188.d0*B**9)
      K=J-1
      IF(J-LL)65,65,70
   65 S(J)=Y
   70 DO 100 I=1,K
      ELP=ELP-1.d0
      J=J-1
      Y=Y-DATAN(E/ELP)
      IF(J-LL)95,95,100
   95 S(J)=Y
  100 CONTINUE
      DEL1=R-TE
      RMAX=.41666667d0*(TF+4.d0*E+3.d0)
      IF(RMAX.LT.10.0d0) RMAX=10.0d0
      DEL=R-RMAX 
      IF(E-5.d0)280,130,130
  130 IF(DABS(DEL1)-DABS (DEL))140,140,280
  140 DEL=DEL1
      IF(DEL)147,145,147
  145 I=2
      GO TO 150
  147 I=1
  150 X=TE
      T1=TF
      T2=T1**2
      T3=E**(2.d0/3.d0)
      T4=T3**2
      T5=T4**2
      T6=T3*T5
      T7=T4*T6
      T8=T3*T7
      T9=E**(1.d0/6.d0)
      Y=1.22340402d0*T9*(1.d0+.495957017d-1/T4-.888888889d-2/T1
     *+.245519918d-2/T6-.910895806d-3/T2+.253468412d-3/T8)
      Z=-.707881773d0/T9*(1.d0-.172826039d0/T3+.317460317d-3/T1
     *-.358121485d-2/T5+.311782468d-3/T2-.907396643d-3/T7)
      GO TO 665
  280 IF(E)285,290,285
  285 IF(DEL)310,290,290
  290 X=R
      I=2
      GO TO 320
  310 X=RMAX
      I=1
  320 T1=TF
      T2=2.d0*X
      T3=X-E*DLOG(T2)+S(1)
      T4=E/T2
      SS=1.d0
      TS=0.d0
      SL=0.d0
      TL=1.d0-E/X
      SSS=1.d0
      STS=0.d0
      SSL=0.d0
      STL=TL
      EN=0.d0
      DO 620 K=1,25
      T5=EN+1.d0
      T6=T5+EN
      T7=EN*T5
      T8=T6*T4/T5
      T9=(T1-T7)/(T2*T5)
      T5=T8*SS-T9*TS
      TS=T8*TS+T9*SS
      SS=T5
      IF(DABS(SS/SSS)-1.0d-10)630,630,540
  540 T5=T8*SL-T9*TL-SS/X
      TL=T8*TL+T9*SL-TS/X
      SL=T5
      SSS=SSS+SS
      STS=STS+TS
      SSL=SSL+SL
      STL=STL+TL
      EN=EN+1.d0
  620 CONTINUE
  630 T8=DSIN(T3)
      T9=DCOS(T3)
      Y=SSS*T9-STS*T8
      Z=SSL*T9-STL*T8
  665 GO TO (670,810),I
  670 M=1
  671 N=DABS(DEL/H)+1.0d0 
      DX=DEL/DFLOAT(N)
      T1=DX/2.d0
      T2=DX/8.d0
      T3=TE
      DO 805 I=1,N
      T4=DX*(T3/X-1.d0)*Y
      X=X+T1
      T5=DX*(T3/X-1.d0)*(Y+T1*Z+T2*T4)
      X=X+T1
      T6=DX*(T3/X-1.d0)*(Y+DX*Z+T1*T5)
      Y=Y+DX*(Z+(T4+2.d0*T5)/6.d0)
      Z=Z+(T4+4.d0*T5+T6)/6.d0
  805 CONTINUE
      GO TO (810,828),M
  810 G(1)=Y
      M=2
      DEL=RP-R
      W=Z
      GO TO 671
  828 GP(1)=Y
      T1=TF
      T8=DSQRT(1.d0+T1)
      G(2)=((1.d0/R+E)*G(1)-W)/T8
      GP(2)=((1.d0/RP+E)*Y-Z)/T8
      T2=1.d0
      T3=2.d0
      DO 910 I=3,LL
      T4=T2+T3
      T5=T2*T3
      T6=T3*DSQRT (T2**2+T1)
      T7=T2*DSQRT (T3**2+T1)
      G(I)=(T4*(E+T5/R )*G(I-1)-T6*G(I-2))/T7
      GP(I)=(T4*(E+T5/RP)*GP(I-1)-T6*GP(I-2))/T7
      T2=T2+1.d0
      T3=T3+1.d0
  910 CONTINUE
      I=L+11
      N=2.d0*R+11.d0
      IF(I-N)960,960,950
  950 N=I
  960 Y=1.0d-20
      YP=Y
      X=Y
      XP=X
      Z=0.d0
      ZP=Z
      T2=N 
 1000 T3=T2+1.d0
      T4=T2+T3
      T5=T2*T3
      T6=T2*DSQRT(T3**2+T1)
      T7=T3*DSQRT(T2**2+T1)
      Y =(T4*(E+T5/R )*Y -T6*Z )/T7
      YP=(T4*(E+T5/RP)*YP-T6*ZP)/T7
      IF(N-LL)1060,1060,1080
 1060 F(N)=Y
      FP(N)=YP
      GO TO 1120
 1080 IF(1.d0-DABS (Y))1090,1090,1120
 1090 CONTINUE 
      Y =Y *1.0d-20 
      YP=YP*1.0d-20
      X =X *1.0d-20
      XP=XP*1.0d-20
 1120 N=N-1
      Z=X
      ZP=XP
      X=Y
      XP=YP
      T2=T2-1.d0
      IF(N)1150,1150,1000
 1150 Y=F(1)*G(2)-F(2)*G(1)
      YP=FP(1)*GP(2)-FP(2)*GP(1)
      Z=1.d0/(Y*T8)
      ZP=1.d0/(YP*T8)
      DO 1180 I=1,LL
      FP(I)=FP(I)*ZP
 1180 F(I)=F(I)*Z
      RETURN
      END

C*****************************************************************************
C*    SUBROUTINE PLEGO: CALCULA OS POLINOMIOS DE LEGENDRE
C*
C*****************************************************************************
      SUBROUTINE PLEGO(TETA,PL)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PL(2001)
      PA=1.0d0
      PL(1)=1.0d0
      X=DCOS(TETA)
      PB=X
      PL(2)=X
      DO 10 N=3,2001
      R=DFLOAT(N-2)
      PC=X*PB*(2.d0*R+1.d0)/(R+1.d0)-PA*R/(R+1.d0)
      PL(N)=PC
      PA=PB
      PB=PC
 10   CONTINUE
      RETURN
      END	
