!------------------------------------------------------------------------------!
!> Calculation of toroidal functions \f$P_{n-1/2}^m \left(z\right)\f$ and 
!! \f$Q_{n-1/2}^m \left(z\right)\f$.
!------------------------------------------------------------------------------!
!> \note Copied and adapted from the DTORH1 routine by Segura and Gil 
!! \cite Segura2000 .
!------------------------------------------------------------------------------!
module dtorh
#include <PB3D_macros.h>
    use num_vars, only: dp, pi, max_str_ln
    use messages
    
    implicit none
    private
    public dtorh1
    
   character(len=max_str_ln) :: err_msg                                         !< error message

contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  INPUT :                                                         !
    !    Z        ARGUMENT OF THE FUNCTIONS                            !
    !    M        ORDER OF THE FUNCTIONS                               !
    !    NMAX     MAXIMUM DEGREE OF THE FUNCTIONS :                    !
    !             WE GET  FUNCTIONS OF ALL THE ORDERS BELOW            !
    !             MIN(NEWN,NMAX). NEWN IS DEFINED BELOW .              !
    !  OUTPUT :                                                        !
    !   *IF MODE IS EQUAL TO 0:                                        !
    !    PL(N)                                                         !
    !             THESE VALUES ARE KEPT IN AN ARRAY                    !
    !    QL(N)                                                         !
    !             THESE VALUES ARE KEPT IN AN ARRAY                    !
    !                                                                  !
    !    NEWN     MAXIMUM ORDER OF FUNCTIONS CALCULATED WHEN           !
    !             PL (NMAX+1) IS LARGER THAN 1/TINY                    !
    !             (OVERFLOW LIMIT = 1/TINY, TINY IS DEFINED BELOW)     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Calculates toroidal harmonics of a fixed order \c m and degrees
    !! up to \c nmax.
    !!
    !! Optionally, the \c mode can be specified to be different from 0 [default]
    !!  - if \c mode=1:
    !!      - The set of functions evaluated is:
    !!        \f[\frac{P_{n-1/2}^m \left(z\right)}{\Gamma(m+1/2)} \ ,\quad 
    !!          \frac{Q_{n-1/2}^m \left(z\right)}{\Gamma(m+1/2)},\f]
    !!        which are respectively stored in the arrays \c pl(n), \c ql(n).
    !!      - newn refers to this new set of functions.
    !!      - Note 1. and note 2. also apply in this case.
    !!  - if \c mode=2:
    !!      - The code performs as for mode 1, but the restriction \f$ z<20 \f$.
    !!      - For the evaluation of the continued fraction is not considered.
    !!
    !! Also, the parameter \c ipre can be used to select a different precision:
    !!  - For \c ipre=1, the precision is \f$10^{-12}\f$, taking
    !!    \f$\epsilon<10^{-12}\f$
    !!  - For \c ipre=2, the precision is \f$10^{-8}\f$, taking
    !!    \f$\epsilon<10^{-8}\f$
    !!
    !! where \f$\epsilon\f$ controls the accuracy.
    !!
    !! \warning Use \c mode 2 only if high \c m 's for \f$z>20\f$ are required.
    !! The evaluation of the cf may fail to converge for too high \c z 's.
    !!
    !! \note
    !!    -# For a precision of \f$10^{-12}\f$, if \f$z>5\f$ and \f$z/m>0.22\f$,
    !!    the code uses a series expansion for \c pl(0).\n
    !!    When \f$z<20\f$ and \f$z/m<0.22\f$ a continued fraction is applied.
    !!    -# For a precision of \f$10^{-8}\f$, if \f$z>5\f$ and \f$z/m>0.12\f$,
    !!    the code uses a series expansion for \c pl(0).\n
    !!    When \f$z<20\f$ and \f$z/m<0.12\f$ a continued fraction is applied.
    !! 
    !! \return ierr
    integer function dtorh1(z,m,nmax,pl,ql,newn,mode,ipre) result(ierr)
        character(*), parameter :: rout_name = 'DTORH1'
        
        ! input / output
        real(dp), intent(in) :: z                                               !< (real) point at which toroidal harmonics are evaluated
        integer, intent(in) :: m                                                !< order of toroidal harmonics (\c m>0)
        integer, intent(in) :: nmax                                             !< maximum degree of toroidal harmonics (\c nmax>0)
        real(dp), intent(inout) :: pl(0:nmax)                                   !< toroidal harmonics of first kind \f$P_{n-1/2}^m\left(z\right)\f$
        real(dp), intent(inout) :: ql(0:nmax)                                   !< toroidal harmonics of second kind \f$Q_{n-1/2}^m\left(z\right)\f$
        integer, intent(inout) :: newn                                          !< maximum order of functions calculated when <tt> pl(nmax+1)>1/tiny </tt> with <tt>tiny=</tt>\f$10^{-290}\f$
        integer, intent(in), optional :: mode                                   !< mode that controls output
        integer, intent(in), optional :: ipre                                   !< precision
        
        ! local variables
        INTEGER NMAXP,MP,NP,N,ICAL,I
        INTEGER :: mode_loc, ipre_loc
        REAL(dp) OVER,TINYSQ,QZ,PISQ,FL,CC,AR,GAMMA,FC,QDC1,QARGU,&
            &ARGU1,DFACQS,FCP,DD,D1,QM0,DFAC3,DFAC4,PL0
        real(dp), allocatable :: PL_loc(:), QL_loc(:)                           ! local copy of pl and ql, possibly larger
        !REAL(dp) DPPI                                                           ! Is never used?
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! THE DIMENSION OF QLMM (INTERNAL ARRAY) MUST BE GREATER THAN M    !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL(dp)  QLMM(0:1001),PR(2)
        real(dp), PARAMETER :: EPS=1.E-14_dp, TINY=1.E-290_dp
        
        ! initialize ierr
        ierr = 0
        
        ! test
        if (NMAX.lt.0) then
            ierr = 1
            err_msg = 'Lowest possible degree is 0. Use recursive properties &
                &of toroidal functions if you need lower.'
            CHCKERR(err_msg)
        end if
        
        ! set optional parameters
        mode_loc = 0
        if (present(mode)) mode_loc = mode
        ipre_loc = 1
        if (present(ipre)) ipre_loc = ipre
        
        OVER=1._dp/TINY
        TINYSQ=DSQRT(TINY)
        IF ((ipre_loc.NE.1).AND.(ipre_loc.NE.2)) THEN
            ierr = 1
            err_msg = 'IPRE MUST BE 1 OR 2'
            CHCKERR(err_msg)
        END IF
        
        PR(1)=.22_dp
        PR(2)=.12_dp
        NMAXP=NMAX-1                                                            ! decremented by one to get more logical ouptut
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   EPS: REQUIRED ACCURACY FOR THE CONTINUED FRACTION     !
        !        (MODIFIED LENTZ)                                 !
        !   TINY: SMALL PARAMETER TO PREVENT OVERFLOWS IN THE CF  !
        !          (CLOSE TO THE UNDERFLOW LIMIT)                 !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF (Z.LE.1._dp) THEN
            ierr = 1
            err_msg = 'IMPROPER ARGUMENT. Z MUST BE GREATER THAN 1'
            CHCKERR(err_msg)
        END IF
        QZ=Z
        PISQ=DSQRT(PI)
        !DPPI=DSQRT(2._dp)/PISQ                                                  ! Is never used?
        FL=M/2._dp
        CC=ABS(FLOAT(INT(FL))-FL)
        IF (CC.LT.0.4_dp) THEN
            AR=1.
        ELSE
            AR=-1.
        END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   WE CHOOSE EXPANSION OR CF FOR PL(0) DEPENDING ON THE VALUES !
        !   OF Z,M AND MODE                                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ICAL=1
        IF (M.NE.0) THEN
            IF ((Z/M).GT.PR(ipre_loc)) ICAL=2
            IF (Z.LT.5._dp) THEN
                ICAL=1
            ELSE IF (Z.GT.20.D0) THEN
                IF ((mode_loc.NE.2).AND.(ICAL.EQ.1)) ICAL=0
            END IF
            IF (ICAL.EQ.0) THEN
                ierr =1
                CHCKERR('YOU MUST CHOOSE MODE=2')
            END IF
        ELSE
            IF (Z.LT.5.D0) THEN
              ICAL=1
            ELSE
              ICAL=2
            END IF
        END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   WE USE THE CODE IF NMAX IS GREATER THAN OR EQUAL TO 2  !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(NMAXP.LT.2) NMAXP=2
        
        ! set up local PL, QL
        allocate(PL_loc(0:nmaxp+1))
        allocate(QL_loc(0:nmaxp+1))
        
        IF (mode_loc.EQ.0) THEN
            GAMMA=GAMMAH(M,OVER)*AR*PISQ  
            IF (ABS(GAMMA).LT.TINY) THEN
                ierr =1
                err_msg = 'M IS TOO LARGE FOR MODE=0, BETTER TRY MODE=1'
                CHCKERR(err_msg)
            END IF              
        ELSE 
            GAMMA=AR
        END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   WE EVALUATE THE CONTINUED FRACTION  USING     !
        !   LENTZ-THOMPSON                                !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ierr = FRAC(Z,M,0,EPS,TINYSQ,FC)
        CHCKERR('')
        QDC1=QZ*QZ-1._dp
        QARGU=QZ/DSQRT(QDC1)
        !DFAC1=DPPI*GAMMA/PI                                                     ! Is never used?
        !DFAC2=GAMMA/DPPI                                                        ! Is never used?
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   WE EVALUATE Q_{-1/2},Q^{1}_{-1/2}               !
        !   USING SLATEC ROUTINES FOR ELLIPTIC FUNCTIONS    !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ARGU1=DSQRT(2._dp/(Z+1._dp))
        QLMM(0)=ARGU1*ELLIP1(ARGU1)
        QLMM(1)=-1._dp/DSQRT(2._dp*(QZ-1._dp))*ELLIP2(ARGU1)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   WE APPLY FORWARD RECURRENCE IN M FOR Q'S  !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        MP=1
        IF (mode_loc.EQ.0) THEN
          1 IF ((MP.LE.M).AND.(ABS(QLMM(MP)).LT.OVER)) THEN
                QLMM(MP+1)=-2._dp*MP*QARGU*QLMM(MP)-&
                    &(MP-0.5_dp)*(MP-0.5_dp)*QLMM(MP-1) 
                MP=MP+1
                GOTO 1 
            ENDIF       
            IF ((MP-1).NE.M) THEN
                ierr = 1
                err_msg = 'M IS TOO LARGE FOR MODE=0, BETTER TRY MODE=1'
                CHCKERR(err_msg)
            END IF              
        ELSE
            QLMM(0)=QLMM(0)/PISQ
            QLMM(1)=QLMM(1)*2._dp/PISQ
          2 IF ((MP.LE.M).AND.(ABS(QLMM(MP)).LT.OVER)) THEN
                D1=MP+0.5_dp
                QLMM(MP+1)=-2._dp*MP*QARGU*QLMM(MP)/D1&
                    &-(MP-0.5_dp)*QLMM(MP-1)/D1               
                MP=MP+1
                GOTO 2
            ENDIF
            IF ((MP-1).NE.M) THEN
                ierr = 1
                err_msg = 'M IS TOO LARGE FOR MODE=1,2'
                CHCKERR(err_msg)
            END IF              
        END IF
        !NMMAX=M                                                                 ! Is never used?
        DFACQS=-(GAMMA/QLMM(M+1))*GAMMA/PI
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   EVALUATION OF PL(0)                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF (ICAL.EQ.1) THEN
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !   WE CALCULATE THE CF FOR P'S                     !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ierr = FRACPS(QARGU,M+1,0,EPS,TINYSQ,FCP) 
            CHCKERR('')
            DD=M+0.5_dp
            IF (mode_loc.NE.0) THEN
                FCP=FCP/DD
                DFACQS=DFACQS/DD
            END IF
            PL_loc(0)=DFACQS/DSQRT(QDC1)/(1._dp-FCP*QLMM(M)/QLMM(M+1))
        ELSE
            CALL EXPAN(Z,mode_loc,ipre_loc,OVER,QARGU,M,PL0)
            PL_loc(0)=PL0
        END IF
        QM0=QLMM(M)
        DFAC3=(GAMMA/QM0)*GAMMA/PI/(0.5_dp-M)
        PL_loc(1)=PL_loc(0)*FC+DFAC3
        NP=1 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   WE USE THE RECURRENCE RELATIONS FOR P'S   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      3 IF ((NP.LE.NMAXP).AND.(ABS((NP-M+0.5_dp)*PL_loc(NP)).LT.OVER)) THEN
            PL_loc(NP+1)=(2._dp*NP*Z*PL_loc(NP)-(NP+M-0.5_dp)*PL_loc(NP-1))/&
                &(NP-M+0.5_dp)
            NP=NP+1
            GOTO 3
        ENDIF      
        NMAXP=NP-1
        NEWN=NMAXP+1                                                            ! incremented by 1 to get more logical output
        DFAC4=(FACTCO(NMAXP,PL_loc(NMAXP+1),M)*GAMMA)*GAMMA/PI/(NMAXP+M+0.5_dp)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   WE EVALUATE THE CONTINUED FRACTION USING LENTZ-THOMPSON  !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ierr = FRAC(Z,M,NMAXP,EPS,TINYSQ,FC)         
        CHCKERR('')
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   EVALUATION OF PL(NMAX+1) AND PL(NMAX) USING   !
        !   THE WRONSKIAN W{PL(NMAX),QL(NMAX)},           !
        !   THE KNOWN VALUES OF PL(NMAX+1) AND PL(NMAX)   !
        !   THE VALUE OF H = QL(NMAX+1)/QL(NMAX)          !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        QL_loc(NMAXP)=DFAC4/(1._dp-FC*PL_loc(NMAXP)/PL_loc(NMAXP+1))
        QL_loc(NMAXP+1)=QL_loc(NMAXP)*FC
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   WE USE THE BACKWARD RECURRENCE RELATION   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO 4 I=1,NMAXP
            NP=NMAXP+1-I
            N=NP-1
            QL_loc(N)=((NP+NP)*Z*QL_loc(NP)-(NP-M+0.5_dp)*QL_loc(NP+1))/&
                &(NP+M-0.5_dp)
      4 CONTINUE
        
        ! copy local PL and QL back
        PL = PL_loc(0:NMAX)
        QL = QL_loc(0:NMAX)
    end function DTORH1
    
    integer function FRAC(Z,M,NMAX,EPS,TINYSQ,FC) result(ierr)
        character(*), parameter :: rout_name = 'FRAC'
        
        INTEGER M,NMAX,N,MM          
        REAL(dp) Z,EPS,TINYSQ,FC,DZ2,DN0,DN1,DN2,DN3,DN4,B,A,C0,D0,DELTA
        
        ! initialize ierr
        ierr = 0
        
        N=NMAX
        MM=0
        DZ2=Z+Z
        DN0=N+M
        DN1=N+1._dp
        DN2=DN0+0.5_dp
        DN3=DN0-0.5_dp
        DN4=1._dp-M-M
        B=2._dp*DN1*Z/DN2 
        A=1._dp
        FC=TINYSQ
        C0=FC
        D0=0._dp
     81 D0=B+A*D0
        IF (DABS(D0).LT.TINYSQ) D0=TINYSQ
        C0=B+A/C0
        IF (DABS(C0).LT.TINYSQ) C0=TINYSQ
        D0=1._dp/D0
        DELTA=C0*D0
        FC=FC*DELTA
        MM=MM+1
        A=-(1.D0+DN4/(DN3+MM))  
        B=DZ2*(DN1+MM)/(DN2+MM)
        IF (MM.LT.10000) THEN
            IF (DABS(DELTA-1.D0).GT.EPS) GOTO 81
        END IF
        IF (MM.EQ.10000) then
            ierr =1
            err_msg = 'CF CONVERGENCE FAILS'
            CHCKERR(err_msg)
        END IF
    end function FRAC

    integer function FRACPS(QZ,M,N,EPS,TINYSQ,FC) result(ierr)
        character(*), parameter :: rout_name = 'FRACPS'
        
        INTEGER M,N,MM         
        REAL(dp) QZ,EPS,TINYSQ,FC,DN2,DZ2,DM,B,A,C0,D0,DELTA
        
        ! initialize ierr
        ierr = 0
        
        MM=0
        DN2=N*N
        DZ2=QZ+QZ
        DM=M-0.5_dp
        B=DZ2*M
        A=DN2-DM*DM
        FC=TINYSQ
        C0=FC
        D0=0._dp
     82 D0=B+A*D0
        IF (DABS(D0).LT.TINYSQ) D0=TINYSQ
        C0=B+A/C0
        IF (DABS(C0).LT.TINYSQ) C0=TINYSQ
        D0=1._dp/D0
        DELTA=C0*D0
        FC=FC*DELTA
        MM=MM+1
        A=DN2-(MM+DM)*(MM+DM)  
        B=DZ2*(MM+M)
        IF (MM.LT.10000) THEN
            IF (ABS(DELTA-1._dp).GT.EPS) GOTO 82
        END IF
        IF (MM.EQ.10000) then
            ierr =1
            err_msg= 'CF CONVERGENCE FAILS'
            CHCKERR(err_msg)
        END IF
    end function FRACPS

    REAL(dp) FUNCTION ELLIP2(AK)
        REAL(dp) AK
        REAL(dp) Q
        
        Q=(1._dp-AK)*(1._dp+AK)
        ELLIP2=DRF(Q)-(AK)**2*DRD(Q)/3.D0
    END FUNCTION ELLIP2
    
    REAL(dp) FUNCTION ELLIP1(AK)
        REAL(dp) AK
        
        ELLIP1=DRF((1._dp-AK)*(1._dp+AK))
    END FUNCTION ELLIP1
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  WE USE SLATEC ROUTINES IN THE EVALUATION OF ELLIPTIC INTEGRALS   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    REAL(dp) FUNCTION DRF (Y)
        REAL(dp) EPSLON, ERRTOL
        REAL(dp) C1, C2, C3, E2, E3, LAMDA
        REAL(dp) MU, S, XN, XNDEV
        REAL(dp) XNROOT, Y, YN, YNDEV, YNROOT, ZN, ZNDEV, ZNROOT
        LOGICAL FIRST
        SAVE ERRTOL,C1,C2,C3,FIRST
        DATA FIRST /.TRUE./
        
        !***FIRST EXECUTABLE STATEMENT  DRF
        IF (FIRST) THEN
            ERRTOL = 1.E-8_dp
            C1 = 1.0_dp/24.0_dp
            C2 = 3.0_dp/44.0_dp
            C3 = 1.0_dp/14.0_dp
        ENDIF
        FIRST = .FALSE.
        DRF = 0.0_dp
        XN = 0._dp
        YN = Y
        ZN = 1._dp
        
     30 MU = (XN+YN+ZN)/3.0_dp
        XNDEV = 2.0_dp - (MU+XN)/MU
        YNDEV = 2.0_dp - (MU+YN)/MU
        ZNDEV = 2.0_dp - (MU+ZN)/MU
        EPSLON = MAX(ABS(XNDEV),ABS(YNDEV),ABS(ZNDEV))
        IF (EPSLON.LT.ERRTOL) GO TO 40
        XNROOT = SQRT(XN)
        YNROOT = SQRT(YN)
        ZNROOT = SQRT(ZN)
        LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
        XN = (XN+LAMDA)*0.250_dp
        YN = (YN+LAMDA)*0.250_dp
        ZN = (ZN+LAMDA)*0.250_dp
        GO TO 30
        
     40 E2 = XNDEV*YNDEV - ZNDEV*ZNDEV
        E3 = XNDEV*YNDEV*ZNDEV
        S  = 1.0_dp + (C1*E2-0.10_dp-C2*E3)*E2 + C3*E3
        DRF = S/SQRT(MU)
    END FUNCTION DRF


    REAL(dp) FUNCTION DRD (Y)
        REAL(dp) EPSLON, ERRTOL
        REAL(dp) C1, C2, C3, C4, EA, EB, EC, ED, EF, LAMDA
        REAL(dp) MU, POWER4, SIGMA, S1, S2, XN, XNDEV
        REAL(dp) XNROOT, Y, YN, YNDEV, YNROOT,  ZN, ZNDEV, ZNROOT    
        LOGICAL FIRST
        SAVE ERRTOL, C1, C2, C3, C4, FIRST
        DATA FIRST /.TRUE./
        
        !***FIRST EXECUTABLE STATEMENT  DRD
        IF (FIRST) THEN
            ERRTOL = 1.E-8_dp
            C1 = 3.0_dp/14.0_dp
            C2 = 1.0_dp/6.0_dp
            C3 = 9.0_dp/22.0_dp
            C4 = 3.0_dp/26.0_dp
        ENDIF
        FIRST = .FALSE.
        
        ! CALL ERROR HANDLER IF NECESSARY.
        
        DRD = 0.0_dp
        
        XN = 0._dp
        YN = Y
        ZN = 1._dp
        SIGMA = 0.0_dp
        POWER4 = 1.0_dp
        
     30 MU = (XN+YN+3.0_dp*ZN)*0.20_dp
        XNDEV = (MU-XN)/MU
        YNDEV = (MU-YN)/MU
        ZNDEV = (MU-ZN)/MU
        EPSLON = MAX(ABS(XNDEV), ABS(YNDEV), ABS(ZNDEV))
        IF (EPSLON.LT.ERRTOL) GO TO 40
        XNROOT = SQRT(XN)
        YNROOT = SQRT(YN)
        ZNROOT = SQRT(ZN)
        LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
        SIGMA = SIGMA + POWER4/(ZNROOT*(ZN+LAMDA))
        POWER4 = POWER4*0.250_dp
        XN = (XN+LAMDA)*0.250_dp
        YN = (YN+LAMDA)*0.250_dp
        ZN = (ZN+LAMDA)*0.250_dp
        GO TO 30
        
     40 EA = XNDEV*YNDEV
        EB = ZNDEV*ZNDEV
        EC = EA - EB
        ED = EA - 6.0_dp*EB
        EF = ED + EC + EC
        S1 = ED*(-C1+0.250_dp*C3*ED-1.50_dp*C4*ZNDEV*EF)
        S2 = ZNDEV*(C2*EF+ZNDEV*(-C3*EC+ZNDEV*C4*EA))
        DRD = 3.0_dp*SIGMA + POWER4*(1.0_dp+S1+S2)/(MU*SQRT(MU))
    END FUNCTION DRD

    REAL(dp) FUNCTION GAMMAH(N,OVER)
        INTEGER N,I,J
        REAL(dp) OVER
        
        I=N
        J=2*I-1
        GAMMAH=1._dp
     85 IF ((J.GE.1).AND.(GAMMAH.LT.OVER)) THEN
            GAMMAH=GAMMAH*J/2._dp
            I=I-1
            J=2*I-1
            GOTO 85
        ENDIF  
        IF (J.GT.1) THEN 
            GAMMAH=0._dp
        END IF   
    END FUNCTION GAMMAH

    REAL(dp) FUNCTION FACTCO(N,PL,M)
        INTEGER N,M,J  
        REAL(dp) PL,X1,X2 
        
        FACTCO=1._dp/PL
        X1=M+0.5_dp
        X2=-M+0.5_dp
        J=N
    86  IF (J.GE.0) THEN
            FACTCO=FACTCO*(J+X1)/(J+X2)
            J=J-1
            GOTO 86
        ENDIF
    END FUNCTION FACTCO

    SUBROUTINE EXPAN(Z,MODE,IPRE,OVER,QARGU,M,PL)
        INTEGER MODE,IPRE,M,K,I 
        REAL(dp) Z,OVER,QARGU,PL,PISQ,DB,FL,CC,DZ,GAMMA,AR,&
            &DFAC,DF1,A0,Z2I,DA1,DA2,DELTA,SUM,DCC,DCCP,DKK
        REAL(dp) PRECI(2)
        
        PRECI(1)=1.E-13_dp
        PRECI(2)=1.E-9_dp
        PISQ=DSQRT(PI)
        DB=2._dp*DLOG(2._dp)
        FL=M/2.
        CC=ABS(FLOAT(INT(FL))-FL)
        IF (CC.LT.0.4_dp) THEN
            AR=1.
        ELSE
            AR=-1.
        END IF      
        DZ=1._dp
        DO 100 I=1,M
            DZ=DZ/QARGU
    100 CONTINUE
        IF (MODE.EQ.0) THEN
            GAMMA=GAMMAH(M,OVER)*AR*PISQ 
        ELSE 
            GAMMA=AR
        END IF
        DFAC=2._dp/PI*DZ*GAMMA/PISQ
        DF1=DLOG(2._dp*Z)
        A0=1._dp/DSQRT(2._dp*Z)
        Z2I=1._dp/(Z*Z)
        DELTA=1._dp
        SUM=0._dp
        K=0
        DA2=FACTOR(M,0)
        DA1=DB+PSI(M,0)
     87 DELTA=(DF1+DA1)*DA2*A0
        SUM=SUM+DELTA
        DCC=0.5_dp+M+2._dp*K
        DCCP=DCC+1._dp
        DKK=K+1._dp
        DA2=DA2*DCCP*DCC/(DKK*DKK)
        DA1=DA1+1._dp/DKK-1._dp/DCCP-1._dp/DCC
        K=K+1
        A0=A0*0.25_dp*Z2I
        IF (ABS(DELTA/SUM).GT.PRECI(IPRE)) GOTO 87
        PL=SUM*DFAC
    END SUBROUTINE EXPAN

    REAL(dp) FUNCTION FACTOR(M,K)
        INTEGER M,K,N,I,J
        REAL(dp) X1
        
        FACTOR=1._dp
        IF (K.GE.1) THEN
        X1=M+0.5_dp
        N=2*K-1
        I=K
        J=N
     88 IF (J.GE.0) THEN
            FACTOR=FACTOR*(J+X1)/I/I
            J=J-1
            I=I-1
            IF (I.EQ.0) I=1
                GOTO 88
            ENDIF 
        END IF  
    END FUNCTION FACTOR

    REAL(dp) FUNCTION PSI(M,K)
        INTEGER M,K,N,J,I
        REAL(dp) FACTR1,FACTR2
        PSI=0._dp
        FACTR1=0._dp
        FACTR2=0._dp
        N=2*K+M
        J=1
        IF (K.GE.1) THEN
         89 IF (J.LE.K) THEN
                FACTR1=FACTR1+1._dp/J
                J=J+1
                GOTO 89
            ENDIF
        END IF
        I=1
     90 IF (I.LE.N) THEN
            FACTR2=FACTR2+1._dp/(2._dp*I-1._dp)
            I=I+1
            GOTO 90
        ENDIF
        IF (N.EQ.0) FACTR2=0._dp
        PSI=FACTR1-2._dp*FACTR2
    END FUNCTION PSI
end module
