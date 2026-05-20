!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!   WRINKLING
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,par,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), par(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      double precision w, wx, wxx, wxxx, r, r2, t
      double precision Delt, Q0, Q1, Q2, Q3
      double precision n, n2, tau, Bo, ep, L

      ep = par(1)
      n = par(2)
      tau = par(3)
      Bo = par(4)
      L = par(5)

      w = U(1)
      wx = U(2)
      wxx = U(3)
      wxxx = U(4)
      t = U(5)
      r = 1.0 + ep**0.25*L*t

      r2 = r**2;
      n2 = n**2;

      Delt = (tau - 1.0)/r2;
      Q0 = Bo + n2/r2*(ep*(n2-4.0)/r2 + (1.0 - Delt))
      Q1 = 1.0/r*(ep*(2.0*n2 + 1.0)/r2 - (1.0 - Delt))/ep**0.25;
      Q2 = -(ep*(2.0*n2 + 1.0)/r2 + (1.0 + Delt))/ep**0.5;
      Q3 = 2.0/r/ep**0.75;

      F(1) = L*wx
      F(2) = L*wxx
      F(3) = L*wxxx
      F(4) = -L*(Q3*wxxx + Q2*wxx + Q1*wx + Q0*w)
      F(5) = 1.0d0

      !print *, U
      !print *, F
      !print *, Q3
      !print *, Q2
      !print *, L
      !print *, n
      !print *, ep
      !call sleep(5)

      END SUBROUTINE FUNC
!---------------------------------------------------------------------- 

      SUBROUTINE STPNT(NDIM,U,PAR,T) 
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),par(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      par(1) = 0.17 
      par(2) = 8.0
      par(3) = 20.0
      par(4) = 0.0
      par(5) = 100.0
      par(6) = 1.238e-4 ! norm

      !U(1)=0.0d0
      !U(2)=0.0d0
      !U(3)=0.0d0
      !U(4)=0.0d0
      !U(5)=T

      END SUBROUTINE STPNT
!---------------------------------------------------------------------- 

      subroutine bcnd(ndim,par,icp,nbc,u0,u1,fb,ijac,dbc)

      implicit none
      integer, intent(in) :: ndim, icp(*), nbc, ijac
      double precision, intent(in) :: par(*), u0(ndim), u1(ndim)
      double precision, intent(out) :: fb(nbc)
      double precision, intent(inout) :: dbc(nbc,*)

      ! Boundary conditions
      fb(1) = u0(1)
      fb(2) = u0(2)
      fb(3) = u1(1)
      fb(4) = u1(2)
      fb(5) = u0(5) ! For t starts at 0
      !fb(5) = u0(3) - 1.0 ! Curvature

      end subroutine bcnd

      SUBROUTINE ICND2 
      END SUBROUTINE ICND2

      SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
      !--------- ----

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
        DOUBLE PRECISION, INTENT(IN) :: PAR(*)
        DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
        DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
        DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)

        FI(1) = U(1)*U(1)-PAR(6) 

      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
