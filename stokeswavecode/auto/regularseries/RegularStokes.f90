!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   Regular series method for Stokes wave 
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

! Evaluates the algebraic equations or ODE right hand side

! Input arguments :
!      NDIM   :   Dimension of the ODE system 
!      U      :   State variables
!      ICP    :   Array indicating the free parameter(s)
!      PAR    :   Equation parameters

! Values to be returned :
!      F      :   ODE right hand side values

! Normally unused Jacobian arguments : IJAC, DFDU, DFDP (see manual)

      IMPLICIT NONE
      INTEGER NDIM, IJAC, ICP(*)
      DOUBLE PRECISION F(NDIM), U(NDIM), PAR(*), DFDU(*), DFDP(*)
      integer :: N
      integer :: i, j, k
      double precision, parameter :: pi = 4.0d0*atan(1.0d0)
      
      complex(kind(1.0D0)) :: dzdf, y, deltabeta0, deltabetaend
      double precision :: phi, absvalf, ep
      double precision :: mu, B
      ! Note that we only use NDIM - 2 components in a
      double precision :: a(NDIM)
      double precision :: W0, Wend

      N = NDIM - 1
      mu = U(N)
      B = U(N + 1)
      ep = par(1)
      
      par(2) = mu
      par(3) = B


      do j = 1, N-1
        a(j) = U(j)
      end do

      do j = 1, N
        phi = 0.5d0*(j-1.0d0)/(N-1.0d0)

        !This is delta + ibeta (without the first term)
        dzdf = 1
        do k = 1, N - 1
            dzdf = dzdf + a(k)*exp(-2.0d0*k*(0,1.0d0)*pi*phi)
        end do

        !The integral term (gravity)
        y = 0.0;
        do k = 1, N - 1
          y = y + 1.0d0/(2.0d0*k*pi)*(cos(2.0d0*k*pi*phi)-1.0d0)*a(k)
        end do

        absvalf = (abs(dzdf))**2
        F(j) = 2.0d0*pi/mu*y + 1.0d0/(2.0d0*absvalf) - B
      end do

      deltabeta0 = 1.0d0
      deltabetaend = 1.0d0
      do k = 1, N - 1
        deltabeta0 = deltabeta0 + a(k)
        deltabetaend = deltabetaend + a(k)*exp(-pi*(0,1.0d0)*k)
      end do
      W0 = 1.0d0/(deltabeta0)
      Wend = 1.0d0/(deltabetaend)

      F(N+1) = 1.0d0 - abs(W0)**2*abs(Wend)**2 - ep**2

      !write(*,*) 'MY M1 = ', M1
      !write(*,*) 'MY M2 = ', M2
      !write(*,*) ' = ', F
      ! print *, F
      !stop
      
      END SUBROUTINE FUNC
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

! Input arguments :
!      NDIM   :   Dimension of the ODE system 

! Values to be returned :
!      U      :   A starting solution vector
!      PAR    :   The corresponding equation-parameter values
!      T      :	  Not used here

      implicit none
      integer :: NDIM, j
      double precision :: U(NDIM), PAR(*), T    
      
      ! Initialize the equation parameters
      par(1) = 0.0500 ! this is epsilon
       
      call getsol(NDIM, U)

      par(2) = 1.0004
      par(3) = 0.4794

      U(NDIM-1) = 1.0004 ! this is mu (Note this is N)
      U(NDIM) = 0.4794 ! this is B (Note this is N+1)
       
      !do j = 1, NDIM
      !  print *, U(j)
      !end do
      !call sleep(10)

      END SUBROUTINE STPNT
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! The following subroutines are not used here,
! but they must be supplied as dummy routines

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      subroutine getsol(N, U)

      implicit none
      integer :: N, j
      double precision :: U(N)

      do j = 1, N
        U(j) = 0.0d0
      end do
      U(1) = 0.0204d0
      U(2) = 8.3295d-4

      end subroutine

