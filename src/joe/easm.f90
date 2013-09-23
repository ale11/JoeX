  subroutine easm(st, wt, wft, bl, as, ar, rey, dmn, cir, ktrmax, bl_rij, ierr)
    implicit none
    integer, parameter  :: dp = selected_real_kind(15, 307) !double precision
    real(dp), dimension(3,3), intent(in)    :: st  !(Rate of strain)*timescale
    real(dp), dimension(3,3), intent(in)    :: wt  !(Mean rotation)*timescale
    real(dp), dimension(3,3), intent(in)    :: wft !(Frame rotation)*timescale
    real(dp), dimension(3,3), intent(in)    :: bl  !wall blocking tensor
    real(dp), dimension(3,3), intent(inout) :: as  !eddy axis tensor from st
    real(dp), dimension(3,3), intent(inout) :: ar  !eddy axis tensor from st, wt

    real(dp), dimension(3,3), intent(out)   :: rey  !reynolds stresses
    real(dp), dimension(3,3), intent(in)    :: dmn  !dimensionality
    real(dp), dimension(3,3), intent(in)    :: cir  !circulicity

    real(dp), dimension(3,3)     :: delta =                          &
         reshape((/1., 0., 0., 0., 1., 0., 0., 0., 1./),(/3,3/))

    integer, intent(in)                     :: ktrmax  !max iters for N-R
    integer, intent(in)                     :: bl_rij  !wall blocking type
    integer, intent(inout)                  :: ierr    !error flag

!   W-J version of LRR:
    real(dp), parameter :: c10 = 3.4_dp
    real(dp), parameter :: c11 = 1.8_dp
    real(dp), parameter :: c2  = 0.36_dp
    real(dp), parameter :: c3  = 1.25_dp
    real(dp), parameter :: c4  = 0.4_dp
    real(dp), parameter :: cs0 = 0._dp

    real(dp), parameter :: gamma0 = 0.5_dp*c11
    real(dp), parameter :: gamma1 = 0.5_dp*c10 + 0.8864_dp
    real(dp), parameter :: a1 = 0.5_dp*(4.0_dp/3.0_dp - c2)
    real(dp), parameter :: a2 = 0.5_dp*(2.0_dp - c4)
    real(dp), parameter :: a3 = 0.5_dp*(2.0_dp - c3)

    real(dp) :: S11,S12,S13,S21,S22,S23,S31,S32,S33
    real(dp) :: W11,W12,W13,W21,W22,W23,W31,W32,W33
    real(dp) :: eta1, eta2, cmu, cmuj, a4, SijSij, WijWij
    real(dp) :: S1kWk1, W1kSk1, S1kSk1
    real(dp) :: S2kWk2, W2kSk2, S2kSk2
    real(dp) :: S3kWk3, W3kSk3, S3kSk3
    real(dp) :: S1kWk2, W1kSk2, S1kSk2
    real(dp) :: S1kWk3, W1kSk3, S1kSk3
    real(dp) :: S2kWk3, W2kSk3, S2kSk3

    S11=St(1,1)
   	S12=St(1,2)
    S13=St(1,3)
    S21=St(2,1)
    S22=St(2,2)
    S23=St(2,3)
    S31=St(3,1)
    S32=St(3,2)
    S33=St(3,3)
    W11=Wt(1,1)
    W12=Wt(1,2)
    W13=Wt(1,3)
    W21=Wt(2,1)
    W22=Wt(2,2)
    W23=Wt(2,3)
    W31=Wt(3,1)
    W32=Wt(3,2)
    W33=Wt(3,3)
    SijSij=S11*S11 + S22*S22 + S33*S33 + 2.*S12*S12 + 2.*S13*S13 + 2.*S23*S23
    WijWij=2.*W12*W12 + 2.*W13*W13 + 2.*W23*W23

    S1kWk1=S11*W11-S12*W12-S13*W13
    W1kSk1=W11*S11+W12*S12+W13*S13
    S1kSk1=S11*S11+S12*S12+S13*S13

    S1kWk2=S11*W12+S12*W22-S13*W13
    W1kSk2=W11*S12+W12*S22+W13*S23
    S1kSk2=S11*S12+S12*S22+S13*S23

    S2kWk2=S12*W12+S22*W22-S23*W23
    W2kSk2=-W12*S12+W22*S22+W23*S23
    S2kSk2=S12*S12+S22*S22+S23*S23

    S3kWk3= S31*W13+S32*W23+S33*W33
    W3kSk3= W31*S13+W32*S23+W33*W33
    S3kSk3= S31*S13+S32*S23+S33*S33

    S1kWk3= S11*W13+S12*W23+S13*W33
    W1kSk3= W11*S13+W12*S23+W13*S33
    S1kSk3= S11*S13+S12*S23+S13*S33

    S2kWk3= S21*W13+S22*W23+S23*W33
    W2kSk3= W21*S13+W22*S23+W23*S33
    S2kSk3= S21*S13+S22*S23+S23*S33

    eta1 = min(SijSij, 1200.0_dp)
    eta2 = min(WijWij, 1200.0_dp)

    call cubic(c10,c11,c2,c3,c4,a1,a2,a3,gamma0,gamma1,eta1,eta2,cmuj,2,cs0)

    cmu = max(-cmuj, 0.0005_dp)

    a4 = 1.0_dp/(gamma1 + 2.0_dp*gamma0*cmu*eta1)

    rey(1,1) = 1.0_dp/3.0_dp - cmu*S11 - a2*a4*cmu*(S1kWk1 - W1kSk1) + 2.0_dp*a3*a4*cmu*(S1kSk1 - SijSij/3.0_dp)
    rey(2,2) = 1.0_dp/3.0_dp - cmu*S22 - a2*a4*cmu*(S2kWk2 - W2kSk2) + 2.0_dp*a3*a4*cmu*(S2kSk2 - SijSij/3.0_dp)
    rey(3,3) = 1.0_dp/3.0_dp - cmu*S33 - a2*a4*cmu*(S3kWk3 - W3kSk3) + 2.0_dp*a3*a4*cmu*(S3kSk3 - SijSij/3.0_dp)
    rey(1,2) = -cmu*S12 - a2*a4*cmu*(S1kWk2 - W1kSk2) + 2.0_dp*a3*a4*cmu*S1kSk2
	rey(2,1) = rey(1,2)
    rey(2,3) = -cmu*S23 - a2*a4*cmu*(S2kWk3 - W2kSk3) + 2.0_dp*a3*a4*cmu*S2kSk3
	rey(3,2) = rey(2,3)
    rey(1,3) = -cmu*S13 - a2*a4*cmu*(S1kWk3 - W1kSk3) + 2.0_dp*a3*a4*cmu*S1kSk3
	rey(3,1) = rey(1,3)

    call blockingr(rey,bl,delta,ierr)

    ierr=0
  end
! -----------------------------------------------------------
!
  subroutine cubic(c10,c11,c2,c3,c4,a1,a2,a3,g0,g1,eta1,eta2,cmuj,igo,cs0)
    implicit none 
    integer, parameter  :: dp = selected_real_kind(15, 307) !double precision

    real (dp) :: c10,c11,c2,c3,c4,a1,a2,a3,g0,g1,eta1,eta2,cmuj,cs0
    integer igo
    real (dp) :: pi,tpi3,eta11,ap,ar,aq,aa,ab,ad,raat,rbbt,raa,rbb,cmub,cmuc,coss,theta
     
    pi   = acos(-1.0_dp)
    tpi3 = 2.0_dp*pi/3.0_dp

!   *******-------------------------------------------*******
!   Solve cubic eqn now (solve for alpha/tau = -cmu):
!   (alpha/tau)**3 + p(alpha/tau)**2 + q(alpha/tau) + r = 0
!   *******-------------------------------------------*******
!   must handle special case when eta -> 0:
    if(eta1.le.1.e-6_dp) then
      cmuj=-g1*a1/(g1*g1 + 2.0_dp*eta2*a2*a2)
    else
!     set up parameters
!     ap is p, coefficient of (alpha1/tau)**2
!     aq is q, coefficient of (alpha1/tau)
!     ar is r
      eta11 = eta1*g0
      ap = -g1/eta11
      ar = g1*a1/(2.0_dp*eta11)**2
      aq = g1*g1 - 2.0_dp*eta11*a1 - 0.666667_dp*eta1*a3*a3 + 2.0_dp*eta2*a2*a2
      aq = aq/(4.0_dp*eta11*eta11)

      aa = (aq - ap*ap/3.0_dp)
      ab = (2.0_dp*ap*ap*ap - 9.0_dp*ap*aq + 27.0_dp*ar)/27.0_dp
      ad = (ab*ab/4.0_dp) + (aa*aa*aa)/27.0_dp
!     if ad > 0,  get one real and 2 imaginary roots:
      if(ad.gt.0.0_dp) then
        raat = -0.5_dp*ab + sqrt(ad)
        rbbt = -0.5_dp*ab - sqrt(ad)
        raa  = (abs(raat))**(1.0_dp/3.0_dp)
        raa  = sign(raa,raat)
        rbb  = (abs(rbbt))**(1.0_dp/3.0_dp)
        rbb  = sign(rbb,rbbt)
        cmuj = -ap/3.0_dp + raa + rbb
!       choose min of real root AND real part of imaj roots:
        cmub = -ap/3.0_dp - 0.5_dp*raa - 0.5_dp*rbb
        cmuj = min(cmuj,cmub)
      else
!     if ad < 0, get 3 real roots:
        coss  = -0.5_dp*ab/sqrt(-aa*aa*aa/27.0_dp)
        theta = acos(coss)
        cmuj  = -ap/3.0_dp + 2.0_dp*sqrt(-aa/3.0_dp)*cos(theta/3.0_dp)
        cmub  = -ap/3.0_dp + 2.0_dp*sqrt(-aa/3.0_dp)*cos(tpi3 + theta/3.0_dp)
        cmuc  = -ap/3.0_dp + 2.0_dp*sqrt(-aa/3.0_dp)*cos(2.0_dp*tpi3 + theta/3.0_dp)
!   choose minimum root:
        cmuj  = min(cmuj,cmub)
        cmuj  = min(cmuj,cmuc)
      end if
    end if
  end


 subroutine blockingr(a,bl,delta,ierr)
    implicit none 
    integer, parameter  :: dp = selected_real_kind(15, 307) !double precision
    
    ! INITIAL DECLARATIONS
    ! routine's inputs and outputs
    real(dp), dimension(3,3), intent(inout) :: a        !eddy axis tensor
    real(dp), dimension(3,3), intent(in)    :: bl       !blocking tensor
    real(dp), dimension(3,3), intent(in)    :: delta    !kronecker delta

    integer, intent(inout)                  :: ierr  !error flag
    
    ! constants
    real(dp)                 :: small = 1.0e-14_dp
    real(dp)                 :: zero = 0.0_dp
    real(dp)                 :: one = 1.0_dp
            
    ! variables
    integer                  :: i,j,k,l    !do loop indeces
    real(dp), dimension(3,3) :: ah         !homogeneous eddy axis tensor
    real(dp)                 :: trace_bl   !bl_ii
    real(dp)                 :: trace_ahbl !ah_in*bl_ni
    real(dp), dimension(3,3) :: h          !partial projection operator
    real(dp)                 :: d2         !normalizing factor
    real(dp)                 :: dinv       !one over sqrt(d2)
    real(dp)                 :: sum

    continue

    trace_bl = bl(1,1) + bl(2,2) + bl(3,3)

    ! return if no blockage
    if (trace_bl < small) return 
       
    ! check for bad data
    if (trace_bl > one) then
      ierr = 7
      return 
    end if

    ah = a

    ! compute normalizing factor
    sum = zero
    do i = 1,3
      do j = 1,3
        sum = sum + ah(i,j)*bl(i,j)
      end do
    end do
    d2 = one - (2.0_dp - trace_bl)*sum
    if (d2 < zero) return
    dinv = one/sqrt(d2);

    ! compute partial projection operator
    do i = 1,3
      do j = 1,3
        h(i,j) = dinv*(delta(i,j) - bl(i,j))
      end do
    end do

    ! apply blockage correction
    a = zero
    do i = 1,3
      do j = i,3
        do k = 1,3
          do l = 1,3
            a(i,j) = a(i,j) + h(i,k)*h(j,l)*ah(k,l)
          end do
        end do
        a(j,i) = a(i,j)
      end do
    end do

  end subroutine blockingr

