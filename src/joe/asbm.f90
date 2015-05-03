!=================================== ASBM ====================================80
!
! Computes the Reynolds stresses using the Algebraic Structure-based Model.
! On return:
!   ierr = 0 ok
!   ierr = 1 norm of error for AS implicit eq. too large
!   ierr = 2 trace_aa not within bounds
!   ierr = 3 trace of blocking vector out of bounds
!   ierr = 4 flattening parameter is NAN
!
!=============================================================================80

  subroutine asbm(s, w, wft, bl, as, ar, rey, dmn, cir, ktrmax, bltype, ierr)
    implicit none
    integer, parameter  :: dp = selected_real_kind(15, 307) !double precision

    ! INITIAL DECLARATIONS
    ! Routine's inputs and outputs
    real(dp), dimension(3,3), intent(in)    :: s   !(Rate of strain)*timescale
    real(dp), dimension(3,3), intent(in)    :: w   !(Mean rotation)*timescale
    real(dp), dimension(3,3), intent(in)    :: wft !(Frame rotation)*timescale
    real(dp), dimension(3,3), intent(in)    :: bl  !wall blocking tensor
    real(dp), dimension(3,3), intent(inout) :: as  !eddy axis tensor from s
    real(dp), dimension(3,3), intent(inout) :: ar  !eddy axis tensor from s, w

    real(dp), dimension(3,3), intent(out)   :: rey  !reynolds stresses
    real(dp), dimension(3,3), intent(out)   :: dmn  !dimensionality
    real(dp), dimension(3,3), intent(out)   :: cir  !circulicity

    integer, intent(in)                     :: ktrmax  !max iters for N-R
    integer, intent(in)                     :: bltype  !wall blocking type 
    integer, intent(inout)                  :: ierr    !error flag

    ! Constants
    real(dp), parameter                     :: a0 = 2.5_dp
    real(dp), parameter                     :: pi = 3.14159265359_dp

    real(dp), parameter                     :: zero = 0.0_dp
    real(dp), parameter                     :: small = 1.0e-06_dp
    real(dp), parameter                     :: fifth = 0.2_dp
    real(dp), parameter                     :: fourth = 0.25_dp
    real(dp), parameter                     :: third = 1.0_dp / 3.0_dp
    real(dp), parameter                     :: half = 0.5_dp
    real(dp), parameter                     :: twoth = 2.0_dp*third
    real(dp), parameter                     :: one = 1.0_dp
    real(dp), parameter                     :: two = 2.0_dp
    real(dp), parameter                     :: three = 3.0_dp
    real(dp), parameter                     :: six = 6.0_dp
    real(dp), parameter, dimension(3,3)     :: delta =                         &
         reshape((/1., 0., 0., 0., 1., 0., 0., 0., 1./),(/3,3/))
    real(dp), parameter, dimension(3,3,3)   :: eps =                           &
         reshape((/0., 0., 0., 0., 0., -1., 0., 1., 0.,                        &
                   0., 0., 1., 0., 0., 0., -1., 0., 0.,                        &
                   0., -1., 0., 1., 0., 0., 0., 0., 0./),(/3,3,3/))
    real(dp), parameter                     :: a_error = 1.0e-10_dp
    
    ! Variables
    integer                  :: i,j,k,l,m,n   !do loop indices
    integer                  :: ktr           !iteration index for N-R
    integer                  :: id, idx       !indeces for N-R vector & matrix
    integer, dimension(3,3)  :: index =                                        &
      reshape((/1, 2, 3, 2, 4, 5, 3, 5, 6/),(/3, 3/))

    logical                  :: strain        !true for strained flows
    logical                  :: rotation      !true for flows with mean rot
    logical                  :: converged     !used for N-R solver
    
    real(dp), dimension(3,3) :: sw            !s_ik*w_kj
    real(dp), dimension(3,3) :: ss            !s_ik*s_kj
    real(dp), dimension(3,3) :: ww            !w_ik*w_kj
    real(dp), dimension(3,3) :: wss           !w_ip*s_pq*s_qj
    real(dp), dimension(3,3) :: wws           !w_ip*w_pq*s_qj
    real(dp), dimension(3,3) :: wsww          !w_ip*s_pq*w_qt*w_tj 
    real(dp), dimension(3,3) :: swss          !s_ip*w_pq*s_qt*s_tj
    real(dp), dimension(3,3) :: wwss          !w_ip*w_pq*s_qt*s_tj
    real(dp), dimension(3,3) :: wssww         !w_ip*s_pq*s_qt*w_tr*s_rj
    real(dp)                 :: trace_ss      !s_ik*s_ki
    real(dp)                 :: trace_ww      !w_ik*w_ki
    real(dp)                 :: trace_sss     !s_ip*s_pq*s_qi
    real(dp)                 :: trace_wws     !w_ip*w_pq*s_qi
    real(dp)                 :: trace_wwss    !w_ip*w_pq*s_qt*s_ti
    real(dp), dimension(3,3) :: basis1        !1st element of integrity basis   
    real(dp), dimension(3,3) :: basis2        !2nd element of integrity basis
    real(dp), dimension(3,3) :: basis3        !3rd element of integrity basis
    real(dp), dimension(3,3) :: basis4        !4th element of integrity basis
    real(dp), dimension(3,3) :: basis5        !5th element of integrity basis
    real(dp), dimension(3,3) :: basis6        !6th element of integrity basis
    real(dp), dimension(3,3) :: basis7        !7th element of integrity basis
    real(dp), dimension(3,3) :: basis8        !8th element of integrity basis
    real(dp), dimension(3,3) :: basis9        !9th element of integrity basis
    real(dp), dimension(3,3) :: basis10       !10th element of integrity basis
    real(dp), dimension(3,3) :: test

    real(dp), dimension(3,3) :: a             !eddy axis tensor
    real(dp)                 :: trace_bl      !bl_ii
    real(dp)                 :: norm_s        !norm of s_ik
    real(dp)                 :: norm_w        !norm of w_ik

    ! variables for the strained-only eddy-axis tensor
    real(dp)                 :: i2, i3        !firs and second invariant of s
    real(dp)                 :: na,nb,nc,nd   !coefficients for N equation
    real(dp)                 :: ndelta        !discriminant for N equation
    real(dp)                 :: ndelta0       !variable for N equation
    real(dp)                 :: ndelta1       !variable for N equation
    real(dp)                 :: p1, p2        !variables for N equation
    real(dp)                 :: nbig          !N (production/dissipation)
    real(dp)                 :: as_alpha      !part of denom. for as coeff's
    real(dp)                 :: as_a          !coeff for as
    real(dp)                 :: as_b          !coeff for as

    ! variables for the rotated eddy-axis tensor
    real(dp)                 :: trace_aa      !a_ij*a_ji
    real(dp)                 :: r_ratio       !w_qr*s_rp*a_pq/s_kn*s_nm*a_mk
    real(dp)                 :: alpha         !2nd coeff for rot matrix H
    real(dp)                 :: beta          !3rd coeff for rot matrix H
    real(dp)                 :: alpha2        !alpha squared
    real(dp)                 :: term          !needed for beta

    ! variables for the structure  scalars
    real(dp)                 :: hat_w         !-a_ij*w_ik*w_kj
    real(dp)                 :: hat_s         !a_ij*s_ik*s_kj
    real(dp)                 :: eta_r         !mean rot over strain parameter
    real(dp)                 :: eta_f         !frame rot over strain parameter
    real(dp)                 :: eta_c1        !used to compute eta_r
    real(dp)                 :: eta_c2        !used to compute eta_f
    real(dp)                 :: oma           !one minus trace_aa
    real(dp)                 :: sqamth        !sqrt(trace_aa - third)

    real(dp)                 :: phi           !jettal scalar
    real(dp)                 :: bet           !correlation scalar
    real(dp)                 :: chi           !flattening scalar
    real(dp)                 :: gam           !helical scalar

    real(dp)                 :: phis          !jettal scalar for any a
    real(dp)                 :: bets          !correlation scalar for any a
    real(dp)                 :: chis          !flattening scalar for any a

    real(dp)                 :: phi1          !jettal scalar for shear
    real(dp)                 :: bet1          !correlation scalar for shear
    real(dp)                 :: chi1          !flattening scalar for shear

    real(dp)                 :: struc_weight  !smoothing parameter
    real(dp)                 :: xp_aa         !extrapolation along trace_aa

    real(dp)                 :: big_phi       !-half + three*half*phi
    real(dp)                 :: big_gam       !gam/sqrt(2)
    real(dp), dimension(10)  :: coeff_p       !phi coeff for integrity basis
    real(dp), dimension(10)  :: coeff_g       !gamma coeff for integrity basis
    real(dp), dimension(10)  :: coeff         !combined coeff_p and coeff_g
    real(dp), dimension(3,3) :: b             !anisotropic component of rij

    continue

    ierr = 0

    ! --------------------------------------------------------------------------
    ! COMPUTE TENSOR BASIS AND INVARIANTS
    ! --------------------------------------------------------------------------
    sw    = zero
    ss    = zero
    ww    = zero
    wss   = zero
    wws   = zero
    wsww  = zero
    swss  = zero
    wwss  = zero
    wssww = zero

    trace_ss   = zero
    trace_ww   = zero
    trace_sss  = zero
    trace_wws  = zero
    trace_wwss = zero

    do i = 1,3
      do j = 1,3
        do k = 1,3
          sw(i,j) = sw(i,j) + s(i,k)*w(k,j)
          ss(i,j) = ss(i,j) + s(i,k)*s(k,j)
          ww(i,j) = ww(i,j) + w(i,k)*w(k,j)
        end do
      end do
      trace_ss = trace_ss + ss(i,i)
      trace_ww = trace_ww + ww(i,i)
    end do

    do i = 1,3
      do j = 1,3
        do k = 1,3
          wss(i,j) = wss(i,j) + w(i,k)*ss(k,j)
          wws(i,j) = wws(i,j) + ww(i,k)*s(k,j)
          swss(i,j) = swss(i,j) + sw(i,k)*ss(k,j)
          wwss(i,j) = wwss(i,j) + ww(i,k)*ss(k,j)
          do l = 1,3
            wsww(i,j) = wsww(i,j) + w(i,k)*s(k,l)*ww(l,j)
            wssww(i,j) = wssww(i,j) + w(i,k)*ss(k,l)*ww(l,j)
          end do
        end do
        trace_sss = trace_sss + ss(i,j)*s(j,i)
      end do
      trace_wws = trace_wws + wws(i,i)
      trace_wwss = trace_wwss + wwss(i,i)
    end do

    ! elements of integrity basis
    basis1  = s
    basis2  = sw + transpose(sw)
    basis3  = ss - third*trace_ss*delta
    basis4  = ww - third*trace_ww*delta
    basis5  = wss + transpose(wss)
    basis6  = wws + transpose(wws) - twoth*trace_wws*delta
    basis7  = wsww + transpose(wsww)
    basis8  = swss + transpose(swss)
    basis9  = wwss + transpose(wwss) - twoth*trace_wwss*delta
    basis10 = wssww + transpose(wssww)

    if (trace_ss > zero) then
      strain = .true.
      norm_s = sqrt(trace_ss)
    else
      strain  = .false.
      norm_s = zero
    end if

    if (trace_ww < zero) then
      rotation = .true.
      norm_w = sqrt(-trace_ww)
    else
      rotation = .false.
      norm_w  = zero
    end if

    ! --------------------------------------------------------------------------
    ! COEFFICIETNS FOR AS
    ! --------------------------------------------------------------------------
    as_a = zero
    as_b = zero

    ! check for strain
    purestrain: if (strain) then
      ! invariants
      i2 = -half*trace_ss
      i3 = third*trace_sss

      ! compute N 
      na = three
      nb = -three*a0
      nc = 12.0_dp*i2
      nd = -(4.0_dp*a0*i2 + 24.0_dp*i3)

      ndelta = 18.0_dp*na*nb*nc*nd - 4.0_dp*nb**three*nd + nb**two*nc**two    &
               - 4.0_dp*na*nc**three - 27.0_dp*na**two*nd**two
      ndelta0 = nb**two - three*na*nc
      ndelta1 = two*nb**three - 9.0_dp*na*nb*nc + 27.0_dp*na**two*nd

      p1 = half*ndelta1
      p2 = half*three*na*sqrt(three*abs(ndelta))

      nbig = nb
      if (ndelta <= zero) then
        nbig = nbig + (p1 + p2)**(third) + ndelta0*(p1 + p2)**(-third)
      else
        nbig = nbig + two*(p1*p1 + p2*p2)**(one/six)*                    &
               cos(third*acos(p1/sqrt(p1*p1 + p2*p2)) + twoth*pi)
      end if
      nbig = -third*nbig/na

      ! compute as
      as_alpha = three*nbig**two + 4.0_dp*i2

      as_a = two*nbig/as_alpha
      as_b = 4.0_dp/as_alpha

    end if purestrain

    ! --------------------------------------------------------------------------
    ! COEFFICIENTS FOR AR
    ! --------------------------------------------------------------------------
    alpha = zero
    beta  = zero

    ! check for rotation
    purerotation: if (rotation) then

      trace_aa = third + as_a**two*trace_ss +                                  &
                 one/six*as_b**two*trace_ss**two + two*as_a*as_b*trace_sss

      r_ratio = -trace_ww/trace_ss*(1.5_dp*(trace_aa - third))**half

      ! compute alpha and beta
      if (r_ratio < one) then
        alpha2 = one - sqrt(one - r_ratio)     !hyperbolic mean flow
      else
        alpha2 = one + sqrt(one - one/r_ratio) !elliptic mean flow
      end if

      if (alpha2 < zero) alpha2 = zero
      alpha = sqrt(alpha2)
        
      term = 4.0_dp - 2.0_dp*alpha2
      if (term < zero) term = zero
      beta = 2.0_dp - sqrt(term)

    end if purerotation

    ! coefficients for eddy-axis tensor
    coeff_p(1)  = as_a*(one - half*(beta**two - two*beta))
    coeff_p(2)  = -one/sqrt(-trace_ww)*as_a*alpha
    coeff_p(3)  = as_b*(one - half*(beta**two - two*beta))
    coeff_p(4)  = - trace_ss/trace_ww*as_b*beta*(beta - two)                   &
                  + trace_wws/(trace_ww**two)*as_a*beta**two                   &
                  + trace_wwss/(trace_ww**two)*as_b*beta**two
    coeff_p(5)  = one/sqrt(-trace_ww)*as_b*alpha
    coeff_p(6)  = one/trace_ww*as_a*(beta**two - three*beta)
    coeff_p(7)  = -one/(sqrt(-trace_ww)*trace_ww)*as_a*alpha*beta
    coeff_p(8)  = zero
    coeff_p(9)  = one/trace_ww*as_b*(beta**two - three*beta)
    coeff_p(10) = -one/(sqrt(-trace_ww)*trace_ww)*as_b*alpha*beta

    ! --------------------------------------------------------------------------
    ! COMPUTE STRUTURE SCALARS
    ! --------------------------------------------------------------------------

    ! bounds check for trace_aa
    if (trace_aa > one) then
      write(*,*) 'trace_aa = ', trace_aa, ' greater than one'
      trace_aa = one
      ierr = 2
    end if
    
    if (trace_aa < third) then
      write(*,*) 'trace_aa = ', trace_aa, ' less than third'
      trace_aa = third
      ierr = 2
    end if

    ! compute eta_r
    hat_w = third*trace_ww                                                     &
          + coeff_p(1)*trace_wws                                               &
          + coeff_p(3)*(trace_wwss - third*trace_ss*trace_ww)                  &
          + coeff_p(4)*one/six*trace_ww**two                                   &
          + coeff_p(6)*third*trace_ww*trace_wws                                &
          + coeff_p(9)*third*trace_ww*trace_wwss

    hat_s = third*trace_ss                                                     &
          + coeff_p(1)*trace_sss                                               &
          + coeff_p(3)*one/six*trace_ss**two                                   &
          + coeff_p(4)*(trace_wwss - third*trace_ss*trace_ww)                  &
          + coeff_p(6)*(third*trace_ss*trace_wws + twoth*trace_ww*trace_sss)   &
          + coeff_p(9)*(third*trace_ss*trace_wwss + twoth*trace_wws*trace_sss)

    eta_c1 = -hat_w/(hat_s + 1.0e-06_dp)

    if ((eta_c1 < zero) .or. (eta_c1 /= eta_c1)) eta_c1 = zero

    eta_r = sqrt(eta_c1)
    eta_f = zero

    ! compute phis, chis, bets
    if (hat_s < zero) then
      ! without strain
      phis = zero
      chis = zero
      bets = one

      if (hat_w > zero) then
        ! with rotation
        phis = third
        chis = zero
        bets = zero
      end if
    else
      ! with strain
      oma = one - trace_aa
      sqamth = sqrt(trace_aa - third)

      if (eta_r <= one) then
        call int_er_lt_one(eta_r,eta_f,oma,sqamth,phis,bets,chis,phi1,bet1,chi1)
      else 
        call int_er_gt_one(eta_r,eta_f,oma,sqamth,phis,bets,chis,phi1,bet1,chi1)
      end if
      struc_weight = exp(-1000.0_dp*abs(eta_r - one)**two)
      !phis = phis*(one - struc_weight) + phi1*(struc_weight)
      bets = bets*(one - struc_weight) + bet1*(struc_weight)
      chis = chis*(one - struc_weight) + chi1*(struc_weight)
    end if
   
    ! compute phi, chi, bet
    xp_aa = 1.5_dp*(trace_aa - third)
    xp_aa = 0.35_dp*xp_aa**2.5_dp + 0.65_dp*xp_aa**half

    phi = phis*xp_aa
    chi = chis*xp_aa
    if (eta_r <= one) then
      bet = bets
    else
      bet = one - max(one - 0.9_dp*(eta_r - one)**0.31_dp, zero)*              &
                  (1.5_dp*(trace_aa - third))**10.0_dp
    end if
    struc_weight = exp(-1000.0_dp*abs(eta_r - one)**two)  !ale
    bet = bet*(one - struc_weight) + bet1*(struc_weight) !ale

    ! compute helical scalar and vector
    gam = two*phi*(one - phi)/(one + chi)
    if (gam < zero) gam = zero
    gam = bet*sqrt(gam)

    ! --------------------------------------------------------------------------
    ! COMPUTE STRUCTURE TENSORS
    ! --------------------------------------------------------------------------
    big_phi = -half + half*three*phi
    big_gam = gam/sqrt(two)

    ! coefficients for the Reynolds stresses
    coeff_g(1)  = as_a*(alpha - alpha*beta)
    coeff_g(2)  = -one/sqrt(-trace_ww)*as_a*(one - half*beta)
    coeff_g(3)  = as_b*(alpha - alpha*beta)
    coeff_g(4)  = trace_ss/trace_ww*two*as_b*(alpha - alpha*beta)              &
                + trace_wws/(trace_ww**two)*two*as_a*alpha*beta                &
                + trace_wwss/(trace_ww**two)*two*as_b*alpha*beta
    coeff_g(5)  = one/sqrt(-trace_ww)*as_b*(one - half*beta)
    coeff_g(6)  = -one/trace_ww*as_a*three*(alpha - twoth*alpha*beta)
    coeff_g(7)  = one/(sqrt(-trace_ww)*trace_ww)*as_a*(beta**two - three*beta)
    coeff_g(8)  = zero
    coeff_g(9)  = -one/trace_ww*three*as_b*(alpha - twoth*alpha*beta)
    coeff_g(10) = one/(sqrt(-trace_ww)*trace_ww)*as_b*(beta**two - three*beta)

    do i = 1,10
      coeff(i) = big_phi*coeff_p(i) + big_gam*coeff_g(i)
    end do

    b = coeff(1)*basis1 + coeff(2)*basis2 + coeff(3)*basis3 +                  &
        coeff(4)*basis4 + coeff(5)*basis5 + coeff(6)*basis6 +                  &
        coeff(7)*basis7 + coeff(8)*basis8 + coeff(9)*basis9 +                  &
        coeff(10)*basis10

    rey = third*delta + b

    ! --------------------------------------------------------------------------
    ! NEAR-WALL CORRECTION
    ! --------------------------------------------------------------------------
    if (bltype == 1) then
      ! blockage correction to Reynolds stress tensor
      ! Note: The other tensors need to be updated
      call blocking(rey,bl,delta,ierr)
    endif

    ! Output some useful data
    dmn(1,1) = coeff(1)
    dmn(1,2) = coeff(2)
    dmn(1,3) = coeff(3)

    dmn(2,1) = coeff(4)
    dmn(2,2) = coeff(5)
    dmn(2,3) = coeff(6)

    dmn(3,1) = coeff(7)
    dmn(3,2) = coeff(9)
    dmn(3,3) = coeff(10)

    cir(1,1) = trace_aa
    cir(1,2) = r_ratio
    cir(1,3) = eta_r

  end subroutine asbm

!=============================== INT_ER_LT_ONE ===============================80
!
! Computes the structure parameters for an arbitrary a_ij by interpolating
! between the plane strain and pure shear states.
!
!=============================================================================80

  subroutine int_er_lt_one(eta_r,eta_f,oma,sqamth,phis,bets,chis,phi1,bet1,chi1)
    implicit none
    integer, parameter  :: dp = selected_real_kind(15, 307) !double precision

    ! INITIAL DECLARATIONS
    ! routine's inputs and outputs
    real(dp), intent(in)    :: eta_r
    real(dp), intent(in)    :: eta_f
    real(dp), intent(in)    :: oma
    real(dp), intent(in)    :: sqamth
    real(dp), intent(inout) :: phis
    real(dp), intent(inout) :: bets
    real(dp), intent(inout) :: chis
    real(dp), intent(inout) :: phi1
    real(dp), intent(inout) :: bet1
    real(dp), intent(inout) :: chi1

    ! constants
    real(dp), parameter     :: zero = 0.0_dp
    real(dp), parameter     :: one = 1.0_dp
    real(dp), parameter     :: two = 2.0_dp
    
    ! variables
    real(dp)                :: eta_f1, eta_f0
    real(dp)                :: phi0
    real(dp)                :: bet0
    real(dp)                :: chi0
    real(dp)                :: param 
    
    continue

    param = sqrt(3.0_dp)/4.0_dp
    if (eta_f < (param*(eta_r - one))) then 
      eta_f1 = eta_f - param*(eta_r - one)
      eta_f0 = eta_f - param*(eta_r - one) - param
    else if (eta_f < zero) then
      eta_f1 = zero
      eta_f0 = eta_f/(one - eta_r)
    else if (eta_f < eta_r) then
      eta_f1 = eta_f/eta_r
      eta_f0 = zero
    else if (eta_f < (eta_r + param*(one - eta_r))) then
      eta_f1 = one
      eta_f0 = one - (one - eta_f)/(one - eta_r)
    else
      eta_f1 = one + eta_f - eta_r - param*(one - eta_r)
      eta_f0 = param + eta_f - eta_r - param*(one - eta_r)
    end if

    ! compute parameters for shear and plane strain states
    call pure_shear(eta_f1,oma,sqamth,phi1,bet1,chi1)
    call plane_strain(eta_f0,oma,sqamth,phi0,bet0,chi0)
    
    ! interpolate along eta_r direction
    phis = phi0 + (phi1 - phi0)*                                               &
                  (0.82_dp*eta_r**two)/(one - (one-0.82_dp)*eta_r**two)
    bets = bet0 + (bet1 - bet0)*eta_r**two
    chis = chi0 + (chi1 - chi0)*eta_r**two
  
  end subroutine int_er_lt_one

!=============================== INT_ER_GT_ONE ===============================80
!
! Computes the structure parameters for an arbitrary a_ij by extrapolating
! for higher values of eta_r than that of pure shear.
!
!=============================================================================80

  subroutine int_er_gt_one(eta_r,eta_f,oma,sqamth,phis,bets,chis,phi1,bet1,chi1)
    implicit none 
    integer, parameter  :: dp = selected_real_kind(15, 307) !double precision

    ! INITIAL DECLARATIONS
    ! routine's inputs and outputs
    real(dp), intent(in)    :: eta_r
    real(dp), intent(in)    :: eta_f
    real(dp), intent(in)    :: oma
    real(dp), intent(in)    :: sqamth
    real(dp), intent(inout) :: phis
    real(dp), intent(inout) :: bets
    real(dp), intent(inout) :: chis
    real(dp), intent(inout) :: phi1
    real(dp), intent(inout) :: bet1
    real(dp), intent(inout) :: chi1

    ! constants
    real(dp), parameter     :: third = 1.0_dp / 3.0_dp
    real(dp), parameter     :: half = 1.0_dp / 2.0_dp
    real(dp), parameter     :: one = 1.0_dp

    ! variables
    real(dp)                :: aux
    real(dp)                :: aux_shift

    continue

    aux = one/(one + half*(eta_r - one)/(oma)**2.5_dp)
    aux_shift = one/(one + half*(eta_r - one)/(oma+third)**2.5_dp)

    ! compute parameters for shear state
    call pure_shear(eta_f,oma,sqamth,phi1,bet1,chi1)
    
    ! extrapolate along eta_r direction
    phis = third + (phi1 - third)*aux_shift
    bets = bet1*aux
    chis = chi1*aux
  
  end subroutine int_er_gt_one

!================================ PURE_SHEAR =================================80
!
! Computes reference structure parameters along the shear line.
!
!=============================================================================80

  subroutine pure_shear(eta_f,oma,sqamth,phi1,bet1,chi1)
    implicit none 
    integer, parameter  :: dp = selected_real_kind(15, 307) !double precision

    ! INITIAL DECLARATIONS
    ! routine's inputs and outputs
    real(dp), intent(in)    :: eta_f
    real(dp), intent(in)    :: oma
    real(dp), intent(in)    :: sqamth
    real(dp), intent(inout) :: phi1
    real(dp), intent(inout) :: bet1
    real(dp), intent(inout) :: chi1

    ! constants
    real(dp), parameter     :: zero = 0.0_dp
    real(dp), parameter     :: fifth = 0.2_dp
    real(dp), parameter     :: one = 1.0_dp

    continue

    if (eta_f < zero) then
      phi1 = (eta_f - one)/(3.0_dp*eta_f - one)
      bet1 = one/(one - eta_f*(one + sqamth)/oma)
      chi1 = fifth*bet1
    else if (eta_f < one) then
      phi1 = one - eta_f
      chi1 = zero
      !chi1 = fifth + (one - fifth)*                                           &
      !       (one - (one - eta_f)**2.0_dp/(one + 3.0_dp*eta_f/oma))
      bet1 = one
    else
      phi1 = (eta_f - one)/(3.0_dp*eta_f - one)
      bet1 = one/(one + 0.8_dp*(eta_f - one)*(eta_f*sqamth)/oma)
      chi1 = one - (one - bet1)*(eta_f - one)/(oma + eta_f - one)
    end if

  end subroutine pure_shear

!============================== PLANE_STRAIN =================================80
!
! Computes reference structure parameters along the shear line.
!
!=============================================================================80

  subroutine plane_strain(eta_f,oma,sqamth,phi0,bet0,chi0)
    implicit none 
    integer, parameter  :: dp = selected_real_kind(15, 307) !double precision

    ! INITIAL DECLARATIONS
    ! routine's inputs and outputs
    real(dp), intent(in)    :: eta_f
    real(dp), intent(in)    :: oma
    real(dp), intent(in)    :: sqamth
    real(dp), intent(inout) :: phi0
    real(dp), intent(inout) :: bet0
    real(dp), intent(inout) :: chi0

    ! constants
    real(dp), parameter     :: zero = 0.0_dp
    real(dp), parameter     :: one = 1.0_dp
    real(dp), parameter     :: two = 2.0_dp
    real(dp), parameter     :: three = 3.0_dp

    ! variables
    real(dp)                :: var1
    real(dp)                :: var2

    continue
    
    var1 = two*eta_f
    var2 = var1*var1

    if ( var2 > 0.75_dp) then
      bet0 = one/(one + (var1 - sqrt(0.75_dp))*(var1*sqamth)/oma)
      phi0 = (one - bet0)/three
      chi0 = -bet0
    else
      phi0 = 0.145_dp*(var2/0.75_dp - (var2/0.75_dp)**9.0_dp)
      chi0 = -(0.342_dp*(var2/0.75_dp) +                                       &
             (one - 0.342_dp)*(var2/0.75_dp)**6.0_dp)
      bet0 = one
    end if
    chi0 = zero

  end subroutine plane_strain

!================================ BLOCKING ===================================80
!
! Modifies the eddy axis tensor to account for wall effects.
! Input:
!        a(i,j)       unblocked (homogeneous) tensor
!        bl(i,j)      blocking tensor
! Output
!        a(i,j)       blocked tensor
!
!=============================================================================80

  subroutine blocking(a,bl,delta,ierr)
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
    real(dp), dimension(3,3) :: h          !partial projection operator
    real(dp)                 :: d2         !normalizing factor
    real(dp)                 :: dinv       !one over sqrt(d2)
    real(dp)                 :: suma

    continue

    trace_bl = bl(1,1) + bl(2,2) + bl(3,3)

    ! return if no blockage
    if (trace_bl < small) return

    ! check for bad data
    if (trace_bl > one) then
      ierr = 3
      return
    end if

    ah = a

    ! compute normalizing factor
    suma = zero
    do i = 1,3
      do j = 1,3
        suma = suma + ah(i,j)*bl(i,j)
      end do
    end do
    d2 = one - (2.0_dp - trace_bl)*suma
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

  end subroutine blocking
