!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   BZX: An interface code for reading boozmn files produced 
!                         by BoozXform with VMEC (in STELOPT)
!                             
!                         by Motoki NAKATA, July 2016 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module BZX_header 

  implicit none

  integer, parameter :: DP = selected_real_kind(14)

! --- constants
  real(kind=DP),    parameter :: pi  = 3.141592653589793_DP, &
                                 twopi = pi * 2._DP

! --- ################## input file name and resoution settings ##################

  character(len=*), parameter :: fname_tag = "vmec_inward" ! file name index of VMEC data

!!!!!!!!!! --- To make field-aligned coord. for GKV: one must set nzeta=0
!!!  integer, parameter :: Ntheta_gkv = 3 ! N_tht value in GKV
!!!  integer, parameter :: nrho  = 501  ! radial grid number in [0 <= rho <= 1]
!!!  integer, parameter :: ntht  = 384  ! poloidal grid number in [-N_theta*pi < theta < N_theta*pi]
!!!  integer, parameter :: nzeta = 0  ! toroidal grid number in [0 <= zeta <= 2*pi] 
!!!
!!!!!!  real(kind=DP), parameter :: alpha_fix  = 0._DP ! field-line label: alpha = zeta - q*theta
!!!  real(kind=DP), parameter :: alpha_fix  = pi/10._DP ! field-line label


!!!!!!! --- To make full 3d coordinates 
  integer, parameter :: Ntheta_gkv = 1 ! N_tht value in GKV
  integer, parameter :: nrho  = 11  ! radial grid number in [0 <= rho <= 1]
  integer, parameter :: ntht  = 64  ! poloidal grid number in [-N_theta*pi < theta < N_theta*pi]
  integer, parameter :: nzeta = 64  ! toroidal grid number in [0 <= zeta <= 2*pi] 
  real(kind=DP), parameter :: alpha_fix  = 0._DP ! field-line label: alpha = zeta - q*theta NOT USED in 3d case


! --- ############################################################################


  ! --- file names 
  character(len=*), parameter :: fname_boozmn = "boozmn."//fname_tag,       &
                                 fname_wout   = "wout_"//fname_tag//".txt", &
                                   wout_txt   = fname_tag

  character(len=*), parameter :: file_log  = "./log_BZX.dat"
  character(len=*), parameter :: file_prf  = "./geom/prof.dat"
  character(len=*), parameter :: file_out1 = "./geom/shape_rtz.dat"
  character(len=*), parameter :: file_out2 = "./geom/shape_trz.dat"
  character(len=*), parameter :: file_out3 = "./geom/shape_tzr.dat"
  character(len=*), parameter :: file_out4 = "./geom/metric_boozer.bin.dat"

  ! --- for reading from boozxform 
  real(kind=DP), dimension(:), allocatable :: iota_b_nu, pres_b_nu,           &
                                              phip_b_nu, phi_b_nu, beta_b_nu, & 
                                              buco_b_nu, bvco_b_nu, qq_nu

  real(kind=DP), dimension(:), allocatable :: iota_b, pres_b, phip_b, phi_b, beta_b,       & 
                                              buco_b, bvco_b, rho, rho_nu, rho2, dphidrho, &
                                              qq, dqdrho, shat, epst, cug, cui, dummy1

  real(kind=DP), dimension(:,:), allocatable :: bmnc_b, rmnc_b, zmns_b, pmns_b, gmnc_b
  real(kind=DP), dimension(:,:), allocatable :: bmns_b, rmns_b, zmnc_b, pmnc_b, gmns_b

  real(kind=DP) :: aspect_b, rmax_b, rmin_b, betaxis_b

  integer :: mnboz_b, mboz_b, nboz_b, nfp_b, ns_b
  integer, dimension(:), allocatable :: ixm_b, ixn_b, jlist
  logical :: lasym_b = .false.

  ! --- for metric calc.
  real(kind=DP), dimension(:,:), allocatable :: bbozc_nu, rbozc_nu, zbozs_nu, pbozs_nu, gbozc_nu 
  real(kind=DP), dimension(:,:), allocatable :: bbozs_nu, rbozs_nu, zbozc_nu, pbozc_nu, gbozs_nu
  real(kind=DP), dimension(:,:), allocatable :: bbozc, rbozc, zbozs, pbozs, gbozc
  real(kind=DP), dimension(:,:), allocatable :: bbozs, rbozs, zbozc, pbozc, gbozs
  real(kind=DP), dimension(:,:), allocatable :: dbbozc, drbozc, dzbozs, dpbozs, dummy2
  real(kind=DP), dimension(:,:), allocatable :: dbbozs, drbozs, dzbozc, dpbozc

  real(kind=DP), dimension(:,:,:), allocatable :: dbb_drho, dbb_dtht, dbb_dzeta
  real(kind=DP), dimension(:,:,:), allocatable :: drr_drho, drr_dtht, drr_dzeta
  real(kind=DP), dimension(:,:,:), allocatable :: dzz_drho, dzz_dtht, dzz_dzeta
  real(kind=DP), dimension(:,:,:), allocatable :: dph_drho, dph_dtht, dph_dzeta

  ! --- for wout read
  real(kind=DP)  :: B0_p, Aminor_p, Rmajor_p, volume_p

  ! --- I/O unit numbers
  integer, parameter :: ibzmn = 100  ! input for boozmn produced by Booz_Xform
  integer, parameter :: iwout = 110  ! input for wout produced by VMEC

  integer, parameter :: olog  = 50   ! for log output
  integer, parameter :: odbg  = 1000 ! for debug 
  integer, parameter :: oprf  = 1001 ! radial profiles
  integer, parameter :: oshp  = 1002 ! (B,R,Z,phi) data
  integer, parameter :: omtr  = 1003 ! metric data
 
  ! --- grids and 3D data
  real(kind=DP), dimension(0:ntht)   :: theta
  real(kind=DP), dimension(0:nzeta)  :: zeta
  real(kind=DP), dimension(:,:,:), allocatable :: bb, rr, zz, ph, ggb
 
  real(kind=DP)  :: phase, cnorm, dummy
  real(kind=DP)  :: Bax, Rax, aa

  ! --- for cubic spline interpolation 
  real(kind=DP), dimension(:), allocatable :: x1_sp, y1_sp, a_sp, b_sp, c_sp
  integer :: nmax


End Module BZX_header 


Program Read_boozmn_from_BOOZXFORM

  use BZX_header

  implicit none 

  integer :: nsval, ierr, jsize, status
  character(len=50) :: version
  character(len=50) :: varname
  character(len=30) :: file_dbg, file_shape

  integer :: imn, it, iz, js

      open(unit=olog, file=file_log)

! --- input from boozmn (binary)
      open(unit=ibzmn, file=fname_boozmn, form="unformatted")
      read(ibzmn, iostat=ierr) nfp_b, ns_b, aspect_b, rmax_b, rmin_b, betaxis_b
                    if (ierr /= 0) print *, "Read error1 !!"
                    if (ierr /= 0) write(olog,*) "Read error1 !!"

      allocate( iota_b_nu(ns_b), pres_b_nu(ns_b), beta_b_nu(ns_b), phip_b_nu(ns_b),  &
                 phi_b_nu(ns_b), bvco_b_nu(ns_b), buco_b_nu(ns_b) )
      iota_b_nu(1) = 0; pres_b_nu(1) = 0; beta_b_nu(1) = 0; phip_b_nu(1) = 0 
       phi_b_nu(1) = 0; bvco_b_nu(1) = 0; buco_b_nu(1) = 0

      do js = 2, ns_b
         read(ibzmn, iostat=ierr) iota_b_nu(js), pres_b_nu(js), beta_b_nu(js), & 
                                  phip_b_nu(js),  phi_b_nu(js), bvco_b_nu(js), buco_b_nu(js)
                  if (ierr /= 0) print *, "Read error2 !!"
                  if (ierr /= 0) write(olog,*) "Read error2 !!"
      end do

      read(ibzmn, iostat=ierr) mboz_b, nboz_b, mnboz_b, jsize
                  if (ierr /= 0) print *, "Read error3 !!"
                  if (ierr /= 0) write(olog,*) "Read error3 !!"
      read(ibzmn, iostat=js) version, lasym_b
      print*, trim(version)
      print*, "# nfp_b, ns_b, aspect_b, rmax_b, rmin_b, betaxis_b = "
      print*, nfp_b, ns_b, aspect_b, rmax_b, rmin_b, betaxis_b
      print*, "# mboz_b, nboz_b, mnboz_b, jsize, lasym_b = "
      print*, mboz_b, nboz_b, mnboz_b, jsize, lasym_b
      print*, "# fname_tag = ", fname_tag
      if (nzeta == 0) then
        print*, "# nrho, ntheta, nzeta, alpha_fix = "
        print*, nrho, ntht, nzeta, alpha_fix
      else 
        print*, "# nrho, ntheta, nzeta = "
        print*, nrho, ntht, nzeta
      end if
   
      write(olog,*) version
      write(olog,*) "# nfp_b, ns_b, aspect_b, rmax_b, rmin_b, betaxis_b = "
      write(olog,*) nfp_b, ns_b, aspect_b, rmax_b, rmin_b, betaxis_b
      write(olog,*) "# mboz_b, nboz_b, mnboz_b, jsize, lasym_b = "
      write(olog,*) mboz_b, nboz_b, mnboz_b, jsize, lasym_b
      write(olog,*) "# fname_tag = ", fname_tag 
      if (nzeta == 0) then
        write(olog,*) "# nrho, ntheta, nzeta, alpha_fix = "
        write(olog,*) nrho, ntht, nzeta, alpha_fix
      else 
        write(olog,*) "# nrho, ntheta, nzeta = "
        write(olog,*) nrho, ntht, nzeta
      end if

                  if (ierr /= 0) print *, "Read error4 !!"
                  if (ierr /= 0) write(olog,*) "Read error4 !!"
      if (lasym_b) then
            print*, " *** up-down asymmetric configuration *** "
            write(olog,*) " *** up-down asymmetric configuration *** "
      else if (.not.lasym_b) then
            write(olog,*) " *** up-down symmetric configuration *** "
      end if

      allocate( bmnc_b(mnboz_b,ns_b-1), rmnc_b(mnboz_b,ns_b-1),             &
                zmns_b(mnboz_b,ns_b-1), pmns_b(mnboz_b,ns_b-1),             &
                gmnc_b(mnboz_b,ns_b-1), ixm_b(mnboz_b), ixn_b(mnboz_b), jlist(ns_b-1) )
      ixm_b = 0; rmnc_b = 0; zmns_b = 0; pmns_b = 0; bmnc_b = 0; gmnc_b = 0; jlist = 0

      if (lasym_b) then
        allocate( bmns_b(mnboz_b,ns_b-1), rmns_b(mnboz_b,ns_b-1),             &
                  zmnc_b(mnboz_b,ns_b-1), pmnc_b(mnboz_b,ns_b-1),             &
                  gmns_b(mnboz_b,ns_b-1), jlist(ns_b-1) )
        rmns_b = 0; zmnc_b = 0; pmnc_b = 0; bmns_b = 0; gmns_b = 0; jlist = 0
       end if


      read(ibzmn,iostat=ierr) ixn_b(1:mnboz_b), ixm_b(1:mnboz_b)
                  if (ierr /= 0) print *, "Read error5 !!"
                  if (ierr /= 0) write(olog,*) "Read error5 !!"
      
      do js = 1, jsize   ! jsize = ns_b - 1
        read(ibzmn, iostat=ierr) jlist(js)
                  if (ierr /= 0) print *, "Read error6 !!"
                  if (ierr /= 0) write(olog,*) "Read error6 !!"

        if ( js > ns_b .OR. js <= 0 ) cycle
          read(ibzmn,iostat=ierr) bmnc_b(1:mnboz_b,js),                        &
                                  rmnc_b(1:mnboz_b,js), zmns_b(1:mnboz_b,js),  &
                                  pmns_b(1:mnboz_b,js), gmnc_b(1:mnboz_b,js)
                  if (ierr /= 0) print *, "Read error7 !!"
                  if (ierr /= 0) write(olog,*) "Read error7 !!"

        if (lasym_b) then
          read(ibzmn,iostat=ierr) bmns_b(1:mnboz_b,js),                        &
                                  rmns_b(1:mnboz_b,js), zmnc_b(1:mnboz_b,js),  &
                                  pmnc_b(1:mnboz_b,js), gmns_b(1:mnboz_b,js)
                  if (ierr /= 0) print *, "Read error8 !!"
                  if (ierr /= 0) write(olog,*) "Read error8 !!"
        end if
      end do

      close (ibzmn)


! --- input from wout_***.txt (ascii)
      open(iwout,file=fname_wout,action='read')

      read(unit=iwout,fmt="(a)") version
      print*, version
      write(olog,*) version
      do while(.true.)
         read(unit=iwout,fmt="(a)",iostat=status) varname
         !!! print*, varname
         if(status==0) then
            if(trim(varname)==trim(wout_txt)) then
               print *, "data point in wout found."
               write(olog,*) "data point in wout found."
               exit
            end if
         else if(status<0) then
            varname='file end'
            print *, varname
            write(olog,*) varname
         end if
      end do

      read(iwout,*) dummy, B0_p, dummy, dummy, dummy, &
                    Aminor_p, Rmajor_p, volume_p  ! note : B0_p = Volume avg. of |B|
      print*, "B0_p, Aminor_p, Rmajor_p, volume_p = "
      print*, B0_p, Aminor_p, Rmajor_p, volume_p 
 
      write(olog,*) "B0_p, Aminor_p, Rmajor_p, volume_p = "
      write(olog,*) B0_p, Aminor_p, Rmajor_p, volume_p 

      close(iwout)
      

! --- extrapolation to the magnetic axis 

      allocate(bbozc_nu(mnboz_b,ns_b),rbozc_nu(mnboz_b,ns_b),zbozs_nu(mnboz_b,ns_b),pbozs_nu(mnboz_b,ns_b))
      allocate(gbozc_nu(mnboz_b,ns_b))
      if (lasym_b) then
        allocate(bbozs_nu(mnboz_b,ns_b),rbozs_nu(mnboz_b,ns_b),zbozc_nu(mnboz_b,ns_b),pbozc_nu(mnboz_b,ns_b))
        allocate(gbozs_nu(mnboz_b,ns_b))
      end if

      do js = 2, ns_b  ! --- copy array    
        do imn = 1, mnboz_b
          bbozc_nu(imn,js) = bmnc_b(imn,js-1)   ! Non-Uniform grid
          rbozc_nu(imn,js) = rmnc_b(imn,js-1)
          zbozs_nu(imn,js) = zmns_b(imn,js-1)
          pbozs_nu(imn,js) = pmns_b(imn,js-1)
          gbozc_nu(imn,js) = gmnc_b(imn,js-1)
        end do
      end do
      if (lasym_b) then
        do js = 2, ns_b
          do imn = 1, mnboz_b
            bbozs_nu(imn,js) = bmns_b(imn,js-1)
            rbozs_nu(imn,js) = rmns_b(imn,js-1)
            zbozc_nu(imn,js) = zmnc_b(imn,js-1)
            pbozc_nu(imn,js) = pmnc_b(imn,js-1)
            gbozs_nu(imn,js) = gmns_b(imn,js-1)
           end do
        end do
      end if


      do imn = 1, mnboz_b   ! extrapolation
        if ( ixm_b(imn) == 0 ) then
          bbozc_nu(imn,1) = 3._DP*bbozc_nu(imn,2) - 3._DP*bbozc_nu(imn,3) + bbozc_nu(imn,4)
          rbozc_nu(imn,1) = 3._DP*rbozc_nu(imn,2) - 3._DP*rbozc_nu(imn,3) + rbozc_nu(imn,4)
          zbozs_nu(imn,1) = 3._DP*zbozs_nu(imn,2) - 3._DP*zbozs_nu(imn,3) + zbozs_nu(imn,4)
          pbozs_nu(imn,1) = 3._DP*pbozs_nu(imn,2) - 3._DP*pbozs_nu(imn,3) + pbozs_nu(imn,4)
          gbozc_nu(imn,1) = 3._DP*gbozc_nu(imn,2) - 3._DP*gbozc_nu(imn,3) + gbozc_nu(imn,4)
        else
          bbozc_nu(imn,1) = 0._DP 
          rbozc_nu(imn,1) = 0._DP
          zbozs_nu(imn,1) = 0._DP
          pbozs_nu(imn,1) = 0._DP
          gbozc_nu(imn,1) = 0._DP
        end if
      end do
      if (lasym_b) then
        do imn = 1, mnboz_b
          if ( ixm_b(imn) == 0 ) then
            bbozs_nu(imn,1) = 3._DP*bbozs_nu(imn,2) - 3._DP*bbozs_nu(imn,3) + bbozs_nu(imn,4)
            rbozs_nu(imn,1) = 3._DP*rbozs_nu(imn,2) - 3._DP*rbozs_nu(imn,3) + rbozs_nu(imn,4)
            zbozc_nu(imn,1) = 3._DP*zbozc_nu(imn,2) - 3._DP*zbozc_nu(imn,3) + zbozc_nu(imn,4)
            pbozc_nu(imn,1) = 3._DP*pbozc_nu(imn,2) - 3._DP*pbozc_nu(imn,3) + pbozc_nu(imn,4)
            gbozs_nu(imn,1) = 3._DP*gbozs_nu(imn,2) - 3._DP*gbozs_nu(imn,3) + gbozs_nu(imn,4)
          else
            bbozs_nu(imn,1) = 0._DP 
            rbozs_nu(imn,1) = 0._DP
            zbozc_nu(imn,1) = 0._DP
            pbozc_nu(imn,1) = 0._DP
            gbozs_nu(imn,1) = 0._DP
          end if
        end do
      end if
!
       phi_b_nu(1) = 0._DP 
      iota_b_nu(1) = 3._DP*iota_b_nu(2) - 3._DP*iota_b_nu(3) + iota_b_nu(4)
      bvco_b_nu(1) = 0._DP
      buco_b_nu(1) = 3._DP*buco_b_nu(2) - 3._DP*buco_b_nu(3) + buco_b_nu(4)
      pres_b_nu(1) = 3._DP*pres_b_nu(2) - 3._DP*pres_b_nu(3) + pres_b_nu(4)
      beta_b_nu(1) = 3._DP*beta_b_nu(2) - 3._DP*beta_b_nu(3) + beta_b_nu(4)
      phip_b_nu(1) = 0._DP
 
       phi_b_nu(:) =  phi_b_nu(:)/twopi
      phip_b_nu(:) = phip_b_nu(:)/twopi

! --- normalization for vmec calculation and (Bax, Rax, a) for GKV

      ! rescale factor --> cnorm*bbozc (not used ordinary)
      cnorm = Rmajor_p*B0_p/abs(bvco_b_nu(ns_b))  

      Bax = bbozc_nu(1,1)     ! Bax = B_{m=0, n=0}(rho=0, theta=0, zeta=0)
      Rax = rbozc_nu(1,1)     ! Rax = B_{m=0, n=0}(rho=0, theta=0, zeta=0) 
      aa = dsqrt(2._DP*abs(phi_b_nu(ns_b))/Bax)    ! aa = (2*Phi_edge/Bax)^{1/2}
      print*, "Bax, Rax, a (for GKV), Phi_edge = ", Bax, Rax, aa, phi_b_nu(ns_b)
      write(olog,*) "Bax, Rax, a (for GKV), Phi_edge = ", Bax, Rax, aa, phi_b_nu(ns_b)

! --- check RHS or LHS: under construction.... 


! --- rho & theta_b & zeta_b grids and q-profile

      allocate(rho(nrho),rho_nu(ns_b),rho2(ns_b),qq_nu(ns_b))
      do js = 1, ns_b
        rho_nu(js) = dsqrt(phi_b_nu(js)/phi_b_nu(ns_b))   ! Non-Uniform rho grid
        rho2(js) = phi_b_nu(js)/phi_b_nu(ns_b)
      end do
      print*, "rho_nu(1), rho_nu(2), iota_bar(1), iota_bar(2), = "
      print*, rho_nu(1), rho_nu(2), iota_b_nu(1), iota_b_nu(2)
      write(olog,*) "rho_nu(1), rho_nu(2), iota_bar(1), iota_bar(2), = "
      write(olog,*) rho_nu(1), rho_nu(2), iota_b_nu(1), iota_b_nu(2)

      do iz = 0, nzeta
        do it = 0, ntht
          ! theta(it) = real(it,kind=DP)*twopi/real(ntht,kind=DP) - pi  ! theta = [-pi:pi]
          theta(it) = real(Ntheta_gkv,kind=DP)*(real(it,kind=DP)*twopi/real(ntht,kind=DP) - pi)  ! theta = [-N*pi:N*pi]
          zeta(iz)  = real(iz,kind=DP)*twopi/real(nzeta,kind=DP)      ! zeta  = [0:2pi]
          ! zeta(iz)  = real(iz,kind=DP)*twopi/real(nzeta,kind=DP)/real(nfp_b,kind=DP)
        end do
      end do

      do js = 1, nrho
        rho(js) = real(js-1,kind=DP)/real(nrho-1,kind=DP)   ! Uniform rho grid
      end do

      do js = 1, ns_b
         ! qq_nu(js) = 1._DP/iota_b_nu(js)
         qq_nu(js) = 1._DP/abs(iota_b_nu(js))
      end do

! --- interpolation to uniform rho-grids

      allocate(qq(nrho),dqdrho(nrho),shat(nrho),epst(nrho))
      allocate(cug(nrho),cui(nrho),dummy1(nrho),phi_b(nrho),dphidrho(nrho))
      allocate(bbozc(mnboz_b,nrho),rbozc(mnboz_b,nrho),zbozs(mnboz_b,nrho),pbozs(mnboz_b,nrho),gbozc(mnboz_b,nrho))
      allocate(dbbozc(mnboz_b,nrho),drbozc(mnboz_b,nrho),dzbozs(mnboz_b,nrho),dpbozs(mnboz_b,nrho),dummy2(mnboz_b,nrho))
      if (lasym_b) then
        allocate(bbozs(mnboz_b,nrho),rbozs(mnboz_b,nrho),zbozc(mnboz_b,nrho),pbozc(mnboz_b,nrho),gbozs(mnboz_b,nrho))
        allocate(dbbozs(mnboz_b,nrho),drbozs(mnboz_b,nrho),dzbozc(mnboz_b,nrho),dpbozc(mnboz_b,nrho))
      end if


      call cubic_spline_pre(rho_nu,qq_nu,ns_b)
      do js = 1, nrho
        call cubic_spline(rho(js),qq(js),dqdrho(js))
        shat(js) = dqdrho(js)*(rho(js)/qq(js))
        epst(js) = rho(js)*aa/Rax
        ! write(2000, "(256ES24.15e3)") rho_nu(js), rho(js), qq_nu(js), qq(js), dqdrho(js), shat(js), epst(js)
      end do

      ! --- cug: B_zeta  (covariant zeta comp. of B, or toroidal current func.)
      call cubic_spline_pre(rho_nu,bvco_b_nu,ns_b)
      do js = 1, nrho
        call cubic_spline(rho(js),cug(js),dummy1(js))
      end do
      ! --- cui: B_theta (covariant theta comp. of B, or poloidal current func.)
      call cubic_spline_pre(rho_nu,buco_b_nu,ns_b)
      do js = 1, nrho
        call cubic_spline(rho(js),cui(js),dummy1(js))
      end do
      call cubic_spline_pre(rho_nu,phi_b_nu,ns_b)
      do js = 1, nrho
        call cubic_spline(rho(js),phi_b(js),dphidrho(js))
      end do

      do imn = 1, mnboz_b
        call cubic_spline_pre(rho_nu,bbozc_nu(imn,:),ns_b)
        do js = 1, nrho
          call cubic_spline(rho(js),bbozc(imn,js),dbbozc(imn,js))
        end do
      end do

      do imn = 1, mnboz_b
        call cubic_spline_pre(rho_nu,rbozc_nu(imn,:),ns_b)
        do js = 1, nrho
          call cubic_spline(rho(js),rbozc(imn,js),drbozc(imn,js))
        end do
      end do
  
      do imn = 1, mnboz_b
        call cubic_spline_pre(rho_nu,zbozs_nu(imn,:),ns_b)
        do js = 1, nrho
          call cubic_spline(rho(js),zbozs(imn,js),dzbozs(imn,js))
        end do
      end do

      do imn = 1, mnboz_b
        call cubic_spline_pre(rho_nu,pbozs_nu(imn,:),ns_b)
        do js = 1, nrho
          call cubic_spline(rho(js),pbozs(imn,js),dpbozs(imn,js))
          ! write(2100, "(256ES24.15e3)") rho_nu(js), rho(js), pbozs_nu(imn,js), pbozs(imn,js), dpbozs(imn,js)
        end do
          ! write(2100,*)
          ! write(2100,*)
      end do

      ! --- gbozc: rootg_boz/dphidrho 
      do imn = 1, mnboz_b
        call cubic_spline_pre(rho_nu,gbozc_nu(imn,:),ns_b)
        do js = 1, nrho
          call cubic_spline(rho(js),gbozc(imn,js),dummy2(imn,js))
        end do
      end do
   

      if (lasym_b) then
        do imn = 1, mnboz_b
          call cubic_spline_pre(rho_nu,bbozs_nu(imn,:),ns_b)
          do js = 1, nrho
            call cubic_spline(rho(js),bbozs(imn,js),dbbozs(imn,js))
          end do
        end do

        do imn = 1, mnboz_b
          call cubic_spline_pre(rho_nu,rbozs_nu(imn,:),ns_b)
          do js = 1, nrho
            call cubic_spline(rho(js),rbozs(imn,js),drbozs(imn,js))
          end do
        end do
  
        do imn = 1, mnboz_b
          call cubic_spline_pre(rho_nu,zbozc_nu(imn,:),ns_b)
          do js = 1, nrho
            call cubic_spline(rho(js),zbozc(imn,js),dzbozc(imn,js))
          end do
        end do

        do imn = 1, mnboz_b
          call cubic_spline_pre(rho_nu,pbozc_nu(imn,:),ns_b)
          do js = 1, nrho
            call cubic_spline(rho(js),pbozc(imn,js),dpbozc(imn,js))
          end do
        end do

        do imn = 1, mnboz_b
          call cubic_spline_pre(rho_nu,gbozs_nu(imn,:),ns_b)
          do js = 1, nrho
            call cubic_spline(rho(js),gbozs(imn,js),dummy2(imn,js))
          end do
        end do
      end if
     

! --- B(rho,theta_b,zeta_b), R(rho,theta_b,zeta_b), Z(rho,theta_b,zeta_b), phi(rho,theta_b,zeta_b)

      allocate( bb(1:nrho,0:ntht,0:nzeta), rr(1:nrho,0:ntht,0:nzeta), & 
                zz(1:nrho,0:ntht,0:nzeta), ph(1:nrho,0:ntht,0:nzeta) )
      allocate( ggb(1:nrho,0:ntht,0:nzeta) )

      allocate( dbb_drho(1:nrho,0:ntht,0:nzeta), drr_drho (1:nrho,0:ntht,0:nzeta), &
                dzz_drho(1:nrho,0:ntht,0:nzeta), dph_drho (1:nrho,0:ntht,0:nzeta), &
                dbb_dtht(1:nrho,0:ntht,0:nzeta), dbb_dzeta(1:nrho,0:ntht,0:nzeta), & 
                drr_dtht(1:nrho,0:ntht,0:nzeta), drr_dzeta(1:nrho,0:ntht,0:nzeta), &
                dzz_dtht(1:nrho,0:ntht,0:nzeta), dzz_dzeta(1:nrho,0:ntht,0:nzeta), &
                dph_dtht(1:nrho,0:ntht,0:nzeta), dph_dzeta(1:nrho,0:ntht,0:nzeta) )

      do js = 1, nrho
        print*, "js-loop num. = ", js
        write(olog,*) "js-loop num. = ", js
        bb(js,:,:) = 0._DP
        rr(js,:,:) = 0._DP
        zz(js,:,:) = 0._DP
        ph(js,:,:) = 0._DP
        ggb(js,:,:) = 0._DP
        dbb_drho(js,:,:) = 0._DP
        drr_drho(js,:,:) = 0._DP
        dzz_drho(js,:,:) = 0._DP
        dph_drho(js,:,:) = 0._DP
        dbb_dtht(js,:,:) = 0._DP
        drr_dtht(js,:,:) = 0._DP
        dzz_dtht(js,:,:) = 0._DP
        dph_dtht(js,:,:) = 0._DP
        dbb_dzeta(js,:,:) = 0._DP
        drr_dzeta(js,:,:) = 0._DP
        dzz_dzeta(js,:,:) = 0._DP
        dph_dzeta(js,:,:) = 0._DP
        do iz = 0, nzeta
          do it = 0, ntht

             do imn = 1, mnboz_b
              if ( nzeta == 0 ) then 
                zeta(iz) = qq(js)*theta(it) + alpha_fix
                phase = -real(ixn_b(imn),kind=DP)*zeta(iz) + real(ixm_b(imn),kind=DP)*theta(it)
               ! phase = real(ixn_b(imn),kind=DP)*zeta(iz) + real(ixm_b(imn),kind=DP)*theta(it)
              else
                phase = -real(ixn_b(imn),kind=DP)*zeta(iz) + real(ixm_b(imn),kind=DP)*theta(it)
               ! phase = real(ixn_b(imn),kind=DP)*zeta(iz) + real(ixm_b(imn),kind=DP)*theta(it)
              end if
              bb(js,it,iz)  =  bb(js,it,iz) + bbozc(imn,js)*dcos(phase)    ! --- only for asymmetric parts
              rr(js,it,iz)  =  rr(js,it,iz) + rbozc(imn,js)*dcos(phase)
              zz(js,it,iz)  =  zz(js,it,iz) + zbozs(imn,js)*dsin(phase)
              ph(js,it,iz)  =  ph(js,it,iz) + pbozs(imn,js)*dsin(phase)
              ggb(js,it,iz) = ggb(js,it,iz) + gbozc(imn,js)*dcos(phase)

              dbb_drho(js,it,iz) = dbb_drho(js,it,iz) + dbbozc(imn,js)*dcos(phase)
              drr_drho(js,it,iz) = drr_drho(js,it,iz) + drbozc(imn,js)*dcos(phase)
              dzz_drho(js,it,iz) = dzz_drho(js,it,iz) + dzbozs(imn,js)*dsin(phase)
              dph_drho(js,it,iz) = dph_drho(js,it,iz) + dpbozs(imn,js)*dsin(phase)

              dbb_dtht(js,it,iz) = dbb_dtht(js,it,iz) - real(ixm_b(imn),kind=DP)*bbozc(imn,js)*dsin(phase)
              drr_dtht(js,it,iz) = drr_dtht(js,it,iz) - real(ixm_b(imn),kind=DP)*rbozc(imn,js)*dsin(phase)
              dzz_dtht(js,it,iz) = dzz_dtht(js,it,iz) + real(ixm_b(imn),kind=DP)*zbozs(imn,js)*dcos(phase)
              dph_dtht(js,it,iz) = dph_dtht(js,it,iz) + real(ixm_b(imn),kind=DP)*pbozs(imn,js)*dcos(phase)

              dbb_dzeta(js,it,iz) = dbb_dzeta(js,it,iz) + real(ixn_b(imn),kind=DP)*bbozc(imn,js)*dsin(phase)
              drr_dzeta(js,it,iz) = drr_dzeta(js,it,iz) + real(ixn_b(imn),kind=DP)*rbozc(imn,js)*dsin(phase)
              dzz_dzeta(js,it,iz) = dzz_dzeta(js,it,iz) - real(ixn_b(imn),kind=DP)*zbozs(imn,js)*dcos(phase)
              dph_dzeta(js,it,iz) = dph_dzeta(js,it,iz) - real(ixn_b(imn),kind=DP)*pbozs(imn,js)*dcos(phase)
            end do

            ph(js,it,iz) = zeta(iz) + ph(js,it,iz)   
            dph_drho(js,it,iz) =  dph_drho(js,it,iz)
            dph_dtht(js,it,iz) =  dph_dtht(js,it,iz)
            dph_dzeta(js,it,iz) = 1._DP + dph_dzeta(js,it,iz)
          end do
        end do
      end do


! --- calculation and output of metric tensor

      call  metric_boozer(nrho) 

! --- output of geom data
      open(unit=oprf, file=file_prf)
      do js = 1, nrho
         write(unit=oprf, fmt="(I5, 256E24.15e3)") js, rho(js),      & 
         qq(js), shat(js), epst(js), cug(js), cui(js), phi_b(js), dphidrho(js)
      end do
      close(oprf)

      open(unit=oshp, file=file_out1)
      do iz = 0, nzeta
      do it = 0, ntht
      do js = 1, nrho
        write(unit=oshp,fmt="(256ES24.15e3)") rho(js), theta(it), zeta(iz), & 
                                              bb(js,it,iz), &
                                              rr(js,it,iz), &
                                              zz(js,it,iz), &
                                              ph(js,it,iz)
      end do
      write(unit=oshp,fmt=*)
      end do
      write(unit=oshp,fmt=*)
      write(unit=oshp,fmt=*)
      end do
      close(oshp)

      open(unit=oshp, file=file_out2)
      do iz = 0, nzeta
      do js = 1, nrho
      do it = 0, ntht
        write(unit=oshp,fmt="(256ES24.15e3)") rho(js), theta(it), zeta(iz), & 
                                              bb(js,it,iz), &
                                              rr(js,it,iz), &
                                              zz(js,it,iz), &
                                              ph(js,it,iz)
      end do
      write(unit=oshp,fmt=*)
      end do
      write(unit=oshp,fmt=*)
      write(unit=oshp,fmt=*)
      end do
      close(oshp)

      open(unit=oshp, file=file_out3)
      do js = 1, nrho
      do iz = 0, nzeta
      do it = 0, ntht
        write(unit=oshp,fmt="(256ES24.15e3)") rho(js), theta(it), zeta(iz), & 
                                              bb(js,it,iz), &
                                              rr(js,it,iz), &
                                              zz(js,it,iz), &
                                              ph(js,it,iz)
      end do
      write(unit=oshp,fmt=*)
      end do
      write(unit=oshp,fmt=*)
      write(unit=oshp,fmt=*)
      end do
      close(oshp)


! --- output for checking the read data

      file_dbg = "./check/check1.dat"
      open(unit=odbg, file=file_dbg)
      do js = 1, ns_b
         write(unit=odbg, fmt="(I5, 256E24.15e3)") js, rho_nu(js),      & 
         iota_b_nu(js), qq_nu(js), pres_b_nu(js), beta_b_nu(js), phip_b_nu(js), phi_b_nu(js),  &
         bvco_b_nu(js), buco_b_nu(js), dsqrt(phi_b_nu(js)/phi_b_nu(ns_b))
      end do
      close(odbg)

      file_dbg = "./check/check2.dat"
      open(unit=odbg, file=file_dbg)
      if (lasym_b) then
        do js = 1, ns_b-1
          do imn = 1, mnboz_b
            write(unit=odbg, fmt="(2I5, 256ES24.15e3)") js, imn, bmnc_b(imn,js), &
                                                 rmnc_b(imn,js), zmns_b(imn,js), &
                                                 pmns_b(imn,js), gmnc_b(imn,js), &
                                                 bmns_b(imn,js),                 &
                                                 rmns_b(imn,js), zmnc_b(imn,js), &
                                                 pmnc_b(imn,js), gmns_b(imn,js)
         
          end do
          write(unit=odbg,fmt=*)
        end do
      else 
        do js = 1, ns_b-1
          do imn = 1, mnboz_b
            write(unit=odbg, fmt="(2I5, 256ES24.15e3)") js, imn, bmnc_b(imn,js), &
                                                 rmnc_b(imn,js), zmns_b(imn,js), &
                                                 pmns_b(imn,js), gmnc_b(imn,js)
          end do
          write(unit=odbg,fmt=*)
        end do
      end if
      close(odbg)

      file_dbg = "./check/check3.dat"
      open(unit=odbg, file=file_dbg)
      do imn = 1, mnboz_b
        write(unit=odbg, fmt="(3I5, ES24.15e3)") imn, ixn_b(imn), ixm_b(imn), bmnc_b(imn,(ns_b-1)/2)
      end do
      close(odbg)


  print *, "DONE!!"
  write(olog,*) "DONE!!"
  stop

End Program Read_boozmn_from_BOOZXFORM


Subroutine metric_boozer (nss) 

  use BZX_header

  implicit none

  integer, intent(in)    :: nss
  real(kind=DP), dimension(:,:,:,:,:), allocatable ::  ggup_boz,  ggdn_boz
  real(kind=DP), dimension(:,:,:),     allocatable ::  ggsq_boz, rootg_boz
  real(kind=DP), dimension(:,:,:),     allocatable ::  rootg_boz0, rootg_boz1 ! another definition
  integer :: js, it, iz, it2, asym_flg
  character(len=30) :: file_dbg

! --- indices 1 to 3 corredpond to rho, theta, zeta
      allocate(ggup_boz(1:nss,0:ntht,0:nzeta,1:3,1:3),ggdn_boz(1:nss,0:ntht,0:nzeta,1:3,1:3)) 
      allocate(ggsq_boz(1:nss,0:ntht,0:nzeta),rootg_boz(1:nss,0:ntht,0:nzeta)) 
      allocate(rootg_boz0(1:nss,0:ntht,0:nzeta),rootg_boz1(1:nss,0:ntht,0:nzeta)) 
  
  
      do iz = 0, nzeta
        do it = 0, ntht
          do js = 1, nss
  
! --- Co-variant metric tensor
            ggdn_boz(js,it,iz,1,1) =   drr_drho(js,it,iz) * drr_drho(js,it,iz)               &
                                     + dzz_drho(js,it,iz) * dzz_drho(js,it,iz)               &
                                     + rr(js,it,iz)**2    * dph_drho(js,it,iz)**2          
            
            ggdn_boz(js,it,iz,1,2) =   drr_drho(js,it,iz) * drr_dtht(js,it,iz)               &
                                     + dzz_drho(js,it,iz) * dzz_dtht(js,it,iz)               &
                                     + rr(js,it,iz)**2    * dph_drho(js,it,iz)               & 
                                                          * dph_dtht(js,it,iz)                  
    
            ggdn_boz(js,it,iz,1,3) =   drr_drho(js,it,iz) * drr_dzeta(js,it,iz)              &
                                     + dzz_drho(js,it,iz) * dzz_dzeta(js,it,iz)              &
                                     + rr(js,it,iz)**2    * dph_drho( js,it,iz)              & 
                                                          * dph_dzeta(js,it,iz)                  
    
            ggdn_boz(js,it,iz,2,1) = ggdn_boz(js,it,iz,1,2)
    
    
            ggdn_boz(js,it,iz,2,2) =   drr_dtht(js,it,iz) * drr_dtht(js,it,iz)               &
                                     + dzz_dtht(js,it,iz) * dzz_dtht(js,it,iz)               &
                                     + rr(js,it,iz)**2    * dph_dtht(js,it,iz)               & 
                                                          * dph_dtht(js,it,iz)                  
    
            ggdn_boz(js,it,iz,2,3) =   drr_dtht(js,it,iz) * drr_dzeta(js,it,iz)              &
                                     + dzz_dtht(js,it,iz) * dzz_dzeta(js,it,iz)              &
                                     + rr(js,it,iz)**2    * dph_dtht( js,it,iz)              & 
                                                          * dph_dzeta(js,it,iz)                  
    
            ggdn_boz(js,it,iz,3,1) = ggdn_boz(js,it,iz,1,3)
    
            ggdn_boz(js,it,iz,3,2) = ggdn_boz(js,it,iz,2,3)
    
            ggdn_boz(js,it,iz,3,3) =   drr_dzeta(js,it,iz) * drr_dzeta(js,it,iz)             &
                                     + dzz_dzeta(js,it,iz) * dzz_dzeta(js,it,iz)             &
                                     + rr(js,it,iz)**2     * dph_dzeta(js,it,iz)             & 
                                                           * dph_dzeta(js,it,iz)     
               
! --- Squared Jacobian
            ggsq_boz(js,it,iz) =  ggdn_boz(js,it,iz,1,1) * ggdn_boz(js,it,iz,2,2) * ggdn_boz(js,it,iz,3,3)   &
                                + ggdn_boz(js,it,iz,1,2) * ggdn_boz(js,it,iz,2,3) * ggdn_boz(js,it,iz,3,1)   &
                                + ggdn_boz(js,it,iz,1,3) * ggdn_boz(js,it,iz,2,1) * ggdn_boz(js,it,iz,3,2)   &
                                - ggdn_boz(js,it,iz,1,3) * ggdn_boz(js,it,iz,2,2) * ggdn_boz(js,it,iz,3,1)   &
                                - ggdn_boz(js,it,iz,1,2) * ggdn_boz(js,it,iz,2,1) * ggdn_boz(js,it,iz,3,3)   &
                                - ggdn_boz(js,it,iz,1,1) * ggdn_boz(js,it,iz,2,3) * ggdn_boz(js,it,iz,3,2)   
    
! --- Jacobian: rootg = sqrt(g)
            rootg_boz(js,it,iz) = dsqrt( ggsq_boz(js,it,iz) )
            rootg_boz0(js,it,iz) = dphidrho(js)*( cug(js) + cui(js)/qq(js) )/bb(js,it,iz)**2
            rootg_boz1(js,it,iz) = dphidrho(js)*ggb(js,it,iz)
            ! --- cug: B_zeta  (covariant zeta comp. of B, or toroidal current func.)
            ! --- cui: B_theta (covariant theta comp. of B, or poloidal current func.)
            ! --- ggb: rootg_boz/dphidrho 

! --- Contra-variant metric tensor
            ggup_boz(js,it,iz,1,1) = (   ggdn_boz(js,it,iz,2,2) * ggdn_boz(js,it,iz,3,3)                     &
                                       - ggdn_boz(js,it,iz,2,3) * ggdn_boz(js,it,iz,3,2) ) / ggsq_boz(js,it,iz)         
            ggup_boz(js,it,iz,1,2) = (   ggdn_boz(js,it,iz,1,3) * ggdn_boz(js,it,iz,3,2)                     &
                                       - ggdn_boz(js,it,iz,1,2) * ggdn_boz(js,it,iz,3,3) ) / ggsq_boz(js,it,iz)         
            ggup_boz(js,it,iz,1,3) = (   ggdn_boz(js,it,iz,1,2) * ggdn_boz(js,it,iz,2,3)                     &
                                       - ggdn_boz(js,it,iz,1,3) * ggdn_boz(js,it,iz,2,2) ) / ggsq_boz(js,it,iz)         
            ggup_boz(js,it,iz,2,1) = ggup_boz(js,it,iz,1,2)
            ggup_boz(js,it,iz,2,2) = (   ggdn_boz(js,it,iz,1,1) * ggdn_boz(js,it,iz,3,3)                     &
                                       - ggdn_boz(js,it,iz,1,3) * ggdn_boz(js,it,iz,3,1) ) / ggsq_boz(js,it,iz)         
            ggup_boz(js,it,iz,2,3) = (   ggdn_boz(js,it,iz,1,3) * ggdn_boz(js,it,iz,2,1)                     &
                                       - ggdn_boz(js,it,iz,1,1) * ggdn_boz(js,it,iz,2,3) ) / ggsq_boz(js,it,iz)         
            ggup_boz(js,it,iz,3,1) = ggup_boz(js,it,iz,1,3)
            ggup_boz(js,it,iz,3,2) = ggup_boz(js,it,iz,2,3)
            ggup_boz(js,it,iz,3,3) = (   ggdn_boz(js,it,iz,1,1) * ggdn_boz(js,it,iz,2,2)                     &
                                       - ggdn_boz(js,it,iz,1,2) * ggdn_boz(js,it,iz,2,1) ) / ggsq_boz(js,it,iz)         
  
          end do
        end do
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  All quantities have not been normalized for GKV here. 
!  Polidal angle range is [-pi:pi], but [0:2pi] for toroidal one. 
!
!  List of output data: 
!            rho: radial coord.,                rho = sqrt(phi/phi_edge) [phi: toloidal flux / 2pi]
!             qq: safety factor,                q(rho)
!           shat: local magnetic shear,         s_hat(rho) = (rho/q)*dq/drho
!           epst: local inverse aspect ratio,   eps(rho) = rho*a/Rax, a is defined as a = sqrt(2phi_edge/Bax) 
!            cug: covariant B-field comp.,      B_zeta(rho)  (or toroidal current function)
!            cui: covariant B-field comp.,      B_theta(rho) (or poloidal current function)
!       dphidrho: radial deriv of phi,          dphi/drho(rho)
!             bb: B-field intensity,            B(rho,theta,zeta) [T]
!             rr: B-field intensity,            R(rho,theta,zeta) [m]
!             zz: B-field intensity,            Z(rho,theta,zeta) [m]
!             ph: B-field intensity,            phi(rho,theta,zeta) in VMEC coord. [radian]
!          dX_dY: derivatives in Booz. coord.,  dX/dY(rho,theta,zeta), X={bb,rr,zz,ph}, Y={rho,tht,zeta}
!      rootg_boz: Jacobian,                     sqrt(g_boozer) [rho,theta,zeta]
!           ggdn: covariant metric comp.,       g_ij [rho,theta,zeta] (i,j) = {rho,theta,zeta}
!           ggup: contravariant metric comp.,   g^ij [rho,theta,zeta] (i,j) = {rho,theta,zeta}
!          nfp_b: stellarator period            nfp_b = 10(LHD), = 5(W-7X), = 4(Heliotron-J) etc.
!           mboz: the number of modes for m     poloidal mode number in Boozer coord. (specified by in_booz)
!           nboz: the number of modes for n     toroidal mode number in Boozer coord. (specified by in_booz)
!          mnboz: the total number of modes     mnboz = nboz+1 + (mboz-1)*(1+2*nboz)
!          ixm_b: poloidal mode number          m(i), i=1, mnboz_b
!          ixn_b: toroidal mode number          n(i), i=1, mnboz_b
!          bbozc: Fourier spectrum of B         B_{n,m}(i), i=1, mnboz_b
!            Bax: |B| at the axis               B_{m=0,n=0}(rho=0,theta=0,zeta=0) [T]
!            Rax: |B| at the axis               R_{m=0,n=0}(rho=0,theta=0,zeta=0) [m]
!             aa: plasma radius                 a = sqrt(2phi_edge/Bax)  [m]
!       asym_flg: flag for up-down asym.        =1: up-down asymmetric fields, =0: up-down symmetric(default)
!      alpha_fix: field line label              alpha_fix = zeta - q*theta
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! --- output of metric data in Boozer coord.
! --- binary output

      if (lasym_b) then 
        asym_flg = 1
      else if (.not.lasym_b) then
        asym_flg = 0
      end if
      
      open(unit=omtr, file=file_out4, form="unformatted")
      write(unit=omtr) nfp_b, nss, ntht, nzeta, mnboz_b, mboz_b, nboz_b, Rax, Bax, aa, volume_p, asym_flg, alpha_fix
      write(unit=omtr) rho, theta, zeta, qq, shat, epst, bb, rootg_boz, rootg_boz0, ggup_boz, & 
                       dbb_drho, dbb_dtht, dbb_dzeta
      write(unit=omtr) rr, zz, ph, bbozc, ixn_b, ixm_b
      if (lasym_b) write(unit=omtr) bbozs
      close(omtr)
  
          ! --- for debug
          file_dbg = "./check/mag_check.dat"
          open(unit=odbg, file=file_dbg)
          do js = 1, nss
          do iz = 0, nzeta
          do it = 0, ntht
             write(odbg,"(256ES24.15e3)") rho(js), theta(it), zeta(iz), qq(js), & 
                                          shat(js), epst(js),                   & 
                                          cui(js), cug(js), dphidrho(js),       &
                                          rootg_boz0(js,it,iz), &    ! 10 
                                          rootg_boz1(js,it,iz), &    ! 11
                                          rootg_boz(js,it,iz), &     ! 12
                                          ggsq_boz(js,it,iz),  &     ! 13
                                          ggdn_boz(js,it,iz,1,1), &  ! 14
                                          ggdn_boz(js,it,iz,1,2), &  ! 15
                                          ggdn_boz(js,it,iz,1,3), &  ! 16
                                          ggdn_boz(js,it,iz,2,1), &  ! 17
                                          ggdn_boz(js,it,iz,2,2), &  ! 18
                                          ggdn_boz(js,it,iz,2,3), &  ! 19
                                          ggdn_boz(js,it,iz,3,1), &  ! 20
                                          ggdn_boz(js,it,iz,3,2), &  ! 21
                                          ggdn_boz(js,it,iz,3,3), &  ! 22 
                                          ggup_boz(js,it,iz,1,1), &  ! 23
                                          ggup_boz(js,it,iz,1,2), &  ! 24
                                          ggup_boz(js,it,iz,1,3), &  ! 25
                                          ggup_boz(js,it,iz,2,1), &  ! 26
                                          ggup_boz(js,it,iz,2,2), &  ! 27
                                          ggup_boz(js,it,iz,2,3), &  ! 28
                                          ggup_boz(js,it,iz,3,1), &  ! 29
                                          ggup_boz(js,it,iz,3,2), &  ! 30
                                          ggup_boz(js,it,iz,3,3), &  ! 31
                                          bb(js,it,iz),           &  ! 32
                                          rr(js,it,iz),           &  ! 33
                                          zz(js,it,iz),           &  ! 34
                                          ph(js,it,iz),           &  ! 35
                                          dbb_drho(js,it,iz),     &  ! 36
                                          drr_drho(js,it,iz),     &  ! 37 
                                          dzz_drho(js,it,iz),     &  ! 38
                                          dph_drho(js,it,iz),     &  ! 39
                                          dbb_dtht(js,it,iz),     &  ! 40 
                                          drr_dtht(js,it,iz),     &  ! 41 
                                          dzz_dtht(js,it,iz),     &  ! 42 
                                          dph_dtht(js,it,iz),     &  ! 43
                                          dbb_dzeta(js,it,iz),    &  ! 44
                                          drr_dzeta(js,it,iz),    &  ! 45
                                          dzz_dzeta(js,it,iz),    &  ! 46
                                          dph_dzeta(js,it,iz)        ! 47
          end do
             write(odbg,*) 
          end do
             write(odbg,*) 
             write(odbg,*) 
          end do
          close(odbg)

    return

End Subroutine metric_boozer


Subroutine cubic_spline_pre(xp,yp,n)

  use BZX_header
  implicit none

  real(kind=DP), dimension (1:n), intent(in) :: xp, yp
  real(kind=DP), dimension(:), allocatable :: ah, bh, ch, dh
  real(kind=DP) :: h1, h2, x0
  integer       :: n, i, i1, i2

        if (allocated(x1_sp)) deallocate ( x1_sp, y1_sp, a_sp, b_sp, c_sp )
        allocate ( x1_sp(1:n), y1_sp(1:n-1), a_sp(1:n-1), b_sp(1:n-1), c_sp(1:n) )
        allocate ( ah(1:n), bh(1:n), ch(1:n), dh(1:n) )

        h1 = xp(2)-xp(1)
        bh(1) = 2*h1
        ch(1) = h1
        dh(1) = 3*(yp(2)-yp(1))

        do i = 2, n-1
          h1 = xp(i)-xp(i-1)
          h2 = xp(i+1)-xp(i)
          ah(i) = h2
          bh(i) = 2*(xp(i+1)-xp(i-1))
          ch(i) = h1
          dh(i) = 3*((yp(i)-yp(i-1))*h2/h1+(yp(i+1)-yp(i))*h1/h2)
        end do

        h1 = xp(n)-xp(n-1)
        ah(n) = h1
        bh(n) = 2*h1
        dh(n) = 3*(yp(n)-yp(n-1))
        call tridiagonal_matrix(ah,bh,ch,dh,n)

        do i = 1, n-1
          h1 = xp(i+1)-xp(i);     h2 = h1*h1
          x1_sp(i) = xp(i)
          y1_sp(i) = yp(i)
          b_sp(i) = (3*(yp(i+1)-yp(i))-(c_sp(i+1)+2*c_sp(i))*h1)/h2
          a_sp(i) = (c_sp(i+1)-c_sp(i)-2*b_sp(i)*h1)/(3*h2)
        end do
        x1_sp(n) = xp(n)
        nmax  = n

        deallocate ( ah, bh, ch, dh )

    return

End Subroutine cubic_spline_pre


Subroutine cubic_spline(x,y,dydx)

  use BZX_header
  implicit none

  real(kind=DP) :: x, y, dydx
  real(kind=DP) :: x0
  integer :: i, i1, i2

      if (x <= x1_sp(2)) then
        i1 = 1
        if (x < x1_sp(1)) print *,'Out of bounds -- cubic_spline',x
      else if (x >= x1_sp(nmax-1)) then
        i1 = nmax-1
        if (x > x1_sp(nmax)) print *,'Out of bounds -- cubic_spline',x
      else
        i1 = 2;  i2 = nmax-1
        do while (i2-i1 > 1)
           i = (i1+i2)/2
           if (x < x1_sp(i)) then
              i2 = i
           else
              i1 = i
           end if
        end do
      end if

      x0 = x - x1_sp(i1)
      y = ((a_sp(i1)*x0 + b_sp(i1))*x0 + c_sp(i1))*x0 + y1_sp(i1)
      dydx = (3.d0*a_sp(i1)*x0 + 2.d0*b_sp(i1))*x0 + c_sp(i1)

  return

End Subroutine cubic_spline


Subroutine tridiagonal_matrix(a1,b1,c1,d1,n)

  use BZX_header
  implicit none

  real(kind=DP), dimension(1:n), intent(in) :: a1, b1, c1, d1
  real(kind=DP), dimension(1:n) :: G, H
  real(kind=DP)                 :: den
  integer :: i, n

      G(1) = -c1(1)/b1(1)
      H(1) = d1(1)/b1(1)

      do i = 2, n
        den  = 1/(b1(i) + a1(i)*G(i-1))
        G(i) = -c1(i)*den
        H(i) = (d1(i) - a1(i)*H(i-1))*den
      end do

      c_sp(n) = H(n)
      do i = n-1, 1, -1
        c_sp(i) = G(i)*c_sp(i+1) + H(i)
      end do

  return

End Subroutine tridiagonal_matrix
