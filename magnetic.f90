!=================================================================================================
! module_magnetic: Contains variables for magnetic forces
!-------------------------------------------------------------------------------------------------
module module_magnetic
  real(8), allocatable, dimension(:,:,:) :: phimag,chip1
!  real(8), parameter :: mu_0 = 1.256637d-6 ! vacuum value
  real(8) :: mu1mag,mu2mag
  character(len=6) :: init_phi_type
  real(8):: init_phi_param

contains

subroutine init_mag
  implicit none
  call ReadMagParameters
  call initialize_ferro
  call phimaginitialize
end subroutine init_mag

!Read parameters
subroutine ReadMagParameters
  use module_grid
  use module_2phase
  use module_IO
  use module_flow
  implicit none
  integer :: ierr
  namelist /magparameters/ mu1mag,mu2mag,init_phi_type,init_phi_param

  open(unit=10, file='inputmagnetic', status='old', action='read', iostat=ierr)
  if (ierr .ne. 0) call err_no_out_dir("ReadMagParameters: error opening 'inputferromag' file --- perhaps it does not exist ?")
  read(10,NML=magparameters)
  close(10)

end subroutine ReadMagParameters

! initialize mag variables
subroutine initialize_ferro
  use module_grid
  use module_flow    
  implicit none
  
  ! Allocate variables
  allocate ( phimag(imin:imax,jmin:jmax,kmin:kmax),&
       chip1(imin:imax,jmin:jmax,kmin:kmax) )

  phimag=0.0d0; chip1 = 0.0d0
end subroutine initialize_ferro

subroutine phimaginitialize

  use module_grid

  implicit none
  !real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(out) :: phimag
  integer :: i,j,k

    if (init_phi_type=='linear') then
      do k=kmin,kmax-1;  do j=jmin,jmax-1;  do i=imax,imax-1
        phimag(i,j,k) = 1.0d0 + init_phi_param*z(k)
      enddo;  enddo;  enddo
    else
     call pariserror('unknown inital phi type')
    endif
   
end subroutine phimaginitialize

!=================================================================================================
! Calculates the magnetic surface force in the momentum equations and adds it to du,dv,dw
!-------------------------------------------------------------------------------------------------
subroutine surfaceMagForce(du,dv,dw,rho)
!  use module_solid
  use module_grid
  use module_vof
  use module_2phase
  use module_surface_tension
  use module_tmpvar
  use module_timer
  
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: du, dv, dw, rho
  real(8) :: deltax
  real(8) ::  Hn_sq, Ht_sq, H(3),Hn(3),Ht(3)
  real(8) :: stencil3x3(-1:1,-1:1,-1:1),mxyz(3),mxyzp(3),mudiff,muratio
  integer :: i,i1,j,k,n,ierr
  integer :: i0,j0,k0
  
  deltax = dx(nx)
  mudiff = mu1mag - mu2mag
  muratio = mu2mag/mu1mag
  
!  call my_timer(7)
  do k=ks,ke;  do j=js,je; do i=is,ieu
     if(abs(cvof(i+1,j,k)-cvof(i,j,k))>EPSC/1d1) then  
       ! there is a non-zero grad H (H Heaviside function)

       !Compute normal vectors and contributions at (i,j,k)
        do k0=-1,1; do j0=-1,1; do i0=-1,1
            stencil3x3(i0,j0,k0) = cvof(i+i0,j+j0,k+k0)
        enddo;enddo;enddo
        call mycs(stencil3x3,mxyz)
        
        H(1) = - 0.5D0*(phimag(i+1,j,k) - phimag(i-1,j,k))/dxh(i)
        H(2) = - 0.5D0*(phimag(i,j+1,k) - phimag(i,j-1,k))/dyh(j)
        H(3) = - 0.5D0*(phimag(i,j,k+1) - phimag(i,j,k-1))/dzh(k)

        !Compute H-normal and H-tangential components
        do i1=1,3
            Hn(i1) = (H(1)*mxyz(1)+H(2)*mxyz(2)+H(3)*mxyz(3))*mxyz(i1)           
            Ht(i1) = (H(i1)-Hn(i1))
        enddo

        !Compute Hnorm and Htan magnitudes
        Hn_sq = Hn(1)*Hn(1)+Hn(2)*Hn(2)+Hn(3)*Hn(3)
        Ht_sq = Ht(1)*Ht(1)+Ht(2)*Ht(2)+Ht(3)*Ht(3)

        !Compute normal vectors and contributions at (i+1,j,k)
        do k0=-1,1; do j0=-1,1; do i0=-1,1
            stencil3x3(i0,j0,k0) = cvof(i+1+i0,j+j0,k+k0)
        enddo;enddo;enddo
        call mycs(stencil3x3,mxyzp)
        
        H(1) = - 0.5D0*(phimag(i+2,j,k) - phimag(i,j,k))/dxh(i+1)
        H(2) = - 0.5D0*(phimag(i+1,j+1,k) - phimag(i+1,j-1,k))/dyh(j)
        H(3) = - 0.5D0*(phimag(i+1,j,k+1) - phimag(i+1,j,k-1))/dzh(k)

        !Compute H-normal and H-tangential components
        do i1=1,3
            Hn(i1) = (H(1)*mxyzp(1)+H(2)*mxyzp(2)+H(3)*mxyzp(3))*mxyzp(i1)           
            Ht(i1) = (H(i1)-Hn(i1))
        enddo

        !Compute Hnorm and Htan magnitudes
        Hn_sq = 0.5D0*(Hn_sq + Hn(1)*Hn(1)+Hn(2)*Hn(2)+Hn(3)*Hn(3))
        Ht_sq = 0.5D0*(Ht_sq + Ht(1)*Ht(1)+Ht(2)*Ht(2)+Ht(3)*Ht(3))

        mxyz = mxyz + mxyzp
        mxyz(1) = mxyz(1)/sqrt(mxyz(1)*mxyz(1) + mxyz(2)*mxyz(2) + mxyz(3)*mxyz(3))
        du(i,j,k) = du(i,j,k) - 0.5*mudiff*(Ht_sq + muratio*Hn_sq)*  &
                    mxyz(1)*(2.0/(rho(i+1,j,k)+rho(i,j,k)))
        
     endif
  enddo; enddo; enddo
  


  do k=ks,ke;  do j=js,jev; do i=is,ie
     if(abs(cvof(i,j+1,k)-cvof(i,j,k))>EPSC/1d1) then  
        ! there is a non-zero grad H (H Heaviside function) 

       !Compute normal vectors and contributions at (i,j,k)
        do k0=-1,1; do j0=-1,1; do i0=-1,1
            stencil3x3(i0,j0,k0) = cvof(i+i0,j+j0,k+k0)
        enddo;enddo;enddo
        call mycs(stencil3x3,mxyz)
        
        H(1) = - 0.5D0*(phimag(i+1,j,k) - phimag(i-1,j,k))/dxh(i)
        H(2) = - 0.5D0*(phimag(i,j+1,k) - phimag(i,j-1,k))/dyh(j)
        H(3) = - 0.5D0*(phimag(i,j,k+1) - phimag(i,j,k-1))/dzh(k)

        !Compute H-normal and H-tangential components
        do i1=1,3
            Hn(i1) = (H(1)*mxyz(1)+H(2)*mxyz(2)+H(3)*mxyz(3))*mxyz(i1)           
            Ht(i1) = (H(i1)-Hn(i1))
        enddo

        !Compute Hnorm and Htan magnitudes
        Hn_sq = Hn(1)*Hn(1)+Hn(2)*Hn(2)+Hn(3)*Hn(3)
        Ht_sq = Ht(1)*Ht(1)+Ht(2)*Ht(2)+Ht(3)*Ht(3)

        !Compute normal vectors and contributions at (i,j+1,k)
        do k0=-1,1; do j0=-1,1; do i0=-1,1
            stencil3x3(i0,j0,k0) = cvof(i+i0,j+1+j0,k+k0)
        enddo;enddo;enddo
        call mycs(stencil3x3,mxyzp)
        
        H(1) = - 0.5D0*(phimag(i+1,j+1,k) - phimag(i-1,j+1,k))/dxh(i)
        H(2) = - 0.5D0*(phimag(i,j+2,k) - phimag(i,j,k))/dyh(j+1)
        H(3) = - 0.5D0*(phimag(i,j+1,k+1) - phimag(i,j+1,k-1))/dzh(k)

        !Compute H-normal and H-tangential components
        do i1=1,3
            Hn(i1) = (H(1)*mxyzp(1)+H(2)*mxyzp(2)+H(3)*mxyzp(3))*mxyzp(i1)           
            Ht(i1) = (H(i1)-Hn(i1))
        enddo

        !Compute Hnorm and Htan magnitudes
        Hn_sq = 0.5D0*(Hn_sq + Hn(1)*Hn(1)+Hn(2)*Hn(2)+Hn(3)*Hn(3))
        Ht_sq = 0.5D0*(Ht_sq + Ht(1)*Ht(1)+Ht(2)*Ht(2)+Ht(3)*Ht(3))

        mxyz = mxyz + mxyzp
        mxyz(2) = mxyz(2)/sqrt(mxyz(1)*mxyz(1) + mxyz(2)*mxyz(2) + mxyz(3)*mxyz(3))

        dv(i,j,k)=dv(i,j,k) -  0.5*mudiff*(Ht_sq + muratio*Hn_sq)*  &
                    mxyz(2)*(2.0/(rho(i,j+1,k)+rho(i,j,k)))

     endif
  enddo; enddo; enddo




  do k=ks,kew;  do j=js,je; do i=is,ie
     if(abs(cvof(i,j,k+1) - cvof(i,j,k))>EPSC/1d1) then  
        ! there is a non-zero grad H (H Heaviside function) 

       !Compute normal vectors and contributions at (i,j,k)
        do k0=-1,1; do j0=-1,1; do i0=-1,1
            stencil3x3(i0,j0,k0) = cvof(i+i0,j+j0,k+k0)
        enddo;enddo;enddo
        call mycs(stencil3x3,mxyz)
        
        H(1) = - 0.5D0*(phimag(i+1,j,k) - phimag(i-1,j,k))/dxh(i)
        H(2) = - 0.5D0*(phimag(i,j+1,k) - phimag(i,j-1,k))/dyh(j)
        H(3) = - 0.5D0*(phimag(i,j,k+1) - phimag(i,j,k-1))/dzh(k)

        !Compute H-normal and H-tangential components
        do i1=1,3
            Hn(i1) = (H(1)*mxyz(1)+H(2)*mxyz(2)+H(3)*mxyz(3))*mxyz(i1)           
            Ht(i1) = (H(i1)-Hn(i1))
        enddo

        !Compute Hnorm and Htan magnitudes
        Hn_sq = Hn(1)*Hn(1)+Hn(2)*Hn(2)+Hn(3)*Hn(3)
        Ht_sq = Ht(1)*Ht(1)+Ht(2)*Ht(2)+Ht(3)*Ht(3)

        !Compute normal vectors and contributions at (i,j,k+1)
        do k0=-1,1; do j0=-1,1; do i0=-1,1
            stencil3x3(i0,j0,k0) = cvof(i+i0,j+j0,k+1+k0)
        enddo;enddo;enddo
        call mycs(stencil3x3,mxyzp)
        
        H(1) = - 0.5D0*(phimag(i+1,j,k+1) - phimag(i-1,j,k+1))/dxh(i)
        H(2) = - 0.5D0*(phimag(i,j+1,k+1) - phimag(i,j-1,k+1))/dyh(j)
        H(3) = - 0.5D0*(phimag(i,j,k+2) - phimag(i,j,k))/dzh(k+1)

        !Compute H-normal and H-tangential components
        do i1=1,3
            Hn(i1) = (H(1)*mxyzp(1)+H(2)*mxyzp(2)+H(3)*mxyzp(3))*mxyzp(i1)           
            Ht(i1) = (H(i1)-Hn(i1))
        enddo

        !Compute Hnorm and Htan magnitudes
        Hn_sq = 0.5D0*(Hn_sq + Hn(1)*Hn(1)+Hn(2)*Hn(2)+Hn(3)*Hn(3))
        Ht_sq = 0.5D0*(Ht_sq + Ht(1)*Ht(1)+Ht(2)*Ht(2)+Ht(3)*Ht(3))

        mxyz = mxyz + mxyzp
        mxyz(3) = mxyz(3)/sqrt(mxyz(1)*mxyz(1) + mxyz(2)*mxyz(2) + mxyz(3)*mxyz(3))

        dw(i,j,k)=dw(i,j,k) -  0.5*mudiff*(Ht_sq + muratio*Hn_sq)*  &
                    mxyz(3)*(2.0/(rho(i,j,k+1)+rho(i,j,k)))

     endif
  enddo; enddo; enddo

end subroutine surfaceMagForce


!=================================================================================================
!=================================================================================================
! The Poisson equation for the magnetic potential, phi, is setup with matrix A as
! A7*Pijk = A1*Pi-1jk + A2*Pi+1jk + A3*Pij-1k + 
!           A4*Pij+1k + A5*Pijk-1 + A6*Pijk+1 + A8
!-------------------------------------------------------------------------------------------------
subroutine SetupMagPoisson(chip1,A,VolumeSource) 
  use module_grid
  use module_BC
  use module_2phase

  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: chip1
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
  real(8), intent(in) :: VolumeSource
  integer :: i,j,k
  
  do k=ks,ke; do j=js,je; do i=is,ie
    !    if(mask(i,j,k))then
    A(i,j,k,1) = 0.5d0*(chip1(i-1,j,k)+chip1(i,j,k))/(dx(i)*dxh(i-1))
    A(i,j,k,2) = 0.5d0*(chip1(i+1,j,k)+chip1(i,j,k))/(dx(i)*dxh(i  ))
    A(i,j,k,3) = 0.5d0*(chip1(i,j-1,k)+chip1(i,j,k))/(dy(j)*dyh(j-1))
    A(i,j,k,4) = 0.5d0*(chip1(i,j+1,k)+chip1(i,j,k))/(dy(j)*dyh(j  ))
    A(i,j,k,5) = 0.5d0*(chip1(i,j,k-1)+chip1(i,j,k))/(dz(k)*dzh(k-1))
    A(i,j,k,6) = 0.5d0*(chip1(i,j,k+1)+chip1(i,j,k))/(dz(k)*dzh(k  ))
    A(i,j,k,7) = sum(A(i,j,k,1:6))
    A(i,j,k,8) = 0.d0 ! no magnetic monopole sources
    !    endif
  enddo; enddo; enddo

!  call MagPoisson_BCs (A)
!  call check_and_debug_Poisson(A,umask,vmask,wmask,pmask,pmask,dt,VolumeSource)

end subroutine SetupMagPoisson





end module module_magnetic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ADD TO paris.f90 main program
!
 ! right after InitConditions call
!  call linfunc(chip1,1.0d0, 1.0d0+muratio,HarmMean)
!
 ! right after SurfaceForce
  ! call SetupMagPoisson
  ! call NewSolver_std
  ! call SurfaceMagForce
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ADD TO paris.f90 in subroutine initialize
!
 ! new variables mumag1,2 and default value=mumag0
 ! new flag test_MagneticDrop
 ! allocate phimag, chip1
 ! phimag, chip1 = 0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ADD TO paris.f90 in subroutine InitConditions
!
!  if (test_MagneticLinearFieldDroplet) ....
!  loop: phimag Initial Conditions: phimag = Az
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

