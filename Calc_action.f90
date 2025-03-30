subroutine Calc_action(at,as,umat,action)

  implicit none
  include 'size.h'
  
  double complex umat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1,1:ndim)
  
  double precision action,kin,pot
  double precision at,as

  integer ix,iy,it
  integer idim,jdim
  integer imat,jmat,kmat

  double complex uu(1:nmat,1:nmat,1:ndim,1:ndim)
  
  kin=0d0
  pot=0d0

  do it=1,nt
     do ix=1,nx
        do iy=1,ny
           uu=(0d0,0d0)
           do imat=1,nmat
              do jmat=1,nmat
                 do kmat=1,nmat
                    uu(imat,jmat,1,2)=uu(imat,jmat,1,2)&
                         &+umat(imat,kmat,it,ix,iy,1)*umat(kmat,jmat,it+1,ix,iy,2)
                    uu(imat,jmat,1,3)=uu(imat,jmat,1,3)&
                         &+umat(imat,kmat,it,ix,iy,1)*umat(kmat,jmat,it+1,ix,iy,3)
                    
                    uu(imat,jmat,2,1)=uu(imat,jmat,2,1)&
                         &+umat(imat,kmat,it,ix,iy,2)*umat(kmat,jmat,it,ix+1,iy,1)
                    uu(imat,jmat,2,3)=uu(imat,jmat,2,3)&
                         &+umat(imat,kmat,it,ix,iy,2)*umat(kmat,jmat,it,ix+1,iy,3)
                    
                    uu(imat,jmat,3,1)=uu(imat,jmat,3,1)&
                         &+umat(imat,kmat,it,ix,iy,3)*umat(kmat,jmat,it,ix,iy+1,1)
                    uu(imat,jmat,3,2)=uu(imat,jmat,3,2)&
                         &+umat(imat,kmat,it,ix,iy,3)*umat(kmat,jmat,it,ix,iy+1,2)
                    
                 end do
              end do
           end do
           
           do imat=1,nmat
              do jmat=1,nmat
                 kin=kin+dble(uu(imat,jmat,1,2)*dconjg(uu(imat,jmat,2,1)))
                 kin=kin+dble(uu(imat,jmat,1,3)*dconjg(uu(imat,jmat,3,1)))
                 !  kin=kin+dble(uu(imat,jmat,2,1)*dconjg(uu(imat,jmat,1,2)))
                 !  kin=kin+dble(uu(imat,jmat,3,1)*dconjg(uu(imat,jmat,1,3)))
                 pot=pot+dble(uu(imat,jmat,2,3)*dconjg(uu(imat,jmat,3,2)))
                 !  pot=pot+dble(uu(imat,jmat,3,2)*dconjg(uu(imat,jmat,2,3)))
              end do
           end do
           
        end do
     end do
  end do

  kin=kin/at*dble(nmat)*(-2d0)
  pot=pot*at/(as*as)*dble(nmat)*(-2d0)
  
  action=kin+pot
  
  return

END subroutine Calc_action
