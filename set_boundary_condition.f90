! umat -> exp(i*P_umat*dtau_umat)*umat
! P_umat -> P_umat - delh_umat*dtau_umat
subroutine set_bc(umat)

  implicit none

  include 'size.h'
  
  double complex umat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1,1:ndim)

  integer imat,jmat,it,ix,iy,idim

  do idim=1,ndim
     do imat=1,nmat
        do jmat=1,nmat
           !do it=1,nt
           do ix=1,nx
              do iy=1,ny
                 umat(imat,jmat,0,ix,iy,idim)=umat(imat,jmat,nt,ix,iy,idim)
                 umat(imat,jmat,nt+1,ix,iy,idim)=umat(imat,jmat,1,ix,iy,idim)
              end do
           end do
           
           do it=1,nt
              !do ix=1,nx
              do iy=1,ny
                 umat(imat,jmat,it,0,iy,idim)=umat(imat,jmat,it,nx,iy,idim)
                 umat(imat,jmat,it,nx+1,iy,idim)=umat(imat,jmat,it,1,iy,idim)
              end do
           end do
           
           do it=1,nt
              do ix=1,nx
                 !do iy=1,ny
                 umat(imat,jmat,it,ix,0,idim)=umat(imat,jmat,it,ix,ny,idim)
                 umat(imat,jmat,it,ix,ny+1,idim)=umat(imat,jmat,it,ix,1,idim)
              end do
           end do

           do it=1,nt
              umat(imat,jmat,it,0,0,idim)=umat(imat,jmat,it,nx,ny,idim)
              umat(imat,jmat,it,0,ny+1,idim)=umat(imat,jmat,it,nx,1,idim)
              umat(imat,jmat,it,nx+1,0,idim)=umat(imat,jmat,it,1,ny,idim)
              umat(imat,jmat,it,nx+1,ny+1,idim)=umat(imat,jmat,it,1,1,idim)
           end do
     
           do ix=1,nx
              umat(imat,jmat,0,ix,0,idim)=umat(imat,jmat,nt,ix,ny,idim)
              umat(imat,jmat,nt+1,ix,0,idim)=umat(imat,jmat,1,ix,ny,idim)
              umat(imat,jmat,0,ix,ny+1,idim)=umat(imat,jmat,nt,ix,1,idim)
              umat(imat,jmat,nt+1,ix,ny+1,idim)=umat(imat,jmat,1,ix,1,idim)
           end do

           do iy=1,ny
              umat(imat,jmat,0,0,iy,idim)=umat(imat,jmat,nt,nx,iy,idim)
              umat(imat,jmat,nt+1,0,iy,idim)=umat(imat,jmat,1,nx,iy,idim)
              umat(imat,jmat,0,nx+1,iy,idim)=umat(imat,jmat,nt,1,iy,idim)
              umat(imat,jmat,nt+1,nx+1,iy,idim)=umat(imat,jmat,1,1,iy,idim)
           end do

           umat(imat,jmat,0,0,0,idim)=umat(imat,jmat,nt,nx,ny,idim)
           umat(imat,jmat,nt+1,0,0,idim)=umat(imat,jmat,1,nx,ny,idim)
           umat(imat,jmat,0,nx+1,0,idim)=umat(imat,jmat,nt,1,ny,idim)
           umat(imat,jmat,0,0,ny+1,idim)=umat(imat,jmat,nt,nx,1,idim)
           umat(imat,jmat,nt+1,nx+1,0,idim)=umat(imat,jmat,1,1,ny,idim)
           umat(imat,jmat,0,nx+1,ny+1,idim)=umat(imat,jmat,nt,1,1,idim)
           umat(imat,jmat,nt+1,0,ny+1,idim)=umat(imat,jmat,1,nx,1,idim)
           umat(imat,jmat,nt+1,nx+1,ny+1,idim)=umat(imat,jmat,1,1,1,idim)
           
        end do
     end do
  end do
  
  return

END subroutine Set_Bc
