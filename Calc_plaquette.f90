subroutine Calc_spatial_plaquette(umat,spatial_plaquette)

  implicit none
  include 'size.h'
  
  double complex umat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1,1:ndim)
  
  double precision action,kin,pot
  double precision at,as

  integer ix,iy,it
  integer idim,jdim
  integer imat,jmat,kmat

  double complex uu(1:nmat,1:nmat,1:ndim)
  double precision spatial_plaquette
  
  spatial_plaquette=0d0

  do it=1,nt
     do ix=1,nx
        do iy=1,ny
           uu=(0d0,0d0)
           do imat=1,nmat
              do jmat=1,nmat
                 do kmat=1,nmat
                    uu(imat,jmat)=uu(imat,jmat)&
                         &+umat(imat,kmat,it,ix,iy,1)*umat(kmat,jmat,it+1,ix,iy,2)
                    uu(imat,jmat,1,3)=uu(imat,jmat,1,3)&
                         &+umat(imat,kmat,it,ix,iy,1)*umat(kmat,jmat,it+1,ix,iy,3)
                    
                 end do
              end do
           end do
           
           do imat=1,nmat
              do jmat=1,nmat
                 spatial_plaquette = spatial_plaquette + dble(uu(imat,jmat)*dconjg(uu(imat,jmat)))
              end do
           end do
           
        end do
     end do
  end do

spatial_plaquette = spatial_plaquette / dble(nx * ny * nt)
  
  return

END subroutine Calc_spatial_plaquette
