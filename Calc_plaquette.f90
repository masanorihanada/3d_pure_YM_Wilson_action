subroutine Calc_plaquette(umat,spatial_plaquette,temporal_plaquette)

  implicit none
  include 'size.h'
  
  double complex umat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1,1:ndim)
  
  double precision action,kin,pot
  double precision at,as

  integer ix,iy,it
  integer idim,jdim
  integer imat,jmat,kmat

  double complex uxuy(1:nmat,1:nmat)
  double complex uxduyd(1:nmat,1:nmat)

  double complex utux(1:nmat,1:nmat)
  double complex utduxd(1:nmat,1:nmat)
  
  double complex utuy(1:nmat,1:nmat)
  double complex utduyd(1:nmat,1:nmat)
  double precision spatial_plaquette,temporal_plaquette
  
  spatial_plaquette=0d0
  temporal_plaquette=0d0

  do it=1,nt
     do ix=1,nx
        do iy=1,ny
           uxuy=(0d0,0d0)
           uxduyd=(0d0,0d0)

           utux=(0d0,0d0)
           utduxd=(0d0,0d0)
           
           utuy=(0d0,0d0)
           utduyd=(0d0,0d0)
           do imat=1,nmat
              do jmat=1,nmat
                 do kmat=1,nmat
                    uxuy(imat,jmat)=uxuy(imat,jmat)&
                         &+umat(imat,kmat,it,ix,iy,2)*umat(kmat,jmat,it,ix+1,iy,3)
                    uxduyd(imat,jmat)=uxduyd(imat,jmat)&
                         &+dconjg(umat(kmat,imat,it,ix,iy+1,2))*dconjg(umat(jmat,kmat,it,ix,iy,3))


                    utux(imat,jmat)=utux(imat,jmat)&
                         &+umat(imat,kmat,it,ix,iy,1)*umat(kmat,jmat,it+1,ix,iy,2)
                    utduxd(imat,jmat)=utduxd(imat,jmat)&
                         &+dconjg(umat(kmat,imat,it,ix+1,iy,1))*dconjg(umat(jmat,kmat,it,ix,iy,2))

                    utuy(imat,jmat)=utuy(imat,jmat)&
                         &+umat(imat,kmat,it,ix,iy,1)*umat(kmat,jmat,it+1,ix,iy,3)
                    utduyd(imat,jmat)=utduyd(imat,jmat)&
                         &+dconjg(umat(kmat,imat,it,ix,iy+1,1))*dconjg(umat(jmat,kmat,it,ix,iy,3))
                    
                 end do
              end do
           end do
           
           do imat=1,nmat
              do jmat=1,nmat
                 spatial_plaquette = spatial_plaquette + dble(uxuy(imat,jmat)*uxduyd(jmat,imat))
                 temporal_plaquette = temporal_plaquette &
                      + dble(utux(imat,jmat)*utduxd(jmat,imat))&
                      + dble(utuy(imat,jmat)*utduyd(jmat,imat))
              end do
           end do
           
        end do
     end do
  end do

  spatial_plaquette = spatial_plaquette / dble(nx * ny * nt)
  temporal_plaquette = temporal_plaquette / dble(2 * nx * ny * nt)
  
  return

END subroutine Calc_plaquette
