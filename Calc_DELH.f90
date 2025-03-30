! umat -> exp(i*P_umat*dtau_umat)*umat
! xmat -> xmat + P_xmat*dtau_xmat
! P_umat -> P_umat - delh_umat*dtau_umat
! P_xmat -> P_xmat - delh_xmat*dtau_xmat
! delh_xmat(imat,jmat)=dS/dxmat(jmat,imat)

SUBROUTINE Calc_DELH(at,as,umat,delh_umat)

  implicit none
  include 'size.h'

  double complex umat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1,1:ndim)
  double complex delh_umat(1:nmat,1:nmat,1:nt,1:nx,1:ny,1:ndim)
  double precision at,as

  integer ix,iy,it
  integer idim,jdim
  integer imat,jmat,kmat

  double complex uu(1:nmat,1:nmat,1:ndim,1:ndim)
  double complex uud(1:nmat,1:nmat,1:ndim,1:ndim)
  double complex udu(1:nmat,1:nmat,1:ndim,1:ndim)
  double complex uuudud(1:nmat,1:nmat,1:ndim,1:ndim)
  double complex uududu(1:nmat,1:nmat,1:ndim,1:ndim)

  delh_umat=(0d0,0d0)
  
  do it=1,nt
     do ix=1,nx
        do iy=1,ny
           uu=(0d0,0d0)
           uud=(0d0,0d0)
           udu=(0d0,0d0)
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
                    !********
                    uud(imat,jmat,1,2)=uud(imat,jmat,1,2)&
                         &+umat(imat,kmat,it,ix,iy,1)*dconjg(umat(jmat,kmat,it+1,ix-1,iy,2))
                    uud(imat,jmat,1,3)=uud(imat,jmat,1,3)&
                         &+umat(imat,kmat,it,ix,iy,1)*dconjg(umat(jmat,kmat,it+1,ix,iy-1,3))

                    uud(imat,jmat,2,1)=uud(imat,jmat,2,1)&
                         &+umat(imat,kmat,it,ix,iy,2)*dconjg(umat(jmat,kmat,it-1,ix+1,iy,1))
                    uud(imat,jmat,2,3)=uud(imat,jmat,2,3)&
                         &+umat(imat,kmat,it,ix,iy,2)*dconjg(umat(jmat,kmat,it,ix+1,iy-1,3))

                    uud(imat,jmat,3,1)=uud(imat,jmat,3,1)&
                         &+umat(imat,kmat,it,ix,iy,3)*dconjg(umat(jmat,kmat,it-1,ix,iy+1,1))
                    uud(imat,jmat,3,2)=uud(imat,jmat,3,2)&
                         &+umat(imat,kmat,it,ix,iy,3)*dconjg(umat(jmat,kmat,it,ix-1,iy+1,2))
                    !********
                    udu(imat,jmat,1,2)=udu(imat,jmat,1,2)&
                         &+dconjg(umat(kmat,imat,it,ix-1,iy,1))*umat(kmat,jmat,it,ix-1,iy,2)
                    udu(imat,jmat,1,3)=udu(imat,jmat,1,3)&
                         &+dconjg(umat(kmat,imat,it,ix,iy-1,1))*umat(kmat,jmat,it,ix,iy-1,3)

                    udu(imat,jmat,2,1)=udu(imat,jmat,2,1)&
                         &+dconjg(umat(kmat,imat,it-1,ix,iy,2))*umat(kmat,jmat,it-1,ix,iy,1)
                    udu(imat,jmat,2,3)=udu(imat,jmat,2,3)&
                         &+dconjg(umat(kmat,imat,it,ix,iy-1,2))*umat(kmat,jmat,it,ix,iy-1,3)

                    udu(imat,jmat,3,1)=udu(imat,jmat,3,1)&
                         &+dconjg(umat(kmat,imat,it-1,ix,iy,3))*umat(kmat,jmat,it-1,ix,iy,1)
                    udu(imat,jmat,3,2)=udu(imat,jmat,3,2)&
                         &+dconjg(umat(kmat,imat,it,ix-1,iy,3))*umat(kmat,jmat,it,ix-1,iy,2)
                    
                 end do
              end do
           end do
           
           uuudud=(0d0,0d0)
           uududu=(0d0,0d0)
           do idim=1,ndim
              do jdim=1,ndim
                 if(idim.NE.jdim)then
                    do imat=1,nmat
                       do jmat=1,nmat
                          do kmat=1,nmat
                             uuudud(imat,jmat,idim,jdim)=uuudud(imat,jmat,idim,jdim)&
                                  &+uu(imat,kmat,idim,jdim)*dconjg(uu(jmat,kmat,jdim,idim))
                             
                             uududu(imat,jmat,idim,jdim)=uududu(imat,jmat,idim,jdim)&
                                  &+uud(imat,kmat,idim,jdim)*udu(kmat,jmat,idim,jdim)
                          end do
                       end do
                    end do
                 end if
              end do
           end do

           do imat=1,nmat
              do jmat=1,nmat
                 delh_umat(imat,jmat,it,ix,iy,1)=&
                      delh_umat(imat,jmat,it,ix,iy,1)&
                      &+(uuudud(imat,jmat,1,2)-dconjg(uuudud(jmat,imat,1,2))&
                      &+uuudud(imat,jmat,1,3)-dconjg(uuudud(jmat,imat,1,3))&
                      &+uududu(imat,jmat,1,2)-dconjg(uududu(jmat,imat,1,2))&
                      &+uududu(imat,jmat,1,3)-dconjg(uududu(jmat,imat,1,3))&
                      )/at*dcmplx(nmat)*(0d0,-1d0)

                 delh_umat(imat,jmat,it,ix,iy,2)=&
                      delh_umat(imat,jmat,it,ix,iy,2)&
                      &+(uuudud(imat,jmat,2,1)-dconjg(uuudud(jmat,imat,2,1))&
                      &+uududu(imat,jmat,2,1)-dconjg(uududu(jmat,imat,2,1)))&
                      &/at*dcmplx(nmat)*(0d0,-1d0)&
                      &+(uuudud(imat,jmat,2,3)-dconjg(uuudud(jmat,imat,2,3))&
                      &+uududu(imat,jmat,2,3)-dconjg(uududu(jmat,imat,2,3)))&
                      &*at/(as*as)*dcmplx(nmat)*(0d0,-1d0)

                 delh_umat(imat,jmat,it,ix,iy,3)=&
                      delh_umat(imat,jmat,it,ix,iy,3)&
                      &+(uuudud(imat,jmat,3,1)-dconjg(uuudud(jmat,imat,3,1))&
                      &+uududu(imat,jmat,3,1)-dconjg(uududu(jmat,imat,3,1)))&
                      &/at*dcmplx(nmat)*(0d0,-1d0)&
                      &+(uuudud(imat,jmat,3,2)-dconjg(uuudud(jmat,imat,3,2))&
                      &+uududu(imat,jmat,3,2)-dconjg(uududu(jmat,imat,3,2)))&
                      &*at/(as*as)*dcmplx(nmat)*(0d0,-1d0)

                 
              end do
           end do
           
        end do
     end do
  end do
  
  return

END SUBROUTINE Calc_DELH
