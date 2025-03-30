subroutine Calc_Ham(at,as,umat,P_umat,ham)

  implicit none
  include 'size.h'
  
  double complex umat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1,1:ndim)
  double complex P_umat(1:nmat,1:nmat,1:nt,1:nx,1:ny,1:ndim)
  double precision ham
  double precision at,as

  integer it,ix,iy
  integer idim
  integer imat,jmat

  call Calc_action(at,as,umat,ham)

  do idim=1,ndim
     do it=1,nt
        do ix=1,nx
           do iy=1,ny
              do imat=1,nmat
                 do jmat=1,nmat
                    ham=ham&
                         +dble(P_umat(imat,jmat,it,ix,iy,idim)&
                         *P_umat(jmat,imat,it,ix,iy,idim))*0.5d0
                 end do
              end do
           end do
        end do
     end do
  end do

  return
  
END subroutine Calc_Ham
