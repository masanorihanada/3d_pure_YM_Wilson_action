! Generate P_umat with Gaussian weight.
SUBROUTINE Generate_P_umat(P_umat)
  
  implicit none
  include 'size.h'
  
  integer it,ix,iy,idim,imat,jmat
  double precision r1,r2,trace

  double complex P_umat(1:nmat,1:nmat,1:nt,1:nx,1:ny,1:ndim)

  do idim=1,ndim
     do it=1,nt
        do ix=1,nx
           do iy=1,ny
              do imat=1,nmat-1
                 do jmat=imat+1,nmat
                    call BoxMuller(r1,r2)
                    P_umat(imat,jmat,it,ix,iy,idim)=&
                         dcmplx(r1/dsqrt(2d0))+dcmplx(r2/dsqrt(2d0))*(0D0,1D0)
                    P_umat(jmat,imat,it,ix,iy,idim)=&
                         dcmplx(r1/dsqrt(2d0))-dcmplx(r2/dsqrt(2d0))*(0D0,1D0)
                 end do
              end do
              trace=0d0
              do imat=1,nmat
                 call BoxMuller(r1,r2)
                 P_umat(imat,imat,it,ix,iy,idim)=dcmplx(r1)
                 trace=trace+r1
              end do
              !traceless condition
              do imat=1,nmat
                 call BoxMuller(r1,r2)
                 P_umat(imat,imat,it,ix,iy,idim)=P_umat(imat,imat,it,ix,iy,idim)-dcmplx(trace/dble(nmat))
              end do
           end do
        end do
     end do
  end do
  
  
     
  return
  
END SUBROUTINE Generate_P_umat
