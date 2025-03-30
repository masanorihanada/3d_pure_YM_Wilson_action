! umat -> exp(i*P_umat*dtau)*umat
! P_umat -> P_umat - delh_umat*dtau

subroutine Molecular_Dynamics(at,as,ntau,dtau_t,dtau_s,umat,P_umat)

  implicit none

  include 'size.h'

  integer ntau
  double complex umat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1,1:ndim)
  
  double precision at,as
  double precision dtau_t,dtau_s,dtau(1:ndim)
  integer ix,iy,it
  integer idim,jdim
  integer imat,jmat,kmat
  integer itau
  
  double complex P_umat(1:nmat,1:nmat,1:nt,1:nx,1:ny,1:ndim)
  double complex delh_umat(1:nmat,1:nmat,1:nt,1:nx,1:ny,1:ndim)
  double complex MAT(1:NMAT,1:NMAT),EXP_iP(1:NMAT,1:NMAT)

  dtau(1)=dtau_t
  dtau(2)=dtau_s
  dtau(3)=dtau_s
  
  ! first step of leap frog
  do idim=1,ndim
     do it=1,nt
        do ix=1,nx
           do iy=1,ny
              do imat=1,nmat
                 do jmat=1,nmat
                    MAT(imat,jmat)=P_umat(imat,jmat,it,ix,iy,idim)*dcmplx(0.5d0*dtau(idim))
                 end do
              end do
              call MATRIX_iEXP(NMAT,MAT,EXP_iP)
              do imat=1,nmat
                 do jmat=1,nmat
                    MAT(imat,jmat)=umat(imat,jmat,it,ix,iy,idim)
                    umat(imat,jmat,it,ix,iy,idim)=(0d0,0d0)
                 end do
              end do
              do imat=1,nmat
                 do jmat=1,nmat
                    do kmat=1,nmat
                       umat(imat,jmat,it,ix,iy,idim)=umat(imat,jmat,it,ix,iy,idim)&
                            +EXP_iP(imat,kmat)*MAT(kmat,jmat) 
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do
  call set_bc(umat)
  ! second,...,Ntau-th   
  do itau=2,ntau
     call Calc_DELH(at,as,umat,delh_umat)
     do idim=1,ndim
        do it=1,nt
           do ix=1,nx
              do iy=1,ny
                 do imat=1,nmat
                    do jmat=1,nmat
                       P_umat(imat,jmat,it,ix,iy,idim)&
                            &=P_umat(imat,jmat,it,ix,iy,idim)&
                            &-delh_umat(imat,jmat,it,ix,iy,idim)*dcmplx(dtau(idim))
                    end do
                 end do
              end do
           end do
        end do
     end do
     
     do idim=1,ndim
        do it=1,nt
           do ix=1,nx
              do iy=1,ny
                 do imat=1,nmat
                    do jmat=1,nmat
                       MAT(imat,jmat)=P_umat(imat,jmat,it,ix,iy,idim)*dcmplx(dtau(idim))
                    end do
                 end do
                 call MATRIX_iEXP(NMAT,MAT,EXP_iP)
                 do imat=1,nmat
                    do jmat=1,nmat
                       MAT(imat,jmat)=umat(imat,jmat,it,ix,iy,idim)
                       umat(imat,jmat,it,ix,iy,idim)=(0d0,0d0)
                    end do
                 end do
                 do imat=1,nmat
                    do jmat=1,nmat
                       do kmat=1,nmat
                          umat(imat,jmat,it,ix,iy,idim)=umat(imat,jmat,it,ix,iy,idim)&
                               +EXP_iP(imat,kmat)*MAT(kmat,jmat) 
                       end do
                    end do
                 end do
                 
              end do
           end do
        end do
     end do
     call set_bc(umat)   
  end do
  
  ! last step
  call Calc_DELH(at,as,umat,delh_umat)
  do idim=1,ndim
     do it=1,nt
        do ix=1,nx
           do iy=1,ny
              do imat=1,nmat
                 do jmat=1,nmat
                    P_umat(imat,jmat,it,ix,iy,idim)&
                         &=P_umat(imat,jmat,it,ix,iy,idim)&
                         &-delh_umat(imat,jmat,it,ix,iy,idim)*dcmplx(dtau(idim))
                 end do
              end do
           end do
        end do
     end do
  end do
  
  do idim=1,ndim
     do it=1,nt
        do ix=1,nx
           do iy=1,ny
              do imat=1,nmat
                 do jmat=1,nmat
                    MAT(imat,jmat)=P_umat(imat,jmat,it,ix,iy,idim)*dcmplx(0.5d0*dtau(idim))
                 end do
              end do
              call MATRIX_iEXP(NMAT,MAT,EXP_iP)
              do imat=1,nmat
                 do jmat=1,nmat
                    MAT(imat,jmat)=umat(imat,jmat,it,ix,iy,idim)
                    umat(imat,jmat,it,ix,iy,idim)=(0d0,0d0)
                 end do
              end do
              do imat=1,nmat
                 do jmat=1,nmat
                    do kmat=1,nmat
                       umat(imat,jmat,it,ix,iy,idim)=umat(imat,jmat,it,ix,iy,idim)&
                            +EXP_iP(imat,kmat)*MAT(kmat,jmat) 
                    end do
                 end do
              end do
              
           end do
        end do
     end do
  end do
  call set_bc(umat) 
  
  return

END subroutine Molecular_Dynamics
