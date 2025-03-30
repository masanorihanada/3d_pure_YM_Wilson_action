!***********************************************************
!***********************************************************
!    Polaakov loop 
SUBROUTINE Calc_Polyakov(umat,Pol_re,Pol_im,pol_phase)
  
  implicit none
  include 'size.h'
  
  double complex umat(1:nmat,1:nmat,0:nt+1,0:nx+1,0:ny+1,1:ndim)
  double complex trace,MAT1(1:nmat,1:nmat),MAT2(1:nmat,1:nmat),eig_P(1:nmat)
  double precision Pol_re,Pol_im,pol_phase(1:nmat,1:nx,1:ny),pi
  double complex pol_pre
  integer imat,jmat,kmat
  integer isite,ix,iy,it

  pi=2d0*dasin(1d0)

  pol_pre=(0d0,0d0)
  do ix=1,nx
     do iy=1,ny
        do imat=1,nmat
           do jmat=1,nmat
              MAT1(imat,jmat)=umat(imat,jmat,1,ix,iy,1)
           end do
        end do
        do it=2,nt
           MAT2=(0d0,0d0)
           do imat=1,nmat
              do jmat=1,nmat
                 do kmat=1,nmat
                    MAT2(imat,jmat)=MAT2(imat,jmat)&
                         +MAT1(imat,kmat)*umat(kmat,jmat,it,ix,iy,1)
                 end do
              end do
           end do
           MAT1=MAT2
        end do
        trace=(0d0,0d0)
        do imat=1,nmat
           trace=trace+MAT2(imat,imat)
        end do
        Pol_pre=Pol_pre+trace

        call Diagonalization(nmat,MAT2,eig_P)
        do imat=1,nmat
           if(dble(eig_P(imat)).GE.0d0)then
              pol_phase(imat,ix,iy)=dasin(dble((0d0,-1d0)*eig_P(imat)))
             
           else if(dble((0d0,-1d0)*eig_P(imat)).GT.0d0)then
              pol_phase(imat,ix,iy)=pi-dasin(dble((0d0,-1d0)*eig_P(imat)))
           else
              pol_phase(imat,ix,iy)=-pi-dasin(dble((0d0,-1d0)*eig_P(imat)))
           end if
           
        end do
     end do
  end do
  Pol_re=dble(Pol_pre)/dble(NMAT*nx*ny)
  Pol_im=dble(Pol_pre*(0d0,-1d0))/dble(NMAT*nx*ny)

  return
  
END SUBROUTINE Calc_Polyakov

!#####################################################
SUBROUTINE Diagonalization(MatSize,MAT,w)

  implicit none

  integer MatSize
  double complex MAT(1:MatSize,1:MatSize)
  integer i
  !lapack
  character jobvl,jobvr
  integer lda,ldvl,ldvr,info,lwork

  double complex w(1:MatSize),VL(1:MatSize),VR(1:MatSize,1:MATSIZE),&
       work(1:2*MatSize),rwork(1:2*MatSize)


  jobvl='N'
  jobvr='V'
  lda=MatSize
  ldvl=1
  ldvr=MatSize
  lwork=2*MatSize

  !MAT2=MAT

  call ZGEEV(jobvl,jobvr,MatSize,MAT,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)


  return

END SUBROUTINE Diagonalization
