!##############################################################################
!######              3d YM, Wilson's plaquette action, HMC            #########
!######                                                               #########
!######                 written by Masanori Hanada                    #########
!######                                                               #########
!####  The definition of g^2 is factor 2 different from the orbifold lattice ##
!##############################################################################
!Mersenne twister.
include 'mt19937.f90'
program D3YM
  
  use mtmod !Mersenne twistor
  implicit none
  include 'size.h'
  include 'include.h'
  
  double precision pol_phase(1:nmat,1:nx,1:ny),Pol_re,Pol_im
  
  open(unit=10,status='OLD',file='input_v1.dat',action='READ')
  read(10,*) input_config
  read(10,*) output_config
  read(10,*) data_output
  read(10,*) output_pol
  read(10,*) init
  read(10,*) at
  read(10,*) as
  read(10,*) ntraj
  read(10,*) nskip
  read(10,*) ntau
  read(10,*) dtau_t
  read(10,*) dtau_s
  read(10,*) mersenne_seed

  close(10)
  !*************************************
  !*** Set the initial configuration ***
  !*************************************

  if(init.EQ.0)then
     !continue from old config
     open(unit=9,status='OLD',file=input_config,action='READ')
     call mtgetu(9)
     read(9,*) itraj
     read(9,*) umat
     close(9)

  else if(init.EQ.1)then
     !initialize random number generator
     call sgrnd(mersenne_seed)
     !new config, cold start
     itraj=1
     umat=(0d0,0d0)
     do idim=1,ndim
        do it=0,nt+1
           do ix=0,nx+1
              do iy=0,ny+1
                 do imat=1,nmat
                    umat(imat,imat,it,ix,iy,idim)=(1d0,0d0)
                 end do
              end do
           end do
        end do
     end do
  end if
  !**************************************************
  !**************************************************
  nacceptance=0 !number of acceptance
  ntrial=0 !number of trial
  
  !************************************
  !************************************
  !     Make the output file
  !************************************
  !************************************
  
  open(unit=10,status='REPLACE',file=data_output,action='WRITE')
  open(unit=20,status='REPLACE',file=output_pol,action='WRITE')
  write(10,*) "#size of the gauge group: nmat=",nmat
  write(10,*) "#ntau=",ntau
  write(10,*) "#dtau for U_t=",Dtau_t
  write(10,*) "#dtau for U_x,U_y=",Dtau_s
  write(10,*) "# traj, ham_fin - ham_init, spatial_plaquette, Re(Pol.), Im(Pol.), acceptance"
  write(10,*) "#------------------------------------------------"
  nacceptance=0
  ntrial=0
  do while (itraj.LE.ntraj)
     
     backup_umat=umat
     call Generate_P_umat(P_umat)
     call Calc_Ham(at,as,umat,P_umat,ham_init)
     call Molecular_Dynamics(at,as,ntau,dtau_t,dtau_s,umat,P_umat)
     call Calc_Ham(at,as,umat,P_umat,ham_fin)
     
     metropolis=grnd()
     ntrial=ntrial+1
     If(dexp(ham_init-ham_fin) > metropolis)THEN
        !accept
        nacceptance=nacceptance+1
     else
        !reject
        umat=backup_umat
     end If
     
     ! write(*,*)itraj,-ham_init+ham_fin,dble(nacceptance)/dble(ntrial)
     
     ! measurements
     if(MOD(itraj,nskip).EQ.0)then
        
        !call Calc_action(at,as,umat,action)
        
        call Calc_spatial_plaquette(umat,spatial_plaquette)
        call Calc_Polyakov(umat,Pol_re,Pol_im,Pol_phase)

        write(10,*)itraj,-ham_init+ham_fin,spatial_plaquette,&
             &Pol_re,Pol_im,dble(nacceptance)/dble(ntrial)
        write(*,*)itraj,-ham_init+ham_fin,spatial_plaquette,&
             &Pol_re,Pol_im,dble(nacceptance)/dble(ntrial)
        do ix=1,nx
           do iy=1,ny
              !do imat=1,nmat
                 !write(20,*)ix,iy,imat,pol_phase(imat,ix,iy)
                 write(20,*)pol_phase(1,ix,iy),pol_phase(2,ix,iy)
              !end do
           end do
        end do
     end if

     itraj=itraj+1
  end do
  !**************************************************
  !**************************************************
  !   End of iteration
  !**************************************************
  !**************************************************

  close(10)
  close(20)



  open(UNIT = 22, File = output_config, STATUS = "REPLACE", ACTION = "WRITE")
  call mtsaveu(22)
  write(22,*) itraj
  write(22,*) umat
  close(22)






end program D3YM

include 'BoxMuller.f90'
include 'MATRIX_iEXP.f90'
include 'Generate_P_umat.f90'
include 'Calc_Ham.f90'
include 'Calc_action.f90'
include 'Molecular_Dynamics.f90'
include 'set_boundary_condition.f90'
include 'Calc_DELH.f90'
include 'Calc_Polyakov.f90'
include 'Calc_plaquette.f90'
