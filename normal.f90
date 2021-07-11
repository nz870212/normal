! A fortran95 program for G95
! By NZ
module gdata
integer, parameter ::  nx=200,ny=200,np=5,Q=200,QQ=400,parent_QQ=15,cores=4
real(kind=8), dimension(nx, ny, np) :: lap_phi,ss,phi_value,phi_value_mid
real(kind=8), dimension(nx, ny) :: mic,mij,Num_p,sum_phi_value,r_c_i,r_c_i_real,N_sub_i
!real(kind=8), dimension(nx, ny) ::S_energy,R_grain,r_sub,r_sub_0,stress_i,nucleation_rate,mis_angle_show
!real(kind=8) :: r_sub_init(QQ+1)
integer, dimension(nx, ny, np) :: phi_name, phi_name_mid
!integer, dimension(nx_0, ny_0, np) :: phi_name_0,parent_phi_name
integer :: T=1000,Qb=21866,Qa=386573,step_mid=100,itimestep
real(kind=8) :: dx, dy, dt, h=1e-6, gb_energy=0.4, gb_energy_sub,gb_energy_real=0.4, gb_width, grain_size,dislocation_0!=1.9e13!1.1229e13
real(kind=8) :: pi=3.1415926, M, M0, M_LAG, R=8.314, dE, timestep_nucleation, ll,stress_y,stress_s,strain_c,stress_c
real(kind=8) :: def_rate, alpha=0.5, b_vector=2.58e-10, miu, K_sub=23, yita, M_1, hard_rate,diameter,stress_0
real(kind=8) :: sigema_0, omega_0
integer,dimension(QQ) :: S
integer :: num_Max_phi_name,timestep
end module gdata
!=============================================================================================================================
program main
use gdata
implicit none
integer :: nnn,MM0
real :: time_begin,time_end
call CPU_TIME(time_begin)
print *, 'start'
call init
do itimestep = 1,200000!itimestep
  call gradient   
  !print*, itimestep
  IF (mod(itimestep,step_mid)==0) THEN
    call microstructure
    call account_S
    CALL printdata(itimestep)
  print *, itimestep, "th iteration is running"
  ENDIF
end do

call CPU_TIME(time_end)
print *, time_end-time_begin
print *, 'stop'
!pause
stop
end program main
!=============================================================================================================================
subroutine init
use gdata
implicit none
integer ::i,j,k,l,mm
real::d,dmin,a,b,e,f,mid,Sst,o,p,G_energy,r_c_h
integer, dimension(1,QQ) :: xx,yy
    phi_name(:,:,1)=QQ+1
    phi_value(:,:,1)=1
    phi_name_mid(:,:,1)=QQ+1
    phi_value_mid(:,:,1)=1
    dx = h
    dy = h
    gb_width=6*h
    M0=1.1713e-7!(pi**2/(8.0*gb_width))*(M0/T)*exp(-Qb/(R*T))
    M_1=(M0)*exp(-Qb/(R*T))
    sigema_0 = 4.0*gb_energy/gb_width
    omega_0 = (2.0/pi)*sqrt(2*gb_energy*gb_width)
    M=(pi**2/(8.0*gb_width))*M_1
    !dt = 0.0075!dx**2/(4.0*omega_0**2*M)
    dt =dx**2/(4.0*((2.0/pi)*sqrt(2*gb_energy*gb_width))**2*(M))
    print*,dt
    dt=3e-5
    
    
    !cal_step=5000/dt
    do i=1,QQ
    call random_number(o)
    call random_number(p)
    xx(1,i)=ceiling(nx*o)
    yy(1,i)=ceiling(ny*p)
    end do
    
   
    do j=1,nx
    do k=1,ny
    a=1.0*j
    b=1.0*k
    dmin=sqrt((a)**2+(b)**2)
    l=1 
    do mm=1,QQ
    e=1.0*xx(1,mm)
    f=1.0*yy(1,mm)
    d=sqrt((e-a)**2+(f-b)**2)
    if (d<=dmin) then
    dmin=d
    l=mm
    phi_name(j,k,1)=l
    phi_value(j,k,1)=1.0
    phi_name_mid(j,k,1)=l
    phi_value_mid(j,k,1)=1.0
    endif
    enddo
    
    enddo
    enddo

	print *, 'intialized'
end subroutine init

!=============================================================================================================================
    
!主要计算步骤，计算每个时间步的相场该变量。
subroutine gradient
use gdata
use omp_lib
implicit none
integer :: i,j,k,ip,jp,im,jm,N,NN=0,ii,jj,kk,num_cal,tmp_BB
real :: tmp_temp_phivalue,sum_temp_phivalue,sum_dis_eng,anti_gb_energy
real :: phis, phin, phiw, phie, phic, phisw, phine, phise, phinw
real,external :: dfphi
integer,dimension(45) :: A=0,B=0
real,dimension(3,3,np) :: temp_value,temp_angle
real :: temp_phi,temp_phi_1,temp_phi_2,evol_phi
integer,dimension(3,3,np) :: temp_name
integer,dimension(:),allocatable :: BB
real,dimension(:),allocatable :: temp_lap
real,dimension(:,:,:),allocatable :: temp_phivalue,temp_anglevalue,sigma_3
real,dimension(:,:),allocatable :: sigma,omiga,a_sigma,a_omiga,delta_angle,anti_MM
real,dimension(5) :: value
integer,dimension(5) :: name

call omp_set_num_threads(cores)
!$OMP PARALLEL PRIVATE(A,B,NN,num_cal,temp_value,temp_name,BB,N,&
& temp_lap,temp_phivalue,temp_anglevalue,temp_angle,sigma,sigma_3,tmp_BB,&
& omiga,i,j,ip,jp,im,jm,kk,ii,jj,k,sum_dis_eng,tmp_temp_phivalue,sum_temp_phivalue,value,name,&
& a_sigma,a_omiga,delta_angle,anti_MM,anti_gb_energy,temp_phi,temp_phi_1,temp_phi_2,evol_phi)
do i=1,nx
    !$OMP DO 
    !schedule(dynamic)
do j=1,ny
A=0
B=0
!periodic boundary condition
ip = mod((i+1),NX)
im = mod((NX+i-1),NX)
jp = mod((j+1),NY)
jm = mod((NY+j-1),NY)
if (im==0) im=NX
if (jm==0) jm=NY
if (ip==0) ip=1
if (jp==0) jp=1

!do ii=1,3
!    do jj=1,3
!        do kk=1,5
!            temp_name(ii,jj,kk)=phi_name(i+(ii-2),j+(jj-2),kk)
!            temp_value(ii,jj,kk)=phi_value(i+(ii-2),j+(jj-2),kk)
!        end do
!    end do
!end do

do kk=1,5
    temp_name(1,1,kk)=phi_name(im,jm,kk)
    temp_value(1,1,kk)=phi_value(im,jm,kk)
end do
do kk=1,5
    temp_name(1,2,kk)=phi_name(im,j,kk)
    temp_value(1,2,kk)=phi_value(im,j,kk)
end do
do kk=1,5
    temp_name(1,3,kk)=phi_name(im,jp,kk)
    temp_value(1,3,kk)=phi_value(im,jp,kk)
end do
do kk=1,5
    temp_name(2,1,kk)=phi_name(i,jm,kk)
    temp_value(2,1,kk)=phi_value(i,jm,kk)
end do
do kk=1,5
    temp_name(2,2,kk)=phi_name(i,j,kk)
    temp_value(2,2,kk)=phi_value(i,j,kk)
end do
do kk=1,5
    temp_name(2,3,kk)=phi_name(i,jp,kk)
    temp_value(2,3,kk)=phi_value(i,jp,kk)
end do
do kk=1,5
    temp_name(3,1,kk)=phi_name(ip,jm,kk)
    temp_value(3,1,kk)=phi_value(ip,jm,kk)
end do
do kk=1,5
    temp_name(3,2,kk)=phi_name(ip,j,kk)
    temp_value(3,2,kk)=phi_value(ip,j,kk)
end do
do kk=1,5
    temp_name(3,3,kk)=phi_name(ip,jp,kk)
    temp_value(3,3,kk)=phi_value(ip,jp,kk)
end do
NN=0
do ii=1,3
    do jj=1,3
        do kk=1,np
            A(kk+3*(jj-1)+15*(ii-1))=temp_name(ii,jj,kk)
        end do
    end do
end do

N=1
B=0
B(N)=A(1)
III:  DO ii = 2,45
  DO jj = 1,N
    IF(B(jj)==A(ii) .or. A(ii)==0) CYCLE III
  ENDDO
  N=N+1
  B(N) = A(ii)
ENDDO III
do ii=1,45
    if (B(ii) .ne. 0) then
        NN=NN+1
    end if
end do
!print*,3
allocate(BB(NN))
!BB:临时储存相场名的数组
BB=B(1:NN)
!print*, BB
num_cal=size(BB)
allocate(temp_lap(NN))
allocate(temp_phivalue(3,3,num_cal))
allocate(sigma(NN,NN))
ALLOCATE(sigma_3(NN,NN,NN))
allocate(omiga(NN,NN))
allocate(anti_MM(NN,NN))
temp_phivalue=0.0
do k=1,num_cal
!
do ii=1,3
do jj=1,3
do kk=1,5
if (temp_name(ii,jj,kk) .eq. BB(k)) then
temp_phivalue(ii,jj,k)=temp_value(ii,jj,kk)
end if
end do
end do
end do
!print*, 'temp_anglevalue =',temp_anglevalue
do ii=1,NN!这里有问题,laplace的值不对。
phie=temp_phivalue(3,2,ii)
phiw=temp_phivalue(1,2,ii)
phis=temp_phivalue(2,3,ii)
phin=temp_phivalue(2,1,ii)
phic=temp_phivalue(2,2,ii)
phine=temp_phivalue(3,3,ii)
phinw=temp_phivalue(1,3,ii)
phise=temp_phivalue(3,1,ii)
phisw=temp_phivalue(1,1,ii)
temp_lap(ii) =(2*(phiw + phie + phis + phin -4*phic)+(phisw + phine &
&+ phise + phinw -4*phic))/(4*h*h)
end do
end do

do k=1,NN

	DO kk=1,NN
		DO ii=1,NN
			DO jj=1,NN
				IF (jj==ii) THEN
					sigma_3(ii,jj,kk)=0.0					
				ELSE
					sigma_3(ii,jj,kk)=sigema_0!omega_0
				END IF
			END DO
		END DO
	END DO
	
do ii=1,NN
do jj=1,NN
if (jj==ii) then
sigma(ii,jj)=0.0
omiga(ii,jj)=0.0
else
sigma(ii,jj)=sigema_0!omega_0
omiga(ii,jj)=omega_0!sigema_0
end if
end do
end do
anti_MM=M
temp_phi_1=0.0
temp_phi=0.0
do ii=1,NN
do jj=1,NN
	temp_phi_2=0.0
			DO kk=1,NN
				temp_phi_2=temp_phi_2+(sigma_3(k,jj,kk)-sigma_3(ii,jj,kk))*temp_phivalue(2,2,jj)*temp_phivalue(2,2,kk)
			END DO
            temp_phi_1=(temp_phi_1+((sigma(k,jj)-sigma(ii,jj))*temp_phivalue(2,2,jj)+(omiga(k,jj)**2-omiga(ii,jj)**2)*temp_lap(jj)/2.0))
end do
    !temp_phi=anti_MM(k,ii)*(temp_phi_1)+temp_phi
end do

sum_dis_eng=0

temp_phi=M*temp_phi_1
if (NN>1) then
evol_phi=(-2.0/NN) * temp_phi
else
evol_phi=0.0
end if
temp_phivalue(2,2,k)=temp_phivalue(2,2,k)+dt*evol_phi
end do
!排序
do ii=1,NN-1
do jj=ii+1,NN
if (temp_phivalue(2,2,jj)>temp_phivalue(2,2,ii)) then
tmp_temp_phivalue=temp_phivalue(2,2,ii)
temp_phivalue(2,2,ii)=temp_phivalue(2,2,jj)
temp_phivalue(2,2,jj)=tmp_temp_phivalue
tmp_BB=BB(ii)
BB(ii)=BB(jj)
BB(jj)=tmp_BB
end if
end do
end do
!删除第五个之后的值并归一化
value=0
name=0
sum_temp_phivalue=0.0
do ii=1,NN
    if(temp_phivalue(2,2,ii)<=0.0)then
        BB(ii)=0
        temp_phivalue(2,2,ii)=0
    end if
end do
do ii=1,NN
sum_temp_phivalue=sum_temp_phivalue+temp_phivalue(2,2,ii)
end do
if (NN<=5) then
value(1:NN)=temp_phivalue(2,2,:)/sum_temp_phivalue
name(1:NN)=BB(:)
else
value=temp_phivalue(2,2,1:5)/sum_temp_phivalue
name=BB(1:5)
end if
do ii=1,5
    if(value(ii)<=0.00001)then
        name(ii)=0
        value(ii)=0
    end if
end do
!赋值
phi_name_mid(i,j,:)=name
phi_value_mid(i,j,:)=value
deallocate(BB)
deallocate(temp_lap)
deallocate(temp_phivalue)
deallocate(sigma)
DEALLOCATE(sigma_3)
deallocate(omiga)
deallocate(anti_MM)

end do
!$OMP END DO
end do
!$OMP END PARALLEL
phi_name=phi_name_mid
phi_value=phi_value_mid

end subroutine gradient
!=============================================================================================================================
!=============================================================================================================================
subroutine microstructure
use gdata
implicit none
integer::i,j,k
mic=0.0
mij=0.0
do k=1,np
    do i=1,nx
        do j=1,ny
            mic(i,j)=mic(i,j) + phi_value(i,j,k)**2
        end do
    end do
end do
end subroutine microstructure
!=============================================================================================================================

    
subroutine account_S
use gdata
implicit none
integer::i,j,k

S=0.0

!do k=1,1
do i=1,nx
do j=1,ny

S(phi_name(i,j,1))=S(phi_name(i,j,1))+1

end do
end do
!end do

end subroutine account_S


!=============================================================================================================================
subroutine printdata(k)
use gdata

  IMPLICIT NONE
  INTEGER :: k,i,j
  CHARACTER(len=80) :: FileName1,FileName2

  WRITE(FileName1,FMT='(A4,I6.6,A4)') "mic_",k,".dat"
  WRITE(FileName2,FMT='(A4,I6.6,A4)') "miS_",k,".dat"

  FileName1 = TRIM(FileName1)
  FileName2 = TRIM(FileName2)

  OPEN(unit=52,file=FileName1,status="new")

  DO j=1,ny
     DO i=1,nx
        WRITE(52,FMT='(F10.4)',ADVANCE='no')mic(i,j)

     ENDDO
        WRITE(52,FMT='(F10.4)',ADVANCE='yes')

  ENDDO

  CLOSE(52)
  
  OPEN(unit=53,file=FileName2,status="new")

  DO i=1,QQ
    WRITE(53,FMT='(I6.6)',ADVANCE='yes')S(i)
  ENDDO

  CLOSE(53)

  PRINT*,k
  RETURN

end subroutine printdata
!==============================================================================================================================