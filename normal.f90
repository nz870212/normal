! A fortran95 PROGRAM for G95
! By NZ
MODULE gdata

	INTEGER, parameter ::  nx_0=612,ny_0=612,nx=512,ny=512,np=5,  Q=200,QQ=1000,parent_QQ=15,cores=24
	REAL(KIND=8), DIMENSION(nx, ny, np) :: lap_phi,ss,phi_value,phi_value_mid,dislocation,nucleation_rate_all,mis_angle,mis_angle_mid,parent_phi_value_mid
	REAL(KIND=8), DIMENSION(nx_0, ny_0, np) :: phi_value_0,dislocation_init,mis_angle_0,parent_phi_value
	REAL(KIND=8), DIMENSION(nx, ny) :: Num_p,sum_phi_value,r_c_i,r_c_i_real,N_sub_i
	REAL(KIND=8), DIMENSION(nx, ny) ::S_energy,R_grain,r_sub,r_sub_0,stress_i,nucleation_rate,mis_angle_show
	REAL(KIND=8) :: r_sub_init(QQ+1)
	INTEGER, DIMENSION(nx, ny, np) :: phi_name, phi_name_mid,parent_phi_name_mid
	INTEGER, DIMENSION(nx_0, ny_0, np) :: phi_name_0,parent_phi_name
	INTEGER :: T=1223,Qb=21866,Qa=386573,step_mid=100,itimestep
	REAL(KIND=8) :: dx, dy, dt, h=0.5e-6, gb_energy=0.4, gb_energy_sub,gb_energy_real=0.4, gb_width, grain_size,dislocation_0!=1.9e13!1.1229e13
	REAL(KIND=8) :: pi=3.1415926, M, M0, M_LAG, R=8.314, dE, timestep_nucleation, ll,stress_y,stress_s,strain_c,stress_c
	REAL(KIND=8) :: def_rate, alpha=0.5, b_vector=2.58e-10, miu, K_sub=23, yita, M_1, hard_rate,diameter,stress_0
	REAL(KIND=8) :: k1, k2, A1=175325.2, A2=0.1778, sigema_0, omega_0, nuclear_rate, dislocation_c,r_c,strain_tol,rate_all
	INTEGER :: num_Max_phi_name, num_may_nucleation, time_nucleation_1, time_nucleation_2,timestep
	!REAL(KIND=8),DIMENSION(20000+QQ*QQ) :: dislocation_average, strain, stress,radiu,volume
	!REAL(KIND=8),DIMENSION(5,20000+QQ*QQ) :: huizong
	REAL(KIND=8),DIMENSION(QQ+1) :: angle_value,M_array
	!REAL(KIND=8),DIMENSION(QQ+1) :: parent_angle_value,nucleation_rate_total
	
END MODULE gdata
!=============================================================================================================================
PROGRAM main

	USE gdata
	IMPLICIT NONE
	INTEGER :: nnn,MM0
	REAL(KIND=8) :: time_begin,time_end
	CALL CPU_TIME(time_begin)
	PRINT *, 'start'
	CALL init
	DO itimestep=1,10!itimestep
		CALL gradient
		CALL dislocation_evol_after		
		IF (itimestep>0) THEN
			PRINT*, itimestep
			IF (mod(itimestep,step_mid)==0) THEN			
				CALL microstructure(itimestep)			
				PRINT *, itimestep, "th iteration is running"
			END IF
		END IF	
	END DO
	CALL CPU_TIME(time_end)
	PRINT *, time_end-time_begin
	PRINT *, 'stop'
	!pause
	STOP
	
END PROGRAM main
!=============================================================================================================================
SUBROUTINE init

	USE gdata
	IMPLICIT NONE
	INTEGER :: i,j,k,l,mm
	REAL(KIND=8) :: d,dmin,a,b,e,f,mid,Sst,o,p,G_energy
	INTEGER, DIMENSION(1,parent_QQ) :: p_xx,p_yy
	INTEGER, DIMENSION(1,QQ) :: xx,yy
	
	dx=h
	dy=h
	gb_width=6*h
    M0=1.17e-7*exp(-Qb/(R*T))!
	
	mid=(3.0*b_vector*ll*yita**2.0*M_1)**(1.0/3.0)
	!修改界面迁移率
	M=(pi**2/(8.0*gb_width))*(M0)
    M_array=M0
	PRINT*, 'M =', M
	sigema_0=4.0*gb_energy/gb_width
	omega_0=(2.0/pi)*sqrt(2*gb_energy*gb_width)
	dt=0.5*dx**2/(4.0*((2.0/pi)*sqrt(2*gb_energy*gb_width))**2*(M))
	itimestep=ceiling(500+50/dt)
	PRINT*, 'timestep =', itimestep
    PRINT*, 'dt =', dt
    
    DO i=1,QQ
        CALL random_number(o)
        CALL random_number(p)        
        xx(1,i)=ceiling(nx_0*o)
        yy(1,i)=ceiling(nx_0*p)        
    END DO
	
    DO j=1,nx_0
        DO k=1,ny_0
            a=1.0*j
            b=1.0*k
            dmin=sqrt((a)**2+(b)**2)
            l=0
            DO mm=1,QQ
                e=1.0*xx(1,mm)
                f=1.0*yy(1,mm)
                d=sqrt(((e-a))**2+((f-b))**2)
                IF (d<=dmin) THEN
                    dmin=d
                    l=mm
                    CALL random_number(o)
                    phi_name_0(j,k,1)=l
                    phi_value_0(j,k,1)=1.0                   
                    mis_angle_0(j,k,:)=90*o
                END IF
            END DO        
        END DO
    END DO
	
    num_Max_phi_name=maxval(phi_name_0)
	
    DO j=1,nx_0
        DO k=1,ny_0
            IF(phi_value_0(j,k,1)==0.0) THEN
                phi_name_0(j,k,1)=num_Max_phi_name+1
                phi_value_0(j,k,1)=1.0               
                mis_angle_0(j,k,:)=15               
            END IF
        END DO
    END DO

    DO k=1,QQ+1
		DO i=1,nx_0
			DO j=1,ny_0
				IF (phi_name_0(i,j,1)==k) angle_value(k)=mis_angle_0(i,j,1)
			END DO
		END DO
    END DO
    
    DO i=1,nx
        DO j=1,ny
            phi_name(i,j,1)=phi_name_0(i+50,j+50,1)            
            phi_value(i,j,1)=phi_value_0(i+50,j+50,1)            
            mis_angle(i,j,:)=mis_angle_0(i+50,j+50,:)            
        END DO
    END DO
    
    phi_value_mid=phi_value
    phi_name_mid=phi_name
    mis_angle_mid=mis_angle
    num_Max_phi_name=maxval(phi_name_0)
	
    PRINT*, num_Max_phi_name
	PRINT *, 'intialized'
	
END SUBROUTINE init


!=============================================================================================================================

!主要计算步骤，计算每个时间步的相场该变量。
SUBROUTINE gradient
	USE gdata
	USE omp_lib
	IMPLICIT NONE
	INTEGER :: i,j,k,ip,jp,im,jm,N,NN=0,ii,jj,kk,num_cal,tmp_BB
	REAL(KIND=8) :: tmp_temp_phivalue,sum_temp_phivalue,sum_dis_eng,anti_gb_energy
	REAL(KIND=8) :: phis, phin, phiw, phie, phic, phisw, phine, phise, phinw
	REAL(KIND=8),external :: dfphi
	INTEGER,DIMENSION(45) :: A=0,B=0
	REAL(KIND=8),DIMENSION(3,3,np) :: temp_value,temp_angle
	REAL(KIND=8) :: temp_phi,temp_phi_1,temp_phi_2,evol_phi
	REAL(KIND=8) :: max_w,min_w,mid_w,w_i_jj_kk,w_ii_jj_kk
	REAL(KIND=8) :: w_12,w_13,w_23,M_ijk,n_M_ijk
	INTEGER,DIMENSION(3,3,np) :: temp_name
	INTEGER,DIMENSION(:),allocatable :: BB
	REAL(KIND=8),DIMENSION(:),allocatable :: temp_lap
	REAL(KIND=8),DIMENSION(:,:,:),allocatable :: temp_phivalue,temp_anglevalue
	REAL(KIND=8),DIMENSION(:,:),allocatable :: sigma,omiga,a_sigma,a_omiga,delta_angle,anti_MM
	REAL(KIND=8),DIMENSION(5) :: value
	INTEGER,DIMENSION(5) :: name
	!REAL,external::dfphi

	CALL omp_set_num_threads(cores)
	!$OMP PARALLEL PRIVATE(A,B,NN,num_cal,temp_value,temp_name,BB,N,&
	& max_w,min_w,mid_w,w_i_jj_kk,w_ii_jj_kk,temp_phi_2,M_ijk,w_12,w_13,w_23,&
	& temp_lap,temp_phivalue,temp_anglevalue,temp_angle,sigma,tmp_BB,&
	& omiga,i,j,ip,jp,im,jm,kk,ii,jj,k,sum_dis_eng,tmp_temp_phivalue,sum_temp_phivalue,value,name,&
	& a_sigma,a_omiga,delta_angle,anti_MM,anti_gb_energy,temp_phi,temp_phi_1,evol_phi)
	!pause
	DO i=1,nx
			!$OMP DO schedule(dynamic)
		DO j=1,ny
		
			A=0
			B=0
			!periodic boundary condition
			ip=mod((i+1),NX)
			im=mod((NX+i-1),NX)
			jp=mod((j+1),NY)
			jm=mod((NY+j-1),NY)
			IF (im==0) im=NX
			IF (jm==0) jm=NY
			IF (ip==0) ip=1
			IF (jp==0) jp=1
			
			DO kk=1,5
				temp_name(1,1,kk)=phi_name(im,jm,kk)
				temp_value(1,1,kk)=phi_value(im,jm,kk)
				temp_angle(1,1,kk)=mis_angle(im,jm,kk)
			END DO
			DO kk=1,5
				temp_name(1,2,kk)=phi_name(im,j,kk)
				temp_value(1,2,kk)=phi_value(im,j,kk)
				temp_angle(1,2,kk)=mis_angle(im,j,kk)
			END DO
			DO kk=1,5
				temp_name(1,3,kk)=phi_name(im,jp,kk)
				temp_value(1,3,kk)=phi_value(im,jp,kk)
				temp_angle(1,3,kk)=mis_angle(im,jp,kk)
			END DO
			DO kk=1,5
				temp_name(2,1,kk)=phi_name(i,jm,kk)
				temp_value(2,1,kk)=phi_value(i,jm,kk)
				temp_angle(2,1,kk)=mis_angle(i,jm,kk)
			END DO
			DO kk=1,5
				temp_name(2,2,kk)=phi_name(i,j,kk)
				temp_value(2,2,kk)=phi_value(i,j,kk)
				temp_angle(2,2,kk)=mis_angle(i,j,kk)
			END DO
			DO kk=1,5
				temp_name(2,3,kk)=phi_name(i,jp,kk)
				temp_value(2,3,kk)=phi_value(i,jp,kk)
				temp_angle(2,3,kk)=mis_angle(i,jp,kk)
			END DO
			DO kk=1,5
				temp_name(3,1,kk)=phi_name(ip,jm,kk)
				temp_value(3,1,kk)=phi_value(ip,jm,kk)
				temp_angle(3,1,kk)=mis_angle(ip,jm,kk)
			END DO
			DO kk=1,5
				temp_name(3,2,kk)=phi_name(ip,j,kk)
				temp_value(3,2,kk)=phi_value(ip,j,kk)
				temp_angle(3,2,kk)=mis_angle(ip,j,kk)
			END DO
			DO kk=1,5
				temp_name(3,3,kk)=phi_name(ip,jp,kk)
				temp_value(3,3,kk)=phi_value(ip,jp,kk)
				temp_angle(3,3,kk)=mis_angle(ip,jp,kk)
			END DO
			
			NN=0
			
			DO ii=1,3
				DO jj=1,3
					DO kk=1,np
						A(kk+3*(jj-1)+15*(ii-1))=temp_name(ii,jj,kk)
					END DO
				END DO
			END DO

			N=1
			B=0
			B(N)=A(1)
			
			III:  DO ii=2,45
			DO jj=1,N
				IF(B(jj)==A(ii) .or. A(ii)==0) CYCLE III
			END DO
			N=N+1
			B(N)=A(ii)
			END DO III
			
			DO ii=1,45
				IF (B(ii) .ne. 0) THEN
					NN=NN+1
				END IF
			END DO
			
			ALLOCATE(BB(NN))
			!BB:临时储存相场名的数组
			BB=B(1:NN)
			num_cal=SIZE(BB)
			
			ALLOCATE(temp_lap(NN))
			ALLOCATE(temp_phivalue(3,3,num_cal))
			ALLOCATE(temp_anglevalue(3,3,num_cal))
			ALLOCATE(sigma(NN,NN))
			ALLOCATE(omiga(NN,NN))
			ALLOCATE(a_omiga(NN,NN))
			ALLOCATE(a_sigma(NN,NN))
			ALLOCATE(delta_angle(NN,NN))
			ALLOCATE(anti_MM(NN,NN))
			
			temp_phivalue=0.0
			
			DO k=1,num_cal			!
				DO ii=1,3
					DO jj=1,3
						DO kk=1,5
							IF (temp_name(ii,jj,kk) .eq. BB(k)) THEN
								temp_phivalue(ii,jj,k)=temp_value(ii,jj,kk)
							END IF
						END DO
					END DO
				END DO

				DO ii=1,NN
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
				END DO
			END DO

			IF (NN == 3) THEN
			
				DO k=1,NN
					DO ii=1,NN
						DO jj=1,NN
							IF (jj==ii) THEN
								sigma(ii,jj)=0.0
								omiga(ii,jj)=0.0
							ELSE
								sigma(ii,jj)=omega_0
								omiga(ii,jj)=sigema_0
							END IF
						END DO
					END DO
					
					DO ii=1,NN
						DO jj=1,NN
							delta_angle(ii,jj)=min(abs(angle_value(BB(ii))-angle_value(BB(jj))),abs(90-(angle_value(BB(ii))-angle_value(BB(jj)))))
						END DO
					END DO
					
					DO ii=1,NN
						DO jj=1,NN
							IF (delta_angle(ii,jj)<14.5 .AND. delta_angle(ii,jj)>3) THEN
								anti_gb_energy=gb_energy*(delta_angle(ii,jj)/15.0)*(1-log(delta_angle(ii,jj)/15.0))
								anti_MM(ii,jj)=M_array(phi_name(i,j,1))*(1-exp(-5*(delta_angle(ii,jj)/15.0)**4))
							ELSE IF (delta_angle(ii,jj)<=3) THEN
								anti_gb_energy=gb_energy*(3.0/15.0)*(1-log(3.0/15.0))
								anti_MM(ii,jj)=M_array(phi_name(i,j,1))*(1-exp(-5*(3.0/15.0)**4))
							ELSE IF (delta_angle(ii,jj)>14.5) THEN
								anti_gb_energy=gb_energy
								anti_MM(ii,jj)=M_array(phi_name(i,j,1))
							END IF
							!IF (itimestep<=100) anti_MM=(pi**2/(8.0*gb_width))*M0
							IF (jj==ii) THEN
								a_sigma(ii,jj)=0.0
								a_omiga(ii,jj)=0.0							
							ELSE
								a_sigma(ii,jj)=2.0*sqrt(2*anti_gb_energy*gb_width)/3.1417926
								a_omiga(ii,jj)=4*anti_gb_energy/gb_width
							END IF
						END DO
					END DO
		
					temp_phi_1=0.0
					temp_phi=0.0
					w_12=0
					w_13=0
					w_23=0
					n_M_ijk=1
					!
					w_12=0.5*(1-(ABS(anti_MM(1,2)-(anti_MM(1,2)+anti_MM(1,3)+anti_MM(2,3))/3))**n_M_ijk &
					& /((ABS(anti_MM(1,2)-(anti_MM(1,2)+anti_MM(1,3)+anti_MM(2,3))/3))**n_M_ijk &
					& +(ABS(anti_MM(1,3)-(anti_MM(1,2)+anti_MM(1,3)+anti_MM(2,3))/3))**n_M_ijk &
					& +(ABS(anti_MM(2,3)-(anti_MM(1,2)+anti_MM(1,3)+anti_MM(2,3))/3))**n_M_ijk))
					w_13=0.5*(1-(ABS(anti_MM(1,3)-(anti_MM(1,2)+anti_MM(1,3)+anti_MM(2,3))/3))**n_M_ijk &
					& /((ABS(anti_MM(1,2)-(anti_MM(1,2)+anti_MM(1,3)+anti_MM(2,3))/3))**n_M_ijk &
					& +(ABS(anti_MM(1,3)-(anti_MM(1,2)+anti_MM(1,3)+anti_MM(2,3))/3))**n_M_ijk &
					& +(ABS(anti_MM(2,3)-(anti_MM(1,2)+anti_MM(1,3)+anti_MM(2,3))/3))**n_M_ijk))
					w_23=0.5*(1-(ABS(anti_MM(2,3)-(anti_MM(1,2)+anti_MM(1,3)+anti_MM(2,3))/3))**n_M_ijk &
					& /((ABS(anti_MM(1,2)-(anti_MM(1,2)+anti_MM(1,3)+anti_MM(2,3))/3))**n_M_ijk &
					& +(ABS(anti_MM(1,3)-(anti_MM(1,2)+anti_MM(1,3)+anti_MM(2,3))/3))**n_M_ijk &
					& +(ABS(anti_MM(2,3)-(anti_MM(1,2)+anti_MM(1,3)+anti_MM(2,3))/3))**n_M_ijk))
					M_ijk=(pi**2/(8.0*gb_width))*(anti_MM(1,2)*w_12+anti_MM(1,3)*w_13+anti_MM(2,3)*w_23)
				
					PRINT *, w_12,w_13,w_23
					
					DO ii=1,NN
						DO jj=1,NN
							temp_phi_2=0
							DO kk=1,NN
								max_w=0
								mid_w=0
								min_w=0
								max_w=MAX(a_sigma(k,jj),a_sigma(k,kk),a_sigma(jj,kk))
								min_w=MIN(a_sigma(k,jj),a_sigma(k,kk),a_sigma(jj,kk))
								mid_w=a_sigma(k,jj)+a_sigma(k,kk)+a_sigma(jj,kk)-max_w-min_w
								IF ((max_w+min_w-2*mid_w) > 0) THEN
									w_i_jj_kk=5*(max_w-(mid_w+min_w)/2)
								ELSE
									w_i_jj_kk=0
								END IF
								IF (min_w == 0 .AND. w_i_jj_kk .NE. 0) PRINT *, max_w,mid_w,min_w,w_i_jj_kk
								
								max_w=0
								mid_w=0
								min_w=0
								max_w=MAX(a_sigma(ii,jj),a_sigma(ii,kk),a_sigma(jj,kk))
								min_w=MIN(a_sigma(ii,jj),a_sigma(ii,kk),a_sigma(jj,kk))
								mid_w=a_sigma(ii,jj)+a_sigma(ii,kk)+a_sigma(jj,kk)-max_w-min_w
								IF ((max_w+min_w-2*mid_w) > 0) THEN
									w_ii_jj_kk=5*(max_w-(mid_w+min_w)/2)
								ELSE
									w_ii_jj_kk=0
								END IF
								IF (min_w == 0 .AND. w_ii_jj_kk .NE. 0) PRINT *, max_w,mid_w,min_w,w_ii_jj_kk
								
								temp_phi_2=temp_phi_2+(w_i_jj_kk-w_ii_jj_kk)*temp_phivalue(2,2,jj)*temp_phivalue(2,2,kk)
								
							END DO
							
							IF (itimestep<=100) THEN
								temp_phi_1=(temp_phi_1+((omiga(k,jj)-omiga(ii,jj))*temp_phivalue(2,2,jj)+(sigma(k,jj)**2-sigma(ii,jj)**2)*temp_lap(jj)/2.0))
							ELSE
								temp_phi_1=temp_phi_2+(temp_phi_1+((a_omiga(k,jj)-a_omiga(ii,jj))*temp_phivalue(2,2,jj)+(a_sigma(k,jj)**2-a_sigma(ii,jj)**2)&
								&*temp_lap(jj)/2.0))
							END IF
						END DO
						IF (itimestep<=100) THEN
							temp_phi=M*(temp_phi_1)+temp_phi
						ELSE
							temp_phi=M_ijk*(temp_phi_1)+temp_phi
						END IF
					END DO			
					
					IF (NN>1 .AND. itimestep<=200) THEN
						evol_phi=(-2.0/NN) * temp_phi
					ELSEIF (NN>1 .AND. itimestep>200) THEN
						evol_phi=(-2.0/NN) * temp_phi
					ELSE
						evol_phi=0.0			
					END IF
						temp_phivalue(2,2,k)=temp_phivalue(2,2,k)+dt*evol_phi
				END DO
				
			ELSE
			
				DO k=1,NN
					DO ii=1,NN
						DO jj=1,NN
							IF (jj==ii) THEN
								sigma(ii,jj)=0.0
								omiga(ii,jj)=0.0
							ELSE
								sigma(ii,jj)=omega_0
								omiga(ii,jj)=sigema_0
							END IF
						END DO
					END DO
					
					DO ii=1,NN
						DO jj=1,NN
							delta_angle(ii,jj)=min(abs(angle_value(BB(ii))-angle_value(BB(jj))),abs(90-(angle_value(BB(ii))-angle_value(BB(jj)))))
						END DO
					END DO
					
					DO ii=1,NN
						DO jj=1,NN
							IF (delta_angle(ii,jj)<14.5 .AND. delta_angle(ii,jj)>3) THEN
								anti_gb_energy=gb_energy*(delta_angle(ii,jj)/15.0)*(1-log(delta_angle(ii,jj)/15.0))
								anti_MM(ii,jj)=M_array(phi_name(i,j,1))*(1-exp(-5*(delta_angle(ii,jj)/15.0)**4))
							ELSE IF (delta_angle(ii,jj)<=3) THEN
								anti_gb_energy=gb_energy*(3.0/15.0)*(1-log(3.0/15.0))
								anti_MM(ii,jj)=M_array(phi_name(i,j,1))*(1-exp(-5*(3.0/15.0)**4))
							ELSE IF (delta_angle(ii,jj)>14.5) THEN
								anti_gb_energy=gb_energy
								anti_MM(ii,jj)=M_array(phi_name(i,j,1))
							END IF
							IF (itimestep<=100) anti_MM=M0
							IF (jj==ii) THEN
								a_sigma(ii,jj)=0.0
								a_omiga(ii,jj)=0.0
							ELSE
								a_sigma(ii,jj)=2.0*sqrt(2*anti_gb_energy*gb_width)/3.1417926
								a_omiga(ii,jj)=4*anti_gb_energy/gb_width
							END IF
						END DO
					END DO
					
					anti_MM=(pi**2/(8.0*gb_width))*anti_MM
					temp_phi_1=0.0
					temp_phi=0.0
					
					DO ii=1,NN
						DO jj=1,NN													
							IF (itimestep<=100) THEN
								temp_phi_1=(temp_phi_1+((omiga(k,jj)-omiga(ii,jj))*temp_phivalue(2,2,jj)+(sigma(k,jj)**2-sigma(ii,jj)**2)*temp_lap(jj)/2.0))
							ELSE
								temp_phi_1=(temp_phi_1+((a_omiga(k,jj)-a_omiga(ii,jj))*temp_phivalue(2,2,jj)+(a_sigma(k,jj)**2-a_sigma(ii,jj)**2)&
								&*temp_lap(jj)/2.0))
							END IF
						END DO
							temp_phi=anti_MM(k,ii)*(temp_phi_1)+temp_phi
					END DO			
					
					IF (NN>1 .AND. itimestep<=200) THEN
						evol_phi=(-2.0/NN) * temp_phi
					ELSEIF (NN>1 .AND. itimestep>200) THEN
						evol_phi=(-2.0/NN) * temp_phi
					ELSE
						evol_phi=0.0			
					END IF
						temp_phivalue(2,2,k)=temp_phivalue(2,2,k)+dt*evol_phi
				END DO
			END IF
			!排序
			DO ii=1,NN-1
				DO jj=ii+1,NN
					IF (temp_phivalue(2,2,jj)>temp_phivalue(2,2,ii)) THEN
						tmp_temp_phivalue=temp_phivalue(2,2,ii)
						temp_phivalue(2,2,ii)=temp_phivalue(2,2,jj)
						temp_phivalue(2,2,jj)=tmp_temp_phivalue
						tmp_BB=BB(ii)
						BB(ii)=BB(jj)
						BB(jj)=tmp_BB
					END IF
				END DO
			END DO
			!删除第五个之后的值并归一化
			value=0
			name=0
			sum_temp_phivalue=0.0
			
			DO ii=1,NN
				IF(temp_phivalue(2,2,ii)<=0.0)THEN
					BB(ii)=0
					temp_phivalue(2,2,ii)=0
				END IF
			END DO
			
			DO ii=1,NN
				sum_temp_phivalue=sum_temp_phivalue+temp_phivalue(2,2,ii)
			END DO
			
			IF (NN<=5) THEN
				value(1:NN)=temp_phivalue(2,2,:)/sum_temp_phivalue
				name(1:NN)=BB(:)
			ELSE
				value=temp_phivalue(2,2,1:5)/sum_temp_phivalue
				name=BB(1:5)
			END IF
			
			DO ii=1,5
				IF(value(ii)<=0.00001) THEN
					name(ii)=0
					value(ii)=0
				END IF
			END DO
			!赋值
			phi_name_mid(i,j,:)=name
			phi_value_mid(i,j,:)=value
			DEALLOCATE(BB)
			DEALLOCATE(temp_lap)
			DEALLOCATE(temp_anglevalue)
			DEALLOCATE(temp_phivalue)
			DEALLOCATE(sigma)
			DEALLOCATE(omiga)
			DEALLOCATE(a_sigma)
			DEALLOCATE(a_omiga)
			DEALLOCATE(delta_angle)
			DEALLOCATE(anti_MM)
		
		END DO
		!$OMP END DO
	END DO
	!$OMP END PARALLEL
	
	phi_name=phi_name_mid
	phi_value=phi_value_mid
	
END SUBROUTINE gradient
!=============================================================================================================================
SUBROUTINE dislocation_evol_after

	USE gdata
	IMPLICIT NONE
	INTEGER :: i,j,k
	
	DO i=1,nx
		DO j=1,ny
			DO k=1,np
				IF (phi_name(i,j,k)>0) THEN
					mis_angle(i,j,k)=angle_value(phi_name(i,j,k))
				ELSE
					mis_angle(i,j,k)=angle_value(phi_name(i,j,1))
				END IF
			END DO
		END DO
	END DO

END SUBROUTINE dislocation_evol_after
!=============================================================================================================================
SUBROUTINE microstructure(l)

	USE gdata
	IMPLICIT NONE
	INTEGER::i,j,k,l
	CHARACTER(len=80) :: FileName1,FileName2
	REAL(KIND=8), DIMENSION(:,:,:), allocatable :: delta_mis_angle
	REAL(KIND=8), DIMENSION(:,:), allocatable :: mic,mij
	
	ALLOCATE(delta_mis_angle(nx,ny,np))
	ALLOCATE(mij(nx,ny))
	ALLOCATE(mic(nx,ny))
	mic=0.0
	mij=0.0
	
	DO j=1,np
		delta_mis_angle(:,:,j)=min(abs(mis_angle(:,:,j)-mis_angle(:,:,1)),abs(90-(mis_angle(:,:,j)-mis_angle(:,:,1))))
	END DO
	
	mij=delta_mis_angle(:,:,2)
	
	DEALLOCATE(delta_mis_angle)
	
	DO k=1,np
		DO i=1,nx
			DO j=1,ny
				mic(i,j)=mic(i,j) + phi_value(i,j,k)**2       
			END DO
		END DO
	END DO
	
	WRITE(FileName1,FMT='(A4,I6.6,A4)') "mic_",l,".dat"
	WRITE(FileName2,FMT='(A4,I6.6,A4)') "mij_",l,".dat"
	FileName1=TRIM(FileName1)
	FileName2=TRIM(FileName2)
	OPEN(unit=52,file=FileName1,status="new")
	
	DO j=1,ny
		DO i=1,nx
			WRITE(52,FMT='(F10.4)',ADVANCE='no')mic(i,j)
		END DO
		WRITE(52,FMT='(F10.4)',ADVANCE='yes')
	END DO
	
	CLOSE(52)
	
	OPEN(unit=53,file=FileName2,status="new")
	
	DO j=1,ny
		DO i=1,nx
			WRITE(53,FMT='(F10.4)',ADVANCE='no')mij(i,j)
		END DO
		WRITE(53,FMT='(F10.4)',ADVANCE='yes')
	END DO
	
	CLOSE(53)
	
	PRINT*,k
	
	RETURN
	
	DEALLOCATE(mij)
	DEALLOCATE(mic)
	
END SUBROUTINE microstructure

