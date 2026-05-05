      program SYDPDZC2
      implicit none
c>>  serial number: Lx
c>>  L: Receptor-Ligand interaction;
c>>  Dissociation rate x: 0 means none (50.0);
c>>                       1 means 0.5
c>>                       2 means 0.05
c>>                       3 means 0.005
c>>                       4 means 0.0005
c>>                       5 means 0.00005

      integer type_num
      parameter (type_num=4)
      real*8 cell_range_x
      parameter (cell_range_x=200.0)
      real*8 cell_range_y
      parameter (cell_range_y=200.0)
      real*8 cell_range_z
      parameter (cell_range_z=100.0)
      real*8 dt
      integer npart_max,npart_tot  
      parameter (npart_max=100)
      integer intranum_max  
      parameter (intranum_max=10)
      parameter (dt=0.02)
      integer nsimu
      parameter (nsimu=10000000)
      integer niter
      parameter (niter=100)
      real*8 iter_tol
      parameter (iter_tol=0.0001)
      real*8 cm_frict_const
      parameter (cm_frict_const=4.0)
      real*8 cc_frict_const
      parameter (cc_frict_const=10.0)
      real*8 spring_length_rp
      parameter (spring_length_rp=20.0)
      real*8 spring_length_bd
      parameter (spring_length_bd=50.0)
      integer complex_nb_ctg_num
      parameter (complex_nb_ctg_num=2) ! the maximal number of subunits in complex
      real*8 pai
      parameter (pai=3.1415926) 
      real*8 reaction_range
      parameter (reaction_range=0.5)
      real*8 reaction_shell
      parameter (reaction_shell=0.5)
      real*8 bond_thetapd,bond_thetapd_cutoff
      parameter (bond_thetapd=180.0,bond_thetapd_cutoff=10.0)
      real*8 bond_thetaot,bond_thetaot_cutoff
      parameter (bond_thetaot=180.0,bond_thetaot_cutoff=10.0)

      real*8 Prob_Dettach,Dettach_Rate,Attach_Prob
      parameter (Dettach_Rate=0.0,Attach_Prob=0.0)

      integer i,j,k,ii,jj
      integer npart(type_num)
      integer intra_num(type_num)
      real*8 radius(type_num)
      real*8 linker(type_num)
      real*8 Pol_radius(type_num)
      real*8 diffu_const(type_num)
      real*8 fluct_const(type_num)
      real*8 rtdiffu_const(type_num)
      real*8 spring_length_plm(type_num)
      integer intra_geo_matrix(intranum_max,intranum_max,type_num)
      real*8 intrabond_length(intranum_max,intranum_max,type_num)
      real*8 x(intranum_max,npart_max,type_num)
      real*8 y(intranum_max,npart_max,type_num)
      real*8 z(intranum_max,npart_max,type_num)
      real*8 temp_x(intranum_max,npart_max,type_num)
      real*8 temp_y(intranum_max,npart_max,type_num)
      real*8 temp_z(intranum_max,npart_max,type_num)
      real*8 x0(intranum_max,npart_max,type_num)
      real*8 y0(intranum_max,npart_max,type_num)
      real*8 z0(intranum_max,npart_max,type_num)
      real*8 f_rp_x(intranum_max,npart_max,type_num)
      real*8 f_rp_y(intranum_max,npart_max,type_num)
      real*8 f_rp_z(intranum_max,npart_max,type_num)
      real*8 f_plm_x(intranum_max,npart_max,type_num)
      real*8 f_plm_y(intranum_max,npart_max,type_num)
      real*8 f_plm_z(intranum_max,npart_max,type_num)
      real*8 f_bd_x(intranum_max,npart_max,type_num)
      real*8 f_bd_y(intranum_max,npart_max,type_num)
      real*8 f_bd_z(intranum_max,npart_max,type_num)
      real*8 f_rdm_x(intranum_max,npart_max,type_num)
      real*8 f_rdm_y(intranum_max,npart_max,type_num)
      real*8 f_rdm_z(intranum_max,npart_max,type_num)
      real*8 f_rdmrt_x(intranum_max,npart_max,type_num)
      real*8 f_rdmrt_y(intranum_max,npart_max,type_num)
      real*8 f_rdmrt_z(intranum_max,npart_max,type_num)
      real*8 f_rdmfl_x(intranum_max,npart_max,type_num)
      real*8 f_rdmfl_y(intranum_max,npart_max,type_num)
      real*8 f_rdmfl_z(intranum_max,npart_max,type_num)
      real*8 vx(intranum_max,npart_max,type_num)
      real*8 vy(intranum_max,npart_max,type_num)
      real*8 vz(intranum_max,npart_max,type_num)
      real*8 vx_old(intranum_max,npart_max,type_num)
      real*8 vy_old(intranum_max,npart_max,type_num)
      real*8 vz_old(intranum_max,npart_max,type_num)
      real*8 vx_new(intranum_max,npart_max,type_num)
      real*8 vy_new(intranum_max,npart_max,type_num)
      real*8 vz_new(intranum_max,npart_max,type_num)
      real*8 ax(intranum_max,npart_max,type_num,
     &     intranum_max,npart_max,type_num)
      real*8 ay(intranum_max,npart_max,type_num,
     &     intranum_max,npart_max,type_num)
      real*8 az(intranum_max,npart_max,type_num,
     &     intranum_max,npart_max,type_num)
      real*8 bx(intranum_max,npart_max,type_num)
      real*8 by(intranum_max,npart_max,type_num)
      real*8 bz(intranum_max,npart_max,type_num)
      integer status(intranum_max,npart_max,type_num)
      integer status_new(intranum_max,npart_max,type_num)
      real*8 bond_dist_cutoff
      real*8 dij
      integer itime
      integer molecule,neighbor
      real*8 RMSD
      real*8 temp_i,temp_j,temp_k
      real*8 dir_x,dir_y,dir_z
      real*8 amp_dir
      real*8 scaled_dir_x,scaled_dir_y,scaled_dir_z
      integer complex_num,complex_num_new ! Number of complex formed in simulation
      integer complex_nb_num(npart_max*intranum_max*type_num)
      integer complex_nb_num_new(npart_max*intranum_max*type_num) ! How many subunits in each complex. Here the max number is 2.
      integer complex_nb_ctg
     &     (npart_max*intranum_max*type_num,complex_nb_ctg_num)
      integer complex_nb_ctg_new
     &     (npart_max*intranum_max*type_num,complex_nb_ctg_num) ! e.g. for complex_nb_ctg(i,j), what is the molecule category of jth subunit of ith complex. (here 1 is RB_A and 2 is RB_B)
      integer complex_nb_ctg_idx
     &     (npart_max*intranum_max*type_num,complex_nb_ctg_num)
      integer complex_nb_ctg_idx_new
     &     (npart_max*intranum_max*type_num,complex_nb_ctg_num) ! e.g. for complex_nb_ctg(i,j), what is the molecule index of jth subunit of ith complex, given its category=complex_nb_ctg(i,j).
      integer complex_nb_ctg_plm_idx
     &     (npart_max*intranum_max*type_num,complex_nb_ctg_num)
      integer complex_nb_ctg_plm_idx_new
     &     (npart_max*intranum_max*type_num,complex_nb_ctg_num) ! e.g. for complex_nb_ctg(i,j), what is the atom index in the polymer of ith complex, given its category=complex_nb_ctg(i,j).
      real*8 Prob_Ass,Prob_Diss
      integer selected_complex
      integer selected_ctg_A,selected_ctg_B
      integer selected_idx_A,selected_idx_B
      integer selected_plm_idx_A,selected_plm_idx_B
      real*8 prob
      real*8 temp
      integer part1,part2,part3
      integer part4,part5,part6
      integer part7,part8,part9
      integer reaction_matrix
     &     (intranum_max,type_num,intranum_max,type_num)
      real*8 Ass_Rate(intranum_max,type_num,intranum_max,type_num)
      real*8 Diss_Rate(intranum_max,type_num,intranum_max,type_num)
      real*8 dist
      real*8 dx,dy,dz
      real*8 theta,phi,psai,phai
      real*8 t(3,3)
      real*8 point_x(3),point_y(3),point_z(3)
      real*8 theta_pd,theta_ot
      integer angle_flag
      real*8 cm_x,cm_y,cm_z
      real*8 amp_rdmrt
      integer dimer_num
      integer check_flag,dimer_nb

      integer oligomer_num
      integer oligomer_subunit_num(npart_max*2)
      integer oligomer_subunit_index(npart_max*2,npart_max*2)
      integer oligomer_state(npart_max*2)
      integer contact_matrix(npart_max*2,npart_max*2)
      integer color_flag(npart_max*2)
      integer oligomer_map(npart_max*2,npart_max*2)
      integer count_flag
      integer complex_mol_idx2(npart_max*intranum_max,2)
      integer mol_tot_num
      integer C2_status(npart_max)
      integer PDZ_status(npart_max)
      integer LR_num,RS_num,SE_num,RR_num

      integer index_mol,Ass_Flag,index_mol2
	  
      real rand3
      double precision r3
      real rand4
      double precision r4
      real rand5
      double precision r5

      r3=5.0
      r4=5.0
      r5=5.0    

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c   construct the initialized position of molecules in 2D
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c>>  customize concentration

      npart(1)=50
      npart(2)=50
      npart(3)=50
      npart(4)=50
	  
      intra_num(1)=2 ! L
      intra_num(2)=4 ! R
      intra_num(3)=4 ! S
      intra_num(4)=4 ! E

      radius(1)=1
      radius(2)=1
      radius(3)=1
      radius(4)=1
	  
      bond_dist_cutoff=2.0+reaction_range

      linker(1)=1
      linker(2)=1
      linker(3)=3
      linker(4)=3

      Pol_radius(1)=10
      Pol_radius(2)=10
      Pol_radius(3)=10
      Pol_radius(4)=10

      diffu_const(1)=5.0
      diffu_const(2)=5.0  
      diffu_const(3)=5.0
      diffu_const(4)=5.0  

      rtdiffu_const(1)=0.5
      rtdiffu_const(2)=0.5
      rtdiffu_const(3)=0.5
      rtdiffu_const(4)=0.5
	  
      fluct_const(1)=5.0
      fluct_const(2)=5.0
      fluct_const(3)=5.0
      fluct_const(4)=5.0
	  
      spring_length_plm(1)=10.0
      spring_length_plm(2)=10.0
      spring_length_plm(3)=10.0
      spring_length_plm(4)=10.0
      
cc>>  set up initial configuration

      do k=1,type_num
         do i=1,npart(k)
 200        continue
            if(k.eq.1)then

               temp_i=rand3(r3)*cell_range_x-cell_range_x/2  
               temp_j=rand3(r3)*cell_range_y-cell_range_y/2  
               temp_k=rand3(r3)*cell_range_z/2       

               do ii=1,k
                  do jj=1,i-1
                     dist=sqrt((temp_i-x(1,jj,ii))**2+
     &                    (temp_j-y(1,jj,ii))**2+
     &                    (temp_k-z(1,jj,ii))**2)
                     if(dist.le.
     &                    (Pol_radius(k)+Pol_radius(ii)))then
                        goto 200
                     endif
                  enddo
               enddo

               x(1,i,k)=temp_i-0.5*(2*radius(k)+linker(k))
               y(1,i,k)=temp_j
               z(1,i,k)=temp_k
               x0(1,i,k)=temp_i-0.5*(2*radius(k)+linker(k))
               y0(1,i,k)=temp_j
               z0(1,i,k)=temp_k
               status(1,i,k)=0

               x(2,i,k)=temp_i+0.5*(2*radius(k)+linker(k))
               y(2,i,k)=temp_j
               z(2,i,k)=temp_k
               x0(2,i,k)=temp_i+0.5*(2*radius(k)+linker(k))
               y0(2,i,k)=temp_j
               z0(2,i,k)=temp_k
               status(2,i,k)=0


            elseif(k.eq.2)then

               temp_i=rand3(r3)*cell_range_x-cell_range_x/2  
               temp_j=rand3(r3)*cell_range_y-cell_range_y/2  
               temp_k=0

               do ii=1,k
                  do jj=1,i-1
                     dist=sqrt((temp_i-x(1,jj,ii))**2+
     &                    (temp_j-y(1,jj,ii))**2+
     &                    (temp_k-z(1,jj,ii))**2)
                     if(dist.le.
     &                    (Pol_radius(k)+Pol_radius(ii)))then
                        goto 200
                     endif
                  enddo
               enddo

               x(1,i,k)=temp_i
               y(1,i,k)=temp_j
               z(1,i,k)=2*(2*radius(k)+linker(k))
               x0(1,i,k)=temp_i
               y0(1,i,k)=temp_j
               z0(1,i,k)=2*(2*radius(k)+linker(k))
               status(1,i,k)=0
               x(2,i,k)=temp_i
               y(2,i,k)=temp_j
               z(2,i,k)=1*(2*radius(k)+linker(k))
               x0(2,i,k)=temp_i
               y0(2,i,k)=temp_j
               z0(2,i,k)=1*(2*radius(k)+linker(k))
               status(2,i,k)=0
               x(3,i,k)=temp_i
               y(3,i,k)=temp_j
               z(3,i,k)=0
               x0(3,i,k)=temp_i
               y0(3,i,k)=temp_j
               z0(3,i,k)=0
               status(3,i,k)=0
               x(4,i,k)=temp_i
               y(4,i,k)=temp_j
               z(4,i,k)=-1*(2*radius(k)+linker(k))
               x0(4,i,k)=temp_i
               y0(4,i,k)=temp_j
               z0(4,i,k)=-1*(2*radius(k)+linker(k))
               status(4,i,k)=0

            elseif(k.eq.3)then

               temp_i=rand3(r3)*cell_range_x-cell_range_x/2  
               temp_j=rand3(r3)*cell_range_y-cell_range_y/2  
               temp_k=-rand3(r3)*cell_range_z/2       

               do ii=1,k
                  do jj=1,i-1
                     dist=sqrt((temp_i-x(1,jj,ii))**2+
     &                    (temp_j-y(1,jj,ii))**2+
     &                    (temp_k-z(1,jj,ii))**2)
                     if(dist.le.
     &                    (Pol_radius(k)+Pol_radius(ii)))then
                        goto 200
                     endif
                  enddo
               enddo

               x(1,i,k)=temp_i
               y(1,i,k)=temp_j
               z(1,i,k)=temp_k+1.5*(2*radius(k)+linker(k))
               x0(1,i,k)=temp_i
               y0(1,i,k)=temp_j
               z0(1,i,k)=temp_k+1.5*(2*radius(k)+linker(k))
               status(1,i,k)=0
               x(2,i,k)=temp_i
               y(2,i,k)=temp_j
               z(2,i,k)=temp_k+0.5*(2*radius(k)+linker(k))
               x0(2,i,k)=temp_i
               y0(2,i,k)=temp_j
               z0(2,i,k)=temp_k+0.5*(2*radius(k)+linker(k))
               status(2,i,k)=0
               x(3,i,k)=temp_i
               y(3,i,k)=temp_j
               z(3,i,k)=temp_k-0.5*(2*radius(k)+linker(k))
               x0(3,i,k)=temp_i
               y0(3,i,k)=temp_j
               z0(3,i,k)=temp_k-0.5*(2*radius(k)+linker(k))
               status(3,i,k)=0
               x(4,i,k)=temp_i
               y(4,i,k)=temp_j
               z(4,i,k)=temp_k-1.5*(2*radius(k)+linker(k))
               x0(4,i,k)=temp_i
               y0(4,i,k)=temp_j
               z0(4,i,k)=temp_k-1.5*(2*radius(k)+linker(k))
               status(4,i,k)=0

            elseif(k.eq.4)then

               temp_i=rand3(r3)*cell_range_x-cell_range_x/2  
               temp_j=rand3(r3)*cell_range_y-cell_range_y/2  
               temp_k=-rand3(r3)*cell_range_z/2       

               do ii=1,k
                  do jj=1,i-1
                     dist=sqrt((temp_i-x(1,jj,ii))**2+
     &                    (temp_j-y(1,jj,ii))**2+
     &                    (temp_k-z(1,jj,ii))**2)
                     if(dist.le.
     &                    (Pol_radius(k)+Pol_radius(ii)))then
                        goto 200
                     endif
                  enddo
               enddo

               x(1,i,k)=temp_i
               y(1,i,k)=temp_j
               z(1,i,k)=temp_k+1.5*(2*radius(k)+linker(k))
               x0(1,i,k)=temp_i
               y0(1,i,k)=temp_j
               z0(1,i,k)=temp_k+1.5*(2*radius(k)+linker(k))
               status(1,i,k)=0
               x(2,i,k)=temp_i
               y(2,i,k)=temp_j
               z(2,i,k)=temp_k+0.5*(2*radius(k)+linker(k))
               x0(2,i,k)=temp_i
               y0(2,i,k)=temp_j
               z0(2,i,k)=temp_k+0.5*(2*radius(k)+linker(k))
               status(2,i,k)=0
               x(3,i,k)=temp_i
               y(3,i,k)=temp_j
               z(3,i,k)=temp_k-0.5*(2*radius(k)+linker(k))
               x0(3,i,k)=temp_i
               y0(3,i,k)=temp_j
               z0(3,i,k)=temp_k-0.5*(2*radius(k)+linker(k))
               status(3,i,k)=0
               x(4,i,k)=temp_i
               y(4,i,k)=temp_j
               z(4,i,k)=temp_k-1.5*(2*radius(k)+linker(k))
               x0(4,i,k)=temp_i
               y0(4,i,k)=temp_j
               z0(4,i,k)=temp_k-1.5*(2*radius(k)+linker(k))
               status(4,i,k)=0
           
            endif


			
c>>>>>   random orientation

            if(k.ne.2)then

               cm_x=0
               cm_y=0
               cm_z=0

               do j=1,intra_num(k)
                  cm_x=cm_x+x0(j,i,k)
                  cm_y=cm_y+y0(j,i,k)
                  cm_z=cm_z+z0(j,i,k)
               enddo
               
               cm_x=cm_x/real(intra_num(k))
               cm_y=cm_y/real(intra_num(k))
               cm_z=cm_z/real(intra_num(k))

               theta=(2*rand3(r3)-1)*pai
               phi=(2*rand3(r3)-1)*pai
               psai=(2*rand3(r3)-1)*pai
            
               t(1,1)=cos(psai)*cos(phi)-cos(theta)*sin(phi)*sin(psai)
               t(1,2)=-sin(psai)*cos(phi)-cos(theta)*sin(phi)*cos(psai)
               t(1,3)=sin(theta)*sin(phi)
               
               t(2,1)=cos(psai)*sin(phi)+cos(theta)*cos(phi)*sin(psai)
               t(2,2)=-sin(psai)*sin(phi)+cos(theta)*cos(phi)*cos(psai)
               t(2,3)=-sin(theta)*cos(phi)
               
               t(3,1)=sin(psai)*sin(theta)
               t(3,2)=cos(psai)*sin(theta)
               t(3,3)=cos(theta)
               
               
               do j=1,intra_num(k)
                  x(j,i,k)=t(1,1)*(x0(j,i,k)-cm_x)+
     &                 t(1,2)*(y0(j,i,k)-cm_y)
     &                 +t(1,3)*(z0(j,i,k)-cm_z)
     &                 +cm_x
                  y(j,i,k)=t(2,1)*(x0(j,i,k)-cm_x)+
     &                 t(2,2)*(y0(j,i,k)-cm_y)
     &                 +t(2,3)*(z0(j,i,k)-cm_z)
     &                 +cm_y
                  z(j,i,k)=t(3,1)*(x0(j,i,k)-cm_x)+
     &                 t(3,2)*(y0(j,i,k)-cm_y)
     &                 +t(3,3)*(z0(j,i,k)-cm_z)
     &                 +cm_z
               enddo
            endif
         enddo

      enddo

cc>> customize geometry matrix

      do i=1,type_num
         do ii=1,intranum_max
            do jj=1,intranum_max
               intra_geo_matrix(ii,jj,i)=0
               intrabond_length(ii,jj,i)=0
            enddo
         enddo
      enddo

      intra_geo_matrix(1,2,1)=1
	  
      intra_geo_matrix(1,2,2)=1
      intra_geo_matrix(2,3,2)=1
      intra_geo_matrix(3,4,2)=1

      intra_geo_matrix(1,2,3)=1
      intra_geo_matrix(2,3,3)=1
      intra_geo_matrix(3,4,3)=1

      intra_geo_matrix(1,2,4)=1
      intra_geo_matrix(2,3,4)=1
      intra_geo_matrix(3,4,4)=1

      intrabond_length(1,2,1)=2*radius(1)+linker(1)

      intrabond_length(1,2,2)=2*radius(2)+linker(2)
      intrabond_length(2,3,2)=2*radius(2)+linker(2)
      intrabond_length(3,4,2)=2*radius(2)+linker(2)

      intrabond_length(1,2,3)=2*radius(3)+linker(3)
      intrabond_length(2,3,3)=2*radius(3)+linker(3)
      intrabond_length(3,4,3)=2*radius(3)+linker(3)

      intrabond_length(1,2,4)=2*radius(4)+linker(4)
      intrabond_length(2,3,4)=2*radius(4)+linker(4)
      intrabond_length(3,4,4)=2*radius(4)+linker(4)
	  
c>>  customize PPI matrix


      complex_num=0
      do i=1,npart_max*type_num
         complex_nb_num(i)=0
         do j=1,complex_nb_ctg_num
            complex_nb_ctg(i,j)=0
            complex_nb_ctg_idx(i,j)=0
         enddo
      enddo

      do i=1,type_num
         do ii=1,intranum_max
            do j=1,type_num
               do jj=1,intranum_max
                  reaction_matrix(ii,i,jj,j)=0
                  Ass_Rate(ii,i,jj,j)=0
                  Diss_Rate(ii,i,jj,j)=0
               enddo
            enddo
         enddo
      enddo



      reaction_matrix(1,1,1,2)=1
      reaction_matrix(2,1,1,2)=1

      reaction_matrix(3,2,3,2)=1

      reaction_matrix(4,2,2,3)=1

      reaction_matrix(1,3,1,4)=1
      reaction_matrix(1,3,2,4)=1
      reaction_matrix(1,3,3,4)=1
      reaction_matrix(1,3,4,4)=1

      reaction_matrix(3,3,1,4)=1
      reaction_matrix(3,3,2,4)=1
      reaction_matrix(3,3,3,4)=1
      reaction_matrix(3,3,4,4)=1	  


      Ass_Rate(1,1,1,2)=50.0
      Ass_Rate(2,1,1,2)=50.0
      Ass_Rate(3,2,3,2)=50.0
      Ass_Rate(4,2,2,3)=50.0
      Ass_Rate(1,3,1,4)=50.0
      Ass_Rate(1,3,2,4)=50.0
      Ass_Rate(1,3,3,4)=50.0
      Ass_Rate(1,3,4,4)=50.0
      Ass_Rate(3,3,1,4)=50.0
      Ass_Rate(3,3,2,4)=50.0
      Ass_Rate(3,3,3,4)=50.0
      Ass_Rate(3,3,4,4)=50.0	

c>>>>>>>>>   Ligand-Receptor
      Diss_Rate(1,1,1,2)=0.00005
      Diss_Rate(2,1,1,2)=0.00005
c>>>>>>>>>   Receptor cis
      Diss_Rate(3,2,3,2)=50.0
c>>>>>>>>>   Receptor-Scaffold
      Diss_Rate(4,2,2,3)=0.0001
c>>>>>>>>>   Scaffold-Signaling
      Diss_Rate(1,3,1,4)=0.0001
      Diss_Rate(1,3,2,4)=0.0001
      Diss_Rate(1,3,3,4)=0.0001
      Diss_Rate(1,3,4,4)=0.0001
      Diss_Rate(3,3,1,4)=0.0001
      Diss_Rate(3,3,2,4)=0.0001
      Diss_Rate(3,3,3,4)=0.0001
      Diss_Rate(3,3,4,4)=0.0001


      do j=1,npart(2)
         C2_status(j)=0  ! the 4th domain in 'S' can attache to membrane
      enddo


c      open (unit=10,file=
c     &     'RTKLD_trj_10042022_001.pdb',
c     &     status='unknown',access='append')
c      do j=1,npart(1)
c         write(10,2100) 'ATOM  ',j,' CA ','ALA', 
c     &        'A',j,x(1,j,1),y(1,j,1),z(1,j,1)
c         write(10,2100) 'ATOM  ',j,' CA ','ILE', 
c     &        'A',j,x(2,j,1),y(2,j,1),z(2,j,1)
c      enddo
c      write(10,2102) 'TER'
c      do j=1,npart(2)
c         write(10,2100) 'ATOM  ',j,' CA ','VAL', 
c     &        'B',j,x(1,j,2),y(1,j,2),z(1,j,2)
c         write(10,2100) 'ATOM  ',j,' CA ','LYS', 
c     &        'B',j,x(2,j,2),y(2,j,2),z(2,j,2)
c         write(10,2100) 'ATOM  ',j,' CA ','ARG', 
c     &        'B',j,x(3,j,2),y(3,j,2),z(3,j,2)
c         write(10,2100) 'ATOM  ',j,' CA ','GLU', 
c     &        'B',j,x(4,j,2),y(4,j,2),z(4,j,2)	 
c      enddo
c      write(10,2102) 'TER'
c      do j=1,npart(3)
c         write(10,2100) 'ATOM  ',j,' CA ','LEU', 
c     &        'C',j,x(1,j,3),y(1,j,3),z(1,j,3)
c         write(10,2100) 'ATOM  ',j,' CA ','ASP', 
c     &        'C',j,x(2,j,3),y(2,j,3),z(2,j,3)
c         write(10,2100) 'ATOM  ',j,' CA ','ASN', 
c     &        'C',j,x(3,j,3),y(3,j,3),z(3,j,3)
c         write(10,2100) 'ATOM  ',j,' CA ','GLN', 
c     &        'C',j,x(4,j,3),y(4,j,3),z(4,j,3)	 
c      enddo
c      write(10,2102) 'TER'
c      do j=1,npart(4)
c         write(10,2100) 'ATOM  ',j,' CA ','GLY', 
c     &        'D',j,x(1,j,4),y(1,j,4),z(1,j,4)
c         write(10,2100) 'ATOM  ',j,' CA ','PRO', 
c     &        'D',j,x(2,j,4),y(2,j,4),z(2,j,4)
c         write(10,2100) 'ATOM  ',j,' CA ','TYR', 
c     &        'D',j,x(3,j,4),y(3,j,4),z(3,j,4)
c         write(10,2100) 'ATOM  ',j,' CA ','HIS', 
c     &        'D',j,x(4,j,4),y(4,j,4),z(4,j,4)	 
c      enddo
c      write(10,2102) 'TER'		
c      write(10,2102) 'END'
c      close(10)
      
     

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c>>>        main simulation of Langevin Dynamics
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      do itime=1,nsimu
c         print*,itime
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         do i=1,type_num
            do molecule=1,npart(i)
               do j=1,intra_num(i)
                  status_new(j,molecule,i)=status(j,molecule,i)
               enddo
            enddo
         enddo
         complex_num_new=complex_num
         do i=1,complex_num
            complex_nb_num_new(i)=complex_nb_num(i)
            do j=1,complex_nb_ctg_num
               complex_nb_ctg_new(i,j)=complex_nb_ctg(i,j)
               complex_nb_ctg_idx_new(i,j)=complex_nb_ctg_idx(i,j)
               complex_nb_ctg_plm_idx_new(i,j)=
     &              complex_nb_ctg_plm_idx(i,j)
            enddo
         enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   perform the diffusions for molecules
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Force Calculation
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c>  calculate the elastic forces between points

         do i=1,type_num
            do molecule=1,npart(i)
               do j=1,intra_num(i)
                  f_rp_x(j,molecule,i)=0
                  f_rp_y(j,molecule,i)=0
                  f_rp_z(j,molecule,i)=0
               enddo
            enddo
         enddo


         do i=1,type_num
            do j=i,type_num
               
               do molecule=1,npart(i)
                  do neighbor=1,npart(j)

                     if((i.eq.j).and.(molecule.eq.neighbor))then

                        do ii=1,intra_num(i)-1
                           do jj=ii+1,intra_num(j)

                              dx=x(jj,neighbor,j)-
     &                             x(ii,molecule,i)
                              dx=dx-cell_range_x*anint(dx/cell_range_x)
                              
                              dy=y(jj,neighbor,j)-
     &                             y(ii,molecule,i)
                              dy=dy-cell_range_y*anint(dy/cell_range_y)
                              
                              dz=z(jj,neighbor,j)-
     &                             z(ii,molecule,i)
                              dz=dz-cell_range_z*anint(dz/cell_range_z)
                              
                              dij=sqrt(dx**2+dy**2+dz**2)
                              
                              if(dij.lt.(radius(i)+radius(j)
     &                             +reaction_range))then
                                 
                                 f_rp_x(ii,molecule,i)=
     &                                f_rp_x(ii,molecule,i)
     &                                +spring_length_rp*
     &                                (1-(radius(i)+radius(j)
     &                                +reaction_range)/dij)*dx
                                 
                                 f_rp_y(ii,molecule,i)=
     &                                f_rp_y(ii,molecule,i)
     &                                +spring_length_rp*
     &                                (1-(radius(i)+radius(j)
     &                                +reaction_range)/dij)*dy
                                 
                                 f_rp_z(ii,molecule,i)=
     &                                f_rp_z(ii,molecule,i)
     &                                +spring_length_rp*
     &                                (1-(radius(i)+radius(j)
     &                                +reaction_range)/dij)*dz
                                 
                                 f_rp_x(jj,neighbor,j)=
     &                                f_rp_x(jj,neighbor,j)
     &                                -spring_length_rp*
     &                                (1-(radius(i)+radius(j)
     &                                +reaction_range)/dij)*dx
                                 
                                 f_rp_y(jj,neighbor,j)=
     &                                f_rp_y(jj,neighbor,j)
     &                                -spring_length_rp*
     &                                (1-(radius(i)+radius(j)
     &                                +reaction_range)/dij)*dy
                                 
                                 f_rp_z(jj,neighbor,j)=
     &                                f_rp_z(jj,neighbor,j)
     &                                -spring_length_rp*
     &                                (1-(radius(i)+radius(j)
     &                                +reaction_range)/dij)*dz

                              endif
                              
                           enddo
                        enddo

                     else

                        do ii=1,intra_num(i)
                           do jj=1,intra_num(j)

                              dx=x(jj,neighbor,j)-
     &                             x(ii,molecule,i)
                              dx=dx-cell_range_x*anint(dx/cell_range_x)
                              
                              dy=y(jj,neighbor,j)-
     &                             y(ii,molecule,i)
                              dy=dy-cell_range_y*anint(dy/cell_range_y)
                              
                              dz=z(jj,neighbor,j)-
     &                             z(ii,molecule,i)
                              dz=dz-cell_range_z*anint(dz/cell_range_z)
                              
                              dij=sqrt(dx**2+dy**2+dz**2)
                              
                              if(dij.lt.(radius(i)+radius(j)
     &                             +reaction_range))then
                                 
                                 f_rp_x(ii,molecule,i)=
     &                                f_rp_x(ii,molecule,i)
     &                                +spring_length_rp*
     &                                (1-(radius(i)+radius(j)
     &                                +reaction_range)/dij)*dx
                                 
                                 f_rp_y(ii,molecule,i)=
     &                                f_rp_y(ii,molecule,i)
     &                                +spring_length_rp*
     &                                (1-(radius(i)+radius(j)
     &                                +reaction_range)/dij)*dy
                                 
                                 f_rp_z(ii,molecule,i)=
     &                                f_rp_z(ii,molecule,i)
     &                                +spring_length_rp*
     &                                (1-(radius(i)+radius(j)
     &                                +reaction_range)/dij)*dz
                                 
                                 f_rp_x(jj,neighbor,j)=
     &                                f_rp_x(jj,neighbor,j)
     &                                -spring_length_rp*
     &                                (1-(radius(i)+radius(j)
     &                                +reaction_range)/dij)*dx
                                 
                                 f_rp_y(jj,neighbor,j)=
     &                                f_rp_y(jj,neighbor,j)
     &                                -spring_length_rp*
     &                                (1-(radius(i)+radius(j)
     &                                +reaction_range)/dij)*dy
                                 
                                 f_rp_z(jj,neighbor,j)=
     &                                f_rp_z(jj,neighbor,j)
     &                                -spring_length_rp*
     &                                (1-(radius(i)+radius(j)
     &                                +reaction_range)/dij)*dz

                              endif

                           enddo
                        enddo

                     endif

                  enddo
               enddo

            enddo
         enddo

c>  calculate intrmolecular forces within each molecule

         do i=1,type_num
            do molecule=1,npart(i)
               do j=1,intra_num(i)
                  f_plm_x(j,molecule,i)=0
                  f_plm_y(j,molecule,i)=0
                  f_plm_z(j,molecule,i)=0
               enddo
            enddo
         enddo

         do i=1,type_num
            do molecule=1,npart(i)

               do ii=1,intra_num(i)-1
                  do jj=ii+1,intra_num(i)

                     if(intra_geo_matrix(ii,jj,i).eq.1)then

                        dx=x(jj,molecule,i)-
     &                       x(ii,molecule,i)
                        dx=dx-cell_range_x*anint(dx/cell_range_x)
                        
                        dy=y(jj,molecule,i)-
     &                       y(ii,molecule,i)
                        dy=dy-cell_range_y*anint(dy/cell_range_y)
                        
                        dz=z(jj,molecule,i)-
     &                       z(ii,molecule,i)
                        dz=dz-cell_range_z*anint(dz/cell_range_z)
                        
                        dij=sqrt(dx**2+dy**2+dz**2)

                        f_plm_x(ii,molecule,i)=
     &                       f_plm_x(ii,molecule,i)
     &                       +spring_length_plm(i)*
     &                       (1-intrabond_length(ii,jj,i)/dij)*dx
                        
                        f_plm_y(ii,molecule,i)=
     &                       f_plm_y(ii,molecule,i)
     &                       +spring_length_plm(i)*
     &                       (1-intrabond_length(ii,jj,i)/dij)*dy
                        
                        f_plm_z(ii,molecule,i)=
     &                       f_plm_z(ii,molecule,i)
     &                       +spring_length_plm(i)*
     &                       (1-intrabond_length(ii,jj,i)/dij)*dz
                        
                        f_plm_x(jj,molecule,i)=
     &                       f_plm_x(jj,molecule,i)
     &                       -spring_length_plm(i)*
     &                       (1-intrabond_length(ii,jj,i)/dij)*dx
                        
                        f_plm_y(jj,molecule,i)=
     &                       f_plm_y(jj,molecule,i)
     &                       -spring_length_plm(i)*
     &                       (1-intrabond_length(ii,jj,i)/dij)*dy
                        
                        f_plm_z(jj,molecule,i)=
     &                       f_plm_z(jj,molecule,i)
     &                       -spring_length_plm(i)*
     &                       (1-intrabond_length(ii,jj,i)/dij)*dz
                        
                     endif

                  enddo
               enddo

            enddo
         enddo

c>  calculate bond forces between subunits of a dimer

         do i=1,type_num
            do molecule=1,npart(i)
               do j=1,intra_num(i)
                  f_bd_x(j,molecule,i)=0
                  f_bd_y(j,molecule,i)=0
                  f_bd_z(j,molecule,i)=0
               enddo
            enddo
         enddo

         do i=1,type_num-1
            do j=1,type_num
               
               do molecule=1,npart(i)
                  do neighbor=1,npart(j)

                     do ii=1,intra_num(i)
                        do jj=1,intra_num(j)

                           do k=1,complex_num_new
                              if(((complex_nb_ctg_new(k,1).eq.i).AND.
     &                             (complex_nb_ctg_idx_new(k,1).eq.
     &                             molecule).AND.
     &                             (complex_nb_ctg_plm_idx_new(k,1).eq.
     &                             ii))
     &                             .AND.
     &                             ((complex_nb_ctg_new(k,2).eq.j).AND.
     &                             (complex_nb_ctg_idx_new(k,2).eq.
     &                             neighbor).AND.
     &                             (complex_nb_ctg_plm_idx_new(k,2).eq.
     &                             jj))
     &                             )then
                                 
                                 dx=x(jj,neighbor,j)-
     &                                x(ii,molecule,i)
                                 dx=dx-cell_range_x*
     &                                anint(dx/cell_range_x)
                                 
                                 dy=y(jj,neighbor,j)-
     &                                y(ii,molecule,i)
                                 dy=dy-cell_range_y*
     &                                anint(dy/cell_range_y)
                                 
                                 dz=z(jj,neighbor,j)-
     &                                z(ii,molecule,i)
                                 dz=dz-cell_range_z*
     &                                anint(dz/cell_range_z)
                                 
                                 dij=sqrt(dx**2+dy**2+dz**2)
                                 
                                 f_bd_x(ii,molecule,i)=
     &                                f_bd_x(ii,molecule,i)
     &                                +spring_length_bd*
     &                                (1-bond_dist_cutoff/dij)*dx
                                 
                                 f_bd_y(ii,molecule,i)=
     &                                f_bd_y(ii,molecule,i)
     &                                +spring_length_bd*
     &                                (1-bond_dist_cutoff/dij)*dy
                                 
                                 f_bd_z(ii,molecule,i)=
     &                                f_bd_z(ii,molecule,i)
     &                                +spring_length_bd*
     &                                (1-bond_dist_cutoff/dij)*dz
                                 
                                 f_bd_x(jj,neighbor,j)=
     &                                f_bd_x(jj,neighbor,j)
     &                                -spring_length_bd*
     &                                (1-bond_dist_cutoff/dij)*dx
                                 
                                 f_bd_y(jj,neighbor,j)=
     &                                f_bd_y(jj,neighbor,j)
     &                                -spring_length_bd*
     &                                (1-bond_dist_cutoff/dij)*dy
                                 
                                 f_bd_z(jj,neighbor,j)=
     &                                f_bd_z(jj,neighbor,j)
     &                                -spring_length_bd*
     &                                (1-bond_dist_cutoff/dij)*dz


                              endif
                           enddo

                        enddo
                     enddo

                  enddo
               enddo

            enddo
         enddo

c>> calculate the stochastic force of each points

c>>>  1) stochastic force for the translational motion of the entire molecules
 
         do i=1,type_num
            do molecule=1,npart(i)
               dir_x=2*rand3(r3)-1
               dir_y=2*rand3(r3)-1
               dir_z=2*rand3(r3)-1
               amp_dir=sqrt(dir_x**2+dir_y**2+dir_z**2)
               scaled_dir_x=dir_x/amp_dir
               scaled_dir_y=dir_y/amp_dir
               scaled_dir_z=dir_z/amp_dir
               do ii=1,intra_num(i)
                  f_rdm_x(ii,molecule,i)=2*sqrt(diffu_const(i))*
     &                 cm_frict_const*scaled_dir_x
                  f_rdm_y(ii,molecule,i)=2*sqrt(diffu_const(i))*
     &                 cm_frict_const*scaled_dir_y
                  f_rdm_z(ii,molecule,i)=2*sqrt(diffu_const(i))*
     &                 cm_frict_const*scaled_dir_z
               enddo
            enddo
         enddo


c>>>  2) stochastic force for the rotational motion of the entire molecules
 
         do i=1,type_num
            do molecule=1,npart(i)

               theta=(2*rand3(r3)-1)*pai
               phi=(2*rand3(r3)-1)*pai
               psai=(2*rand3(r3)-1)*pai
               
               t(1,1)=cos(psai)*cos(phi)-cos(theta)*sin(phi)*sin(psai)
               t(1,2)=-sin(psai)*cos(phi)-cos(theta)*sin(phi)*cos(psai)
               t(1,3)=sin(theta)*sin(phi)
               
               t(2,1)=cos(psai)*sin(phi)+cos(theta)*cos(phi)*sin(psai)
               t(2,2)=-sin(psai)*sin(phi)+cos(theta)*cos(phi)*cos(psai)
               t(2,3)=-sin(theta)*cos(phi)
               
               t(3,1)=sin(psai)*sin(theta)
               t(3,2)=cos(psai)*sin(theta)
               t(3,3)=cos(theta)

               cm_x=0
               cm_y=0
               cm_z=0
               do j=1,intra_num(i)
                  if(x(j,molecule,i)-x(1,molecule,i)
     &                 .gt.cell_range_x*0.5)then
                     temp_x(j,molecule,i)=x(j,molecule,i)-cell_range_x
                  elseif(x(j,molecule,i)-x(1,molecule,i)
     &                 .lt.-cell_range_x*0.5)then
                     temp_x(j,molecule,i)=x(j,molecule,i)+cell_range_x
                  else
                     temp_x(j,molecule,i)=x(j,molecule,i)
                  endif
                  if(y(j,molecule,i)-y(1,molecule,i)
     &                 .gt.cell_range_y*0.5)then
                     temp_y(j,molecule,i)=y(j,molecule,i)-cell_range_y
                  elseif(y(j,molecule,i)-y(1,molecule,i)
     &                 .lt.-cell_range_y*0.5)then
                     temp_y(j,molecule,i)=y(j,molecule,i)+cell_range_y
                  else
                     temp_y(j,molecule,i)=y(j,molecule,i)
                  endif
c                  if(z(j,molecule,i)-z(1,molecule,i)
c     &                 .gt.cell_range_z*0.5)then
c                     temp_z(j,molecule,i)=z(j,molecule,i)-cell_range_z
c                  elseif(z(j,molecule,i)-z(1,molecule,i)
c     &                 .lt.-cell_range_z*0.5)then
c                     temp_z(j,molecule,i)=z(j,molecule,i)+cell_range_z
c                  else
                     temp_z(j,molecule,i)=z(j,molecule,i)
c                  endif
               enddo
               do j=1,intra_num(i)
                  cm_x=cm_x+temp_x(j,molecule,i)
                  cm_y=cm_y+temp_y(j,molecule,i)
                  cm_z=cm_z+temp_z(j,molecule,i)
               enddo
               cm_x=cm_x/real(intra_num(i))
               cm_y=cm_y/real(intra_num(i))
               cm_z=cm_z/real(intra_num(i))

               amp_rdmrt=2*sqrt(rtdiffu_const(i))*cm_frict_const

               do j=1,intra_num(i)
                  f_rdmrt_x(j,molecule,i)=amp_rdmrt*(
     &                 t(1,1)*(temp_x(j,molecule,i)-cm_x)+
     &                 t(1,2)*(temp_y(j,molecule,i)-cm_y)+
     &                 t(1,3)*(temp_z(j,molecule,i)-cm_z))
                  f_rdmrt_y(j,molecule,i)=amp_rdmrt*(
     &                 t(2,1)*(temp_x(j,molecule,i)-cm_x)+
     &                 t(2,2)*(temp_y(j,molecule,i)-cm_y)+
     &                 t(2,3)*(temp_z(j,molecule,i)-cm_z))
                  f_rdmrt_z(j,molecule,i)=amp_rdmrt*(
     &                 t(3,1)*(temp_x(j,molecule,i)-cm_x)+
     &                 t(3,2)*(temp_y(j,molecule,i)-cm_y)+
     &                 t(3,3)*(temp_z(j,molecule,i)-cm_z))
               enddo

            enddo
         enddo

c>>>  3) stochastic force for the internal fluctuations of the entire molecules

         do i=1,type_num
            do molecule=1,npart(i)
               do ii=1,intra_num(i)
                  dir_x=2*rand4(r4)-1
                  dir_y=2*rand4(r4)-1
                  dir_z=2*rand5(r5)-1
                  amp_dir=sqrt(dir_x**2+dir_y**2+dir_z**2)
                  scaled_dir_x=dir_x/amp_dir
                  scaled_dir_y=dir_y/amp_dir
                  scaled_dir_z=dir_z/amp_dir                  
                  f_rdmfl_x(ii,molecule,i)=2*sqrt(fluct_const(i))*
     &                 cm_frict_const*scaled_dir_x
                  f_rdmfl_y(ii,molecule,i)=2*sqrt(fluct_const(i))*
     &                 cm_frict_const*scaled_dir_y
                  f_rdmfl_z(ii,molecule,i)=2*sqrt(fluct_const(i))*
     &                 cm_frict_const*scaled_dir_z
               enddo
            enddo
         enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>> calculate the velocity of each particle
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         do i=1,type_num
            do molecule=1,npart(i)
               do ii=1,intra_num(i)
                  vx_old(ii,molecule,i)=0
                  vy_old(ii,molecule,i)=0
                  vz_old(ii,molecule,i)=0
                  do j=1,type_num
                     do neighbor=1,npart(j)
                        do jj=1,intra_num(j)
                           ax(ii,molecule,i,jj,neighbor,j)=0
                           ay(ii,molecule,i,jj,neighbor,j)=0
                           az(ii,molecule,i,jj,neighbor,j)=0
                        enddo
                     enddo
                  enddo
                  bx(ii,molecule,i)=0
                  by(ii,molecule,i)=0
                  bz(ii,molecule,i)=0
               enddo
            enddo
         enddo
 
         do i=1,type_num
            do molecule=1,npart(i)
               do ii=1,intra_num(i)
                  ax(ii,molecule,i,ii,molecule,i)=
     &                 ax(ii,molecule,i,ii,molecule,i)+cm_frict_const
                  ay(ii,molecule,i,ii,molecule,i)=
     &                 ay(ii,molecule,i,ii,molecule,i)+cm_frict_const
                  az(ii,molecule,i,ii,molecule,i)=
     &                 az(ii,molecule,i,ii,molecule,i)+cm_frict_const
                  do j=1,type_num
                     do neighbor=1,npart(j)
                        do jj=1,intra_num(j)
                           if((i.ne.j).OR.((i.eq.j)
     &                          .AND.(molecule.ne.neighbor)).OR.
     &                          ((i.eq.j).AND.(molecule.ne.neighbor)
     &                          .AND.(ii.ne.jj)))then
                              dij=sqrt((x(ii,molecule,i)-
     &                             x(jj,neighbor,j))**2+
     &                             (y(ii,molecule,i)-
     &                             y(jj,neighbor,j))**2+
     &                             (z(ii,molecule,i)-
     &                             z(jj,neighbor,j))**2)
                              if(dij.lt.(radius(i)+radius(j)))then
                                 ax(ii,molecule,i,ii,molecule,i)=
     &                                ax(ii,molecule,i,ii,molecule,i)+
     &                                cc_frict_const
                                 ax(ii,molecule,i,jj,neighbor,j)=
     &                                ax(ii,molecule,i,jj,neighbor,j)+
     &                                cc_frict_const
                                 ay(ii,molecule,i,ii,molecule,i)=
     &                                ay(ii,molecule,i,ii,molecule,i)+
     &                                cc_frict_const
                                 ay(ii,molecule,i,jj,neighbor,j)=
     &                                ay(ii,molecule,i,jj,neighbor,j)+
     &                                cc_frict_const
                                 az(ii,molecule,i,ii,molecule,i)=
     &                                az(ii,molecule,i,ii,molecule,i)+
     &                                cc_frict_const
                                 az(ii,molecule,i,jj,neighbor,j)=
     &                                az(ii,molecule,i,jj,neighbor,j)+
     &                                cc_frict_const
                              endif
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo

         do i=1,type_num
            do molecule=1,npart(i)
               do ii=1,intra_num(i)
                  bx(ii,molecule,i)=
     &                 f_rp_x(ii,molecule,i)
     &                 +f_plm_x(ii,molecule,i)
     &                 +f_bd_x(ii,molecule,i)
     &                 +f_rdm_x(ii,molecule,i)
     &                 +f_rdmrt_x(ii,molecule,i)
     &                 +f_rdmfl_x(ii,molecule,i)
                  by(ii,molecule,i)=
     &                 f_rp_y(ii,molecule,i)
     &                 +f_plm_y(ii,molecule,i)
     &                 +f_bd_y(ii,molecule,i)
     &                 +f_rdm_y(ii,molecule,i)
     &                 +f_rdmrt_y(ii,molecule,i)
     &                 +f_rdmfl_y(ii,molecule,i)
                  bz(ii,molecule,i)=
     &                 f_rp_z(ii,molecule,i)
     &                 +f_plm_z(ii,molecule,i)
     &                 +f_bd_z(ii,molecule,i)
     &                 +f_rdm_z(ii,molecule,i)
     &                 +f_rdmrt_z(ii,molecule,i)
     &                 +f_rdmfl_z(ii,molecule,i)
               enddo
            enddo
         enddo

         do k=1,niter

            do i=1,type_num
               do molecule=1,npart(i)
                  do ii=1,intra_num(i)
                     vx_new(ii,molecule,i)=bx(ii,molecule,i)
                     vy_new(ii,molecule,i)=by(ii,molecule,i)
                     vz_new(ii,molecule,i)=bz(ii,molecule,i)
                     do j=1,type_num
                        do neighbor=1,npart(j)
                           do jj=1,intra_num(j)
                              if((i.ne.j).OR.((i.eq.j)
     &                             .AND.(molecule.ne.neighbor)).OR.
     &                             ((i.eq.j).AND.(molecule.ne.neighbor)
     &                             .AND.(ii.ne.jj)))then
                                 vx_new(ii,molecule,i)=
     &                                vx_new(ii,molecule,i)+
     &                                ax(ii,molecule,i,jj,neighbor,j)*
     &                                vx_old(jj,neighbor,j)
                                 vy_new(ii,molecule,i)=
     &                                vy_new(ii,molecule,i)+
     &                                ay(ii,molecule,i,jj,neighbor,j)*
     &                                vy_old(jj,neighbor,j)
                                 vz_new(ii,molecule,i)=
     &                                vz_new(ii,molecule,i)+
     &                                az(ii,molecule,i,jj,neighbor,j)*
     &                                vz_old(jj,neighbor,j)
                              endif
                           enddo
                        enddo
                     enddo
                     vx_new(ii,molecule,i)=vx_new(ii,molecule,i)
     &                    /ax(ii,molecule,i,ii,molecule,i)
                     vy_new(ii,molecule,i)=vy_new(ii,molecule,i)
     &                    /ay(ii,molecule,i,ii,molecule,i)
                     vz_new(ii,molecule,i)=vz_new(ii,molecule,i)
     &                    /az(ii,molecule,i,ii,molecule,i)
                  enddo
               enddo
            enddo

            RMSD=0
            npart_tot=0
            do i=1,type_num
               do molecule=1,npart(i)
                  RMSD=RMSD+(vx_new(ii,molecule,i)-
     &                 vx_old(ii,molecule,i))**2+
     &                 (vy_new(ii,molecule,i)-
     &                 vy_old(ii,molecule,i))**2+
     &                 (vz_new(ii,molecule,i)-
     &                 vz_old(ii,molecule,i))**2
               enddo
               npart_tot=npart_tot+npart(i)
            enddo

            RMSD=sqrt(RMSD/real(npart_tot))

            if(RMSD.lt.iter_tol)then
               goto 100
            endif
            do i=1,type_num
               do molecule=1,npart(i)
                  do ii=1,intra_num(i)
                     vx_old(ii,molecule,i)=vx_new(ii,molecule,i)
                     vy_old(ii,molecule,i)=vy_new(ii,molecule,i)
                     vz_old(ii,molecule,i)=vz_new(ii,molecule,i)
                  enddo
               enddo
            enddo

         enddo
         
 100     continue

         do i=1,type_num
            do molecule=1,npart(i)
               do ii=1,intra_num(i)
                  vx(ii,molecule,i)=vx_new(ii,molecule,i)
                  vy(ii,molecule,i)=vy_new(ii,molecule,i)
                  vz(ii,molecule,i)=vz_new(ii,molecule,i)
               enddo
            enddo
         enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>>>  update the position of each cell
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         do i=1,type_num
            do j=1,npart(i)
               do ii=1,intra_num(i)
                  x(ii,j,i)=x(ii,j,i)+vx(ii,j,i)*dt
                  y(ii,j,i)=y(ii,j,i)+vy(ii,j,i)*dt
                  if((i.eq.2).and.(ii.eq.3))then ! the transmembrane domain of R has to be confined in the plasma membrane
                     z(ii,j,i)=0  ! the z coordinate of plasma membrane is zero
                  elseif((i.eq.3).and.(ii.eq.4))then ! the C2 domain of S has the probability to attach below the membrane surface
                     if(C2_status(j).eq.1)then
                        z(ii,j,i)=0	
                     elseif(C2_status(j).eq.0)then	
                        z(ii,j,i)=z(ii,j,i)+vz(ii,j,i)*dt
                     endif									 
                  else
                     z(ii,j,i)=z(ii,j,i)+vz(ii,j,i)*dt
                  endif
cc>>>>>>>>>>>>>>>>>   boundary condition: periodic boundary along x and y direction, reflection along z direction above and below the membrane surface, as well as at the top and bottom of the simulation box 

c>>>    x and y direction follow period boundary condition
                  x(ii,j,i)=x(ii,j,i)-
     &                 cell_range_x*anint(x(ii,j,i)/cell_range_x)
                  y(ii,j,i)=y(ii,j,i)-
     &                 cell_range_y*anint(y(ii,j,i)/cell_range_y)

c>>>  along z direction can not move outside of the simulation box at the top and bottom

                  if(z(ii,j,i).gt.cell_range_z/2)then
                     z(ii,j,i)=z(ii,j,i)-abs(vz(ii,j,i)*dt*2)
                  elseif(z(ii,j,i).lt.-cell_range_z/2)then
                     z(ii,j,i)=z(ii,j,i)+abs(vz(ii,j,i)*dt*2)
                  endif

c>>>  L can not move below the membrane, S and E can not move above the membrane
c>>>  the C2 domain of S also can attach to the membrane 

                  if(i.eq.1)then
                     if(z(ii,j,i).lt.0)then
                        z(ii,j,i)=z(ii,j,i)+abs(vz(ii,j,i)*dt*2)
                     endif
                  elseif(i.eq.4)then
                     if(z(ii,j,i).gt.0)then
                        z(ii,j,i)=z(ii,j,i)-abs(vz(ii,j,i)*dt*2)
                     endif
                  elseif(i.eq.3)then
                     if(ii.ne.4)then
                        if(z(ii,j,i).gt.0)then
                           z(ii,j,i)=z(ii,j,i)-abs(vz(ii,j,i)*dt*2)
                        endif
                     elseif(ii.eq.4)then
                        if(z(ii,j,i).gt.0)then
                           temp=rand5(r5)
                           if((C2_status(j).eq.0).AND.
     &                          (temp.le.Attach_Prob))then
                              z(ii,j,i)=cell_range_z/2.0
                              C2_status(j)=1
                           else
                              z(ii,j,i)=z(ii,j,i)-abs(vz(ii,j,i)*dt*2)
                           endif	
                        endif
                     endif
                  endif

               enddo
            enddo
         enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   perform the reaction for molecules
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  1) two molecules associate into a complex
         do i=1,type_num-1
            do j=1,type_num
               do ii=1,intra_num(i)
                  do jj=1,intra_num(j)
               
                     if(reaction_matrix(ii,i,jj,j).eq.1)then

                        do molecule=1,npart(i)
                           do neighbor=1,npart(j)
                        
                              if(((i.eq.j).and.(ii.eq.jj).and.
     &                             (molecule.ne.neighbor))
     &                             .or.(i.ne.j))then
	 
                                 if((status_new(ii,molecule,i).eq.0)
     &                                .AND.
     &                                (status_new(jj,neighbor,j)
     &                                .eq.0)
     &                                )then

c>>>>>>>>>>>>>>>>>>>>  set up the rules for different specific interaction as phosphorylate regulation

                                    Ass_Flag=0

cccccccccccccccccc   the first rule is that only after R form L induced dimer, which will trigger the phosphorylation of intracellular domain in R to interact with S


                                    if((i.eq.2).AND.(j.eq.3)
     &                                   .AND.(ii.eq.4).AND.(jj.eq.2)
     &                                   )then ! the interaction between R and S
                                       
                                       if(status_new(1,molecule,2)
     &                                      .eq.1)then ! the R has to be ligand bound
                                          
                                          index_mol=0
                                          
                                          do k=1,complex_num_new ! found out which ligand bind to the receptor
                                             if((complex_nb_ctg_new
     &                                            (k,2).eq.2).AND.
     &                                            (complex_nb_ctg
     &                                            _idx_new(k,2).eq.
     &                                            molecule).AND.
     &                                            (complex_nb_ctg
     &                                            _plm_idx_new(k,2)
     &                                            .eq.1))then
                                                
                                                index_mol=
     &                                               complex_nb_ctg
     &                                               _idx_new(k,1)
                                                
                                             endif
                                          enddo
     &                                            
                                          if((index_mol.ne.0).AND.
     &                                         (status_new(1,
     &                                         index_mol,1)
     &                                         .eq.1).AND.
     &                                         (status_new(2,
     &                                         index_mol,1)
     &                                         .eq.1))then ! both subunits of the L dimer need to bind to R
                                             
                                             Ass_Flag=1
                                             
                                          endif
                                          
                                       endif
                                          

cccccccccccccccccccc  the second rule is that only after the S bind to R, which trigger the phosphorylation of S to further recruit E 

                                    elseif((i.eq.3).AND.(j.eq.4)
     &                                      )then ! the interaction between S and E
                                       
                                       if(status_new(2,molecule,3)
     &                                      .eq.1)then ! the S has to be bound to R
                                          
                                          Ass_Flag=1
                                          
                                       endif
                               

cccccccccccccccccccc  the third rule is that only after the S bind to R, which trigger the cis association between two R 

                                    elseif((i.eq.2).AND.(j.eq.2)
     &                                   .AND.(ii.eq.3).AND.(jj.eq.3)
     &                                   )then ! the lateral interaction between two R

                                       if((status_new(1,molecule,2)
     &                                      .eq.1).AND.
     &                                      (status_new(1,neighbor,2)
     &                                      .eq.1))then ! the R has to be ligand bound

                                          index_mol=0
                                          
                                          do k=1,complex_num_new ! found out which ligand bind to the first receptor 
                                             if((complex_nb_ctg_new
     &                                            (k,2).eq.2).AND.
     &                                            (complex_nb_ctg
     &                                            _idx_new(k,2).eq.
     &                                            molecule).AND.
     &                                            (complex_nb_ctg
     &                                            _plm_idx_new(k,2)
     &                                            .eq.1))then
                                                
                                                index_mol=
     &                                               complex_nb_ctg
     &                                               _idx_new(k,1)
                                                
                                             endif
                                          enddo


                                          index_mol2=0
                                          
                                          do k=1,complex_num_new ! found out which ligand bind to the second receptor 
                                             if((complex_nb_ctg_new
     &                                            (k,2).eq.2).AND.
     &                                            (complex_nb_ctg
     &                                            _idx_new(k,2).eq.
     &                                            neighbor).AND.
     &                                            (complex_nb_ctg
     &                                            _plm_idx_new(k,2)
     &                                            .eq.1))then
                                                
                                                index_mol2=
     &                                               complex_nb_ctg
     &                                               _idx_new(k,1)
                                                
                                             endif
                                          enddo         


                                          if((index_mol.ne.0)
     &                                         .AND.(index_mol2.ne.0)
     &                                         .AND.(index_mol.ne.
     &                                         index_mol2))then ! the two R need to bind to different L in order to form a cis-interaction
                                             
                                             Ass_Flag=1
                                             
                                          endif
                                          


                                       endif

cc>>>>>>>>>>>>>>>   other type of interaction

                                    elseif((i.eq.1).AND.(j.eq.2)
     &                                      )then ! the ligand receptor interaction

                                       Ass_Flag=1
                                      

                                    endif

                               
                                    
                                    if(Ass_Flag.eq.1)then ! the rules are statisfied

c>>>>>>>>>>>>>calculate distance of two reaction sites between RB_A and RB_B
                                       dist=sqrt((x(jj,neighbor,j)-
     &                                      x(ii,molecule,i))**2+
     &                                      (y(jj,neighbor,j)-
     &                                      y(ii,molecule,i))**2+
     &                                      (z(jj,neighbor,j)-
     &                                      z(ii,molecule,i))**2)
                                       
                                       if(dist.lt.bond_dist_cutoff
     &                                      +reaction_shell)then
                                          
                                          
                                          
                                          Prob_Ass=Ass_Rate(ii,i,jj,j)
     &                                         *dt ! ready to change
                                          prob=rand3(r3)
                                          if((prob.lt.Prob_Ass).AND.
     &                                         (prob.gt.0.000))then
                                             status_new(ii,molecule,i)=1
                                             status_new(jj,neighbor,j)=1
                                             complex_num_new=
     &                                            complex_num_new+1
                                             complex_nb_num_new
     &                                            (complex_num_new)=2
                                             complex_nb_ctg_new
     &                                            (complex_num_new,1)=i
                                             complex_nb_ctg_new
     &                                            (complex_num_new,2)=j
                                             complex_nb_ctg_idx_new
     &                                            (complex_num_new,1)=
     &                                            molecule
                                             complex_nb_ctg_idx_new
     &                                            (complex_num_new,2)=
     &                                            neighbor
                                             complex_nb_ctg_plm_idx_new
     &                                            (complex_num_new,1)=ii
                                             complex_nb_ctg_plm_idx_new
     &                                            (complex_num_new,2)=jj
                                             
                                             
                                          endif
                                       endif
                                       
                                    endif

                                   

                                 endif

                              endif
                              
                           enddo
                        enddo

                     endif
                  enddo
               enddo

            enddo
         enddo



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  2)  complex dissociate into monomers

         do i=1,type_num
            do molecule=1,npart(i)
               do ii=1,intra_num(i)
                  if(status_new(ii,molecule,i).eq.1)then
                     do j=1,complex_num_new
                        if(((complex_nb_ctg_new(j,1).eq.i).AND.
     &                       (complex_nb_ctg_idx_new(j,1).eq.molecule)
     &                       .AND.
     &                       (complex_nb_ctg_plm_idx_new(j,1).eq.ii))
     &                      .or.
     &                      ((complex_nb_ctg_new(j,2).eq.i).AND.
     &                       (complex_nb_ctg_idx_new(j,2).eq.molecule)
     &                       .AND.
     &                       (complex_nb_ctg_plm_idx_new(j,2).eq.ii))
     &                       )then
                           selected_complex=j
                           selected_ctg_A=complex_nb_ctg_new(j,1)
                           selected_ctg_B=complex_nb_ctg_new(j,2)
                           selected_idx_A=complex_nb_ctg_idx_new(j,1)
                           selected_idx_B=complex_nb_ctg_idx_new(j,2)
                           selected_plm_idx_A=
     &                          complex_nb_ctg_plm_idx_new(j,1)
                           selected_plm_idx_B=
     &                          complex_nb_ctg_plm_idx_new(j,2)
                        endif
                     enddo
                     Prob_Diss=Diss_Rate
     &                    (selected_plm_idx_A,selected_ctg_A,
     &                    selected_plm_idx_B,selected_ctg_B)*dt
                  
                     part1=int(rand3(r3)*10)
                     part2=int(rand5(r5)*10)
                     part3=int(rand5(r5)*10)
                     temp=rand3(r3)
                     part4=int(rand3(r3)*10)
                     part5=int(rand4(r4)*10)
                     do j=1,int(rand4(r4)*10)
                        temp=rand5(r5)
                        temp=rand5(r5)
                     enddo
                     part6=int(rand5(r5)*10)
                     part7=int(rand5(r5)*10)
                     part8=int(rand4(r4)*10)
                     part9=int(rand5(r5)*10)
                     
                     prob=real(part1)/10.0+real(part2)/100.0
     &                    +real(part3)/1000.0+real(part4)/10000.0
     &                    +real(part5)/100000.0
     &                    +real(part6)/1000000.0
     &                    +real(part7)/10000000.0
     &                    +real(part8)/100000000.0
     &                    +real(part9)/1000000000.0
                     if((prob.lt.Prob_Diss).AND.(prob.gt.0.000))then
                        status_new(selected_plm_idx_A,
     &                       selected_idx_A,selected_ctg_A)=0
                        status_new(selected_plm_idx_B,
     &                       selected_idx_B,selected_ctg_B)=0
                        do k=selected_complex+1,complex_num_new
                           complex_nb_num_new(k-1)=complex_nb_num_new(k)
                           complex_nb_ctg_idx_new(k-1,1)=
     &                          complex_nb_ctg_idx_new(k,1)
                           complex_nb_ctg_idx_new(k-1,2)=
     &                          complex_nb_ctg_idx_new(k,2)
                           complex_nb_ctg_new(k-1,1)=
     &                          complex_nb_ctg_new(k,1)
                           complex_nb_ctg_new(k-1,2)=
     &                          complex_nb_ctg_new(k,2)
                           complex_nb_ctg_plm_idx_new(k-1,1)=
     &                          complex_nb_ctg_plm_idx_new(k,1)
                           complex_nb_ctg_plm_idx_new(k-1,2)=
     &                          complex_nb_ctg_plm_idx_new(k,2)
                        enddo
                        complex_num_new=
     &                       complex_num_new-1
                     endif
                  endif
               enddo
            enddo
         enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  3)  C2 domain dettach from membrane

         do i=1,npart(2)
            if(C2_status(i).eq.1)then
               Prob_Dettach=Dettach_Rate*dt		
               part1=int(rand3(r3)*10)
               part2=int(rand5(r5)*10)
               part3=int(rand5(r5)*10)
               temp=rand3(r3)
               part4=int(rand3(r3)*10)
               part5=int(rand4(r4)*10)
               do j=1,int(rand4(r4)*10)
                  temp=rand5(r5)
                  temp=rand5(r5)
               enddo
               part6=int(rand5(r5)*10)
               part7=int(rand5(r5)*10)
               part8=int(rand4(r4)*10)
               part9=int(rand5(r5)*10)
               
               prob=real(part1)/10.0+real(part2)/100.0
     &              +real(part3)/1000.0+real(part4)/10000.0
     &              +real(part5)/100000.0
     &              +real(part6)/1000000.0
     &              +real(part7)/10000000.0
     &              +real(part8)/100000000.0
     &              +real(part9)/1000000000.0
               if((prob.lt.Prob_Dettach).AND.(prob.gt.0.000))then
                  C2_status(i)=0
               endif
            endif		
         enddo


         
ccccccccccccccccccccccccccccccccccccc



         do i=1,type_num
            do molecule=1,npart(i)
               do ii=1,intra_num(i)
                  status(ii,molecule,i)=status_new(ii,molecule,i)
               enddo
            enddo
         enddo
         complex_num=complex_num_new
         do i=1,complex_num
            complex_nb_num(i)=complex_nb_num_new(i)
            do j=1,complex_nb_ctg_num
               complex_nb_ctg(i,j)=complex_nb_ctg_new(i,j)
               complex_nb_ctg_idx(i,j)=complex_nb_ctg_idx_new(i,j)
               complex_nb_ctg_plm_idx(i,j)=
     &              complex_nb_ctg_plm_idx_new(i,j)
            enddo
         enddo



cccccccccccccccccccccccccccccccccccccccccccccc
c>>>   data outputx
cccccccccccccccccccccccccccccccccccccccccccccc


         if(MOD(itime,100).eq.0)then

            open (unit=10,file=
     &           'RTKLD_rec_Rcis_L5.dat',
     &           status='unknown',access='append')
            LR_num=0
            do i=1,complex_num
               if((complex_nb_ctg(i,1).eq.1).and.
     &              (complex_nb_ctg(i,2).eq.2))then
                  LR_num=LR_num+1
               endif
            enddo
            RS_num=0
            do i=1,complex_num
               if((complex_nb_ctg(i,1).eq.2).and.
     &              (complex_nb_ctg(i,2).eq.3))then
                  RS_num=RS_num+1
               endif
            enddo
            SE_num=0
            do i=1,complex_num
               if((complex_nb_ctg(i,1).eq.3).and.
     &              (complex_nb_ctg(i,2).eq.4))then
                  SE_num=SE_num+1
               endif
            enddo
            RR_num=0
            do i=1,complex_num
               if((complex_nb_ctg(i,1).eq.2).and.
     &              (complex_nb_ctg(i,2).eq.2))then
                  RR_num=RR_num+1
               endif
            enddo
            write(10,2200) 'index',itime,complex_num,LR_num,
     &           RS_num,SE_num,RR_num
            close(10)

         endif


         if(MOD(itime,10000).eq.0)then

            open (unit=10,file=
     &           'RTKLD_trj_Rcis_L5.pdb',
     &           status='unknown',access='append')

            do j=1,npart(1)
               write(10,2100) 'ATOM  ',j,' CA ','ALA', 
     &              'A',j,x(1,j,1),y(1,j,1),z(1,j,1)
               write(10,2100) 'ATOM  ',j,' CA ','ILE', 
     &              'A',j,x(2,j,1),y(2,j,1),z(2,j,1)
            enddo
            write(10,2102) 'TER'
            do j=1,npart(2)
               write(10,2100) 'ATOM  ',j,' CA ','VAL', 
     &              'B',j,x(1,j,2),y(1,j,2),z(1,j,2)
               write(10,2100) 'ATOM  ',j,' CA ','LYS', 
     &              'B',j,x(2,j,2),y(2,j,2),z(2,j,2)
               write(10,2100) 'ATOM  ',j,' CA ','ARG', 
     &              'B',j,x(3,j,2),y(3,j,2),z(3,j,2)
               write(10,2100) 'ATOM  ',j,' CA ','GLU', 
     &              'B',j,x(4,j,2),y(4,j,2),z(4,j,2)	 
            enddo
            write(10,2102) 'TER'
            do j=1,npart(3)
               write(10,2100) 'ATOM  ',j,' CA ','LEU', 
     &              'C',j,x(1,j,3),y(1,j,3),z(1,j,3)
               write(10,2100) 'ATOM  ',j,' CA ','ASP', 
     &              'C',j,x(2,j,3),y(2,j,3),z(2,j,3)
               write(10,2100) 'ATOM  ',j,' CA ','ASN', 
     &              'C',j,x(3,j,3),y(3,j,3),z(3,j,3)
               write(10,2100) 'ATOM  ',j,' CA ','GLN', 
     &              'C',j,x(4,j,3),y(4,j,3),z(4,j,3)	 
            enddo
            write(10,2102) 'TER'
            do j=1,npart(4)
               write(10,2100) 'ATOM  ',j,' CA ','GLY', 
     &              'D',j,x(1,j,4),y(1,j,4),z(1,j,4)
               write(10,2100) 'ATOM  ',j,' CA ','PRO', 
     &              'D',j,x(2,j,4),y(2,j,4),z(2,j,4)
               write(10,2100) 'ATOM  ',j,' CA ','TYR', 
     &              'D',j,x(3,j,4),y(3,j,4),z(3,j,4)
               write(10,2100) 'ATOM  ',j,' CA ','HIS', 
     &              'D',j,x(4,j,4),y(4,j,4),z(4,j,4)	 
            enddo
            write(10,2102) 'TER'		
            write(10,2102) 'END'
            

            close(10)

         endif

cccccccccccccccccccccccccccccccccccccccccccccc

      enddo

ccccccccccccccccccccccc


 2102 format(A3)   
 2101 format(I10)               
 2100 format(A6,I5,1x,A4,1x,A3,1x,A1,I4,4x,3F8.3)
 2200 format(A5,I10,I5,I5,I5,I5,I5)
 2300 format(I5,6I3)
 2500 format(I3,1x,A15,I5)

cccccccccccccccccccccc

      stop
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand3(r3)
      double precision s,u,v,r3
      s=65536.0
      u=2053.0
      v=13849.0
      m=r3/s
      r3=r3-m*s
      r3=u*r3+v
      m=r3/s
      r3=r3-m*s
      rand3=r3/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand4(r4)
      double precision s,u,v,r4
      s=65536.0
      u=2053.0
      v=13849.0
      m=r4/s
      r4=r4-m*s
      r4=u*r4+v
      m=r4/s
      r4=r4-m*s
      rand4=r4/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc                      
ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand5(r5)
      double precision s,u,v,r5
      s=65536.0
      u=2053.0
      v=13849.0
      m=r5/s
      r5=r5-m*s
      r5=u*r5+v
      m=r5/s
      r5=r5-m*s
      rand5=r5/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc                      
ccccccccccccccccccccccccccccccccccccccccccccccccc            
            
      subroutine gettheta(point_x,point_y,point_z,test_theta)

      implicit none
      real*8 point_x(3),point_y(3),point_z(3),test_theta
ccc
      real*8 lx(2),ly(2),lz(2),lr(2)
      real*8 conv,doth1,doth2,conv2
      integer in

ccccccccccccccccccc
c>>  bond lengt
ccccccccccccccccccc

      do in=1,2
         lx(in)=0
         ly(in)=0
         lz(in)=0
         lr(in)=0
      enddo

      lx(1)=point_x(2)-point_x(1)
      ly(1)=point_y(2)-point_y(1)
      lz(1)=point_z(2)-point_z(1)
      lr(1)=sqrt(lx(1)**2+ly(1)**2+lz(1)**2)
      
      lx(2)=point_x(3)-point_x(2)
      ly(2)=point_y(3)-point_y(2)
      lz(2)=point_z(3)-point_z(2)
      lr(2)=sqrt(lx(2)**2+ly(2)**2+lz(2)**2)
 
cccccccccccccccccccc
c>>  theta value
cccccccccccccccccccc
      
      test_theta=0

      conv=180/3.14159 
      do in=1,1
         doth1=-1*(lx(in+1)*lx(in)+ly(in+1)*ly(in)+lz(in+1)*lz(in))
         doth2=doth1/(lr(in+1)*lr(in))
         if(doth2.gt.1)then
            doth2=1
         endif
         if(doth2.lt.-1)then
            doth2=-1
         endif
         test_theta=acos(doth2)*conv
      enddo

cccccccccccccccc

      return
      end

cccccccccccccccccccccccccccccccccccc
                        
