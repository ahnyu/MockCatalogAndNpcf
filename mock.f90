module Precision
       implicit none
       integer, parameter :: dl = KIND(1.d0)
       integer, parameter :: sp = KIND(1.0)
end module Precision
module Pub
    use Precision
    implicit none
    real(dl), parameter :: tol=1.d-6
    real(dl), parameter :: PI=3.1415926d0
    real(dl) :: rs_sat,c_sat,rand2
!    integer, parameter :: nhalo=1000000
    integer, parameter :: nhalo=176330227
    
    contains
       subroutine init_random_seed()

       INTEGER :: i, n, clock
       INTEGER, DIMENSION(:), ALLOCATABLE :: seed

       CALL RANDOM_SEED(size = n)
       ALLOCATE(seed(n))

       CALL SYSTEM_CLOCK(COUNT=clock)

       seed = clock + 37 * (/ (i - 1, i = 1, n) /)
       CALL RANDOM_SEED(PUT = seed)

       DEALLOCATE(seed)
       end


end module Pub
Subroutine bisection(f,x1,x2,eps,Root,flag)
!============================================================
! Solutions of equation f(x)=0 on [x1,x2] interval
! Method: Bisectional (closed domain) (a single root)
! Alex G. January 2010
!------------------------------------------------------------
! input ...
! f   - function - evaluates f(x) for any x in [x1,x2]
! x1  - left endpoint of initial interval
! x2  - right endpoint of initial interval
! eps - desired uncertainity of the root as |b-a|<eps
! output ...
! Root  - root of the equation f(x)=0
! flag  - indicator of success
!         >0 - a single root found, flag=number of iterations
!          0 - no solutions for the bisectional method
!         <0 - not a root but singularity, flag=number of iterations
!
! Comments: Function f(x) has to change sign between x1 and x2
!           Max number of iterations - 200 (accuracy (b-a)/2**200)
!====================================================================
use Precision
implicit none
real(dl):: f, x1, x2, eps, Root
real(dl):: a, b, c
integer i, flag
integer, parameter:: iter=200

!* check the bisection condition

if(f(x1)*f(x2)>0.0) then
  flag = 0
  return
end if

!* initialize calculations
a=x1
b=x2

!* Iterative refining the solution 
do i=1,iter
  c=(b+a)/2.0
  if(f(a)*f(c).le.0.0) then
      b = c
    else
      a=c
  end if
! condition(s) to stop iterations)
  if(abs(b-a)<= eps) exit  
end do
Root=(b+a)/2.0

!* check if it is a root or singularity
if (abs(f(Root)) < 1.0) then
  flag=i
  else
  flag = -i
end if
end subroutine bisection


function nfw(r)
use Precision
use Pub
implicit none
real(dl) nfw,r
nfw=(dlog((rs_sat+r)/rs_sat)-r/(r+rs_sat))/(dlog(1.d0+c_sat)-c_sat/(1+c_sat))-rand2
end function nfw

      
    
program hod_selection
      use Precision
      use random
      use Pub
      use libnpcf
      implicit none
      type,bind(C) :: float3
                real(sp) :: x,y,z
      end type

      type(npcf) :: corr
      type(float3), dimension(:),allocatable :: pos, triangles, temp
      integer :: s1, timesRans, numShells, numTriangles, dummy
      real(dl) :: r_max, r_min, V_box
      real(dl), dimension(:), allocatable :: threePoint, twoPoint, shells 
      real(dl), dimension(:), allocatable :: mass, x, y, z, rvir, rs 
      real(dl) :: m_cut, sigma, kappa, m1, alpha !baseline hod parameters
      real(sp) :: N_cent, N_sat
      real(dl) :: rand, rand_mu, rand_phi
      real(dl) :: r,dx,dy,dz
      real(dl) :: costheta,sintheta,cosphi,sinphi
      real(dl) :: g_x, g_y, g_z
      real(dl) :: starttime1, starttime2,stoptime1,stoptime2,stoptime3,stoptime4
      real(dl), external :: nfw 
      integer, external :: ZBQLPOI
      character(4) :: mock_index
      
      integer :: nmock,i,j,count_all,flag,count_sat,imock
      
      m_cut=10.d0**13.08d0
      m1=10.d0**14.06d0
      sigma=0.98d0
      kappa=1.13d0
      alpha=0.9d0
      timesRans = 10
      numShells = 20
      V_box = 2500.d0**3
      r_max=20.d0
      r_min=0.d0
      nmock=600
      
      corr = npcf(timesRans, numShells, V_box, r_min, r_max)
      numTriangles = corr%getNumTriangles()

      allocate(shells(numShells))
      allocate(twoPoint(numShells))
      allocate(triangles(numTriangles))
      allocate(threePoint(numTriangles))
      dummy = corr%getShells(shells)
      dummy = corr%getTriangles(triangles)

      call cpu_time(starttime1) 
      allocate(mass(nhalo))
      allocate(x(nhalo))
      allocate(y(nhalo))
      allocate(z(nhalo))
      allocate(rvir(nhalo))
      allocate(rs(nhalo))
      open(unit=100, file='/home/hyzhang/Documents/data/hod/MDhalos.txt')
     
      do i=1,nhalo

          read(100,*) mass(i),rvir(i),rs(i),x(i),y(i),z(i)
          rvir(i)=rvir(i)/1000.d0
          rs(i)=rs(i)/1000.d0
      end do
      close(100)
      call cpu_time(stoptime1)
      write(*,*) 'Reading Halo File Complete',' reading time =',stoptime1-starttime1
      write(*,*) '5-parameter HoD model galaxy selection start'
      do imock=1,nmock
      call cpu_time(starttime2)
      count_all=0
      count_sat=0
      write(mock_index,'(I4.4)') imock
!      open(unit=200, file='/home/hyzhang/Documents/data/hod/mock_MD/gal_catalog/central_output_'//trim(mock_index)//'.dat')
!      open(unit=300, file='/home/hyzhang/Documents/data/hod/mock_MD/gal_catalog/satellite_output_'//trim(mock_index)//'.dat')
      open(unit=400, file='/home/hyzhang/Documents/data/hod/mock_MD/gal_catalog/all_output_'//trim(mock_index)//'.dat')
      call init_random_seed() 
      allocate(temp(nhalo))
          
      do i=1,nhalo
          N_cent = 1.d0/2.d0*erfc(dlog10(m_cut/mass(i))/(dsqrt(2.d0)*sigma))
          N_sat = ((mass(i)-kappa*m_cut)/m1)**alpha
          call RANDOM_NUMBER(rand)
!          if(random_Poisson(N_cent,.true.).gt.0) then
          if(rand.lt.N_cent) then
              if(x(i).le.2500.d0.and.y(i).le.2500.d0.and.z(i).le.2500.d0) then
              count_all=count_all+1
              temp(count_all)%x=x(i)
              temp(count_all)%y=y(i)
              temp(count_all)%z=z(i)
!              write(200,"(3F9.3)") temp(count_all)
              write(400,"(3F9.3)") temp(count_all)
              end if

          end if
          if(N_sat.ge.1.d0) then
               !write(*,*) N_sat
               do j=1,random_Poisson(N_sat,.true.)
               count_all=count_all+1
               count_sat=count_sat+1
               call RANDOM_NUMBER(rand2)
               rs_sat=rs(i)
               c_sat=rvir(i)/rs(i)
               call bisection(nfw,0.d0,rvir(i),1.d-4,r,flag)
               call RANDOM_NUMBER(rand_mu)
               costheta=rand_mu*2.d0-1.d0
               sintheta=dsqrt(1-rand_mu**2)
               call RANDOM_NUMBER(rand_phi)
               cosphi=dcos(rand_phi*2*PI)
               sinphi=dsin(rand_phi*2*PI)
               dx=r*sintheta*cosphi
               dy=r*sintheta*sinphi
               dz=r*costheta
               temp(count_all)%x=x(i)+dx               
               temp(count_all)%y=y(i)+dy               
               temp(count_all)%z=z(i)+dz               
               if(temp(count_all)%x.gt.2500.d0) temp(count_all)%x=temp(count_all)%x-2500.d0
               if(temp(count_all)%y.gt.2500.d0) temp(count_all)%y=temp(count_all)%y-2500.d0
               if(temp(count_all)%z.gt.2500.d0) temp(count_all)%z=temp(count_all)%z-2500.d0
               if(temp(count_all)%x.lt.0.d0) temp(count_all)%x=temp(count_all)%x+2500.d0
               if(temp(count_all)%y.lt.0.d0) temp(count_all)%y=temp(count_all)%y+2500.d0
               if(temp(count_all)%z.lt.0.d0) temp(count_all)%z=temp(count_all)%z+2500.d0
!               write(300,"(3F9.3)") temp(count_all)
               write(400,"(3F9.3)") temp(count_all)
               end do               
           end if
              
      end do
      call cpu_time(stoptime2)
      write(*,*) 'galaxy selection and galaxy catalog file writing complete',' time = ',stoptime2-starttime2
      write(*,*) 'mock index = ',imock,' satellite number = ',count_sat
      write(*,*) 'mock index = ',imock,' total number = ',count_all
      write(*,*) 'D.P. npcf library based 2pcf and 3pcf calculation start'
      allocate(pos(count_all))
      pos(1:count_all)=temp(1:count_all)

      dummy=corr%setNumParticles(count_all)
      dummy=corr%calculateCorrelations(pos)
      dummy=corr%get2pt(twoPoint)
      dummy=corr%get3pt(threePoint)
      deallocate(pos)
      deallocate(temp)
      
      open(unit=500, file='/home/hyzhang/Documents/data/hod/mock_MD/2pcf/twopoint_'//trim(mock_index)//'.dat')
      open(unit=600, file='/home/hyzhang/Documents/data/hod/mock_MD/3pcf/threepoint_'//trim(mock_index)//'.dat')
      do i=1,numshells
                write(500,*) shells(i),twoPoint(i)
      end do
      do i=1,numTriangles
                write(600,*) triangles(i), threePoint(i)
      end do 
      call cpu_time(stoptime3)
      write(*,*) '2pcf and 3pcf calculation and file writing complete', 'time = ', stoptime3-stoptime2
!      close(200)
!      close(300)
      close(400)
      close(500)
      close(600)
      write(*,*) 'time for mock index ',imock,' = ',stoptime3-starttime2 
      end do!nmock
      call cpu_time(stoptime4)
      write(*,*) nmock,' mocks calculation complete, total time = ',stoptime4-starttime1 
      deallocate(mass)
      deallocate(x)
      deallocate(y)
      deallocate(z)
      deallocate(rs)
      deallocate(rvir)
      deallocate(shells)
      deallocate(triangles)
      deallocate(twoPoint)
      deallocate(threePoint)

    end program hod_selection
    
              

