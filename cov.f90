program cf_data_analysis
implicit none

integer, parameter :: dl=kind(1.d0)
integer ibin, iarea, imock, il, iline, jline, mock_used, iw 
integer, parameter :: Nmock=600, Nw=1, Nline=20, Nl=2
integer ::  ierr

character*512 :: filename, cfroot
character(4) :: x, y
character*16 :: w
character*512 :: outroot
real(dl) :: s(Nline), cftmp(Nline,Nmock), scftmp(Nline)
real(dl) :: cf(Nline), scf(Nline)
real(dl) :: cfall(Nmock,Nline*Nw), cfall_mean(Nline*Nw), cf_mean(Nline)
real(dl) :: cov(Nline*Nw,Nline*Nw), corr(Nline*Nw,Nline*Nw)
real(dl) :: err_cf(Nline),dummy 


outroot = '/home/hyzhang/Documents/data/hod/mock_MD_run1/cov/'

 cftmp = 0.d0; scftmp = 0.d0
 s=0.d0; cf=0.d0; scf=0.d0
 ierr = 0
 cfall =0.d0; cfall_mean=0.d0
 cov=0.d0; corr=0.d0

!  do iw=1, Nw !!redshift loop!!


	 mock_used=0
      do imock=1, Nmock
   
       write (x,'(I4.4)') imock

	   cfroot = '/home/hyzhang/Documents/data/hod/mock_MD_run1/2pcf/'
	   filename=trim(cfroot)//'twopoint_'//trim(x)//'.dat'
	  
	  open(unit=50, file=filename)
	  	  
	   do iline=1, Nline
	  
	   read(50,*,iostat=ierr) s(iline), cftmp(iline,imock)
           scftmp(iline) = s(iline)*cftmp(iline,imock) 
           cf(iline)=cf(iline)+cftmp(iline,imock) 
           scf(iline)=scf(iline)+scftmp(iline) 
           


	   end do !iline
	  
	  close(50)
	  
	  if(ierr==0) then 
	    mock_used = mock_used+1
	  else
	    exit
	  end if 	
	 end do !imock
 	 
	 cf_mean(:)  = cf(:)/dble(mock_used)
	 scf(:) = scf(:)/dble(mock_used)
	 
	
	write(*,*) mock_used,'mocks are used.'
		
  
  cfall_mean(1: Nline)            =  cf_mean(1:Nline)
  
  
       do imock=1, mock_used

          cfall(imock,1:Nline)            =  cftmp(1:Nline,imock)

       end do
!   end do ! iw	   
  
  cov=0.d0; corr=0.d0
  
 
     do imock=1, mock_used
   	  	  	 
	   do iline=1, Nline
	    do jline = 1, Nline
		
	      cov(iline, jline)  =  cov(iline, jline)  + &
		                      (cfall(imock,iline)-cfall_mean(iline)) * &
		                      (cfall(imock,jline)-cfall_mean(jline))/dble(mock_used-1)
		end do !jline
	   end do  !iline
	    
     end do !imock
	 

	   do iline=1, Nline
	    do jline = 1, Nline
			  
	      corr(iline, jline)  =  cov(iline, jline) / &
		                         sqrt(cov(iline,iline))/sqrt(cov(jline,jline))

	    end do !jline
	   end do  !iline
	
	filename =   trim(outroot)//'2pcf_cov.dat' 
	
	open(unit=50,file=filename)
	 do iline=1, Nline
	   write(50,'(20e20.10)') cov(iline,:)
	 end do
	close(50)
	
  
	filename =   trim(outroot)//'2pcf_err.dat'	
	open(unit=50,file=filename)
	
	 do iline=1, Nline	
		  err_cf(iline) = sqrt(cov(iline, iline))
	 end do
	
	 do iline=1, Nline
	 
	    write(50,'(3e20.10)') s(iline), &
		                         s(iline)*cf_mean(iline), s(iline)*err_cf(iline)
	 
	 end do
	 
 	close(50)
	
	
   filename =   trim(outroot)//'2pcf_corr.dat' 
	
	open(unit=50,file=filename)
	 do iline=1, Nline
	   write(50,'(20e20.10)') corr(iline,:)
	 end do
	close(50)
	   


end program cf_data_analysis
