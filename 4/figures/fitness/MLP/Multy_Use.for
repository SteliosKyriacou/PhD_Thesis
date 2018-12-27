c     --------------------------------------------------------------
      program Multyuse       ! 08.05.2008
c     --------------------------------------------------------------
c     Creates the fitted curve 
c	Only for kfree=1 (for now) 
c     --------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
      
	inum=20
      open(3,file='curve.dat')
	do i=1,inum
	do jj=1,inum
		x=0.0 +(i-1)*(1.d0/inum)*(0.45- 0.0)
		y=0.25 +(jj-1)*(1.d0/inum)*(0.75-0.25)
		j=system('del task.dat')
		j=system('del task.res')
		open(1,file='task.dat')
		write(1,*) 2
		write(1,*) x
		write(1,*) y
		close(1)
		j=system('./use_plainMLP.exe 1')
		open(2,file='task.res')
		read(2,*) res
		close(2)
		write(3,*) x,y,res 
	enddo
	write(3,*) 
       enddo
      close(3)
      end
