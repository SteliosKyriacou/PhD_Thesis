c     --------------------------------------------------------------
      program Multyuse       ! 08.05.2008
c     --------------------------------------------------------------
c     Creates the fitted curve 
c	Only for kfree=1 (for now) 
c     --------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
      
	inum=100
      open(3,file='curve.dat')
	do i=1,inum
		x=-0.25 +(i-1)*(1.d0/inum)*(0.25+0.25)
		j=system('del task.dat')
		j=system('del task.res')
		open(1,file='task.dat')
		write(1,*) 1
		write(1,*) x
		close(1)
		j=system('./use_plainMLP.exe 1')
		open(2,file='task.res')
		read(2,*) res
		close(2)
		write(3,*) x,y,res 
	enddo
      close(3)
      end
