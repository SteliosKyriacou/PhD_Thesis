      implicit double precision (a-h,o-z)
      dimension x(40),A(40),B(40,40),y(40)
c
      pi=4.d0*datan(1.d0)
      open(1,file='task.dat')
      open(2,file='A.dat')
      open(3,file='B.dat')
      read(1,*) kfreeall
      do i=1,kfreeall
	    read(1,*)x(i)
            read(2,*)A(i)
            read(3,*) (B(i,j),j=1,kfreeall)
      enddo

	close(1)
	close(2)
	close(3)
c     Rotate
      do i=1,kfreeall
        y(i)=0
        do j=1,kfreeall
          y(i)=y(i)+B(i,j)*x(j)
        enddo
         !y(i)=x(i)
      enddo
c
	sum1=0
	do i=1,kfreeall
	   !sum1=sum1+A(i)*y(i)**2
	   !sum1=sum1+y(i)**2
	   !sum1=sum1+A(i)*x(i)**2
	   sum1=sum1+1000**((i-1)/(kfreeall-1))*x(i)**2!inria sep
	   !sum1=sum1+1000**((i-1)/(kfreeall-1))*y(i)**2!inria non-sep
	enddo
	f1=sum1
c
      open(2,file='task.res')
        write(2,*) f1
      close(2)
c
      end
