      program welded_beam
      implicit double precision (a-h, o-z)
!
!     Read task.dat
!     -------------      
      open(1,file='task.dat')
         read(1,*) ibid
         read(1,*) x1
         read(1,*) x2
         read(1,*) x3  
         x4=x3  
      close(1)
!
!     Hardcoded data
!     --------------
      pp   = 6000.d0    ! [lbs]
      el   = 14.d0
      dmax = 0.25d0     ! [in]
      ee   = 30.d6      ! [psi]
      gg   = 12.d6      ! [psi]
      tmax = 13600.d0   ! [psi]
      smax = 30000.d0   ! [psi]
!
!     Cost function
!     -------------      
      fx = 1.10471d0*x1*x1*x2 + 0.04811d0*x3*x4*(14.d0 + x2)
!
!     Constraints
!     -----------
!     preparation
      t1 = pp/(dsqrt(2.d0)*x1*x2)     
      am = pp*(el + 0.5d0*x2)
      term0 = 0.25d0*(x1+x3)*(x1+x3)
      ar = dsqrt( 0.25d0*x2*x2 + term0 )
      aj = 2.d0*dsqrt(2.d0)* x1*x2* (x2*x2/12.d0 + term0) 
      t2 = am*ar/aj
      tx = dsqrt( t1*t1 + t1*t2*x2/ar + t2*t2 )
!
      sx = 6.d0*pp*el/(x4*x3*x3)
      dx = 4.d0*pp* el**3.d0/(ee*x4*x3**3.d0)
      term1 = 4.013d0*ee/el/el*dsqrt(x3*x3*x4**6.d0/36.d0)
      term2 = 1.d0 - 0.5d0*x3/el*dsqrt(0.25d0*el/gg)
      pcx = term1*term2
!
      g1 = tx - tmax
      g2 = sx - smax
      g3 = x1 - x4
      g4 = 1.10471d0*x1*x1 + 0.04811d0*x3*x4*(14.d0 + x2) - 50.d0
      g5 = 0.125d0 - x1
      g6 = dx - dmax
      g7 = pp - pcx
!
!     Erase previous files (task.res, task.cns)
!     -----------------------------------------
      open (1,file='task.res')
      close(1,status='delete')
      open (1,file='task.cns')
      close(1,status='delete')
!      
!     Write task.res
!     --------------     
      open (1,file='task.res')!minimization of cost and deflection
         write(1,*) fx
!	 write(1,*) sx 
	 write(1,*) dx 
      close(1)
!
!     Write task.cns
!     --------------      
      open(1,file='task.cns')
         write(1,*) g1
         write(1,*) g2
         write(1,*) g3
         write(1,*) g4
         write(1,*) g5
         write(1,*) g6
         write(1,*) g7
      close(1)
!
      stop
      end program welded_beam
