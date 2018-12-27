c     --------------------------------------------------------------
      program USE_Plain_MultiLayer_Perceptron       ! 11.03.2008
c     --------------------------------------------------------------
c     Only for responses ; Not for the gradient
c     --------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
      parameter(kparam=21,mobj=1,kparamfix=100)
      parameter(maxlay=17,maxunit=200,maxwei=10000)
      dimension  bound1(kparam), bound2(kparam)
      dimension xfstring(kparam+mobj)
      dimension layunits(maxlay),laypo(maxlay+1),listwei(maxunit+1)
      dimension unitin(maxunit),unitout(maxunit),unitd(maxunit)
      dimension wei(maxwei),biasout(maxlay)
      dimension martyp(kparamfix)  !<0>:fixed,<1>:variable
      dimension valfix(kparamfix)
      integer iarray(3)
      character filen1*120,aaa
c
c-Read Basic Integer Quantities from the Training Data
c-----------------------------------------------------
      open(3,file='samples_nd.dat')
        read(3,*)kpatterns    !   useless
        read(3,*)kfree
        read(3,*)nobj
        write(*,*)' Number of training patterns = ',kpatterns
        write(*,*)' Number of free variables    = ',kfree
        write(*,*)' Number of responses         = ',nobj
      close(3)
      kfree2=kfree
      
c
c-Reads dofS lower & upper bounds
c--------------------------------
      call getarg(1,aaa)
      if (aaa.eq.'') then
        write(*,*)' Enter <1> if the data are allredy De-scaled '
        read(*,*) iDescaledOnOff
      else
        read (aaa,*) iDescaledOnOff
      endif
      if(iDescaledOnOff.ne.iDescaledOnOff) then
            write(*,*)' Enter filename with dof bounds '
            read(*,*)filen1
            open(1,file=filen1)
            read(1,*) kfreeall   !  total number of parameters
            kfree = 0
            do i=1,kfreeall
                  read(1,*) kdigi,bound1i,bound2i
                  if ( kdigi.eq.0 ) then
                        martyp(i) = 0      !  fixed
                        valfix(i) = bound1i
                  else
                        martyp(i) = 1      !  variable
                        kfree = kfree+1
                        bound1(kfree)=bound1i
                        bound2(kfree)=bound2i
                  endif
            enddo
            write (*,*)' #DOFS = ',kfree,' from ',kfreeall
            close(1)
            if(kfree.ne.kfree2) stop ' check kfree !!! '
      else
          kfreeall=kfree
          kfree=0
            do i=1,kfreeall
                        martyp(i) = 1      !  variable
                        kfree = kfree+1
                        bound1(kfree)=0
                        bound2(kfree)=1
            enddo
            write (*,*)' #DOFS = ',kfree,' from ',kfreeall
      endif
c-Build Network
c--------------
      call build_perceptron(maxlay,layunits,laypo,kfree,nobj,
     1           nunit,maxunit,unitout,unitin,unitd,wei,listwei,
     2           maxwei,nwei,maxLiter,eta,nlay,kk,monoff,kpatterns)
c
c-Read Computed Weights
c----------------------
      call Read_ANNstate(nlay,nunit,nwei,wei,maxwei,kter)
   
c-Read Prediction point (dimensional)
c------------------------------------
      open(1,file='task.dat')
        read(1,*)nfreeall
        if(nfreeall.ne.kfreeall) stop ' check nfreeall !!! '
        k2=0
        do i=1,nfreeall
         if(martyp(i).eq.0)then
           read(1,*)
         else
           k2=k2+1
           read(1,*)xfst
ccc        xfstring(k2)=xfst            !!!!!  kirk
           xfstring(k2)=(xfst-bound1(k2))/(bound2(k2)-bound1(k2))
         endif
        enddo
      close(1)
c
c-Guess Responses
c----------------
       call use_perceptron(maxlay,layunits,laypo,kfree,nobj,
     1           nunit,maxunit,unitout,unitin,unitd,wei,listwei,
     2           maxwei,nwei,nlay,xfstring,kparam,mobj,kk)
c
      open(1,file='task.res')
        do i=1,nobj
          write(1,*) xfstring(kfree+i)
        enddo
      close(1)

      end
c
c
c
c **********************************************************************
      subroutine build_perceptron(maxlay,layunits,laypo,kfree,nobj,
     1           nunit,maxunit,unitout,unitin,unitd,wei,listwei,
     2           maxwei,nwei,maxLiter,eta,nlay,kk,monoff,kpatterns)
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension layunits(maxlay),laypo(maxlay+1),biasout(maxlay)
      dimension unitin(maxunit),unitout(maxunit),unitd(maxunit)
      dimension wei(maxwei)
      dimension listwei(maxunit+1)
c
      open(1,file='ann.ini')
      read(1,*)  ! skip 
      read(1,*)  ! skip 
c-maxLiter: Max. Learning iterations
      read(1,*) maxLiter
c-eta: Learning rate
      read(1,*) eta
c-monoff: Sequential (on-line) or batch (off-line) training
      read(1,*) monoff ! on=1, off=2
      if(monoff.eq.1) then
        write(*,*)' #####  On-line or SEQUENTIAL Training  '
      elseif(monoff.eq.2) then
        write(*,*)' #####  Off-line or BATCH Training  '
      else
        stop ' Correct MONOFF (1 or 2) in ANN.INI'
      endif
c-nlay: Number of Network Levels (ilay=1:input, ilay=nlay: output)
      read(1,*) nlay
      if(nlay.gt.maxlay) stop 'InitANN: Increase maxlay'
c-kk: <0>=no bias units, <1>=with bias units
      read(1,*) kk
c-biasout: Output value of each bias unit (if zero, ignored)
      do ilay =1,nlay
        read(1,*) layunits(ilay),biasout(ilay)
        if(ilay.eq.1   .and.kfree.gt.0)  layunits(1)   =kfree
        if(ilay.eq.nlay.and.nobj .gt.0)  layunits(nlay)=nobj
        if(ilay.ne.nlay) layunits(ilay) = layunits(ilay)+kk ! plus bias
      enddo
      close(1)
c
c     Construct data-structures, lists,...
c     ------------------------------------
c-laypo: Pointer to the first unit of each layer - including bias
c-unitout: Output from each unit (after activation is applied)
c-nunit: Total number of units - including bias
c-unitin: Input per unit (sum)
c-unitd: Deltas per unit
      laypo(1)=1
      do ilay=1,nlay
        laypo(ilay+1)=laypo(ilay)+layunits(ilay)
        if(kk.eq.1.and.ilay.ne.nlay)
     1         unitout(laypo(ilay+1)-1)=biasout(ilay) ! last=bias
      enddo
      nunit = laypo(nlay+1)-1
      if(nunit.gt.maxunit) stop 'InitANN: Increase maxunit'
c
c-listwei: Pointers (list) to the wei array (for each unit points
c               to the receiving weights)
      do i=laypo(1),laypo(2)-1 ! no incoming synapses for layer 1
        listwei(i)=0
      enddo
      
      listwei(laypo(2))=1
      do ilay=2,nlay   !   for all but the input layer
        kunitlay=layunits(ilay-1) ! #links to previous layer (incl_bias)
        if(kk.eq.0.or.ilay.eq.nlay)then
          lastun=laypo(ilay+1)-1  !  last of current layer
        else
          lastun=laypo(ilay+1)-2  !  last of current layer
        endif
        do iunit=laypo(ilay),lastun  ! all units of this layer
          listwei(iunit+1)=listwei(iunit)+kunitlay
        enddo
        if(kk.eq.1.and.ilay.lt.nlay) listwei(lastun+2)=listwei(lastun+1)
      enddo
c-nwei: total number of weights
      nwei = listwei(nunit+1)-1
      if(nwei.gt.maxwei) stop 'InitANN: Increase maxwei'
c
c     Initialize Weights: small number
c     --------------------------------
      do i=1,nwei
        wei(i)= (-1)**i*(1+dfloat(i)/dfloat(nwei))*1.d-02    !  kirk
      enddo
c
      return
      end
c
c **********************************************************************
      subroutine Read_ANNstate(nlay,nunit,nwei,wei,maxwei,kter)
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension wei(maxwei)
      write(*,*)' ##########   Reads file WEIGHTS.DAT'
      open(1,file='weights.dat')
      read(1,'(10i10)') nlayb,nunitb,nweib,kter
      write(*,*)' of the ',kter,' epoch !!! '
      if(nlayb.ne.nlay.or.nunitb.ne.nunit.or.nweib.ne.nwei)
     1    stop    ' STOP:  Incompatible Continuation File'
      read(1,*) (wei(i),i=1,nwei)
      close(1)
      return
      end
c
c **********************************************************************
      subroutine use_perceptron(maxlay,layunits,laypo,kfree,nobj,
     1           nunit,maxunit,unitout,unitin,unitd,wei,listwei,
     2           maxwei,nwei,nlay,xfstring,kparam,mobj,kk)
c **********************************************************************
      implicit double precision (a-h,o-z)
c-Logistic
      act(hv)=1.d0/(1.d0+1.d0*dexp(-hv))
c-Hyperbolic Tangent
c     parameter(act_a=1.7159,act_b=0.6666d0)
c     daxp(hv)=dexp(dmax1(dmin1(1.d+2,hv),-1.d+2))
c     daxp(hv)=dexp(hv)
c     dhtan(hv)=(daxp(2.d0*hv)-1.d0)/(daxp(2.d0*hv)+1.d0)
c     act(hv)=act_a*dhtan(act_b*hv)
c-Dummy (for output unit)
      actL(hv)=hv
c-Dummy (for output unit)
c     actL(hv)=act(hv)
c
      dimension layunits(maxlay),laypo(maxlay+1),listwei(maxunit+1)
      dimension unitin(maxunit),unitout(maxunit),unitd(maxunit)
      dimension xfstring(kparam+mobj),wei(maxwei)
c
c-Presents the new input pattern to input units
      do j=1,kfree
        unitout(j) = xfstring(j)
      enddo
c-Propagates the signal forward, up to the exit
      do 15 ilay=2,nlay ! for each layer (except input)
          lastun = laypo(ilay+1)-1-kk  ! its last unit (excl. bias)
          if(ilay.eq.nlay)lastun = laypo(ilay+1)-1
          do irunit=laypo(ilay),lastun ! layer's receiving units
            sum=0.d0
            do idum=1,layunits(ilay-1) ! signals from prev. layer
              iunit = laypo(ilay-1) + (idum-1)
              iwei  = listwei(irunit) + (idum-1)
              sum   = sum + wei(iwei)*unitout(iunit)
            enddo    !  previous layer
            unitin(irunit)  = sum
            if(ilay.eq.nlay)then
              unitout(irunit) = actL(sum)
            else
              unitout(irunit) = act(sum)
            endif
          enddo    !  current layer
 15   continue
c-Restores responses
      do i=1,nobj
        xfstring(kfree+i)=unitout(laypo(nlay)-1+i)
      enddo
c
      return
      end
