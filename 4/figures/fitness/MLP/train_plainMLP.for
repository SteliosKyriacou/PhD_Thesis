c     --------------------------------------------------------------
      program Plain_MultiLayer_Perceptron_for_testing    ! 15.03.2008
c     --------------------------------------------------------------
c     Only for responses ; Not for the gradient
c     --------------------------------------------------------------
c     Check the activation, there is likely a multiplier (2*act, 1*dact)
c     Sequential training: Not ready yet. Needs 2D matrices for storage
c     Trained on known GRAD, too!
c     Works only with one OUTPUT
c     Attention:  bounds should be [0-1], see DENOM (possible mistake)
c     Logistic function of Hyperbolic tangent (comment in & out)
c     Last layer without activation !!!!!  (actL)
c     --------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
      parameter(kparam=21,mobj=1,kpoplme=600)
      parameter(maxdbindiv=1000,maxlay=7,maxunit=200,maxwei=10000)
      parameter(maxdb=maxdbindiv*(kparam+mobj))
      parameter(maxdbder=maxdbindiv*kparam*mobj)
      dimension  bound1(kparam), bound2(kparam)
      dimension iwksp(kpoplme)
      dimension xfstring(kparam+mobj),dbder(0:maxdbder)
      dimension db(0:maxdb),dweipat(maxwei),apold(maxwei)
      dimension layunits(maxlay),laypo(maxlay+1),listwei(maxunit+1)
      dimension unitin(maxunit),unitout(maxunit),unitd(maxunit)
      dimension wei(maxwei),biasout(maxlay),dwei(maxwei),dweiold(maxwei)
      dimension weiold(maxwei)
      integer iarray(3)
c
c-Read Data
c----------
      open(3,file='samples_nd.dat')
      read(3,*)kpatterns
      read(3,*)kfree
      read(3,*)nobj
      write(*,*)' Number of training patterns = ',kpatterns
      write(*,*)' Number of free variables    = ',kfree
      write(*,*)' Number of responses         = ',nobj
      if(nobj.gt.mobj) stop 'Increase mobj'
      if(kpatterns.gt.maxdbindiv) stop 'Increase maxdbindiv!'
      do i=1,kpatterns
        read(3,*) (xfstring(ip),ip=1,kfree+nobj)
        call fill_in_DB(maxdb,db,mdb,xfstring,kparam,mobj,kfree,nobj)
      enddo
      close(3)

c-Network Training
c-----------------
      do i=1,kpatterns
        iwksp(i)=i  !   take all of them - skip distance criterion
      enddo

      call build_perceptron(maxlay,layunits,laypo,kfree,nobj,
     1           nunit,maxunit,unitout,unitin,unitd,wei,listwei,
     2           maxwei,nwei,maxLiter,eta,nlay,kk,monoff,kpatterns)
      interactive=1
      do iphase=0,0    !  kirk08 ; Only for responses; not derivatives
      call train_perceptron(maxlay,layunits,laypo,kfree,nobj,
     1           nunit,maxunit,unitout,unitin,unitd,wei,listwei,
     2           maxwei,nwei,maxLiter,eta,nlay,interactive,monoff,
     3           kpatterns,iwksp,kpoplme,db,maxdb,kk,dwei,dweiold,
     4           dweipat,apold,weiold,xfstring,kparam,mobj,
     5           dbder,maxdbder,bound1,bound2,iphase)
      enddo
c
      end
c
c
c
c **********************************************************************
      subroutine fill_in_DB
     1    (maxdb,db,mdb,xfstring,kparam,mobj,kfree,nobj)
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension xfstring(kparam+mobj),db(0:maxdb)
c
      kur_end=mdb*(kfree+nobj)
      new_end=kur_end+kfree+nobj
      mdb=mdb+1
      if (new_end.gt.maxdb) stop  ' Increase MAXDB '
      do i=1,kfree+nobj
         db(kur_end+i) = xfstring(i)
      enddo
c
      return
      end
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
        wei(i)= (-1)**(i+1)*(1+dfloat(i+5)/dfloat(nwei))*1.d-02 !  kirk
      enddo
c
c     Printouts about the net structure
c     ---------------------------------
      write(*,*) ' #####  Number of units - including bias = ',nunit
      write(*,*) ' #####  Number of unknown weights        = ',nwei
      write(*,*) ' #####  Number of Training Patterns      = ',kpatterns
      kpa=kpatterns*(1+kfree)*nobj
      write(*,*) ' #####  Number of Conditions             = ',kpa
c     if (kpa.gt.nwei) then
c       write(*,*)' >> '
c       write(*,*)' >>  More conditions than unknowns!!!!'
c       write(*,*)' >>  Shall I continue (press ENTER)??? ',char(72)
c       write(*,*)' >> --->>>  '
c       pause
c     endif
      goto 50   !     return      !   kirk attention
      write(*,*) ' *** LAYUNITS( ) : '
      do ilay =1,nlay
        write(*,*) ilay,layunits(ilay)
      enddo
      write(*,*) ' *** LAYPO( ) : '
      do ilay=1,nlay+1
        write(*,*) ilay,laypo(ilay)
      enddo
      write(*,*) ' *** LISTWEI( ) : '
      do i=1,nunit+1
        write(*,*) i,listwei(i)
      enddo
c
c     Possible Continuation
c     ---------------------
 50   write(*,*)' Enter <1> for continuation of weights '
      read(*,*)konti
      if(konti.eq.1)then
         call Read_ANNstate(nlay,nunit,nwei,wei,maxwei,kter)
      endif
c
      return
      end
c
c **********************************************************************
      subroutine train_perceptron(maxlay,layunits,laypo,kfree,nobj,
     1           nunit,maxunit,unitout,unitin,unitd,wei,listwei,
     2           maxwei,nwei,maxLiter,eta,nlay,interactive,monoff,
     3           kpatterns,iwksp,kpoplme,db,maxdb,kk,dwei,dweiold,
     4           dweipat,apold,weiold,xfstring,kparam,mobj,
     5           dbder,maxdbder,bound1,bound2,iphase)
c **********************************************************************
c-IPHASE=0 (training on responses only)
c-IPHASE=1 (training on responses and derivatives)
      implicit double precision (a-h,o-z)
c-Logistic Activation Function
      act(hv)=1.d0/(1.d0+1.d0*dexp(-hv))
      dact(hv)=1.d0*(1.d0-act(hv))*act(hv)
c-Hyperbolic Tangent
      parameter(act_a=1.7159,act_b=0.6666d0)
c     daxp(hv)=dexp(dmax1(dmin1(1.d+2,hv),-1.d+2))
c     daxp(hv)=dexp(hv)
c     dhtan(hv)=(daxp(2.d0*hv)-1.d0)/(daxp(2.d0*hv)+1.d0)
c     act(hv)=act_a*dhtan(act_b*hv)
c     dact(hv)=act_a*4.d0*act_b*daxp(2.d0*act_b*hv)/
c    1         (daxp(2.d0*act_b*hv)+1.d0)**2
c-Dummy (for output unit)
      actL(hv)=hv
      dactL(hv)=1.d0
c-Dummy (for output unit)
c     actL(hv)=act(hv)
c     dactL(hv)=dact(hv)
c
      dimension layunits(maxlay),laypo(maxlay+1)
      dimension unitin(maxunit),unitout(maxunit),unitd(maxunit)
      dimension wei(maxwei),db(0:maxdb),dwei(maxwei),dweiold(maxwei)
      dimension dweipat(maxwei),apold(maxwei),weiold(maxwei)
      dimension listwei(maxunit+1),iwksp(kpoplme)
      dimension xfstring(kparam+mobj)
      dimension dbder(0:maxdbder)
      dimension unitd1(maxunit,kparam),unitd2(maxunit,kparam)
      dimension  bound1(kparam), bound2(kparam)
      dimension xsample(100),ysample(100)
      dimension asample(3,3),rsample(3),ssample(3)
c
  56  write(*,*)' >> '
      write(*,*)' >>  Choose Solver, ENTER : '
      write(*,*)' >>  <1> for Steepest-Descent              '
      write(*,*)' >>  <2> for Polak-Ribiere                 '
      write(*,*)' >>  <3> for Fletcher-Reeves               '
      write(*,*)' >> --->>>  '
      read(*,*) ksolve
      if(ksolve.lt.0.or.ksolve.gt.3) goto 56
c
      open(12,file='conv')
      maxiter = maxLiter
      if(interactive.ne.0) then
        write(*,*)' >> '
        write(*,*)' >>  ENTER training iterations and eta ' 
        write(*,*)' >>  (-ve Iterations means Continuation) : '
        write(*,*)' >> --->>>  '
        read(*,*) maxiter, eta
        if(maxiter.lt.0) then
          call Read_ANNstate(nlay,nunit,nwei,wei,maxwei,kter)
          maxiter=-maxiter
          do iwei=1,nwei
            weiold(iwei)=wei(iwei)
          enddo
        endif
      endif
      kter = 0        !   always
c
c-Initializes-Zeroes old weight corrections
      do iwei=1,nwei
        dweiold(iwei)=0.d0
      enddo
c
c-Loop on Epochs starts here
 100  do 101 iter=1,maxiter
      kter = kter+1           ! Epoch counter
      cost = 0.d0             ! Initial epoch error
c-Re-orders training patterns (shifts current iwksp)
      do i=1,kpatterns
        iwksp(i)=iwksp(i)+1   !  sort of rotation among patterns
        if(iwksp(i).gt.kpatterns)iwksp(i)=iwksp(i)-kpatterns !close loop
      enddo
c-Zeroes DWEI for the new epoch, needed in batch training
      do iwei=1,nwei
        dwei(iwei)=0.d0
      enddo
c-Presents one pattern to the input units, at a time
ccc   randompat=ran3(idum) !  kirk, will multiply dwei upon addition ...
      do 66 ikp=1,kpatterns
        npat=iwksp(ikp)    !  the presented pattern
        idpat= (npat-1)*(kfree+nobj)
c-Presentation of a pattern at the first layer units
        do j=1,kfree
          unitout(j) = db(idpat+j) ! input=output at the first layer
        enddo
c-Propagates the signal forward, up to the exit unit(s)
        do 15 ilay=2,nlay ! for each layer (except the input)
          lastun = laypo(ilay+1)-1-kk  ! its last unit (excl. bias)
          if(ilay.eq.nlay)lastun = laypo(ilay+1)-1
          do irunit=laypo(ilay),lastun ! layer's receiving units
            sum=0.d0
            do idum=1,layunits(ilay-1) ! signals from prev. layer
              iunit = laypo(ilay-1) + idum-1  ! prev. layer's units
              iwei  = listwei(irunit) + idum-1 ! corresponding weight
              sum   = sum + wei(iwei)*unitout(iunit)
            enddo
            unitin(irunit)  = sum  ! input to activation f.
            if(ilay.eq.nlay)then  ! this is the guessed response
               unitout(irunit) = actL(sum)  !output from activation f.
            else
               unitout(irunit) = act(sum)  !output from activation f.
            endif
        enddo
 15   continue


c-SFHNA-PRINTOUT   ------
      if(iter.eq.maxiter)then
         open(22,file='unitout.dat')
           do ilay=2,nlay ! for each layer (except the input)
             lastun = laypo(ilay+1)-1-kk  ! its last unit (excl. bias)
             if(ilay.eq.nlay)lastun = laypo(ilay+1)-1
             do irunit=laypo(ilay),lastun ! layer's receiving units
               write(22,*)ilay,irunit, unitout(irunit)
             enddo
           enddo
         close(22)
      endif
c-SFHNA   ------
      
 
c-Computes deltas for the output layer
      i=0
      do iunit= laypo(nlay),laypo(nlay+1)-1 ! output layer's units
        i=i+1
        deviation=db(idpat+kfree+i)-unitout(iunit)
        unitd(iunit) = -dactL(unitin(iunit))*deviation
      enddo    !  iunit
c-Sweeps backwards and computes deltas for the hidden layers
      do ilay=nlay-1,2,-1   !  current layer = ILAY
        do icunit=laypo(ilay),laypo(ilay+1)-1 ! units of layer ILAY
          sum=0.d0
          indx=icunit-laypo(ilay)   !  first value=0
          nnext = layunits(ilay+1)  !  number of next layer units
          if(ilay+1.ne.nlay) nnext=nnext-kk  ! except bias, if any
          sump=dact(unitin(icunit))
          do idum=1,nnext ! next layer units, except bias
            iunit= laypo(ilay+1) + idum-1  ! unit on upper level
            iwei = listwei(iunit)+ indx
            sum  = sum+wei(iwei)*unitd(iunit)
          enddo     !  idum
          unitd(icunit)= dact(unitin(icunit))*sum ! d(E)/d(w_ij)
        enddo   !  icunit
      enddo     !  ilay
c-Computes DWEIPAT=(grad(f)/dw) DWEIPAT, per pattern
      do ilay=2,nlay
        lastun = laypo(ilay+1)-1 ! last of the current layer, excl.bias
        if(ilay.ne.nlay) lastun=lastun-kk
        do icunit=laypo(ilay),lastun  ! current layer units
          do idum=1,layunits(ilay-1) ! previous layer units
            iunit = laypo(ilay-1)  + (idum-1)
            iwei  = listwei(icunit)+ (idum-1)
            dweipat(iwei) = unitd(icunit)*unitout(iunit)
          enddo   ! previous layer units
        enddo   ! current layer units
      enddo   ! current layer
c-Updates Weights, in ON-LINE training 
      if(monoff.eq.1) then
        if (kter.eq.1.or.ksolve.eq.1) then ! epoch_1, Steepest-Descent
            bfr=0.d0  ! for Fletcher-Reeves
            do iwei=1,nwei
               apold(iwei)=0.d0   !  old derivatives (zeroing needed?)
            enddo
        else  !    Polak-Ribiere  or  Fletcher-Reeves
            stop  ' NOT-READY ! '
            sumar_fr=0.d0
            sumar_pr=0.d0
            sumpar=0.d0
            do iwei=1,nwei
               sumar_pr=sumar_pr+dweipat(iwei)*
     :                  (dweipat(iwei)-dweiold(iwei))
               sumar_fr=sumar_fr+dweipat(iwei)*dweipat(iwei)
               sumpar=sumpar+dweiold(iwei)*dweiold(iwei)
            enddo
            sumpar=sumpar+1.d-10
            if(ksolve.eq.2) then
               bfr=sumar_pr/sumpar   !Polak Ribiere
c-Limiter     !bfr=dmax1(bfr,0.d0)   !constraint Haikin, page 240
            elseif(ksolve.eq.3) then
               bfr=sumar_fr/sumpar   !Fletcher Reeves
            endif
        endif
        do iwei=1,nwei
cccc      randompat=ran3(idum)   !  kirk, SOS
          ap=-dweipat(iwei)+bfr*apold(iwei)
          dweiold(iwei)=dweipat(iwei)
          apold(iwei)=ap
          wei(iwei)=wei(iwei)+eta*ap !   *randompat    !  kirk
        enddo
c-Collects Delta_Weights DWEI from all patterns, in OFF-LINE training 
      else   !  monoff.eq.2
        do iwei=1,nwei
          dwei(iwei)=dwei(iwei)+dweipat(iwei)
        enddo
      endif   !  monoff
c-Checks for convergence. Attention: After weight correction??
      dist    = 0.d0    !   cost for this pattern
      i=0
      do iunit=laypo(nlay), laypo(nlay+1)-1 ! except bias
        i = i+1
        deviation=db(idpat+kfree+i)-unitout(iunit)
        dist = dist+deviation*deviation
      enddo
        cost    = cost+0.5d0*dist   !  cost for all epoch's patterns
  66  continue        !   go for next DB pattern
c-Updates weights for Off-line correction, once for all patterns
      if(monoff.eq.2) then
        if (kter.eq.1.or.ksolve.eq.1) then ! epoch_1, Steepest-Descent
            bfr=0.d0   !!! Steepest descent for the first step!
            do iap=1,nwei
               apold(iap)=0.d0
            enddo
        else
            bfr=0.d0
            sumar_fr=0.d0
            sumar_pr=0.d0
            sumpar=0.d0
            do iwei=1,nwei
               sumar_pr=sumar_pr+dwei(iwei)*(dwei(iwei)-dweiold(iwei))
               sumar_fr=sumar_fr+dwei(iwei)*dwei(iwei)
               sumpar  =sumpar  +dweiold(iwei)*dweiold(iwei)
            enddo
      !!!      sumpar=sumpar+1.d-20
            if(ksolve.eq.2) then
               bfr=sumar_pr/sumpar   !Polak Ribiere
               bfr=dmax1(bfr,0.d0)   !constraint Haikin, page 240
            elseif(ksolve.eq.3) then
               bfr=sumar_fr/sumpar   !Fletcher Reeves
            endif
        endif
c-stores current ... during the search of optimal ETA
        do iwei=1,nwei
          weiold(iwei)=wei(iwei)
          apold(iwei)=-dwei(iwei)+bfr*apold(iwei)
          dweiold(iwei)=dwei(iwei)
        enddo
c
cc      write(*,*)kter,' bfr=',bfr
        costmin=1.d12
        ksample=0
        do intertest=10,60,5
          fw= ( dfloat(intertest)/10.d0  -  4.d0 ) / 2.d0
          etatest=eta*10**fw
          do iwei=1,nwei
            wei(iwei)=weiold(iwei)+etatest*apold(iwei)
          enddo
          costtest=0.d0   !  for this ETATEST value
          do npat=1,kpatterns
           idpat= (npat-1)*(kfree+nobj)
           do j=1,kfree
             xfstring(j) = db(idpat+j)
           enddo
           call use_perceptron(maxlay,layunits,laypo,kfree,nobj,
     1           nunit,maxunit,unitout,unitin,unitd,wei,listwei,
     2           maxwei,nwei,nlay,xfstring,kparam,mobj,kk)
           dist    = 0.d0    !   cost for this pattern
           i=0
           do iunit=laypo(nlay), laypo(nlay+1)-1 ! except bias
             i = i+1
             deviation=db(idpat+kfree+i)-unitout(iunit)
             dist = dist+deviation*deviation
           enddo
           costtest = costtest+0.5d0*dist
          enddo    !     npat
          ksample=ksample+1
          xsample(ksample)=etatest
          ysample(ksample)=costtest
          if(costtest.lt.costmin)then
            costmin=costtest
            etabest=etatest
            kksample=ksample
          endif
        enddo    !     intertest
c-find best eta, by fitting polynome
        goto 667  !  skip
        if(kksample.eq.1)          kksample=2
        if(kksample.eq.ksample)    kksample=ksample-1
        asample(1,1)=xsample(kksample-1)**2
        asample(1,2)=xsample(kksample-1)
        asample(1,3)=1.d0
          rsample(1)=ysample(kksample-1)
        asample(2,1)=xsample(kksample)**2
        asample(2,2)=xsample(kksample)
        asample(2,3)=1.d0
          rsample(2)=ysample(kksample)
        asample(3,1)=xsample(kksample+1)**2
        asample(3,2)=xsample(kksample+1)
        asample(3,3)=1.d0
          rsample(3)=ysample(kksample+1)
        call gauss(3,3,asample,rsample,ssample)
        etabest=-ssample(2)/ssample(1)/2.d0
        if(etabest.lt.1.d-2*eta) etabest=1.d-2*eta
        if(etabest.gt.1.d+2*eta) etabest=1.d+2*eta
cc      write(*,*)' #####  ETA_best        ',etabest
  667   continue     !  skip
        
c-repeat computation for best ETA
          do iwei=1,nwei
            wei(iwei)=weiold(iwei)+etabest*apold(iwei)
          enddo
c-kirk    eta=etabest       !   SOS  kirk,  too risky !!!!
      endif
c-Computes and prints out cost, etc (cost prior to this correction)
      cost    = dsqrt(cost   /dfloat(kpatterns))
      if(dlog10(cost).lt.-10.d0) then
        write(*,'(i10,2x,f10.5)') kter,dlog10(cost+1.d-16)
        write(12,'(i10,2(2x,f10.5))') kter,dlog10(cost+1.d-16),eta
        call Save_ANNstate(nlay,nunit,nwei,wei,maxwei,kter)
        write(*,*)
        write(*,*)' ################################################# '
        write(*,*)' #####  Training converged to machine accuracy ! '
        write(*,*)' #####  log10(cost) = ',dlog10(cost)
        write(*,*)' #####  after ',kter,' epochs '
        call mess_ksolve(ksolve)
        write(*,*)' #####  with ',kpatterns,' patterns '
        write(*,*)' #####  and  ',kfree,' design variables '
        write(*,*)' ################################################# '
        return
      endif
c     if(kter.lt.10000000.or.iter.eq.maxiter) goto 13
c     if(mod(kter,10)   .eq.0.and.kter.le.100   ) goto 13
      if(mod(kter,50)   .eq.0                   ) goto 13
c     if(mod(kter,100)  .eq.0.and.kter.le.1000  ) goto 13
c     if(mod(kter,500)  .eq.0.and.kter.le.10000 ) goto 13
c     if(mod(kter,5000) .eq.0) goto 13
      goto 101   !  skips writing convergence
  13  write(* ,'(i10,3(2x,f10.5))') kter,dlog10(cost+1.d-16),eta
      write(12,'(i10,3(2x,f10.5))') kter,dlog10(cost+1.d-16),eta

c-Change a weight
      goto 101
 102  if(iter.eq.maxiter)then
         write(*,*)' Enter a NON-ZERO number to randomize a weight'
         read(*,*)irawei
         if(irawei.ne.0)then
            randompar=ran3(idum)   !  kirk, SOS
            iwei=1+(nwei-1)*randompar
            bbb      =ran3(idum)*.1*(-1)**iwei
            write(*,*)' Change',iwei,bbb,wei(iwei)
            wei(iwei)=bbb
         endif
      endif

      
 101  continue        !    go for next cycle
c
      if(interactive.ne.0) then
        write(*,'(a,i10)')'  #####  Epoch      =',kter
        write(*,'(a,f12.5,a,i3)')'  #####  Current eta=',eta,' bias=',kk
        call mess_ksolve(ksolve)
        call Save_ANNstate(nlay,nunit,nwei,wei,maxwei,kter)
        write(*,*)' >> '
        write(*,*)' >>  ENTER Number of extra iterations &  eta  : '
        write(*,*)'  (-ve number of iterations being shuffling)'
        write(*,*)' >> --->>>  '
        read(*,*) maxiter,eta 
        if(maxiter.lt.0)then
            randompar=ran3(idum)   !  kirk, SOS
            iwei=1+(nwei-1)*randompar
            bbb      =ran3(idum)*.1*(-1)**iwei
            write(*,*)' Change',iwei,bbb,wei(iwei)
            wei(iwei)=bbb
            maxiter=-maxiter
        endif
        if(maxiter.gt.0) goto 100
      endif
      close(12)
c
      return
      end
c
c **********************************************************************
      subroutine mess_ksolve(ksolve)
c **********************************************************************
      if(ksolve.eq.1)then
         write(*,*)' #####  Running with Steepest-Descent   '
      elseif(ksolve.eq.2)then
         write(*,*)' #####  Running with Polak-Ribiere      '
      elseif(ksolve.eq.3)then
         write(*,*)' #####  Running with Fletcher-Reeves    '
      endif
      return
      end
c
c **********************************************************************
      subroutine Save_ANNstate(nlay,nunit,nwei,wei,maxwei,kter)
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension wei(maxwei)
      write(*,*)' ##########   Writes file WEIGHTS.DAT'
      open(1,file='weights.dat')
      write(1,'(10i10)') nlay,nunit,nwei,kter
      write(1,*) (wei(i),i=1,nwei)
      close(1)
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
c
c **********************************************************************
      subroutine gauss(kdim,n,a,b,x)
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension a(kdim,kdim),b(kdim),x(kdim)
      eps=1.d-14
      do k=1,n-1
        if (dabs(a(k,k)).lt.eps)   stop   'Divide by Zero !!!!'
        do i=k+1,n 
           factor=a(i,k)/a(k,k)
           do j=k+1,n
           a(i,j)=a(i,j)-factor*a(k,j)
           enddo
           b(i)  =b(i)  -factor*b(k)
        enddo
      enddo
c
c-back substitution
      x(n)=b(n)/a(n,n)
      do i=n-1,1,-1
         sum=0.d0
         do j=i+1,n
           sum=sum+a(i,j)*x(j)
         enddo
         x(i)=(b(i)-sum)/a(i,i)
      enddo
c
      return
      end

      include 'recipes_ga.for'
