C---A Fortran code for the maximum likelihood estimation (MLE) for the two star ERG model of a bipartite network  
C---Bipartite network: Example bank firm network
C-- Input Parameters....
C--N: No of firms, M: No of banks, L: Total number of links, ni= number
C of parameter to be estimated. ni=3 for two star model
C-- nconf: No of ensamble avarage------   
C---iseeds are the intital random number seed....

C---output will be written in "theta_2Star.dat". col1: no of
C-- ensamble avg, col2: no of converge ensamble avg, col3: estimated
C----parameter theta, col4: <theta^2>, col5: std of theta 


      implicit double precision(a-h,o-z)
      parameter (N=10,M=5,L1=17,kmax=M,kmax1=N)
      parameter (ni=3,nconf=1000,maxrun1=7+3*ni)
      parameter (iseed2=1449942249,iseed1=76438761,iseed3=387489371)
      integer  Pneb(N,kmax),Pdegree(N),Pneb1(N,kmax),Pdegree1(N)
      integer  Adegree(M),Adegree1(M),Aneb(M,kmax1),Aneb1(M,kmax1)
      dimension theta(ni),z(ni),u(ni),zo(ni),zsum(ni),zsqsum(ni)
      dimension  D1inv(ni),f3(ni),thetasum(ni),thetasqsum(ni)
      dimension zf(ni)

      iran1=iseed1
      iran2=iseed2
      iran3=iseed3
     
      Pdegree=0
      Adegree=0
         
C-----------Observed network--------------------------------------
C-----------Read the linklist data of observed network and creat the degree array and neigbour list of each node--------------------------------------
C----Link list data does not cotain any parallel edges-------
C---'linklist.dat' is a two columns data file. col1: firm index{1,N}, colo2: bank index {1,M}-------
     
      open(1,file="linklist.dat")
      do 11 iL=1,L
      read(1,*)i,j
      Pdegree(i)=Pdegree(i)+1
      Adegree(j)=Adegree(j)+1
      Pneb(i,Pdegree(i))=j
      Aneb(j,Adegree(j))=i
 11   enddo
      close(1)
C----------Calculate network statistics in observed network zo----

c-----intialize
      zo=0
      
C----Edge zo1
      zo(1)=L
       
C--------Two star in Bank node set
      do i=1,M
      k=Adegree(i)
      zo(2)=zo(2)+ (dfloat(k*(k-1)))/2
      enddo
C--------Two star in firm node set
      do i=1,N
      k=Pdegree(i)
      zo(3)=zo(3)+ (dfloat(k*(k-1)))/2
      enddo
      


C-------copy of the observed network---------------------------------------
      Pneb1=Pneb
      Pdegree1=Pdegree
      Aneb1=Aneb
      Adegree1=Adegree
      z=zo




      thetasum=0.0d00
      thetasqsum=0.0d00
      kount=0
      zf=0.0d00
     
      
      do 5000 iconf=1,nconf
     
c-------Phase 1------------------------------------------

c--To determine Scaling Matrix D0 and initial parameters

      theta=0.0d00
      zsum=0.d00
      zsqsum=0.0d00
      a=0.10d00


      do irun1=1,maxrun1

C-----Starting from the observed network generate maxrun1=(7+3*ni) new network using MCMC------ 
      Pneb=Pneb1
      Pdegree=Pdegree1
      Aneb=Aneb1
      Adegree=Adegree1
      z=zo


      call MCMC(N,M,ni,Pneb,Pdegree,Aneb,Adegree,z,theta,
     .iran1,iran2,iran3)

C-----average stattics of the generated network------------------
      do ir=1,ni
      zsum(ir)=zsum(ir)+z(ir)
      zsqsum(ir)=zsqsum(ir)+z(ir)*z(ir)
      enddo       
      
      enddo

      count=maxrun1

C----Scaling matrix D1-------------------
      do ir=1,ni
      av=zsum(ir)/count
      avsq=zsqsum(ir)/count
      D1=avsq-av*av
      D1inv(ir)=1.0d00/D1
      theta(ir)=theta(ir)-a*D1inv(ir)*(av-zo(ir))
      enddo
C------------Phase 2---------------------------
      
      itime=0
      isub=0
           
 17   it=0
      isub=isub+1
      a=0.40d00
      zsum=0.d00
      zsqsum=0.0d00
      count=0.0d00

 18   it=it+1
      imt=(7+ni)*(2**(1.333d00*it))+200
      do 1000 im=1,imt
C------Start from the observed network-----------------
      itime=itime+1
      Pneb=Pneb1
      Pdegree=Pdegree1
      Aneb=Aneb1
      Adegree=Adegree1
      z=zo
C-----Generating one graph and statistics Z from the model

      call MCMC(N,M,ni,Pneb,Pdegree,Aneb,Adegree,z,theta,
     .iran1,iran2,iran3)

c---------------------------------
C------Optimization: Newton Raphson minimization-----
      do ir=1,ni
      theta(ir)=theta(ir)-a*D1inv(ir)*(z(ir)-zo(ir)) 
      enddo
     



      if(it.gt.4)then
      do ir=1,ni
      zsum(ir)=zsum(ir)+z(ir)
      zsqsum(ir)=zsqsum(ir)+z(ir)*z(ir)
      enddo
  
      
      
      count=count+1
      ic=0
      do ir=1,ni
      f1=zsum(ir)/count
      f2=dsqrt(zsqsum(ir)/count-f1*f1)
      f3(ir)=(f1-zo(ir))/f2
      if(dabs(f3(ir)).gt.0.1d00)ic=1
      enddo

      if(ic.lt.1)then
      kount=kount+1
 
      open(1,file='theta_2Star.dat')
      do ir=1,ni

      thetasum(ir)=thetasum(ir)+theta(ir)
      thetasqsum(ir)=thetasqsum(ir)+theta(ir)*theta(ir)
      std=dsqrt((thetasqsum(ir)/kount)-((thetasum(ir)/kount)
     .*(thetasum(ir)/kount)))
 
      zf(ir)=zf(ir)+zsum(ir)/count

      write(1,*)iconf,kount,thetasum(ir)/kount,thetasqsum(ir)/kount,std
c     .zo(ir),zf(ir)/kount
      enddo
      close(1)
      goto 5000
      endif
      endif

 1000 enddo
       a=a/2.0d00
       if(isub.gt.50.and.it.eq.6)goto 5000
       if(it.lt.6)goto 18
       if(it.gt.5)goto 17
 5000 enddo
      end

      subroutine MCMC(N,M,ni,Pneb,Pdegree,Aneb,Adegree,z,theta,
     .iran1,iran2,iran3)
      implicit double precision(a-h,o-z)
C----set maxrun= 10*N*M---------

      parameter (maxrun=10*10*5)
      integer  Pneb(N,M),Pdegree(N)
      integer  Adegree(M),Aneb(M,N) 
      dimension u(ni),theta(ni),z(ni)
      dimension Fr1(N),Fr2(N),Fr3(N),Fr4(N),Fr5(N)
      dimension Fr6(N),Fr7(N),Fr8(N),Fr9(N)
      dimension Br1(M),Br2(M),Br3(M),Br4(M),Br5(M)


      

    
      do 3000 irun=1,maxrun

 20   i=1+ranf1(iran1)*N
      if(i.gt.N)goto 20
 21   j=1+ranf2(iran2)*M
      if(j.gt.M)goto 21
     
      kj=Adegree(j)
      ki=Pdegree(i)
      do k=1,ki
      if(Pneb(i,k).eq.j)goto 333
      enddo
C-----calculation of change statistics if link added------
      u=0.0d00
      u(1)=1
      u(2)=kj
      u(3)=ki
      
      Hr=0.0d00
      do ir=1,ni
      Hr=Hr+theta(ir)*u(ir)
      enddo

      if(ranf3(iran3).lt.dexp(Hr))then
      do ir=1,ni
      z(ir)=z(ir)+u(ir)
      enddo
      Pdegree(i)=Pdegree(i)+1
      Adegree(j)=Adegree(j)+1
      Pneb(i,Pdegree(i))=j
      Aneb(j,Adegree(j))=i
c------------------------------------
      endif
      goto 50
C-----calculation of change statistics if link deleted------

 333  u=0.0d00
      u(1)=-1.0d00
      u(2)=1.0d00-kj
      u(3)=1.0d00-ki


      Hr=0.0d00
      do ir=1,ni
      Hr=Hr+theta(ir)*u(ir)
      enddo


      if(ranf3(iran3).lt.dexp(Hr))then
      isuc=isuc+1
      do ir=1,ni
      z(ir)=z(ir)+u(ir)
      enddo
      i1=k
      do k=1,kj
      if(Aneb(j,k).eq.i)goto 444
      enddo
 444  j1=k
      ineb=Pneb(i,ki)
      Pdegree(i)=Pdegree(i)-1
      if(j.ne.ineb)Pneb(i,i1)=ineb

      jneb=Aneb(j,kj)
      Adegree(j)=Adegree(j)-1
      if(i.ne.jneb)Aneb(j,j1)=jneb

c------------------------------------
   

      endif

 50   continue

 3000 enddo
      
      end
      
      
      
      


*----------------------------------------------------------------------------------
*     RANDOM NUMBER GENERATOR 1
*----------------------------------------------------------------------------------

      double precision function ranf1(iran1)
      iran1=iran1*1566083941
      if(iran1.lt.0)iran1=iran1+2147483647+1
      iran1=iran1*1566083941
      if(iran1.lt.0)iran1=iran1+2147483647+1
      ranf1=iran1*4.6566128752458D-10
      return
      end

*----------------------------------------------------------------------------------
*     RANDOM NUMBER GENERATOR 2
*----------------------------------------------------------------------------------

      double precision function ranf2(iran2)
      iran2=iran2*1664525
      if(iran2.lt.0)iran2=iran2+2147483647+1
      iran2=iran2*1664525
      if(iran2.lt.0)iran2=iran2+2147483647+1
      ranf2=iran2*4.6566128752458D-10
      return
      end

*----------------------------------------------------------------------------------
*     RANDOM NUMBER GENERATOR 3
*----------------------------------------------------------------------------------

      double precision function ranf3(iran3)
      iran3=iran3*16807
      if(iran3.lt.0)iran3=iran3+2147483647+1
      iran3=iran3*16807
      if(iran3.lt.0)iran3=iran3+2147483647+1
      ranf3=iran3*4.6566128752458D-10
      return
      end
