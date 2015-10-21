	program Force

	implicit none  
	real*8 x(300000),y(300000),z(300000),NaveF(100),NaveN(100)
	real*8 alpha(300000),delta(300000)
	real*8 r(300000),rmax,rmin,dr,rr(100) 
	real*8 dx,dy,dz,pi,V(1000) 
	real*8 N(100,300000),M(100,300000),NF(100,300000),NN(100,300000)
	real*8 Nave(100),Nvar(100),NvarN(100),NvarF(100)
	real*8 r1,r2,r3,r4,r5,r6,pp0
	real*8 rborder1(300000),rborder2(300000)
	real*8 xfrac,qq,rborder2A,rborder2B
	real*8 q,r_min,r_max,ar_min,ar_max,dec_min,dec_max
	real*8 rkmin,rkmax,dk,rk(0:50),rlambda,drk,rrmin,rrminnew 
	real*8 rd,rp,tstar,h1,h2(10000)
	real*8 xA,yA,zA,deltaA,alphaA,rA,tstarA
	real*8 xB,yB,zB,deltaB,alphaB,rB,tstarB
	
	integer i,nd,imax,j,k,nc,ichoice,ichoice2,JJ,kmax,ncN,ncF
	integer ichoice3,kmax2,l,ik 
	
	character*40 filein1,filein2,filein3
	character*40 fileout1,fileout2,fileout3,fileout4



	

	rkmax=0
	rlambda=1. 
	rmax=-10000
	rmin=10000
	pi=2.*asin(1.)

c	rborder_max=-1000
c	write(*,*)'imax='
c	read(*,*)imax

	write(*,*)'random sampling 0<p<1 ='
	read(*,*)pp0

	write(*,*)'radius ='
	read(*,*)xfrac

	imax=50

	open(30,file='fileinput.info',status='unknown') 
	open(99,file='adr.dat',status='unknown') 

c	open(98,file='rborder.dat',status='unknown') 
	
	read(30,*)r_min
	read(30,*)r_max
	read(30,*)ar_min
	read(30,*)ar_max
	read(30,*)dec_min
	read(30,*)dec_max
	read(30,*)filein1
	read(30,*)filein2
	read(30,*)filein3
	read(30,*)fileout1
	read(30,*)fileout2
	read(30,*)fileout3
	read(30,*)fileout4
	
	write(*,*)r_min
	write(*,*)r_max
	write(*,*)ar_min
	write(*,*)ar_max
	write(*,*)dec_min
	write(*,*)dec_max
	write(*,*)filein1
c	write(*,*)filein2
c	write(*,*)filein3
c	write(*,*)fileout1
c	write(*,*)fileout2
c	write(*,*)fileout3
c	write(*,*)fileout4
	
	open(10,file=filein1,status='old') 
	open(20,file=fileout1,status='unknown') 
	

	   ar_min=ar_min*pi/180.
	   ar_max=ar_max*pi/180.
	   dec_min=dec_min*pi/180.
           dec_max=dec_max*pi/180.

c**************************************************************
c       READ FILE
c**************************************************************

	   write(20,*)'# r(1),Nave(2),Nvar(3),nave(4),nave(5),nvar(6)
	1,nc(7),NaveF(8),Navan(9),delta(10),vardelta(11)'

c**************************************************************


	   i=1 
 12	   continue 
 15	   read(10,*,end=11)x(i),y(i),z(i) 
	   if(rand(0).gt.pp0)goto 77

c*********************Transformation into spherical coordinates *******
	
	r(i)=sqrt(x(i)**2.+y(i)**2.+z(i)**2.)
	if(r(i).lt.r_min)goto 77
	if(r(i).gt.r_max)goto 77
	
	if(r(i).lt.rmin)rmin=r(i)
	if(r(i).gt.rmax)rmax=r(i)

         IF(X(i).EQ.0.) THEN
            IF(Y(i).GT.0)ALPHA(i)=PI/2.
            IF(Y(i).LT.0)ALPHA(i)=3*PI/2.
            GOTO 111
         END IF

         IF(Y(i).EQ.0)THEN
            IF(X(i).GT.0)ALPHA(i)=0
            IF(X(i).LT.0)ALPHA(i)=PI
            GOTO 111
         END IF

         IF(X(i).GT.0.AND.Y(i).GT.0)Q=0
         IF(X(i).LT.0.AND.Y(i).GT.0)Q=2
         IF(X(i).LT.0.AND.Y(i).LT.0)Q=2
         IF(X(i).GT.0.AND.Y(i).LT.0)Q=4
         
         alpha(i)=atan(y(i)/x(i)) + PI/2. *Q 
 111      CONTINUE

         IF(X(i).EQ.0.AND.Y(i).EQ.0)THEN
            IF(Z(i).GT.0)DELTA(i)=PI/2.
            IF(Z(i).LT.0)DELTA(i)=-PI/2.
            GOTO 112
            END IF
            

         delta(i)=asin(z(i)/r(i))
 112      CONTINUE 


c******************************************************************	  
	  if(alpha(i).gt.pi)alpha(i)=alpha(i)-2.*pi

	  if(delta(i).lt.dec_min)goto 77
	  if(delta(i).gt.dec_max)goto 77
	  if(alpha(i).lt.ar_min)goto 77
	  if(alpha(i).gt.ar_max)goto 77

	write(99,*)alpha(i)*180./pi,delta(i)*180./pi,r(i)
c	write(98,*)i,rborder1(i),rborder2(i) 


	   i=i+1 
 77	   continue
	   goto 12 

 11	   close(10) 
	   close(99)
	   nd=i-1 
	   write(*,*)'nd,rmax,rmin =',nd,rmax,rmin 
c	   read(*,*)
	   



c**************************************************************
c       COMPUTING GAMMA 
c**************************************************************

c******** parameters of the r-range

	   rmin=0.1 
	   rmax=200. 

	   qq=10.**(log(rmax/rmin)/log(10.)*1./imax)
	   write(*,*)'imax =',imax,qq

c	   read(*,*)
		 
C*******DO LOOP ON THE r VALUES TO COMPUTE GAMMA 

c       rr(i) is the hight of the cylinder

c       h2 is the cylinder base radius


	   do i=1,imax 
	      rr(i)=rmin*qq**i 
	      write(*,*)rr(i),i
	      h2(i)=xfrac 
c	      h2(i)=rr(i)/xfrac 
c	      h1=sqrt(rr(i)**2.+h2(i)**2.)

C*******1ST    DO LOOP ON # OF POINTS
	      
	      do j=1,nd 
		 if(mod(j,5000).eq.0) write(*,*)j
c		 write(*,*)j,rr(i),rborder(j) 
		 M(i,j)=0.

c*****  BOUNDARY CONDITIONS CYLINDER

c***************** center 	  
    
	  r1=abs(r(j)*cos(delta(j))*sin(alpha(j)-ar_min))
	  r2=abs(r(j)*cos(delta(j))*sin(ar_max-alpha(j)))
	  r3=abs(r(j)*sin(delta(j)-dec_min))
	  r4=abs(r(j)*sin(dec_max-delta(j)))
 	  
	  r5=r(j)-r_min
	  r6=r_max-r(j)

	  rborder1(j)=r6
	  if(r5.lt.rborder1(j))rborder1(j)=r5
	  
	  rborder2(j)=r1
	  if(r2.lt.rborder2(j))rborder2(j)=r2
	  if(r3.lt.rborder2(j))rborder2(j)=r3
	  if(r4.lt.rborder2(j))rborder2(j)=r4

c************************* POINT A

		 tstarA=1.+rr(i)/r(j)
		 xA=x(j)*tstarA
		 yA=y(j)*tstarA
		 zA=z(j)*tstarA


	rA=sqrt(xA**2.+yA**2.+zA**2.)


         IF(XA.EQ.0.) THEN
            IF(YA.GT.0)ALPHAA=PI/2.
            IF(YA.LT.0)ALPHAA=3*PI/2.
            GOTO 111
         END IF

         IF(YA.EQ.0)THEN
            IF(XA.GT.0)ALPHAA=0
            IF(XA.LT.0)ALPHAA=PI
            GOTO 611
         END IF

         IF(XA.GT.0.AND.YA.GT.0)Q=0
         IF(XA.LT.0.AND.YA.GT.0)Q=2
         IF(XA.LT.0.AND.YA.LT.0)Q=2
         IF(XA.GT.0.AND.YA.LT.0)Q=4
         
         alphaA=atan(yA/xA) + PI/2. *Q 
 611      CONTINUE

         IF(XA.EQ.0.AND.YA.EQ.0)THEN
            IF(ZA.GT.0)DELTAA=PI/2.
            IF(ZA.LT.0)DELTAA=-PI/2.
            GOTO 612
            END IF
            

         deltaA=asin(zA/rA)
 612      CONTINUE 

	  if(alphaA.gt.pi)alphaA=alphaA-2.*pi

	  if(alphaA.gt.ar_max)goto 13
	  if(alphaA.lt.ar_min)goto 13
	  if(deltaA.gt.dec_max)goto 13
	  if(deltaA.lt.dec_min)goto 13
	  if(rA.gt.r_max)goto 13
	  if(rA.lt.r_min)goto 13


	  r1=abs(rA*cos(deltaA)*sin(alphaA-ar_min))
	  r2=abs(rA*cos(deltaA)*sin(ar_max-alphaA))
	  r3=abs(rA*sin(deltaA-dec_min))
	  r4=abs(rA*sin(dec_max-deltaA))
 	  	  
	  rborder2A=r1
	  if(r2.lt.rborder2A)rborder2A=r2
	  if(r3.lt.rborder2A)rborder2A=r3
	  if(r4.lt.rborder2A)rborder2A=r4




c************************* POINT B

		 tstarB=1.-rr(i)/r(j)
		 xB=x(j)*tstarB
		 yB=y(j)*tstarB
		 zB=z(j)*tstarB


	rB=sqrt(xB**2.+yB**2.+zB**2.)


         IF(XB.EQ.0.) THEN
            IF(YB.GT.0)ALPHAB=PI/2.
            IF(YB.LT.0)ALPHAB=3*PI/2.
            GOTO 711
         END IF

         IF(YB.EQ.0)THEN
            IF(XB.GT.0)ALPHAB=0
            IF(XB.LT.0)ALPHAB=PI
            GOTO 711
         END IF

         IF(XB.GT.0.AND.YB.GT.0)Q=0
         IF(XB.LT.0.AND.YB.GT.0)Q=2
         IF(XB.LT.0.AND.YB.LT.0)Q=2
         IF(XB.GT.0.AND.YB.LT.0)Q=4
         
         alphaB=atan(yB/xB) + PI/2. *Q 
 711      CONTINUE

         IF(XB.EQ.0.AND.YB.EQ.0)THEN
            IF(ZB.GT.0)DELTAB=PI/2.
            IF(ZB.LT.0)DELTAB=-PI/2.
            GOTO 712
            END IF
            

         deltaB=asin(zB/rB)
 712      CONTINUE 

 	  if(alphaB.gt.pi)alphaB=alphaB-2.*pi

	  if(alphaB.gt.ar_max)goto 13
	  if(alphaB.lt.ar_min)goto 13
	  if(deltaB.gt.dec_max)goto 13
	  if(deltaB.lt.dec_min)goto 13
	  if(rB.gt.r_max)goto 13
	  if(rB.lt.r_min)goto 13

	  r1=abs(rB*cos(deltaB)*sin(alphaB-ar_min))
	  r2=abs(rB*cos(deltaB)*sin(ar_max-alphaB))
	  r3=abs(rB*sin(deltaB-dec_min))
	  r4=abs(rB*sin(dec_max-deltaB))
 	  	  
	  rborder2B=r1
	  if(r2.lt.rborder2B)rborder2B=r2
	  if(r3.lt.rborder2B)rborder2B=r3
	  if(r4.lt.rborder2B)rborder2B=r4



c****************** bc for the three points **** begin 

	      if(rborder1(j).lt.rr(i))goto 13
	      if(rborder2(j).lt.h2(i))goto 13
	      if(rborder2A.lt.h2(i))goto 13
	      if(rborder2B.lt.h2(i))goto 13

c****************** bc for the three points  **** end 

		 M(i,j)=1.	!LABEL FOR THE JTH POINT AT THE SCALE ITH
		 N(i,j)=0.
		 NF(i,j)=0.
		 NN(i,j)=0.


C*******2ND    DO LOOP ON # OF POINTS
		 
		 do k=1,nd
		    if(k.eq.j)goto 88
		    dx=x(j)-x(k)
		    dy=y(j)-y(k)
		    dz=z(j)-z(k)
		    rd=sqrt(dx**2.+dy**2.+dz**2.)

		    tstar=(x(j)*x(k)+y(j)*y(k)+z(j)*z(k))/
     & (x(j)**2.+y(j)**2.+z(j)**2.)
		    rp=sqrt( (x(k)-tstar*x(j))**2. + 
     & (y(k)-tstar*y(j))**2. + (z(k)-tstar*z(j))**2. ) 

		    if(rp.le.h2(i))then
		       if(sqrt(rd**2-rp**2).le.rr(i))then
			  N(i,j)=N(i,j)+1. 
			  if(r(k).gt.r(j))NF(i,j)=NF(i,j)+1. 
			  if(r(k).lt.r(j))NN(i,j)=NN(i,j)+1. 
c			  if(j.eq.10)write(100+i,*)x(k),y(k),z(k)
		       end if
		    end if
 88		    continue 
		 end do		!2nd (k)
 13		 continue 
	      end do		!1st (j) 

	   end do		!scales (i) 

	   write(*,*)'do loops on r, and Np2 finished'

c*************************************************************
c*******Gamma and its variance 

	   do i=1,imax

	      Nave(i)=0.
	      Nvar(i)=0.

	      NaveF(i)=0.
	      NaveN(i)=0. 

	      V(i)=(2.*rr(i))*pi*h2(i)**2. 
c       (4.*pi/3.)*rr(i)**3. 
	      nc=0
	      ncN=0
	      ncF=0

	      do j=1,nd
		 if(M(i,j).eq.1)then
		    Nave(i)=Nave(i)+N(i,j)
                    if(NF(i,j).gt.0)then
		       NaveF(i)=NaveF(i)+NF(i,j)
		       ncF=ncF+1
		       end if
		       if(NN(i,j).gt.0)then
			  NaveN(i)=NaveN(i)+NN(i,j)
			  ncN=NcN+1
		       end if
		    nc=nc+1
		 end if
	      end do 

	      Nave(i)=Nave(i)/nc
	      NaveN(i)=NaveN(i)/nc
	      NaveF(i)=NaveF(i)/nc


	      do j=1,nd
		 if(M(i,j).eq.1)then
		       Nvar(i)=Nvar(i)+(N(i,j)-Nave(i))**2.
		       NvarF(i)=NvarF(i)+(NF(i,j)-NaveF(i))**2.
		       NvarN(i)=NvarN(i)+(NN(i,j)-NaveN(i))**2.
		 end if
	      end do
	      


		 if(nc.gt.1)then
	      Nvar(i)=sqrt(Nvar(i)/(float(nc)*(float(nc)-1.)))
	      NvarF(i)=sqrt(NvarF(i)/(float(nc)*(float(nc)-1.)))
	      NvarN(i)=sqrt(NvarN(i)/(float(nc)*(float(nc)-1.)))

	      if(ncF.gt.0.and.ncN.gt.0)then
		 write(20,'(11(e12.5,2x))')rr(i),Nave(i),Nvar(i),
	1	      Nave(i)/V(i),Nave(i)/V(i),Nvar(i)/V(i),float(nc)
	1	      ,NaveF(i) ,NaveN(i),(NaveN(i)-NaveF(i)) /(NaveN(i)
	2	      +NaveF(i)), sqrt(NvarN(i)**2+NvarF(i)**2)
	2	      /(NaveN(i)+NaveF(i))
	      end if
	      
		    rkmax=rr(i)
		 end if
 17		 continue
	   end do
	   close(20)
c***************************************************************	   
c***************************************************************	   
c       SL 
	   jj=0
	   write(*,*)'rkmax =',rkmax 

c***** scale at which writes the outputs

	   rk(0)=0
	   rk(1)=2.
	   rk(2)=5.
	   rk(3)=10.
	   rk(4)=20.
	   rk(5)=30.
	   rk(6)=40.
	   rk(7)=50.
	   rk(8)=60.
	   rk(9)=70.
	   rk(10)=80.
	   rk(11)=90.
	   rk(12)=100.
	   rk(13)=110.
	   rk(14)=120.
	   rk(15)=130.
	   rk(16)=140.
	   rk(17)=150.
	   rk(18)=160.
	   rk(19)=170.
	   rk(20)=180.
	   rk(21)=190.
	   rk(22)=200.
	   rk(23)=250.
	   rk(24)=300.
	   rk(25)=350.
	   rk(26)=400.
	   rk(27)=450.
	   rk(28)=500.
	   rk(29)=550.
	   rk(30)=600.
	   rk(31)=650.


	   kmax=31 

	   do k=1,kmax
	      if(rk(k).le.rkmax)kmax2=k
	   end do

	   write(*,*)'kmax2, max scale =',kmax2,rk(kmax2)

c	   kmax2=kmax 



c******* writes the output files

	   
	   do k=1,kmax2
	      
      if(k.eq.1) open(15,file='SL_CP002.dat',status='unknown') 
      if(k.eq.2) open(15,file='SL_CP005.dat',status='unknown') 
      if(k.eq.3) open(15,file='SL_CP010.dat',status='unknown') 
      if(k.eq.4) open(15,file='SL_CP020.dat',status='unknown') 
      if(k.eq.5) open(15,file='SL_CP030.dat',status='unknown') 
      if(k.eq.6) open(15,file='SL_CP040.dat',status='unknown') 
      if(k.eq.7) open(15,file='SL_CP050.dat',status='unknown') 
      if(k.eq.8) open(15,file='SL_CP060.dat',status='unknown') 
      if(k.eq.9) open(15,file='SL_CP070.dat',status='unknown') 
      if(k.eq.10)open(15,file='SL_CP080.dat',status='unknown') 
      if(k.eq.11)open(15,file='SL_CP090.dat',status='unknown') 
      if(k.eq.12)open(15,file='SL_CP100.dat',status='unknown') 
      if(k.eq.13)open(15,file='SL_CP110.dat',status='unknown') 
      if(k.eq.14)open(15,file='SL_CP120.dat',status='unknown') 
      if(k.eq.15)open(15,file='SL_CP130.dat',status='unknown') 
      if(k.eq.16)open(15,file='SL_CP140.dat',status='unknown') 
      if(k.eq.17)open(15,file='SL_CP150.dat',status='unknown') 
      if(k.eq.18)open(15,file='SL_CP160.dat',status='unknown') 
      if(k.eq.19)open(15,file='SL_CP170.dat',status='unknown') 
      if(k.eq.20)open(15,file='SL_CP180.dat',status='unknown') 
      if(k.eq.21)open(15,file='SL_CP190.dat',status='unknown') 
      if(k.eq.22)open(15,file='SL_CP200.dat',status='unknown') 
      if(k.eq.23)open(15,file='SL_CP250.dat',status='unknown') 
      if(k.eq.24)open(15,file='SL_CP300.dat',status='unknown') 	      
      if(k.eq.25)open(15,file='SL_CP350.dat',status='unknown') 
      if(k.eq.26)open(15,file='SL_CP400.dat',status='unknown') 
      if(k.eq.27)open(15,file='SL_CP450.dat',status='unknown') 
      if(k.eq.28)open(15,file='SL_CP500.dat',status='unknown') 
      if(k.eq.29)open(15,file='SL_CP550.dat',status='unknown') 
      if(k.eq.30)open(15,file='SL_CP600.dat',status='unknown') 
      if(k.eq.31)open(15,file='SL_CP650.dat',status='unknown') 


      write(15,*)'#r(1),N(2),x(3),y(4),z(5),rr(5),al(7),del(8)'

	      
	      rrmin=100000.
	      do i=1,imax
		 rrminnew=abs(rr(i)-rk(k))
		 if(rrminnew.lt.rrmin)then
		    rrmin=rrminnew
		    l=i
		 end if
	      end do
	      write(*,*)rk(k),rr(l)
	      
	      
	      
	      do j=1,nd
		 if(M(l,j).eq.1)then
		    write(15,'(9(e12.5,2x))')r(j),N(l,j),
     %    x(j),y(j),z(j),rr(l),alpha(j),delta(j)
		 end if
	      end do		!points
	      
	      
 68	      continue
	   end do		!scales where writes on files 
	   

c********************************************************
 1111	   continue 

	   end  
	
 
