	program Force

	implicit none  
	real*8 x(300000),y(300000),z(300000),NaveF(100),NaveN(100)
	real*8 alpha(300000),delta(300000)
	real*8 r(300000),rmax,rmin,dr,rr(100) 
	real*8 dx,dy,dz,pi,V(1000) 
	real*8 N(100,300000),M(100,300000),NF(100,300000),NN(100,300000)
	real*8 Fx(100,300000),Fy(100,300000),Fz(100,300000)
	real*8 Nave(100),Nvar(100),F(100,300000),pp0,NvarN(100)
	1    ,NvarF(100)
	real*8 rrmax,qq 
	real*8 r1,r2,r3,r4,r5,r6,rborder_max

	real*8 rb3P(300000),rb2P(300000),rb1P(300000)
	real*8 rb1P3(300000),rb1P2(300000)
	real*8 rb2P3(300000),rb2P2(300000)

	real*8 q,r_min,r_max,ar_min,ar_max,dec_min,dec_max
	real*8 rkmin,rkmax,dk,rk(0:50),rlambda,drk,rrmin,rrminnew 
	real*8 rd,rp,tstar,h1,h2(10000),omega,xfrac
	integer i,nd,imax,j,k,nc,ichoice,ichoice2,JJ,kmax
	integer ichoice3,kmax2,l,ik,jjk,ii,id,Nkt,ncN,ncF



	real*8 x1,y1,z1,x0,y0,z0,x2,y2,z2,x3,y3,z3,x4,y4,z4
	real*8 api,bpi,cpi,a1,b1,api2,bpi2,cpi2,api3,bpi3,cpi3
	real*8 api4,bpi4,cpi4,rt
	real*8 xp2,yp2,zp2,xp3,yp3,zp3,dd1,dd2,dd3,tstar2,tstar3
	real*8 alphap2,deltap2,rp2,alphap3,deltap3,rp3
	real*8 dp3pi1,dp3pi4,dp2pi1,dp2pi4
	real*8 xp,yp,zp,xt,yt,zt  

	REAL*8 XPA,YPA,ZPA,ALPHAPA,DELTAPA,RPA,tstarA,dpa
	REAL*8 XPB,YPB,ZPB,ALPHAPB,DELTAPB,RPB,tstarB,dpb

	character*40 filein1,filein2,filein3
	character*40 fileout1,fileout2,fileout3,fileout4



	
	jjk=1
	rkmax=0
	rlambda=1. 
	rmax=-10000
	rmin=10000
	pi=2.*asin(1.)

c	rborder_max=-1000
c	write(*,*)'imax='
c	read(*,*)imax

	write(*,*)'random samplig 0<p<1  ='
	read(*,*)pp0

	write(*,*)'radius  ='
	read(*,*)xfrac

	imax=50

	open(30,file='fileinput.info',status='unknown') 
	open(99,file='adr.dat',status='unknown') 


c******************** test test test 
	open(11,file='geometry.dat',status='unknown') 
	open(12,file='P1.dat',status='unknown') 
	open(13,file='P2.dat',status='unknown') 
	open(14,file='P3.dat',status='unknown')
 	open(15,file='P4.dat',status='unknown') 

c 	open(16,file='Point.dat',status='unknown') 
c 	open(17,file='retta_t1.dat',status='unknown')
c	open(18,file='retta_s1.dat',status='unknown') 
c	open(19,file='Pc2.dat',status='unknown')
c 	open(21,file='Pc3.dat',status='unknown') 
c	open(22,file='Points.dat',status='unknown') 
c	open(27,file='Points_cyl.dat',status='unknown') 
c******************** test test test 



	
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


c*********** test test test
c        ar_min=-36.0 
c        ar_max=36.0
c        dec_min=-48.
c        dec_max=48
c*********** test test test
	

	ar_min=ar_min*pi/180.
	ar_max=ar_max*pi/180.
	dec_min=dec_min*pi/180.
	dec_max=dec_max*pi/180.
	
	omega=(ar_max-ar_min)*(sin(dec_max)-sin(dec_min))
	write(*,*)'solid angle=',omega 


c************************* coordinates of the extremes 


	   x4=1*r_max*cos(ar_min)*cos(dec_min)
	   y4=1*r_max*sin(ar_min)*cos(dec_min)
	   z4=1*r_max*sin(dec_min)
	   
	   x3=1*r_max*cos(ar_min)*cos(dec_max)
 7	   y3=1*r_max*sin(ar_min)*cos(dec_max)
	   z3=1*r_max*sin(dec_max)

	   
	   x2=1*r_max*cos(ar_max)*cos(dec_min)
	   y2=1*r_max*sin(ar_max)*cos(dec_min)
	   z2=1*r_max*sin(dec_min)
	   
	   x1=1*r_max*cos(ar_max)*cos(dec_max)
	   y1=1*r_max*sin(ar_max)*cos(dec_max)
	   z1=1*r_max*sin(dec_max)
	   
	   
	   x0=0
	   y0=0
	   z0=0
	   
	   write(11,*)x0,y0,z0
	   write(11,*)x4,y4,z4
	   write(11,*)x2,y2,z2
	   write(11,*)x0,y0,z0
	   write(11,*)x3,y3,z3
	   write(11,*)x1,y1,z1
	   write(11,*)x0,y0,z0
	   write(11,*)x2,y2,z2
	   write(11,*)x1,y1,z1
	   write(11,*)x0,y0,z0
	   write(11,*)x3,y3,z3
	   write(11,*)x4,y4,z4
	   CLOSE(11) 

	 
c******************************* Piano P1 [P0,P1,P4]
	   api=(y2*z4-z2*y4)
	   bpi=-(x2*z4-z2*x4)
	   cpi=(x2*y4-y2*x4)

	   do i=1,1000
	      xp=200*rand(0)
	      yp=150*2*(rand(0)-0.5)
	      zp=-(xp*api+yp*bpi)/cpi
	      write(12,*)xp,yp,zp 
	   end do
c******************************* Piano P2 [P0,P3,P4]
	   api2=(y3*z4-z3*y4)
	   bpi2=-(x3*z4-z3*x4)
	   cpi2=(x3*y4-y3*x4)

	   do i=1,1000
	      xp=200*rand(0)
	      zp=150*2*(rand(0)-0.5)
	      yp=-(xp*api2+zp*cpi2)/bpi2
	      write(13,*)xp,yp,zp 
	   end do

c******************************* Piano P3 [P0,P1,P2]

	   api3=(y1*z2-z1*y2)
	   bpi3=-(x1*z2-z1*x2)
	   cpi3=(x1*y2-y1*x2)

	   do i=1,1000
	      xp=200*rand(0)
	      zp=150*2*(rand(0)-0.5)
	      yp=-(xp*api3+zp*cpi3)/bpi3 
	      write(14,*)xp,yp,zp 
      end do

c******************************* Piano P4 [P0,P1,P3]
	   api4=(y1*z3-z1*y3)
	   bpi4=-(x1*z3-z1*x3)
	   cpi4=(x1*y3-y1*x3)

	   do i=1,1000
	      xp=200*rand(0)
	      yp=150*2*(rand(0)-0.5)
	      zp=-(xp*api4+yp*bpi4)/cpi4
	      write(15,*)xp,yp,zp 
	   end do


c**************************************************************
c       READ FILE
c**************************************************************



	   write(20,*)'# r(1),Nave(2),Nvar(3),nave(4),nave(5),nvar(6)
	1,nc(7),NaveF(8),Navan(9),delta(10),vardelta(11)'


	   i=1 
 12	   continue 
 15	   read(10,*,end=11)x(i),y(i),z(i) 
	   if(rand(0).gt.pp0)goto 77

c*********************Transformation into spherical coordinates *******
	
	r(i)=sqrt(x(i)**2.+y(i)**2.+z(i)**2.)
	if(r(i).lt.r_min)goto 77
	if(r(i).gt.r_max)goto 77
	
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

	   i=i+1 
 77	   continue
	   goto 12 

 11	   close(10) 
	   close(99)
	   nd=i-1 
	   write(*,*)'nd=',nd

	   write(*,*)'ar_min,ar_max=',ar_min,ar_max
	   write(*,*)'dec_min,dec_max=',dec_min,dec_max



c**************************************************************
c       COMPUTING GAMMA 
c**************************************************************

c******** parameters of the r-range

	   rrmin=0.1
	   rrmax=200. 

	   qq=10.**(log(rrmax/rrmin)/log(10.)*1./imax)
	   write(*,*)'imax,qq =',imax,qq
	   do i=1,imax 
	      rr(i)=rrmin*qq**i 
	      h2(i)=xfrac 
c	      h2(i)=rr(i)/xfrac 

c	      write(*,*)'i,rr(i),rrmin,rrmax,qq,h2(i)=',
c     & i,rr(i),rrmin,rrmax,qq,h2(i) 
	      end do


C*******DO LOOP ON THE r VALUES TO COMPUTE GAMMA 
	      do i=1,imax 
		 write(*,*)'i,rr(i),h2(i)=',i,rr(i),h2(i) 

C*******1ST    DO LOOP ON # OF POINTS
	      
	      do j=1,nd
		 if(mod(j,5000).eq.0) write(*,*)j
		 M(i,j)=0.


c*****  BOUNDARY CONDITIONS CYLINDER

		 rb1P(j)=r(j)-r_min
		 if(r_max-r(j).lt.rb1P(j))rb1P(j)=r_max-r(j)		 
		 r1=abs(r(j)*cos(delta(j))*sin(alpha(j)-ar_min))
		 r2=abs(r(j)*cos(delta(j))*sin(ar_max-alpha(j)))
		 r3=abs(r(j)*sin(delta(j)-dec_min))
		 r4=abs(r(j)*sin(dec_max-delta(j)))
		 rb2P(j)=r1
		 if(r2.lt.rb2P(j))rb2P(j)=r2
		 if(r3.lt.rb2P(j))rb2P(j)=r3
		 if(r4.lt.rb2P(j))rb2P(j)=r4


c*************************************       retta t1 

		 b1=(api*z(j)/x(j)-cpi)/(bpi-api*y(j)/x(j))
		 a1=-(y(j)*b1+z(j))/x(j)

c******************************************* Punto PA

		 tstarA= sqrt(rr(i)**2./(a1**2.+b1**2.+1.))
		 xpA=x(j)+a1*tstarA
 		 ypA=y(j)+b1*tstarA
		 zpA=z(j)+tstarA

	
		 rPA =sqrt(xpA**2.+yPA**2.+zPA**2.)

		 		 
		 IF(XPA.EQ.0.) THEN
		    IF(YPA.GT.0)ALPHAPA=PI/2.
		    IF(YPA.LT.0)ALPHAPA=3*PI/2.
		    GOTO 611
		 END IF
		 
		 IF(YPA.EQ.0)THEN
		    IF(XPA.GT.0)ALPHAPA=0
		    IF(XPA.LT.0)ALPHAPA=PI
		    GOTO 611
		 END IF
		 
		 IF(XPA.GT.0.AND.YPA.GT.0)Q=0
		 IF(XPA.LT.0.AND.YPA.GT.0)Q=2
		 IF(XPA.LT.0.AND.YPA.LT.0)Q=2
		 IF(XPA.GT.0.AND.YPA.LT.0)Q=4
		 
		 alphaPA=atan(yPA/xPA) + PI/2. *Q 
 611		 CONTINUE
		 
		 IF(XPA.EQ.0.AND.YPA.EQ.0)THEN
		    IF(ZPA.GT.0)DELTAPA=PI/2.
		    IF(ZPA.LT.0)DELTAPA=-PI/2.
		    GOTO 612
		 END IF
		 
		 
		 deltaPA=asin(zPA/rPA)
 612		 CONTINUE 
		 if(alphaPA.gt.pi)alphaPA=alphaPA-2.*pi

		 IF(ALPHAPA.LT.AR_MIN)GOTO 13
		 IF(ALPHAPA.GT.AR_MAX)GOTO 13
		 IF(DELTAPA.LT.DEC_MIN)GOTO 13
		 IF(DELTAPA.GT.DEC_MAX)GOTO 13
		 if(rPA.lt.r_min)goto 13
		 if(rPA.gt.r_max)goto 13

c******************** Punto PB

		 tstarB=- sqrt(rr(i)**2./(a1**2.+b1**2.+1.))
		 xpB=x(j)+a1*tstarB
 		 ypB=y(j)+b1*tstarB
		 zpB=z(j)+tstarB

	
		 rPB =sqrt(xPB**2.+yPB**2.+zPB**2.)
		 if(rPB.lt.r_min)goto 13
		 if(rPB.gt.r_max)goto 13
		 		 
		 IF(XPB.EQ.0.) THEN
		    IF(YPB.GT.0)ALPHAPB=PI/2.
		    IF(YPB.LT.0)ALPHAPB=3*PI/2.
		    GOTO 711
		 END IF
		 
		 IF(YPB.EQ.0)THEN
		    IF(XPB.GT.0)ALPHAPB=0
		    IF(XPB.LT.0)ALPHAPB=PI
		    GOTO 711
		 END IF
		 
		 IF(XPB.GT.0.AND.YPB.GT.0)Q=0
		 IF(XPB.LT.0.AND.YPB.GT.0)Q=2
		 IF(XPB.LT.0.AND.YPB.LT.0)Q=2
		 IF(XPB.GT.0.AND.YPB.LT.0)Q=4
		 
		 alphaPB=atan(yPB/xPB) + PI/2. *Q 
 711		 CONTINUE
		 
		 IF(XPB.EQ.0.AND.YPB.EQ.0)THEN
		    IF(ZPB.GT.0)DELTAPB=PI/2.
		    IF(ZPB.LT.0)DELTAPB=-PI/2.
		    GOTO 712
		 END IF
		 
		 
		 deltaPB=asin(zPB/rPB)
 712		 CONTINUE 
		 if(alphaPB.gt.pi)alphaPB=alphaPB-2.*pi


		 IF(ALPHAPB.LT.AR_MIN)GOTO 13
		 IF(ALPHAPB.GT.AR_MAX)GOTO 13
		 IF(DELTAPB.LT.DEC_MIN)GOTO 13
		 IF(DELTAPB.GT.DEC_MAX)GOTO 13


c****************** bc for the CENTRAL POINT ************************
	      if(rb1P(j).lt.h2(i))goto 13
	      if(rb2P(j).lt.h2(i))goto 13
c****************** bc for the CENTRAL POINT ************************

		 M(i,j)=1.	!LABEL FOR THE Jth POINT AT THE SCALE Ith 
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

		    tstar=(a1*x(k)+b1*y(k)+z(k)-a1*x(j)-b1*y(j)-z(j))
     & /(a1*a1+b1*b1+1)
		    rp =abs(sqrt((x(J)+a1*tstar-x(k))**2.
     & +(y(j)+b1*tstar-y(k))**2. +(z(j)+tstar-z(k))**2.))

		    rt=abs(sqrt(rd**2-rp**2))
		    if(rp.le.h2(i))then
		       if(rt.le.rr(i))then
			  N(i,j)=N(i,j)+1. 

			  if(alpha(k).gt.alpha(j))NF(i,j)=NF(i,j)+1. 
			  if(alpha(k).lt.alpha(j))NN(i,j)=NN(i,j)+1. 

		       end if
		    end if
 88		    continue 
		 end do		!2nd (k)
 13		 continue 
	      end do		!1st (j) 
 99	      continue 
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

c	      if(ncF.gt.0.and.ncN.gt.0)then
		 write(20,'(11(e12.5,2x))')rr(i),Nave(i),Nvar(i),
	1	      Nave(i)/V(i),Nave(i)/V(i),Nvar(i)/V(i),float(nc)
	1	      ,NaveF(i) ,NaveN(i),(NaveN(i)-NaveF(i)) /(NaveN(i)
	2	      +NaveF(i)), sqrt(NvarN(i)**2+NvarF(i)**2)
	2	      /(NaveN(i)+NaveF(i))
c	      end if



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
	      
      if(k.eq.1) open(15,file='SL_CC_002.dat',status='unknown') 
      if(k.eq.2) open(15,file='SL_CC_005.dat',status='unknown') 
      if(k.eq.3) open(15,file='SL_CC_010.dat',status='unknown') 
      if(k.eq.4) open(15,file='SL_CC_020.dat',status='unknown') 
      if(k.eq.5) open(15,file='SL_CC_030.dat',status='unknown') 
      if(k.eq.6) open(15,file='SL_CC_040.dat',status='unknown') 
      if(k.eq.7) open(15,file='SL_CC_050.dat',status='unknown') 
      if(k.eq.8) open(15,file='SL_CC_060.dat',status='unknown') 
      if(k.eq.9) open(15,file='SL_CC_070.dat',status='unknown') 
      if(k.eq.10)open(15,file='SL_CC_080.dat',status='unknown') 
      if(k.eq.11)open(15,file='SL_CC_090.dat',status='unknown') 
      if(k.eq.12)open(15,file='SL_CC_100.dat',status='unknown') 
      if(k.eq.13)open(15,file='SL_CC_110.dat',status='unknown') 
      if(k.eq.14)open(15,file='SL_CC_120.dat',status='unknown') 
      if(k.eq.15)open(15,file='SL_CC_130.dat',status='unknown') 
      if(k.eq.16)open(15,file='SL_CC_140.dat',status='unknown') 
      if(k.eq.17)open(15,file='SL_CC_150.dat',status='unknown') 
      if(k.eq.18)open(15,file='SL_CC_160.dat',status='unknown') 
      if(k.eq.19)open(15,file='SL_CC_170.dat',status='unknown') 
      if(k.eq.20)open(15,file='SL_CC_180.dat',status='unknown') 
      if(k.eq.21)open(15,file='SL_CC_190.dat',status='unknown') 
      if(k.eq.22)open(15,file='SL_CC_200.dat',status='unknown') 
      if(k.eq.23)open(15,file='SL_CC_250.dat',status='unknown') 
      if(k.eq.24)open(15,file='SL_CC_300.dat',status='unknown') 	      
      if(k.eq.25)open(15,file='SL_CC_350.dat',status='unknown') 
      if(k.eq.26)open(15,file='SL_CC_400.dat',status='unknown') 
      if(k.eq.27)open(15,file='SL_CC_450.dat',status='unknown') 
      if(k.eq.28)open(15,file='SL_CC_500.dat',status='unknown') 
      if(k.eq.29)open(15,file='SL_CC_550.dat',status='unknown') 
      if(k.eq.30)open(15,file='SL_CC_600.dat',status='unknown') 
      if(k.eq.31)open(15,file='SL_CC_650.dat',status='unknown') 


      write(15,*)'#r(1),N(2),F(3),x(4),y(5),z(6),rr(7),al(8),del(9)'

	      
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
     %  x(j),y(j),z(j),rr(l),alpha(j),delta(j)
c     %   F(l,j),x(j),y(j),z(j),rr(l),alpha(j),delta(j)
		 end if
	      end do		!points
	      
	      
 68	      continue
	   end do		!scales where writes on files 
	   

c********************************************************
 1111	   continue 

	   end  
	
 
