c  Prrogram calculates time of the Planet cooling starting from temperature T_0
c
c
	implicit real*8 (a-h,o-z) 
c
	dimension tem(1000),rad(1000),vol(1000),ulay(1000),xmas(1000)
	dimension q(1000),esurf(1000),den(1000)
	dimension rr(600),rho(600)
	logical myf

c  Parameters
c
	pi=3.14159265359d0
	Gn=6.67408d-11  ! Gravitational constant, m^3/kg/s^2
c
	RadEar=6371.0*1.0e+3 ! m  Earth mean radius
c	Pmass=5.97237e+24 ! kg  Earth
c
	Radius=1737.1*1.0e+3 ! m Moon
	Pmass=7.342e+22      ! kg Moon
c
	delta=RadEar-Radius  ! difference in radii of the Planet and Earth
	board=1500.0*1.0d+3  ! m, boarder for r/active isotopes from surface
	crust=20.0*1.0d+3    ! m, crust depth from surface
	zero=-273.15         ! K, zero temperature in Celsium
	dh=100               ! m, standard expantion of layer
c
	Stef_B=5.670367e-8 ! Stefan Boltzmsnn constant W m^-2 K^-4
c	condi=1.5d0 ! thermal conductivity of Planet W m^-1 K^-1 (Planet ground)
	condi=1.0e+2 ! thermal conductivity of Planet W m^-1 K^-1 (iron)
	conds=1.5 ! thermal conductivity of Planet W m^-1 K^-1 (iron)
	capacy=0.5*1.0e+3 ! thermal capacity of Planet J kg^-1 K^-1 (most materials 0.1 - 0.9)
	tau=6.8e+3 ! My, roughly averaged live time for U, Th, K
c	Qth=5.0e-11 ! W/kg heat release in volume unit (Earth = 4.0e-11)
c
	tins=700.-zero ! temperature inside the Planet
	tout=0.         ! K, temperature outside the Planet (Cosmos)
	print *,'----------------------------------------'
	print '(a8,e12.4,a2)','Tstart= ',tins,' K'
	print *,'----------------------------------------'
c
	xln2=alog(2.)
	tk40=1.277d+3/xln2 ! tau in My
	tu38=4.468d+3/xln2 ! tau in My
	th32=1.405d+4/xln2 ! tau in My
c
	qk40=3.67693d-11   ! W/kg now if suppose F_tot = 244 TW
	qu38=2.12646d-12   ! W/kg now if suppose F_tot = 244 TW
	qt32=2.04274d-12   ! W/kg now if suppose F_tot = 244 TW
c
	qk40s=qk40*exp(4.5d+3/tk40)   ! W/kg at start if suppose F_tot = 244 TW
	qu38s=qu38*exp(4.5d+3/tu38)   ! W/kg at start if suppose F_tot = 244 TW
	qt32s=qt32*exp(4.5d+3/th32)   ! W/kg at start if suppose F_tot = 244 TW
c
	Qth=qk40s+qu38s+qt32s ! W/kg heat release in volume unit (Earth = 4.0e-11)
c	print '(e15.7,a5)', Qth,' W/kg'
c
	n=1000
	step_r=Radius/n ! step in radius, m
	dt=3.1556952d+7  ! step in time, 1 year (in sec)
	tcur=0.0
c
	inquire(file='planet.txt', Exist=myf)
c	print *,myf
c
	if(myf) then
	 open(1,file='planet.txt')
	 read(1,'(2e15.7)') Radius,Pmass
	 read(1,'(2e15.7)') condi,conds,capacy
	 read(1,'(4e15.7)') qk40s,qu38s,qt32s,tcur
	 read(1,'(i5)') n
	 do i=1,n
	  read(1,*) rad(i),tem(i)
c	  print '(e15.5,e15.7)',rad(i),tem(i)
	 enddo
	 close(1)
	 open(2,file='planet.txt')
	 write(2,'(2e15.7)') Radius,Pmass
	 write(2,'(2e15.7)') condi,conds,capacy
	else
	 open(2,file='planet.txt')
	 write(2,'(2e15.7)') Radius,Pmass
	 write(2,'(2e15.7)') condi,conds,capacy
	 do i=1,n
	  tem(i)=tins                 ! K
	 enddo
	endif
c
	print '(a3,e12.4,a2)','R= ',Radius,' m'
	print '(a3,e12.4,a3)','M= ',Pmass,' kg'
	Pvolum=4.*pi*Radius**3/3.
	pldens=Pmass/Pvolum
	print '(a6,e12.4,a8)','dens= ',pldens,' kg/m^3'
	print *,'----------------------------------------'
	print '(e12.4,a12)',Stef_B,' W m^-2 K^-4'
	print '(e12.4,a12)',condi,' W m^-1 K^-1'
	print '(e12.4,a12)',conds,' W m^-1 K^-1'
	print '(e12.4,a12)',capacy,' J g^-1 K^-1'
	print '(e15.7,a5)', Qth,' W/kg'
	print '(e15.7,a5)', tcur,' year'
c
	open(1,file='prem.txt') ! density profile real in g/cm^3
	read(1,*) nden
	do i=1,nden
	 read(1,*) rr(i),rho(i)
c	print '(e12.4,a3,e12.4,a7)',rr(i),' m ',rho(i),' g/cm^3'
	enddo
c	print *,nden
	close(1)
c
	do i=1,n
	 rad(i)=Radius-step_r*(i-1)     ! m
	 esurf(i)=4.*pi*rad(i)*rad(i)    ! m^2
	do j=1,nden
	 if(rad(i)+delta.le.rr(j)*1.0d3) then
	  den(i)=rho(j)*1.0d3 ! kg m^-3
	  goto 5
	 endif
	enddo
 5	continue
c	print '(e12.4,a3,e12.4,a7)',rad(i),' m ',den(i),' kg/m^3'
	enddo
c
c   Starting conditions
c
	do i=1,n
	 if(i.ge.n) then
	  vol(i)=4.*pi*rad(i)**3/3.   ! m^3
	 else
	  vol(i)=4.*pi*(rad(i)**3-rad(i+1)**3)/3. ! m^3
	 endif
	 xmas(i)=den(i)*vol(i)         ! kg
	 ulay(i)=capacy*xmas(i)*tem(i) ! J
	 if(rad(i).ge.Radius-board) then
	  q(i)=Qth*den(i)              ! W/m^3
	 else
	  q(i)=0.0                  ! W/m^3
	 endif
c	print '(i5,e12.4)',i,rad(i)
c	print '(e12.4,a4,e12.4,a7)',vol(i),' m^3',den(i),' kg/m^3'
c	print '(i5,a5,e12.4,a3)',i,' E_l= ',ulay(i),' J '
	enddo
c
	nsec=1000000  ! cycle length in years
	curt=tins
c
	print *,'----------------------------------------'
c
c	  en_emit=0.2*Stef_B*curt*curt*curt*curt*esurf(1) ! W
c
	do k=1,200  ! cycle on epoches ( 1 epoche = 10^6 years )
c
	 tcur=tcur+1.0
	 print '(i10,a3)',k,' My'
c
	qk40=qk40s*exp(-tcur/tk40)   ! W/kg at start if suppose F_tot = 244 TW
	qu38=qu38s*exp(-tcur/tu38)   ! W/kg at start if suppose F_tot = 244 TW
	qt32=qt32s*exp(-tcur/th32)   ! W/kg at start if suppose F_tot = 244 TW
c
	Qth=qk40+qu38+qt32 ! W/kg heat release in volume unit (Earth = 4.0e-11)
c	 print '(e15.7,a5)', Qth,' W/kg'
c
	 do j=1,nsec  ! cycle on years ( 1 year = 3.1556952d+7 sec )
c
	  do i=1,n ! cycle on n layers of Planet radius
c
c	print '(i5)',i
c	print '(a5,e12.4,3a)',' E_l= ',ulay(i),' J '
c	print '(e12.4,3a)',tem(i),' K '
c
	   if(i.le.1) then
	    emit=0.01*Stef_B*curt*curt*curt*curt*esurf(i) ! emitted Energy, W
	   if(emit.gt.ulay(i)) emit=ulay(i)/dt
	    erad=q(i)*vol(i)             ! Energy released by radioactivity, W
	    grT=0.                       ! temperature gradient, K/m 
	    flT=conds*grT                ! Energy flux density, W/m^2
	    transf=0.
	    en_emit=emit
	   else
	    emit=0.
	    erad=q(i)*vol(i)             ! Energy released by radioactivity, W
	    if(tem(i-1).eq.tem(i)) then
	    grT=0.                       ! temperature gradient, K/m 
	    else
	    grT=(tem(i-1)-tem(i))/step_r ! temperature gradient, K/m 
	    endif
	     if(rad(i).ge.Radius-crust) then
	      flT=conds*grT                ! Energy flux density, W/m^2
	     else
	      flT=condi*grT                ! Energy flux density, W/m^2
	     endif
	      transf=flT*esurf(i)          ! Energy transfered between layers (per sec), W
	   endif
c	print '(3e12.4)',erad,emit,transf
c	print '(4e12.4)',xmas(i),capacy,grT,condi
	if(i.gt.1) then
	   ulay(i)=ulay(i)+(erad-emit+transf)*dt ! change of i-th layer energy
	   ulay(i-1)=ulay(i-1)-transf*dt  ! change of (i-1)-th layer energy on transf energy from i-th
	   tem(i)=ulay(i)/xmas(i)/capacy  ! temperature of i-th layer, K
	   tem(i-1)=ulay(i-1)/xmas(i-1)/capacy  ! temperature of (i-1)-th layer, K
	else
	   ulay(i)=ulay(i)+(erad-emit+transf)*dt ! change of i-th layer energy
	   tem(i)=ulay(i)/xmas(i)/capacy  ! temperature of i-th layer, K
	endif
c	print '(a5,e12.4,a2)',' E_l= ',ulay(i),' J'
c	print '(e12.4,a2)',tem(i),' K'
c	tcrit=1.0d+3+2.0*i
c	if(tem(i).gt.tcrit) then
c	print '(2e12.4,a2)',tcrit,tem(i),' K'
c	 sma=0.
c	 do m=i+1,n
c	 sma=sma+xmas(m)
c	 enddo
c	 du=sma*Gn*xmas(i)/rad(i)/rad(i)*dh
c	 du=0.1*ulay(i)
c	 ulay(i)=ulay(i)-du             ! decrease of i-th layer energy
c	 tem(i)=ulay(i)/xmas(i)/capacy  ! temperature of i-th layer, K
c	print '(i5)',i
c	print '(a5,e12.4,a2)',' E_l= ',ulay(i),' J'
c	print '(e12.4)',du
c	print '(e12.4,a2)',tem(i),' K'
c	endif
c
c
	   q(i)=Qth*den(i)
c
	  enddo ! cycle on n layers
c
	  curt=tem(1)
	  if(tem(2).le.10.) goto 3
c
	 enddo  ! cycle on years
c
	print '(e15.7,a5)', Qth,' W/kg'
	print '(e15.7,a5)', tcur,' year'
	print *,'----------------------------------------'
 	print '(a7,e12.5,a2)',' Emitt ',en_emit*dt,' W'
	print '(e15.5,e15.7,a2)',rad(1),tem(1),' K'
	print '(e15.5,e15.7,a2)',rad(2),tem(2),' K'
	print '(e15.5,e15.7,a2)',rad(n/4),tem(n/4),' K'
	print '(e15.5,e15.7,a2)',rad(n/2),tem(n/2),' K'
	print '(e15.5,e15.7,a2)',rad(n-n/4),tem(n-n/4),' K'
	print '(e15.5,e15.7,a2)',rad(n),tem(n),' K'
c
	enddo ! cycle on epoches
c
 3	print *,'----------------------------------------'
c
	write(2,'(4e15.7)') qk40s,qu38s,qt32s,tcur
	write(2,'(i5)') n
	do i=1,n
c	print '(e12.4,a1,e15.7,3a)',rad(i),' ',tem(i),' K '
	write(2,'(e12.4,a1,e15.7)') rad(i),' ',tem(i)
	enddo
	close(2)
c
	stop
	end

