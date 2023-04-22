!--------------------------------
!
!
!!   ע���޸�di
!
!
!-------------------------------

! usage : gfortran Read_DAT-new.f90 -o Read_DAT_new

program cygread  
parameter (mx=101,my=101,mz=101)
integer i,j,k,m,n1,n2,n3
dimension x(mx*my*mz,12),xx(mx),yy(my),zz(mz),fac(mx*my*mz)
real xmax,ymax,zmax,xmin,ymin,zmin,dx,dy,dz
real gama,di
character*7 input
character*16 output
character*3 cn
write(*,*)'Please input the n1(start) and n2 (end) and step'
	read(*,*)n1
	read(*,*)n2
	read(*,*)n3

do  n=n1,n2,n3
write(*,*)'Data is being read...',n
    input='m3ds'//cn(n)
	output(1:12)='cygp12var'//cn(n)
	output(13:16)='.dat' 
	gama=1.66666666667    
  	xmax=15
	ymax=50
	zmax=50
!      xmin=-xmax
	xmin=-xmax
	ymin=-ymax
	zmin=-zmax
	dx=(xmax-xmin)/(mx-1)
	dy=(ymax-ymin)/(my-1)
	dz=(zmax-zmin)/(mz-1)
	di=0.0

!---------------
    do i=1,mx
	xx(i)=xmin+(i-1)*dx
	end do
!--------------
	do  j=1,my
	yy(j)=ymin+(j-1)*dy
	end do
!--------------
	do  k=1,mz
	zz(k)=zmin+(k-1)*dz
	enddo
!--------------

    open(unit=12,file=input,status="unknown",form="formatted")
    read(12,10)((x(i,j),j=1,11) ,i=1,mx*my*mz)
	close(12)
10	format(11(1x,e10.4))

    write(*,*)'****1'
    do i=1,mx*my*mz

	x(i,2)=x(i,2)/x(i,1)
	x(i,3)=x(i,3)/x(i,1)
	x(i,4)=x(i,4)/x(i,1)    

!	x(i,12)=x(i,2)-x(i,9)*di/x(i,1)
!	x(i,13)=x(i,3)-x(i,10)*di/x(i,1)
!	x(i,14)=x(i,4)-x(i,11)*di/x(i,1)
    x(i,8)=(gama-1)*((x(i,8)-x(i,1)*(x(i,2)*x(i,2)+x(i,3)*x(i,3)+x(i,4)*x(i,4))/2-(x(i,5)*x(i,5)+x(i,6)*x(i,6)+x(i,7)*x(i,7))/2))!ѹǿ

	x(i,12)=x(i,5)*x(i,5)+x(i,6)*x(i,6)+x(i,7)*x(i,7)
	enddo

     write(*,*)'****2'
   open(unit=13,file=output,status='unknown',form='formatted')
	write(13,*)'TITLE=magnetic_field'
	write(13,*)'VARIABLES="x" "y" "z" "rho" "Vx" "Vy" "Vz" "bx" "by" "bz" "Jx" "Jy" "Jz" "Bt2" "pr" '!    "fac" '   !���з��ò���
!	write(13,*)'VARIABLES="x" "y" "z" "rho" "Vx" "Vy" "Vz" "bx" "by" "bz" "Jx" "Jy" "Jz" "Vex" "Vey" "Vez" "pr" '!    "fac" '   !���з��ò���
	write(13,200)mx,my,mz

   write(13,20) (((xx(i),yy(j),zz(k),(x(my*mx*(k-1)+mx*(j-1)+i,m),m=1,7),(x(my*mx*(k-1)+mx*(j-1)+i,m),m=9,11),&
   (x(my*mx*(k-1)+mx*(j-1)+i,m),m=12,12),(x(my*mx*(k-1)+mx*(j-1)+i,m),m=8,8),i=1,mx),j=1,my),k=1,mz)

20	format(15(1x,e10.4))
200	format('ZONE I=',i3,' J=',i3,' K=',i3,' F=POINT ')

    write(*,*) 'end'

enddo
end


character*3 function cn(n)
!-----assume that n is no greater than 999
!-----separate the digits
      n1=n/100
      n2=(n-100*n1)/10
      n3=(n-100*n1-10*n2)
!-----stick together cn using char function
      n1=n1+48
      n2=n2+48
      n3=n3+48
      cn(1:1)=char(n1)
      cn(2:2)=char(n2)
      cn(3:3)=char(n3)
      return
      end
