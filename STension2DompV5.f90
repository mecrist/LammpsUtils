PROGRAM stension_final
implicit none

character(50)            ::inputfile
integer                  ::nconf,j,l,k,d,o,i,tid
integer                  ::natoms
integer                  ::nbins,m
real*8,parameter         ::pi=3.14159
integer                  ::ntype,b,dummyint,nnumber,jump
character(100)           ::dummychar,id
real*8                   ::t1,delta,vx,vy,vz
real*8                   ::hix,lox,temp,sigmabig,temp2,bla
real*8                   ::hiy,loy,hiz,loz,std,ste
real*8                   ::x,y,z,rxij,ryij,rzij
real*8,dimension(:),allocatable :: Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,temp1
real*8                   ::dp,rmax,dr,low,high,r,dV,norm,theta
real*8                   ::xh,xl,yh,yl,zh,zl
real*8                   ::dummy
real*8,dimension(:),allocatable          ::xx,yy,zz,Ptot,pwater,dpdr
real*8,dimension(:),allocatable          ::Ptot1,Ptot2,Ptot3,Ptot4,psili,sigma
real*8,dimension(:,:),allocatable        ::soma
integer,dimension(:),allocatable       ::nntype
integer,dimension(:),allocatable  ::idd
!ccccccccccccccccccccccccccccccccccccccccc

!Setting de walltime timer
t1=secnds(0.0)

!Read the input files
nconf=800  !Number of steps saved in the dump, take into account the writing frequency
nbins=210  !Number of bins to divide the cylinder, to calculate separately the pressure on each bin then derive
jump=1

!Open the necessary files
open(9,file='positions.lammpstrj')

!Read initial data from the dump file, to get natoms

nconf=nconf-jump

do i=1,jump
  print*,'Jumping ',i
  read(9,*)dummychar
  read(9,*)dummyint
  read(9,*)dummychar
  read(9,*)natoms
  read(9,*)dummychar
  read(9,*)lox,hix
  read(9,*)loy,hiy
  read(9,*)loz,hiz
  read(9,*)dummychar
  do j=1,natoms
    read(9,*)nnumber,ntype,dummyint,x,y,z,Pxx,Pyy,Pzz !IMPORTANT!! That is the format that the dump file should be: Atom Number, Atom Type, x, y, z, Pxx, Pyy, Pzz, Pxy, Pxz, Pyz
  enddo
enddo

!Calculataing the thickness of each bin
dr=(hiz-loz)/nbins
print*,'bin size = ',dr

!Allocating the variables now that the natoms is known
allocate(xx(1:natoms),yy(1:natoms),zz(1:natoms),Ptot(1:natoms),soma(1:natoms,1:nconf),temp1(1:nbins))
allocate(Pxx(1:natoms),Pyy(1:natoms),Pzz(1:natoms),Pxy(1:natoms),Pxz(1:natoms),Pyz(1:natoms))
allocate(pwater(1:natoms),Ptot1(1:natoms),Ptot2(1:natoms),sigma(1:nconf),dpdr(1:nconf))
allocate(Ptot3(1:natoms),Ptot4(1:natoms),psili(1:natoms),nntype(1:natoms),idd(1:natoms))

!Reseting the counters
sigmabig=0.0
o=0
m=0

do b=1,nconf !loop on step number
   m=m+1 ! Incrementing counter m, this counter is used if you want to skip some steps (oversized dump file)
         !to change the number of skipped steps change the ifs in lines 82,99,120,147 from 1 to the desired value.
         ! You can also set the m at line 65 to a negative value, if you want to skip some initial equilibration steps.
!Reading dump head
   read(9,*)dummychar
   read(9,*)dummyint
   read(9,*)dummychar
   read(9,*)dummyint
   read(9,*)dummychar
   read(9,*)lox,hix
   read(9,*)loy,hiy
   read(9,*)loz,hiz
   read(9,*)dummychar

   do j=1,natoms !Loop over atom number
     read(9,*)nnumber,ntype,dummyint,xx(nnumber),yy(nnumber),zz(nnumber),Pxx(nnumber),Pyy(nnumber),Pzz(nnumber) !Reading individual atom data
     zz(nnumber)=zz(nnumber)-loz
     idd(j)=nnumber !Setting the data from the temp variables
     nntype(nnumber)=ntype !Setting the data from the temp variables
     Ptot(nnumber)=(Pzz(nnumber)-((Pxx(nnumber)+Pyy(nnumber))/2))!/(3*(hix-lox)*(hiy-loy)*(hiz-loz))
   enddo

   do l=1,nbins
     soma(l,b)=0.0 !Reseting the pressure counter 
     low=(l-1)*dr !Defining the bin start
     high=l*dr    !Defining the bin end
     do k=1,natoms
       d=zz(k)
       if ((d.gt.low).and.(d.le.high)) then !testing if it is inside the bin
         soma(l,b)=soma(l,b)+Ptot(k) !-(Pzz(k)/dr)+(((Pxx(k)+Pyy(k))/2)/(2*(hix-lox)*(hiy-loy))) !Ptot(k) !Accumulate on pressure counter
       endif
     enddo
     soma(l,b)=-soma(l,b)/(3*(hix-lox)*(hiy-loy)*dr)
   enddo

   !Reseting Counters 
   temp=0.0
   temp1=0.0
   temp2=0.0
   sigma(b)=0.0
   dpdr(b)=0

 if(m==1)then
   do i = 2,nbins !Loop over selected bins, those who has derivatives
     !dV=(hiz-loz)*(hiy-loy)*dr !Volume of the bin
     !soma(i,b)=soma(i,b)/(3*(hix-lox)*(hiy-loy)*dr)
     dpdr(b)=dpdr(b)+(dr*(soma(i,b)+soma(i-1,b))/2) 
   enddo
 endif

 if(m==1)then
   sigma(b)=(dpdr(b))*1.01325E-5*0.5
   sigmabig=sigmabig+(sigma(b)) !Summing up all the integral values
   o=o+1 !incrementing the counter, to get the average after
   m=0
 endif
 
 bla=(sigmabig/o) !Calculate the average interfacial tension up to this step
 print *,b,sigma(b),bla !print the step number, the interfacial tension of this step, the average interfacial tension up to this step
enddo

sigmabig=(sigmabig/o) !calculate the final average for the interfacial tension
std=0 !reset counter

do i=1,nconf !Loop over all steps
  std=std+(sigma(i)-sigmabig)**2 !Calculating the standard deviation
enddo 

std=sqrt(std/DBLE(nconf)) !Calculating the standard deviation
ste=std/(sqrt(DBLE(nconf))) !Calculating the standard deviation

delta=secnds(t1) !Getting the walltime counter

temp=0
temp1=0

! do i=1,nconf
!  temp=temp+dpdr(i)
! enddo

open(47,file='pressures-james.dat')

do i=1,nbins
 do j=1,nconf
  temp1(i)=temp1(i)+soma(i,j)
 enddo
 temp1(i)=temp1(i)/nconf
 write(47,*)(i-1)*dr,temp1(i)
enddo
 
temp=temp/nconf

!Write the final values
write(*,*)'Execution Time = ',delta
write(*,'(3A15)')'Surf. Tension','Std. Dev.','Std. Err.'
write(*,'(3F15.5)')sigmabig,std,ste
! print*,temp

! close(11)
! close(10)
close(9)


 END PROGRAM stension_final
