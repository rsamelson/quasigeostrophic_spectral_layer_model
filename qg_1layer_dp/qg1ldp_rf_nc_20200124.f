c
      program qg1l_dp
c
c -------------------------------
c  Public distribution copy of quasigeostrophic model code used for the calculations described by:
c
c       Samelson, R. M., D. B. Chelton, and M. G. Schlax, 2019.  The ocean mesoscale regime
c         of the reduced-gravity quasi-geostrophic model.  J. Phys. Oceanogr., 49, 2469–2498,
c         DOI: 10.1175/JPO-D-18-0260.1.
c
c   This version adds code for netcdf output files.   RMS  24 Jan 2020
c   Example compile command:
c    gfortran -I/usr/local/include -L/usr/local/lib -lnetcdf -lnetcdff qg1ldp_rf_nc.f ../../fft_src/fft99f.f -o qg1ldp_rf_nc.x
c
c  Copyright 2019, 2020 Roger M. Samelson
c
c Permission is hereby granted, free of charge, to any person obtaining a copy of this software
c  and associated documentation files (the "Software"), to deal in the Software without restriction,
c  including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
c  and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
c  subject to the following conditions:
c
c The above copyright notice and this permission notice shall be included in all copies or substantial
c  portions of the Software.
c
c THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
c  LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
c  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
c  WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
c  SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
c ----------------------
c
c     This code uses the fft99f.f two-dimensional Fast Fourier Transform routine.  The version
c      that has been used for the published calculations with this code was obtained from
c      dsl@ncar.ucar.edu on Monday 9 September 1991 and contains the following history comments:
C      PACKAGE FFT99F HISTORY:
C              THE PACKAGE WAS WRITTEN BY CLIVE TEMPERTON AT ECMWF IN
C              NOVEMBER, 1978.  IT WAS MODIFIED, DOCUMENTED, AND TESTED
C              FOR NCAR BY RUSS REW IN SEPTEMBER, 1980.
c
c -------------------------------
c
c   Code history
c
c    Recent modifications for calculations described in Samelson et al. (2019)
c
c     RMS Dec 2017:  convert to doubly periodic from qg1lc_rf.f
c
c     RMS Oct 2017:  add random forcing
c                    convert to 1-layer model from qg2lc_rf.f
c
c    Original source and previous modifications
c
c     The two-layer quasi-geostrophic model code that is the original source for this code
c      was written by R. M. Samelson to perform the calculations described by:
c       Samelson, R. M., and J. Pedlosky, 1990. Local baroclinic instability of flow over
c         variable topography. Journal of Fluid Mechanics, 221, 411-436.
c       Oh, S. P., J. Pedlosky, and R. M. Samelson, 1993. Linear and finite-amplitude localized
c         baroclinic instability. Journal of the Atmospheric Sciences, 50(16), 2772-2784.
c
c     This version has been modified by Karl Helfrich from
c     the original code written by Roger Samelson.
c     KRH 10/92
c
c	  The Roberts scheme for controlling the saw-tooth
c	  computational instability has been added.
c	  KRH 11-12-92
c
c     This version has friction in the from of
c     a hyperviscosity = -r del^6 psi_n
c     KRH 12-23-93
c
c  N - x, M - y, XK, XL wavenumber constants (= 2 pi/Lx or pi/Ly)
c
      parameter (N=256,M=128)
      real a0(N+2,M+2),a1(N+2,M+2),ap(N+2,M+2)
      real p1(N+2,M+2)
      real fcap0(N/2+1,M/2+1),phi_a(N/2+1,M/2+1,2)
      real force_a(N+2,M+2)
c netcdf file and variable ids
      integer ncid,psi_varid,Fcap0_varid,nrec
      common /dumdif/ ddum((2*N+2)*(2*M+2)+(N+2)*(M+2))
      common /dumjac/ dumj(4*(2*N+2)*(2*M+2))
      common /workft/ work((2*N+1)*2*M),trigs(3*N+1),ifax(13)
      common /workft2/ work2((2*M+1)*2*(N+2)),trigs2(3*M+1),ifax2(13)
c r: del^6; rek: del^2, rpsi: del^0 
      common /frict/ ifrict,r,rek,rpsi
      common /parms/ xk,xl,dt,bb,f1,rob,twopi,dtn,c0_rf
c
c     if a continuation run set t=time of restart and icont=1
c
      data nt/5/,np/1/,ki/0/,io/1/,t/0./,icont/0/
c      data nt/50/,np/1/,ki/0/,io/1/,t/0./,icont/0/
c      data nt/400/,np/1/,ki/0/,io/1/,t/0./,icont/0/
c      data nt/2500/,np/1/,ki/0/,io/1/,t/0./,icont/0/
c
c
c --
      dt=0.002
c      j2m=int(1.0/dt+0.05)
      j2m=int(0.5/dt+0.05)
c      j2m=1
c
c  beta parameter
      bb=0.616
c  F1=1 (L=L_R)
      f1=1
      xk=1./6.
      xl=2*1./6
c  Robert filter parameter for leapfrog time-stepping
c      rob=0.001
c  Friction parameters
      ifrict=2
      rpsi=0.0216
      rek=0.0
      r0_del6=0.00865
      r=400000*r0_del6/(xk**2*(N/2)**2+xl**2*(M/2)**2)**3
c      r=200000*rpsi/(xk**2*(N/2)**2+xl**2*(M/2)**2)**3
c  Forcing autocorrelation timescale, wavenumbers, and time-step
      tau=0.925
      f_amp=1./sqrt(tau)
      dxkf=1./3.
      xkf=1.0
      c0_rf=0.
      ndt=10
      dtn=dt/ndt
      rcap=(tau-0.5*dtn)/(tau+0.5*dtn)
c
c  Parameter rescaling following Appendix a. Scaling of Samelson et al. (2019).
c  No longer required because normalization of forcing (force_a,fcap0)
c   has been corrected in this version of the code.
c  Default parameter values in code were computed by this rescaling
c   with F_stddev=0.89.
c
c     Acap_F=F_stddev**(2./3.)
c     bb=0.57/Acap_F
c     tau=1.00*Acap_F
c     rpsi=0.020/Acap_F
c     r0_del6=0.008/Acap_F
c
c --
c
      call init(ki,N+2,M+2,a0,a1,ap,icont)
      call rforc_init(N+2,M+2,N/2+1,M/2+1,xkf,dxkf,fcap0,force_a)
c
c
c      write(20,'(1x,2i5,8f10.5)') N,M,DT,XK,XL
c      write(20,'(1x,8e15.5)') BB,F1,rpsi,rek,r
c      write(20,'(1x,4i5,f7.2)') nt,ki,io,np,t
c
      write(22,'(1x,2i5,8f10.5)') N,M,DT,XK,XL
      write(22,'(1x,8e15.5)') BB,F1,rpsi,rek,r
      write(22,'(1x,8e15.5)') f_amp,tau,dxkf,xkf,c0_rf
      write(22,'(1x,4i5,f7.2)') nt,ki,io,np,t
c
c  initialize netcdf file and write initial data
      call out_nc_init(ncid
     .,  bb,f1,xk,xl,dt,ifrict,rpsi,rek,r0_del6,r,tau,f_amp
     .,  dxkf,xkf,c0_rf,ndt,dtn,rcap
     .,  N+2,M+2,psi_varid,Fcap0_varid,a0,force_a,p1)
c   netcdf file output record counter
      nrec=1
c
c      if(io.eq.1) call out(N+2,M+2,a0,p1,20)
c      if(io.eq.1) call out(N+2,M+2,force_a,p1,21)
c ***
      call intclc(N+2,M+2,a0,force_a,ep,ek,eq,ef)
      write(22,'(1x,e15.5,8(1x,e15.5))') t,ep,ek,eq,ef
      write( 6,'(1x,e15.5,8(1x,e15.5))') t,ep,ek,eq,ef
c ***
c  Start integration loop
c
c
      do 500 j1=1,nt
       do 400 j2=1,j2m
c  Runge-Kutta step
         call rkstep((N+2)*(M+2),t,a0,a1,ap)
c  Leapfrog time-step
c         call lfstep((N+2)*(M+2),t,a0,a1,ap)
c  Random forcing step
         call rforc_step(N+2,M+2,N/2+1,M/2+1,ndt,rcap,fcap0,phi_a
     .                   ,force_a,f_amp,a0)
400    continue
c ***
       call intclc(N+2,M+2,a0,force_a,ep,ek,eq,ef)
       write(22,'(1x,e15.5,8(1x,e15.5))') t,ep,ek,eq,ef
       write( 6,'(1x,e15.5,8(1x,e15.5))') t,ep,ek,eq,ef
c ***
c  write netcdf file record
      if(io.eq.1.and.((j1/np)*np).eq.j1) nrec=nrec+1
      if(io.eq.1.and.((j1/np)*np).eq.j1) call out_nc(ncid,nrec
     .,  N+2,M+2,psi_varid,Fcap0_varid,a0,force_a,p1)
c
c      if(io.eq.1.and.((j1/np)*np).eq.j1)call out(N+2,M+2,a0,p1,20)
c      if(io.eq.1.and.((j1/np)*np).eq.j1)call out(N+2,M+2,force_a,p1,21)
500   continue
c
      close(20)
      close(21)
      close(22)
      stop
      end

      subroutine rforc_init(nn,mm,n2,m2,xkf,dxkf,fcap0,force_a)
c      call with n2=nn/2,m2=mm/2
      real fcap0(n2,m2)
      real force_a(nn,mm)
      common /parms/ xk,xl,dt,bb,f1,rob,twopi,dtn,c0_rf
c
      fcap0sum=0.
      do l=0,m2-1
       k=0
       xkcap2=XK**2*k**2+XL**2*l**2
       fcap0(k+1,l+1)=
     .         max(0.,((xkf+dxkf)**2-xkcap2)*(xkcap2-(xkf-dxkf)**2)) 
       fcap0sum=fcap0sum+fcap0(k+1,l+1)**2
       k=n2-1
       xkcap2=XK**2*k**2+XL**2*l**2
       fcap0(k+1,l+1)=
     .         max(0.,((xkf+dxkf)**2-xkcap2)*(xkcap2-(xkf-dxkf)**2)) 
       fcap0sum=fcap0sum+fcap0(k+1,l+1)**2
      enddo
      do k=0,n2-1
       l=0
       xkcap2=XK**2*k**2+XL**2*l**2
       fcap0(k+1,l+1)=
     .         max(0.,((xkf+dxkf)**2-xkcap2)*(xkcap2-(xkf-dxkf)**2)) 
       fcap0sum=fcap0sum+fcap0(k+1,l+1)**2
       l=m2-1
       xkcap2=XK**2*k**2+XL**2*l**2
       fcap0(k+1,l+1)=
     .         max(0.,((xkf+dxkf)**2-xkcap2)*(xkcap2-(xkf-dxkf)**2)) 
       fcap0sum=fcap0sum+fcap0(k+1,l+1)**2
      enddo
      fcap0sum=2*fcap0sum
      do l=1,m2-2
        do k=1,n2-2
          xkcap2=XK**2*k**2+XL**2*l**2
          fcap0(k+1,l+1)=
     .         max(0.,((xkf+dxkf)**2-xkcap2)*(xkcap2-(xkf-dxkf)**2)) 
          fcap0sum=fcap0sum+4.*fcap0(k+1,l+1)**2
        enddo
      enddo
      fcap0norm=1./sqrt(fcap0sum)
      write(6,*) fcap0norm
      do l=1,m2
        do k=1,n2
          fcap0(k,l)=fcap0norm*fcap0(k,l)
        enddo
      enddo
      do l=1,mm
        do k=1,nn
         force_a(k,l)=0.
        enddo
      enddo
      return
      end


      subroutine rforc_step(nn,mm,n2,m2,ndt,rcap,fcap0,phi_a
     .                     ,force_a,f_amp,a0)
c      call with n2=nn/2,m2=mm/2
      real a0(nn,mm)
      real fcap0(n2,m2),phi_a(n2,m2,2)
      real force_a(nn,mm)
      common /parms/ xk,xl,dt,bb,f1,rob,twopi,dtn,c0_rf
c
      do jdt=1,ndt
       force_a(1,1)=0.
       force_a(2,1)=0.
       force_a(1,2)=0.
       force_a(2,2)=0.
       call random_number(phi_a)
       do k=1,n2-1
         l=0
         cp1=cos(twopi*phi_a(k+1,l+1,1))
         sp1=sin(twopi*phi_a(k+1,l+1,1))
         cp2=1.
         force_a(2*k+1,2*l+1)=rcap*force_a(2*k+1,2*l+1)
     .+fcap0(k+1,l+1)*sqrt(1-rcap**2)*cp1*cp2
         force_a(2*k+2,2*l+1)=rcap*force_a(2*k+2,2*l+1)
     .+fcap0(k+1,l+1)*sqrt(1-rcap**2)*sp1*cp2
         l=m2-1
         cp1=cos(twopi*phi_a(k+1,l+1,1))
         sp1=sin(twopi*phi_a(k+1,l+1,1))
         cp2=1.
         force_a(2*k+1,2*l+1)=rcap*force_a(2*k+1,2*l+1)
     .+fcap0(k+1,l+1)*sqrt(1-rcap**2)*cp1*cp2
         force_a(2*k+2,2*l+1)=rcap*force_a(2*k+2,2*l+1)
     .+fcap0(k+1,l+1)*sqrt(1-rcap**2)*sp1*cp2
       enddo
       do l=1,m2-1
         k=0
         cp1=1.
         cp2=cos(twopi*phi_a(k+1,l+1,2))
         sp2=sin(twopi*phi_a(k+1,l+1,2))
         force_a(2*k+1,2*l+1)=rcap*force_a(2*k+1,2*l+1)
     .+fcap0(k+1,l+1)*sqrt(1-rcap**2)*cp1*cp2
         force_a(2*k+1,2*l+2)=rcap*force_a(2*k+1,2*l+2)
     .+fcap0(k+1,l+1)*sqrt(1-rcap**2)*cp1*sp2
         k=n2-1
         cp1=1.
         cp2=cos(twopi*phi_a(k+1,l+1,2))
         sp2=sin(twopi*phi_a(k+1,l+1,2))
         force_a(2*k+1,2*l+1)=rcap*force_a(2*k+1,2*l+1)
     .+fcap0(k+1,l+1)*sqrt(1-rcap**2)*cp1*cp2
         force_a(2*k+1,2*l+2)=rcap*force_a(2*k+1,2*l+2)
     .+fcap0(k+1,l+1)*sqrt(1-rcap**2)*cp1*sp2
       enddo
       do l=1,m2-2
         do k=1,n2-2
            cp1=cos(twopi*phi_a(k+1,l+1,1))
            sp1=sin(twopi*phi_a(k+1,l+1,1))
            cp2=cos(twopi*phi_a(k+1,l+1,2))
            sp2=sin(twopi*phi_a(k+1,l+1,2))
            force_a(2*k+1,2*l+1)=rcap*force_a(2*k+1,2*l+1)
     .+fcap0(k+1,l+1)*sqrt(1-rcap**2)*cp1*cp2
            force_a(2*k+2,2*l+1)=rcap*force_a(2*k+2,2*l+1)
     .+fcap0(k+1,l+1)*sqrt(1-rcap**2)*sp1*cp2
            force_a(2*k+1,2*l+2)=rcap*force_a(2*k+1,2*l+2)
     .+fcap0(k+1,l+1)*sqrt(1-rcap**2)*cp1*sp2
            force_a(2*k+2,2*l+2)=rcap*force_a(2*k+2,2*l+2)
     .+fcap0(k+1,l+1)*sqrt(1-rcap**2)*sp1*sp2
         enddo
       enddo
      enddo

c c  zonal propagation
c        do l=1,m2
c          do k=1,n2
c             cp1=cos(xk*k*c0_rf*dtn)
c             sp1=sin(xk*k*c0_rf*dtn)
c             fa_temp=
c      .    cp1*force_a(2*k+1,2*l+1)+sp1*force_a(2*k+2,2*l+1)
c             force_a(2*k+2,2*l+1)=
c      .   -sp1*force_a(2*k+1,2*l+1)+cp1*force_a(2*k+2,2*l+1)
c             force_a(2*k+1,2*l+1)=fa_temp
c             fa_temp=
c      .    cp1*force_a(2*k+1,2*l+2)+sp1*force_a(2*k+2,2*l+2)
c             force_a(2*k+2,2*l+2)=
c      .   -sp1*force_a(2*k+1,2*l+2)+cp1*force_a(2*k+2,2*l+2)
c             force_a(2*k+1,2*l+2)=fa_temp
c          enddo
c        enddo
c       enddo
c

      do l=0,m2-1
        do k=0,n2-1
          a0(2*k+1,2*l+1)=a0(2*k+1,2*l+1)
     .     -dt*f_amp*force_a(2*k+1,2*l+1)/(XK**2*k**2+XL**2*l**2+F1)
          a0(2*k+2,2*l+1)=a0(2*k+2,2*l+1)
     .     -dt*f_amp*force_a(2*k+2,2*l+1)/(XK**2*k**2+XL**2*l**2+F1)
          a0(2*k+1,2*l+2)=a0(2*k+1,2*l+2)
     .     -dt*f_amp*force_a(2*k+1,2*l+2)/(XK**2*k**2+XL**2*l**2+F1)
          a0(2*k+2,2*l+2)=a0(2*k+2,2*l+2)
     .     -dt*f_amp*force_a(2*k+2,2*l+2)/(XK**2*k**2+XL**2*l**2+F1)
       enddo
      enddo
c
      return
      end
 
 
      subroutine lfstep(nm,t,a0,a1,ap)
c  One leapfrog time-step
      parameter (N=256,M=128)
c
      common /parms/ xk,xl,dt,bb,f1,rob,twopi,dtn,c0_rf
c
      real a0(nm),a1(nm),ap(nm)
      call diff(N+2,M+2,t,a1,a0,ap)
      do 40 i=1,(N+2)*(M+2)
       savea=a1(i)
       a1(i)=a0(i)+ap(i)*2.*DT
       a0(i)=savea+rob*(a0(i)+a1(i)-2.*savea)
40    continue
      t=t+DT
      return
      end


      subroutine diff(nn,mm,t,a,a0,ap)
c  Calculate time derivatives of spectral amplitudes
      parameter (N=256,M=128)
      common /dumdif/ c(2*N+2,2*M+2),d(N+2,M+2)
      common /frict/ ifrict,r,rek,rpsi
c
      common /parms/ xk,xl,dt,bb,f1,rob,twopi,dtn,c0_rf
c
      real a(nn,mm),a0(nn,mm)
      real ap(nn,mm)
c
c
c  Single-layer Jacobian
       call jacobc(nn,mm,a,F1,c)
c
c  Friction
c
      if(ifrict.ge.1)then 
c  ( dissip = - r*del^6(phi) ) hyperviscosity
        do 500 l=0,M/2
          do 490 k=0,N/2
             d(2*k+1,2*l+1)=r*(XK**2*k**2+XL**2*l**2)**3*a0(2*k+1,2*l+1)
             d(2*k+2,2*l+1)=r*(XK**2*k**2+XL**2*l**2)**3*a0(2*k+2,2*l+1)
             d(2*k+1,2*l+2)=r*(XK**2*k**2+XL**2*l**2)**3*a0(2*k+1,2*l+2)
             d(2*k+2,2*l+2)=r*(XK**2*k**2+XL**2*l**2)**3*a0(2*k+2,2*l+2)
  490     continue
  500   continue             
      if(ifrict.eq.2)then
c  add Ekman friction (equal in both layers) and linear psi damping
        do l=0,M/2
          do k=0,N/2
             d(2*k+1,2*l+1)=d(2*k+1,2*l+1)
     .        +(rek*(XK**2*k**2+XL**2*l**2)+rpsi)*a0(2*k+1,2*l+1)
             d(2*k+2,2*l+1)=d(2*k+2,2*l+1)
     .        +(rek*(XK**2*k**2+XL**2*l**2)+rpsi)*a0(2*k+2,2*l+1)
             d(2*k+1,2*l+2)=d(2*k+1,2*l+2)
     .        +(rek*(XK**2*k**2+XL**2*l**2)+rpsi)*a0(2*k+1,2*l+2)
             d(2*k+2,2*l+2)=d(2*k+2,2*l+2)
     .        +(rek*(XK**2*k**2+XL**2*l**2)+rpsi)*a0(2*k+2,2*l+2)
          enddo
	enddo 
       endif
      else
c  no friction
        do 540 l=1,M+2
        do 540 k=1,N+2
          d(k,l)=0.0
540     continue
      endif
c
c  Linear parts explicit
c
      do 200 l=0,M/2
       do 100 k=0,N/2
        ap(2*k+1,2*l+1)=
     .(1./(XK**2*k**2+XL**2*l**2+F1))
     .*( -BB*XK*k*a(2*k+2,2*l+1)
     .   +c(2*k+1,2*l+1)-d(2*k+1,2*l+1) )
        ap(2*k+2,2*l+1)=
     .(1./(XK**2*k**2+XL**2*l**2+F1))
     .*( BB*XK*k*a(2*k+1,2*l+1)
     .   +c(2*k+2,2*l+1)-d(2*k+2,2*l+1) )
        ap(2*k+1,2*l+2)=
     .(1./(XK**2*k**2+XL**2*l**2+F1))
     .*( -BB*XK*k*a(2*k+2,2*l+2)
     .   +c(2*k+1,2*l+2)-d(2*k+1,2*l+2) )
        ap(2*k+2,2*l+2)=
     .(1./(XK**2*k**2+XL**2*l**2+F1))
     .*( BB*XK*k*a(2*k+1,2*l+2)
     .   +c(2*k+2,2*l+2)-d(2*k+2,2*l+2) )
100    continue
200   continue

      return
      end
c
      subroutine jacobc(nn,mm,a,fcap,c)
c  Jacobian for psi=sum(Akl*exp(ikx+ily))
      parameter (N=256,M=128)
      real a(nn,mm),c(2*nn-2,2*mm-2)
      common /dumjac/ u(2*N+2,2*M+2),v(2*N+2,2*M+2),qx(2*N+2,2*M+2)
     .,               qy(2*N+2,2*M+2)
      common /workft/ work((2*N+1)*2*M),trigs(3*N+1),ifax(13)
      common /workft2/ work2((2*M+1)*2*(N+2)),trigs2(3*M+1),ifax2(13)
c
      common /parms/ xk,xl,dt,bb,f1,rob,twopi,dtn,c0_rf
c
c  Get u,v,qx,qy wave terms at (k,l)
      do l=0,M/2
       do k=0,N/2
        u(2*k+1,2*l+1)=+XL*l*a(2*k+1,2*l+2)
        u(2*k+1,2*l+2)=-XL*l*a(2*k+1,2*l+1)
        u(2*k+2,2*l+1)=+XL*l*a(2*k+2,2*l+2)
        u(2*k+2,2*l+2)=-XL*l*a(2*k+2,2*l+1)
        v(2*k+1,2*l+1)=-XK*k*a(2*k+2,2*l+1)
        v(2*k+2,2*l+1)=+XK*k*a(2*k+1,2*l+1)
        v(2*k+1,2*l+2)=-XK*k*a(2*k+2,2*l+2)
        v(2*k+2,2*l+2)=+XK*k*a(2*k+1,2*l+2)
        qx(2*k+1,2*l+1)=
     . -XK*k*(-(XK**2*k**2+XL**2*l**2+fcap)*a(2*k+2,2*l+1))
        qx(2*k+2,2*l+1)=
     . +XK*k*(-(XK**2*k**2+XL**2*l**2+fcap)*a(2*k+1,2*l+1))
        qx(2*k+1,2*l+2)=
     . -XK*k*(-(XK**2*k**2+XL**2*l**2+fcap)*a(2*k+2,2*l+2))
        qx(2*k+2,2*l+2)=
     . +XK*k*(-(XK**2*k**2+XL**2*l**2+fcap)*a(2*k+1,2*l+2))
        qy(2*k+1,2*l+1)=
     . -XL*l*(-(XK**2*k**2+XL**2*l**2+fcap)*a(2*k+1,2*l+2))
        qy(2*k+1,2*l+2)=
     . +XL*l*(-(XK**2*k**2+XL**2*l**2+fcap)*a(2*k+1,2*l+1))
        qy(2*k+2,2*l+1)=
     . -XL*l*(-(XK**2*k**2+XL**2*l**2+fcap)*a(2*k+2,2*l+2))
        qy(2*k+2,2*l+2)=
     . +XL*l*(-(XK**2*k**2+XL**2*l**2+fcap)*a(2*k+2,2*l+1))
       enddo
      enddo
      do l=1,M+2
       do k=N+3,2*N+2
        u(k,l)=0.
        v(k,l)=0.
        qx(k,l)=0.
        qy(k,l)=0.
       enddo
      enddo
      do l=M+3,2*M+2
       do k=1,2*N+2
        u(k,l)=0.
        v(k,l)=0.
        qx(k,l)=0.
        qy(k,l)=0.
       enddo
      enddo 
c
c
c  Transform from (k,l) to (x,l) and then to (x,y)
      inc=1
      jump=2*N+2
      call fft991(u,work,trigs,ifax,inc,jump,2*N,M+2,+1)
      call fft991(v,work,trigs,ifax,inc,jump,2*N,M+2,+1)
      call fft991(qx,work,trigs,ifax,inc,jump,2*N,M+2,+1)
      call fft991(qy,work,trigs,ifax,inc,jump,2*N,M+2,+1)
      inc=2*N+2
      jump=1
      call fft991(u,work2,trigs2,ifax2,inc,jump,2*M,2*N,+1)
      call fft991(v,work2,trigs2,ifax2,inc,jump,2*M,2*N,+1)
      call fft991(qx,work2,trigs2,ifax2,inc,jump,2*M,2*N,+1)
      call fft991(qy,work2,trigs2,ifax2,inc,jump,2*M,2*N,+1)
c  Compute jacobian on doubled grid
      do l=1,2*M
       do k=1,2*N
        c(k,l)=u(k,l)*qx(k,l)+v(k,l)*qy(k,l)
       enddo
      enddo
c  Transform from (x,y) to (x,l) and then to (k,l)
c   Could repack to avoid transforming zeros at k=2
      inc=2*N+2
      jump=1
      call fft991(c,work2,trigs2,ifax2,inc,jump,2*M,2*N,-1)
      inc=1
      jump=2*N+2
      call fft991(c,work,trigs,ifax,inc,jump,2*N,M+2,-1)

c  Done
      return
      end
c
c
      subroutine out(nn,mm,a,p1,iun)
      parameter (N=256,M=128)
      real a(nn,mm),p1(nn,mm)
      real u(2*N+2,M),v(2*N+2,M),q(2*N+2,M)
      common /workft/ work((2*N+1)*2*M),trigs(3*N+1),ifax(13)
      common /workft2/ work2((2*M+1)*2*(N+2)),trigs2(3*M+1),ifax2(13)
c
      common /parms/ xk,xl,dt,bb,f1,rob,twopi,dtn,c0_rf
c
      do 20 l=1,M+2
       do 10 k=1,N+2
        p1(k,l)=a(k,l)
10     continue
20    continue
c
      call set99(trigs,ifax,N)
      call set99(trigs2,ifax2,M)
      inc=1
      jump=N+2
      call fft991(p1,work,trigs,ifax,inc,jump,N,M,+1)
      inc=N+2
      jump=1
      call fft991(p1,work2,trigs2,ifax2,inc,jump,M,N,+1)
      call set99(trigs,ifax,2*N)
      call set99(trigs2,ifax2,2*M)
c
      Mout=M/16
c      Mout=min(M,22)
      do 40 k=1,N
       do lm=1,Mout
        la=1+(lm-1)*16
	lb=lm*16
        write(iun,'(1x,i5,20(1x,e10.4))') k,(p1(k,l),l=la,lb)
       enddo
40    continue 
      return
      end


      subroutine intclc(nn,mm,a,f_a,ep,ek,eq,ef)
c  Calculate integral quantities
c   for psi=sum(Akl*exp(ikx)*sinly)+sum*(A0l*cosly)
c      ep = psi^2
c      ek = K^2 psi^2
c      eq = K^4 psi^2
c      ef = force_a^2
c    
      parameter (N=256,M=128)
c
      common /parms/ xk,xl,dt,bb,f1,rob,twopi,dtn,c0_rf
c
      real a(nn,mm),f_a(nn,mm)
c  Calculate spectral sums
      ep=0.
      ek=0.
      eq=0.
      ef=0.
c  Get u,v,qx,qy wave terms at (k,l)
      do l=0,M/2
       do k=0,N/2
        ep=ep+a(2*k+1,2*l+1)**2+a(2*k+1,2*l+2)**2
     .       +a(2*k+2,2*l+1)**2+a(2*k+2,2*l+2)**2
        ek=ek+(XK**2*k**2+XL**2*l**2)
     .        *(a(2*k+1,2*l+1)**2+a(2*k+1,2*l+2)**2
     .         +a(2*k+2,2*l+1)**2+a(2*k+2,2*l+2)**2)
        eq=eq+(XK**2*k**2+XL**2*l**2)**2
     .        *(a(2*k+1,2*l+1)**2+a(2*k+1,2*l+2)**2
     .         +a(2*k+2,2*l+1)**2+a(2*k+2,2*l+2)**2)
        ef=ef+f_a(2*k+1,2*l+1)**2+f_a(2*k+1,2*l+2)**2
     .       +f_a(2*k+2,2*l+1)**2+f_a(2*k+2,2*l+2)**2
       enddo
      enddo

      return
      end


      subroutine init(ki,nn,mm,a0,a1,ap,icont)
c  Initialization
      parameter (N=256,M=128)
      real a0(nn,mm),a1(nn,mm)
      real ap(nn,mm)
      real ca(M),cb(M)
      common /dumdif/ ddum((2*N+2)*(2*M+2)+(N+2)*(M+2))
      common /dumjac/ dumj(4*(2*N+2)*(2*M+2))
      common /workft/ work((2*N+1)*2*M),trigs(3*N+1),ifax(13)
      common /workft2/ work2((2*M+1)*2*(N+2)),trigs2(3*M+1),ifax2(13)
c
      common /parms/ xk,xl,dt,bb,f1,rob,twopi,dtn,c0_rf
c
c
      twopi=8*atan(1.0)
      write(6,*) twopi
c  zero initial stream function arrays
       do l=1,M+2
        do k=1,N+2
         a0(k,l)=0.
         a1(k,l)=0.
        enddo
       enddo
c  single waves
c       a0(2*6+1,2*1+1)=0.05
c       a0(2*6+1,2*3+1)=0.05
c
c  Read input file in (k,l) space
c       Mout=M/16
c      Mout=min(M,22)
c       do k=1,N
c        do lm=1,Mout
c         la=1+(lm-1)*16
c         lb=lm*11
c         read(5,'(1x,i5,20(1x,e10.4))') kk,(a0(k,l),l=la,lb)
c        enddo
c       enddo
c       
c
c  Zero arrays
      do i=1,(2*N+2)*(2*M+2)+(N+2)*(M+2)
       ddum(i)=0.
      enddo
      do i=1,4*(2*N+2)*(2*M+2)
       dumj(i)=0.
      enddo
c
c
      call set99(trigs,ifax,2*N)
      call set99(trigs2,ifax2,2*M)
c
c  Disturbance at step -1 for leapfrog
c
c      call diff(N+2,M+2,t,a0,a0,ap)
c      do 90 l=1,M+2
c       do 80 k=1,N+1
c        a1(k,l)=a0(k,l)-ap(k,l)*DT
c80     continue
c      a1(N+2,l)=0.
c90    continue

      return
      end


! ----------------------------------------------------------------
!
      subroutine rkstep(nm,t,a0,a1,ap)
c rkstep:  2nd-order Runge-Kutta step
c   Adapt code from pgom for rk step, using eustpr euler step routine
c  One R-K time-step
      parameter (N=256,M=128)
c
      common /parms/ xk,xl,dt,bb,f1,rob,twopi,dtn,c0_rf
c
      real a0(nm),a1(nm)
      real ap(nm)
c      real al0(nn,mm),al1(nn,mm)
c      real alp(nn,mm)

c Predictor:
      call diff(N+2,M+2,t,a0,a0,ap)
c  predictor step
      call eustep((N+2)*(M+2),dt,a0,a1,ap)

c Corrector:
      call diff(N+2,M+2,t,a1,a0,ap)
c  corrector step
      call eustep((N+2)*(M+2),dt,a0,a0,ap)

c Predictor and corrector steps completed. Now average results:
      do i=1,(N+2)*(M+2)
       a0(i)=0.5*(a0(i)+a1(i))
      enddo

      t=t+dt
      return
      end

      subroutine eustep(nm,dt,a0,a1,ap)
! This routine does a simple Euler timestep
! Input variables: nm,dt,a0,ap
! Variable changed on output: a1
      real a0(nm),a1(nm)
      real ap(nm)
!  Euler time step
      do i=1,nm
       a1(i)=a0(i)+ap(i)*dt
      enddo
      return
      end
!^M
! ----------------------------------------------------------------^M


      subroutine out_nc(ncid,nrec,nn,mm,psi_varid,Fcap0_varid
     .,                 a0,force_a,p1)
      implicit none
c   QG spectral model dimensions
      integer N,M
      parameter (N=256,M=128)
      character*(*) FILE_NAME
      parameter (FILE_NAME = 'qg1ldp_rf_out.nc')
      integer ncid,nrec
      integer psi_varid, Fcap0_varid
c
      integer nn,mm
      real a0(nn,mm),force_a(nn,mm),p1(nn,mm)
      common /workft/ work((2*N+1)*2*M),trigs(3*N+1),ifax(13)
      common /workft2/ work2((2*M+1)*2*(N+2)),trigs2(3*M+1),ifax2(13)
      real work,trigs,work2,trigs2
      integer ifax,ifax2
      integer k,l,inc,jump

C   Output is 3D data: an MxN lat-lon grid, for one timestep
      integer NDIMS
      parameter (NDIMS = 3)

      include 'netcdf.inc'
      integer NLATS, NLONS
      parameter (NLATS = M, NLONS = N)

C   The start and count arrays will tell the netCDF library where to
C   write our data.
      integer start(NDIMS), count(NDIMS)

C   Error handling.
      integer retval

c -----------
      retval = nf_open(FILE_NAME, nf_write, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)


C   These settings tell netcdf to write one timestep of data. (The
C   setting of start(3) inside the loop below tells netCDF which
C   timestep to write.)
      count(1) = NLONS
      count(2) = NLATS
      count(3) = 1
      start(1) = 1
      start(2) = 1
      start(3) = nrec

C   Write the stream function and forcing arrays.
c   The arrays hold one time step of data.
c
      do l=1,M+2
       do k=1,N+2
        p1(k,l)=a0(k,l)
       enddo
      enddo
c
      call set99(trigs,ifax,N)
      call set99(trigs2,ifax2,M)
      inc=1
      jump=N+2
      call fft991(p1,work,trigs,ifax,inc,jump,N,M,+1)
      inc=N+2
      jump=1
      call fft991(p1,work2,trigs2,ifax2,inc,jump,M,N,+1)
      call set99(trigs,ifax,2*N)
      call set99(trigs2,ifax2,2*M)
c
      retval = nf_put_vara_real(ncid, psi_varid, start, count, p1)
      if (retval .ne. nf_noerr) call handle_err(retval)
c
      do l=1,M+2
       do k=1,N+2
        p1(k,l)=force_a(k,l)
       enddo
      enddo
c
      call set99(trigs,ifax,N)
      call set99(trigs2,ifax2,M)
      inc=1
      jump=N+2
      call fft991(p1,work,trigs,ifax,inc,jump,N,M,+1)
      inc=N+2
      jump=1
      call fft991(p1,work2,trigs2,ifax2,inc,jump,M,N,+1)
      call set99(trigs,ifax,2*N)
      call set99(trigs2,ifax2,2*M)
c
      retval = nf_put_vara_real(ncid, Fcap0_varid, start, count, p1)
      if (retval .ne. nf_noerr) call handle_err(retval)

C   Close the file. This causes netCDF to flush all buffers and make
C   sure your data are really written to disk.
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      print *,'*** SUCCESS writing state file ',nrec

      end

! ----------------------------------------------------------------^M



      subroutine out_nc_init(ncid
     .,  bb,f1,xk,xl,dt,ifrict,rpsi,rek,r0_del6,r,tau,f_amp
     .,  dxkf,xkf,c0_rf,ndt,dtn,rcap
     .,  nn,mm,psi_varid,Fcap0_varid,a0,force_a,p1)
      implicit none
c   QG spectral model dimensions
      integer N,M
      parameter (N=256,M=128)
C   Specify output file name
      character*(*) FILE_NAME
      parameter (FILE_NAME = 'qg1ldp_rf_out.nc')
      integer ncid
c
      integer nn,mm
      real a0(nn,mm),force_a(nn,mm),p1(nn,mm)
      common /workft/ work((2*N+1)*2*M),trigs(3*N+1),ifax(13)
      common /workft2/ work2((2*M+1)*2*(N+2)),trigs2(3*M+1),ifax2(13)
      real work,trigs,work2,trigs2
      integer ifax,ifax2
      integer k,l,inc,jump

C   Output is 3D data: an MxN lat-lon grid, with 1
C   timestep of initial data; also numerous model parameters.
      integer NDIMS, NRECS
      parameter (NDIMS = 3, NRECS = 1)

      include 'netcdf.inc'
      integer NLATS, NLONS
      parameter (NLATS = M, NLONS = N)
      character*(*) LAT_NAME, LON_NAME, REC_NAME
      parameter (LAT_NAME = 'latitude', LON_NAME = 'longitude')
      parameter (REC_NAME = 'time')
      integer lon_dimid, lat_dimid, rec_dimid

C   The start and count arrays will tell the netCDF library where to
C   write our data.
      integer start(NDIMS), count(NDIMS)

C   These program variables hold the latitudes and longitudes.
      real lats(NLATS), lons(NLONS)
      integer lon_varid, lat_varid
      real twopi

C   Create netCDF variables for stream function psi and
C   forcing Fcap0 fields.
      character*(*) psi_NAME, Fcap0_NAME
      parameter (psi_NAME='psi')
      parameter (Fcap0_NAME='Fcap0')
      integer psi_varid, Fcap0_varid
      integer dimids(NDIMS)
c   QG spectral model parameters
      real bb,f1,xk,xl,dt,rpsi,rek,r0_del6,r,tau,f_amp
     .,  dxkf,xkf,c0_rf,dtn,rcap
      integer ifrict,ndt
C   Create a netCDF variable for each model parameter
      character*(*) bb_NAME, f1_NAME, xk_NAME, xl_NAME, dt_NAME
     ., rpsi_NAME, rek_NAME, r0_del6_NAME, r_NAME, tau_NAME, f_amp_NAME
     ., dxkf_NAME, xkf_NAME, c0_rf_NAME, ndt_NAME, dtn_NAME, rcap_NAME
     ., ifrict_NAME
      parameter (bb_NAME='bb',f1_NAME='f1',xk_NAME='xk',xl_NAME='xl')
      parameter (dt_NAME='dt',rpsi_NAME='rpsi',rek_NAME='rek'
     .,          r0_del6_NAME='r0_del6')
      parameter (r_NAME='r',tau_NAME='tau',f_amp_NAME='f_amp'
     .,          dxkf_NAME='dxkf')
      parameter (xkf_NAME='xkf',c0_rf_NAME='c0_rf',ndt_NAME='ndt'
     .,          dtn_NAME='dtn')
      parameter (rcap_NAME='rcap',ifrict_NAME='ifrict')
      integer bb_varid, f1_varid, xk_varid, xl_varid, dt_varid
     .,       rpsi_varid
     .,       rek_varid, r0_del6_varid, r_varid, tau_varid, f_amp_varid
     .,       dxkf_varid, xkf_varid, c0_rf_varid, ndt_varid, dtn_varid
     .,       rcap_varid, ifrict_varid

C   Each variable carries a "units" attribute; all are dimensionless.
      character*(*) UNITS
      parameter (UNITS = 'Units')
      character*(*) UNITS_nd
      parameter (UNITS_nd = 'dimensionless')

C   Use these to construct latitude and longitude data
      integer START_LAT, START_LON
      parameter (START_LAT = 0, START_LON = 0)

C   Loop indices.
      integer lat, lon, rec, i

C   Error handling.
      integer retval

c -----------

C   Create the file. 
      write(6,*) FILE_NAME
      retval = nf_create(FILE_NAME, nf_clobber, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) FILE_NAME

c   Define parameter variables, assign units attributes.
      write(6,*) bb_NAME 
      retval = nf_def_var(ncid, bb_NAME, NF_REAL, 0, 0, bb_varid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) bb_NAME 
      retval = nf_put_att_text(ncid, bb_varid, UNITS, len(UNITS_nd)
     .,UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) bb_NAME 
      write(6,*) f1_NAME 
      retval = nf_def_var(ncid, f1_NAME, NF_REAL, 0, 0, f1_varid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, f1_varid, UNITS, len(UNITS_nd)
     .,UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) xk_NAME 
      retval = nf_def_var(ncid, xk_NAME, NF_REAL, 0, 0, xk_varid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, xk_varid, UNITS, len(UNITS_nd)
     .,UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) xl_NAME 
      retval = nf_def_var(ncid, xl_NAME, NF_REAL, 0, 0, xl_varid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, xl_varid, UNITS, len(UNITS_nd)
     .,UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) dt_NAME 
      retval = nf_def_var(ncid, dt_NAME, NF_REAL, 0, 0, dt_varid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, dt_varid, UNITS, len(UNITS_nd)
     .,UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) rpsi_NAME 
      retval = nf_def_var(ncid, rpsi_NAME, NF_REAL, 0, 0, rpsi_varid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, rpsi_varid, UNITS, len(UNITS_nd)
     .,UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) rek_NAME 
      retval = nf_def_var(ncid, rek_NAME, NF_REAL, 0, 0, rek_varid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, rek_varid, UNITS, len(UNITS_nd)
     .,UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) r0_del6_NAME 
      retval = nf_def_var(ncid, r0_del6_NAME, NF_REAL, 0, 0
     ., r0_del6_varid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, r0_del6_varid, UNITS, len(UNITS_nd)
     .,UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) r_NAME 
      retval = nf_def_var(ncid, r_NAME, NF_REAL, 0, 0, r_varid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, r_varid, UNITS, len(UNITS_nd)
     .,UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) tau_NAME 
      retval = nf_def_var(ncid, tau_NAME, NF_REAL, 0, 0, tau_varid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, tau_varid, UNITS, len(UNITS_nd)
     .,UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) f_amp_NAME 
      retval = nf_def_var(ncid, f_amp_NAME, NF_REAL, 0, 0, f_amp_varid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, f_amp_varid, UNITS, len(UNITS_nd)
     .,UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) dxkf_NAME 
      retval = nf_def_var(ncid, dxkf_NAME, NF_REAL, 0, 0, dxkf_varid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, dxkf_varid, UNITS, len(UNITS_nd)
     .,UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) xkf_NAME 
      retval = nf_def_var(ncid, xkf_NAME, NF_REAL, 0, 0, xkf_varid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, xkf_varid, UNITS, len(UNITS_nd)
     .,UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) c0_rf_NAME 
      retval = nf_def_var(ncid, c0_rf_NAME, NF_REAL, 0, 0, c0_rf_varid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, c0_rf_varid, UNITS, len(UNITS_nd)
     .,UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) ndt_NAME 
      retval = nf_def_var(ncid, ndt_NAME, NF_INT, 0, 0, ndt_varid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, ndt_varid, UNITS, len(UNITS_nd)
     .,UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) dtn_NAME 
      retval = nf_def_var(ncid, dtn_NAME, NF_REAL, 0, 0, dtn_varid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, dtn_varid, UNITS, len(UNITS_nd)
     .,UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) rcap_NAME 
      retval = nf_def_var(ncid, rcap_NAME, NF_REAL, 0, 0, rcap_varid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, rcap_varid, UNITS, len(UNITS_nd)
     .,UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(6,*) ifrict_NAME 
      retval = nf_def_var(ncid, ifrict_NAME, NF_INT, 0, 0, ifrict_varid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, ifrict_varid, UNITS, len(UNITS_nd)
     .,UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)

C   Define the dimensions for the spectral variables.
c   The record dimension is defined to have
C   unlimited length - it can grow as needed. In this case it is
C   the time dimension.
      retval = nf_def_dim(ncid, LAT_NAME, NLATS, lat_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, LON_NAME, NLONS, lon_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, REC_NAME, NF_UNLIMITED, rec_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C   Define the coordinate variables. We will only define coordinate
C   variables for lat and lon.  Ordinarily we would need to provide
C   an array of dimension IDs for each variable's dimensions, but
C   since coordinate variables only have one dimension, we can
C   simply provide the address of that dimension ID (lat_dimid) and
C   similarly for (lon_dimid).
      retval = nf_def_var(ncid, LAT_NAME, NF_REAL, 1, lat_dimid 
     .,    lat_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, LON_NAME, NF_REAL, 1, lon_dimid 
     .,    lon_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C   Assign units attributes to coordinate variables.
      retval = nf_put_att_text(ncid, lat_varid, UNITS, len(UNITS_nd) 
     .,   UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, lon_varid, UNITS, len(UNITS_nd) 
     .,   UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      if (retval .ne. nf_noerr) call handle_err(retval)

C   Create dimensionless lat-lon grid values
      twopi=8*atan(1.0)
      do lat = 1, NLATS
         lats(lat) = START_LAT + (lat - 1) * twopi/(xl*M)
      end do
      do lon = 1, NLONS
         lons(lon) = START_LON + (lon - 1) * twopi/(xk*N)
      end do

C   The dimids array is used to pass the dimids of the dimensions of
C   the netCDF variables. Both of the netCDF variables we are creating
C   share the same three dimensions. In Fortran, the unlimited
C   dimension must come last on the list of dimids.
      dimids(1) = lon_dimid
      dimids(2) = lat_dimid
      dimids(3) = rec_dimid

C   Define the netCDF variables for the stream function and forcing data.
      retval = nf_def_var(ncid, psi_NAME, NF_REAL, NDIMS, dimids 
     .,    psi_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, Fcap0_NAME, NF_REAL, NDIMS, dimids 
     .,    Fcap0_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C   Assign units attributes to the stream function and forcing variables.
      retval = nf_put_att_text(ncid, psi_varid, UNITS, len(UNITS_nd) 
     .,    UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, Fcap0_varid, UNITS, len(UNITS_nd) 
     .,    UNITS_nd)
      if (retval .ne. nf_noerr) call handle_err(retval)

C   End define mode.
      write(6,*) FILE_NAME
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

c  write parameter values
      retval = nf_put_var_real(ncid, bb_varid, bb)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, f1_varid, f1)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, xk_varid, xk)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, xl_varid, xl)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, dt_varid, dt)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, rpsi_varid, rpsi)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, rek_varid, rek)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, r0_del6_varid, r0_del6)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, r_varid, r)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, tau_varid, tau)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, f_amp_varid, f_amp)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, dxkf_varid, dxkf)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, xkf_varid, xkf)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, c0_rf_varid, c0_rf)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_int(ncid, ndt_varid, ndt)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, dtn_varid, dtn)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, rcap_varid, rcap)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_int(ncid, ifrict_varid, ifrict)

C   Write the coordinate variable data. This will put the latitudes
C   and longitudes of our data grid into the netCDF file.
      retval = nf_put_var_real(ncid, lat_varid, lats)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, lon_varid, lons)
      if (retval .ne. nf_noerr) call handle_err(retval)

C   These settings tell netcdf to write one timestep of data. (The
C   setting of start(4) inside the loop below tells netCDF which
C   timestep to write.)
      count(1) = NLONS
      count(2) = NLATS
      count(3) = 1
      start(1) = 1
      start(2) = 1

C   Write the stream function and forcing arrays.
c   The arrays only hold the initial data.
c
      do l=1,M+2
       do k=1,N+2
        p1(k,l)=a0(k,l)
       enddo
      enddo
c
      call set99(trigs,ifax,N)
      call set99(trigs2,ifax2,M)
      inc=1
      jump=N+2
      call fft991(p1,work,trigs,ifax,inc,jump,N,M,+1)
      inc=N+2
      jump=1
      call fft991(p1,work2,trigs2,ifax2,inc,jump,M,N,+1)
      call set99(trigs,ifax,2*N)
      call set99(trigs2,ifax2,2*M)
c
      do rec = 1, NRECS
         start(3) = rec
         retval = nf_put_vara_real(ncid, psi_varid, start, count, p1)
         if (retval .ne. nf_noerr) call handle_err(retval)
      end do
c
      do l=1,M+2
       do k=1,N+2
        p1(k,l)=force_a(k,l)
       enddo
      enddo
c
      call set99(trigs,ifax,N)
      call set99(trigs2,ifax2,M)
      inc=1
      jump=N+2
      call fft991(p1,work,trigs,ifax,inc,jump,N,M,+1)
      inc=N+2
      jump=1
      call fft991(p1,work2,trigs2,ifax2,inc,jump,M,N,+1)
      call set99(trigs,ifax,2*N)
      call set99(trigs2,ifax2,2*M)
c
      do rec = 1, NRECS
         start(3) = rec
         retval = nf_put_vara_real(ncid, Fcap0_varid, start, count, p1)
         if (retval .ne. nf_noerr) call handle_err(retval)
      end do

C   Close the file. This causes netCDF to flush all buffers and make
C   sure your data are really written to disk.
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
   
      print *,'*** SUCCESS writing initial state file: '
     ., FILE_NAME
      end

      subroutine handle_err(errcode)
      implicit none
      include 'netcdf.inc'
      integer errcode

      print *, 'Error: ', nf_strerror(errcode)
      stop 2
      end
