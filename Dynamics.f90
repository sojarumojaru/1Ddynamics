
Program qdensity
 

 implicit none
 include 'fftw3.f'
 include 'mkl.fi'
 integer, parameter :: WP = selected_real_kind(15,307)
 Integer *8  plan, planb, counter, i,j, ierr, flag, loopcount, flag2, leftfluxi, rightfluxi
 Complex(WP), allocatable :: U(:,:,:), CCstore(:,:,:)
 !Integer *8, allocatable :: plan(:), planb(:)
 complex(WP), allocatable :: T(:)
 Integer nstates
 complex(WP), allocatable :: Psi(:,:), Psif(:,:), dummypsi(:), dummypsif(:)
 complex(WP) II, leftlow, lefthigh, rightlow, righthigh, temp1
 Real*8, allocatable :: CC(:,:),CC1(:,:),V(:,:)
 Real*8, dimension(26):: plist
 Real*8, allocatable :: gtx(:), gtk(:), d(:,:,:)
 Complex(WP), allocatable :: Vd(:,:), Vp(:,:), Vpd(:,:)
 Real*8, allocatable :: l(:)
 Complex*16, allocatable :: density4SVD(:,:), density4SVDc(:,:), gamma_e(:,:)
 complex(WP), allocatable :: wvadia(:)
 Real *8 xmin, xmax, hx, ht, p, q, e1, e2, m, sigma, time, tmax, pi, &
  & hk, ck, k, kin, Xx, centre, b, u0, alpha, gmma, cleftlow, clefthigh, crightlow, crighthigh,  sums, scndd,&
  & leftflux, rightflux, centralnorm, cutoff, readvariable
 Real *8, allocatable :: momentum(:), Energy(:), normk(:), normx(:), psisq(:), expx(:),&
  & norml(:), normr(:), decay(:)
 !**************The external call function variables (i.e. things required for zdgemm and zgesvd
 Integer N4svd, M4svd, INFO, LWORK, LDA, LDU, LDVT, ngrid, statecount
 Real*8, Allocatable :: S4svd(:), RWORK(:)
 complex(WP), allocatable :: U4svd(:,:), VT4SVD(:,:), work(:)
 Character fdummy
 !** Look at https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/zgesvd_ex.f.htm
 Complex*16 alphadgemm, betadgemm
 Integer ldadgemm, ldbdgemm, ldcdgemm, Ndgemm, Mdgemm, kdgemm


           !**********************************************************************************************!
 Open(2, File = 'Input.in')
 Open(3, File = 'norm') 
 Open(7, File= 'SV.dat')  
 Open(8, File = 'SVV.dat')
 Open(9, File = 'Wavefunc.dat')
 Open(6, File = 'rho_12.dat')
 !*************Read inputs***************
 read(2,*) fdummy, readvariable
 nstates = readvariable

 read(2,*) fdummy, readvariable
 ngrid = readvariable

 read(2,*) fdummy, readvariable
 cutoff = readvariable

 read(2,*) fdummy, readvariable
 u0 = readvariable
!This is for tha absorbing boundary condition
 read(2,*) fdummy, readvariable
 alpha = readvariable
!This is for the absorbing boundary condition
 read(2,*) fdummy, readvariable
 kin = readvariable

 read(2,*) fdummy, readvariable
 m = readvariable

 read(2,*) fdummy, readvariable
 ht = readvariable
 !******************set parameters and initialize**********************
 write(*,*) 'Successfully read the file'
 II= (0,1.0)
 pi = 3.1415926535
 N4SVD=nstates
 M4SVD=ngrid
 LDA = M4SVD
 LDU = M4SVD
 LDVT =N4SVD
 LWORK = -1
 allocate(wvadia(nstates), Vd(nstates,nstates), Vp(nstates, nstates), Vpd(nstates,nstates))
 allocate(CC(nstates, nstates), CC1(nstates,nstates), V(nstates,nstates), l(nstates))

 allocate(U4SVD(ngrid,ngrid), VT4SVD(nstates,nstates),WORK(1000),S4svd(nstates),RWORK(5*nstates))
 allocate(momentum(nstates), Energy(nstates), normk(nstates), normx(nstates), psisq(nstates), expx(nstates),&
  & norml(nstates), normr(nstates), decay(nstates))
 allocate(Psi(ngrid,nstates),Psif(ngrid,nstates),U(nstates,nstates,ngrid),d(nstates,nstates,ngrid),T(ngrid),gtx(ngrid),gtk(ngrid),&
  & CCstore(ngrid,nstates,nstates), density4SVD(ngrid,nstates), density4SVDc(ngrid,nstates), gamma_e(nstates,nstates),&
  &dummypsi(ngrid), dummypsif(ngrid))
 if (ierr.ne.0) write(*,*) 'Bless the lord!'
 gamma_e=0.0 
 Ndgemm=nstates
 Mdgemm=nstates
 kdgemm=ngrid
 ldadgemm = kdgemm
 ldbdgemm = kdgemm
 ldcdgemm = Mdgemm
 alphadgemm =1.0
 betadgemm=0.0
 !gamma_e is the c, a = density4SVD^t, b = density4SVD
 !do j = 1,26 

 !kin=plist(j)
 plan =0
 planb=0
 flag=0
 flag2=0
 psi=0.0
 psif=0.0
 centralnorm=0
 sigma =20.0/abs(kin)

 normx=0
 expx = 0
 momentum= 0
 Energy = 0
 
 norml = 0
 normr = 0

 k=0

 xmax= cutoff+min(7.0*sigma,30.0)+30.0
 xmin=-8.0*sigma-cutoff-10.0

 hx=(xmax-xmin)/ngrid
 hk =2.0*pi/(xmax-xmin)

 centre = -4.0*sigma-cutoff
 cleftlow=0
 clefthigh=0
 crightlow=0 
 crighthigh=0
 leftflux = centre-3*sigma
 rightflux= xmax-sigma
 
 leftfluxi = int((leftflux-xmin)/hx)
 rightfluxi = int((rightflux-xmin)/hx)
 write(*,*) 'Beginning the propagation for k = ', kin

 !counter = 0


 !**The U and T matrix being initialised as well as the mapping of grid to space **!
 !**This section is specific to this Model ONLY** ideally should be a subroutine***!
 do i = 1, ngrid
        gtx(i) = xmin+((i-1)*hx)                                                   !            
        !write(*,*) gtx(i)                                                         !
        if (i<=((ngrid/2.0)+1)) then                                               ! 
              gtk(i)=dble(i-1)*hk                                                  !        
        else if(i>((ngrid/2.0)+1)) then                                         
              gtk(i)=dble(i-ngrid-1)*hk
        end if                                                                     !         
        T(i) = exp(-ht*II*gtk(i)**2/(2*m))                                         !  
 end do                                                                            !          
                                                                                   !     
 do counter = 1,ngrid                                                              !                      
        !counter=counter+1                                                         !       
        call Potential(gtx(counter),V)                                             !                 
        e1= V(1,1)+V(2,2)                                                          !             
        e2 = (V(1,1)-V(2,2))**2+4*V(2,1)*V(1,2)                                    !    
        l(2) = (e1+Sqrt(e2))/2 !The two eigen values are evaluated                 ! 
        l(1) = (e1-Sqrt(e2))/2                                                     !
        Vd=0.0                                                                     !    
        Vd(1,1) = l(1)                                                             !                
        Vd(2,2)=l(2)                                                               !      
        call createEigenVector(l(1), l(2), CC, V)!The eigen vector matrix is formed!   
        call inverse2(CC,CC1)                                                      !              
        Vp(1,1) = exp(-0.5*II*l(1)*ht)                                             !       
        CCstore(counter, 1 , 1) =CC(1,1)                                           !           
        CCstore(counter, 1 , 2) =CC(1,2)                                           !             
        CCstore(counter, 2 , 1) =CC(2,1)                                           !      
        CCstore(counter, 2 , 2) =CC(2,2)                                           !
        Vp(2,2) = exp(-0.5*II*l(2)*ht)                                             !         
        Vpd = matmul(CC1,matmul(Vp,CC))                                            !            
        U(1,1,counter) = Vpd(1,1)                                                  !                  
        U(1,2,counter) = Vpd(1,2)                                                  !                      
        U(2,1,counter) = Vpd(2,1)                                                  !                                                        
        U(2,2,counter) = Vpd(2,2)                                                  !                         
        !write(7,*) q, V, Vpd                                                      !                     
                                                                                   !             
 end do                                                                            !    
 k=0                                                                               !
 !*****Saswata, I hope you are smart enough to read this portion and realise that this 
 !*****portion needs a subroutine call. Prefereably externally************************

 !***************************Initialising the populations**************************!
 normx=0         
 counter = 0 
 do counter = 1,ngrid
     
        Xx = gtx(counter)-centre
        b=2.0/(sigma*sigma)
        psi(counter,1) = exp(-Xx*Xx/(sigma*sigma))*exp(II*kin*Xx)*sqrt(sqrt(b/pi))
        psi(counter,2)= 0 
        do statecount = 1,nstates
              psisq(statecount) = (abs(psi(counter,statecount))**2)
              normx(statecount) =normx(statecount)+psisq(statecount)*hx
              expx(statecount) = expx(statecount)+gtx(counter)*psisq(statecount)*hx
        end do
        write(9,'(4E15.6)') gtx(counter), abs(psi(counter,1))**2, abs(psi(counter,2))**2, 0!<--.
        !                                                                                      :
        !This step is specific for this Model ONLY****************_____________________________;
        !write(8,*) counter, ck , exp(ck)
  end do

 write(9,*) ''

 expx=0
 normx=0

!************************************************************************************!

 call dfftw_plan_dft_1d(plan, ngrid, dummypsi, dummypsif, FFTW_FORWARD,FFTW_ESTIMATE)
 call dfftw_plan_dft_1d(planb, ngrid, dummypsif, dummypsi, FFTW_BACKWARD,FFTW_ESTIMATE)

!**********************************Propagate*****************************************!

 norml = 0
 normr = 0
 !norml2 = 0
 !normr2 = 0
   
 loopcount=0

 time = 0
 !tmax = 5
 tmax =4.0*abs(centre*m/kin)
 write(*,*) 'tmax = ', tmax
 !do while(time<tmax)
 do while(flag2.ne.1)
  
!************************************Again this Operation is specific to a two state system*****************%*    
        do counter =1,ngrid
              temp1 = U(1,1, counter)*psi(counter,1) + U(1,2, counter)*psi(counter,2)
              psi(counter,2) = U(2,1, counter)*psi(counter,1) + U(2,2, counter)*psi(counter,2)
              psi(counter,1) = temp1
        end do
!***********************************************************************************************************%    
        do statecount = 1, nstates
              do counter = 1,ngrid
                    dummypsi(counter)=psi(counter, statecount)
              end do
              call dfftw_execute_dft(plan,dummypsi, dummypsif)

              do counter = 1,ngrid
                    psif(counter,statecount)=dummypsif(counter)
              end do
        end do
              !call dfftw_execute_dft(plan2,psi2, psi2f)

        do counter = 1, ngrid
              do statecount = 1,nstates

                    psif(counter,statecount) = T(counter)*psif(counter,statecount)
                    k=0            
                    k=gtk(counter)
              
                    psisq(statecount) = abs(psif(counter,statecount)**2)
                    momentum(statecount) = momentum(statecount)+k*psisq(statecount)*hk
                    Energy(statecount) = Energy(statecount) + (((k)**2)/(2*m))*psisq(statecount)*hk
                    normk(statecount) = normk(statecount)+psisq(statecount)*hk
              end do 
        end do
        
        do statecount = 1, nstates
              do counter = 1,ngrid
                    dummypsif(counter)=psif(counter, statecount)
              end do
              call dfftw_execute_dft(planb,dummypsif, dummypsi)

              do counter = 1,ngrid
                    psi(counter,statecount)=dummypsi(counter)
              end do
        end do
        !call dfftw_execute_dft(planb, psi1f, psi1)
        !call dfftw_execute_dft(planb2, psi2f, psi2)
          
       
        do counter =1,ngrid
              temp1 = (U(1,1, counter)*psi(counter,1) + U(1,2, counter)*psi(counter,2))/dble(ngrid)
              psi(counter,2) = (U(2,1, counter)*psi(counter,1) + U(2,2, counter)*psi(counter,2))/dble(ngrid)
              psi(counter,1) = temp1

              gmma = u0/(cosh(alpha*abs(gtx(counter)-xmax)))**2+ u0/(cosh(alpha*abs(gtx(counter)-xmin)))**2
              do statecount = 1,nstates
                    decay(statecount) = gmma*ht*psi(counter,statecount)
                    !decay2 = gmma*ht*psi2(counter)

                    psi(counter,statecount) = (1-gmma*ht)*psi(counter,statecount)
                    !psi2(counter) = (1-gmma*ht)*psi2(counter)
              

                    psisq(statecount) = (abs(psi(counter,statecount))**2)
                    normx(statecount)=normx(statecount)+psisq(statecount)*hx
                    expx(statecount) = expx(statecount)+gtx(counter)*psisq(statecount)*hx
              
                    !psisq2 = (abs(psi2(counter))**2)


                    !normx2=normx2+psisq2*hx

                    !expx2 = expx2+gtx(counter)*psisq2*hx


                    if((gtx(counter)>=-cutoff).AND.gtx(counter)<=cutoff) then
                
                           centralnorm=centralnorm +(psisq(statecount))*hx

                    end if
              end do
        end do

        if(mod(loopcount,500)==0) then
              flag=1
              do counter = 1, ngrid
                    do statecount = 1,nstates
                           wvadia(statecount) =psi(counter,statecount)
                    end do
                    do i = 1,nstates
                           do j = 1,nstates 
                                  CC(i,j) = CCstore(counter, i , j)
                           end do
                    end do
                    wvadia = matmul(CC,wvadia)
                  
                    do statecount = 1,nstates
                           density4SVD(counter,statecount) =wvadia(statecount)
                           density4SVDc(counter,statecount) =wvadia(statecount)
                    end do

!*******************************This is another statement specific to my system
                    write(9,'(4E15.6)') gtx(counter) , abs(wvadia(1))**2, abs(wvadia(2))**2, centralnorm
                    write(3,'(3E15.6)') time, normx(1)+normx(2), normx(1)
              end do
        end if

        LWORK = 1000
        if( flag==1) then
 
              call ZGESVD('A','A', ngrid, nstates, density4SVDc,ngrid,S4SVD,U4SVD, ngrid, VT4SVD,nstates,&
               &WORK,LWORK,RWORK,INFO)

              if(INFO.ne.0) write(*,*) 'Did not converge'

              write(7,'(3E15.6)') time, S4SVD(1)*sqrt(hx), S4SVD(2)*sqrt(hx)


              call ZGEMM('C','N', nstates,nstates, ngrid, alphadgemm, density4SVD, ngrid, density4SVD, &
               &ngrid, betadgemm, gamma_e, nstates)
              gamma_e = gamma_e*hx
     
              write(6,'(5E15.6)') time, abs(gamma_e(1,1)), dreal(gamma_e(1,2)), dimag(gamma_e(2,1)), abs(gamma_e(2,2)) 
              gamma_e=0.0
              do counter = 1,ngrid
                    write(8, '(4E15.6)') time, gtx(counter), abs(U4SVD(counter,1)), abs(U4SVD(counter,2))
              end do
        end if
        momentum = momentum/normk
        Energy = Energy/normk
        expx = expx/normx
!************************Another Section specific to this 2 state model*******!!!!!        
        leftlow =((psi(leftfluxi,1)*(conjg(psi(leftfluxi,1))-conjg(psi(leftfluxi+1,1))))&
         &-(conjg(psi(leftfluxi,1))*(psi(leftfluxi,1)-psi(leftfluxi+1,1))))/hx!*scndd
        lefthigh =((psi(leftfluxi,2)*(conjg(psi(leftfluxi,2))-conjg(psi(leftfluxi+1,2))))&
         &-(conjg(psi(leftfluxi,2))*(psi(leftfluxi,2)-psi(leftfluxi+1,2))))/hx!*scndd
        rightlow =((psi(rightfluxi,2)*(conjg(psi(rightfluxi,2))-conjg(psi(rightfluxi+1,2))))&
         &-(conjg(psi(rightfluxi,2))*(psi(rightfluxi,2)-psi(rightfluxi+1,2))))/hx!*scndd
        righthigh =((psi(rightfluxi,1)*(conjg(psi(rightfluxi,1))-conjg(psi(rightfluxi+1,1))))&
         &-(conjg(psi(rightfluxi,1))*(psi(rightfluxi,1)-psi(rightfluxi+1,1))))/hx!*scndd

        cleftlow = cleftlow + abs(Real(leftlow*ht/(2*m*II)))
        clefthigh = clefthigh +abs(Real(lefthigh*ht/(2*m*II)))
        crightlow = crightlow + abs(Real(rightlow*ht/(2*m*II)))
        crighthigh = crighthigh + abs(Real(righthigh*ht/(2*m*II)))

      
        time=time+ht
        loopcount = loopcount+1

        if((time>=tmax).AND.(centralnorm<1E-4)) then
              flag2=1
              write(*,*) 'Terminating condition satisfied'
        end if

        momentum=0
        Energy = 0
        normx=0
        normk=0
        expx=0
        centralnorm=0       
        if(flag==1) then
              write(9,*) ''
              write(8,*) ''
              flag=0
        end if
 end do !The time loop
 write(*,'(4E15.4)') crighthigh, crightlow, cleftlow, clefthigh

 !***************************************************************************************************************
 do counter = 1,ngrid
        if(gtx(counter)<-cutoff) then 
        
              do statecount = 1, nstates 
                    psisq(statecount) = (abs(psi(counter,statecount))**2)

                    norml(statecount) = norml(statecount)+psisq(statecount)*hx
              end do 
        else if(gtx(counter)> cutoff) then
              do statecount = 1, nstates 
                    psisq(statecount) = (abs(psi(counter,statecount))**2)
              end do
              do statecount = 1,nstates
                    normr(statecount) = norml(statecount)+psisq(nstates-statecount+1)*hx
              end do 
      
        end if
 end do
   

 call dfftw_destroy_plan(plan) 
 call dfftw_destroy_plan(planb) 
 !end do ! The Kin loop
 close(7)
 close(6)
 close(8)
 close(9)

 deallocate(psi,psif, stat=ierr)
end

!**************************Subroutines**************!!

subroutine Potential(x,V)

 IMPLICIT none
!
!		declarations
 Real*8, intent(in) :: x
 Real*8, intent(inout),dimension(2,2) ::V

call PotentialSAC(x,V)
end subroutine
Subroutine Potentialdual(x, V)
    IMPLICIT none
!
!		de!larations
    REAL*8, intent(in) :: x
    Real*8, intent(inout),dimension(2,2) ::V
    Real*8 A, B, C, D, E0
    A= 0.10
    B=0.28
    C=0.015
    D=0.06
    E0 = 0.05
!
    
    V(1,1) = 0.0
    
    V(2,2) = -A*(exp(-B*x*x))+E0
    
    V(1,2) = C*(exp(-(D*x**2)))
    V(2,1)=V(1,2)
    

END subroutine Potentialdual 

Subroutine PotentialEX(x, V)
    IMPLICIT none
!
!		declarations
    REAL*8, intent(in) :: x
    Real*8, intent(inout),dimension(2,2) ::V
    Real*8 A, B,C
    A = 0.0006
    B = 0.10
    C = 0.90
!
    if(x<0) then 
      V(1,2) = B*(exp(C*x))
    else
      V(1,2) = B*(2- exp(-C*x))
    end if
    V(1,1) = A
    V(2,1)=V(1,2)
    V(2,2) = -V(1,1)

END subroutine PotentialEX 

Subroutine PotentialSAC(x, V)
 IMPLICIT none
!
!		declarations
 Real*8, intent(in) :: x
 Real*8, intent(inout),dimension(2,2) ::V
!
 if(x>0) then 
        V(1,1) = 0.01*(1- exp(-1.6*x))
 else
        V(1,1) = -0.01*(1- exp(1.6*x))
 end if
 V(1,2) = 0.005*(exp(-(x**2)))
 V(2,1)=V(1,2)
 V(2,2) = -V(1,1)

END subroutine PotentialSAC 

Subroutine flatPotential(x, V)
 IMPLICIT none
!
!		declarations
 REAL*8, intent(in) :: x
 Real*8, intent(inout),dimension(2,2) :: V
 V=0
end subroutine flatPotential

!********************************************************************************!

!********************************************************************************!

Subroutine qPotential(x, V)
 IMPLICIT none
!
!		declarations
 REAL*8, intent(in) :: x
 Real*8, intent(inout),dimension(2,2) :: V
 V=0
 V(1,1) =0.01*(x-2)*(x-2)
end subroutine qPotential

!*********************************************************************************!

Subroutine Potentialp(x, Vp)
 IMPLICIT none
!
!		de!larations
 REAL*8, intent(in) :: x
 Real*8, intent(inout),dimension(2,2) ::Vp
!
 if(x>0) then 
        Vp(1,1) = 0.01* exp(-1.6*x)*1.6
 else
        Vp(1,1) = 0.01*(1.6*exp(1.6*x))
 end if
 Vp(1,2) = -0.005*(exp(-(x**2)))*2*x
 Vp(2,1)=Vp(1,2)
 Vp(2,2) = -Vp(1,1)

END subroutine Potentialp 

!*********************************************************************************!
subroutine deriv(psi, pos, hx, ngrid, vali, xmin)
  implicit none
  integer, parameter :: WP = selected_real_kind(15,307)
  complex(WP), intent(in):: psi(:)
  Real*8, intent(in):: hx, pos, xmin
  Real*8, intent(inout) :: vali
  integer, intent(in):: ngrid
  integer posgrid

 ! allocate(psi(ngrid))
  posgrid =int((pos-xmin)/hx)
  vali = (-2*abs(psi(posgrid))**2+abs(psi(posgrid+1))**2+abs(psi(posgrid-1))**2)/(hx**2)
end subroutine deriv
  
subroutine createEigenVector(l1, l2, CC, V)
 implicit none
 Real*8 k1,k2,d1,d2
 Real*8, intent(in):: l1,l2
 Real*8, intent(inout), dimension(2,2):: CC, V
 if(abs(v(1,2))>1d-15) then
        k1 = (l1-V(1,1))/(V(1,2))
        k2 = (l2-V(1,1))/(V(1,2))
        d1 = Sqrt(1.0+k1**2)
        d2 = Sqrt(1.0+k2**2)

        CC(1,1) = 1.0/d1
        CC(1,2) = 1.0/d2
        CC(2,1) = k1/d1
        CC(2,2) = k2/d2
 else if (V(1,1)<V(2,2)) then
        CC(1,1) = 1.0
        CC(1,2) = 0.0
        CC(2,1) = 0.0
        CC(2,2) = 1.0
 else
        CC(1,1) = 0.0
        CC(1,2) = 1.0
        CC(2,1) = 1.0
        CC(2,2) = 0.0
    
 end if  
end subroutine createEigenVector

subroutine inverse2(A, B)
    implicit none

    Real*8, dimension(2,2), intent(in)::A
    Real*8, dimension(2,2), intent(inout)::B
    
    B= Transpose(A)    
end subroutine inverse2
