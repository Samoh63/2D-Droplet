        subroutine init()
        implicit none
        include 'head.inc'
        integer ix, iy, k, x, y
        real*8 a(lx,ly), b(lx,ly), feq(0:8)
        real*8 rho1, Re1
        !integer xc_sp,yc_sp ! Locations of sphere centers in domain
        !real*8 rpoint       ! distance of points from centers of spheres
        !integer ns,R          ! number of Spheres
        ! Read main parameters from file params.in
        open(1,file='params.in')
        read(1,*) ix, iy
        read(1,*) t_max ! Maximum time step
        read(1,*) Nwri ! Dump data per ’Nwri’ steps
        read(1,*) rho_h, rho_l ! Maximum and minimum density
        ! Relaxation time for high-density and low-density fluid,
        ! respectively
        read(1,*) tau_h, tau_l
        read(1,*) thick ! Thickness of the interface, several
        ! lattice units
        read(1,*) sigma ! Surface tension
        read(1,*) RR ! Droplet or bubble’s radii
        read(1,*) UU ! Velocity of droplet or bubble
        read(1,*) Mob ! Mobility
        read(1,*) BD ! if there is any boundary condition ****************************
        read(1,*) rh_w   ! Radius of spheres located at the bottom of domain
        close(1)
        DD = 2.d0 *RR
        LmG = rho_h - rho_l
        ! If domain size (’ix’ and ’iy’ in the params.in) is different from
        ! that in head.inc,
        ! please check params.in again and rerun the program.326 Chapter 9
        if(ix .ne. lx .or. iy .ne. ly) then
            write(*,*) "error in lx , ly input! "
            stop
        endif
        ! \beta and \kappa are main parameters in the free-energy expression
        beta = 12.d0*sigma/(thick*( (rho_h - rho_l)**4 ))
        kappa = 1.5d0*sigma*thick/( (rho_h - rho_l)**2 )
        Re1 = UU *DD*3.d0/(tau_h-0.5)
        We = UU*UU*DD*rho_h/sigma
        ! Write down the main parameters in the simulation (record in
        ! params.dat)
        open(2,file='./out/params.dat')
        write(2,'("lx=", i5, 3X, "ly=", i5)') ix, iy
        write(2,'("t_max=", i10)') t_max
        write(2,'("rho_h=", f12.5,4X, "rho_l=", f12.5)') rho_h, rho_l
        write(2,'("tau_h=", f12.5,4X, "tau_l=", f12.5)') tau_h, tau_l
        write(2,'("thickness=", f12.7)') thick
        write(2,'("beta=", f12.6)') beta
        write(2,'("Nwri=", i9)') Nwri
        write(2,'("RR=", f12.5)') RR
        write(2,'("sigma=",f12.8)') sigma
        write(2,'("kappa=",f12.8)') kappa
        write(2,'("UU=",f12.8)') UU
        write(2,'("Eo=",f12.8)') Eo
        write(2,'("Mo=",f12.8)') Mo
        write(2,'("We=",f12.8)') We
        write(2,'("Re1=",f12.8)') Re1
        write(2,'("Mob=",f12.8)') Mob
        close(2)
        !--------------
        rhoa = 0.5d0*(rho_h+rho_l)
        rhom = 0.5d0*(rho_h-rho_l)
        !--------------
        x_cen= (1+lx)/2
        y_cen= RR
        do 21 y = 1, ly
            do 31 x = 1, lx
                u_x(x,y) = 0.d0
                u_y(x,y) = 0.d0
                if(sqrt( dble(x-x_cen)**2+dble(y-y_cen)**2) .lt. RR) then
                    u_x(x,y) = UU
                endif
                rho(x,y) = rhoa + rhom&
                         *tanh( (RR - dsqrt( dble(x-x_cen)**2&
                         + dble(y-y_cen)**2 ) ) *2.d0/thick)
                rh(x,y) = (rho(x,y)-rho_l)/(rho_h-rho_l)
   31       continue
   21   continue
        do 25 y = 1, ly
            do 15 x = 1, lx
                obst(x,y) = 0
                rhdx(x,y) = 0.d0
                rhdy(x,y) = 0.d0
                brhdx(x,y) = 0.d0
                brhdy(x,y) = 0.d0
                mudx(x,y) = 0.d0
                mudy(x,y) = 0.d0
                bmudx(x,y) = 0.d0
                bmudy(x,y) = 0.d0
                lap(x,y) = 0.d0
                pb(x,y) = 0.d0
   15       continue
25      continue
        !initial values for solid boundary in Domain
         do 22 x = 1 , lx
                obst(x,1) = 1
                obst(x,ly) = 1
                rh(x,1) = rh_w
                rh(x,ly) = rh_w
22      continue
        !-------------------------
        ! Laplacian operator
        call laplace(rh,lap)
        do 40 y = 1, ly
            do 30 x = 1, lx
                mu(x,y) = 4.d0 * beta*(rh(x,y) -0.5d0)&
                        *(rh(x,y) -1.d0)*(rh(x,y)-0.d0)&
                        - kappa* lap(x,y)
   30       continue
   40   continue
        !-------------------------
        call Grad(rh, rhdx, rhdy, brhdx, brhdy, 1)
        call Grad(mu, mudx, mudy, bmudx, bmudy, 1)
        !-------------------------
        do 80 y = 1, ly
            do 70 x = 1, lx
                call get_feq(u_x(x,y), u_y(x,y), feq)
                do k=0, 8
                    ff(k,x,y) = rh(x,y)*feq(k) - 0.5d0 *feq(k)/c_squ&
                              *( (xc(k)-u_x(x,y) )*(rhdx(x,y)*c_squ - rh(x,y)*mudx(x,y))&
                              +(yc(k)-u_y(x,y) )*(rhdy(x,y)*c_squ - rh(x,y)*mudy(x,y)))
            enddo
   70      continue
80      continue
        !********************************************************************
!                        obst(x,y) = 1
!                        rh(x,y) = 0.041
!                    endif
!        ns = lx/R
!        xc_sp = R
!        yc_sp = R
!        do 12 k = 1, ns
!            do 22 y = 2, 2*R
!                do 32 x = 2, lx
!                    rpoint = nint(sqrt(real(x-xc_sp)**2+real(y-yc_sp)**2))
!                    if (rpoint .le. R) then
!                        obst(x,y) = 1
!                        rh(x,y) = 0.041
!                    endif
!32              continue
!22          continue
!            xc_sp=xc_sp + 2*R + 2
!12      continue
        end
