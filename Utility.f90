      !   Calculate the convergence of the numerical solution
        subroutine calc_error(time)
        implicit none
        include 'head.inc'
        integer time,i,j!,n_free
        real*8 err_a,err_b!,a_vel
        err_a = 0.d0 
        err_b = 0.d0
        do 30 j = 2, ly - 1
            do 30 i = 1, lx
                if(obst(i,j) .eq. 0 ) then
                    err_a = err_a + ( (u_x(i,j) - up(i,j))*(u_x(i,j) - up(i,j)) &
                          +(u_y(i,j) - vp(i,j))*(u_y(i,j) - vp(i,j)) )
                    err_b = err_b + (u_x(i,j)*u_x(i,j) + u_y(i,j)*u_y(i,j))
                endif
   30   continue
        error = dsqrt(err_a/err_b)
        write(*,*) error
        write(39,*) time, error
        return
        end
        !-----------------------
        subroutine get_feq(ux, uy, fequ)
        implicit none
        include 'head.inc'
        integer i
        real*8 ux,uy,fequ(0:8),u_n(0:8),u_squ,p1, tmp1
        u_squ = ux * ux + uy * uy
        tmp1 = u_squ *1.5d0
        fequ(0) = t_k(0)* ( 1.d0 - tmp1 )
        do i = 1, 8
            u_n(i) = ux *xc(i) + uy* yc(i)
            fequ(i) = t_k(i)* (1.d0+ u_n(i) / c_squ&
                    + u_n(i) * u_n(i) /con2&
                    - tmp1 )
        enddo
        end
        !---------------------------------------------------------
        ! Get gradient values for macro-variable ’p’
        subroutine Grad(p,dx,dy,bdx,bdy, id)
        implicit none
        include 'head.inc'
        integer x,y,xp, yp, xp2, yp2, xn, xn2, yn, yn2, id
        real*8 p(lx,ly), dx(lx,ly), dy(lx,ly)
        real*8 bdx(lx,ly), bdy(lx,ly)
        do 10 y = 1, ly
            do 20 x = 1, lx
                ! Periodic boundary condition for left and right boundaries
                !***********************************************************
                !For top and bottom boundaries when the points (xs + e_a dt) and (xs + 2*e_a dt) are located outside 
                !the computational domain while (xs - e_a dt) and (xs - 2*e_a dt) are not, any unknown variable / outside
                !the fluid domain is approximated by (Lee 2010)
                !(xs + e_a dt)   = (xs - e_a dt)
                !(xs + 2*e_a dt) = (xs - 2*e_a dt)
                !***********************************************************
                !The rh_w in solid nodes will normally be kept at its initial specified value. This means that the calculation
                !of rh is only applicable to the fluid nodes. However, the calculation of Sai(rh) and sai(rho) should be applied
                !to all nodes, including solid nodes. 
                !Sai(rho) = p − c_s^2*rho
                !Sai(rh) = p_th − c_s^2*rh
                xp = x + 1
                yp = y + 1
                xn = x - 1
                yn = y - 1
                xp2 = x + 2
                xn2 = x - 2
                yp2 = y + 2
                yn2 = y-2
                if(x   .eq. lx) xp  = 1
                if(x   .eq. 1 ) xn  = lx
                if(xp2 .gt. lx) xp2 = xp2 - lx
                if(xn2 .lt. 1 ) xn2 = xn2 + lx
                !*****************top and bottom boundaries*****************
                if(y .eq. 1 ) then                    
                    yn = yp
                    xn = xp
                end if
                if(yn2 .lt. 1 ) then 
                    yn2 = yp2
                    xn2 = xp2
                end if
                if(y  .eq. ly) then
                    yp = yn
                    xp = xn
                end if
                if(yp2 .gt. ly) then
                    yp2 = yn2
                    xp2 = xn2
                end if
                !----------------------------------------------------
                dx(x,y) = (p(xp,y)-p(xn,y))/3.d0 +(p(xp,yp)-p(xn,yn))/12.d0&
                        + (p(xp,yn)-p(xn,yp))/12.d0
                dy(x,y) = (p(x,yp)-p(x,yn))/3.d0 +(p(xp,yp)-p(xn,yn))/12.d0&
                        + (p(xn,yp)-p(xp,yn))/12.d0
                ! Biased difference
                if(id .eq. 1) then
                    bdx(x,y) =(-p(xp2,y)+4.d0*p(xp,y)-(-p(xn2,y) +4.d0*p(xn,y)) )/6.d0&
                             +(-p(xp2,yp2)+4.d0*p(xp,yp) -(-p(xn2,yn2)+4.d0*p(xn,yn))&
                             +(-p(xp2,yn2)+4.d0*p(xp,yn))-(-p(xn2,yp2)+4.d0*p(xn,yp)))/24.d0
                    bdy(x,y) = (-p(x,yp2)+4.d0*p(x,yp)-(-p(x,yn2) +4.d0*p(x,yn)))/6.d0&
                             +(-p(xp2,yp2)+4.d0*p(xp,yp) -(-p(xn2,yn2)+4.d0*p(xn,yn))&
                             +(-p(xn2,yp2)+4.d0*p(xn,yp))-(-p(xp2,yn2)+4.d0*p(xp,yn)) )/24.d0
                endif
   20       continue
10      continue
        end
        !---------------------------------------
        ! Get directional derivatives at position (x,y)
        subroutine xi_difference(x,y, p, xid, xidc)
        implicit none
        include 'head.inc'
        real*8 p(lx,ly), xid(0:8), xidc(0:8)
        integer x, y, xp, yp, xn, yn, xp2, yp2, xn2, yn2
        xp = x + 1
        yp = y + 1
        xn = x - 1
        yn = y - 1
        xp2 = x + 2
        xn2 = x - 2
        yp2 = y + 2
        yn2 = y-2
        if(x   .eq. lx) xp  = 1
        if(x   .eq. 1 ) xn  = lx
        if(xp2 .gt. lx) xp2 = xp2 - lx
        if(xn2 .lt. 1 ) xn2 = xn2 + lx
        !*****************top and bottom boundaries*****************
        if(y .eq. 1 ) then                    
            yn = yp
            xn = xp
        end if
        if(yn2 .lt. 1 ) then 
            yn2 = yp2
            xn2 = xp2
        end if
        if(y  .eq. ly) then
            yp = yn
            xp = xn
        end if
        if(yp2 .gt. ly) then
            yp2 = yn2
            xp2 = xn2
        end if
        
        !----------------------
        ! Half of biased and central finite differences
        xid(0) = 0.d0
        xid(1) =0.25d0* (-p(xp2,y)+5.d0*p(xp,y) -3.d0*p(x,y) -p(xn,y))
        xid(2) =0.25d0* (-p(x,yp2)+5.d0*p(x,yp) -3.d0*p(x,y) -p(x,yn)) 
        xid(3) =0.25d0* (-p(xn2,y)+5.d0*p(xn,y) -3.d0*p(x,y) -p(xp,y))
        xid(4) =0.25d0* (-p(x,yn2)+5.d0*p(x,yn) -3.d0*p(x,y) -p(x,yp))
        xid(5) =0.25d0*(-p(xp2,yp2)+5.d0*p(xp,yp)-3.d0*p(x,y) -p(xn,yn))
        xid(6) =0.25d0*(-p(xn2,yp2)+5.d0*p(xn,yp)-3.d0*p(x,y) -p(xp,yn))
        xid(7) =0.25d0*(-p(xn2,yn2)+5.d0*p(xn,yn)-3.d0*p(x,y) -p(xp,yp))
        xid(8) =0.25d0*(-p(xp2,yn2)+5.d0*p(xp,yn)-3.d0*p(x,y) -p(xn,yp))
        ! Central difference
        xidc(0) =0.d0
        xidc(1) =0.5d0* (p(xp,y) -p(xn,y))
        xidc(2) =0.5d0* (p(x,yp) -p(x,yn))
        xidc(3) =0.5d0* (p(xn,y) -p(xp,y))
        xidc(4) =0.5d0* (p(x,yn) -p(x,yp))
        xidc(5) =0.5d0*(p(xp,yp) -p(xn,yn))
        xidc(6) =0.5d0*(p(xn,yp) -p(xp,yn))
        xidc(7) =0.5d0*(p(xn,yn) -p(xp,yp))
        xidc(8) =0.5d0*(p(xp,yn) -p(xn,yp))
   10   continue
        end
        !-------------------------------------------
        ! Laplacian operator
        subroutine laplace(p,nab)
        implicit none
        include 'head.inc'
        real*8 p(lx,ly)
        real*8 nab(lx,ly)
        integer x,y,xp,xn,yp,yn
        do 10 y = 1, ly
            do 20 x = 1, lx
                yp = mod(y,ly) + 1
                xp = mod(x,lx) + 1
                yn = ly - mod(ly + 1 - y, ly)
                xn = lx - mod(lx + 1 - x, lx)
                nab(x,y) = ( p(xp,yp) +p(xn,yp) +p(xp,yn) +p(xn,yn)&
                         + 4.d0* (p(xp,y)+ p(xn,y) +p(x,yp) +p(x,yn) )&
                         - 20.d0* p(x,y) )/6.d0
   20       continue
   10   continue
        return
        end
