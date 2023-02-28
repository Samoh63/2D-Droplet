        subroutine stream()
        implicit none
        include "head.inc"
        integer k
        real*8 f_hlp(0:8,lx,ly)
        integer x,y,x_e,x_w,y_n,y_s,l,n
        ! Streaming step
        do 10 y = 2, ly-1!**********************
            do 10 x = 1, lx
                y_n = y + 1!************************
                x_e = mod(x,lx) + 1
                y_s = y - 1!************************
                x_w = lx - mod(lx + 1 - x, lx)
                f_hlp(1,x_e,y ) = ff(1,x,y)
                f_hlp(2,x ,y_n) = ff(2,x,y)
                f_hlp(3,x_w,y ) = ff(3,x,y)
                f_hlp(4,x ,y_s) = ff(4,x,y)
                f_hlp(5,x_e,y_n) = ff(5,x,y)
                f_hlp(6,x_w,y_n) = ff(6,x,y)
                f_hlp(7,x_w,y_s) = ff(7,x,y)
                f_hlp(8,x_e,y_s) = ff(8,x,y)
   10       continue
        do 19 y = 1, ly
                do 20 x = 1, lx
                        do 21 k = 1, 8
                            ff(k,x,y) = f_hlp(k,x,y)
   21                   continue
   20           continue
   19   continue
        do 30 y = 2, ly - 1!*******
            do 30 x = 1, lx
                y_n = y + 1!************
                x_e = mod(x,lx) + 1
                y_s = y - 1!************
                x_w = lx - mod(lx + 1 - x, lx)
                f_hlp(1,x_e,y ) = gp(1,x,y)
                f_hlp(2,x ,y_n) = gp(2,x,y)
                f_hlp(3,x_w,y ) = gp(3,x,y)
                f_hlp(4,x ,y_s) = gp(4,x,y)
                f_hlp(5,x_e,y_n) = gp(5,x,y)
                f_hlp(6,x_w,y_n) = gp(6,x,y)
                f_hlp(7,x_w,y_s) = gp(7,x,y)
                f_hlp(8,x_e,y_s) = gp(8,x,y)
   30       continue
        do 39 y = 2, ly-1
            do 40 x = 1, lx
                do 31 k = 1, 8
                    gp(k,x,y) = f_hlp(k,x,y)
   31           continue
   40       continue
   39   continue
    end
        !-----------------------------------------------
        ! Collision steps
        !The bounce-back boundary condition is applied instead of the collision step on all lattice nodes 
        !that are solid walls.
        subroutine collision()
        implicit none
        include 'head.inc'
        integer x, y , k
        real*8 feq2(0:8), geq2(0:8), feq(0:8), ximud(0:8)
        real*8 xirhod(0:8), xirhodc(0:8), ximudc(0:8)
        real*8 xirhpd(0:8), xirhpdc(0:8)
        real*8 tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
        do 10 y= 2, ly-1 !*************************
            do 20 x =1, lx
                tau = tau_l + rh(x,y)*(tau_h - tau_l)
                call get_feq(u_x(x,y), u_y(x,y), feq)
                call xi_difference(x,y,rhp, xirhpd, xirhpdc)
                call xi_difference(x,y, mu, ximud, ximudc)
                call xi_difference(x,y, rho, xirhod, xirhodc)
                tmp1 = (rhdx(x,y)+brhdx(x,y))*0.5d0*LmG
                tmp2 = (rhdy(x,y)+brhdy(x,y))*0.5d0*LmG
                tmp3 = (mudx(x,y)+bmudx(x,y))*0.5d0*rh(x,y)
                tmp4 = (mudy(x,y)+bmudy(x,y))*0.5d0*rh(x,y)
                tmp5 = (rhpdx(x,y)+brhpdx(x,y))*0.5d0
                tmp6 = (rhpdy(x,y)+brhpdy(x,y))*0.5d0
                do 30 k=0, 8
                    feq2(k) = rh(x,y)*feq(k) - 0.5d0 *feq(k)/c_squ&
                            *( xirhpdc(k) -rh(x,y)*ximudc(k)&
                            +(-u_x(x,y) )*(rhpdx(x,y) - rh(x,y)*mudx(x,y))&
                            +(-u_y(x,y) )*(rhpdy(x,y) - rh(x,y)*mudy(x,y))&
                            + Mob*lapmu(x,y)*c_squ)
                    geq2(k) = pb(x,y)*t_k(k) + rho(x,y)*c_squ*(feq(k)-t_k(k))&
                            - 0.5d0*c_squ*(feq(k)-t_k(k))&
                            *(xirhodc(k)&
                            -u_x(x,y)*rhdx(x,y)*LmG&
                            -u_y(x,y) *rhdy(x,y)*LmG)&
                            - 0.5d0 *feq(k)&
                            *( -rh(x,y)*ximudc(k)&
                            -u_x(x,y)*(-rh(x,y)*mudx(x,y))&
                            -u_y(x,y)*(-rh(x,y)*mudy(x,y)))
                    !--------------------
                    ff(k,x,y) = ff(k,x,y)-(ff(k,x,y)- feq2(k))/tau&
                              + feq(k)/c_squ *(xirhpd(k)-rh(x,y)*ximud(k)&
                              +(-u_x(x,y))*(tmp5-tmp3)&
                              +(-u_y(x,y))*(tmp6-tmp4)&
                              + Mob*lapmu(x,y)*c_squ)
                    gp(k,x,y) = gp(k,x,y)-(gp(k,x,y)-geq2(k))/tau&
                              + c_squ*(feq(k)-t_k(k))&
                              *(xirhod(k)&
                              -u_x(x,y)*tmp1&
                              -u_y(x,y)*tmp2)&
                              + feq(k)&
                              *( -rh(x,y)*ximud(k)&
                              -u_x(x,y)*(-tmp3)&
                              -u_y(x,y)*(-tmp4))
   30           continue
   20       continue
   10   continue
        end
        !********************************************************
!        !Bounce Back Steps
        subroutine bounceback()
        implicit none
        include "head.inc"
        integer x, y, k
        do 21 x = 1, lx
            ff(2,x,1 ) = ff(4,x,1 )
            ff(5,x,1 ) = ff(7,x,1 )
            ff(6,x,1 ) = ff(8,x,1 )
            gp(2,x,1 ) = gp(4,x,1 )
            gp(5,x,1 ) = gp(7,x,1 ) 
            gp(6,x,1 ) = gp(8,x,1 )
            ff(4,x,ly) = ff(2,x,ly)
            ff(7,x,ly) = ff(5,x,ly)
            ff(8,x,ly) = ff(6,x,ly)
            gp(4,x,ly) = gp(2,x,ly)
            gp(7,x,ly) = gp(5,x,ly) 
            gp(8,x,ly) = gp(6,x,ly)
21      continue
        end
