        subroutine getuv( )
        implicit none
        include "head.inc"
        integer x,y
        do 10 y = 1, ly
            do 10 x = 1, lx
                if(obst(x,y) .eq. 0) then
                    rh(x,y) = ff(0,x,y) +ff(1,x,y) +ff(2,x,y)&
                            + ff(3,x,y) +ff(4,x,y) +ff(5,x,y)&
                            + ff(6,x,y) +ff(7,x,y) +ff(8,x,y)&
                            + 0.5d0* Mob *lapmu(x,y)
                    rho(x,y) =rho_h* rh(x,y) + rho_l*(1.d0- rh(x,y) )
                    u_x(x,y) = (gp(1,x,y) + gp(5,x,y) + gp(8,x,y)&
                             -(gp(3,x,y) + gp(6,x,y) + gp(7,x,y))) /rho(x,y)/c_squ
                    u_y(x,y) = (gp(2,x,y) + gp(5,x,y) + gp(6,x,y)&
                             -(gp(4,x,y) + gp(7,x,y) + gp(8,x,y))) /rho(x,y)/c_squ
                !else
                !    rh(x,y) = 0.041
                !    u_x(x,y) = 0.0
                !    u_y(x,y) = 0.0
                endif
   10   continue
        !-------------------------
        ! Get density gradients (both central difference and biased
        ! difference)
        call Grad(rh, rhdx, rhdy, brhdx, brhdy, 1)
        ! Get Laplacian operator for ‘rh’
        call laplace(rh,lap)
        do 40 y = 2, ly-1
            do 30 x = 1, lx
                mu(x,y) = 4.d0 * beta*(rh(x,y) -0.5d0)&
                        *(rh(x,y) -1.d0)*(rh(x,y)-0.d0)&
                        - kappa* lap(x,y)
   30       continue
   40   continue
        ! Laplacian operator
        call laplace(mu, lapmu)
        ! Get gradients for \mu (both central difference and biased
        ! difference)
        call Grad(mu, mudx, mudy, bmudx, bmudy, 1)
        !---------------
        do 80 y = 2, ly-1
            do 70 x = 1, lx
                ! Velocity components in x and y directions
                u_x(x,y) = u_x(x,y) -0.5d0*rh(x,y)*mudx(x,y)/rho(x,y)
                u_y(x,y) = u_y(x,y) -0.5d0*rh(x,y)*mudy(x,y)/rho(x,y)
                ! Pressure
                pb(x,y) = (gp(0,x,y) +gp(1,x,y) +gp(2,x,y)&
                        +gp(3,x,y) +gp(4,x,y) +gp(5,x,y)&
                        +gp(6,x,y) +gp(7,x,y) +gp(8,x,y))&
                        + 0.5d0 *( u_x(x,y)* rhdx(x,y)&
                        + u_y(x,y)* rhdy(x,y) )*c_squ
                rhp(x,y) =rh(x,y)*c_squ -pb(x,y)
   70       continue
   80   continue
        ! Get gradients for \rhp
        call Grad(rhp, rhpdx, rhpdy, brhpdx, brhpdy, 1)
        end
