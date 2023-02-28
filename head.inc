        ! 1. This program is aimed to simulate a droplet or a bubble immersed
        ! in another fluid.
        ! 2. The fortran files: head.inc, Main.for, Macro.for,
        ! StreamCollision.for, Initial.for,
        ! Utility.for, Output.for, params.in should be put in an identical
        ! folder and compile
        ! together.
        ! 3. The main parameters control the flow are in params.in.
        ! “params.in” begin with
        ! the first line “101 101 !lx ly"
        integer lx,ly,BD
        PARAMETER(lx=201,ly=61)
        
        common /fg/ ff(0:8,lx,ly),gp(0:8,lx,ly),pb(lx,ly),u_x(lx,ly),u_y(lx,ly),t_k(0:8)
        real*8 ff, gp, pb, u_x, u_y, t_k
        
        common /rh/ rho(lx,ly),rh(lx,ly),rhp(lx,ly),up(lx,ly),vp(lx,ly)&
        ,xc(0:8),yc(0:8),ex(0:8),ey(0:8),opp(0:8)
        real*8 rho, rh, rhp, xc, yc, up, vp
        integer ex, ey, opp
        
        common/cons/ cc, c_squ, con2, rho_h, rho_l, RR, thick, kappa
        real*8 cc, c_squ, con2, rho_h, rho_l, RR, thick, kappa
        
        common/aa/ beta, tau_h, tau_l, tau, rhoa, rhom, rh_w, error
        real*8 beta, tau_h, tau_l, tau, rhoa, rhom, rh_w, error
        
        common/bb/ LmG, Eo, Mo, DD, UU, We, Mob, sigma
        real*8 LmG, Eo, Mo, DD, UU, We, Mob, sigma
        
        common/dd/ rhdx(lx,ly),rhdy(lx,ly),mudx(lx,ly),mudy(lx,ly)&
        ,rhpdx(lx,ly), rhpdy(lx,ly), brhpdx(lx,ly),brhpdy(lx,ly)
        real*8 rhdx, rhdy, mudx, mudy,rhpdx,rhpdy, brhpdx,brhpdy

        common/lap/ lap(lx,ly), mu(lx,ly)
        real*8 lap, mu
        
        common/tt/ t_max, Nwri, x_cen, y_cen
        integer t_max, Nwri, x_cen, y_cen
        
        common /ob/ R, ns, xc_sp, yc_sp, rpoint, obst(lx,ly)    
        integer R, ns, xc_sp, yc_sp, obst
        real*8 rpoint 
        
        common/mu/ bmudx(lx,ly), bmudy(lx,ly), brhdx(lx,ly)&
        ,brhdy(lx,ly),lapmu(lx,ly)
        
        real*8 bmudx, bmudy, brhdx, brhdy, lapmu

        common/boubac/ xtemp, ytemp, xtemp_rev, ytemp_rev
        integer xtemp, ytemp, xtemp_rev, ytemp_rev
