      program main
      implicit none
      include 'head.inc'
      integer k, x, y, time
      real*8 t_0, t_1, t_2, feq(0:8), feq2(0:8)
      integer*4 now(3)
      !---------------------------------------
      ! Author: Haibo Huang, huanghb@ustc.edu.cn
      !---------------------------------------
      ! Discrete velocity model (D2Q9),
      ! components e_ix, e_iy, where i=0,1,2,3,...8
      data xc/0.d0,1.d0,0.d0, -1.d0, 0.d0, 1.d0, -1.d0, -1.d0, 1.d0/,&
           yc/0.d0,0.d0,1.d0, 0.d0, -1.d0, 1.d0, 1.d0, -1.d0, -1.d0/,&
           ex/0, 1, 0, -1, 0, 1, -1, -1, 1/,&
           ey/0, 0, 1, 0, -1, 1, 1, -1, -1/!,&
           !opp/3, 4, 1, 2, 7, 8, 5, 6/!***********************************
      cc = 1.d0
      c_squ = cc *cc / 3.d0
      ! Weighting coefficient in the equilibrium distribution functions.
      t_0 = 4.d0 / 9.d0
      t_1 = 1.d0 / 9.d0
      t_2 = 1.d0 / 36.d0
      con2=(2.d0 * c_squ *c_squ)
      t_k(0) = t_0
      do 1 k =1,4
          t_k(k) = t_1
    1 continue
      do 2 k =5,8
          t_k(k) = t_2
    2 continue
      error= 1.d0
      ! Initialization
      call init()
      !------------------------------------------------
      open(39,file='out/residue.dat')
      do 100 time = 0, t_max
          if(mod(time,Nwri) .eq. 0) then
              ! Output the time
              write(*,*) time
              call itime(now)
              write(*,"(i2.2, ':', i2.2, ':', i2.2)") now
              ! Criterion for convergency (calculate velocity field error)
          if(time .gt. 5) call calc_error(time)
              do y= 1, ly
                  do x= 1, lx
                      up(x,y) = u_x(x,y)
                      vp(x,y) = u_y(x,y)
                  enddo
              enddo
          if(error .lt. 1.0e-10) goto 101
          endif
          !----------
          ! Collision step
          call collision()
          ! Streaming step
          call stream()
          ! Boundary Condition 
          call bounceback()
          ! Update macro-variables (density, velocity, and pressure)
          call getuv()
          !----------------------------------------
          if(mod(time ,Nwri) .eq. 0) then
              ! Output Tecplot data file
              call write_results(u_x,u_y, time)
          endif
      !--------------------------------------------
100   continue
101   close(39)
      end
