        ! Write binary Tecplot file
        subroutine write_results(upx,upy, n)
        implicit none
        include "head.inc"
        integer x,y,i,j,n, k
        real*8 upx(lx,ly),upy(lx,ly)
        REAL*4 ZONEMARKER, EOHMARKER
        integer*2 ii, jj
        integer*1 obob
        integer, parameter:: kmax=1
        character*40 Title,var
        character*40 V1,V2,V3,V4, V5, V6, V7, V8
        character*40 Zonename1
        character filename*46, B*7, C*4
        write(B,'(i7.7)') n
        filename='out/_'//B//'.plt'
        open(41,file=filename, form='binary')
        !--------------------------------------
        ZONEMARKER= 299.0
        EOHMARKER = 357.0
        !--------------------------------------
        write(41) "#!TDV101"
        write(41) 1
        Title="LeeCF10"
        call dumpstring(Title)
        ! Number of variables (NumVar) in the datafile.
        write(41) 7
        V1='X'
        call dumpstring(V1)
        V2='Y'
        call dumpstring(V2)
        V3='u'
        call dumpstring(V3)
        V4='v'
        call dumpstring(V4)
        V5='rho'
        call dumpstring(V5)
        V7='obst'
        call dumpstring(V7)
        V8='pb'
        call dumpstring(V8)
        ! Zones---------------------
        ! Zone marker. Value = 299.0
        write(41) ZONEMARKER
        Zonename1="ZONE 001"
        call dumpstring(Zonename1)
        ! Zone Color
        write(41) -1
        ! ZoneType
        write(41) 0
        ! DataPacking 0=Block, 1=Point
        write(41) 1
        ! Specify Var Location.
        ! Don’t specify=0, all data is located at the nodes. Specify=1
        write(41) 0
        ! Number of user defined face neighbor connections (value >= 0)
        write(41) 0
        ! IMax,JMax,KMax
        write(41) lx
        write(41) ly
        write(41) kmax 
        ! Auxiliary name/value pair to follow=1 No more Auxiliary
        ! name/value pairs=1.
        write(41) 0
        ! HEADER OVER-------
        ! EOHMARKER, value=357.0
        write(41) EOHMARKER
        ! II. Data section+++++++++++++++++++++++++
        write(41) Zonemarker
        ! Data format, 1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte,
        ! 6=Bit
        write(41) 4
        write(41) 4
        write(41) 1
        write(41) 1
        write(41) 1
        write(41) 5
        write(41) 1
        ! Has variable sharing 0 = no, 1 = yes.
        write(41) 0
        ! Zone number to share connectivity list with (-1 = no sharing).
        write(41) -1
        ! Zone Data. Each variable is in data format as specified above.
        do k=1,kmax
            do j=1,ly
                do i=1,lx
                    if(obst(i,j).ne. 0) then
                    else
                    endif
                    ii = i
                    jj = j
                    obob = obst(i,j)
                    write(41) ii
                    write(41) jj
                    write(41) real(upx(i,j))
                    write(41) real(upy(i,j))
                    write(41) real(rho(i,j))
                    write(41) obob
                    write(41) real(pb(i,j))
                end do
            end do
        end do
        close(41)
        end
        !--------------------------------
        subroutine dumpstring(instring)
        character(40) instring
        integer len
        len=LEN_TRIM(instring)
        do ii=1,len
          I=ICHAR(instring(ii:ii))
          write(41) I
        end do
        write(41) 0
        return
        end
