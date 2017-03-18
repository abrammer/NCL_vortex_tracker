

C   background_flow
C       -  Take box average at each point.
C       -  Doesn't account for cyclic or poles, just leaves border
C       - Last Edited 2014/8/14   - By A. Brammer
C   NCLFORTSTART
        subroutine background_flow(indata, bkgrd, nx,ny, dx, dy)
        real indata(nx,ny), bkgrd(nx,ny)
        integer nx, ny
        integer rx, ry
        real dx, dy
C NCLEND
        print *," * * * * * * * * * * * * "
        print *,"* * * * * * * * * * * * "
        print *," * * * * * * * * * * * * "
        print *,"  Box Average Routine "
        print *," * * * * * * * * * * * * "
        print *," Does not work near the boundaries"
        print *," Input Parameters: "
        print *," nlon =", nx ," nlat =",ny
        print *," width =",dx, "height =",dy

        rx = (dx*nx)/360
        print *, rx
        ry = (dy*(ny-1))/180
        print *, dx," = ",rx,"grid points"

        do i=(rx+1), (nx-rx)
            do j=ry, (ny-ry)
               bkgrd(i,j) = sum( indata( (i-rx):(i+rx),(j-ry):(j+ry) ) )
            enddo
        enddo

        bkgrd = bkgrd /( ((2*rx)+1)*((2*ry)+1) )

        return
        end


C  circle_avg
C       - Very Simple implementation
C       - Doesn't account for lat just takes circle on a equidistant grid
C       - Last edited 2014/8/14   - A. Brammer
C   NCLFORTSTART
        subroutine circle_avg(indata, bkgrd, nx,ny,dx)
        real indata(nx,ny), bkgrd(nx,ny), tempv
        integer nx, ny,cx ,i,j,y, ii, iii, dxi
        integer rx, ry, tot_cx, divider
        real dx
C NCLEND
        dxi = int(dx) ! need an integer version of dx

        tot_cx = 0
        do y=1, dxi
        tot_cx = tot_cx + int( sqrt( dx**2 - y**2 )+0.5)
        end do
        divider = (tot_cx*4) -1

        do i=1,nx
            do j=(dxi+1),(ny-dxi-1)
            tempv = 0.0
            y=1
            cx =   int( sqrt( dx**2 - y**2 )+0.5)
            tempv = tempv + sum( indata( (i-cx):(i+cx) , (j+(y-1))  ))
                do y=2,dxi
                cx =   int( sqrt( dx**2 - y**2 )+0.5)
                    do ii=(i-cx), (i+cx)
                        iii = ii
                        if(iii.lt.1)then
                            iii = iii+nx
                        end if
                        if(iii.gt.nx)then
                            iii = iii-nx
                        end if
        tempv = tempv + indata( iii , (j+(y-1)) )
        tempv = tempv + indata( iii , (j-(y-1)) )
                    end do
                enddo
       bkgrd(i,j) = tempv / divider
        tempv = 0.0
            enddo
        enddo
        bkgrd(:,1:dxi) = 9.96921e+36
        bkgrd(:,(ny-dxi):ny) = 9.96921e+36  !This if NCLs Float Missing Val.
        return
        end



C       circle_avg_m
C           - Accounted for Latitide and loops over cyclic
C           - Expects cyclic data, without overlap point
C           - Uses haversine to calc distances, some error but less than a grid space
C           - Last edited 2013/03/22  A Brammer.
C   NCLFORTSTART
        subroutine circle_avg_m(indata, bkgrd, inlat, nx,ny,dxm)
        real indata(nx,ny), bkgrd(nx,ny), tempv, R
        real lat(ny), inlat(ny), divider, lat1, lat2, eps, fillvalue
        integer nx, ny,cx ,i,j,y, ii, iii,dx, dxi,dxmi
        integer rx, ry, tot_cx
        real dxm, clat(ny)
C NCLEND
         eps = 2.0E-07
         fillvalue = 9.96921e+36
         R=6371. ! Earth radius km.
         pi=(355./113.) ! pi Approx good to 6 decimal places

         dxmi = int(dxm) ! need an integer version of dx
c       cx defines the longitudinal grid points for each lat.

          dx  = int(  dxmi/(2*pi*(R/nx)) ) ! expects consistent spaced
          dxi = int(dx)  ! not sure why this is needed but indexing fails without it.

c      Convert latitude to radians and cos(lat)
         lat = inlat   ! work with temp variable so as to not overwrite.
         lat = lat*(pi/180.)   ! convert to radians
         clat = cos(lat)         ! Take cos of radians

        do i=1,nx  ! Loop across longitudes
            do j=(dx+2),(ny-dx-2)  ! loop lats excluding radius of circle
            tempv = 0.0
            y=0
            divider = 0.0
            do y=-dx,dx ! work way up circle
              lat1 = lat(j)  !  centre of circle
              lat2 = lat(j+y) ! vertical distance from cirlce centre.
              xr=acos(-((sin(lat1)*sin(lat2))
     +         -cos(dxmi/R))/(cos(lat1)*cos(lat2)))         !! haversine trig
              xd = int( xr/(pi/180.)/(360./nx) )
              cx = xd ! conv. to grid points.
                 do ii=(i-cx), (i+cx)  ! work across circle
                        iii = ii
                        if(iii.lt.1)then ! loop over break assumes global w/o cyclic
                            iii = iii+nx
                        end if
                        if(iii.gt.nx)then
                            iii = iii-nx
                        end if
                  if( abs(indata( iii , (j+y))-fillvalue).gt.eps )then    !! don't include missing values
                     tempv = tempv + (clat(j+y))*indata( iii , (j+y) )
                     divider = divider + (clat(j+y))
                  end if
C                 tempv = tempv + (clat(j+y))*indata( iii , (j+y) )
C                 divider = divider + (clat(j+y))
                    end do
             end do
             if(divider .lt. eps)then            !! if any non-missing calc average
                 bkgrd(i,j) = fillvalue
              else
                 bkgrd(i,j) =  tempv/divider
              end if
             tempv = 0.0      !  not sure why this is here, it's reset at the top of the loop
            enddo  ! end j loop
        enddo      ! end i loop
        bkgrd(:,1:(dxi+2)) = indata(:,1:(dxi+2))  ! Fill poles with unsmoothed data.
        bkgrd(:,(ny-dxi-2):ny) = indata(:,(ny-dxi-2):ny)
        return
        end


C   NCLFORTSTART
        subroutine circle_avg_m_point(indata, bkgrd, inlat,inlon, nx,ny,
     +    dxm, rlat, rlon, rlat1,fillvalue  )
        integer rx, ry, tot_cx, rlat1
        real indata(nx,ny), bkgrd(rlat1), tempv, R
        real lat(ny), inlat(ny), inlon(nx), divider, lat1, lat2
        integer nx, ny,cx ,i,j,y, ii, iii, dx,dxi,dxmi
        real dxm, clat(ny), rlat(rlat1), rlon(rlat1), eps
C NCLEND
          eps = 2.0E-07
c
c       Last updated: April 22nd 2013
c       Considers latitude in size/shape of circle and weigths for average.
c       Potential error in that I use simple Trig for circle distances.
c
         R=6371. ! Earth radius km.
         pi=(355./113.) ! pi Approx good to 6 decimal places

         dxmi = int(dxm) ! need an integer version of dx
c       cx defines the longitudinal grid points for each lat.

          dx  = int(  dxmi/(2*pi*(R/nx)) )
          dxi = int(dx)  ! not sure why this is needed but indexing fails without it.

c      Convert latitude to radians and cos(lat)
         lat = inlat   ! work with temp variable so as to not overwrite.
         lat = lat*(pi/180.)   ! convert to radians
         clat = cos(lat)         ! Take cos of radians

        do k=1,rlat1
            j = minloc(abs(inlat - rlat(k) ) ,1 )
            i = minloc(abs(inlon - rlon(k) ) ,1 )
            tempv = 0.0
            y=0
            divider = 0.0
            do y=-dx,dx ! work way up circle
             if( j+y .ge. ny-5 .OR. j+y .le. 5 )then ! stuff gets funky at the poles
               cycle
             end if
              lat1 = lat(j)  !  centre of circle
              lat2 = lat(j+y) ! vertical distance from cirlce centre.
              xr=acos(-((sin(lat1)*sin(lat2))
     +         -cos(dxmi/R))/(cos(lat1)*cos(lat2)))
              xd = int( xr/(pi/180.)/(360./nx) )
              cx = xd ! conv. to grid points.
                 do ii=(i-cx), (i+cx)  ! work across circle
                        iii = ii
                        if(iii.lt.1)then ! loop over break assumes global w/o cyclic
                            iii = iii+nx
                        end if
                        if(iii.gt.nx)then
                            iii = iii-nx
                        end if
                  if( abs(indata( iii , (j+y))-fillvalue).gt.eps )then    !! don't include missing values
                  tempv = tempv + (clat(j+y))*indata( iii , (j+y) )
                  divider = divider + (clat(j+y))
                        end if
                    end do
             end do
            if(divider .lt. eps)then
              bkgrd(k) = fillvalue
            else
              bkgrd(k) =  tempv/divider
            end if
        end do
        return
        end


C       circle_avg_m
C           - Accounted for Latitide and loops over cyclic
C           - Expects cyclic data, without overlap point
C           - Uses haversine to calc distances, some error but less than a grid space
C           - Last edited 2013/03/22  A Brammer.
C   NCLFORTSTART
        subroutine circle_inout_m_point(bkgrd, inlat,inlon,nx,ny,
     +    dxm, rlat, rlon )
        real indata(nx,ny), bkgrd(nx,ny), tempv, R
        real lat(ny),inlat(ny),inlon(nx), rlat,rlon, divider, lat1, lat2
        integer nx, ny,cx ,i,j,y, ii, iii,dx, dxi,dxmi
        integer rx, ry, tot_cx
        real dxm, clat(ny)
C NCLEND

         R=6371. ! Earth radius km.
         pi=(355./113.) ! pi Approx good to 6 decimal places

         dxmi = int(dxm) ! need an integer version of dx
c       cx defines the longitudinal grid points for each lat.

          dx  = int(  dxmi/(2*pi*(R/nx)) )
          dxi = int(dx)  ! not sure why this is needed but indexing fails without it.

c      Convert latitude to radians and cos(lat)
         lat = inlat   ! work with temp variable so as to not overwrite.
         lat = lat*(pi/180.)   ! convert to radians
         clat = cos(lat)         ! Take cos of radians
            j = minloc(abs(inlat - rlat ) ,1 )
            i = minloc(abs(inlon - rlon ) ,1 )
            do y=-dx,dx ! work way up circle
              lat1 = lat(j)  !  centre of circle
              lat2 = lat(j+y) ! vertical distance from cirlce centre.
              xr=acos(-((sin(lat1)*sin(lat2))
     +         -cos(dxmi/R))/(cos(lat1)*cos(lat2)))
              xd = int( xr/(pi/180.)/(360./nx) )
              cx = xd ! conv. to grid points.
                 do ii=(i-cx), (i+cx)  ! work across circle
                        iii = ii
                        if(iii.lt.1)then ! loop over break assumes global w/o cyclic
                            iii = iii+nx
                        end if
                        if(iii.gt.nx)then
                            iii = iii-nx
                        end if
                        bkgrd( iii , (j+y) ) = 1
                  end do
             end do
        return
        end


C   NCLFORTSTART
        subroutine circle_inout_m_points(bkgrd, inlat,inlon,nx,ny,
     +    dxm, rlat, rlon, nr)
        real indata(nx,ny), bkgrd(nx,ny), tempv, R
        real lat(ny),inlat(ny),inlon(nx), rlat(nr),rlon(nr), divider
        real lat1, lat2
        integer nx, ny,cx ,i,j,y, ii, iii,dx, dxi,dxmi
        integer rx, ry, tot_cx, nr,l
        real dxm, clat(ny)
C NCLEND

         R=6371. ! Earth radius km.
         pi=(355./113.) ! pi Approx good to 6 decimal places

         dxmi = int(dxm) ! need an integer version of dx
c       cx defines the longitudinal grid points for each lat.

         dx  = int(  dxmi/(2*pi*(R/nx)) )
         dxi = int(dx)  ! not sure why this is needed but indexing fails without it.

c      Convert latitude to radians and cos(lat)
         lat = inlat   ! work with temp variable so as to not overwrite.
         lat = lat*(pi/180.)   ! convert to radians
         clat = cos(lat)         ! Take cos of radians
          do l = 1, nr
            j = minloc(abs(inlat - rlat(l) ) ,1 )
            i = minloc(abs(inlon - rlon(l) ) ,1 )
            do y=-dx,dx ! work way up circle
              lat1 = lat(j)  !  centre of circle
              lat2 = lat(j+y) ! vertical distance from cirlce centre.
              xr=acos(-((sin(lat1)*sin(lat2))
     +         -cos(dxmi/R))/(cos(lat1)*cos(lat2)))
              xd = int( xr/(pi/180.)/(360./nx) )
              cx = xd ! conv. to grid points.
                 do ii=(i-cx), (i+cx)  ! work across circle
                        iii = ii
                        if(iii.lt.1)then ! loop over break assumes global w/o cyclic
                            iii = iii+nx
                        end if
                        if(iii.gt.nx)then
                            iii = iii-nx
                        end if
                        bkgrd( iii , (j+y) ) = bkgrd( iii , (j+y) )+ 1
                  end do
            end do
          end do
        return
        end


C   NCLFORTSTART
        subroutine circle_stddev_m_point(indata, bkgrd,ave, inlat,inlon,
     +    nx, ny, dxm, rlat, rlon, rlat1,fillvalue  )
        integer rx, ry, tot_cx, rlat1
        real indata(nx,ny), bkgrd(rlat1), ave(rlat1), tempv, R,fillvalue
        real lat(ny), inlat(ny), inlon(nx), divider, lat1, lat2
        integer nx, ny,cx ,i,j,y, ii, iii,dx, dxi,dxmi
        real dxm, clat(ny), rlat(rlat1), rlon(rlat1),eps
C NCLEND
c
c       Last updated: April 22nd 2013
c       Considers latitude in size/shape of circle and weigths for average.
c       Potential error in that I use simple Trig for circle distances.
c
c
c
          eps = 2.0E-07

         R=6371. ! Earth radius km.
         pi=(355./113.) ! pi Approx good to 6 decimal places

         dxmi = int(dxm) ! need an integer version of dx
c       cx defines the longitudinal grid points for each lat.

          dx  = int(  dxmi/(2*pi*(R/nx)) )
          dxi = int(dx)  ! not sure why this is needed but indexing fails without it.

c      Convert latitude to radians and cos(lat)
         lat = inlat   ! work with temp variable so as to not overwrite.
         lat = lat*(pi/180.)   ! convert to radians
         clat = cos(lat)         ! Take cos of radians

        do k=1,rlat1
         if(abs(ave(k)-fillvalue).lt.eps)then
                 bkgrd(k) = fillvalue
            cycle
         end if
         j = minloc(abs(inlat - rlat(k) ) ,1 )
         i = minloc(abs(inlon - rlon(k) ) ,1 )

            tempv = 0.0
            y=0
            divider = 0.0
            do y=-dx,dx ! work way up circle
             if( j+y .ge. ny .OR. j+y .le. 0 )then
               cycle
             end if
              lat1 = lat(j)  !  centre of circle
              lat2 = lat(j+y) ! vertical distance from cirlce centre.
              xr=acos(-((sin(lat1)*sin(lat2))
     +         -cos(dxmi/R))/(cos(lat1)*cos(lat2)))
              xd = int( xr/(pi/180.)/(360./nx) )
              cx = xd ! conv. to grid points.
                 do ii=(i-cx), (i+cx)  ! work across circle
                        iii = ii
                        if(iii.lt.1)then ! loop over break assumes global w/o cyclic
                            iii = iii+nx
                        end if
                        if(iii.gt.nx)then
                            iii = iii-nx
                        end if
                  if( abs(indata( iii , (j+y))-fillvalue).gt.eps )then    !! don't include missing values
            tempv = tempv + (clat(j+y))*((indata(iii,(j+y) )-ave(k))**2)
                  divider = divider + (clat(j+y))
                  end if
                    end do
             end do
        if(divider .lt. eps)then
        bkgrd(k) = fillvalue
        else
        bkgrd(k) =  sqrt( tempv/divider )
        end if
        end do
        return
        end



C   NCLFORTSTART
        subroutine circle_lthresh_m_point(indata,bkgrd,thresh,
     +    inlat,inlon,nx, ny, dxm, rlat, rlon, rlat1,fillvalue  )
        integer rx, ry, tot_cx, rlat1
        real indata(nx,ny), bkgrd(rlat1), ave(rlat1), tempv, R,fillvalue
        real lat(ny), inlat(ny), inlon(nx), divider, lat1, lat2
        integer nx, ny,cx ,i,j,y, ii, iii,dx, dxi,dxmi
        real dxm, clat(ny), rlat(rlat1), rlon(rlat1),eps
C NCLEND
c
c       Last updated: April 22nd 2013
c       Considers latitude in size/shape of circle and weigths for average.
c       Potential error in that I use simple Trig for circle distances.
c
c
c
          eps = 2.0E-07

         R=6371. ! Earth radius km.
         pi=(355./113.) ! pi Approx good to 6 decimal places

         dxmi = int(dxm) ! need an integer version of dx
c       cx defines the longitudinal grid points for each lat.

          dx  = int(  dxmi/(2*pi*(R/nx)) )
          dxi = int(dx)  ! not sure why this is needed but indexing fails without it.

c      Convert latitude to radians and cos(lat)
         lat = inlat   ! work with temp variable so as to not overwrite.
         lat = lat*(pi/180.)   ! convert to radians
         clat = cos(lat)         ! Take cos of radians

        do k=1,rlat1
         j = minloc(abs(inlat - rlat(k) ) ,1 )
         i = minloc(abs(inlon - rlon(k) ) ,1 )

            tempv = 0.0
            y=0
            divider = 0.0
            do y=-dx,dx ! work way up circle
             if( j+y .ge. ny .OR. j+y .le. 0 )then
               cycle
             end if
              lat1 = lat(j)  !  centre of circle
              lat2 = lat(j+y) ! vertical distance from cirlce centre.
              xr=acos(-((sin(lat1)*sin(lat2))
     +         -cos(dxmi/R))/(cos(lat1)*cos(lat2)))
              xd = int( xr/(pi/180.)/(360./nx) )
              cx = xd ! conv. to grid points.
                 do ii=(i-cx), (i+cx)  ! work across circle
                        iii = ii
                        if(iii.lt.1)then ! loop over break assumes global w/o cyclic
                            iii = iii+nx
                        end if
                        if(iii.gt.nx)then
                            iii = iii-nx
                        end if
                  if( abs(indata( iii , (j+y))-fillvalue).gt.eps )then  !! don't include missing values
                    if(indata( iii , (j+y)).lt. thresh)then
                      tempv = tempv + (clat(j+y))
                    end if
                  end if
                 divider = divider + (clat(j+y))
                end do
             end do
        if(divider .lt. eps)then
        bkgrd(k) = fillvalue
        else
        bkgrd(k) = tempv/divider
        end if
        end do
        return
        end

C   NCLFORTSTART
        subroutine circle_gthresh_m_point(indata,bkgrd,thresh,
     +    inlat,inlon,nx, ny, dxm, rlat, rlon, rlat1,fillvalue  )
        integer rx, ry, tot_cx, rlat1
        real indata(nx,ny), bkgrd(rlat1), ave(rlat1), tempv, R,fillvalue
        real lat(ny), inlat(ny), inlon(nx), divider, lat1, lat2
        integer nx, ny,cx ,i,j,y, ii, iii,dx, dxi,dxmi
        real dxm, clat(ny), rlat(rlat1), rlon(rlat1),eps
C NCLEND
c
c       Last updated: April 22nd 2013
c       Considers latitude in size/shape of circle and weigths for average.
c       Potential error in that I use simple Trig for circle distances.
c
c
c
          eps = 2.0E-07

         R=6371. ! Earth radius km.
         pi=(355./113.) ! pi Approx good to 6 decimal places

         dxmi = int(dxm) ! need an integer version of dx
c       cx defines the longitudinal grid points for each lat.

          dx  = int(  dxmi/(2*pi*(R/nx)) )
          dxi = int(dx)  ! not sure why this is needed but indexing fails without it.

c      Convert latitude to radians and cos(lat)
         lat = inlat   ! work with temp variable so as to not overwrite.
         lat = lat*(pi/180.)   ! convert to radians
         clat = cos(lat)         ! Take cos of radians

        do k=1,rlat1
         j = minloc(abs(inlat - rlat(k) ) ,1 )
         i = minloc(abs(inlon - rlon(k) ) ,1 )

            tempv = 0.0
            y=0
            divider = 0.0
            do y=-dx,dx ! work way up circle
             if( j+y .ge. ny .OR. j+y .le. 0 )then
               cycle
             end if

              lat1 = lat(j)  !  centre of circle
              lat2 = lat(j+y) ! vertical distance from cirlce centre.
              xr=acos(-((sin(lat1)*sin(lat2))
     +         -cos(dxmi/R))/(cos(lat1)*cos(lat2)))
              xd = int( xr/(pi/180.)/(360./nx) )
              cx = xd ! conv. to grid points.
                 do ii=(i-cx), (i+cx)  ! work across circle
                        iii = ii
                        if(iii.lt.1)then ! loop over break assumes global w/o cyclic
                            iii = iii+nx
                        end if
                        if(iii.gt.nx)then
                            iii = iii-nx
                        end if
                  if( abs(indata( iii , (j+y))-fillvalue).gt.eps )then  !! don't include missing values
                    if(indata( iii , (j+y)) .gt. thresh)then
                      tempv = tempv + (clat(j+y))
                    end if
                  end if
                 divider = divider + (clat(j+y))
                end do
             end do
        if(divider .lt. eps)then
        bkgrd(k) = fillvalue
        else
        bkgrd(k) = tempv/divider
        end if
        end do
        return
        end




C   NCLFORTSTART
        subroutine circle_missing_m_point(indata,bkgrd,
     +    inlat,inlon,nx, ny, dxm, rlat, rlon, rlat1,fillvalue  )
        integer rx, ry, tot_cx, rlat1
        integer nx, ny,cx ,i,j,y, ii, iii,dx, dxi,dxmi
        real indata(nx,ny), bkgrd(rlat1), ave(rlat1), tempv, R,fillvalue
        real lat(ny), inlat(ny), inlon(nx), divider, lat1, lat2
        real dxm, clat(ny), rlat(rlat1), rlon(rlat1),eps
C NCLEND
c
c       Last updated: April 22nd 2013
c       Considers latitude in size/shape of circle and weigths for average.
c       Potential error in that I use simple Trig for circle distances.
c
c
        eps = 2.0E-07

         R=6371. ! Earth radius km.
         pi=(355./113.) ! pi Approx good to 6 decimal places

         dxmi = int(dxm) ! need an integer version of dx
c       cx defines the longitudinal grid points for each lat.

          dx  = int(  dxmi/(2*pi*(R/nx)) )
          dxi = int(dx)  ! not sure why this is needed but indexing fails without it.

c      Convert latitude to radians and cos(lat)
         lat = inlat   ! work with temp variable so as to not overwrite.
         lat = lat*(pi/180.)   ! convert to radians
         clat = cos(lat)         ! Take cos of radians

        do k=1,rlat1
         j = minloc(abs(inlat - rlat(k) ) ,1 )
         i = minloc(abs(inlon - rlon(k) ) ,1 )

            tempv = 0.0
            y=0
            divider = 0.0
            do y=-dx,dx ! work way up circle
             if( j+y .ge. ny .OR. j+y .le. 0 )then
               cycle
             end if
              lat1 = lat(j)  !  centre of circle
              lat2 = lat(j+y) ! vertical distance from cirlce centre.
              xr=acos(-((sin(lat1)*sin(lat2))
     +         -cos(dxmi/R))/(cos(lat1)*cos(lat2)))
              xd = int( xr/(pi/180.)/(360./nx) )
              cx = xd ! conv. to grid points.
                 do ii=(i-cx), (i+cx)  ! work across circle
                        iii = ii
                        if(iii.lt.1)then ! loop over break assumes global w/o cyclic
                            iii = iii+nx
                        end if
                        if(iii.gt.nx)then
                            iii = iii-nx
                        end if
                  if( abs(indata( iii , (j+y))-fillvalue).lt.eps )then  !! don't include missing values
                      tempv = tempv + (clat(j+y))
                  end if
                 divider = divider + (clat(j+y))
                end do
             end do
        bkgrd(k) = tempv/divider
        end do
        return
        end


