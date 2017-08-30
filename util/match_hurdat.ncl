undef("match_hurdat")
procedure match_hurdat(rtracks, stormNames,type)
local years, genlocs, lat, lon, gtimes, sdata, slat, g, tr
begin
    ut_time = cd_calendar(rtracks&time ,0)
    years = toint(unique(ut_time(:,0)))

    genlocs = HURDAT_gen_info(years)
    if(all(genlocs.eq.False)) then return end if
    lat := string_latlon_to_numeric(genlocs(3,:))
    lon := string_latlon_to_numeric(genlocs(4,:))

    gtimes = cd_convert( genlocs&time, rtracks&time@units)

    do g=0, dimsizes(gtimes)-1
        sdata := HURDAT_storm_info(years, genlocs(0,g) )
        slat := string_latlon_to_numeric(sdata(:,4))
        slon := string_latlon_to_numeric(sdata(:,5))
        stime := sdata&time
        si = min( intersect_ind(sdata&time, rtracks&time) )
        if(ismissing(si)) continue end if
        ni = min( (/ max(intersect_ind(sdata&time,rtracks&time)),si+10 /) )
        print(si+"   "+ni)
        if(si .eq. ni) continue end if
        glats := conform( rtracks(0,:,{stime(si:ni)}), slat(si:ni) ,1)
        glons := conform( rtracks(0,:,{stime(si:ni)}), slon(si:ni) ,1)
        dists :=  gc_latlon( glats, glons,  rtracks(0,:,{stime(si:ni)}),  rtracks(1,:,{stime(si:ni)}), 0, 4)
        tr = minind(dim_avg_n(dists,1))
        if( avg(dists(tr,:)).lt. 500 ) 
            print("Match "+genlocs(0,g)+","+genlocs(1,g))
            ci := intersect_ind(sdata&time, rtracks&time)
            rtracks(0,tr,{sdata&time(ci)}) = (/ slat(ci) /)
            rtracks(1,tr,{sdata&time(ci)}) = (/ slon(ci) /)
            type(tr, {sdata&time(ci)} ) = (/sdata(ci,3)/)
            rtracks(:,tr, {max(sdata&time(ci)):} ) = rtracks@_FillValue
            stormNames(tr) = genlocs(0,g)
        end if
    end do
end
