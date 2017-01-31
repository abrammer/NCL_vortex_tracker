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
;    if(.not.any(rtracks&time .eq. gtimes(g) ))
;        continue ; if time is out of bounds move on.
;    end if
    sdata := HURDAT_storm_info(years, genlocs(0,g) )
    slat := string_latlon_to_numeric(sdata(:,4))
    slon := string_latlon_to_numeric(sdata(:,5))
    stime := sdata&time
    if(max(stime).lt.min(rtracks&time) .or. min(stime).gt.max(rtracks&time))
      continue ; if time is out of bounds move on.
    end if
    si = min( intersect_ind(sdata&time, rtracks&time) )
    ni = min( (/si+8, dimsizes(slat)-1, max(intersect_ind(sdata&time, rtracks&time))  /) )
    do tr=0, dimsizes(rtracks(0,:,0))-1
      if( .not.all(ismissing( rtracks(0,tr,{stime(si:ni)}) ) )  )
        dists := gc_latlon( slat(si:ni), slon(si:ni),  rtracks(0,tr,{stime(si:ni)}),  rtracks(1,tr,{stime(si:ni)}), 10, 4)
        if( any( dists .le. 250  ))
          print("Match "+tr+"   "+genlocs(1,g)+" "+min(dists) )
          ci := intersect_ind(sdata&time, rtracks&time)
          rtracks(0,tr,{sdata&time(ci)}) = (/ slat(ci) /)
          rtracks(1,tr,{sdata&time(ci)}) = (/ slon(ci) /)
          type(tr, {sdata&time(ci)} ) = (/sdata(ci,3)/)
          rtracks(:,tr, {max(sdata&time(ci)):} ) = rtracks@_FillValue
          stormNames(tr) = genlocs(1,g)
        end if
      end if
    end do
end do
end
