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
    do tr=0, dimsizes(rtracks(0,:,0))-1
      if(.not.any(rtracks&time .eq. gtimes(g) ))
          continue ; if time is out of bounds move on.
      end if
      if( .not.ismissing( rtracks(0,tr,{gtimes(g)}) ) )       
          sdata := HURDAT_storm_info(years, genlocs(0,g) )
          slat := string_latlon_to_numeric(sdata(:,4))
          slon := string_latlon_to_numeric(sdata(:,5))
          ni = min( (/8, dimsizes(slat)-1/) )
          stime := sdata&time
          si = min( intersect_ind(sdata&time, rtracks&time) )
        if( any( gc_latlon( slat(si:ni), slon(si:ni),  rtracks(0,tr,{stime(si:ni)}),  rtracks(1,tr,{stime(si:ni)}), 10, 4) .le. 500  ))
          print("Match "+genlocs(1,g))
          ;print( (/ gc_latlon( lat(g), lon(g),  tracks(0, {gtimes(g)} ),  tracks(1, {gtimes(g)} ), 10, 4) /) )
;          sdata := HURDAT_storm_info(years, genlocs(0,g) )
;          slat := string_latlon_to_numeric(sdata(:,4))
;          slon := string_latlon_to_numeric(sdata(:,5))
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
