undef("make_lon_monotonic")
function make_lon_monotonic(lon)
local lon, dl, nl,isM,disc,diff
begin

dl = lon(1) - lon(0)
nl = dimsizes(lon)

isM = isMonotonic(lon)

if(isM.ne.0)
    return(lon)
end if

disc :=  ind( abs(lon(1:) - lon(:nl-2)) .gt. 2*dl ) 

if(ismissing(disc))
;    print("make_lon_monotonic Error: No discontiuity found.")
;    print("Returned without adjustment")
    return(lon)
end if
if(dimsizes(disc).gt.1)
;    print("Can't Handle multiple discontinuities yet")
;    print("Returned without adjustment")
;    return(lon)
    disc := disc(0)
end if

diff = lon(disc+1) - lon(disc)
if(diff.lt.0) 
    lon(disc+1:) = lon(disc+1:)+360
else
   lon(disc+1:) = lon(disc+1:)-360
end if

if(all(lon .gt. 360 ) )
    lon = lon - 360
end if
if(all(lon .lt. -360 ) )
    lon = lon + 360
end if

isM = isMonotonic(lon)
if(isM.eq.0)
  print("unravelling multiple discontinuities")
   lon = make_lon_monotonic(lon)
end if
        
return(lon)
end



undef("make_lon_monotonic_dummy")
procedure make_lon_monotonic_dummy()
local lon
begin

lon = lonGlobeF(360, "lon", "longitude", "degrees_east")
lon = lonPivot(lon, 73)  ;;  Lets pivot at 73. 
lon&lon = where(lon&lon.gt.360, lon&lon-360, lon&lon) ;; Wrap values so they are within <360 
printVarSummary(lon)
;;;;   Recreation of lon values from chi.nc

print(isMonotonic(lon&lon))  
print("Adjust lons so that they are monotonic")
lon&lon = make_lon_monotonic(lon&lon)  ;;  Adjust lon so that it is monotonic. order not chagned, just add/subtract 360

printVarSummary(lon)
print(isMonotonic(lon&lon))  ;; Success hopefully. 
printMinMax(lon&lon,0)

end


