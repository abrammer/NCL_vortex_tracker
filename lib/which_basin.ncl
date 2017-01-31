
undef("which_basin")
function which_basin(lat, lon)
local basin_lats, basin_lons, lat, lon
begin
retval = new( dimsizes(lat), "string")
retval = "un"
;; epac
basin_lats :=  (/0,0,60,60, 20, 15,12,10,8.3,9.25,0/)
basin_lons := (/-72,-140,-140,-100,-100, -90,-85,-84,-81.5,-79,-72/)
retval = where(gc_inout( lat, lon, basin_lats, basin_lons), "ep",retval)
;; atl
basin_lons :=  (/-72,0,0,-100,-100, -90,-85,-84,-81.5,-79,-72/)
basin_lats :=  (/0,0,60,60, 20, 15,12,10,8.3,9.25,0/)
retval = where(gc_inout( lat, lon, basin_lats, basin_lons), "al",retval)

;; cpac
basin_lons := (/180, 180, 220, 220, 180/)
basin_lats := (/0,60,60,0,0/)
retval = where(gc_inout( lat, lon, basin_lats, basin_lons), "cp",retval)

;; wpac
basin_lons := (/180, 180, 100, 100, 180/)
basin_lats := (/0,60,60,0,0/)
retval = where(gc_inout( lat, lon, basin_lats, basin_lons), "wp",retval)

;; io
basin_lons := (/40, 40, 100, 100, 40/)
basin_lats := (/0,60,60,0,0/)
retval = where(gc_inout( lat, lon, basin_lats, basin_lons), "io",retval)

;; io ?
basin_lons := (/20, 20, 120, 120, 20/)
basin_lats := (/0,-60,-60,0,0/)
retval = where(gc_inout( lat, lon, basin_lats, basin_lons), "io",retval)

;if( gc_inout( lat, lon, basin_lats, basin_lons))
;return("ep")
;end if
;basin_lons :=  (/-72,0,0,-100,-100, -90,-85,-84,-81.5,-79,-72/)
;basin_lats :=  (/0,0,60,60, 20, 15,12,10,8.3,9.25,0/)
;if( gc_inout( lat, lon, basin_lats, basin_lons))
;return("al")
;end if
;basin_lons := (/180, 180, 220, 220, 180/)
;basin_lats := (/0,60,60,0,0/)
;if( gc_inout( lat, lon, basin_lats, basin_lons))
;return("cp")
;end if
;basin_lons := (/180, 180, 100, 100, 180/)
;basin_lats := (/0,60,60,0,0/)
;if( gc_inout( lat, lon, basin_lats, basin_lons))
;return("wp")
;end if
;basin_lons := (/40, 40, 100, 100, 40/)
;basin_lats := (/0,60,60,0,0/)
;if( gc_inout( lat, lon, basin_lats, basin_lons))
;return("io")
;end if
;basin_lons := (/20, 20, 120, 120, 20/)
;basin_lats := (/0,-60,-60,0,0/)
;if( gc_inout( lat, lon, basin_lats, basin_lons))
;return("io")
;end if
return(retval)
end 


