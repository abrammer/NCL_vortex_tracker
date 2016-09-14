load "/home11/grad/2010/abrammer/Desktop/research/noaa_postdoc/lib/make_lon_monotonic.ncl"
load "/home11/grad/2010/abrammer/Desktop/research/noaa_postdoc/lib/check_closed_contour.ncl"


undef("unique_locations")
function unique_locations(data,rkm)
local dims, n, rpi, data, rkm, nmi, retdata
begin

dims = dimsizes(data)
do n=0, dims(0)-2
  if(ismissing( data(n,0))) continue end if
   rpi := ind( gc_latlon( data(n,0), data(n,1), data(n+1:,0), data(n+1:,1),0,4) .lt. rkm )
   if(.not.all(ismissing(rpi)))
     rpi = rpi+n+1
     data(n,:) = dim_avg_n(data(array_append_record(n, rpi,0),:),0)
     data(rpi,:) = data@_FillValue
   end if
end do
    return data
end

undef("unique_locations_compare")
function unique_locations_compare(data1, data,rkm)
local dims, n, rpi, data, rkm, nmi, retdata, dataout, dataout_
begin
if(dimsizes(dimsizes(data)).eq.1)       ;;;; make sure 2 subscripts works below.
  dataout = new( (/1, dimsizes(data)/), float)
  dataout(0,:) = data
else
   dataout = data
end if

dims = dimsizes(data1)
do n=0, dims(0)-2
 if(.not.ismissing(data1(n,0)))
    ;;; find data locations close to data1 locations
   rpi := ind( gc_latlon( data1(n,0), data1(n,1), dataout(:,0), dataout(:,1),0,4) .lt. rkm )
   if(.not.all(ismissing(rpi)))
    ;;;;; Just delete data locations
     dataout(rpi,:) = dataout@_FillValue
    ;;;;; Strip empty values and downsize array. 
     nmi := ind(.not.ismissing(dataout(:,0)))
     if(all(ismissing(nmi)))        ;; if no new locations left return
      return dataout
     end if
     if(dimsizes(nmi).eq.1)         ;; make sure singleton is returned as 2d array.
       dataout_ = new( (/1, dims(1)/), float)
       dataout_(0,:) = dataout(nmi,:)
       dataout := dataout_
     else
       dataout := dataout(nmi,:)  
     end if
      ;;;;;  end strip
   end if
  end if
end do

return dataout
end

undef("shift_grid")
function shift_grid(data, pLat, pLon,dy, dx)
local data, pLon, pLat, dy, dx, sub_locs, data_
begin
    
  ;;;;;   shift grid so we're not at the edge of a domain
  if(pLon-dx .lt. data&lon(0)  .or.  pLon+dx .gt. last(data&lon) )
      if( pLon-dx .lt. data&lon(0)  )
        data_ = lonPivot(data, pLon+180)
        data_&lon = make_lon_monotonic(data_&lon)
        data_&lon = data_&lon - 360.
      else
        data_ = lonPivot(data, pLon-180)
        data_&lon = make_lon_monotonic(data_&lon)
      end if
      sublocs := data_({pLat-dy:pLat+dy}, {pLon-dx:pLon+dx})
  else
  sublocs := data({pLat-dy:pLat+dy}, {pLon-dx:pLon+dx})
  end if

    return sublocs
end


undef("loc_max_dist")
function loc_max_dist(x, cyclic, dist, locs)
local x, cyclic, dist, locs,\
 maxs, dimnames, mlons, mlats, dists, dist_ind, nmax, new_max
begin
  maxs = local_max(x, cyclic, 0.)
  if(all(maxs.eq.0))
     return(0)
  end if
  dimnames = getVarDimNames(x)
  mlons =  x&$dimnames(1)$(maxs@xi)
  mlats =  x&$dimnames(0)$(maxs@xi)
  dists = gc_latlon(locs(0), locs(1), mlats, mlons,0,4)
  dist_ind = ind(dists .le. dist)
  if(all(ismissing(dist_ind)))
    return(0)
  end if
  nmax = dist_ind( maxind(maxs@maxval(dist_ind)))  
  new_max = 1
  new_max@val = maxs@maxval(nmax)
  new_max@$dimnames(1)$ = x&$dimnames(1)$(maxs@xi(nmax))
  new_max@$dimnames(0)$ = x&$dimnames(0)$(maxs@yi(nmax))
  return(new_max)
end



undef("loc_max")
function loc_max(x, cyclic, thres)
local maxs, new_max, max_thres
begin
  maxs = local_max(x, cyclic, 0.)
  if(all(maxs.eq.0))
     return(0)
  end if
  max_thres = ind(maxs@maxval .ge. thres )
  dimnames = getVarDimNames(x)
  if(all(ismissing(max_thres)))
      new_max = 0.
  else
    new_max = dimsizes(max_thres)
    new_max@val = maxs@maxval(max_thres)
    new_max@$dimnames(1)$ = x&$dimnames(1)$(maxs@xi(max_thres))
    new_max@$dimnames(0)$ = x&$dimnames(0)$(maxs@yi(max_thres))
  end if
  return(new_max)
end

undef("loc_min")
function loc_min(x, cyclic, thres)
local mins, new_min, min_thres
begin
  mins = local_min(x, cyclic, 0.)
  if(all(mins.eq.0))
     return(0)
  end if
  min_thres = ind(mins@minval .le. thres )
  dimnames = getVarDimNames(x)
  if(all(ismissing(min_thres)))
      new_min = 0.
  else
    new_min = dimsizes(min_thres)
    new_min@val = mins@minval(min_thres)
    new_min@$dimnames(1)$ = x&$dimnames(1)$(mins@xi(min_thres))
    new_min@$dimnames(0)$ = x&$dimnames(0)$(mins@yi(min_thres))
  end if
  return(new_min)
end



undef("getlatlon_points")
function getlatlon_points(grid, lat, lon)
local retval, lat_i, lon_i, nmi
begin
  retval = new( dimsizes(lat), typeof(grid) )
  nmi = ind(.not.ismissing(lat))
  if(.not.all(ismissing(nmi)))
    lon_i = ind_nearest_coord(lon(nmi), grid&lon, 0)
    lat_i = ind_nearest_coord(lat(nmi), grid&lat, 0)
    do i=0,dimsizes(nmi)-1
      retval(nmi(i)) = (/ grid(lat_i(i), lon_i(i) ) /)
    end do
  end if
  return(retval)
end


undef("subgrid_location_km")
function subgrid_location_km(data[*][*]:float, lat[1]:numeric, lon[1]:numeric, minima[1]:logical, km[1]:numeric)
local lati, loni, lap, lam, lop, lom, subgrid, lats, lons, alat, alon, aval, dists
begin

;delta = 2   ;  grid points

;lati = ind_nearest_coord(lat,data&lat, 0)
;loni = ind_nearest_coord(lon, data&lon, 0)
;lap = lati+delta < dimsizes(data&lat)-1
;lam = lati-delta > 0
;lop = loni+delta < dimsizes(data&lon)-1
;lom = loni-delta > 0

subgrid = data;(lam:lap, lom:lop) 
if(minima)  subgrid = subgrid*-1  end if
if(any(subgrid.lt.0)) subgrid = subgrid+abs(min(subgrid)) end if

lats = conform(subgrid, subgrid&lat, 0)
lons = conform(subgrid, subgrid&lon, 1)
dists = gc_latlon(lat, lon, lats, lons, 0, 4)
dists = ((km - dists) > 0)/ (0.5*km) < 1.
subgrid = subgrid * dists

alat = dim_avg_wgt(ndtooned(lats), ndtooned(subgrid), 1)
alon = dim_avg_wgt(ndtooned(lons), ndtooned(subgrid), 1)
aval = dim_avg_wgt(ndtooned(data), ndtooned(subgrid), 1)

return( (/alat, alon, aval/) )
end




undef("subgrid_location")
function subgrid_location(data[*][*]:float, lat[1]:numeric, lon[1]:numeric, minima[1]:logical, delta[1]:numeric)
local lati, loni, lap, lam, lop, lom, subgrid, lats, lons, alat, alon, aval
begin

;delta = 2   ;  grid points

lati = ind_nearest_coord(lat,data&lat, 0)
loni = ind_nearest_coord(lon, data&lon, 0)
lap = lati+delta < dimsizes(data&lat)-1
lam = lati-delta > 0
lop = loni+delta < dimsizes(data&lon)-1
lom = loni-delta > 0

subgrid = data(lam:lap, lom:lop) 
if(minima)  subgrid = subgrid*-1  end if
if(any(subgrid.lt.0)) subgrid = subgrid+abs(min(subgrid)) end if

lats = conform(subgrid, subgrid&lat, 0)
lons = conform(subgrid, subgrid&lon, 1)


alat = dim_avg_wgt(ndtooned(lats), ndtooned(subgrid), 1)
alon = dim_avg_wgt(ndtooned(lons), ndtooned(subgrid), 1)
aval = dim_avg_wgt(ndtooned(data(lam:lap, lom:lop)), ndtooned(subgrid), 1)

return( (/alat, alon, aval/) )
end


debug = False

undef("get_mass_center")
function get_mass_center(data, locs)
local minima, rkm,dimz, rankz, retval, pLat, pLon, dx, sublocs, lats, lons, val1d, dist
begin
dx = 10;
minima :=  get_res_value_keep(data, "minima", False) 
rkm :=  get_res_value_keep(data, "rkm", 400.) 

dimz = dimsizes(locs)
rankz = dimsizes(dimz)
if(rankz.eq.2)
  retval = new( (/dimz(0), dimz(1)>3/), float)
  do r=0, dimz(0)-1
   retval(r,:) = get_mass_center( data, locs(r,:) )
  end do
  return retval
else
  pLat = locs(0)
  pLon = locs(1)
  pVal = 0.
  retval = new( (/1, dimz(0)>3/), float)

  if(any(ismissing(locs(:1))))
                                    if(debug) print("Missing locs") end if
   return(retval)
  end if
    
  ;;;;;   shift grid so we're not at the edge of a domain
  if(pLon-dx .lt. data&lon(0)  .or.  pLon+dx .gt. last(data&lon) )
      if( pLon-dx .lt. data&lon(0)  )
        data_ = lonPivot(data, pLon+180)
        data_&lon = make_lon_monotonic(data_&lon)
        data_&lon = data_&lon - 360.
      else
        data_ = lonPivot(data, pLon-180)
        data_&lon = make_lon_monotonic(data_&lon)
      end if
      sublocs := data_({pLat-dx:pLat+dx}, {pLon-dx:pLon+dx})
  else
  sublocs := data({pLat-dx:pLat+dx}, {pLon-dx:pLon+dx})
  end if


  ;;;;;   Find mass center over large area first
;large_scale :=  subgrid_location(sublocs, pLat, pLon, minima, 7)
;large_scale :=  subgrid_location_km(sublocs, pLat, pLon, minima, 400)
;print(""+large_scale)
large_scale :=  subgrid_location_km(sublocs,  pLat, pLon,minima, rkm)
;large_scale :=  subgrid_location_km(sublocs, large_scale(0), large_scale(1), minima, 300)
sub_scale :=  subgrid_location_km(sublocs, large_scale(0), large_scale(1), minima, rkm/2)
;
; if(any(ismissing(large_scale)))
;                                   if(debug) print("Missing large scale") end if
;  return(retval)
; end if
;;
;;large_scale = locs
;  ;;;;;  Find center within smaller area; make sure it's within rkm distance. 
;  ;;;;;  if multiple take maxima/ minima 
;subsub := sublocs({large_scale(0)-6:large_scale(0)+6},{large_scale(1)-6:large_scale(1)+6})
;if( get_res_value_keep(data, "minima", False) )
;  if(debug) print("Searching for minima") end if
;
;  mini := loc_min(subsub, False, percentile(ndtooned(subsub),85.) )
;  if(debug) print("Found "+mini+" potential minima") end if
;
;  if(mini.eq.0) return(retval) end if
;  dist := gc_latlon(large_scale(0), large_scale(1), mini@lat, mini@lon,0,4)
;;  dist = gc_latlon(pLat, pLon, mini@lat, mini@lon,0,4)
;  mi :=ind(dist .le. rkm)
; if(all(ismissing(mi))) return(retval) end if
;  mi := mi( minind(mini@val(mi)) )
;  do r=150,rkm,150
;    if(any(dist .lt. r))
;      di := ind(dist .lt. r)
;      mi := di( minind(mini@val(di)))
;    end if
;  end do
;else
;  mini := loc_max(subsub, False, percentile(ndtooned(subsub),85.))
;  if(debug) print("Found "+mini+" potential maxima") end if
;  if(mini.eq.0) return(retval) end if
;;  dist = gc_latlon(pLat(0), locs(1), mini@lat, mini@lon,0,4)
;  dist := gc_latlon(large_scale(0), large_scale(1), mini@lat, mini@lon,0,4)
;  mi :=ind( dist .le. rkm)
;  if(all(ismissing(mi))) return(retval) end if
;  do rk=150,rkm,150
;;  print(r)
;  if(any(dist .lt. rk))
;    di := ind(dist .lt. rk)
;    mi := di( maxind(mini@val(di)))
;;    print(mini@lat(mi))
;    break
;  end if
;  end do
;end if
;mark = gsn_add_polymarker(wks,  map, mini@lon(mi(0)),  mini@lat(mi(0)), False)
;mark = gsn_add_polymarker(wks,  map, large_scale(1),  large_scale, False)
;df(wks, map)
;;print(mini@lat(mi))
;;;  find mass center over small region. 
; sub_scale = subgrid_location(sublocs, mini@lat(mi(0)), mini@lon(mi(0)), minima, 2)
 retval(0,:2)  = sub_scale

  return(   retval )
end if
 
end





undef("is_minima")
function is_minima(mslp, lat,lon, dist)
local var_sub, var, loc_min, lmi, min_lat,min_lon, ret
begin
mslp@minmima = True
retval = get_mass_center(mslp, (/lat, lon/) )
min_lat = retval(0,0)
min_lon = retval(0,1)

if( gc_latlon(lat, lon, min_lat, min_lon, 5, 4).lt. dist)
  ret = True
  ret@lat = min_lat
  ret@lon = min_lon
  ret@val = retval(0,2)
  return ret
else
  return False
end if
end

undef("is_maxima")
function is_maxima(mslp, lat,lon, dist)
local var_sub, loc_max, lmi, locs, min_lat, ret, mslp, lat, lon, dist,retval
begin
;print(lat+"   "+lon)
retval = get_mass_center(mslp, (/lat, lon/) )
min_lat = retval(0,0)
min_lon = retval(0,1)

if( gc_latlon(lat, lon, min_lat, min_lon, 5, 4).lt. dist .and. retval(0,2) .gt. 0. )
  ret = True
  ret@lat = min_lat
  ret@lon = min_lon
  ret@val = retval(0,2)
  return ret
else
  return False
end if
end


undef("is_env_maxima")
function is_env_maxima(ght, lat,lon)
local ght, gh_env, gh_core,ght_sub, loc_min,lmi, min_lat,min_lon, ret
begin
ght!0 = "lat"
ght!1 = "lon"
gh_env = circle_avg_m_point(ght, 1000., lat, lon)
gh_core= circle_avg_m_point(ght, 100., lat, lon)

if(gh_core .gt. gh_env ) ;return False end if
  ret = True
  ret@val = gh_core - gh_env
  return ret
else
  return False
end if
end



;undef("check_closed_contour")
;function check_closed_contour(data, in_locs, critvalue)
;local fake_wks, fakeres, fakeplot, dimz, rankz, nloc, retval, center_loc,\
;    crits, loop, interval, cVal, iso, j,b,e,x,y, critvalue, data, closed
;begin
;
;fake_wks = create "fake" imageWorkstationClass noparent
;"wkHeight":2000
;"wkWidth":2000
;end create
;
;;fake_wks = gsn_open_wks("png","")
;fakeres= True
;fakeres@gsnDraw = False
;fakeres@gsnFrame = False
;fakeplot = gsn_csm_contour(fake_wks, data, fakeres)
;
;dimz = dimsizes(in_locs)
;rankz = dimsizes(dimz)
;if(rankz .eq.2)
;  nloc = dimz(0)
;  retval = new(nloc, logical)
;  center_loc = in_locs
;else
;  retval = new( 1, logical)
;  nloc = 1
;  center_loc = new((/1, dimz(0)/), float)
;  center_loc(0,:) = in_locs
;end if
;
;crits = new(nloc, float)
;crits = tofloat(critvalue)
;
;loop = get_res_value_keep(critvalue, "loop", False)
;interval = get_res_value_keep(critvalue, "interval", 1.)
;max_size = get_res_value_keep(critvalue, "max_size", 5e6)
;min_size = get_res_value_keep(critvalue, "min_size", 1e5)
;
;area = new(dimsizes(retval), float)
;keep_start = new(dimsizes(retval), integer)
;keep_size = new(dimsizes(retval), integer)
;keep_y = 0.
;keep_x = 0.
;
;do l=0, nloc-1
; if(ismissing(center_loc(l,0))) continue end if
;
;  cVal = crits(l)
;  ;closed = False
;  do while( (cVal .le. crits(l)) .or. (loop .and. closed)   )
;   iso := get_isolines(fakeplot, cVal)
;   do j = 0, iso@segment_count -1
;    b = iso@start_point(j)
;    e = b + iso@n_points(j) - 1
;    if(iso(0,b).eq.iso(0,e) .and. iso(1,b).eq.iso(1,e))
;      y :=iso(0,b:e)
;      x := iso(1,b:e)
;      closed = gc_inout(center_loc(l,0), center_loc(l,1), y, x)
;      if(closed)
;        if( .not. gc_clkwise(y, x) )
;          y = y(::-1)
;          x = x(::-1)
;        end if
;       area(l) = area_poly_sphere(y, x, 6371)
;        if( (area(l) .gt. min_size) .and. (area(l) .lt. max_size))
;          retval(l) = True
;          keep_start(l) = dimsizes(keep_y)
;          keep_size(l) = dimsizes(y)
;          keep_y := array_append_record(keep_y, y,0)
;          keep_x := array_append_record(keep_x, x,0)
;   ;      retval@lats := y
;   ;      retval@lons := x
;   ;      retval@pres = cVal
;          break
;       end if
;      end if
;    end if
;   end do
;   cVal := cVal + interval
;  end do
;end do
;retval@area = area
;retval@start_point = keep_start
;retval@n_points = keep_size
;retval@y = keep_y
;retval@x = keep_x
;;delete(fake_wks)
;return( retval )
;end



undef("get_mand_positions")
;;;;;;;;;;;;; Get 6 center locations using same variables as Marchok GFDL
;;;;;;;;        ;; 925 used instead of 10m. 
;;;;;;;; 
;;;;;;;; 
function get_mand_positions(in_time, in_lat, in_lon, rkm)
local vars, opts, levs, minima, mins, mand_data, rkm, temp_min, opt,data, locs, in_time, in_lat, in_lon, rkm
begin

vars = (/"vr", "pmsl","g","g","vr","vr", "uv", "uv"/)
opts := (/True,False,  True,True, True, True, True, True/)
levs = (/ 850, -999.,  850,  700, 700,   925, 850, 700/)
minima=(/False,True,  True, True, False, False, True, True/)
rkms =(/ 500., 350, 350, 350, 350, 350, 125, 125/)

dimz = dimsizes(in_lat)
mand_data = new( (/dimz(0), 3,dimsizes(vars)/), float)
mand_data@_FillValue = default_fillvalue("float")
if(.not. all(ismissing(in_lat)))
do v=0, dimsizes(vars)-1
  opt = opts(v)
  opt@levs = levs(v)
  data := grab_data(in_time, vars(v), opt)  ;; read in global grid from file
  data@minima = minima(v)
  if(in_lat .lt. 0 .and. vars(v).eq."vr" )
    data@minima = True      ; Things spin the wrong way. 
  end if
  data@rkm = rkms(v)

    locs = new( (/dimz(0), 2/), float)
    locs(:,0) = in_lat
    locs(:,1) = in_lon
    mand_data(:,:2,v) = get_mass_center(data, locs)    
    if(v.eq.0 .and. any(.not.ismissing( mand_data(:,0,0))) )
     nmi := ind(.not.ismissing( mand_data(:,0,0)))
     in_lat(nmi) = dim_avg_n((/in_lat(nmi), mand_data(nmi,0,0)/),0)
     in_lon(nmi) = dim_avg_n((/in_lon(nmi), mand_data(nmi,1,0)/),0)
    end if

;  end if
end do
end if
mand_data!0 = "systems"
mand_data!1 = "data"
mand_data!2= "variables"
mand_data&data = (/"lat","lon","val"/)
mand_data&variables = vars


return (mand_data)
end 



undef("get_wave_positions")
;;;;;;;;;;;;; Get extra 6 center locations using same curv. vort and OW
;;;;;;;;      
;;;;;;;; 
;;;;;;;; 
function get_wave_positions(in_time, in_lat, in_lon, rkm)
local vars, opts, levs, minima, mins, mand_data, rkm, temp_min, opt,v
begin

vars = (/"cv",  "cv",   "cv",   "g",    "vr",    "vr"/)
opts := (/ True, True,  True,    True,   True,     True/)
levs = (/ 850,   700,   500,    850,    850,     700/)
r = (/ 350.,     350,   350,    350,    250,     250/)
minima =(/False, False, False, True,    False,  False/)

mand_data = new( (/dimsizes(in_lat), 3, dimsizes(vars)/), float)

if(.not. all(ismissing(in_lat)))
  do v=0, dimsizes(vars)-1
    opt = opts(v)
    opt@levs = levs(v)
    copy_VarAtts(rkm, opt)
    data := grab_data(in_time, vars(v), opt)  ;; read in global grid from file
    data!0 = "lat"
    data!1 = "lon"
;    pvs(data)
;    data = circle_avg_m(data, r(v))
    do l=0, dimsizes(in_lat)-1
      if(.not. ( ismissing(in_lat(l)) .and. ismissing(in_lon(l))) )
        if(minima(v))
          temp_min := is_minima(data, in_lat(l), in_lon(l), rkm) 
        else
          temp_min := is_maxima(data, in_lat(l), in_lon(l), rkm)
        end if
        mand_data(l,0,v)  = get_res_value(temp_min, "lat", mand_data@_FillValue)
        mand_data(l,1,v)  = get_res_value(temp_min, "lon", mand_data@_FillValue)
        mand_data(l,2,v)  = get_res_value(temp_min, "val", mand_data@_FillValue)  
      end if
    end do
  end do
end if

mand_data!0 = "systems"
mand_data!1= "data"
mand_data!2= "variables"
mand_data&data = (/"lat","lon","val"/)
mand_data&variables = vars

return (mand_data)
end 







undef("get_starting_targets")
function get_starting_targets(time, lat_bds[2], lon_bds[2], opt)
local min_thresh, cv, uwnd, vwnd, maxes, maxes_mass, easterly, tmp_wave_locs
begin
rkm = get_res_value_keep(opt, "rkm",300)
min_thresh = get_res_value_keep(opt, "min_threshold", 0.1)

cv := grab_data(time, "ow", opt)
cv = circle_avg_m(cv, rkm)
uwnd := grab_data(time, "u", opt)
vwnd := grab_data(time, "v", opt)


maxes := loc_max(cv({lat_bds(0):lat_bds(1)},{lon_bds(0):lon_bds(1)}), True, min_thresh)
maxes_mass := new( (/maxes, 3/), float)
do m=0, maxes-1
 maxes_mass(m,:) =  get_mass_center(cv, (/maxes@lat(m), maxes@lon(m)/) )
end do

;;;;;;;;;;  Remove duplicate locations
maxes_mass := unique_locations(maxes_mass, 100.)
maxes_mass := maxes_mass(ind(.not.ismissing(maxes_mass(:,0))),:)        
;;;;;;;;;;;;;;;;;;;;;

maxes_uwnd := circle_avg_m_point(uwnd, 500, maxes_mass(:,0), maxes_mass(:,1) )
maxes_vwnd := circle_avg_m_point(vwnd, 500, maxes_mass(:,0), maxes_mass(:,1) )
easterly := ind(maxes_uwnd .lt. 1)     ;; easterly with a bit of wiggle room
maxes_mass := maxes_mass(easterly,:)

;;;;;;;;;;;;;   Recenter potential targets on combined wave location  ;;; too slow to do here. 
;tmp_wave_locs := get_wave_positions(time, maxes_mass(:,0), maxes_mass(:,1), 600)
;maxes_mass := combine_mand_positions( tmp_wave_locs )
;maxes_mass := unique_locations(maxes_mass, 250.)

return maxes_mass
end


undef("zScore")
;;;;;;;;;;;;;  Modified z-score to detect outliers
;;;;;;;; 
;Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and Handle Outliers", The ASQC
;Basic References in Quality Control: Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
;;;;;;;; 
;;;;;;;; 
function zScore(lat, lon)
local  MAD, Mi, lat, lon, dist, nmi
begin
nmi = ind(.not. ismissing(lat))
if(all(ismissing(nmi))) return default_fillvalue("float") end if
    
Mi = new(dimsizes(lat), float, default_fillvalue("float"))
dist = new(dimsizes(lat), float, default_fillvalue("float"))

med_lat = dim_median(lat(nmi))
med_lon = dim_median(lon(nmi))

dist(nmi) = gc_latlon(med_lat, med_lon, lat(nmi), lon(nmi),0,4)        ;; absolute deviation from median center
MAD =  dim_median(dist)  > 5.  <100.           ;; Median Absolute Deviation ; min deviation of 10km so up-to 50km deviation allowed
Mi = (0.6745*dist)/MAD      ;; > 3.5 Considered outlier. 

;print(dist+"    "+Mi)

return(Mi)
end


undef("combine_mand_positions")
function combine_mand_positions(mand_data)
local mand_data, z_score, avg_loc, rankz, dimz, retval
begin
dimz = dimsizes(mand_data)
rankz = dimsizes(dimz)

if(rankz.eq.3)
  retval = new( (/dimz(0), 2/), float)
  do r=0, dimz(0)-1
   retval(r,:) = combine_mand_positions( mand_data(r,:,:) )
  end do
  return retval
else
if( (num(.not.ismissing(mand_data(0,:5))).le.3 ).or. (all(ismissing(mand_data(0,:1)))) )
    return( (/mand_data@_FillValue, mand_data@_FillValue/) )
 end if 

    zinds := ind( gc_latlon(avg(mand_data(0,:1)), avg(mand_data(1,:1)) , mand_data(0,:), mand_data(1,:), 0, 4) .lt. 350)
    if(all(ismissing(zinds)))
	  zinds := ind( gc_latlon( mand_data(0,0), mand_data(1,0) , mand_data(0,:), mand_data(1,:), 0, 4) .lt. 350)
 	end if 
    if(dimsizes(zinds).lt.3)
        return( (/mand_data@_FillValue, mand_data@_FillValue/) )
    end if 
    avg_loc = dim_avg_n( mand_data(:1,zinds),1)
    dists := gc_latlon(avg_loc(0), avg_loc(1), mand_data(0,6:7), mand_data(1,6:7), 0,4)
    dinds := ind(dists.le.125)
    if(.not.all(ismissing(dinds)))
        zinds := array_append_record(zinds, dinds, 0)
        avg_loc = dim_avg_n( mand_data(:1,zinds),1)
    end if
    return(avg_loc)
;  end if
  
  
end if
end 

undef("combine_wave_positions")
function combine_wave_positions(mand_data)
local mand_data, z_score, avg_loc, rankz, dimz, retval
begin
dimz = dimsizes(mand_data)
rankz = dimsizes(dimz)

if(rankz.eq.3)
  retval = new( (/dimz(0), 2/), float)
  do r=0, dimz(0)-1
   retval(r,:) = combine_wave_positions( mand_data(r,:,:) )
  end do
  return retval
else
if( num(.not.ismissing(mand_data(0,(/0,1,2,4,5/)))).lt.2 ) 
    return( (/mand_data@_FillValue, mand_data@_FillValue/) )
 end if 
 zinds := ind( gc_latlon(avg(mand_data(0,:)), avg(mand_data(1,:)) , mand_data(0,:), mand_data(1,:), 0, 4) .lt. 350)
 if(all(ismissing(zinds)))
   zinds := ind( gc_latlon( mand_data(0,0), mand_data(1,0) , mand_data(0,:), mand_data(1,:), 0, 4) .lt. 350)
 end if 
 if(dimsizes(zinds).lt.2)
     return( (/mand_data@_FillValue, mand_data@_FillValue/) )
 end if 
 avg_loc = dim_avg_n( mand_data(:1,zinds),1)
 return(avg_loc)
;  end if
end if
end 




undef("get_wind_speed")
function get_wind_speed(in_time, in_lat, in_lon)
local opt, u, v, uv, temp_min, maxwnd, in_time, in_lat, in_lon
begin
    opt := True
    opt@levs = (/925/)
    opt@levSelection = "Explicit"
    uv := grab_data(in_time, "uv", opt)
    temp_min := is_maxima(uv, in_lat, in_lon, 250.)
return(temp_min)
end


undef("get_aux_variables")
function get_aux_variables(in_time, in_lat, in_lon)
local opt, hgt, hgtd
begin

;;;;;   Thickness Anomaly   ;; is Warm core? ;;;;; ;;;;; ;;;;; 
    opt := True
    opt@levs = (/200,850/)
    opt@levSelection = "Explicit"
    hgt := grab_data(in_time, "g", opt)
    hgtd := hgt(0,:,:) - hgt(1,:,:)
    copy_VarCoords(hgt(0,:,:), hgtd)
    thick_max := is_env_maxima(hgtd, in_lat,in_lon)
;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; 

;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; 
    data := grab_data(in_time, "pmsl", False)
    temp_min := is_minima(data, in_lat, in_lon, 150.) 
    if(temp_min)
    crit := temp_min@val
    crit = crit+200.
    crit@loop = True
    crit@interval = 200.
    mslp_closed := check_closed_contour(data, (/ in_lat, in_lon/),crit)
    else
      mslp_closed = False
    end if
;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; 


;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; 
    opt := True
    opt@levs = (/600,900/)
    opt@levSelection = "Explicit"
    hgt := grab_data(in_time, "g", opt)
    hgtd := hgt(0,:,:) - hgt(1,:,:)
    copy_VarCoords(hgt(0,:,:), hgtd)
    vtl := is_env_maxima(hgtd, in_lat,in_lon)
;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; 

    opt := True
    opt@levs = (/300,600/)
    opt@levSelection = "Explicit"
    hgt := grab_data(in_time, "g", opt)
    hgtd := hgt(0,:,:) - hgt(1,:,:)
    copy_VarCoords(hgt(0,:,:), hgtd)
    vtu := is_env_maxima(hgtd, in_lat,in_lon)
;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; 

aux_data = new( 4, float)
aux_data(0) = get_res_value(thick_max, "val", aux_data@_FillValue)
aux_data(1) = get_res_value(mslp_closed, "pres", aux_data@_FillValue)
aux_data(2) = get_res_value(vtl, "val", aux_data@_FillValue)
aux_data(3) = get_res_value(vtu, "val", aux_data@_FillValue)

return(aux_data)
end 

undef("advect_location")
function advect_location(mand_data, ups, vps, dt)
local dims, ndims, new_data, lons, mand_data
begin
dims = dimsizes(mand_data)
ndims = dimsizes(dims)

new_data = mand_data
if( ndims .eq. 1)
  new_data(0) = new_data(0) + ((vps*dt) /111120.)
  new_data(1) = new_data(1)  + ((ups*dt) /111120.)*cos(torad(new_data(0) ))
else
  new_data(:,0) = new_data(:,0) + ((vps*dt) /111120.)
  new_data(:,1) = new_data(:,1)  + ((ups*dt) /111120.)*cos(torad(new_data(:,0) ))
end if

return( new_data )
end

