undef("ab_debug")
procedure ab_debug(var)
begin
if(isdefined("debug_") .and. debug_)
  if( isstring(var) )
    print((/var/))
  end if
  if( isnumeric(var) )
    if( product(dimsizes(var)) .lt. 20 )
      print(var)
    else
      printVarSummary(var)
    end if
  end if
end if
end 



undef("first")
function first(x)
local nmi
begin
 nmi = ind(.not.ismissing(x))
 return(x(nmi(0)))
end

undef("last")
function last(x)
local dims,nmi
 begin
 nmi = ind(.not.ismissing(x))
 dims = dimsizes(nmi)
 return(x(nmi(dims(0)-1)))
end


undef("dum_cumcount")
function dum_cumcount(x[*][*]:logical)
local x, nx, y, i
begin
  nx := dimsizes(x)
  y := new( nx, integer)
  y = 0
  do i=0,nx(1)-1
    y(:,i) = where( x(:,i), y(:,(i-1>0))+1, 0)
  end do
  copy_VarCoords(x, y)
return y
end




undef("unwrap_lon")
function unwrap_lon(lon, thresh)
local lon, dl, nl,isM,disc,diff,nmi
begin

nl = dimsizes(lon)
nmi := ind(.not.ismissing(lon))

do d=0,dimsizes(nmi)-2

  if( lon(nmi(d+1)) - lon(nmi(d)).lt. thresh)
    lon(nmi(d+1))  = lon(nmi(d+1))+360
  end if
  if(lon(nmi(d+1)) - lon(nmi(d)).gt. thresh)
    lon(nmi(d+1))  = lon(nmi(d+1))-360
  end if
end do

if(any(lon.lt.0))
 lon = lon + 360
end if

return(lon)
end

undef("wrap_lon")
function wrap_lon(lon)
local lon, dl, nl,isM,disc,diff,nmi
begin

do while(any(lon.gt.180))
lon = where(lon.gt.180, lon-360, lon)
end do

do while(any(lon.lt.-180))
lon = where(lon.lt.-180, lon+360, lon)
end do


return(lon)
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


undef("loc_max_dist")
function loc_max_dist(x, cyclic, dist, inlocs)
local x, cyclic, dist, locs,\
 maxs, dimnames, mlons, mlats, dists, dist_ind, nmax, new_max,inlocs
begin

  maxs = local_max(x, cyclic, 0.)
  if(all(maxs.eq.0))
     return(0)
  end if
  dimnames = getVarDimNames(x)
  mlons =  x&$dimnames(1)$(maxs@xi)
  mlats =  x&$dimnames(0)$(maxs@yi)
  dimz = dimsizes(inlocs)
  rankz = dimsizes(dimz)
  if(rankz.eq.1)
  n=1
  locs = new( (/2, 1/), float)
  locs(:,0) = inlocs
  else
  n=dimz(1)
  locs = inlocs
  end if
  retval := new(n, logical)
  retval = False
  retval@val = new(n, float)
  retval@$dimnames(0)$ =  new(n, float)
  retval@$dimnames(1)$ =  new(n, float)
  do i=0,n-1
    dists = gc_latlon(locs(0,i), locs(1,i), mlats, mlons,0,4)
    dist_ind := ind(dists .le. dist)
     if(all(ismissing(dist_ind)))
       continue
     end if
    nmax := dist_ind( maxind(maxs@maxval(dist_ind)))  
    retval(i) = True
    retval@val(i) = maxs@maxval(nmax)
    retval@$dimnames(1)$(i) = x&$dimnames(1)$(maxs@xi(nmax))
    retval@$dimnames(0)$(i) = x&$dimnames(0)$(maxs@yi(nmax))
  end do
  
  return(retval)
end

undef("loc_min_dist")
function loc_min_dist(data, cyclic, dist, inlocs)
local data, cyclic, dist, locs,\
 maxs, dimnames, mlons, mlats, dists, dist_ind, nmax, new_max,inlocs
begin

  maxs = local_min(data, cyclic, 0.)
  maxs@val = maxs@minval
  if(all(maxs.eq.0))
     return(0)
  end if
  dimnames = getVarDimNames(data)
  mlons :=  data&$dimnames(1)$(maxs@xi)
  mlats :=  data&$dimnames(0)$(maxs@yi)
  dimz = dimsizes(inlocs)
  rankz = dimsizes(dimz)
  if(rankz.eq.1)
  n=1
  locs = new( (/2, 1/), float)
  locs(:,0) = inlocs
  else
  n=dimz(1)
  locs = inlocs
  end if
  retval := new(n, logical)
  retval = False
  retval@val = new(n, float)
  retval@$dimnames(0)$ =  new(n, float)
  retval@$dimnames(1)$ =  new(n, float)
  do i=0,n-1
    dists := gc_latlon(locs(0,i), locs(1,i), mlats, mlons,0,4)
    dist_ind := ind(dists .le. dist)
     if(all(ismissing(dist_ind)))
       continue
     end if
    nmax := dist_ind( minind(maxs@val(dist_ind)))  
    retval(i) = True
    retval@val(i) = maxs@val(nmax)
    retval@$dimnames(1)$(i) = data&$dimnames(1)$(maxs@xi(nmax))
    retval@$dimnames(0)$(i) = data&$dimnames(0)$(maxs@yi(nmax))
  end do
  
  return(retval)
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
     new_max = 0
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
local mins, new_min, min_thres, dimnames
begin
 mins = local_min(x, cyclic, 0.)
 if(all(mins.eq.0))
    return(0)
 end if
 min_thres = ind(mins@minval .le. thres )
 dimnames = getVarDimNames(x)
 if(all(ismissing(min_thres)))
     new_min = 0
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

subgrid = data;(lam:lap, lom:lop) 
if(minima)  subgrid = subgrid*-1  end if
if(any(subgrid.lt.0)) subgrid = subgrid+abs(min(subgrid)) end if

lats = conform(subgrid, subgrid&lat, 0)
lons = conform(subgrid, subgrid&lon, 1)
dists := gc_latlon(lat, lon, lats, lons, 0, 4)
subgrid = where(dists .gt. (km > min(dists)), 0, subgrid^2)

alat = dim_avg_wgt(ndtooned(lats), ndtooned(subgrid), 1)
alon = dim_avg_wgt(ndtooned(lons), ndtooned(subgrid), 1)
if(.not. any( ismissing( (/alat,alon/) ) ))
  aval = data({alat}, {alon})
else
  aval = new(1, float)
end if
return( (/alat, alon, aval/) )
end



;
;undef("subgrid_location")
;function subgrid_location(data[*][*]:float, lat[1]:numeric, lon[1]:numeric, minima[1]:logical, delta[1]:numeric)
;local lati, loni, lap, lam, lop, lom, subgrid, lats, lons, alat, alon, aval
;begin
;
;;delta = 2   ;  grid points
;
;lati = ind_nearest_coord(lat,data&lat, 0)
;loni = ind_nearest_coord(lon, data&lon, 0)
;lap = lati+delta < dimsizes(data&lat)-1
;lam = lati-delta > 0
;lop = loni+delta < dimsizes(data&lon)-1
;lom = loni-delta > 0
;
;subgrid = data(lam:lap, lom:lop) 
;if(minima)  subgrid = subgrid*-1  end if
;if(any(subgrid.lt.0)) subgrid = subgrid+abs(min(subgrid)) end if
;
;lats = conform(subgrid, subgrid&lat, 0)
;lons = conform(subgrid, subgrid&lon, 1)
;
;
;alat = dim_avg_wgt(ndtooned(lats), ndtooned(subgrid), 1)
;alon = dim_avg_wgt(ndtooned(lons), ndtooned(subgrid), 1)
;aval = dim_avg_wgt(ndtooned(data(lam:lap, lom:lop)), ndtooned(subgrid), 1)
;
;return( (/alat, alon, aval/) )
;end


 
 
undef("get_mass_center")
function get_mass_center(data, locs)
local minima, rkm,dimz, rankz, retval, pLat, pLon, dx, sublocs, lats, lons, val1d, dist
begin

minima :=  get_res_value_keep(data, "minima", False) 
rkm :=  get_res_value_keep(data, "rkm", 400.) 
min_rkm :=  get_res_value_keep(data, "min_rkm", 100.) 

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
   ab_debug("Missing locs")
   return(retval)
  end if
    
  ;;;;;   shift grid so we're not at the edge of a domain
  dx = 10;
  sublocs := shift_grid(data, pLat, pLon,dx, dx)

  ;;;;;   Find mass center over large area first
  n=1
  sub_scale :=  subgrid_location_km(sublocs,  pLat, pLon,minima, rkm)  
        ab_debug(sprintf("%6.4f",sub_scale(0))+"     "+sprintf("%6.4f",sub_scale(1)))
  sub_scale :=  subgrid_location_km(sublocs, sub_scale(0), sub_scale(1), minima, rkm/n)
        ab_debug(sprintf("%6.4f",sub_scale(0))+"     "+sprintf("%6.4f",sub_scale(1)))
; check_centerfix(wks, sublocs, sub_scale )

  n=2 
  do while( rkm/(n*1.) .gt. min_rkm )
        ab_debug(""+ rkm/(n-1))
   sub_scale :=  subgrid_location_km(sublocs, sub_scale(0), sub_scale(1), minima, rkm/n)
        ab_debug(sprintf("%6.4f",sub_scale(0))+"     "+sprintf("%6.4f",sub_scale(1)))
   sub_scale :=  subgrid_location_km(sublocs, sub_scale(0), sub_scale(1), minima, rkm/n)
        ab_debug(sprintf("%6.4f",sub_scale(0))+"     "+sprintf("%6.4f",sub_scale(1)))
;    check_centerfix(wks, sublocs, sub_scale )
   n=n+1
  end do

 sub_scale :=  subgrid_location_km(sublocs, sub_scale(0), sub_scale(1), minima, min_rkm)
        ab_debug(sprintf("%6.4f",sub_scale(0))+"     "+sprintf("%6.4f",sub_scale(1)))
 sub_scale :=  subgrid_location_km(sublocs, sub_scale(0), sub_scale(1), minima, min_rkm)
        ab_debug(sprintf("%6.4f",sub_scale(0))+"     "+sprintf("%6.4f",sub_scale(1)))
  
  if( sub_scale(1) .lt. min(data&lon))
       sub_scale(1) =  sub_scale(1)+360.
  end if
  if( sub_scale(1) .gt. max(data&lon))
      sub_scale(1) =  sub_scale(1)-360.
  end if

  retval(0,:2)  = sub_scale

  return(   retval )
end if
 
end





undef("is_minima")
function is_minima(data, in_lat,in_lon, dist)
local var_sub, var, loc_min, lmi, min_lat,min_lon, retval, ret,data, lat, lon
begin
data@minima = True
retval :=  loc_min_dist(data, True, dist, (/in_lat, in_lon/) )
if(any(retval))

  data@rkm = dist 
  retval := get_mass_center(data, (/retval@lat, retval@lon/) )
  min_lat = retval(0,0)
  min_lon = retval(0,1)

  if(.not.any(ismissing(retval(0,:))) .and.  gc_latlon(in_lat, in_lon, min_lat, min_lon, 5, 4).lt. dist)
   ret = True
   ret@lat = min_lat
   ret@lon = min_lon
   ret@val = retval(0,2)
   return ret
  end if
end if
return False
end

undef("is_maxima")
function is_maxima(data, in_lat,in_lon, dist)
local var_sub, loc_max, lmi, locs, min_lat, ret, mslp, lat, lon, dist,retval
begin

retval =  loc_max_dist(data, True, dist, (/in_lat, in_lon/) )
if(any(retval))
 data@rkm = dist 
 retval := get_mass_center(data, (/retval@lat, retval@lon/) )
 min_lat = retval(0,0)
 min_lon = retval(0,1)

  if(.not.any(ismissing(retval(0,:))) .and.  gc_latlon(in_lat, in_lon, min_lat, min_lon, 5, 4).lt. dist)
   ret = True
   ret@lat = min_lat
   ret@lon = min_lon
   ret@val = retval(0,2)
   return ret
  end if
 end if
 return False
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

undef("get_mand_positions")
;;;;;;;;;;;;; Get 6 center locations using same variables as Marchok GFDL
;;;;;;;;        ;; 925 used instead of 10m. 
;;;;;;;; 
;;;;;;;; 
function get_mand_positions(in_time, in_lat, in_lon, rkm)
local vars, opts, levs, minima, mins, mand_data, rkm, temp_min, opt,data, locs, in_time, in_lat, in_lon, rkm
begin

vars = (/  "pmsl","vr",  "g","g","vr","vr", "uv", "uv"/)
opts := (/False,   True,True,True, True, True, True, True/)
levs = (/  -999., 850,  850,  700, 700,   925, 850, 700/)
minima=(/ True,   False,True, True, False, False, True, True/)
rkms =(/  300,    300., 300, 300, 300, 300, 125, 125/)

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

    do l=0, dimsizes(in_lat)-1
      if(.not. ( ismissing(in_lat(l)) .and. ismissing(in_lon(l))) )
        if(  data@minima )
          temp_min := is_minima(data, in_lat(l), in_lon(l), data@rkm) 
        else
          temp_min := is_maxima(data, in_lat(l), in_lon(l), data@rkm)
        end if
;      if(.not.ismissing(temp_min@lat))
;         check_centerfix(wks, data, (/temp_min@lat, temp_min@lon/))
;      else
;        check_centerfix(wks, data, (/in_lat(l), in_lon(l)/))
;     end if
        mand_data(l,0,v)  = get_res_value(temp_min, "lat", mand_data@_FillValue)
        mand_data(l,1,v)  = get_res_value(temp_min, "lon", mand_data@_FillValue)
        mand_data(l,2,v)  = get_res_value(temp_min, "val", mand_data@_FillValue)  
      end if
    end do
  end do
end if
mand_data!0 = "systems"
mand_data!1 = "data"
mand_data!2= "variables"
mand_data&data = (/"lat","lon","val"/)
mand_data&variables = vars


return (mand_data)
end 



undef("combine_wave_positions")
function combine_wave_positions(mand_data)
local mand_data, z_score, avg_loc, rankz, dimz, retval
begin
dimz = dimsizes(mand_data)
rankz = dimsizes(dimz)
ab_debug(" ---- combine_wave_positions   ----")
if(rankz.eq.3)
  retval = new( (/dimz(0), 2/), float)
  do r=0, dimz(0)-1
;    print("   "+r)
   retval(r,:) = combine_wave_positions( mand_data(r,:,:) )
  end do
  return retval
else
if( num(.not.ismissing(mand_data(0,(/0,1,2,4,5/)))).le.2 ) 
    ab_debug("Less than 2 center fixes found")
    return( (/mand_data@_FillValue, mand_data@_FillValue/) )
 end if 
 mand_data(1,:) = unwrap_lon(mand_data(1,:), 50.) ;; fix dateline issue.
 zinds := ind( gc_latlon(avg(mand_data(0,:)), avg(mand_data(1,:)) , mand_data(0,:), mand_data(1,:), 0, 4) .lt. 350)
;print(zinds +"  ")
 if(all(ismissing(zinds)))
       ab_debug("No centers within 350 km")
   zinds := ind( gc_latlon( mand_data(0,0), mand_data(1,0) , mand_data(0,:), mand_data(1,:), 0, 4) .lt. 350)
 end if 
 if(dimsizes(zinds).le.2)
        ab_debug("Only 1 center within 350 km")
     return( (/mand_data@_FillValue, mand_data@_FillValue/) )
 end if 
 avg_loc = dim_avg_n( mand_data(:1,zinds),1)
 return(avg_loc)
;  end if
end if
end 



undef("get_wave_positions")
;;;;;;;;;;;;; Get extra 6 center locations using same curv. vort and OW
;;;;;;;;      
;;;;;;;; 
;;;;;;;; 
function get_wave_positions(in_time, in_lat, in_lon, options)
local vars, opts, levs, minima, mins, mand_data, rkm, temp_min, opt,v,nmi
begin

vars = (/"cv",  "cv",   "cv",   "g",    "vr",    "vr"/)
opts := (/ True, True,  True,    True,   True,     True/)
levs = (/ 850,   700,   500,    850,    850,     700/)
rkms = (/ 300.,     300,   300,    200,    200,     200/)
minima =(/False, False, False, True,    False,  False/)

dimz = dimsizes(in_lat)
mand_data = new( (/dimz(0), 3,dimsizes(vars)/), float)
mand_data@_FillValue = default_fillvalue("float")

if(.not. all(ismissing(in_lat)))
  do v=0, dimsizes(vars)-1
    opt := opts(v)
    copy_VarAtts(options, opt)
    opt@levs := levs(v)
    data := grab_data(in_time, vars(v), opt)  ;; read in global grid from file
    opt@rkm = get_res_value_keep(options, "rkm", rkms(v))
    data!0 = "lat"
    data!1 = "lon"
    data := data({50:-20},:)
    data@min_rkm = get_res_value_keep(options, "min_rkm", 55)
    data@rkm = rkms(v)
    data = circle_avg_m(data,  data@rkm)
    do l=0, dimsizes(in_lat)-1
      if(.not. ( ismissing(in_lat(l)) .and. ismissing(in_lon(l))) )
        if(minima(v))
          temp_min := is_minima(data, in_lat(l), in_lon(l), opt@rkm) 
        else
          temp_min := is_maxima(data, in_lat(l), in_lon(l), opt@rkm)
        end if
        mand_data(l,0,v)  = get_res_value(temp_min, "lat", mand_data@_FillValue)
        mand_data(l,1,v)  = get_res_value(temp_min, "lon", mand_data@_FillValue)
        mand_data(l,2,v)  = get_res_value(temp_min, "val", mand_data@_FillValue)  
      end if
;      if(ismissing(mand_data(l,1,v)))
;      print("missing")
;         check_centerfix(wks, data, (/in_lat, in_lon/) )
;      else
;         check_centerfix(wks, data, mand_data(l,:1,v))
;      end if
    end do
  end do
end if

mand_data!0 = "systems"
mand_data!1= "data"
mand_data!2= "variables"
mand_data&data = (/"lat","lon","val"/)
mand_data&variables = vars


retval = mand_data(:,:,0)
retval = retval@_FillValue

retval(:,:1) = (/ combine_wave_positions(mand_data) /)
nmi = ind(.not.ismissing(retval(:,0)) .and. .not.ismissing(mand_data(:,2,4)) )
if(.not.all(ismissing(nmi)))
  retval(nmi,2) = (/ mand_data(nmi,2,4) /)
end if
retval(:,1) = where(retval(:,1).gt.180, retval(:,1)-360., retval(:,1) ) 

return (retval)
end 







undef("get_starting_targets")
function get_starting_targets(time, lat_bds[2], lon_bds[2], opt)
local min_thresh, cv, uwnd, vwnd, maxes, maxes_mass, easterly, tmp_wave_locs, retval
begin
rkm = get_res_value_keep(opt, "rkm",300)
min_thresh = get_res_value_keep(opt, "min_threshold", 0.1e-5)

cv := grab_data(time, "cv", opt)
cv = circle_avg_m(cv, rkm*1.5)

maxes := loc_max(cv({lat_bds(0):lat_bds(1)},{lon_bds(0):lon_bds(1)}), True, min_thresh)
if(maxes .eq. 0)
  maxes_mass := new( (/1, 3/), float)
  return maxes_mass
end if

maxes_mass := new( (/maxes, 3/), float)
do m=0, maxes-1
 maxes_mass(m,:) =  get_mass_center(cv, (/maxes@lat(m), maxes@lon(m)/) )
end do

;;;;;;;;;;  Remove duplicate locations
maxes_mass := unique_locations(maxes_mass, rkm)
keep_inds = ind(.not.ismissing(maxes_mass(:,0)))
if(dimsizes(keep_inds).eq.1)
 retval := new( (/1,3/), float)
 if(.not.all(ismissing(keep_inds)))
   retval(0,:) = maxes_mass(keep_inds,:)
 end if
else
 retval := maxes_mass(keep_inds,:)        
end if

;;;;;;;;;;;;;;;;;;;;;

;maxes_uwnd := circle_avg_m_point(uwnd, 500, maxes_mass(:,0), maxes_mass(:,1) )
;maxes_vwnd := circle_avg_m_point(vwnd, 500, maxes_mass(:,0), maxes_mass(:,1) )
;easterly := ind(maxes_uwnd .lt. 1)     ;; easterly with a bit of wiggle room
;maxes_mass := maxes_mass(easterly,:)
;if(dimsizes(easterly).eq.1)
; retval := new( (/1,3/), float)
; retval(0,:) = maxes_mass
; maxes_mass := retval
;end if

;;;;;;;;;;;;;   Recenter potential targets on combined wave location  ;;; too slow to do here. 
;tmp_wave_locs := get_wave_positions(time, maxes_mass(:,0), maxes_mass(:,1), 600)
;maxes_mass := combine_mand_positions( tmp_wave_locs )
;maxes_mass := unique_locations(maxes_mass, 250.)

return retval
end


;undef("zScore")
;;;;;;;;;;;;;;  Modified z-score to detect outliers
;;;;;;;;; 
;;Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and Handle Outliers", The ASQC
;;Basic References in Quality Control: Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
;;;;;;;;; 
;;;;;;;;; 
;function zScore(lat, lon)
;local  MAD, Mi, lat, lon, dist, nmi
;begin
;nmi = ind(.not. ismissing(lat))
;if(all(ismissing(nmi))) return default_fillvalue("float") end if
;    
;Mi = new(dimsizes(lat), float, default_fillvalue("float"))
;dist = new(dimsizes(lat), float, default_fillvalue("float"))
;
;med_lat = dim_median(lat(nmi))
;med_lon = dim_median(lon(nmi))
;
;dist(nmi) = gc_latlon(med_lat, med_lon, lat(nmi), lon(nmi),0,4)        ;; absolute deviation from median center
;MAD =  dim_median(dist)  > 5.  <100.           ;; Median Absolute Deviation ; min deviation of 10km so up-to 50km deviation allowed
;Mi = (0.6745*dist)/MAD      ;; > 3.5 Considered outlier. 
;
;;print(dist+"    "+Mi)
;
;return(Mi)
;end


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
    mand_data(1,:) = unwrap_lon(mand_data(1,:), 50.) ;; fix dateline issue.
    zinds := ind( gc_latlon(avg(mand_data(0,:1)), avg(mand_data(1,:1)) , mand_data(0,:5), mand_data(1,:5), 0, 4) .lt. 350)
    if(all(ismissing(zinds)))
	  zinds := ind( gc_latlon( mand_data(0,0), mand_data(1,0) , mand_data(0,:5), mand_data(1,:5), 0, 4) .lt. 350)
 	end if 
    if(dimsizes(zinds).lt.3)
        return( (/mand_data@_FillValue, mand_data@_FillValue/) )
    end if 
    avg_loc = dim_avg_n( mand_data(:1,zinds),1)
    dists := gc_latlon(avg_loc(0), avg_loc(1), mand_data(0,6:7), mand_data(1,6:7), 0,4)
    dinds := ind(dists.le.125)
    if(.not.all(ismissing(dinds)))
        dinds = dinds + 6
        zinds := array_append_record(zinds, dinds, 0)
        avg_loc = dim_avg_n( mand_data(:1,zinds),1)
    end if
    return(avg_loc)
;  end if
  
  
end if
end 


;undef("get_wind_speed")
;function get_wind_speed(in_time, in_lat, in_lon)
;local opt, u, v, uv, temp_min, maxwnd, in_time, in_lat, in_lon
;begin
;    opt := True
;    opt@levs = (/925/)
;    opt@levSelection = "Explicit"
;    uv := grab_data(in_time, "uv", opt)
;    temp_min := is_maxima(uv, in_lat, in_lon, 250.)
;return(temp_min)
;end

undef("check_TC_characteristics")
function check_TC_characteristics(in_time, in_lat, in_lon, opt)
local opt, hgt, hgtd
begin

;;;;;   Thickness Anomaly   ;; is Warm core? ;;;;; ;;;;; ;;;;; 
;    opt := True
;    opt@levs = (/200,850/)
;    opt@levSelection = "Explicit"
;    hgt := grab_data(in_time, "g", opt)
;    hgtd := hgt(0,:,:) - hgt(1,:,:)
;    copy_VarCoords(hgt(0,:,:), hgtd)
;    thick_max := is_env_maxima(hgtd, in_lat,in_lon)
;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; 

;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; 
    data := grab_data(in_time, "mslp", opt)
    locs = (/in_lat, in_lon/)
    temp_min :=  loc_min_dist(data, True, 200, locs); (data, in_lat, in_lon, 200.) 

    loc_data := new( (/ dimsizes(in_lat),7 /), float)
    loc_data(:,0) = in_lat
    loc_data(:,1) = in_lon

    crit := temp_min@val
    crit = crit+100.
    crit = crit > 98000
    crit@loop = True
    crit@interval = 100.
    mslp_closed := check_closed_contour(data, loc_data ,crit)
    loc_data(:,2) = temp_min@val
    loc_data(:,3) = mslp_closed@val
    loc_data(:,4) = mslp_closed@area
    
;    else
;      mslp_closed = False
;    end if
;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; 


;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; 
   g_opt := True
    copy_VarAtts(opt, g_opt)
   g_opt@levs := (/200,850/)
   g_opt@levSelection := "Explicit"
   hgt := grab_data(in_time, "g", g_opt)
   hgtd := hgt(0,:,:) - hgt(1,:,:)
   copy_VarCoords(hgt(0,:,:), hgtd)
   hgt_max := loc_max_dist(hgtd,True, 200, locs)
   loc_data(:,5) =  hgt_max@val

;   vtu := is_env_maxima(hgtd, locs(0),locs(1))
;   loc_data(:,6) =  get_res_value_keep(vtu, "val",loc_data(:,6) )
;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; 



  u := grab_data(in_time, "u", g_opt)
  v := grab_data(in_time, "v", g_opt)
  mag = sqrt( u^2 + v^2)
  totshear = abs(mag(0,:,:) - mag(1,:,:))
  copy_VarCoords_not_n(v,totshear,0)
  loc_data(:,6) = (/circle_avg_m_point(totshear, 500,  in_lat,in_lon)/)


;    opt := True
;    opt@levs = (/300,600/)
;    opt@levSelection = "Explicit"
;    hgt := grab_data(in_time, "g", opt)
;    hgtd := hgt(0,:,:) - hgt(1,:,:)
;    copy_VarCoords(hgt(0,:,:), hgtd)
;    vtu := is_env_maxima(hgtd, in_lat,in_lon)
;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; 

;aux_data = new( 4, float)
;;aux_data(0) = get_res_value(thick_max, "val", aux_data@_FillValue)
;aux_data(1) = get_res_value(mslp_closed, "pres", aux_data@_FillValue)
;;aux_data(2) = get_res_value(vtl, "val", aux_data@_FillValue)
;aux_data(3) = get_res_value(vtu, "val", aux_data@_FillValue)

return(loc_data(:,2:))
end 


;undef("get_aux_variables")
;function get_aux_variables(in_time, in_lat, in_lon)
;local opt, hgt, hgtd
;begin
;
;;;;;;   Thickness Anomaly   ;; is Warm core? ;;;;; ;;;;; ;;;;; 
;;    opt := True
;;    opt@levs = (/200,850/)
;;    opt@levSelection = "Explicit"
;;    hgt := grab_data(in_time, "g", opt)
;;    hgtd := hgt(0,:,:) - hgt(1,:,:)
;;    copy_VarCoords(hgt(0,:,:), hgtd)
;;    thick_max := is_env_maxima(hgtd, in_lat,in_lon)
;;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; 
;
;;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; 
;    data := grab_data(in_time, "pmsl", False)
;    temp_min := is_minima(data, in_lat, in_lon, 200.) 
;    if(temp_min)
;    crit := temp_min@val
;    crit = crit+200.
;    crit@loop = True
;    crit@interval = 200.
;    mslp_closed := check_closed_contour(data, loc_data,crit)
;    else
;      mslp_closed = False
;    end if
;;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; 
;
;
;;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; 
;;    opt := True
;;    opt@levs = (/600,900/)
;;    opt@levSelection = "Explicit"
;;    hgt := grab_data(in_time, "g", opt)
;;    hgtd := hgt(0,:,:) - hgt(1,:,:)
;;    copy_VarCoords(hgt(0,:,:), hgtd)
;;    vtl := is_env_maxima(hgtd, in_lat,in_lon)
;;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; 
;
;;    opt := True
;;    opt@levs = (/300,600/)
;;    opt@levSelection = "Explicit"
;;    hgt := grab_data(in_time, "g", opt)
;;    hgtd := hgt(0,:,:) - hgt(1,:,:)
;;    copy_VarCoords(hgt(0,:,:), hgtd)
;;    vtu := is_env_maxima(hgtd, in_lat,in_lon)
;;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; ;;;;; 
;
;aux_data = new( 4, float)
;;aux_data(0) = get_res_value(thick_max, "val", aux_data@_FillValue)
;aux_data(1) = get_res_value(mslp_closed, "pres", aux_data@_FillValue)
;;aux_data(2) = get_res_value(vtl, "val", aux_data@_FillValue)
;;aux_data(3) = get_res_value(vtu, "val", aux_data@_FillValue)
;
;return(aux_data)
;end 

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




undef("nudge_lons")
function nudge_lons( lons, obs_lon )
begin
    if( obs_lon - lons(0) .lt. -90 )
        obs_lon =obs_lon+360.
    end if
    if( obs_lon - lons(0) .gt. 90 )
        obs_lon =obs_lon-360.
    end if
lons(1:) = lons(1:) + (obs_lon - lons(0))/2.
return wrap_lon(lons)
end 



undef("nudge_lats")
function nudge_lats( lats, obs_lat )
begin
lats(1:) = lats(1:) + (obs_lat - lats(0))/2.
return (lats)
end 



undef("filter_tracks")
procedure filter_tracks(track_locs,t)
local n_miss, n_track, kill
begin
    ;; optional procedure to kill tracks that are too weak to track \
    ;; filters based on number of missing vorticity values within last 8 track positions
    ;; if track is longer than 4 times and 50% are missing get rid of it
  n_miss :=  dim_num_n(.not.ismissing(track_locs(:,0,(t-8)>0:t)) .and. ismissing(track_locs(:,2,(t-8)>0:t)),1 ) 
  n_track := dim_num_n(.not.ismissing(track_locs(:,0,(t-8)>0:t)),1) 
  kill := ind( n_track.gt.4 .and. n_miss .gt. (n_track/2.) )
  if(.not.all(ismissing(kill)))
   print("kill kill")
   track_locs(kill,:,t) =  track_locs@_FillValue
  end if
end


undef("advect_tracks")
function advect_tracks(time, track_loc)
local nmi, rkm,time,track_loc, steering, uwnd, vwnd, wave_ups, wave_vps
begin
  nmi := ind(.not.ismissing(track_loc(:,0)))
  retvals = new( dimsizes(track_loc(:,:1)), float)
  if(all(ismissing(nmi)))  return   end if
  dt = get_res_value_keep(time, "dt", 6)
  rkm = 350
  steering := True
  copy_VarAtts(time, steering)
  steering@levs := (/800, 300/)
  steering@levSelection = "Range"
  uwnd := dim_avg_n_Wrap(grab_data(time, "u", steering),0)
  vwnd := dim_avg_n_Wrap(grab_data(time, "v", steering),0)
  wave_ups := circle_avg_m_point(uwnd, rkm, track_loc(nmi,0), track_loc(nmi,1) )
  wave_vps := circle_avg_m_point(vwnd, rkm, track_loc(nmi,0), track_loc(nmi,1) )
  retvals(nmi,:) = advect_location(track_loc(nmi,:1), wave_ups, wave_vps, 60.*60.*dt)
 return(  retvals  )
end



undef("is_closed_circulation")
function is_closed_circulation( u, v, lat, lon)
local DEPT_THRESH, PI \
    ,radii, nradii,center, circles, u_circle, v_circle, wind_angles, grad, dept \
    , max_wspd, retval
begin
;  '''  Grab u and v wind components from circles around center.
;       Calculate wind angles and the gradient of angles around circle
;       Sum of the angles should = 360
;       if any single gradient is greater than threshold, set sum to that value
;       Can play with that threshold later.
load "~/NCL/functions/circles.ncl"
load "~/NCL/functions/interpolation/interp.ncl"


    DEPT_THRESH = 5.
    PI = 355/113.

    radii = ispan(50, 250, 50)
    nradii = dimsizes(radii)
    center = (/lat, lon/)

    circles := new( (/nradii,2,360/), float)
    do r=0,nradii-1
      circles(r,:,:) = circle_km(center(0), center(1), radii(r) )
    end do

    u_circle := bicubic_interp(u,  circles(:,0,:), circles(:,1,:) )
    v_circle := bicubic_interp(v,  circles(:,0,:), circles(:,1,:) )

    wind_angles := atan2( u_circle, v_circle)
    grad := center_finite_diff_n(todeg(wind_angles),1, True,0,1)
    grad := abs(where(abs(180-grad).lt.abs(grad), 180-grad, grad))
    dept = dim_sum_n(grad,1)
    dept = abs(where( abs(360-dept).lt.dept,360-dept , dept))
    dept = where( dim_max_n( grad ,1 ) .gt. DEPT_THRESH, dim_max_n( grad ,1 ) , dept)

    retval = new( (/360, 4/), float)
    if( num(dept .le. DEPT_THRESH) .ge. 2 )
        max_wspd =  dim_max_n( sqrt( u_circle^2+v_circle^2),1)
        max_wspd = where(dept .gt. DEPT_THRESH, 0, max_wspd)
        retind = maxind(max_wspd)
;        retind = max( ind( dept .le. DEPT_THRESH))
        retval(:,0) = circles(retind,0,:)
        retval(:,1) = circles(retind,1,:)
        retval(:,2) = u_circle(retind,:)
        retval(:,3) = v_circle(retind,:)
        retval@radius = radii(retind)
        retval@area   = PI*radii(retind)^2
        retval@departure = dept(retind)
    else
        retval@area = 0
    end if
return retval
end




