load "~/NCL/functions/remove_loops/remove_loops.ncl"


undef("remove_line_loops")
function remove_line_loops(cone_in)
local nc, ii, md, nmin, cone, cone_in
begin
cone = cone_in
nc = dimsizes(cone(0,:))-1
do ii=0,nc-2
  md := sqrt( (cone(0,ii) - cone(0,ii+2:nc))^2 + (cone(1,ii) - cone(1,ii+2:nc))^2) 
  nmin := local_min_1d(md, False, 0.,1)
  if( all(.not.ismissing(nmin)) .and. any(md(nmin).lt.0.5) )
    nmin := nmin(minind(md(nmin)))
    if( (nmin .lt. nc/5) )
      cone(0,ii:(ii+nmin)) = cone@_FillValue
      nmi := ind(.not.ismissing(cone(0,:)))
      cone := cone(:, nmi)
    end if
  end if
  nc = dimsizes(cone(0,:))- 1
end do

return(cone)
end



undef("unwrap_degrees")
;;;;  Check for discontinuities in array greater than thresh
;;;;  unwrap array so it's monotonic
function unwrap_degrees(lon, thresh)
local lon, thresh, nl, nmi, d
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

if(all(lon .gt. 360 ) )
    lon = lon - 360
end if
if(all(lon .lt. -360 ) )
    lon = lon + 360
end if

return(lon)
end

undef("unwrap_radians")
;;;;  Check for discontinuities in array greater than thresh
;;;;  unwrap array so it's monotonic
function unwrap_radians(lon, thresh)
local lon, thresh, nl, nmi, d
begin
pi   = 355/113.

nl = dimsizes(lon)
nmi := ind(.not.ismissing(lon))

do d=0,dimsizes(nmi)-2

  if( lon(nmi(d+1)) - lon(nmi(d)).lt. thresh)
    lon(nmi(d+1))  = lon(nmi(d+1))+(2*pi)
  end if
  if(lon(nmi(d+1)) - lon(nmi(d)).gt. thresh)
    lon(nmi(d+1))  = lon(nmi(d+1))-(2*pi)
  end if
end do

if(all(lon .gt. (pi) ) )
    lon = lon - (2*pi)
end if
if(all(lon .lt. -(pi) ) )
    lon = lon + (2*pi)
end if

return(lon)
end



undef("along_across_error")
function along_across_error(obs_lat, obs_lon, fcst_lat, fcst_lon)
local pi, pi2d,deg2km,  fxdir, fydir, dbx, dby, vlen, alnorm, axnorm, faldir, faxdir\
    , obs_lat, obs_lon, fcst_lat, fcst_lon
begin
pi   = 355/113.
p2d = (pi/180.)
deg2km = 111.

dimz = dimsizes(obs_lat)
rankz = dimsizes(dimz)

fdimz = dimsizes(fcst_lat)
frankz = dimsizes(fdimz)

if(rankz.eq.1 .and. (frankz.eq.2) .and. (dimz(0).eq.fdimz(1)) )
tobs_lat = conform(fcst_lat, obs_lat, 1)
tobs_lon = conform(fcst_lon, obs_lon, 1)
else
tobs_lat = obs_lat
tobs_lon = obs_lon
end if
dimz := dimsizes(tobs_lat)
rankz := dimsizes(dimz)

fxdir = (fcst_lon - tobs_lon) * deg2km * cos(0.5*(tobs_lat+fcst_lat)*p2d)
fydir = (fcst_lat - tobs_lat) * deg2km

dbx = center_finite_diff_n(tobs_lon, 1, False, 0, rankz-1)* deg2km * cos(0.5*(tobs_lat)*p2d)
dby = center_finite_diff_n(tobs_lat, 1, False, 0, rankz-1)* deg2km 

;dbx = (tobs_lon(t+1) - tobs_lon(t-1)) * deg2km * cos(0.5*(tobs_lat(t)+fcst_lat(1))*p2d)
;dby = (tobs_lat(t+1) - tobs_lat(t-1)) * deg2km

vlen = sqrt( dbx * dbx + dby * dby ) > 0.0001

;alnorm = (/ dbx / vlen  ,  dby / vlen /)

alnormx  =  dbx / vlen
alnormy  =  dby / vlen 

;axnorm = (/ alnorm(1,:)  , -alnorm(0,:) /) 

faldir = alnormx* fxdir + alnormy * fydir
faxdir = alnormy * fxdir + -alnormx * fydir

return( (/faldir, faxdir/) )
end


undef("track_dir_dist")
function track_dir_dist(tobs_lat, tobs_lon)
local pi, p2d, deg2km, dimz, rankz, dbx, dby, dir, vlen
begin
pi   = 355/113.
p2d = (pi/180.)
deg2km = 111.

dimz := dimsizes(tobs_lat)
rankz := dimsizes(dimz)

dbx = center_finite_diff_n(tobs_lon, 1, False, 0, rankz-1)* deg2km * cos(0.5*(tobs_lat)*p2d)
dby = center_finite_diff_n(tobs_lat, 1, False, 0, rankz-1)* deg2km 
dir = atan2(dby, dbx)
vlen = sqrt( dbx * dbx + dby * dby ) > 0.0001

return( (/dir, vlen/) )
end

undef("advect_by_km_radians")
function advect_by_km_radians(s_lat,s_lon, dist, angle)
local pi, p2d, deg2km, s_lat, s_lon, dist, angle, fxdist, fydist, fxdeg, fydeg, r_lat, r_lon
begin
;           
;-pi | pi   +     0
;         -pi/2 
pi   = 355/113.
p2d = (pi/180.)
deg2km = 111.

fxdist = cos(angle)*dist
fydist = sin(angle)*dist

fxdeg = fxdist / deg2km
fydeg = fydist / (deg2km*cos(0.5*(s_lat)*p2d) )

r_lat = s_lat + fydeg
r_lon = s_lon + fxdeg
return( (/r_lat, r_lon/) )
end 



undef("make_cone")
function make_cone(tcdata, rmaxhr) 
local nens, nmi, mean_locs, wgt, ti, mean_lat, mean_lon, mean_track, ens_spread,\
     lcone, rcone, end_cone, start_cone, cone_lats, cone_lons, cone, t,r, maxhr, rmaxhr
begin
cone = new( (/2, 1001/), float)

nens = get_res_value_keep(tcdata, "nens", 21)

if(any(getvardims(tcdata).eq."tau"))
  tvar = "tau"
else
  tvar="time"
end if
;ct = dim_num_n(ismissing(tcdata(0,:,:)),1)
 wgt := tcdata(0,:,0)
 wgt = 1
 
 do e=0, dimsizes(tcdata(1,:,0))-1
   tcdata(1,e,:) = unwrap_degrees(tcdata(1,e,:), 20)
end do
 
 mean_locs := dim_avg_wgt_n_Wrap(tcdata, wgt, nens/3, 1)

 ;wgt := gc_latlon(mean_locs(0,{72}),mean_locs(1,{72}), tcdata(0,:,{72}),tcdata(1,:,{72}),0,4) > 25.
 ;wgt@_FillValue = default_fillvalue(typeof(wgt))
 ;wgt = where(ismissing(tcdata(1,:,{72})), wgt@_FillValue, wgt)
 ;ti = ind_nearest_coord(72, mean_locs&$tvar$, 0)
 ;do t=ti+1, dimsizes(mean_locs&$tvar$)-1
 ;  mean_locs(:,t) = dim_avg_wgt_n(tcdata(:,:,t), 1/wgt, nens/3, 1)
 ;end do
; mean_locs_hr = linint1_n_Wrap(mean_locs&time, mean_locs, False, fspan(min(mean_locs&time), max(mean_locs&time), 100), 0, 1)
 mean_locs = runave_n(mean_locs, 4, 1, 1)
 mean_locs = runave_n(mean_locs, 4, 1, 1)
; mean_locs = (/ mean_locs_hr(:,{mean_locs&time}) /)

 nmi = max(ind(.not.ismissing(mean_locs(0,:))))
 if(ismissing(nmi)) return cone end if
 maxhr = totype(min( (/rmaxhr, mean_locs&$tvar$(nmi)/) ), typeof(rmaxhr) )
 mean_lat := mean_locs(0,{:maxhr})
 mean_lon := mean_locs(1,{:maxhr})
;print( mean_lat+"  "+mean_lon)

 mean_lat =  cubicmsg(mean_lat, False)
 mean_lon =  cubicmsg(mean_lon, False)
if(dimsizes(mean_lat&$tvar$).lt.6) return cone end if 

 mean_track := track_dir_dist(mean_lat({:maxhr}), mean_lon({:maxhr}))
 mean_track(0,:) = unwrap_radians(mean_track(0,:), pi/3.)

 mean_track = runave_n(mean_track,3,1,1)
 mean_track@_FillValue = 0
 mean_track = linmsg_n(mean_track,-1, 1)
; mean_track(0,:) = where(mean_track(0,:) .gt. (2*pi+ pi/2.), mean_track(0,:) - pi, mean_track(0,:) )
 
 
 errors := along_across_error(mean_lat, mean_lon, tcdata(0,:,{mean_lat&$tvar$}), tcdata(1,:,{mean_lat&$tvar$}))
 ens_spread := dim_stddev_n(errors(1,:,:),0)
 ens_spread = runave_n(ens_spread,3,1,0)

 do t=0, dimsizes(ens_spread)-1
  ens_spread(t) = dim_max(ens_spread(:t))
 end do
 ens_spread = tofloat( ens_spread  + (mean_lon&$tvar$) ) > 50. < 750.
;    print(ens_spread)

 lcone := new( (/2,dimsizes(ens_spread)/), float)
 rcone := new((/2,dimsizes(ens_spread)/), float)
 do t=0, dimsizes(ens_spread)-1
  if(.not. ismissing(mean_lat(t)) .and. .not. ismissing(mean_track(0,t)) )
;    print(  mean_track(0,t) +","+ (mean_track(0,t)-1.57 ))
     lcone(:,t) = advect_by_km_radians(mean_lat(t), mean_lon(t), ens_spread(t), mean_track(0,t)-1.57)
     rcone(:,t) = advect_by_km_radians(mean_lat(t), mean_lon(t), ens_spread(t), mean_track(0,t)+1.57)
   end if
 end do
;exit
end_cone:= new( (/2, 35/), float)
start_cone:= new( (/2, 35/), float)
t= dimsizes(ens_spread)-1
do r=0,34
     rad = -1.55 + (pi/35.)*r
    end_cone(:,r) = advect_by_km_radians(mean_lat(t), mean_lon(t), ens_spread(t), mean_track(0,t)+rad)
    start_cone(:,r) = advect_by_km_radians(mean_lat(0), mean_lon(0), ens_spread(0), mean_track(0,0)+(pi-rad))
end do

cone = new( (/2, 1001/), float)
lasti = max(ind(.not.ismissing(rcone(0,:))))
if(ismissing(lasti))
    return cone
end if
cone_dist = (/ gc_latlon( end_cone(0,0), end_cone(1,0),rcone(0,lasti), rcone(1,lasti),0,4), \
                gc_latlon( end_cone(0,34), end_cone(1,34),rcone(0,lasti), rcone(1,lasti),0,4) /)
if(cone_dist(0) .gt. cone_dist(1))
    end_cone = end_cone(:,::-1)
end if
cone_dist = (/ gc_latlon( start_cone(0,0), start_cone(1,0),lcone(0,0), lcone(1,0),0,4), \
                gc_latlon( start_cone(0,34), start_cone(1,34),lcone(0,0),lcone(1,0),0,4) /)
if(cone_dist(0) .gt. cone_dist(1))
    start_cone = start_cone(:,::-1)
end if




cone_lats := array_append_record( rcone(0,:), end_cone(0,:),0)
cone_lons := array_append_record( rcone(1,:), end_cone(1,:),0)
cone_lats := array_append_record( cone_lats, lcone(0,::-1),0)
cone_lons := array_append_record(cone_lons, lcone(1,::-1),0)
cone_lats := array_append_record( cone_lats, start_cone(0,:),0)
cone_lons := array_append_record(cone_lons, start_cone(1,:),0)
cone_lats := array_append_record( cone_lats, cone_lats(0),0)
cone_lons := array_append_record(cone_lons, cone_lons(0),0)

nsmth = 10
cone(0,:) = runave( linint1(fspan(0,100, dimsizes(cone_lats)), cone_lats, True, fspan(0,100,1001), 0), nsmth, -1)
cone(1,:) = runave( linint1(fspan(0,100, dimsizes(cone_lats)), cone_lons, True, fspan(0,100,1001), 0), nsmth, -1)

coned :=    remove_loops(cone) 
npts = dimsizes(coned(0,:))
cone(0,:) = runave( linint1(fspan(0,100, npts), coned(0,:), True, fspan(0,100,1001), 0), nsmth, -1)
cone(1,:) = runave( linint1(fspan(0,100, npts), coned(1,:), True, fspan(0,100,1001), 0), nsmth, -1)


return cone
end 



undef("plot_cone")
procedure plot_cone(wks, map, xx, yy, resell)
local my, mx, ellcol, resell, xell, yell, Pb, rho, sigmax, sigmay, fac, rdex, radians, rdistance, xloc, yloc, dist, prob ,x, y, c, res,t, resell
begin

mean_lat := dim_avg_n(yy, 0)
mean_lon := dim_avg_n(xx, 0)
mean_track := track_dir_dist(mean_lat, mean_lon)
mean_track = runave_n(mean_track, 4, 1, 1)
mean_track(0,:) = where(mean_track(0,:).lt.0, mean_track(0,:)+(2*pi), mean_track(0,:))
dims = dimsizes(xx)

cone_edges = new( (/2,2, dims(1)/), float)
cone_lats := new((dims(1)*2)+38, float)
cone_lons := new((dims(1)*2)+38, float)

do t=0, dims(1)-1
  x := xx(:,t)
  y := yy(:,t)

  my = mean_lat(t)
  mx = mean_lon(t)
    if(all(ismissing(x))) continue end if

  Pb= new( (/2,2/), float)
  Pb(:,:) = 0.0

  Pb(0,0) = dim_sum( (x(:)-mx) * (x(:)-mx) )
  Pb(1,1) = dim_sum( (y(:)-my) * (y(:)-my) )
  Pb(1,0) = dim_sum( (y(:)-my) * (x(:)-mx) )

  Pb(0,1) = Pb(1,0)

  if( Pb(0,0) .eq.0 .or. Pb(1,1) .eq. 0)
    continue
  end if
  Pb(:,:) = Pb(:,:) / int2flt(num(.not.ismissing(y))-1)
  rho     = Pb(1,0) / (sqrt(Pb(0,0)) * sqrt(Pb(1,1)))
  if(abs(rho) .eq. 1) 
     rho = rho * 0.95
;     continue 
  end if
  sigmax  = sqrt(Pb(0,0))
  sigmay  = sqrt(Pb(1,1))
  fac     = 1.0 / (2.0 * (1 - rho * rho))

  xell := new(21, float)
  yell := new(21, float)
  rdistance =ispan(200,1000,1)
  do  rdex = 0, 1
    radians  = mean_track(0,t)+((rdex*3>1)*pi/2.)
    xstart = cos(radians)
    ystart = sin(radians)
    xloc = xstart*rdistance/80.0
    yloc = ystart*rdistance/80.0
    dist = sqrt(xloc * xloc + yloc * yloc)
    prob = exp(-1.0 * fac * ((xloc/sigmax)^2 + (yloc/sigmay)^2 - 2.0 * rho * (xloc/sigmax)*(yloc/sigmay)))
    prob_i = min(ind(prob.lt.0.256))
    if(all(ismissing(prob_i)))
         prob_i = dimsizes(xloc)-1
    end if
    if(rdex.eq.0)
       cone_lats(20+t) =  yloc(prob_i) + my
       cone_lons(20+t) =  xloc(prob_i) + mx
    end if
    if(rdex.eq.1)
       cone_lats(dimsizes(cone_lats)-1-t) =  yloc(prob_i) + my
       cone_lons(dimsizes(cone_lats)-1-t) =  xloc(prob_i) + mx
    end if
  end do

  if(t.eq.dims(1)-1)
    rdex = 0
    do radians  = mean_track(0,t)+(pi/2.), mean_track(0,t)+(3*pi/2.), pi/21.
      xstart = cos(radians)
      ystart = sin(radians)
      xloc = xstart*rdistance/80.0
      yloc = ystart*rdistance/80.0
      dist = sqrt(xloc * xloc + yloc * yloc)
      prob = exp(-1.0 * fac * ((xloc/sigmax)^2 + (yloc/sigmay)^2 - 2.0 * rho * (xloc/sigmax)*(yloc/sigmay)))
      prob_i = min(ind(prob.lt.0.256))
      if(all(ismissing(prob_i)))
           prob_i = dimsizes(xloc)-1
      end if
      cone_lons(20+dims(1)+rdex) = xloc(prob_i) + mx
      cone_lats(20+dims(1)+rdex) = yloc(prob_i) + my
      rdex = rdex +1
    end do
    end if
  if(t.eq.0)
    rdex = 0
    do radians  = mean_track(0,t)-(pi/2.) , mean_track(0,t)+(pi/2.), pi/21.
      xstart = cos(radians)
      ystart = sin(radians)
      xloc = xstart*rdistance/80.0
      yloc = ystart*rdistance/80.0
      dist = sqrt(xloc * xloc + yloc * yloc)
      prob = exp(-1.0 * fac * ((xloc/sigmax)^2 + (yloc/sigmay)^2 - 2.0 * rho * (xloc/sigmax)*(yloc/sigmay)))
      prob_i = min(ind(prob.lt.0.256))
      if(all(ismissing(prob_i)))
           prob_i = dimsizes(xloc)-1
      end if
      cone_lons(rdex) = xloc(prob_i) + mx
      cone_lats(rdex) = yloc(prob_i) + my
      rdex = rdex +1
    end do
    end if


end do

cts = str_sub_str(str_match(getvaratts(map),"cone"),"cone","" )
if(all(ismissing(cts)))
   ct = 0
else
   ct = max(toint(cts))+1
end if 
;
restemp = resell
if(.not.isatt(restemp,"gsLineThicknessF"))
 restemp@gsLineThicknessF := 1.
end if
restemp@gsLineThicknessF := 1.15*  restemp@gsLineThicknessF
;
;do i=0, iso@segment_count-1
;    if(iso@n_points(i).lt. 20) continue end if
;    j = iso@start_point(i)
;    k = iso@start_point(i)+iso@n_points(i)-1
nlons = dimsizes(cone_lons)
cone_lons(nlons-1) = cone_lons(0)
cone_lats(nlons-1) = cone_lats(0)
cone_lons = linmsg(cone_lons, -1)
cone_lats = linmsg(cone_lats, -1)
;print(any(ismissing(cone_lons)))
smth = toint(dimsizes(cone_lons)*0.15)
cone_lons := runave(cone_lons, smth, -1)
cone_lats := runave(cone_lats, smth, -1)
cone_lons(nlons-1) = cone_lons(0)
cone_lats(nlons-1) = cone_lats(0)
;print(cone_lons(nlons-1)+"   "+cone_lons(0))
map@$tostring("cone"+ct)$ = gsn_add_polyline(wks,map,  cone_lons, cone_lats, restemp)
;end do
end

 undef("add_traj_line")
 procedure add_traj_line( wks, map, Lat, Lon, Value, levels, colors,gsres)
 local nmi, lat, lon, value, newlon, newlat, newval, lres
 begin
 nmi = ind( (.not. ismissing(Lat)) .and. (.not. ismissing(Lon)) .and. (.not. ismissing(Value)) )
 if(all(ismissing(nmi))) return end if
 lon =  Lon(nmi)
 lat = Lat(nmi)
 value = Value(nmi)

 newlon = new(dimsizes(lon)*2, typeof(lon))
 newlon(::2) = lon
 newlon(1::2) = lon
 newlat = new(dimsizes(lat)*2, typeof(lat))
 newlat(::2) = lat
 newlat(1::2) = lat
 newval = new(dimsizes(value)*2, typeof(value))
 newval(::2) = value
 newval(1::2) = value

 lres = True
 lres@gsLineThicknessF = get_res_value_keep(gsres,"gsLineThicknessF", 2.)
 lres@gsSegments = ispan(1, dimsizes(newlat)-2,2)
 lres@gsEdgesOn   = get_res_value_keep(gsres, "gsEdgesOn", True )   ; Turn on edges
 lres@gsEdgeColor = get_res_value_keep(gsres, "gsEdgeColor", "white")
 lres@gsFillIndex = get_res_value_keep(gsres, "gsFillIndex", 0)      ; Solid fill, the default


 icolors = get_color_indexes( colors, levels, newval)
; print(icolors)
 if(isnumeric(colors))
     cmap = colors
 else
     cmap = read_colormap_file(colors)
 end if
;  printVarSummary(cmap)
 lres@gsColors = cmap(icolors,:)

 cts = str_sub_str(str_match(getvaratts(map),"traj"),"traj","" )
 if(all(ismissing(cts)))
     ct = 0
 else
     ct = max(stringtoint(cts))+1
 end if 

 map@$tostring("traj"+ct)$ = gsn_add_polyline(wks, map, newlon,  newlat, lres)

 end


undef("which_basin")
function which_basin(lat, lon)
local basin_lats, basin_lons, lat, lon,basins
begin
basins = new( dimsizes(lat), "string", "un")

lat1d = ndtooned(lat)
lon1d = ndtooned(lon)
basin1d = ndtooned(basins)

basin_lats :=  (/0,0,60,60, 20, 15,12,10,8.3,9.25,0/)
basin_lons := (/-72,-140,-140,-100,-100, -90,-85,-84,-81.5,-79,-72/)
nmi := ind( ismissing(basin1d) .and. .not. ismissing(lat1d) )
if(any(.not.ismissing(nmi)))
basin1d(nmi) =  where( gc_inout( lat1d(nmi), lon1d(nmi), basin_lats, basin_lons), "ep", basin1d(nmi) )
end if
;basins =  where( gc_inout( lat, lon, basin_lats, basin_lons), "ep", basins )

basin_lons :=  (/-72,0,0,-100,-100, -90,-85,-84,-81.5,-79,-72/)
basin_lats :=  (/0,0,60,60, 20, 15,12,10,8.3,9.25,0/)
;basins =  where( gc_inout( lat, lon, basin_lats, basin_lons), "al", basins )
nmi :=  ind( ismissing(basin1d) .and. .not. ismissing(lat1d) )
if(any(.not.ismissing(nmi)))
    basin1d(nmi) =  where( gc_inout( lat1d(nmi), lon1d(nmi), basin_lats, basin_lons), "al", basin1d(nmi) )
end if

basin_lons := (/180, 180, 220, 220, 180/)
basin_lats := (/0,60,60,0,0/)
;basins =  where( gc_inout( lat, lon, basin_lats, basin_lons), "cp", basins )
nmi := ind( ismissing(basin1d) .and. .not. ismissing(lat1d) )
if(any(.not.ismissing(nmi)))
basin1d(nmi) =  where( gc_inout( lat1d(nmi), lon1d(nmi), basin_lats, basin_lons), "cp", basin1d(nmi) )
end if

basin_lons := (/180, 180, 100, 100, 180/)
basin_lats := (/0,60,60,0,0/)
;basins =  where( gc_inout( lat, lon, basin_lats, basin_lons), "wp", basins )
nmi :=  ind( ismissing(basin1d) .and. .not. ismissing(lat1d) )
if(any(.not.ismissing(nmi)))
basin1d(nmi) =  where( gc_inout( lat1d(nmi), lon1d(nmi), basin_lats, basin_lons), "wp", basin1d(nmi) )
end if


basin_lons := (/40, 40, 100, 100, 40/)
basin_lats := (/0,60,60,0,0/)
;basins =  where( gc_inout( lat, lon, basin_lats, basin_lons), "io", basins )
nmi :=  ind( ismissing(basin1d) .and. .not. ismissing(lat1d) )
if(any(.not.ismissing(nmi)))
basin1d(nmi) =  where( gc_inout( lat1d(nmi), lon1d(nmi), basin_lats, basin_lons), "io", basin1d(nmi) )
end if


basin_lons := (/20, 20, 120, 120, 20/)
basin_lats := (/0,-60,-60,0,0/)
;basins =  where( gc_inout( lat, lon, basin_lats, basin_lons), "io", basins )
nmi :=  ind( ismissing(basin1d) .and. .not. ismissing(lat1d) )
if(any(.not.ismissing(nmi)))
basin1d(nmi) =  where( gc_inout( lat1d(nmi), lon1d(nmi), basin_lats, basin_lons), "io", basin1d(nmi) )
end if

basins = onedtond(basin1d, dimsizes(basins))

return basins
end 



undef("directional_track_error")
function directional_track_error(inlat, inlon, clat, clon, angle)
local alat, alon, rangle, major_dlat, major_dlon, a, retval
begin
    alat = inlat - clat
    alon = inlon - clon
    rangle = torad(angle)

    major_dlat = torad( alat )
    major_dlon = torad( alon )
    
    a = sin(major_dlat/2.)^2 + (cos(avg(torad(inlat)))^2 * sin( major_dlon/2.)^2)
    retval= 6371. *2.* atan2( sqrt(a), sqrt(1-a))

    major_rangle := atan2(major_dlat, major_dlon)
    retval = retval * cos( rangle -major_rangle )
    return retval
end 


undef("calculate_ellipse_points")
function calculate_ellipse_points( clat, clon, major, minor, rotation)
local th, pi, scale, rot, xpts1, ypts1, xpts2, ypts2, ypts, xpts
begin
  th = ispan(0,360,1)
  pi = acos(-1.0)
  rot = rotation

  xpts1 = major*cos(th*pi/180)
  ypts1 = minor*sin(th*pi/180)

  xpts2 = xpts1*cos(rot*pi/180) - ypts1*sin(rot*pi/180)
  ypts2 = xpts1*sin(rot*pi/180) + ypts1*cos(rot*pi/180)

  ypts =ypts2 + clat
  xpts =xpts2 + clon

return( (/ypts, xpts/) )
end

undef("calculate_ellipse_axis")
function calculate_ellipse_axis( in_lat, in_lon)
local nmi, lat, lon, eof, rangle, angle, retval, ell_points, alon, alat\
 , major_dlat, major_dlon, minor_dlat, minor_dlon,a
begin
    pi = 355/113.

    nmi := ind( .not. ismissing( in_lat))
    if(all(ismissing(nmi))) return new(3, "float", default_fillvalue("float"))  end if
    lat = in_lat(nmi)
    lon = in_lon(nmi)
    if(max(lat).eq.min(lat)) 
       lat = lat + random_normal(0, 0.15, dimsizes(nmi))
    end if
    if(max(lon).eq.min(lon)) 
       lon = lon + random_normal(0, 0.15, dimsizes(nmi))
    end if
    eof = eofunc((/lat, lon/),2, False)
    rangle = atan2(eof(0,0), eof(0,1))    
    if(all(ismissing(rangle))) return new(3, "float", default_fillvalue("float"))  end if
    if(rangle .lt. 0 )
         rangle = rangle + pi
    end if

    angle = todeg(rangle)
    eigen1 = sqrt(eof@eval(0))
    eigen2 = sqrt(eof@eval(1))
    retval = (/eigen1, eigen2, angle/)

    retval@mlat = avg(lat)
    retval@mlon = avg(lon)


    xpts1 =  torad(1*eigen1*cos(rangle)) ;+ retval@mlon
    ypts1 =  torad(1*eigen1*sin(rangle)) ;+ retval@mlat

    xpts2 =  torad(1*eigen2*cos(rangle+pi/2.)) ;+ retval@mlon
    ypts2 =  torad(1*eigen2*sin(rangle+pi/2.)) ;+ retval@mlat

    eigen_a = sin(ypts1/2.)^2 + (cos(avg(torad(lat)))^2 * sin( xpts1/2.)^2)
    eigen1_km =  6371. * 2.* atan2( sqrt(eigen_a), sqrt(1-eigen_a))

    eigen_a = sin(ypts2/2.)^2 + (cos(avg(torad(lat)))^2 * sin( xpts2/2.)^2)
    eigen2_km =  6371. * 2.* atan2( sqrt(eigen_a), sqrt(1-eigen_a))
    retval@eigen_km= (/ eigen1_km,  eigen2_km /)


    ell_points = calculate_ellipse_points( avg(lat), avg(lon), retval(0),  retval(1), retval(2) )
    retval@ellipse_2x = ell_points(1,:)
    retval@ellipse_2y = ell_points(0,:)

    alon := lon - avg(lon)
    alat := lat - avg(lat)

    retval@major = alon*cos( rangle )   + alat*sin( rangle )
    retval@minor = alon*cos( rangle+pi/2. )+ alat*sin( rangle+pi/2. )

    retval@major_km = directional_track_error(lat, lon, avg(lat), avg(lon), angle)
    retval@minor_km = directional_track_error(lat, lon, avg(lat), avg(lon), angle-90)

    retval@tech = in_lat&tech(nmi)

return ( retval )
end


undef("calculate_along_axis")
function calculate_along_axis(in_lat, in_lon)
local pi, pi2d,deg2km,  mlat, mlon, dby, dbx, rangle, angle, in_lat, in_lon
begin

    pi   = 355/113.
    p2d = (pi/180.)
    deg2km = 111.

    mlat = dim_avg_n_Wrap(in_lat,0)
    mlon = dim_avg_n_Wrap(in_lon,0)
    dby = center_finite_diff(mlat, 1, False, 0)* deg2km     
    dbx = center_finite_diff(mlon, 1, False, 0)* deg2km * cos(0.5*(mlat)*p2d)
    rangle = atan2(dby, dbx)
    if(all(ismissing(rangle))) return new(dimsizes(in_lat), "float", default_fillvalue("float"))  end if
    rangle = where(rangle.lt.0, rangle+pi, rangle)
    angle = todeg(rangle)
    copy_VarCoords_not_n(in_lat, angle,0)
    return( angle )
end



undef("plot_elipse")
function plot_elipse(wks, plot, x, y, resell, restext)
local my, mx, ellcol, resell, xell, yell, Pb, rho, sigmax, sigmay, fac, rdex, radians, rdistance, xloc, yloc, dist, prob ,x, y, c, res
begin
ell_axis = calculate_ellipse_axis(y, x)

xell = ell_axis@ellipse_2x
yell = ell_axis@ellipse_2y
mx = ell_axis@mlon
my = ell_axis@mlat
dums = new(3, graphic)
if(all(ismissing(xell)))
return dums
end if

restemp = resell
if(.not.isatt(restemp,"gsLineThicknessF"))
  restemp@gsLineThicknessF = 1.
end if
restemp@gsLineThicknessF = 1.15*  restemp@gsLineThicknessF
restemp@gsLineColor = "White"

dums(0) = gsn_add_polyline(wks,plot,xell,yell,restemp)
dums(1) = gsn_add_polyline(wks,plot,xell,yell,resell)
if(restext .and. isatt(restext, "text") )
 dums(2) = gsn_add_text(wks,plot,restext@text,mx,my,restext)
end if
return(dums)
end



 
undef("data_to_adeck")
function data_to_adeck( id,basin, tcdata)
local id, basin, tcdata, lats, lons, rlons, subsets, vmax, mslp, time, storm_id, date, date_num, year\
    , lat_s, lon_s, vmax_s,mslp_s, type, id_s, outtext, delim, smin
begin
  lats = tcdata(0,:,:)
  lons = tcdata(1,:,:)
  rlons = lons
  subsets = dimsizes(tcdata&ens_mem)
  do i=0, subsets-1
    lons(i,:) = wrap_lon(lons(i,:))
  end do
  vmax = tcdata(2,:,:)
  mslp = tcdata(3,:,:)
  time = tcdata&time
  storm_id = id
  date := tcdata@initial_time
  date_num = cd_inv_string(date, "%Y%N%D")
  year = cd_string(date_num, "%Y")

  lats@_FillValue = default_fillvalue("integer")
  lons@_FillValue = default_fillvalue("integer")
  
;  print(tcdata@tech)
  lat_s := where(lats.lt.0, sprinti("%3i",toint( abs(lats)*10 ) )+"S", sprinti("%3i",toint( abs(lats)*10  ))+"N")
  lon_s := where(lons.lt.0, sprinti("%4i",toint( abs(lons)*10 ) )+"W", sprinti("%4i",toint( abs(lons)*10 ) )+"E")
  vmax_s := where(ismissing(vmax), "    ", sprintf("%3.0f",vmax))
  do while(avg(mslp).gt.10000)
  print("reducing mslp by 10x")
  print("Average mslp was "+avg(mslp))
  mslp = mslp / 10.
  end do
  mslp_s := where(ismissing(mslp),"     ",sprintf("%4.0f",mslp))
  type = "XX"
  
  id_s = sprinti("%02i",toint(storm_id))
  outtext := ""
  delim = ", "
  smin = 0
  do s=smin, subsets-1
    nmi := ind(.not.ismissing(lats(s,:)))
    if(all(ismissing(nmi)))   continue end if
    tech := tcdata@tech(s)
    text := basin+delim\                           ;;;   BASIN        2           WP, EP, CP, IO, SH, AL, LS as defined in basin.dat
     + id_s+delim\                                  ;;;   CY           2          01 to 99, 01 to 49 are real storms, 80 to 89 are test storms, 90 to 99 are INVESTS
     + date+delim\                                  ;;;   YYYYMMDDHH  10          only valid DTGs
     + "03"+delim\                                 ;;;    TECHNUM/MIN  2          00 <= TECHNUM <= 99,  00 <= MIN < 60
     + tech+ delim\                                 ;;;   TECH         4          up to four alphanumeric characters
     + sprintf("%3.0f",time(nmi))+delim\       ;;;    TAU          3          -24 <= TAU <= 240 hours
     + lat_s(s,nmi)+delim\                        ;;;     LatN/S       4          0 <= Lat <= 900  N/S, in tenths of degrees
     + lon_s(s,nmi)+delim\                        ;;;     LonE/W       5          0 <= Lon <=1800  E/W, in tenths of degrees
     + vmax_s(s,nmi)+delim\                      ;;;     VMAX         3          10 <= VMAX <= 250 kt
     + mslp_s(s,nmi)+delim\                      ;;;     MSLP         4           850 to 1050 mb
     + "XX"                                         ;;;     TY           2            as defined in tcdevel.dat
    outtext := array_append_record(outtext, text,0)
  end do
  if(dimsizes(outtext).eq.1 .and. outtext .eq. "") return end if
  outtext@filename = str_lower("A"+basin+sprinti("%02i",toint(storm_id))+year)
  outtext@storm_id = storm_id
return( outtext(1:) ) ;; first line is blank
end

 
undef("adeck_to_data")
function adeck_to_data(adeck, time)
local time_s, adeck_string, dataread, techs, storms, utechs, ustorms, ntech, nstorm, fhrs, nfhrs, afhr, retdata,\
       lats, lons, vmax, mslp, ct, s, t, i
begin

  if(str_get_cols(adeck, -2, -1).eq."gz")
  grep = "zgrep"
  else
  grep = "grep"
  end if
  tech_regex = get_res_value_keep(time, "tech_regex", ".*")
  
  time_s = cd_string(time, "%Y%N%D%H"  )
  tech_regex_str = time_s+str_join(tech_regex, time_s)
  ;print("grep '"+tech_regex_str+"' "+adeck)
  ;print(tech_regex_str)
  adeck_string := systemfunc(grep+" '"+tech_regex_str+"' "+adeck)
  if(all(ismissing(adeck_string)))
       retdata=  new( (/4, 1, 1/), float)
       return retdata
  end if
;  asciiwrite("./test_string", adeck_string)
  ;print(adeck_string)
  csv_string := str_split_csv(adeck_string,",",0)
  dataread := (csv_string)
  delete(adeck_string)
  ;print(dataread)
  ;print(dataread(:,4))
  techs :=  str_strip(dataread(:,4))
  storms := str_strip(dataread(:,1))
  utechs = get_unique_values(techs)
  ustorms = get_unique_values(storms)
  ntech := dimsizes(utechs)
  nstorm := dimsizes(ustorms)
  
  fhrs  = toint(dataread(:,5))
  nfhrs := (max(fhrs)-min(fhrs))/6
  if(nfhrs .lt. 1) return new((/4,1,1/), float) end if
  afhr := fspan(min(fhrs), max(fhrs), nfhrs+1)
  retdata = new( (/4, ntech*nstorm, nfhrs+1/), float)
  retdata!2 = "time"
  retdata&time  = afhr
  retdata!1 = "ens_mem"
  retdata&ens_mem = new( (/ntech*nstorm/), "string")
  retdata@storm_id = new( (/ntech*nstorm/), "string")
  lats = string_latlon_to_numeric(dataread(:,6))
  lons = string_latlon_to_numeric(dataread(:,7)) 
  vmax = tofloat(dataread(:,8))  
  mslp = tofloat(dataread(:,9))  
  
  lats@_FillValue = default_fillvalue(typeof(lats))
  lons@_FillValue = default_fillvalue(typeof(lons))
  lats = where(lats.eq.0 .and. lons.eq.0 , lats@_FillValue, lats)
  lons = where(lats.eq.0 .and. lons.eq.0 , lons@_FillValue, lons)  
  ct=0     
  print("preventing a segfault --- n="+afhr)

  do s=0, nstorm-1
   do t=0, ntech-1
     i := ind(techs.eq.utechs(t) .and. storms.eq.ustorms(s)  .and. fhrs%6 .eq. 0 )
     retdata&ens_mem(ct) = utechs(t)
     retdata@storm_id(ct) = ustorms(s)
;     if( .not. any(afhr .eq. (/fhrs(i)/) ) ) continue end if
     retdata(0,ct, {fhrs(i)}) =(/ lats(i) /)
     retdata(1,ct, {fhrs(i)}) =(/ lons(i) /)
     retdata(2,ct, {fhrs(i)}) =(/ vmax(i) /)
     retdata(3,ct, {fhrs(i)}) =(/ mslp(i) /)
     ct = ct+1
   end do
  end do
retdata@init_time   = time
return retdata
end
   
undef("bdeck_to_data")
function bdeck_to_data(bdeck)
local time_s, adeck_string, dataread, techs, storms, utechs, ustorms, ntech, nstorm, fhrs, nfhrs, afhr, retdata,\
       lats, lons, vmax, mslp, ct, s, t, i
begin

  if(str_get_cols(bdeck, -2, -1).eq."gz")
  cat = "zcat"
  else
  cat = "cat"
  end if
; tech_regex = get_res_value_keep(time, "tech_regex", ".*")
  
;  time_s = cd_string(time, "%Y%N%D%H"  )
;  print(grep+" '"+time_s+tech_regex+"' "+adeck)
  adeck_string := systemfunc(cat+" "+bdeck)
  if(all(ismissing(adeck_string)))
       retdata=  new( (/4, 1, 1/), float)
       return retdata
  end if
  dataread = str_left_strip(str_split_csv(adeck_string,",",0))
  techs :=  dataread(:,4)
  storms := dataread(:,1)
  
  utechs = get_unique_values(techs)
  ustorms = get_unique_values(storms)
  
  ntech := dimsizes(utechs)
  nstorm := dimsizes(ustorms)
  
  times  := ( cd_inv_string(dataread(:,2), "%Y%N%D%H") )
  ntimes = dimsizes(get_unique_values(times))
  retdata = new( (/ntimes, 4, ntech*nstorm/), float)
  retdata!0 = "time"
  retdata&time  = get_unique_values(times)
  retdata&time@units = times@units
  retdata!2 = "ens_mem"
  retdata&ens_mem = new( (/ntech*nstorm/), "string")
  retdata@storm_id = new( (/ntech*nstorm/), "string")
  
  lats = string_latlon_to_numeric(dataread(:,6))
  lons = string_latlon_to_numeric(dataread(:,7)) 
  vmax = tofloat(dataread(:,8))  
  mslp = tofloat(dataread(:,9))  
  
  lats@_FillValue = default_fillvalue(typeof(lats))
  lons@_FillValue = default_fillvalue(typeof(lons))
  lats = where(lats.eq.0 .and. lons.eq.0 , lats@_FillValue, lats)
  lons = where(lats.eq.0 .and. lons.eq.0 , lons@_FillValue, lons)  
  
  ct=0
  do s=0, nstorm-1
   do t=0, ntech-1
     i := ind(dataread(:,4).eq.utechs(t) .and. dataread(:,1).eq.ustorms(s) )
     retdata&ens_mem(ct) = utechs(t)
     retdata@storm_id(ct) = ustorms(s)
     do ii=0, dimsizes(i)-1
     retdata({times(i(ii))},0,ct) =(/ lats(i(ii)) /)
     retdata({times(i(ii))},1,ct) =(/ lons(i(ii)) /)
     retdata({times(i(ii))},2,ct) =(/vmax(i(ii)) /)
     retdata({times(i(ii))},3,ct) =(/ mslp(i(ii)) /)
     end do
     ct = ct+1
   end do
  end do
   
return retdata
end
   



undef("plot_ensemble_track_map") 
function plot_ensemble_track_map(wks, lats, lons, plot_var)
begin
    rlat = (abs(max(lats(:,:)) - min(lats(:,:)))+2)/2.
    rlon = (abs(max(lons(:,:)) - min(lons(:,:)))+2)/2.
    map_dim = min( (/ max( (/rlat, rlon, 12.5/) ), 50/) )
    mlat = max(lats(:,:)) - rlat
    mlon = max(lons(:,:)) - rlon


    ;;; double check some time coord stuff, get fhr 
    curr_time = get_res_value_keep(lons, "init_time", min(lons&time))
    curr_time@units = "hours since 1800-01-01 00:00:00"
	if( min(lons&time) .gt. 1e4)
      fhr = lons&time - curr_time
      adjust_time = curr_time
    else
      fhr = lons&time
      adjust_time = 0
    end if
    minfhr = min(fhr) >0
    minfhr = minfhr + (24 - toint(minfhr)%24)
    maxfhr = max(fhr)
    maxfhr = maxfhr - (toint(maxfhr)%24)

    ;; get analysis hours
    minahr = min(fhr) 
    maxahr =0

    dims = dimsizes(lats)


   cnMin = get_res_value_keep(plot_var, "cnMinLevelValF", min(plot_var) )
   cnMax = get_res_value_keep(plot_var, "cnMaxLevelValF", max(plot_var) )
   cnLevels = get_res_value_keep(plot_var, "cnLevels", fspan(cnMin, cnMax, 21) )

   title = get_res_value_keep(plot_var, "long_name", "")
   units = get_res_value_keep(plot_var, "units", "")
   cnFillPalette = get_res_value_keep(plot_var, "cnFillPalette", "MPL_StepSeq")


    mpres := True
    mpres@mpMaxLatF = mlat + map_dim
    mpres@mpMinLatF = mlat - map_dim;min(lats)-5
    mpres@mpMaxLonF = mlon + map_dim;max(lons)+5
    mpres@mpMinLonF = mlon - map_dim;min(lons)-5
    mpres@mpCenterLonF =mlon
    mpres@tmXBTickSpacingF = 10
    mpres@tmYLTickSpacingF = 10
    mpres@tmXBOn = False
    mpres@tmYLOn = False
    mpres@tmYROn = False
    mpres@tmXTOn = False
    mpres@mpInlandWaterFillColor = "dodgerblue"
    mpres@mpOceanFillColor = "dodgerblue"
    mpres@mpLandFillColor = "darkgoldenrod1"
    mpres@gsnDraw = False
    mpres@gsnFrame = False
    mpres@mpDataSetName          = "Earth..4"
    mpres@mpDataBaseVersion      = "MediumRes"      ; Use GMT coastline
    mpres@mpGeophysicalLineThicknessF =0.5
    mpres@mpNationalLineThicknessF =1.
    mpres@mpOutlineSpecifiers = "North America:states"
    mpres@mpOutlineDrawOrder = "PostDraw"
    mpres@mpOutlineBoundarySets  ="AllBoundaries"
    mpres@mpUSStateLineColor = "darkgoldenrod3" 
    mpres@mpNationalLineColor = "darkgoldenrod4"  
    mpres@mpGeophysicalLineColor = "darkgoldenrod4"  
    mpres@mpOutlineOn = True
    mpres@mpGridAndLimbOn = True
    mpres@mpGridLonSpacingF = 10.
    mpres@mpGridLatSpacingF = 10.
    mpres@mpGridMaskMode = "MaskNotOcean"
    mpres@mpGridLineDashPattern = 12
    mpres@mpGridLineColor = "dodgerblue4"
    mpres@tmXBLabelFontHeightF = 0.008
    mpres@txBackgroundFillColor = (/255.,255.,255.,125./)/255.

    sstres = True
    sstres@cnLinesOn = True
    sstres@gsnDraw = False
    sstres@gsnFrame = False
    sstres@cnLevels := (/24,25,26,27,28,29,30/)
    sstres@cnLevelSelectionMode = "ExplicitLevels"
    sstres@cnLineLabelAngleF = 0
    sstres@cnLineLabelFontHeightF = 0.008
    sstres@cnLineLabelBackgroundColor = -1
    sstres@cnLabelMasking = True
    sstres@cnLineThicknesses := (/0.2, 0.5, 2.5, 1., 1.25, 1.5,1.75/)
    sstres@cnMonoLineThickness = False
    sstres@gsnRightString = ""
    sstres@gsnLeftString  = ""
    sstres@cnInfoLabelOn = False
    sstres@cnLineColor = "dodgerblue4"
    sstres@cnLineLabelFontColor ="dodgerblue4"

    plres = True
    plres@gsLineThicknessF = 2
    plres@gsLineThicknessF := 7.
    plres@gsLineColor = "floralwhite"
    plres@gsLineOpacityF = 0.6
    plres@gsLineThicknessF = 9.

    lb_res := True
    lb_res@lbOrientation = "Horizontal"
    lb_res@lbFormat = "%2.0f"
    lb_res@lbLabelFontHeightF = 0.008
    lb_res@vpWidthF = 0.5
    lb_res@vpHeightF :=lb_res@lbLabelFontHeightF  *2.5
    lb_res@lbLeftMarginF = 0
    lb_res@lbRightMarginF = 0
    lb_res@lbTopMarginF = 0
    lb_res@lbBottomMarginF = 0
    lb_res@lbBoxMajorExtentF =1.
    lb_res@lbBoxMinorExtentF =1.
    lb_res@lbLabelPosition  = "Center"             ; label position
    lb_res@lbLabelJust  = "CenterCenter"             ; label position
    lb_res@lbLabelAlignment = "BoxCenters"
    lb_res@lbTitleOffsetF = 0.
    lb_res@lbTitleString = units
    lb_res@lbTitlePosition = "Right"
    lb_res@lbTitleDirection = "Across"
    lb_res@lbPerimOn = True
    lb_res@lbPerimFill = 0
    lb_res@lbPerimFillColor = "White"
    lb_res@lbTitleFontHeightF = lb_res@lbLabelFontHeightF 

    annres := True
    annres@amOrthogonalPosF=-0.5
    annres@amParallelPosF=-0.5
    annres@amJust = "TopLeft"

    ellres = True
    ellres@gsLineThicknessF = 4.5

    elltxt = True
    elltxt@txFontHeightF = 0.012
    elltxt@txFontColor = "White"

    linetxt = True
    linetxt@txFontHeightF = 0.005
    linetxt@txFontColor = "Black"

    gsMark = True
    gsMark@gsMarkerIndex = 16
    gsMark@gsMarkerSizeF = 0.008

    time_res = True
    time_res@txFontHeightF = lb_res@lbLabelFontHeightF 
    time_res@amJust = "BottomLeft"
    time_res@txPerimOn = True
    time_res@txBackgroundFillColor = "White"
    time_res@tfPolyDrawOrder = "PostDraw"
    
    map := gsn_csm_map(wks, mpres)
    interior_tickmarks(wks, map,mpres)
    do i=0,dims(0)-1
     add_traj_line(wks, map, lats(i,:), lons(i,:), plot_var(i,:), cnLevels, cnFillPalette,plres)
     mxnmi = max(ind(.not.ismissing(lats(i,:))))
     if(.not.ismissing(mxnmi))
       map@$tostring("line_label"+i)$ = gsn_add_text(wks, map,plot_var&tech(i), lons(i,mxnmi), lats(i,mxnmi), linetxt)
     end if
    end do
    map@time_label = gsn_add_text(wks, map, cd_string(curr_time, "%H:%M %D %c %Y"), mpres@mpMinLonF, mpres@mpMinLatF, time_res)
    map@labelbar = apb_create_labelbar(wks, cnLevels(::2), cnFillPalette, lb_res)
;    map@labelbar = apb_create_labelbar(wks, fspan(1,18,18), "MPL_StepSeq", lb_res)
    map@annoid1 = gsn_add_annotation(map,map@labelbar,annres)


    do f=minfhr,maxfhr,24
     elltxt@text = f
     if(num(ismissing(lons(:,{f+adjust_time}))).gt.10) continue end if
     map@$tostring("ell"+f)$ = plot_elipse(wks, map, lons(:,{f+adjust_time}), lats(:,{f+adjust_time}), ellres, elltxt)
    end do

    do f=minahr, maxahr,6
     map@$tostring("anl_"+f)$ = gsn_add_polymarker(wks, map, avg(lons(:,{f+adjust_time})), avg(lats(:,{f+adjust_time})), gsMark)
    end do
    if(isfilepresent(cd_string(curr_time, "~/data/SST/AVHRR/%Y.nc")))
    sstf = addfile(cd_string(curr_time, "~/data/SST/AVHRR/%Y.nc"), "r")
    sst_time = cd_convert(curr_time, sstf->time@units)
    sst_time = floor(sst_time)-2
    sst_time_inds := get1Dindex(sstf->time, (/sst_time, sst_time+1, sst_time+2/) )
    nmi = max(ind(.not.ismissing(sst_time_inds)))
    if(.not.ismissing(nmi))
      sst := short2flt(sstf->sst(sst_time_inds(nmi),:,:))
      wrf_smooth_2d(sst,8)
      sst_map = gsn_csm_contour(wks, sst, sstres)
      overlay(map, sst_map)
    end if
    end if
return map
end


undef("plot_ensemble_intensity_metogram") 
function plot_ensemble_intensity_metogram(wks, plot_var)
begin

plot_var = (/ cubicmsg_n(plot_var,0,0) /)

cnMin = get_res_value_keep(plot_var, "cnMinLevelValF", min(plot_var) )
cnMax = get_res_value_keep(plot_var, "cnMaxLevelValF", max(plot_var) )
cnLevels = get_res_value_keep(plot_var, "cnLevels", fspan(cnMin, cnMax, 21) )
lbLevels = get_res_value_keep(plot_var, "lbLevels", cnLevels )
if(dimsizes(lbLevels) .gt. 10  .and. .not. isatt(plot_var, "lbLevels") )
    lbSpacer = dimsizes(lbLevels)/10
    lbLevels := lbLevels(::lbSpacer)
end if
title = get_res_value_keep(plot_var, "long_name", "")
units = get_res_value_keep(plot_var, "units", "")
cnFillPalette = get_res_value_keep(plot_var, "cnFillPalette", "MPL_StepSeq")

nens =  dimsizes(plot_var(:,0))
techs = plot_var&tech

;print(nens)

xarray = new( dimsizes(plot_var&time)+1, typeof(plot_var&time) )
xarray(:dimsizes(xarray)-2) = plot_var&time
xarray(dimsizes(xarray)-1) = xarray(dimsizes(xarray)-2)+ (xarray(dimsizes(xarray)-2)-xarray(dimsizes(xarray)-3))
;pvs(xarray)
 tm_res := True
 tm_res@txBackgroundFillColor := -1
 tm_res@txPosYF := 0.5
 tm_res@txPerimSpaceF := 0.1
 tm_res@txPerimColor = "white"
 tm_res@txPerimThicknessF = 0
 tm_res@txFontHeightF = 0.013
 res := True
 res@cnFillBackgroundColor = "grey88"
 res@cnMissingValFillColor = "grey88"

 res@gsnRightString = ""
 res@gsnLeftString = ""
 res@gsnCenterString = ""
 res@vpHeightF = 0.3
 res@vpWidthF  = 0.6
 res@gsnDraw = False
 res@gsnFrame = False
 res@cnFillMode = "CellFill"
 res@cnFillOn = True
 res@cnLinesOn = False
 res@cnFillPalette = cnFillPalette
 res@cnLevelSelectionMode = "ExplicitLevels"
 res@cnLevels = cnLevels 
 res@tmYLLabelFontHeightF = (res@vpHeightF/nens)*0.8
 res@tmYLLabels = techs
 res@tmYLValues := ispan(0,nens-1,1)+0.5
 res@tmYLMode = "Explicit"
 xbInds = ind( (mod(plot_var&time,24) .eq. 0 ) .and. (plot_var&time .ne. min(plot_var&time)) .and. (plot_var&time .ne. max(plot_var&time)) )
 if(all(ismissing(xbInds)))
         xbInds := ind( ( mod(plot_var&time,12) .eq. 0 ) .and. (plot_var&time .ne. min(plot_var&time)) .and. (plot_var&time .ne. max(plot_var&time)) )
 end if
 if(all(ismissing(xbInds)))
         xbInds := ind( (mod(plot_var&time,6) .eq. 0 ) .and. (plot_var&time .ne. min(plot_var&time)) .and. (plot_var&time .ne. max(plot_var&time)) )
 end if
 if(all(ismissing(xbInds)))
         xbInds := ind(  (plot_var&time .ne. min(plot_var&time)) .and. (plot_var&time .ne. max(plot_var&time)) )
 end if

 res@tmXBValues := plot_var&time(xbInds)
 res@tmXBLabels := res@tmXBValues+"hr"
 res@tmXBMode = "Explicit"
 res@tmYLMajorLengthF = 0.
 res@sfYArray := ispan(0,nens,1)
 res@sfXArray := xarray
 res@trYMinF = 0
 res@trYMaxF = nens+1.5
 res@trXMinF = max( (/-12,min( (/0, min(plot_var&time) /) ) /) )
 res@trXMaxF = max( (/120, max(plot_var&time)/) )
 res@lbLabelBarOn = False
 res@cnLineLabelsOn = False
 res@cnInfoLabelOn = False
 res@tmXBMajorLengthF = res@vpHeightF
 res@tmXBMajorThicknessF = 0.07
 res@cnConstFLabelString = "No Data Found ~C~ System likely too weak"
 res@cnNoDataLabelString = "No Data Found ~C~ System likely too weak"

 lb_res := True
 lb_res@lbOrientation = "Horizontal"
 lb_res@lbFormat = "%2.0f"
 lb_res@lbLabelFontHeightF = 0.008
; lb_res@lbMinLabelSpacingF =0
 lb_res@vpWidthF = res@vpWidthF/1.8
 lb_res@vpHeightF := tm_res@txFontHeightF*1.2
 lb_res@lbLeftMarginF = 0
 lb_res@lbRightMarginF = 0
 lb_res@lbTopMarginF = 0
 lb_res@lbBottomMarginF = 0
 lb_res@lbBoxMajorExtentF =1.
 lb_res@lbBoxMinorExtentF =1.
 lb_res@lbLabelPosition  = "Center"             ; label position
 lb_res@lbLabelAlignment = "BoxCenters"
 lb_res@lbLabelJust = "CenterCenter"
 lb_res@lbTitleOffsetF = 0.
 lb_res@lbTitleString = units
 lb_res@lbTitlePosition = "Right"
 lb_res@lbTitleDirection = "Across"
 lb_res@lbPerimOn = True
 lb_res@lbPerimFill = 0
 lb_res@lbPerimFillColor = "White"
 lb_res@lbTitleFontHeightF = lb_res@lbLabelFontHeightF *0.8
 lb_res@lbLabelStride = 2
; lb_res@txPerimSpaceF = 0.

 annres := True
 annres@amOrthogonalPosF=-0.5
 annres@amParallelPosF=-0.5
 annres@amJust = "TopLeft"


 plot = gsn_csm_contour(wks, plot_var(:, :), res)
 interior_tickmarks(wks, plot, tm_res)
 interior_subtitles(wks, plot, "", "", title+" | Ens Mem x Time",  tm_res)
 plot@labelbar = apb_create_labelbar(wks, lbLevels, cnFillPalette, lb_res)
 plot@annoid1 = gsn_add_annotation(plot,plot@labelbar,annres)
 return plot
end 



