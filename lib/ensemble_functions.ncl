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

if(rankz.eq.1 .and. dimz(0) .eq. dimsizes(fcst_lat(0,:)))
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

alnormx =  dbx / vlen
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
nens = get_res_value_keep(tcdata, "nens", 21)

; nmi := ind(dim_num_n(.not.ismissing(tcdata(0,:,:)),0) .gt. nens/10. ) 
 mean_locs:=dim_avg_n_Wrap(tcdata(:,:,:),1)
 wgt := gc_latlon(mean_locs(0,{72}),mean_locs(1,{72}), tcdata(0,:,{72}),tcdata(1,:,{72}),0,4) > 25.
 wgt@_FillValue = default_fillvalue(typeof(wgt))
 wgt = where(ismissing(tcdata(1,:,{72})), wgt@_FillValue, wgt)
 ti = ind_nearest_coord(72, mean_locs&time, 0)
 do t=ti+1, dimsizes(mean_locs&time)-1
   mean_locs(:,t) = dim_avg_wgt_n(tcdata(:,:,t), 1/wgt, nens/6, 1)
 end do
; mean_locs_hr = linint1_n_Wrap(mean_locs&time, mean_locs, False, fspan(min(mean_locs&time), max(mean_locs&time), 100), 0, 1)
 mean_locs = runave_n(mean_locs, 4, 1, 1)
; mean_locs = (/ mean_locs_hr(:,{mean_locs&time}) /)
 
 nmi = max(ind(.not.ismissing(mean_locs(0,:))))
 maxhr = totype(min( (/rmaxhr, mean_locs&time(nmi)/) ), typeof(rmaxhr) )
 mean_lat := mean_locs(0,{:maxhr})
 mean_lon := mean_locs(1,{:maxhr})

 mean_lat =  cubicmsg(mean_lat, False)
 mean_lon =  cubicmsg(mean_lon, False)

 mean_track := track_dir_dist(mean_lat({:maxhr}), mean_lon({:maxhr}))
 errors := along_across_error(mean_lat, mean_lon, tcdata(0,:,{mean_lat&time}), tcdata(1,:,{mean_lat&time}))
 ens_spread := dim_stddev_n(errors(1,:,:),0)
 do t=0, dimsizes(ens_spread)-1
  ens_spread(t) = dim_max(ens_spread(:t))
 end do
 ens_spread = tofloat( ens_spread  + (mean_lon&time) ) > 50.

 lcone := new( (/2,dimsizes(ens_spread)/), float)
 rcone := new((/2,dimsizes(ens_spread)/), float)
 do t=0, dimsizes(ens_spread)-1
  if(.not. ismissing(mean_lat(t)) .and. .not. ismissing(mean_track(0,t)) )
     lcone(:,t) = advect_by_km_radians(mean_lat(t), mean_lon(t), ens_spread(t), mean_track(0,t)-1.57)
     rcone(:,t) = advect_by_km_radians(mean_lat(t), mean_lon(t), ens_spread(t), mean_track(0,t)+1.57)
   end if
 end do
;print(mean_track(0,:))
;mean_track(0,:) = where(mean_track(0,:).gt. pi/2., mean_track(0,:)-pi, mean_track(0,:))
;print(mean_track(0,:))

 end_cone:= new( (/2, 35/), float)
 start_cone:= new( (/2, 35/), float)
 t= dimsizes(ens_spread)-1
 do r=0,34
     rad = -1.55 + (pi/35.)*r
    end_cone(:,r) = advect_by_km_radians(mean_lat(t), mean_lon(t), ens_spread(t), mean_track(0,t)-rad)
    start_cone(:,r) = advect_by_km_radians(mean_lat(0), mean_lon(0), ens_spread(0), mean_track(0,0)+(pi-rad))
end do

cone_lats := array_append_record( rcone(0,:), end_cone(0,:),0)
cone_lons := array_append_record( rcone(1,:), end_cone(1,:),0)
cone_lats := array_append_record( cone_lats, lcone(0,::-1),0)
cone_lons := array_append_record(cone_lons, lcone(1,::-1),0)
cone_lats := array_append_record( cone_lats, start_cone(0,:),0)
cone_lons := array_append_record(cone_lons, start_cone(1,:),0)
cone_lats := array_append_record( cone_lats, cone_lats(0),0)
cone_lons := array_append_record(cone_lons, cone_lons(0),0)

cone = new( (/2, 1001/), float)
cone(0,:) = runave( linint1(fspan(0,100, dimsizes(cone_lats)), cone_lats, True, fspan(0,100,1001), 0), 40, -1)
cone(1,:) = runave( linint1(fspan(0,100, dimsizes(cone_lats)), cone_lons, True, fspan(0,100,1001), 0), 40, -1)

coned :=    remove_loops(cone) 
npts = dimsizes(coned(0,:))
cone(0,:) = runave( linint1(fspan(0,100, npts), coned(0,:), True, fspan(0,100,1001), 0), 40, -1)
cone(1,:) = runave( linint1(fspan(0,100, npts), coned(1,:), True, fspan(0,100,1001), 0), 40, -1)


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
print(any(ismissing(cone_lons)))
    smth = toint(dimsizes(cone_lons)*0.15)
    cone_lons := runave(cone_lons, smth, -1)
    cone_lats := runave(cone_lats, smth, -1)
    cone_lons(nlons-1) = cone_lons(0)
    cone_lats(nlons-1) = cone_lats(0)
    print(cone_lons(nlons-1)+"   "+cone_lons(0))
    map@$tostring("cone"+ct)$ = gsn_add_polyline(wks,map,  cone_lons, cone_lats, restemp)
;end do
end

