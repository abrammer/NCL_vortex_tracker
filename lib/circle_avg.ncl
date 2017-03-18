
undef("circle_avg")
function circle_avg(x,r)
local dims, ndims, i
begin
external bx01 "/home/abrammer/Desktop/research/track_cyclones/lib/circle_avg.so"
dims = dimsizes(x)
ndims = dimsizes(dims)

ret = x
ret = default_fillvalue(typeof(x))

if ndims.eq.2 then
  bx01::circle_avg(x, ret, dims(ndims-1), dims(ndims-2),r)
  else if ndims.eq.3 then
    do i=0, dims(0)-1
      bx01::circle_avg(x(i,:,:), ret(i,:,:), dims(ndims-1), dims(ndims-2),toint(r))
    end do
  else
    print("only supports 2 or 3 dimension arrays")
    ret = x
  end if
end if
ret@smoothed = "Smoothed over radius of "+r+" grid cells"
ret@_FillValue = default_fillvalue(typeof(x))
return(ret)
end



undef("circle_avg_m")
function circle_avg_m(x,r)
local dims, ndims, i,r
begin
external bx01 "/home/abrammer/Desktop/research/track_cyclones/lib/circle_avg.so"
dims = dimsizes(x)
ndims = dimsizes(dims)

ret = x
ret = default_fillvalue(typeof(x))
ret!(ndims-2) = "lat"
ret!(ndims-1) = "lon"
if ndims.eq.2 then
  bx01::circle_avg_m(x, ret, ret&lat, dims(ndims-1), dims(ndims-2),tofloat(r))
  else if ndims.eq.3 then
    do i=0, dims(0)-1
      bx01::circle_avg_m(x(i,:,:), ret(i,:,:), ret&lat, dims(ndims-1), dims(ndims-2), tofloat(r) )
    end do
  else
    print("only supports 2 or 3 dimension arrays")
    ret = x
  end if
end if
ret@smoothed = "Smoothed over radius of "+r+" km"
ret@_FillValue = default_fillvalue(typeof(x))
return(ret)
end



undef("circle_avg_m_point")
function circle_avg_m_point(x,r,lat,lon)
local dims, ndims, i,r
begin
external bx01 "/home/abrammer/Desktop/research/track_cyclones/lib/circle_avg.so"
dims = dimsizes(x)
ndims = dimsizes(dims)

lon = where(lon.gt.max(x&lon), lon-360., lon)
lon = where(lon.lt.min(x&lon), lon+360., lon)

if(ndims.eq.2)
 ret = new( dimsizes(lat), typeof(x) )
end if
if(ndims.eq.3)
 ret = new( (/ dims(0), dimsizes(lat)/) , typeof(x) )
end if
if(ndims.gt.3 .or. ndims.lt.2)
 print("Only 2 or 3 Dimension arrays supported")
 return(default_fillvalue(typeof(x) ))
end if

if(isatt(x,"_FillValue"))
 ret@_FillValue = x@_FillValue
else
 x@_FillValue = ret@_FillValue
end if

if(ndims.eq.2)
  nmi = ind(.not.ismissing(lat))
  tmpret = ret(nmi)
  bx01::circle_avg_m_point(x, tmpret, x&lat,x&lon, dims(ndims-1), dims(ndims-2),tofloat(r), lat(nmi), lon(nmi), dimsizes(nmi),ret@_FillValue)
  ret(nmi) = tmpret
  else if ndims.eq.3 then
    do i=0, dims(0)-1
        nmi := ind(.not.ismissing(lat))
      bx01::circle_avg_m_point(x(i,:,:), ret(i,nmi), x&lat,x&lon, dims(ndims-1), dims(ndims-2), tofloat(r), lat(nmi), lon(nmi), dimsizes(nmi),ret@_FillValue)
    end do
  else
    print("only supports 2 or 3 dimension arrays")
  end if
end if
ret@smoothed = "Smoothed over radius of "+r+" km"
ret@_FillValue = default_fillvalue(typeof(x))

ret := rm_single_dims(ret)
return(ret)
end

undef("circle_inout_m_point")
function circle_inout_m_point(x,r,lat,lon)
local dims, ndims, i,r
begin
external bx01 "/home/abrammer/Desktop/research/track_cyclones/lib/circle_avg.so"
dims = dimsizes(x)
ndims = dimsizes(dims)

ret = x
ret = default_fillvalue(typeof(x))
ret!(ndims-2) = "lat"
ret!(ndims-1) = "lon"
ret = 0
if ndims.eq.2 then
  bx01::circle_inout_m_point( ret, ret&lat, ret&lon, dims(ndims-1), dims(ndims-2),tofloat(r),lat, lon)
  else if ndims.eq.3 then
    do i=0, dims(0)-1
      bx01::circle_inout_m_point(ret(i,:,:), ret&lat,ret&lon, dims(ndims-1), dims(ndims-2), tofloat(r),lat, lon )
    end do
  else
    print("only supports 2 or 3 dimension arrays")
    ret = x
  end if
end if
ret@smoothed = "Grid cells within "+r+" km"
ret@_FillValue = default_fillvalue(typeof(x))
return(ret)
end

undef("circle_inout_m_points")
function circle_inout_m_points(x,r,lat,lon)
local dims, ndims, i,r
begin
external bx01 "/home/abrammer/Desktop/research/track_cyclones/lib/circle_avg.so"
dims = dimsizes(x)
ndims = dimsizes(dims)

ret = x
ret = default_fillvalue(typeof(x))
ret!(ndims-2) = "lat"
ret!(ndims-1) = "lon"
nlat = dimsizes(lat)

if(dimsizes(lon).ne.nlat)
  print("Error circle_inout_m_points -- lat and lon must be equal dimensions")
  return(0)
end if

ret = 0
if ndims.eq.2 then
  bx01::circle_inout_m_points( ret, ret&lat, ret&lon, dims(ndims-1), dims(ndims-2),tofloat(r),lat, lon, nlat)
  else if ndims.eq.3 then
    do i=0, dims(0)-1
      bx01::circle_inout_m_points(ret(i,:,:), ret&lat,ret&lon, dims(ndims-1), dims(ndims-2), tofloat(r),lat, lon,nlat)
    end do
  else
    print("only supports 2 or 3 dimension arrays")
    ret = x
  end if
end if
ret@smoothed = "Grid cells within "+r+" km"
ret@_FillValue = default_fillvalue(typeof(x))
return(ret)
end




undef("circle_stddev_m_point")
function circle_stddev_m_point(x,r,lat,lon)
local dims, ndims, i,r, bx01
begin
external bx01 "/home/abrammer/Desktop/research/track_cyclones/lib/circle_avg.so"
dims = dimsizes(x)
ndims = dimsizes(dims)

if(ndims.eq.2)
 ret = new( dimsizes(lat), typeof(x) )
 ave = new( dimsizes(lat), typeof(x) )

end if
if(ndims.eq.3)
 ret = new( (/ dims(0), dimsizes(lat)/) , typeof(x) )
 ave = new( (/ dims(0), dimsizes(lat)/) , typeof(x) )

end if
if(ndims.gt.3 .or. ndims.lt.2)
 print("Only 2 or 3 Dimension arrays supported")
 return(default_fillvalue(typeof(x) ))
end if

if(isatt(x,"_FillValue"))
ret@_FillValue = x@_FillValue
ave@_FillValue = x@_FillValue
else
x@_FillValue = ret@_FillValue
end if

if(ndims.eq.2)
  bx01::circle_avg_m_point(x, ave, x&lat,x&lon, dims(ndims-1), dims(ndims-2),tofloat(r), lat, lon, dimsizes(lat),ret@_FillValue)
  bx01::circle_stddev_m_point(x, ret, ave, x&lat,x&lon, dims(ndims-1), dims(ndims-2),tofloat(r), lat, lon, dimsizes(lat),ret@_FillValue)

  else if ndims.eq.3 then
    do i=0, dims(0)-1
      bx01::circle_avg_m_point(x(i,:,:), ave(i,:), x&lat,x&lon, dims(ndims-1), dims(ndims-2), tofloat(r), lat, lon, dimsizes(lat),ret@_FillValue)
      bx01::circle_stddev_m_point(x(i,:,:), ret(i,:), ave(i,:), x&lat,x&lon, dims(ndims-1), dims(ndims-2),tofloat(r), lat, lon, dimsizes(lat), ret@_FillValue)
    end do
  else
    print("only supports 2 or 3 dimension arrays")
  end if
end if
ret@smoothed = "Smoothed over radius of "+r+" km"
ret@_FillValue = default_fillvalue(typeof(x))

ret := rm_single_dims(ret)
return(ret)
end


undef("circle_gthresh_m_point")
function circle_gthresh_m_point(x,r,thresh,lat,lon)
local dims, ndims, i,r, bx01,thresh
begin
external bx01 "/home/abrammer/Desktop/research/track_cyclones/lib/circle_avg.so"
dims = dimsizes(x)
ndims = dimsizes(dims)

if(ndims.eq.2)
 ret = new( dimsizes(lat), typeof(x) )

end if
if(ndims.eq.3)
 ret = new( (/ dims(0), dimsizes(lat)/) , typeof(x) )
end if
if(ndims.gt.3 .or. ndims.lt.2)
 print("Only 2 or 3 Dimension arrays supported")
 return(default_fillvalue(typeof(x) ))
end if

if(isatt(x,"_FillValue"))
ret@_FillValue = x@_FillValue
else
x@_FillValue = ret@_FillValue
end if


if(ndims.eq.2)
  bx01::circle_gthresh_m_point(x, ret,thresh, x&lat,x&lon, dims(ndims-1), dims(ndims-2),tofloat(r), lat, lon, dimsizes(lat),ret@_FillValue)
  else if ndims.eq.3 then
    do i=0, dims(0)-1
      bx01::circle_gthresh_m_point(x(i,:,:), ret(i,:),thresh, x&lat,x&lon, dims(ndims-1), dims(ndims-2), tofloat(r), lat, lon, dimsizes(lat),ret@_FillValue)
    end do
  else
    print("only supports 2 or 3 dimension arrays")
  end if
end if
ret@smoothed = "Smoothed over radius of "+r+" km"
ret@_FillValue = default_fillvalue(typeof(x))

ret := rm_single_dims(ret)
return(ret)
end



undef("circle_lthresh_m_point")
function circle_lthresh_m_point(x,r,thresh,lat,lon)
local dims, ndims, i,r, bx01,thresh
begin
external bx01 "/home/abrammer/Desktop/research/track_cyclones/lib/circle_avg.so"
dims = dimsizes(x)
ndims = dimsizes(dims)

if(ndims.eq.2)
 ret = new( dimsizes(lat), typeof(x) )

end if
if(ndims.eq.3)
 ret = new( (/ dims(0), dimsizes(lat)/) , typeof(x) )
end if
if(ndims.gt.3 .or. ndims.lt.2)
 print("Only 2 or 3 Dimension arrays supported")
 return(default_fillvalue(typeof(x) ))
end if

if(isatt(x,"_FillValue"))
ret@_FillValue = x@_FillValue
else
x@_FillValue = ret@_FillValue
end if


if(ndims.eq.2)
  bx01::circle_lthresh_m_point(x, ret,thresh, x&lat,x&lon, dims(ndims-1), dims(ndims-2),tofloat(r), lat, lon, dimsizes(lat),ret@_FillValue)
  else if ndims.eq.3 then
    do i=0, dims(0)-1
      bx01::circle_lthresh_m_point(x(i,:,:), ret(i,:),thresh, x&lat,x&lon, dims(ndims-1), dims(ndims-2), tofloat(r), lat, lon, dimsizes(lat),ret@_FillValue)
    end do
  else
    print("only supports 2 or 3 dimension arrays")
  end if
end if
ret@smoothed = "Smoothed over radius of "+r+" km"
ret@_FillValue = default_fillvalue(typeof(x))

ret := rm_single_dims(ret)
return(ret)
end

undef("circle_missing_m_point")
function circle_missing_m_point(x,r,lat,lon)
local dims, ndims, i,r
begin
external bx01 "/home/abrammer/Desktop/research/track_cyclones/lib/circle_avg.so"
dims = dimsizes(x)
ndims = dimsizes(dims)

if(ndims.eq.2)
 ret = new( dimsizes(lat), typeof(x) )
end if
if(ndims.eq.3)
 ret = new( (/ dims(0), dimsizes(lat)/) , typeof(x) )
end if
if(ndims.gt.3 .or. ndims.lt.2)
 print("Only 2 or 3 Dimension arrays supported")
 return(default_fillvalue(typeof(x) ))
end if

if(isatt(x,"_FillValue"))
ret@_FillValue = x@_FillValue
else
x@_FillValue = ret@_FillValue
end if


if(ndims.eq.2)
  bx01::circle_missing_m_point(x, ret, x&lat,x&lon, dims(ndims-1), dims(ndims-2),tofloat(r), lat, lon, dimsizes(lat),ret@_FillValue)
  else if ndims.eq.3 then
    do i=0, dims(0)-1
      bx01::circle_missing_m_point(x(i,:,:), ret(i,:), x&lat,x&lon, dims(ndims-1), dims(ndims-2), tofloat(r), lat, lon, dimsizes(lat),ret@_FillValue)
    end do
  else
    print("only supports 2 or 3 dimension arrays")
  end if
end if
ret@smoothed = "Smoothed over radius of "+r+" km"
ret@_FillValue = default_fillvalue(typeof(x))

ret := rm_single_dims(ret)
return(ret)
end



;;;; change lat degrees to radians
undef("torad")
function torad(x)
local pi, rx
begin
    pi = 3.14159265358979325
    rx = x * pi/180.
    return(rx)
end
;;;; Expects equal distances between grid points.
undef("deg_2_m")
function deg_2_m(lat)
local R, pi, dy,dx,dlat
begin
R = 6378.1
pi = 3.14159265358979325
dy = new( dimsizes(lat), float)
dx = new( dimsizes(lat), float)
dlat = torad(abs( lat(1)-lat(0) ))
dy =  (dlat*R)
dx =  cos( torad(lat) )*(dlat*R)
dy = dy * 1000
dx = dx * 1000  ; convert to metres !!
dx = where(dx.lt.0,1.0,dx)
return( (/dy, dx/) )
end


undef("test_function")
function test_function()
begin
;load "~/NCL/loadFiles.ncl"
wks = gsn_open_wks("x11","Test")
res = True
res@cnFillOn = True
;res@mpProjection = "Satellite"
res@cnLinesOn = False
res@cnFillMode = "RasterFill"
res@mpMinLatF = 0
 res@mpMaxLatF=45
 res@mpMaxLonF = 22.5
 res@mpMinLonF = -22.5

external bx01 "/home/abrammer/Desktop/research/track_cyclones/lib/circle_avg.so"
tf = addfile("/cfsr/data/1999/t.1999.0p5.anl.nc","r")
t = tf->t(0,0,:,:)
t2 = t
bx01::circle_inout_m_point( t, t&lat,t&lon, dimsizes(t&lon), dimsizes(t&lat),500., 15, -15)

print( t({::13},0))
plot= gsn_csm_xy(wks, (/t(:,0),t2(:,0),t3(:,0)/),t&lat, False)



dxdy = deg_2_m(t&lat)
rlat = conform(t,  dxdy(0,:), 0 )
rlon = conform(t,  dxdy(1,:), 0 )
copy_VarCoords(t, rlon)
area = rlat/1000 * rlon/1000
bx01::circle_avg_m(area, t, t&lat, dimsizes(t&lon), dimsizes(t&lat),550)
copy_VarCoords(t,area)
plot= gsn_csm_xy(wks, t(:,0),t&lat, False)

t=0
t({60},{0}) = 100
tout = t
tout_old = t
bx01::circle_avg(t, tout_old, dimsizes(t&lon), dimsizes(t&lat),10)
;tout_old@_FillValue = default_fillvalue(typeof(tout))
dat = systemfunc("date")
bx01::circle_avg_m(t, tout, t&lat, dimsizes(t&lon), dimsizes(t&lat),550)
tout@_FillValue = default_fillvalue(typeof(tout))
tout  =where(tout.gt.1000, tout@_FillValue, tout)
wallClockElapseTime(dat, "" , 0)
replace_ieeenan(tout, tout@_FillValue, 0)
printMinMax(tout,0)
plot_old = gsn_csm_contour_map(wks, t2,res)


;plot = gsn_csm_contour_map(wks, tout_old,res)
plot_old = gsn_csm_contour_map(wks, tout,res)
end



