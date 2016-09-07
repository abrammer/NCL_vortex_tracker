load "~/NCL/loadFiles.ncl"
load "~/NCL/functions/HURDAT_storm_info.ncl"
load "/cfsr/code/get_cfsr_climo.ncl"
load "~/NCL/functions/grid_to_cyl/grid_to_cyl.ncl"
load "~/Desktop/research/noaa_postdoc/hurricane_track_functions.ncl"
load "~/NCL/functions/boxplot.ncl"


begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;   ibTracs                   ;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ncin = addfile("~/data/ibtracs/Allstorms.ibtracs_all.v03r03.nc","r")
lat := short2flt(ncin->source_lat(:,:,9) )
lon := short2flt(ncin->source_lon(:,:,9) )

nmi := ind( dim_num_n(.not.ismissing(lat),1) .ne. 0)
times := ncin->source_time(nmi,0)
year = cd_calendar(times, 0)
storms = nmi( ind(year(:,0).ge.1979))

lat := lat(storms, :190)
lon := lon(storms, :190)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;   erai                      ;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


out_filename = "erai_jtwc.nc"
inf = addfile(out_filename,"r")

sdata := inf->storm_data(:,:190,:,:)
center_locs := sdata(:,:,:,0)
do s=0, dimsizes(sdata&storm)-1
center_locs(s,:,:1) = combine_mand_positions(sdata(s,:,:,:))
end do

slat = sdata(:,:190,0,:)
slon = sdata(:,:190,1,:)
sval = sdata(:,:190,2,:)

dims = dimsizes(center_locs)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


wks = gsn_open_wks("X11","")
tdist :=  new( (/2, dims(0), dims(1)/), float)
do s=0, dimsizes(sdata&storm)-1
tdist(:,s,:) = gc_dist_angle(lat(s,:), lon(s,:), center_locs(s,:,0), center_locs(s,:,1) )
end do

tdist@_FillValue = default_fillvalue(typeof(tdist))
tdist(0,:,:) = where(ismissing(lat), tdist@_FillValue, tdist(0,:,:))
tdist(0,:,:) = where(ismissing(center_locs(:,:,0)), tdist@_FillValue, tdist(0,:,:)
tdist(1,:,:) = where(ismissing(tdist(0,:,:)) , tdist@_FillValue, tdist(1,:,:) )


lat_error = slat - conform(slat, center_locs(:,:,0), (/0,1/) )
lon_error = slon - conform(slon, center_locs(:,:,1), (/0,1/) )

levs = (/ 850, -999.,  850,  700, 700,   925, 850, 700/)
levs@_FillValue = -999
names =(/"vr", "pmsl","g","g","vr","vr", "uv", "uv"/) + levs
names(1) = "pmsl"

res = True
res@gsnDraw = False
res@gsnFrame = False
res@trXMaxF := 3.5
res@trXMinF := -3.5
res@trYMinF = 0
res@trYMaxF = 8
res@tmYLMode = "Explicit"
res@tmYLLabels = names
res@tmYLValues = fspan(0.5, 7.5, dimsizes(names))
boxres = True
boxres@box_width = 0.5

plot = gsn_csm_blank_plot(wks, res)
do v=0, 12
if(.not.all(ismissing(slat(:,:,v))))
add_boxplot_horizontal(wks, plot, ndtooned(lat_error(:,:,v)), v+0.5, boxres)
end if
end do
df(wks, plot)

plot = gsn_csm_blank_plot(wks, res)
boxres = True
boxres@box_width = 0.5
do v=0, 12
if(.not.all(ismissing(slat(:,:,v))))
add_boxplot_horizontal(wks, plot, ndtooned(lon_error(:,:,v)),  v+0.5, boxres)
end if
end do
df(wks, plot)




;;;; spatial bias per variable -- map
grid := new( (/51, 51/), float)
grid!0 = "lat"
grid!1 = "lon"
grid&lat = fspan(-5,5,51)
grid&lon = fspan(-5,5,51)

cres = True
cres@cnFillOn = True
cres@cnLevels = fibonacci(ispan(2, 11,1))
cres@cnLevelSelectionMode = "ExplicitLevels"
grid =0
count := toint(grid)
do v=0,6
grid =0
count := toint(grid)
bin_sum(grid, count, grid&lat, grid&lon, ndtooned(lat_error(:,:,v)), ndtooned(lon_error(:,:,v)),  (ndtooned(lat_error(:,:,v))*0)+1)
plot = gsn_csm_contour(wks, grid, cres)
end do


dims = dimsizes(sdata)
da := new( (/dims(0), dims(1), 2/), float)
loc := new( (/dims(0), dims(1), 2/), float)
wspd := new( (/dims(0), dims(1), 2/), float)

do s=0, dims(0)-1
progress_bar(s, dims(0))

  id = s_ids(s)/10000
  year = s_ids(s) - id*10000

  storm_info :=  HURDAT_storm_info(year,id)
  storm_info := storm_info(ind(toint(storm_info&time)%6.eq.0),:)
  lat := string_latlon_to_numeric(storm_info(:,4))
  lon := string_latlon_to_numeric(storm_info(:,5))
  lat!0 = "time"
  lat&time = storm_info&time

  do t=0,dimsizes(lat)-1
    if(.not.all(ismissing(sdata(s, t, :, :5))))
       cfsr_loc = combine_mand_positions(sdata(s, t, :, :5))
       hur_loc = (/lat(t), lon(t)/)
       wspd(s,t,0) = tofloat( storm_info(t, 6) )
       da(s,t,:) = gc_dist_angle(lat(t), lon(t), cfsr_loc(0), cfsr_loc(1))
       wspd(s,t,1) = get_wind_speed(lat&time(t), cfsr_loc(0), cfsr_loc(1) )
       loc(s,t,:) =  cfsr_loc
    end if
  end do
end do
wspd = where(wspd.lt.0, wspd@_FillValue, wspd)
wspd(:,:,1) = wspd(:,:,1)*1.94384449244





wks = gsn_open_wks("X11","")


opt = True
opt@levs = 1000
grid := grab_cfsr( cd_inv_calendar(1979,01,01,00,00,00,"hours since 1800-01-01",0),"u", opt)
grid := grid(::3,::3)

grid = 0
count := toint(grid)
bin_sum(grid, count,grid&lon, grid&lat,  ndtooned(loc(:,:,1)),ndtooned(loc(:,:,0)),  ndtooned(da(:,:,0)))
count = where(count.le.3, count@_FillValue, count)
grid = grid/count

res = True
res@cnFillOn = True
res@cnLinesOn = False
res@mpMaxLatF = 60
res@mpMinLatF = -10
res@mpMaxLonF = 20
res@mpMinLonF = -110
res@cnLevelSpacingF = 25.
res@cnMinLevelValF = 0.
res@cnMaxLevelValF = 350.
res@cnLevelSelectionMode = "ManualLevels"
res@cnFillMode = "RasterFill"
res@cnRasterSmoothingOn = True
plot = gsn_csm_contour_map(wks, grid, res)

delete(cnLevels)
getvalues plot@contour
    "cnLevels" :cnLevels
end getvalues
hist_res := True
hist_res@gsnHistogramBinIntervals := cnLevels
hist_res@gsnHistogramBarWidthPercent=100
hist_res@tmXBLabelStride = 2
hist_res@tiXAxisString = "Distance [km]"
hist_res@vpHeightF = 0.55
hist_res@vpWidthF = 0.75

hist = gsn_histogram(wks, ndtooned(da(:,:,0)), hist_res)


;annoRes = True
;annoRes@amParallelPosF = 0.5
;annoRes@amOrthogonalPosF = -0.5
;annoRes@amJust = "TopRight"
;dum = gsn_add_annotation(plot, hist, annoRes)
;df(wks, plot)

end
