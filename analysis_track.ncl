load "~/Desktop/research/track_cyclones/lib/circle_avg.ncl"
load "~/Desktop/research/track_cyclones/util/progress_bar.ncl"

load "~/Desktop/research/track_cyclones/lib/get_cfsr_climo.ncl"  ;; this script needs to be changed
load "~/Desktop/research/track_cyclones/lib/tracking_functions.ncl"




year=2016

begin

stime = cd_inv_calendar(year, 07,0,00,00,00,"hours since 1800-01-01 00:00:00",0)
etime = cd_inv_calendar(year, 11,00,00,00,00,"hours since 1800-01-01 00:00:00",0) ; trick to get last day of prev month
times := fspan(stime, etime, toint((etime-stime)/6)+1)
times@units = stime@units

track_locs := new( (/1000,  8, dimsizes(times)/) , float)
track_locs@ct = 0
track_locs!2 = "time"
track_locs&time = times


plot = False
if(plot)
wkst = "x11"
wkst@wkPause = False
wks = gsn_open_wks(wkst,"tracking_images")
end if

greys = span_color_indexes("gsltod", 7)
mpres = True
mpres@mpMinLatF = -10
mpres@mpMaxLatF = 50
mpres@mpMinLonF = -180
mpres@mpMaxLonF = 180
mpres@gsnDraw = False
mpres@gsnFrame =  False
shade_res = True
shade_res@gsnDraw = False
shade_res@gsnFrame = False
shade_res@tmXTOn = False
shade_res@tmYROn = False
shade_res@gsnRightString =""; "Valid: "+cd_calendar(validate,-3)
shade_res@gsnLeftString = "";""+date+" Forecast: "+(t*6)+"hr"
shade_res@gsnAddCyclic = True
shade_res@cnFillPalette = "sunshine_diff_20lev"
shade_res@cnLevelSelectionMode = "ExplicitLevels"
shade_res@cnLevels = (/ -75,-55,-35,-15,0,5,7.5,10,12.5,15,20,25,30,35,40,45,50,75,100,125,150/)
shade_res@cnInfoLabelOn = False
shade_res@cnLineLabelsOn = False
shade_res@cnFillOn = True
shade_res@cnLinesOn = False
shade_res@cnFillMode = "RasterFill"
shade_res@cnRasterSmoothingOn = True
shade_res@lbLabelBarOn = True
shade_res@lbLabelAngleF = 45
shade_res@lbOrientation = "horizontal"
shade_res@lbLabelFontHeightF = 0.006
shade_res@lbLabelAutoStride = False
copy_VarAtts(mpres, shade_res)

opt := True
opt@levs = (/700/)
opt@min_threshold = 0.5e-5
opt@rkm := 400.

ntimes = dimsizes(times)-1
do t=0, ntimes
  progress_bar(t, ntimes)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;  Find new starting points
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       ;; lat bounds and lon bounds for starting targets
  maxes_mass := get_starting_targets(times(t), (/5,35/), (/-150,40/), opt)
  track_locs(track_locs@ct:track_locs@ct+dimsizes(maxes_mass(:,0))-1,:2,t) = (/ maxes_mass /)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;  Combine potential locations into tracked locations
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  track_locs(:,:2,t) = get_multi_positions(times(t), track_locs(:,0,t), track_locs(:,1,t), opt)
  track_locs(:,:,t) = unique_locations(track_locs(:,:,t), opt@rkm)
  filter_tracks(track_locs, t, opt)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;  Propagate/ Extrapolate systems
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  if( all(ismissing(track_locs(:,0,t))) .or. t.eq. dimsizes(times)-1)    continue   end if
  track_locs(:,:1,t+1) = (/ advect_tracks(times(t), track_locs(:,:,t))  /); 
  
  
  if(plot)
     vr = grab_data(times(t), "vr", opt)
     vr = circle_avg_m(vr, 200.)*10^6
     map = gsn_csm_contour_map(wks, vr, shade_res)
     do tt=t-6,t
       if(tt.lt.0)continue end if
       mres = True
       mres@gsMarkerIndex = 16
       mres@gsMarkerColor = greys(t-tt)
       map@$tostring("points"+tt)$ = gsn_add_polymarker(wks, map, track_locs(:,1,tt),  track_locs(:,0,tt),mres)
     end do
     draw(map)
     frame(wks)
   end if
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end do


track_locs(:,1,:) = where(track_locs(:,1,:).gt.180., track_locs(:,1,:)-360., track_locs(:,1,:))
track_locs(:,1,:) = where(track_locs(:,1,:).lt.-180., track_locs(:,1,:)+360., track_locs(:,1,:))

nmi := ind(  (  (dim_num_n(.not.ismissing(track_locs(:,0,:)),1) .ge. 8 ) .and.  (dim_sum_n(  sqrt(  (track_locs(:,1,1:) - track_locs(:,1,:dimsizes(times)-2))^2 + ( track_locs(:,1,1:) - track_locs(:,1,:dimsizes(times)-2))^2 ),1 ).gt. 10 ) ))
ct := dimsizes(nmi)
track_locs(:ct-1,:,:) = track_locs(nmi,:,:)
track_locs(ct:,:,:) = track_locs@_FillValue
track_locs@ct = ct

nmi := ind(  dim_max_n(track_locs(:,0,:),1) .gt. 5 )
ct := dimsizes(nmi)
track_locs(:ct-1,:,:) = track_locs(nmi,:,:)
track_locs(ct:,:,:) = track_locs@_FillValue
track_locs@ct = ct

nmi := ind(  dim_min_n(track_locs(:,1,:),1) .lt. -10 )
ct := dimsizes(nmi)
track_locs(:ct-1,:,:) = track_locs(nmi,:,:)
track_locs(ct:,:,:) = track_locs@_FillValue
track_locs@ct = ct


track_locs := track_locs(:ct-1,:,:)

trash = get_res_value( track_locs, "vals", 0)
trash := get_res_value( track_locs, "ct", 0)

track_locs!0 = "id"
track_locs&id = ispan(0, dimsizes(track_locs(:,0,0))-1, 1)

track_locs!2 = "time"
track_locs&time = times

ncout = addfile(""+year+"_cfsr_pregen_tracks.nc","c")
ncout->track_lat = track_locs(:,0,:)
ncout->track_lon = track_locs(:,1,:)
ncout->track_val = track_locs(:,2,:)
ncout->track_val@units = "10e-5"
ncout->track_val@longname = "200km 850 hPa Rel. Vort."
ncout->all_track_data = track_locs


status_exit(0)
end
status_exit(1)
;
;opt@levs = 850
;cv = grab_data(times(t), "vr", opt)
;cv = circle_avg_m(cv, 500)
;res@gsnCenterString = "fhr:"+t
;plot = gsn_csm_contour_map(wks, cv({13:23},{-39:-29}), res)
;;dum := gsn_add_text(wks, plot, sprinti("%3i",ispan(1,maxes,1)), maxes@lon, maxes@lat, txres)
;;dum := gsn_add_text(wks, plot, sprinti("%3i",ispan(1,dimsizes(maxes_mass(:,0)),1)), maxes_mass(:,1), maxes_mass(:,0), txres)
;;dum1 := gsn_add_text(wks, plot, sprinti("%3i",ispan(1,dimsizes(track_locs(:,1,t)),1)), track_locs(:,1,t), track_locs(:,0,t), txres)
;dum :=new(1000, graphic)
;do i=0,0;999
;if(.not.all(ismissing(track_locs(i,1,:))))
;dum(i) = gsn_add_polyline(wks, plot, track_locs(i,1,:t), track_locs(i,0,:t), False)
;end if
;end do
;;dum := gsn_add_text(wks, plot, sprinti("%3i",ispan(1,dimsizes(tmp_wave_locs(0,0,:)),1)), tmp_wave_locs(0,1,:), tmp_wave_locs(0,0,:), txres)
;
;;dum2 := gsn_add_text(wks, plot, "-", new_locs(:,1), new_locs(:,0), txres)
;
;df(wks, plot)
;;end if
;
;end do
;
;
;dist = new(200, float)
;do i=0,999
;if(.not.all(ismissing(track_locs(i,1,:))))
;dist(i) = dim_sum_n(  sqrt(  (track_locs(:,1,1:) - track_locs(:,1,:dimsizes(times)-2))^2 + ( track_locs(:,1,1:) - track_locs(:,1,:dimsizes(times)-2))^2 ),1 ) 
;end if
;end do
;
;
;
;tmp_wave_locs := get_wave_positions(cv&time(1), new_locs(:,0), new_locs(:,1), rkm)
;new_locs_new := combine_mand_positions( tmp_wave_locs )
;
;
;plot = gsn_csm_contour_map(wks, cv(t,:,:), res)
;;dum := gsn_add_text(wks, plot, sprinti("%3i",ispan(1,maxes,1)), maxes@lon, maxes@lat, txres)
;;dum := gsn_add_text(wks, plot, sprinti("%3i",ispan(1,dimsizes(maxes_mass(:,0)),1)), maxes_mass(:,1), maxes_mass(:,0), txres)
;dum := gsn_add_text(wks, plot, sprinti("%3i",ispan(1,dimsizes(wave_locs(:,0)),1)), wave_locs(:,1), wave_locs(:,0), txres)
;dum := gsn_add_text(wks, plot, sprinti("%3i",ispan(1,dimsizes(wave_locs(:,0)),1)), new_locs(:,1), new_locs(:,0), txres)
;dum1 := gsn_add_text(wks, plot, "_"+sprinti("%3i",ispan(1,dimsizes(wave_locs(:,0)),1)), new_locs_new(:,1), new_locs_new(:,0), txres)
;
;;dum1 := gsn_add_text(wks, plot, sprintf("%2.0f",maxes_vwnd), maxes_mass(:,1), maxes_mass(:,0), txres)
;;dum := gsn_add_text(wks, plot, sprintf("%2.0f",wave_locs(:,0)), wave_locs(:,1), wave_locs(:,0), txres)
;
;
;df(wks, plot)
;
;
;
;
;settings = True
;settings@lon_bounds = (/-35,5/)
;settings@min_time = 8
;settings@min_lon = 10
;tracks1d :=  track_hov(cva(:,{-60:20}), settings) 
;
;tracks2d := make_2d_track(cv(:,{70000},:,:), track_hov(cva(:,:), settings) )
;
;
;
;wkstype = "X11"
;wkstype@wkPause = False
;wks =gsn_open_wks(wkstype,"")
;res = True
;res@mpMaxLatF = 40
;res@mpMinLatF = 0
;res@mpMaxLonF = 20
;res@mpMinLonF = -120
;res@cnFillOn = True
;res@cnFillPalette = "wamwg256"
;res@cnFillMode = "RasterFill"
;res@cnRasterSmoothingOn = True
;res@cnMinLevelValF = 0.25
;res@cnMaxLevelValF = 6.
;res@cnLevelSpacingF  := 0.5
;res@cnLevelSelectionMode = "ManualLevels"
;res@cnLinesOn = False
;res@gsnDraw = False
;res@gsnFrame = False
;res@trYReverse = False
;mres = True
;mres@gsMarkerIndex = 16
;txres = True
;txres@txFontHeightF = 0.008
;
;
;plot = gsn_csm_contour(wks, cva, res)
;lres = True
;lres@gsLineColor = "White"
;lres@gsLineThicknessF = 3.
;do t=0, dimsizes(tracks1d(:,0))-1
;  plot@$unique_string("s")$ = gsn_add_polyline(wks, plot, tracks1d(t,:), tracks1d&time, lres)
;end do
;df(wks, plot)
;
;
;
