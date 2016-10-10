undef("df")
procedure df(wks, plot)
begin
draw(plot)
frame(wks)
end 
;check_location(wks, in_time, (/in_lat, in_lon/), sub_scale)

undef("check_location")
procedure check_location(wks, in_time, obs, locs)
local obs, model, time, res,gres, pres, vr_res, ores, mres, ares, map, plot, glot, mslp,vr850
begin


res = True
res@gsnMaximize = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True
res@vpWidthF = 0.6
;res@cnLinesOn = True
pres = True
pres@gsnDraw = False
pres@gsnFrame = False
pres@cnLevelSpacingF = 2.
pres@cnLineLabelsOn = True
pres@cnInfoLabelOn = True
;pres@cnHighLabelsOn  = True
pres@cnLowLabelsOn = True
pres@cnFillOn = False
pres@cnLabelMasking = True
pres@cnLineLabelBackgroundColor = -1
pres@gsnAddCyclic = True
pres@gsnRightString = ""
pres@gsnLeftString = ""
pres@cnInfoLabelOn= False

;
gres = True
gres@gsnDraw = False
gres@gsnFrame = False
;gres@cnLevelSpacingF = 10.
gres@cnLineColor = "Brown"
gres@cnLineThicknessF =3.
gres@cnLineLabelsOn = False
gres@cnInfoLabelOn = False
;gres@cnHighLabelsOn  = True
;gres@cnLowLabelsOn = True
gres@cnLabelMasking = True
gres@cnLineLabelBackgroundColor = -1
gres@gsnAddCyclic = False
gres@gsnRightString = ""
gres@gsnLeftString = ""
gres@cnMinLevelValF = 0
gres@cnMaxLevelValF = 30
gres@cnLevelSpacingF = 5
gres@cnLevelSelectionMode = "ManualLevels"

vr_res = True

vr_res@gsnAddCyclic = True
vr_res@gsnDraw = False
vr_res@gsnFrame = False
vr_res@cnFillOn = True
vr_res@cnFillMode = "RasterFill"
vr_res@cnFillPalette = "sunshine_diff_20lev"
vr_res@cnLevelSelectionMode = "ExplicitLevels"
vr_res@cnLevels = (/ -75,-55,-35,-15,0,5,7.5,10,12.5,15,20,25,30,35,40,45,50,75,100,125,150/)*(10^-5)
vr_res@cnLinesOn = False
vr_res@cnRasterSmoothingOn = True
vr_res@lbOrientation = "vertical"
vr_res@pmLabelBarOrthogonalPosF = -0.15
vr_res@lbLabelBarOn = False
vr_res@cnLineLabelsOn = False
vr_res@cnInfoLabelOn = False
vr_res@gsnRightString = ""
vr_res@gsnLeftString = ""



ores = True
ores@gsMarkerIndex = 16
ores@gsMarkerSizeF = 0.008

mres = True
mres@gsMarkerIndex = 10
ores@gsMarkerSizeF = 0.015

ares = True
ares@gsMarkerIndex = 0
ares@txFontHeightF = 0.012

mslp  =grab_data(in_time, "pmsl", False)
mslp = mslp/100.
opt = True
opt@levs = 700
;g850  =grab_data(in_time, "cv", opt)
opt@levs = 850
vr850  =grab_data(in_time, "vr", opt)
;vr850 = circle_avg_m(vr850, 150.)

dx = abs(obs(1) - locs(1)) > 5
dy = abs(obs(0) - locs(0)) > 5
dxy = max( (/dx, dy/) )

res@mpMinLatF = avg((/obs(0),locs(0)/))-dxy
res@mpMaxLatF = avg((/obs(0),locs(0)/))+dxy

res@mpMinLonF = avg((/obs(1),locs(1)/))-dxy
res@mpMaxLonF = avg((/obs(1),locs(1)/))+dxy

map = gsn_csm_map(wks,res)
vlot = gsn_csm_contour(wks, vr850, vr_res)
plot = gsn_csm_contour(wks, mslp, pres)
;glot = gsn_csm_contour(wks, g850({res@mpMinLatF:res@mpMaxLatF},{res@mpMinLonF:res@mpMaxLonF}), gres)

overlay(map, vlot)

overlay(map, plot)
;overlay(map, glot)
;


map@obs = gsn_add_polymarker(wks, map, obs(1), obs(0), ores)
map@mod1 = gsn_add_polymarker(wks, map, locs(1), locs(0), mres)
if(isatt(locs, "lat"))
map@mod2 = gsn_add_text(wks, map,tostring(ispan(1,dimsizes(locs@lat),1)), locs@lon, locs@lat, ares)
end if

df(wks, map)

end


undef("check_centerfix")
procedure check_centerfix(wks, data, locs)
local obs, model, time, res,gres, pres, vr_res, ores, mres, ares, map, plot, glot, mslp,vr850
begin


res = True
res@gsnMaximize = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True
res@vpWidthF = 0.6
;res@cnLinesOn = True
pres = True
pres@gsnDraw = False
pres@gsnFrame = False
pres@cnLevelSpacingF = 2.
pres@cnLineLabelsOn = True
pres@cnInfoLabelOn = True
;pres@cnHighLabelsOn  = True
pres@cnLowLabelsOn = True
pres@cnFillOn = False
pres@cnLabelMasking = True
pres@cnLineLabelBackgroundColor = -1
pres@gsnAddCyclic = True
pres@gsnRightString = ""
pres@gsnLeftString = ""
pres@cnInfoLabelOn= False

;
gres = True
gres@gsnDraw = False
gres@gsnFrame = False
;gres@cnLevelSpacingF = 10.
gres@cnLineColor = "Brown"
gres@cnLineThicknessF =3.
gres@cnLineLabelsOn = False
gres@cnInfoLabelOn = False
;gres@cnHighLabelsOn  = True
;gres@cnLowLabelsOn = True
gres@cnLabelMasking = True
gres@cnLineLabelBackgroundColor = -1
gres@gsnAddCyclic = False
gres@gsnRightString = ""
gres@gsnLeftString = ""
gres@cnMinLevelValF = 0
gres@cnMaxLevelValF = 30
gres@cnLevelSpacingF = 5
gres@cnLevelSelectionMode = "ManualLevels"

vr_res = True

vr_res@gsnAddCyclic = False
vr_res@gsnDraw = False
vr_res@gsnFrame = False
vr_res@cnFillOn = True
vr_res@cnFillMode = "RasterFill"
;vr_res@cnFillPalette = "sunshine_diff_20lev"
;vr_res@cnLevelSelectionMode = "ExplicitLevels"
;vr_res@cnLevels = (/ -75,-55,-35,-15,0,5,7.5,10,12.5,15,20,25,30,35,40,45,50,75,100,125,150/)*(10^-5)
vr_res@cnLinesOn = False
vr_res@cnRasterSmoothingOn = True
vr_res@lbOrientation = "vertical"
vr_res@pmLabelBarOrthogonalPosF = -0.15
vr_res@lbLabelBarOn = False
vr_res@cnLineLabelsOn = False
vr_res@cnInfoLabelOn = False
vr_res@gsnRightString = ""
vr_res@gsnLeftString = ""



ores = True
ores@gsMarkerIndex = 16
ores@gsMarkerSizeF = 0.008

mres = True
mres@gsMarkerIndex = 10
ores@gsMarkerSizeF = 0.015

ares = True
ares@gsMarkerIndex = 0
ares@txFontHeightF = 0.012

;mslp  =grab_data(in_time, "pmsl", False)
;mslp = mslp/100.
;opt = True
;opt@levs = 700
;;g850  =grab_data(in_time, "cv", opt)
;opt@levs = 850
;vr850  =grab_data(in_time, "vr", opt)
;vr850 = circle_avg_m(vr850, 150.)

;dx = abs(obs(1) - locs(1)) > 5
;dy = abs(obs(0) - locs(0)) > 5
;dxy = max( (/dx, dy/) )
dxy = 5
res@mpMinLatF = locs(0) -dxy
res@mpMaxLatF = locs(0) +dxy
res@mpMinLonF = locs(1)-dxy
res@mpMaxLonF = locs(1)+dxy

map = gsn_csm_map(wks,res)
;vlot = gsn_csm_contour(wks, vr850, vr_res)
plot = gsn_csm_contour(wks, data({res@mpMinLatF-1:res@mpMaxLatF+1},{res@mpMinLonF-1:res@mpMaxLonF+1}), vr_res)
;glot = gsn_csm_contour(wks, g850({res@mpMinLatF:res@mpMaxLatF},{res@mpMinLonF:res@mpMaxLonF}), gres)

;overlay(map, vlot)

overlay(map, plot)
;overlay(map, glot)
;


;map@obs = gsn_add_polymarker(wks, map, obs(1), obs(0), ores)
map@mod1 = gsn_add_polymarker(wks, map, locs(1), locs(0), mres)
;if(isatt(locs, "lat"))
;map@mod2 = gsn_add_text(wks, map,tostring(ispan(1,dimsizes(locs@lat),1)), locs@lon, locs@lat, ares)
;end if

df(wks, map)

end


;print(outdata(t,:,2))
;newdata = (get_mass_center(g850, obs))
;
;ores = True
;ores@gsMarkerIndex = 16
;mres = True
;mres@gsMarkerIndex = 10
;ares = True
;ares@gsMarkerIndex = 0
;ares@txFontHeightF = 0.012
;
;map@obs = gsn_add_polymarker(wks, map, obs(1), obs(0), ores)
;map@model = gsn_add_polymarker(wks, map, model(1), model(0), mres)
;map@vort := gsn_add_polymarker(wks, map, outdata(t,1,0), outdata(t,0,0), ores)
;
;map@allmodel := gsn_add_text(wks, map, ispan(1,8,1), outdata(t,1,:7), outdata(t,0,:7), ares)
;map@newmodel = gsn_add_polymarker(wks, map, newdata(0,1), newdata(0,0), mres)
;
;hdt_line = gsn_add_polyline(wks, plot, lon, lat, False)
;cfs_line = gsn_add_polyline(wks,map, dim_avg_n(outdata(:,1,:), 1), dim_avg_n(outdata(:,0,:), 1), False)
;
;df(wks, map)
;
;
;
