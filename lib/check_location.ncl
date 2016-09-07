
;check_location(wks, in_time, (/in_lat, in_lon/), sub_scale)

undef("check_location")
procedure check_location(wks, in_time, obs, locs)
local obs, model, time, res,gres, pres, vr_res, ores, mres, ares, map, plot, glot
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
pres@cnLevelSpacingF = 50.
pres@cnLineLabelsOn = True
pres@cnInfoLabelOn = True
;pres@cnHighLabelsOn  = True
;pres@cnLowLabelsOn = True
pres@cnFillOn = False
pres@cnLabelMasking = True
pres@cnLineLabelBackgroundColor = -1
pres@gsnAddCyclic = True

;
gres = True
gres@gsnDraw = False
gres@gsnFrame = False
;gres@cnLevelSpacingF = 10.
gres@cnLineColor = "Brown"
gres@cnLineThicknessF =3.
gres@cnLineLabelsOn = True
gres@cnInfoLabelOn = False
;gres@cnHighLabelsOn  = True
;gres@cnLowLabelsOn = True
gres@cnLabelMasking = True
gres@cnLineLabelBackgroundColor = -1
gres@gsnAddCyclic = False
vr_res = True

vr_res@gsnAddCyclic = True
vr_res@gsnDraw = False
vr_res@gsnFrame = False
vr_res@cnFillOn = True
vr_res@cnFillMode = "RasterFill"
vr_res@cnLevelSelectionMode = "ManualLevels"
vr_res@cnMinLevelValF = 0.5*10^-5
vr_res@cnMaxLevelValF = 15*10^-1
vr_res@cnLevelSpacingF = 15*10^-3
vr_res@cnLinesOn = False
vr_res@cnRasterSmoothingOn = True
vr_res@lbOrientation = "vertical"
vr_res@pmLabelBarOrthogonalPosF = -0.15
ores = True
ores@gsMarkerIndex = 16

mres = True
mres@gsMarkerIndex = 10

ares = True
ares@gsMarkerIndex = 0
ares@txFontHeightF = 0.012

mslp  =grab_data(in_time, "pmsl", False)
opt = True
opt@levs = 700
g850  =grab_data(in_time, "cv", opt)
opt@levs = 850
vr850  =grab_data(in_time, "cv", opt)
vr850 = circle_avg_m(vr850, 300.)

res@mpMinLatF = obs(0)-10
res@mpMaxLatF = obs(0)+10
res@mpMinLonF = obs(1)-10
res@mpMaxLonF = obs(1)+10

map = gsn_csm_map(wks,res)
vlot = gsn_csm_contour(wks, vr850, vr_res)
plot = gsn_csm_contour(wks, mslp, pres)
glot = gsn_csm_contour(wks, g850({res@mpMinLatF:res@mpMaxLatF},{res@mpMinLonF:res@mpMaxLonF}), gres)

overlay(map, vlot)

overlay(map, plot)
overlay(map, glot)
;


map@obs = gsn_add_polymarker(wks, map, obs(1), obs(0), ores)
map@mod1 = gsn_add_polymarker(wks, map, locs(1), locs(0), mres)
;map@allmodel := gsn_add_text(wks, map, ispan(1,8,1), outdata(t,1,:7), outdata(t,0,:7), ares)

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
