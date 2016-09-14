load "~/NCL/loadFiles.ncl"
load "~/NCL/functions/print_r.ncl"
load "~/NCL/functions/percentile.ncl"
load "lib/get_cfsr_climo.ncl"
load "lib/make_lon_monotonic.ncl"
load "lib/hurricane_track_functions.ncl"
load "lib/check_location.ncl"



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





begin

c_inds = (/0,3, 7,9,10,13,15,16,18/)
centers = (/  "hurdat_atl","td9636","reunion","atcf","mlc_natl","ds824_sh","ds824_ni","bom","ds824_au","jtwc_sh","jtwc_wp","td9635","ds824_wp","jtwc_io","cma","hurdat_epa","jtwc_ep","ds824_ep","jtwc_cp","tokyo","neumann","hko","cphc","wellington","newdelhi","nadi"/)


out_filename = "data/cfsr_ibtracs_v03r08.nc"

if(.not.fileexists(out_filename))
     setfileoption("nc","Format", "NetCDF4") 
     ncdf = addfile(out_filename ,"c")  ; open output netCDF file

     ;===================================================================
     ; create global attributes of the file (optional)
     ;===================================================================
     fAtt               = True            ; assign file attributes
     fAtt@title         = "ERAI Storm Data recentered from HURDAT"
     fAtt@Conventions   = "None"
     fAtt@creation_date = systemfunc ("date")
     fileattdef( ncdf, fAtt )            ; copy file attributes

     ;===================================================================
     ; make time an UNLIMITED dimension; recommended  for most applications
     ;===================================================================
     filedimdef(ncdf,"storm",-1,True) 

     outdata := new( (/1,191, 3, 13/), float)
     outdata!0 = "storm"
     outdata!1 = "time"
     outdata!2 = "data"
     outdata!3 = "variables"

     time = ispan(0,190*6, 6)
     time!0 = "time"
     time&time = time

     storm_time = new( (/1, 191/), double)
     storm_time!0 = "storm"
     storm_time&storm = 0
     storm_time!1 = "time"
     storm_time&time = time
     storm_time@units = "hours since 1800-01-01 00:00:00"

     storm = 0
     storm!0 = "storm"

     storm_id = new(1, integer)
     storm_id!0 = "storm"
     storm_id&storm = storm


     variables = ispan(0,12,1)
     variables!0= "variables"
     variables&variables = variables

     data      = ispan(0,2,1)
     data!0 = "data"
     data&data = data

     outdata&storm = storm
     outdata&time = time
     outdata&data = data
     outdata&variables = variables

     storm_pos =  outdata(:,:,0,0)
     delete(storm_pos@data)
     delete(storm_pos@variables)

     storm_sn = new( (/1, 13/), character)
     storm_sn!0 = "storm"
     storm_sn&storm = storm
     storm_sn!1 = "ncharsn"
     storm_sn&ncharsn = ispan(0,12,1) 
     storm_sn@_FillValue = default_fillvalue("character")

     c_inds!0 = "center"
     c_inds&center = c_inds

     source_pos = outdata(:,:,0,:dimsizes(c_inds)-1)
     source_pos!2 = "center"
     source_pos&center = c_inds
     delete(source_pos@data)
     source_pos@center = centers(c_inds)

     ncdf->data      = data        
     ncdf->variables = variables 
     ncdf->time = time
     ncdf->center = c_inds
     ncdf->storm_data  = outdata
     ncdf->storm_time = storm_time
     ncdf->storm = storm_id
     ncdf->storm_lat = storm_pos
     ncdf->storm_lon =storm_pos
     ncdf->storm_sn =  storm_sn 
end if

;setfileoption("nc","Format", "NetCDF4") 
;ncdf = addfile(out_filename, "w")

ncin = addfile("~/data/ibtracs/Allstorms.ibtracs_all.v03r08.nc","r")

;source_lat := short2flt(ncin->source_lat(:,:,c_inds))
;source_lat@souce = centers(c_inds)
;source_lon := short2flt(ncin->source_lon(:,:,c_inds))
;source_lon@souce = centers(c_inds)
;lat := dim_avg_n_Wrap( source_lat,2)
;lon := dim_avg_n_Wrap( source_lon, 2)

source_lat := short2flt(ncin->lat_for_mapping(:,:))
source_lon := short2flt(ncin->lon_for_mapping(:,:))
lat := source_lat
lon := source_lon



nmi := ind( dim_num_n(.not.ismissing(lat),1) .ne. 0)
times := ncin->source_time(nmi,0)
year := cd_calendar(times, 0)
storms := nmi( ind(year(:,0).ge.1979))

storm_id := ncin->storm_sn(storms,:)
printVarSummary(storm_id)

rlon = lon
rlat = lat

;out_ct = dimsizes(ncdf->storm_id)
;if(out_ct .eq. 1)
;   out_ct = 0
;end if

;ss = dimsizes(ncdf->storm_id)
do s=0,dimsizes(storms)-1
progress_bar(s, dimsizes(storms))

time := cd_convert(ncin->source_time(storms(s),:), "hours since 1800-01-01 00:00:00")

outdata := new( (/dimsizes(time),3,13/) , float)
outdata!0 = "time"
outdata!1 = "data"
outdata!2 = "variables"
outlocs := new( (/dimsizes(time),2/) , float)
outlocs!0 = "time"
outlocs!1 = "data"

lon(storms(s), :) = unwrap_lon(lon(storms(s),:), 6)
lon(storms(s), :)  = linmsg(lon(storms(s), :) , 0)
lat(storms(s), :)  = linmsg(lat(storms(s), :) , 0)
lon(storms(s), :) = wrap_lon(lon(storms(s),:))

do t=0, dimsizes(time)-1
  in_lat := lat(storms(s),t)
  in_lon := lon(storms(s),t)
  in_time := time(t)
;  print( in_lat +"   "+in_lon)
 if(.not.ismissing(in_lat))
    new_loc := rm_single_dims(get_mand_positions(in_time, in_lat, in_lon, 350.))
    outdata(t,:,:7) = (/ new_loc /)
    c_loc = combine_mand_positions(new_loc)
    c_loc = where(abs(c_loc) .gt. 720, c_loc@_FillValue, c_loc)
;    print( c_loc(0)+"   "+c_loc(1) )
   if(.not.ismissing(c_loc(0)))
      outlocs(t,:) = (/c_loc/)
;    print(  lat(storms(s),t+1)+"   "+lon(storms(s),t+1)  )
    if(t.ne.dimsizes(time)-1)
     lat(storms(s),t+1:) = lat(storms(s),t+1:) + ( c_loc(0) - lat(storms(s),t) )
     lon(storms(s),t+1:) = lon(storms(s),t+1:) + ( c_loc(1) - lon(storms(s),t) )
;    print(  lat(storms(s),t+1)+"   "+lon(storms(s),t+1)  )
    end if
    end if
;   check_location(wks, in_time, (/rlat(storms(s),t),rlon(storms(s),t)/), outlocs(t,:) )
  else
  if(all(ismissing(  lat(storms(s),t:))))
   break
  end if
  end if
end do

ncdf->storm(s) = (/storms(s)/)
ncdf->storm_data(s,:dimsizes(time)-1,:,:) = (/outdata/)
ncdf->storm_time(s,:dimsizes(time)-1) = (/time/)
ncdf->storm_lat(s,:dimsizes(time)-1) = (/outlocs(:,0)/)
ncdf->storm_lon(s,:dimsizes(time)-1) = (/outlocs(:,1)/)
ncdf->storm_sn(s,:) = (/ storm_id(s,:) /)
end do
end





undef("combine_track_points")
function combine_track_points()
local inf, sdata, alllats, alllons, nmi, dist, table, dists, v, vv
begin

inf = addfile("test.nc","r")

sdata = inf->storm_data
storm_id  = inf->storm_id
storm_no = storm_id/10000
copy_VarCoords(storm_id, storm_no)
storm_year = storm_id - storm_no*10000
copy_VarCoords(storm_id, storm_year)

storm_time = inf->storm_time
storm_time@units = "hours since 1800-01-01 00:00:00"

outlocs := sdata(:,:,:1,0)
outlocs = outlocs@_FillValue
outlocs!0 = "storm"
outlocs!1 = "time"
outlocs&time = sdata&time
do s=0, dimsizes(sdata&storm)-1
   outlocs(s,:,:) = (/ combine_mand_positions(sdata(s,:,:,:7)) /)
end do
tmp := get_res_value(outlocs, "variables", False)   ; little trick to strip attributes. 
tmp := get_res_value(outlocs, "vals", False)
delete(tmp)

hurlocs = outlocs
hurlocs = hurlocs@_FillValue
do s=0, dimsizes(hurlocs&storm)-1
  progress_bar(s, dimsizes(hurlocs&storm) )
  si := HURDAT_storm_info(storm_year(s), storm_no(s))
  nsi = dimsizes(si(:,0))
  hurlocs(s,:nsi-1,0) = string_latlon_to_numeric(si(:,4))
  hurlocs(s,:nsi-1,1) = string_latlon_to_numeric(si(:,5))
end do

hurlocs!0 = "storm"
hurlocs!1 = "time"


IPE = (/ outlocs(:,:,0) /)
do s=0, dimsizes(hurlocs&storm)-1
  IPE(s,:) = gc_latlon(hurlocs(s,:,0), hurlocs(s,:,1), outlocs(s,:,0), outlocs(s,:,1), 0, 4)
end do

IPE = where(ismissing(hurlocs(:,:,0)), IPE@_FillValue, IPE)
IPE = where(ismissing(outlocs(:,:,0)), IPE@_FillValue, IPE)

IPE!0 = "storm"
IPE!1 = "time"
delete([/IPE@spacing, IPE@gclat, IPE@gclon/])


outlocs&storm = ispan(0, dimsizes(outlocs&storm)-1, 1)
ncout = addfile("erai_storm_locations.nc","c")
ncout->cfsr_lat = outlocs(:,:,0)
ncout->cfsr_lon = outlocs(:,:,1)
ncout->hurdat_lat = hurlocs(:,:,0)
ncout->hurdat_lon = hurlocs(:,:,1)
ncout->ipe = IPE
ncout->storm_id = storm_id
ncout->storm_no = storm_no
ncout->storm_year = storm_year
ncout->storm_time = storm_time

end



undef("extra_code")
procedure extra_code()
begin

hist1D(IPE,0)

IPE_max = dim_max_n(IPE,1)

p_i = ind(IPE_max .gt. 500)

wks = gsn_open_wks("X11","")

mpres = True
mpres@gsnDraw =False
mpres@gsnFrame = False
do p=0, dimsizes(p_i)-1
  i = p_i(p)
  mpres@mpMaxLatF := max(outlocs(i,:,0))+5
  mpres@mpMinLatF := min(outlocs(i,:,0))-5
  mpres@mpMaxLonF := max(outlocs(i,:,1))+5
  mpres@mpMinLonF := min(outlocs(i,:,1))-5
  mpres@gsnCenterString = storm_id(i)+"    "+IPE_max(i)
  map = gsn_csm_map(wks, mpres)
  track_hur = gsn_add_polyline(wks, map, hurlocs(i,:,1),hurlocs(i,:,0), False)
  track_re = gsn_add_polyline(wks, map, outlocs(i,:,1),outlocs(i,:,0), False)
  df(wks, map)
end do



alllats := reshape(sdata(:,:,0,:7), (/556*200,8/) )
alllons := reshape(sdata(:,:,1,:7), (/556*200,8/) )

nmi := ind(dim_num_n(ismissing(alllats),1).ne.6)

alllats := alllats(nmi,:)
alllons := alllons(nmi,:)


dist = gc_latlon(alllats(:,0), alllons(:,0), alllats(:,1), alllons(:,1),0,4)

table := new( (/8,8/), float)
dists := new( (/8,8/), float)

do v=0,7
  nmi := ind(.not.ismissing(alllats(:,v)))
  do vv=v,7
  nnmi := ind(.not.ismissing(alllons(nmi,vv)))
  dist := gc_latlon( alllats(nmi(nnmi),v), alllons(nmi(nnmi),v),  alllats(nmi(nnmi),vv), alllons(nmi(nnmi),vv),0, 4)
  dists(v,vv) = avg(dist)
  end do
  table(v,:) = avg( dim_num_n(ismissing(alllons(nmi,:)),1) )
end do

dists@_FillValue = -1

print_r(dists)
print(table(:,0))
end
;
;data := combine_mand_positions(outdata)
;
;
;
;
;obs = (/lat(storms(s),t),  lon(storms(s),t) /)
;
;wks = gsn_open_wks("X11","")
;res = True
;res@gsnMaximize = True
;res@gsnDraw = False
;res@gsnFrame = False
;res@cnLinesOn = True
;pres = True
;pres@gsnDraw = False
;pres@gsnFrame = False
;pres@cnLevelSpacingF = 50.
;pres@cnLineLabelsOn = True
;pres@cnInfoLabelOn = True
;pres@cnHighLabelsOn  = True
;pres@cnLowLabelsOn = True
;pres@cnFillOn = False
;pres@cnLabelMasking = True
;pres@cnLineLabelBackgroundColor = -1
;pres@gsnAddCyclic = True
;
;;
;gres = True
;gres@gsnDraw = False
;gres@gsnFrame = False
;gres@cnLevelSpacingF = 10.
;gres@cnLineColor = "Brown"
;gres@cnLineLabelsOn = True
;gres@cnInfoLabelOn = False
;gres@cnHighLabelsOn  = True
;gres@cnLowLabelsOn = True
;gres@cnLabelMasking = True
;gres@cnLineLabelBackgroundColor = -1
;gres@gsnAddCyclic = True
;vr_res = True
;
;vr_res@gsnAddCyclic = True
;vr_res@gsnDraw = False
;vr_res@gsnFrame = False
;vr_res@cnFillOn = True
;vr_res@cnFillMode = "RasterFill"
;vr_res@cnLevelSelectionMode = "ManualLevels"
;vr_res@cnMinLevelValF = 0.5*10^-5
;vr_res@cnMaxLevelValF = 5*10^-5
;vr_res@cnLevelSpacingF = 0.5*10^-5
;vr_res@cnLinesOn = False
;vr_res@cnRasterSmoothingOn = True
;vr_res@lbOrientation = "vertical"
;ores = True
;ores@gsMarkerIndex = 16
;mres = True
;mres@gsMarkerIndex = 10
;ares = True
;ares@gsMarkerIndex = 0
;ares@txFontHeightF = 0.012
;
;t=25
;obs = (/lat(storms(s),t),  lon(storms(s),t) /)
;;model = (/data(t,0:1) /)
;in_time = time(t)
;;
;;;undef("check_location")
;;;procedure check_location(time, obs, model)
;;;begin
;;;
;mslp  =grab_data(in_time, "pmsl", False)
;opt = True
;opt@levs = 850
;g850  =grab_data(in_time, "g", opt)
;vr850  =grab_data(in_time, "vr", opt)
;vr850 = vr850*-1
;res@mpMinLatF = obs(0)-10
;res@mpMaxLatF = obs(0)+10
;res@mpMinLonF = obs(1)-10
;res@mpMaxLonF = obs(1)+10
;
;map = gsn_csm_map(wks,res)
;vlot = gsn_csm_contour(wks, vr850, vr_res)
;plot = gsn_csm_contour(wks, mslp, pres)
;glot = gsn_csm_contour(wks, g850, gres)
;
;overlay(map, vlot)
;
;overlay(map, plot)
;overlay(map, glot)
;;
;
;
;map@obs = gsn_add_polymarker(wks, map, obs(1), obs(0), ores)
;map@mod = gsn_add_polymarker(wks, map, locs(t,1), locs(t,0), mres)
;map@allmodel := gsn_add_text(wks, map, ispan(1,8,1), outdata(t,1,:7), outdata(t,0,:7), ares)
;
;df(wks, map)
;
;;print(outdata(t,:,2))
;;newdata = (get_mass_center(g850, obs))
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
;;map@newmodel = gsn_add_polymarker(wks, map, newdata(0,1), newdata(0,0), mres)
;
;;hdt_line = gsn_add_polyline(wks, plot, lon, lat, False)
;;cfs_line = gsn_add_polyline(wks,map, dim_avg_n(outdata(:,1,:), 1), dim_avg_n(outdata(:,0,:), 1), False)
;;
;df(wks, map)
;;
;
;
;
;
;if(isdefined("vals"))
;  delete([/vals, names/] )
;end if
;do year=2012, 2014
;  print(year)
;  storms := HURDAT_gen_info(year)
;  if(isdefined("vals"))
;  cnt = dimsizes(vals(:,0,0))
;  vals := array_append_record(vals,  new((/dimsizes(storms&time),7,5/), float),0 )
;  names := array_append_record(names, new(dimsizes(storms&time), "string"),0 )
;  names!0 = "time"
;  else
;  cnt = 0
;  vals = new((/dimsizes(storms&time),7,5/), float)
;  names := new(dimsizes(storms&time), "string")
;  names!0 = "time"
;  end if
;  names(cnt:) = storms(1,:)
;  names&time(cnt:) = storms&time
;
;  do s=0, dimsizes(storms&time)-1
;    storm_info := HURDAT_storm_info(year, storms(1,s))
;    ts= min(ind(storm_info(:,3).eq."TS"))
;    if(ismissing(ts))
;    continue
;    end if
;;    do tt=0,6
;    tt=3
;    t = (tt-3)+ts
;    if(t .lt. 0 .or. t.ge. dimsizes(storm_info&time))
;      continue
;    end if
;
;    lat := string_latlon_to_numeric(storm_info(:,4))
;    lon := string_latlon_to_numeric(storm_info(:,5))
;
;    pres = grab_cfsr(storm_info&time(t), "pmsl", False)
;    pres = pres/100.
;    
;    opt = True
;    opt@levs = (/200,850/)
;    opt@levSelection = "Explicit"
;    hgt = grab_cfsr(storm_info&time(t), "g", opt)
;    hgtd = hgt(0,:,:) - hgt(1,:,:)
;    copy_VarCoords(hgt(0,:,:), hgtd)
;
;    opt = True
;    opt@levs = (/700,850/)
;    opt@levSelection = "Range"
;    u := dim_avg_n_Wrap(grab_cfsr(storm_info&time(t), "u", opt),0)
;    v := dim_avg_n_Wrap(grab_cfsr(storm_info&time(t), "v", opt),0)
;    vr = uv2vr_cfd(u,v, u&lat, u&lon,1)
;    copy_VarCoords(u,vr)
;    vr = circle_avg_m(vr, 200.)
;    vr = vr*10^5
;    
;    vr_max := is_maxima(vr, lat(t),lon(t), 250.)
;
;    mslp_min := is_minima(pres, lat(t),lon(t), 250.)
;    if(mslp_min)
;    crit := mslp_min@val
;    crit = crit+2.
;    crit@loop = True
;    crit@interval = 2.
;    mslp_closed := check_closed_contour(pres, (/lat(t), lon(t)/),crit)
;    else
;    mslp_closed := False
;    end if
;    thick_max := is_maxima(hgtd, lat(t), lon(t), 250.)
;    thick_max := is_env_maxima(hgtd, lat(t),lon(t), 250.)
;
;    
;;    mslp_min = (/ mslp_closed /)
;  ;  if( vr_max@val*10^5 .lt. 4) vr_max = False end if
;  ;  if(thick_max@val .lt. 20) thick_max = False end if
;  ;  if(mslp_min@val/100. .gt. 1020) mslp_min = False end if 
;  
;    lat_pos = (/get_res_value_keep(mslp_min, "lat", -999), \
;                get_res_value_keep(vr_max, "lat", -999),\
;                get_res_value_keep(thick_max, "lat", -999) /)
;    lon_pos = (/get_res_value_keep(mslp_min, "lon", -999), \
;                get_res_value_keep(vr_max, "lon", -999),\
;                get_res_value_keep(thick_max, "lon", -999) /)
;              
;    lat_pos@_FillValue = -999
;    lon_pos@_FillValue = -999   
;    avg_lat =  avg(lat_pos)
;    avg_lon = avg(lon_pos)
;  
;    vals(cnt+s,tt,:) = (/get_res_value_keep(mslp_min, "val", vals@_FillValue), \
;                get_res_value_keep(vr_max, "val", vals@_FillValue),\
;                get_res_value_keep(thick_max, "val", vals@_FillValue),\
;                get_res_value_keep(mslp_closed, "pres", vals@_FillValue),\
;                gc_latlon(lat(t), lon(t), avg_lat, avg_lon,0, 4)  /)
;    end do
;  end do
;end do
;
;ni := ind( dim_num_n( (/ vals(:,3,0).lt. 1012,.not.(ismissing(vals(:,3,3))), vals(:,3,1).ge.8, vals(:,3,2).gt.45 /) ,0).lt.2 )
;print(dim_num_n(dim_num_n( (/ vals(:,6,0).lt. 1012.,.not.(ismissing(vals(:,6,3))), vals(:,6,1).ge.8, vals(:,6,2).gt.45 /) ,0).lt.3 ,0))
;
;
;end 
