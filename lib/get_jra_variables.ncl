load "/home/abrammer/NCL/functions/uv2cvF_Wrap.ncl"
load "/home/abrammer/NCL/functions/uv2owF_Wrap.ncl"

undef("JRAhyb_to_plev")
function JRAhyb_to_plev(Ps)
local A, B, retval
begin
 A = (/ 0, 0, 0, 0, 0, 0, 0, 0, 133.0510113, 364.9041489, 634.6027164, 959.7971673, 1347.680042, 1790.907396, 2294.84169, 2847.484778, 3468.871488, 4162.956463, 4891.880833, 5671.82424, 6476.712996, 7297.469895, 8122.159791, 8914.082201, 9656.181911, 10329.43618, 10912.63844, 11369.64783, 11695.3716, 11861.25309, 11855.43432, 11663.35537, 11285.40406, 10729.94941, 10014.61505, 9167.247036, 8226.244908, 7201.56898, 6088.673009, 4950, 4000, 3230, 2610, 2105, 1700, 1370, 1105, 893, 720, 581, 469, 377, 301, 237, 182, 136, 97, 65, 39, 20, 0 /) 
 B = (/ 1, 0.997, 0.994, 0.989, 0.982, 0.972, 0.96, 0.946, 0.92666949, 0.904350959, 0.879653973, 0.851402028, 0.8195232, 0.785090926, 0.748051583, 0.709525152, 0.668311285, 0.624370435, 0.580081192, 0.534281758, 0.48823287, 0.442025301, 0.395778402, 0.350859178, 0.307438181, 0.265705638, 0.225873616, 0.189303522, 0.155046284, 0.124387469, 0.096445657, 0.072366446, 0.052145959, 0.035700506, 0.022853849, 0.01332753, 0.006737551, 0.00248431, 0.00011327, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)

retval = pres_hybrid_jra55(Ps, A, B)
copy_VarCoords(Ps, retval(0,:,:))
retval!0 = "lev"
retval&lev = ispan(1,60,1)
 return(retval)
end




undef("grab_data")
function grab_data(time_in, var, opt)
local cf, climo_dir, climo_suff, doy, dims, varin, time, var, opt,u,v,file_dir,file_suff,cdtime,ff, u, v 
begin
time = time_in
sort_ind = dim_pqsort(time,2)
unsort_ind = get1Dindex(sort_ind, ispan(0,dimsizes(sort_ind)-1,1))


if(var.eq."vr")
 u = grab_data(time_in, "u", opt)
 v = grab_data(time_in, "v", opt)
 varin = uv2vr_cfd(u,v,u&lat,u&lon,3)
 copy_VarMeta(u, varin)
 varin@units = "s~S~-1"
 varin@long_name = "Relative Vorticity"
 return(varin)
end if
if(var.eq."uv")
  u = grab_data(time_in, "u", opt)
  v = grab_data(time_in, "v", opt)
  varin = sqrt(u^2+v^2)
  copy_VarMeta(u, varin)
  varin@long_name = "Wind Speed"
  return(varin)
end if
if(var.eq."cv")
     u = grab_data(time_in, "u", opt)
     v = grab_data(time_in, "v", opt)
     varin = uv2cvF_Wrap(u,v)
     return(varin)
end if
if(var.eq."ow")
     u = grab_data(time_in, "u", opt)
     v = grab_data(time_in, "v", opt)
     varin = uv2owF_Wrap(u,v)
     return(varin)
end if
if(var.eq."mslp")
    varin = grab_data(time_in, "pmsl",opt)
    return(varin)
end if

gribNames = (/"u-component_of_wind_hybrid", "v-component_of_wind_hybrid", "SLP","Specific_humidity_hybrid","Geopotential_height_hybrid", "Pressure_surface","Temperature_hybrid"/)
vars      = (/"u", "v", "mslp","q", "g", "p_sfc","t"/)
v_ind = ind(var.eq.vars)
if(ismissing(v_ind))
    print("Variable not found, try --")
    print(vars)
end if
var = gribNames(v_ind)


base_url="http://rda.ucar.edu/thredds/dodsC/aggregations/JRA55/ds628.0/"
mdl_dir="7/TP"
sfc_dir="12/TP"
sfc_f = addfile(base_url+sfc_dir,"r")
;time = get_current_time()
;time = time-1000.
time := cd_convert(time, sfc_f->time@units)
time_ind = ind_nearest_coord(time, sfc_f->time, 0)
if( sfc_f->time(time_ind).ne. time)
 print(" Using Nearest time -- check inputs ")
 print(cd_string(time, "%Y%N%D %H%M")+"  --- "+cd_string(sfc_f->time(time_ind), "%Y%N%D %H%M"))
end if


var_f = addfile(base_url+mdl_dir,"r")
varin = var_f->$var$(time_ind,:,:,:)
varin@_FillValue =  default_fillvalue("float")

if(opt .and. isatt(opt, "levs"))
 sfc_pressure = sfc_f->Pressure_surface(time_ind,:,:)
 sfc_pressure@_FillValue = default_fillvalue("float")
 plevs = JRAhyb_to_plev(sfc_pressure)
 replace_ieeenan(plevs, default_fillvalue("float"), 0)
 varin := int2p_n_Wrap(plevs, varin, opt@levs*100., 2,0)
end if
varin := rm_single_dims(varin)
varin = lonFlip(varin)

ndims = dimsizes(dimsizes(varin))
if(dimsizes(time).gt.1)
  if(ndims.eq.3)
      return(varin(unsort_ind,:,:))
  end if
  if(ndims.eq.4)
      return(varin(unsort_ind,:,:,:))
  end if
else
  return(varin)
end if

end

undef("cd_climo_time")
function cd_climo_time(x)
local ut, jan1st, retx
begin
 ut = cd_calendar(x, -5)
 ut(:,1:2) = 1
 ut(:,3:) = 0
 jan1st = cd_inv_calendar(ut(:,0), ut(:,1), ut(:,2), ut(:,3), ut(:,4), ut(:,5),x@units,0)
 retx = x - jan1st
 retx@units = "hours since 0-1-1 00:00:00"
 return(retx)
end



undef("grab_climo")
function grab_climo(time_in, var, opt)
local cf, climo_dir, climo_suff, doy, dims, varin, time, var, opt,u,v,file_dir,file_suff,cdtime,ff, u, v, years
begin
time = time_in
sort_ind = dim_pqsort(time,2)
unsort_ind = get1Dindex(sort_ind, ispan(0,dimsizes(sort_ind)-1,1))

if(var.eq."vr")
 u = grab_climo(time_in, "u", opt)
 v = grab_climo(time_in, "v", opt)
 varin = uv2vr_cfd(u,v,u&lat,u&lon,3)
 copy_VarMeta(u, varin)
 varin@units = "s~S~-1"
 varin@long_name = "Relative Vorticity"
 return(varin)
end if
if(var.eq."uv")
  u = grab_climo(time_in, "u", opt)
  v = grab_climo(time_in, "v", opt)
  varin = sqrt(u^2+v^2)
  copy_VarMeta(u, varin)
  varin@long_name = "Wind Speed"
  return(varin)
end if
if(var.eq."cv")
     u = grab_climo(time_in, "u", opt)
     v = grab_climo(time_in, "v", opt)
     varin = uv2cvF_Wrap(u,v)
     return(varin)
end if
if(var.eq."ow")
     u = grab_climo(time_in, "u", opt)
     v = grab_climo(time_in, "v", opt)
     varin = uv2owF_Wrap(u,v)
     return(varin)
end if
if(var.eq."mslp")
    varin = grab_climo(time_in, "pmsl",opt)
    return(varin)
end if

file_dir = "/erai/"
file_suff = ".nc"
years = "mean"
ff = addfiles(file_dir+years+"/"+var+"."+years+file_suff, "r")

time := cd_climo_time(time)
nt = dimsizes(time)-1

dims = getfilevardims(ff[0],var)
if(any(dims.eq."lev"))
   if(isatt(opt, "levs"))
        if(isatt(opt,"levSelection") .and.(opt@levSelection .eq. "Range" ))
           varin = ff[:]->$var$({time}, {opt@levs(0):opt@levs(1)}, :, :)
        else
           varin = ff[:]->$var$({time}, {opt@levs}, :, :)
        end if
 else
   varin = ff[:]->$var$({time}, :, :, :)
 end if
else
varin = ff[:]->$var$({time}, :, :)
end if

varin = lonFlip(varin)

ndims = dimsizes(dimsizes(varin))
if(dimsizes(time).gt.1)
  if(ndims.eq.3)
      return(varin(unsort_ind,:,:))
  end if
  if(ndims.eq.4)
      return(varin(unsort_ind,:,:,:))
  end if
else
  return(varin)
end if

end



