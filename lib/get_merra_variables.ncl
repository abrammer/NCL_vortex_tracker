load "/home/abrammer/NCL/functions/uv2cvF_Wrap.ncl"
load "/home/abrammer/NCL/functions/uv2owF_Wrap.ncl"



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
    print("calculating CV")
     u = grab_data(time_in, "u", opt)
     v = grab_data(time_in, "v", opt)
     pvs(u)
     pvs(v)
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

gribNames = (/"U", "V", "SLP","QV","H", "PS","T"/)
vars      = (/"u", "v", "mslp","q", "g", "p_sfc","t"/)
v_ind = ind(var.eq.vars)
if(ismissing(v_ind))
    print("Variable not found, try --")
    print(vars)
end if
var_ = gribNames(v_ind)


file_dir = "https://goldsmr5.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I6NPANA.5.12.4"
streams = (/100, 200, 300, 400/)
stream_starts = (/ \t
    cd_inv_calendar(1979,01,01,00,00,00,time@units,0),   \
    cd_inv_calendar(1992,01,01,00,00,00,time@units,0),   \
    cd_inv_calendar(2002,01,01,00,00,00,time@units,0),   \
    cd_inv_calendar(2011,01,01,00,00,00,time@units,0) /)
stream_ends = (/ \
    cd_inv_calendar(1992,01,01,00,00,00,time@units,0),   \
    cd_inv_calendar(2002,01,01,00,00,00,time@units,0),   \
    cd_inv_calendar(2011,01,01,00,00,00,time@units,0),   \
    cd_inv_calendar(2020,01,01,00,00,00,time@units,0) /) \
     -6

s_ind = ind( time .ge. stream_starts .and. time.le.stream_ends )    
if(ismissing(s_ind))
  print("Stream not found -- is it 2020?")
  print("update the script to support stream 500")
end if
stream = streams(s_ind)    
fname = cd_string(time, "/%Y/%N/MERRA2_"+stream+".inst6_3d_ana_Np.%Y%N%D.nc4")

url = file_dir+fname
ff = addfile(file_dir+fname, "r")
;print(ff)

nt = dimsizes(time)-1


time = cd_convert(time, ff->time@units)

dims = getfilevardims(ff,var_)
if(any(dims.eq."lev"))
   if(isatt(opt, "levs"))
        if(isatt(opt,"levSelection") .and.(opt@levSelection .eq. "Range" ))
           varin = ff->$var_$({time}, {opt@levs(0):opt@levs(1)}, :, :)
        else
           varin = ff->$var_$({time}, {opt@levs}, :, :)
        end if
 else
   varin = ff->$var_$({time}, :, :, :)
 end if
else
 varin = ff->$var_$({time}, :, :)
end if

varin&lat := tofloat( varin&lat)
varin&lat@units = "degrees_north"
varin&lon := tofloat( varin&lon)
varin&lon@units = "degrees_east"

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



