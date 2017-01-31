undef("conditional_lonFlip")
function conditional_lonFlip(data)
local dims, ndims, data, retdata
begin
  dims = dimsizes(data)
  ndims= dimsizes(dims)
  data!(ndims-1) = "lon"
  if(min(data&lon).gt.0)
    data&lon = data&lon - 360
  end if
  cyclic = (data&lon(dims(ndims-1)-1) - 360 ) .eq. data&lon(0) 
  if(cyclic)
    if(ndims.eq.4)
     retdata = data(:,:,:,:dims(3)-2)
    end if
    if(ndims.eq.3)
     retdata = data(:,:,:dims(2)-2)
    end if
    if(ndims.eq.2)
     retdata = data(:,:dims(1)-2)
    end if
    if(retdata&lon(0).eq.0) 
      retdata = lonFlip(retdata)
    end if
    return(retdata)
  else
    if(data&lon(0).eq.0) 
      data = lonFlip(retdata)
    end if
    return(data)
  end if
end

undef("make_global")
function make_global(data)
local dims, names, nd, inlats, inlons, retdata, lat, lon, data
begin

dims = dimsizes(data)
names = getvardims(data)
nd = dimsizes(dims)
inlats = data&$names(nd-2)$
inlons = data&$names(nd-1)$

dims(nd-2) = (dims(nd-1)/2)+1

retdata = new(dims, typeof(data))
lat = fspan(-90,90,dims(nd-2)) 
lat@units = "degrees_north"
lon = data&$names(nd-1)$
lon@units = "degrees_east"
retdata!(nd-2) = "lat"
retdata!(nd-1) = "lon"
retdata&lat = lat
retdata&lon = lon

if(nd.eq.2)
   retdata({inlats},:) = (/ data /) 
else
   retdata(:,{inlats},:) = (/data/)
end if
if(min(retdata&lon).gt.0)
retdata&lon = retdata&lon - 360
end if
copy_VarAtts(data, retdata)
copy_VarCoords_not_n(data, retdata, (/nd-2, nd-1/) )
return(retdata)
end


undef("grab_data")
function grab_data(times, var, opt)
local times, var, opt, var_short, gribnames, vars, dims, ens, fhr, lev, tmp, ret_dims, retval, filli, inf
begin

;lat = "lat_0"
;lon = "lon_0"
;
if(var.eq."vr")
 u = grab_data(times, "u", opt)
 v = grab_data(times, "v", opt)
 varin = uv2vr_cfd(u,v,u&lat,u&lon,3)
 copy_VarMeta(u, varin)
 varin@units = "s~S~-1"
 varin@long_name = "Relative Vorticity"
 return(varin)
end if
if(var.eq."uv")
  u = grab_data(times, "u", opt)
  v = grab_data(times, "v", opt)
  varin = sqrt(u^2+v^2)
  copy_VarMeta(u, varin)
  varin@long_name = "Wind Speed"
  return(varin)
end if
if(var.eq."cv")
     u = grab_data(times, "u", opt)
     v = grab_data(times, "v", opt)
     varin = uv2cvF_Wrap(u,v)
     return(varin)
end if
if(var.eq."ow")
     u = grab_data(times, "u", opt)
     v = grab_data(times, "v", opt)
     varin = uv2owF_Wrap(u,v)
     return(varin)
end if

gribnames = (/"T_GDS0_ISBL","GH_GDS0_ISBL", "TP_GDS0_SFC", "MSL_GDS0_SFC", "U_GDS0_ISBL", "V_GDS0_ISBL"/)
vars = (/"tmp_pres","hgt_pres",  "apcp_sfc" ,"pres_msl", "ugrd_pres", "vgrd_pres"/)
dims= (/4,4,3,3,4,4/)
var_short = (/"t", "g","precip","mslp","u","v"/)
vi = ind(var_short.eq.var)
if(ismissing(vi))
  print("variable name not found options are:")
  print("     "+var_short)
  return(default_fillvalue("float"))
end if


ens = get_res_value_keep(opt, "ens", "e00")
fhr = get_res_value_keep(opt, "fhr", (/0/) )
lev = get_res_value_keep(opt, "levs", (/100,1000/) )

if(dimsizes(times).gt.1)
  do t=0, dimsizes(times)-1
    tmp = grab_data(times(t), var, opt)
    if(t.eq.0)
     ret_dims = array_append_record(dimsizes(times), dimsizes(tmp),0)
     retval := new( ret_dims, typeof(tmp) )
     retval!0 = "time"
     retval&time = times

     if(dimsizes(ret_dims).eq.5)  copy_VarMeta(tmp, retval(0,:,:,:,:) ) end if
     if(dimsizes(ret_dims).eq.4)  copy_VarMeta(tmp, retval(0,:,:,:)  ) end if
     if(dimsizes(ret_dims).eq.3)  copy_VarMeta(tmp, retval(0,:,:)  )  end if
    end if
      
     if(dimsizes(ret_dims).eq.5)  retval(t,:,:,:,:) = (/tmp /)  end if
     if(dimsizes(ret_dims).eq.4)  retval(t,:,:,:)   = (/tmp /)  end if
     if(dimsizes(ret_dims).eq.3)  retval(t,:,:)     = (/tmp /)  end if
  end do
else

etype = str_capital(str_get_cols(ens,0,0))
if(etype.eq."C")
  etype = "E"
end if
  eind = toint(str_get_cols(ens,1,2))
  suffix:= where(etype.eq."D".and. fhr.eq.0,"01", "00")
  vtime = times
  vtime = vtime + fhr
  filli = cd_string(times, "/brammerlab_rit/ecmwf/realtime/%Y/E1"+etype+"%N%D%H00")+cd_string(vtime, "%N%D%H")+suffix
  inf := addfile(filli+"1.grib1", "r")
  if(dims(vi).eq.4) 
    if(etype.eq."D")
  retval = inf->$gribnames(vi)$( {min(lev):max(lev)}, ::-1, :)
  retval = conditional_lonFlip(retval)
    else
  retval = conditional_lonFlip(inf->$gribnames(vi)$(eind, {min(lev):max(lev)}, ::-1, :))
    end if
  end if
  if(dims(vi).eq.3) 
    if(etype.eq."D")
  retval = conditional_lonFlip(inf->$gribnames(vi)$( ::-1, :))
    else
  retval = conditional_lonFlip(inf->$gribnames(vi)$(eind, ::-1, :))
    end if
  end if
  retval := rm_single_dims(retval)
  ndims = dimsizes(dimsizes(retval))
  retval!(ndims-2) = "lat"
  retval!(ndims-1) = "lon"
end if

retval := make_global(retval)
return(retval)
end

