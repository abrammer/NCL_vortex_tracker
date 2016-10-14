


undef("grab_data")
function grab_data(times, var, opt)
local times, var, opt, var_short, gribnames, vars, dims, ens, fhr, lev, tmp, ret_dims, retval, filli, inf
begin

lat = "lat_0"
lon = "lon_0"

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

var_short = (/"g","precip","mslp","u","v"/)
vi = ind(var_short.eq.var)
if(ismissing(vi))
  print("variable name not found options are:")
  print("     "+var_short)
  return(default_fillvalue("float"))
end if

gribnames = (/"HGT_P1_L100_GLL0", "APCP_P11_L1_GLL0_acc6h", "PRMSL_P1_L101_GLL0", "UGRD_P1_L100_GLL0", "VGRD_P1_L100_GLL0"/)
vars = (/"hgt_pres",  "apcp_sfc" ,"pres_msl", "ugrd_pres", "vgrd_pres"/)
dims= (/4,3,3,4,4/)-1

ens = get_res_value_keep(opt, "ens", "c00")
fhr = get_res_value_keep(opt, "fhr", 0 )
lev = get_res_value_keep(opt, "levs", (/100,1000/) )
lev = lev*100

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


dir = "/cas2/unidata/GRIB/gefs/"
filli = dir+cd_string(times, "GFS_%y%N%D_"+ens+"_%H_")+sprinti("%03i",toint(fhr))+".grib2"

;  filli = cd_string(times, "/brammerlab_rit/gefs/%Y%N/%Y%N%D%H/"+ens+"/"+vars(vi)+"_%Y%N%D%H_"+ens+".grib2")
  inf = addfile(filli, "r")
  if(dims(vi).eq.3) 
  retval = lonFlip(inf->$gribnames(vi)$({min(lev):max(lev)}, ::-1, :))
  end if
  if(dims(vi).eq.2) 
  retval = lonFlip(inf->$gribnames(vi)$( ::-1, :))
  end if
  retval := rm_single_dims(retval)
  ndims = dimsizes(dimsizes(retval))
  retval!(ndims-2) = "lat"
  retval!(ndims-1) = "lon"
end if


return(retval)
end

