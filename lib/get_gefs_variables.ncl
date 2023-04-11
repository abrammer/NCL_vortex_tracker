load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "/home/abrammer/NCL/functions/uv2cvF_Wrap.ncl"
load "/home/abrammer/NCL/functions/uv2owF_Wrap.ncl"

undef("grab_data")
function grab_data(times, var, opt)
local times, var, opt, var_short, gribnames, vars, dims, ens, fhr, lev, tmp, ret_dims, retval, filli, inf,\
      t, e, vi, lat_name, lon_name, lev_name,u,v
begin

lat_name = "lat_0"
lon_name = "lon_0"
lev_name = "lv_ISBL0"

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

var_short = (/"g","precip","mslp","u","v","q"/)
vi = ind(var_short.eq.var)
if(ismissing(vi))
  print("variable name not found options are:")
  print("     "+var_short)
  return(default_fillvalue("float"))
end if

gribnames = (/"HGT_P1_L100_GLL0", "APCP_P11_L1_GLL0_acc6h", "PRES_P1_L101_GLL0", "UGRD_P1_L100_GLL0", "VGRD_P1_L100_GLL0","SPFH_P1_L100_GLL0"/)
vars = (/"hgt_pres",  "apcp_sfc" ,"pres_msl", "ugrd_pres", "vgrd_pres", "spfh_pres"/)
dims= (/4,3,3,4,4,4/)

ens = get_res_value_keep(opt, "ens", 0)
fhr = get_res_value_keep(opt, "fhr", (/0,120/) )
lev = get_res_value_keep(opt, "levs", (/100,1000/) )
lev = lev*100

ens := sprinti("p%02i", ens)
ens = where(ens .eq. "p00", "c00",ens)

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
   do e=0,dimsizes(ens)-1
    filli = cd_string(times, "/brammerlab_rit/gefs/%Y%N/%Y%N%D%H/"+ens(e)+"/"+vars(vi)+"_%Y%N%D%H_"+ens(e)+".grib2")
    inf = addfile(filli, "r")
    if(dims(vi).eq.4) 
      tmp = lonFlip(inf->$gribnames(vi)$({min(fhr):max(fhr)},{min(lev):max(lev)}, ::-1, :))
    end if
    if(dims(vi).eq.3) 
      tmp = lonFlip(inf->$gribnames(vi)$({min(fhr):max(fhr)}, ::-1, :))
    end if
    if(e .eq. 0)
      ret_dims = array_append_record(dimsizes(ens), dimsizes(tmp),0)
      retval := new( ret_dims, typeof(tmp) )
      copy_VarAtts(tmp, retval  )
      retval!0 = "ensemble"
      retval&ensemble = ens
      dName = getvardims(tmp)
      oName = dName
      do r=0,dimsizes(dName)-1
        if(dName(r) .eq. lat_name) oName(r) = "lat" end if
        if(dName(r) .eq. lon_name) oName(r) = "lon" end if
        if(dName(r) .eq. lev_name) oName(r) = "lev" end if
        retval!(r+1) = oName(r)
        retval&$oName(r)$  = tmp&$dName(r)$ 
      end do
    end if
    if(dimsizes(ret_dims).eq.5)  retval(e,:,:,:,:)   = (/tmp /)  end if
    if(dimsizes(ret_dims).eq.4)  retval(e,:,:,:)   = (/tmp /)  end if
    if(dimsizes(ret_dims).eq.3)  retval(e,:,:)     = (/tmp /)  end if
   end do
   retval := rm_single_dims(retval)
   ndims = dimsizes(dimsizes(retval))
   retval!(ndims-2) = "lat"
   retval!(ndims-1) = "lon"
end if


return(retval)
end

