undef("time_since_change")
function time_since_change(filename)
local nowtime, filetime
begin
if(.not. fileexists(filename) )
    return 0
else    
    nowtime = toint(systemfunc("date +%s"))
    filetime = toint(systemfunc("stat -c%Y "+filename))
    return nowtime - filetime
end if
end

undef("isfilevarpresent")
function isfilevarpresent(filli, var, ndim,lev)
local inf, vars
begin

  if( time_since_change(filli) .lt. 120) return False end if
  if( .not. isfilepresent(filli) )  return False end if
  inf = addfile(filli+".grib2", "r")
  vars = getfilevarnames(inf)
  if( .not. any( vars .eq. var) )  return False end if
  ndims = dimsizes(getfilevardimsizes(inf, var)) 
  if(ndims .ne. ndim )
   return False
  end if
  if(ndims .eq. 3)
  	vardims = getfilevardims(inf, var)
  	levs = inf->$vardims(0)$
  	overlap = get1Dindex(levs, lev)
  	if(any(ismissing(overlap)))
  	return False
  	end if
  end if
  ;  If you've got here it should be true? 
  return True
end




undef("grab_data")
function grab_data(times, var, opt)
local times, var, opt, var_short, gribnames, vars, dims, ens, fhr, lev, tmp, ret_dims, retval, filli, inf
begin

lat = "lat_0"
lon = "lon_0"

if(var.eq."dv")
 u = grab_data(times, "u", opt)
 v = grab_data(times, "v", opt)
 varin = uv2dv_cfd(u,v,u&lat,u&lon,3)
 copy_VarMeta(u, varin)
 varin@units = "s~S~-1"
 varin@long_name = "Divergence"
 return(varin)
end if
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

var_short = (/"g","precip","mslp","u","v","tpw","q", "t"/)
vi = ind(var_short.eq.var)
if(ismissing(vi))
  print("variable name not found options are:")
  print("     "+var_short)
  return(default_fillvalue("float"))
end if


ens = get_res_value_keep(opt, "ens", "c00")
fhr = get_res_value_keep(opt, "fhr", 0 )
lev = get_res_value_keep(opt, "levs", (/100,1000/) )
lev = lev*100

lev = where(lev .eq. 80000, 85000, lev)
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
  if(ens .eq. "d00" )
      dir = "/cas2/unidata/GRIB/gfs1p0/"
      filli = dir+cd_string(times, "GFS_%H_%y%N%D_")+sprinti("%03i",toint(fhr))
      gribnames = (/"HGT_P0_L100_GLL0", "APCP_P11_L1_GLL0_acc6h", "PRMSL_P0_L101_GLL0", "UGRD_P0_L100_GLL0", "VGRD_P0_L100_GLL0","PWAT_P0_L200_GLL0","RH_P1_L100_GLL0_1", "TMP_P0_L100_GLL0"/)
  else
      dir = "/cas2/unidata/GRIB/gefs/"
      filli = dir+cd_string(times, "GFS_%y%N%D_"+ens+"_%H_")+sprinti("%03i",toint(fhr))
      gribnames = (/"HGT_P1_L100_GLL0", "APCP_P11_L1_GLL0_acc6h", "PRMSL_P1_L101_GLL0", "UGRD_P1_L100_GLL0", "VGRD_P1_L100_GLL0","PWAT_P1_L200_GLL0","RH_P1_L100_GLL0_1", "TMP_P1_L100_GLL0"/)

  end if
    ;print(filli)
  vars = (/"hgt_pres",  "apcp_sfc" ,"pres_msl", "ugrd_pres", "vgrd_pres", "tpw","t"/)
  dims= (/4,3,3,4,4,3,4,4/)-1

  ct = 0
  do while( .not. isfilevarpresent(filli, gribnames(vi), dims(vi),lev )  )      ; we've been looping for an hour.  Give up and accept defeat. 
        print(filli)
        print(addfile(filli,"r"))
        print("sleeping")
        sleep(60)
        ct = ct+1
        if( ct .gt. 60)
            print("ERROR 1:  File hasn't been found for an hour. ")
            return retval
        end if
  end do 
  

  inf = addfile(filli+".grib2", "r")
  
  if(dims(vi).eq.3) 
       if(isatt(opt,"levSelection") .and.(opt@levSelection .eq. "Range" ))
          retval = lonFlip(inf->$gribnames(vi)$({min(lev):max(lev)}, ::-1, :))
       else
          retval = lonFlip(inf->$gribnames(vi)$({lev}, ::-1, :))
       end if
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

