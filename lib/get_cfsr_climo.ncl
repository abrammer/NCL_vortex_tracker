load "/home/abrammer/NCL/functions/uv2cvF_Wrap.ncl"
load "/home/abrammer/NCL/functions/uv2owF_Wrap.ncl"



undef("grab_data")
function grab_data(time_in, var, opt)
local cf, climo_dir, climo_suff, doy, dims, varin, time, var, opt,u,v,file_dir,file_suff,cdtime,ff, u, v, years
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

file_dir = "/cfsr/data/"
file_suff = ".0p5.anl.nc"
cdtime = cd_calendar(time, 0)
years = get_unique_values(cdtime(:,0))
ff = addfiles(file_dir+years+"/"+var+"."+years+file_suff, "r")

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



