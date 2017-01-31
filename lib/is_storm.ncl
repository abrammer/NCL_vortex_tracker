
undef("is_storm")
function is_storm( data, thres_values )
local dimz, ndimz, _data, data, thres_values, logic_data, count,i,j,k
begin
dimz =  dimsizes(data)
ndimz = dimsizes(dimz)
if(ndimz .eq. 4)
  rdimz := (/ dimz(0),dimz(1), dimz(ndimz-1) /) 
  retval = new( rdimz, float)
  do i=0, dimz(0)-1
    if(all(ismissing(data(i,:,:,:)))) continue end if
    do j=0, dimz(1)-1
      if(all(ismissing(data(i,j,:,:)))) continue end if
      retval(i,j,:) = is_storm( data(i,j,:,:), thres_values)
    end do
  end do
  return retval
end if
if(ndimz .eq. 3)
  rdimz := (/ dimz(0), dimz(ndimz-1) /) 
  retval = new( rdimz, float)
  do i=0, dimz(0)-1
    if(all(ismissing(data(i,:,:)))) continue end if
    retval(i,:) = is_storm( data(i,:,:), thres_values)
  end do
  return retval
end if

if(ndimz.eq.1)
 dimz := array_append_record(1, dimz, 0)
 _data = data
 data := new(dimz, typeof(_data))
 data(0,:) = _data
 delete(_data)
end if
logic_data = new(dimz, integer)
do i=0, dimz(0)-1
  logic_data(i,:) = where( data(i,:) .ge. thres_values(i), 1, 0) 
end do
logic_data = where(ismissing(logic_data), 0, logic_data)
count = dim_acumrun_n(dim_min_n(logic_data,0),4,1,0)
return count
end
