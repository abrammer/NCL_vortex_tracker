undef("progress_bar")
procedure progress_bar(x, n)
local xf
begin

str = get_res_value_keep(n, "str", "")
xf =x + 1.   ;; convert to float and offset for 0 index
system("tput cuu1 && echo '"+sprintf("% 5.1f", 100*xf/n)+"% Completed" +" "+str+"'")

end

