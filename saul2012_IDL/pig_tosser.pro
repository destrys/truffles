pro pig_tosser,piggys,pig_out

restore,piggys

goods=fltarr(n_elements(pig_sp))
for i=0,n_elements(pig_sp)-1 do begin
flyingv=where(pig_sp[i].spm EQ 1,ct)
if ((ct GT 1) and (max(flyingv) LT 2040) and (min(flyingv) GT 5)) then goods[i]=1
endfor

pig_out=pig_sp[where(goods EQ 1)]

stop

end

