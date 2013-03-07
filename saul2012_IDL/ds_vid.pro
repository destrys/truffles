pro ds_vid,cv
for i=0,2000 do begin
tv,bytscl(cv[*,*,i],min=0,max=30)
wait,0.05
endfor


end
