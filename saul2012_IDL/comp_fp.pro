q =0
for n=3, 4 do begin
for m=19, 28 do begin
    if not( (m ge 27) and (n eq 4)) then begin 
        restore, '~goldston/cv/pig_fit_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'
        q = q+n_elements(pig_fit)
    endif
    endfor
endfor

pf = replicate({gfit:fltarr(6), sigs:fltarr(6), twodfit:fltarr(9), kp:0.}, q)
q=0
for n=3, 4 do begin
for m=19, 28 do begin
    if not( (m ge 27) and (n eq 4)) then begin 
        restore, '~goldston/cv/pig_fit_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'
        np = n_elements(pig_fit)
        pf[q:q+np-1].gfit = pig_fit.gfit
        pf[q:q+np-1].sigs = pig_fit.sigs
        pf[q:q+np-1].twodfit = pig_fit.twodfit
        if ((n eq 3) and (m eq 24)) then pf[q:q+2].kp = 1.
        q = q+n_elements(pig_fit)
        endif
    endfor
endfor

save, pf, f='pf.sav'

end
