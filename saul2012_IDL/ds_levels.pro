pro ds_levels,cvfile

;
; Plot mean value of each frequency slice
;
;
;
;
;
;
;
;
;

; Restore CV File
restore,cvfile

; Build Arrays
sz=size(cv)
cv_rms=fltarr(sz[3])
cv_stddev=cv_rms
cv_mean=cv_rms
; Loop, booo.
for i=0,sz[3]-1 do begin
cv_rms[i]=sqrt(mean((cv[*,*,i])^2,/nan))
cv_stddev[i]=stddev(cv[*,*,i],/nan)
cv_mean[i]=mean(cv[*,*,i],/nan)
endfor

!p.multi=[0,1,3]
plot,cv_mean,tit='Mean of each slice'
plot,cv_rms,tit='RMS of each slice'
plot,cv_stddev,tit='STDDEV of each slice'

stop





end
