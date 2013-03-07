pro cv_plots,infile,outfile

;
; Make plots of the convolved cubes' properties as a function of v
;
;

restore,infile
sz=size(cv)
stds=fltarr(sz[3])
means=fltarr(sz[3])
meds=median(median(cv,dim=1),dim=1)
maxs=max(max(cv,dim=1,/nan),dim=1,/nan)
for i=0,sz[3]-1 do begin
stds[i]=stddev(cv[*,*,i],/nan)
means[i]=mean(cv[*,*,i],/nan)
endfor

save,meds,maxs,means,stds,file=outfile
end
