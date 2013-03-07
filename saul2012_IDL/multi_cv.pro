pro multi_cv,i,j,dataf,cv=cv,outfile=outfile,nosave=nosave,ra=ra,dec=dec

;
; MULTI_CV
;
; Convolves DAta Cube with multiple kernels, and combines them
;
;
;
;
;

;;;; Define ;;;;
five_three=2.16263
ten_three=8.6476
fif_three=19.2923

extract_coords,dataf,ra,dec,vels,fits

best_cv=fltarr(512,512,2048)

cv_all_2d,i,j,dataf=dataf,fwhm=3.4,ra=ra,dec=dec,cv=cv,fits=fits,vels=vels,/nosave

best_cv=cv


cv_all_2d,i,j,dataf=dataf,fwhm=5,ra=ra,dec=dec,cv=cv,fits=fits,vels=vels,/nosave

cv/=five_three

best_cv=best_cv>cv

cv_all_2d,i,j,dataf=dataf,fwhm=10,ra=ra,dec=dec,cv=cv,fits=fits,vels=vels,/nosave

cv/=ten_three

best_cv=best_cv>cv


cv_all_2d,i,j,dataf=dataf,fwhm=15,ra=ra,dec=dec,cv=cv,fits=fits,vels=vels,/nosave


cv/=fif_three

cv=best_cv>cv


if ~keyword_set(nosave) then begin
if ~keyword_set(outfile) then begin
save, dec, ra, cv, f= '~/cv/cv_' + string(i, f='(I2.2)') +'_' + string(j, f='(I2.2)')+ '.sav'
endif else begin
save,dec,ra,cv,f=outfile
endelse
endif

end
