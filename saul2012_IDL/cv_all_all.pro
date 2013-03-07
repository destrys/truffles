pro cv_all_all

; A wrapper for cv_all_2d to do a ton of convolutions



;cv_all_2d, 15, 3,dataf='~/truffles/cubes/GALFA_HI_RA+DEC_124.00+26.35_W_SPAM_0.fits',fwhm=3.4,outfile='~/cv/cv_15_03_0_3.4.sav'

;stop

;cv_all_2d, 15, 3,dataf='~/truffles/cubes/GALFA_HI_RA+DEC_124.00+26.35_W_SPAM_1.fits',fwhm=3.4,outfile='~/cv/cv_15_03_1_3.4.sav'

;cv_all_2d, 15, 3,dataf='~/truffles/cubes/GALFA_HI_RA+DEC_124.00+26.35_W_SPAM_2.fits',fwhm=3.4,outfile='~/cv/cv_15_03_2_3.4.sav'

;cv_all_2d, 15, 3,dataf='~/truffles/cubes/GALFA_HI_RA+DEC_124.00+26.35_W_SPAM_3.fits',fwhm=3.4,outfile='~/cv/cv_15_03_3_3.4.sav'






for i=1,3 do begin

outf='~/cv_15_03_0_'+strtrim(i*5,1)+'.fits'

cv_all_2d, 15, 3,dataf='~/truffles/cubes/GALFA_HI_RA+DEC_124.00+26.35_W_SPAM_0.fits',fwhm=i*5,outfile=outf

endfor

for i=1,3 do begin

outf='~/cv_15_03_1_'+strtrim(i*5,1)+'.fits'

cv_all_2d, 15, 3,dataf='~/truffles/cubes/GALFA_HI_RA+DEC_124.00+26.35_W_SPAM_1.fits',fwhm=i*5,outfile=outf

endfor

for i=1,3 do begin

outf='~/cv_15_03_2_'+strtrim(i*5,1)+'.fits'

cv_all_2d, 15, 3,dataf='~/truffles/cubes/GALFA_HI_RA+DEC_124.00+26.35_W_SPAM_2.fits',fwhm=i*5,outfile=outf

endfor

for i=1,3 do begin

outf='~/cv_15_03_3_'+strtrim(i*5,1)+'.fits'

cv_all_2d, 15, 3,dataf='~/truffles/cubes/GALFA_HI_RA+DEC_124.00+26.35_W_SPAM_3.fits',fwhm=i*5,outfile=outf

endfor


end

