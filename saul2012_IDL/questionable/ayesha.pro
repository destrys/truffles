pro f2test

restore,'spores_ayesha.sav',/verb
;oldspore=fits

tm=systime(/sec)
;spore,fits,chin,dataf='GALFA_HI_RA+DEC_332.00+18.35_W.fits',outfile='spores_ayesha.sav',ra=ra,dec=dec,vels=vels
print, (systime(/sec)-tm)/60., 'SPORE runtime (minutes)'


tm=systime(/sec)
cv_all_2d, 15, 3,ra=ra,dec=dec,chin=chin,fwhm=7,v_fwhm=5,vels=vels,cv=cv,fits=fits,outfile='cv_ayesha_7.5.sav',wtout='wt_ayesha_7.5.sav'
print, (systime(/sec)-tm)/60., 'CV_ALL runtime (minutes)'


tm=systime(/sec)
cv_all_2d, 15,3, ra=ra,dec=dec,chin=chin,fwhm=18,v_fwhm=5,vels=vels,cv=cv,fits=fits,outfile='cv_ayesha_18.5.sav',wtout='wt_ayesha_18.5.sav'
print, (systime(/sec)-tm)/60., 'CV_ALL runtime (minutes)'


tm=systime(/sec)
cv_all_2d, 15, 3,ra=ra,dec=dec,chin=chin,fwhm=7,v_fwhm=15,vels=vels,cv=cv,fits=fits,outfile='cv_ayesha_7.15.sav',wtout='wt_ayesha_7.15.sav'
print, (systime(/sec)-tm)/60., 'CV_ALL runtime (minutes)'


tm=systime(/sec)
cv_all_2d, 15,3, ra=ra,dec=dec,chin=chin,fwhm=18,v_fwhm=15,vels=vels,cv=cv,fits=fits,outfile='cv_ayesha_18.15.sav',wtout='wt_ayesha_18.15.sav'
print, (systime(/sec)-tm)/60., 'CV_ALL runtime (minutes)'




tm=systime(/sec)
find_reg_new, 15,3,datacv='cv_ayesha_7.5.sav',fits=fits,outfile='rmask_ayesha_7.5.sav'
print, (systime(/sec)-tm)/60., 'Find_reg runtime (minutes)'



tm=systime(/sec)
find_reg_new, 15,3,datacv='cv_ayesha_18.5.sav',fits=fits,outfile='rmask_ayesha_18.5.sav'
print, (systime(/sec)-tm)/60., 'Find_reg runtime (minutes)'



tm=systime(/sec)
find_reg_new, 15,3,datacv='cv_ayesha_7.15.sav',fits=fits,outfile='rmask_ayesha_7.15.sav'
print, (systime(/sec)-tm)/60., 'Find_reg runtime (minutes)'



tm=systime(/sec)
find_reg_new, 15,3,datacv='cv_ayesha_18.15.sav',fits=fits,outfile='rmask_ayesha_18.15.sav'
print, (systime(/sec)-tm)/60., 'Find_reg runtime (minutes)'









stop
end
