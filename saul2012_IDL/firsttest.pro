pro firsttest

restore,'spores_15_03_0.sav',/verb
oldspore=fits

tm=systime(/sec)
;spore,fits,chin,outfile='spores_15_03_0.sav',ra=ra,dec=dec,vels=vels,/nosmooth
print, (systime(/sec)-tm)/60., 'SPORE runtime (minutes)'


tm=systime(/sec)
;cv_all_2d, 15, 3,ra=ra,dec=dec,chin=chin,fwhm=7,v_fwhm=5,vels=vels,cv=cv,fits=fits,outfile='cv_15_03_0_7.5.sav'
print, (systime(/sec)-tm)/60., 'CV_ALL runtime (minutes)'


tm=systime(/sec)
;cv_all_2d, 15,3, ra=ra,dec=dec,chin=chin,fwhm=18,v_fwhm=5,vels=vels,cv=cv,fits=fits,outfile='cv_15_03_0_18.5.sav'
print, (systime(/sec)-tm)/60., 'CV_ALL runtime (minutes)'


tm=systime(/sec)

find_reg, 15,3,datacv='cv_15_03_0_7.5.sav',fits=fits,outfile='rmask_15_03_0_7.5.sav'

print, (systime(/sec)-tm)/60., 'Find_reg runtime (minutes)'



tm=systime(/sec)

find_reg, 15,3,datacv='cv_15_03_0_18.5.sav',fits=fits,outfile='rmask_15_03_0_18.5.sav'

print, (systime(/sec)-tm)/60., 'Find_reg runtime (minutes)'






stop
end
