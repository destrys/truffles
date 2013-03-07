pro truffles_newer
m=20
n=20

mimsy=fltarr(4)

;;;;;; Set up string stuff ;;;;;;
dataf='GALFA_HI_RA+DEC_268.00+10.35_W.fits'

outstring='_268+10'

sporeout='spores'+outstring+'.sav'
spout='pig_sp'+outstring+'.sav'

cvs=['cv'+outstring+'_7.5.sav',$
'cv'+outstring+'_18.5.sav',$
'cv'+outstring+'_7.15.sav',$
'cv'+outstring+'_18.15.sav']
rmasks=['rmask'+outstring+'_7.5.sav',$
'rmask'+outstring+'_18.5.sav',$
'rmask'+outstring+'_7.15.sav',$
'rmask'+outstring+'_18.15.sav']
wts=['wt'+outstring+'_7.5.sav',$
'wt'+outstring+'_7.15.sav',$
'wt'+outstring+'_18.5.sav',$
'wt'+outstring+'_18.15.sav']

save,dataf,outstring,sporeout,spout,cvs,rmasks,wts,f='names'+outstring+'.sav'


;restore,sporeout,/verb

tm=systime(/sec)
;spore,fits,chin,dataf=dataf,outfile=sporeout,ra=ra,dec=dec,vels=vels
print, (systime(/sec)-tm)/60., 'SPORE runtime (minutes)'

mimsy=[[mimsy],[memory()]]
;fits=float(fits)
tm=systime(/sec)
;cv_all_2d,ra=ra,dec=dec,chin=chin,fwhm=7,v_fwhm=5,vels=vels,cv=cv,fits=fits,outfile=cvs[0],wtout=wts[0]
print, (systime(/sec)-tm)/60., 'CV_ALL runtime (minutes)'

mimsy=[[mimsy],[memory()]]


;restore,sporeout,/verb
;fits=float(fits)

tm=systime(/sec)
;cv_all_2d, ra=ra,dec=dec,chin=chin,fwhm=18,v_fwhm=5,vels=vels,cv=cv,fits=fits,outfile=cvs[1],wtout=wts[1]
print, (systime(/sec)-tm)/60., 'CV_ALL runtime (minutes)'

mimsy=[[mimsy],[memory()]]

;restore,sporeout,/verb
;fits=float(fits)

tm=systime(/sec)
;cv_all_2d,ra=ra,dec=dec,chin=chin,fwhm=7,v_fwhm=15,vels=vels,cv=cv,fits=fits,outfile=cvs[2],wtout=wts[2]
print, (systime(/sec)-tm)/60., 'CV_ALL runtime (minutes)'

mimsy=[[mimsy],[memory()]]

;restore,sporeout,/verb
;fits=float(fits)

tm=systime(/sec)
;cv_all_2d, ra=ra,dec=dec,chin=chin,fwhm=18,v_fwhm=15,vels=vels,cv=cv,fits=fits,outfile=cvs[3],wtout=wts[3]
print, (systime(/sec)-tm)/60., 'CV_ALL runtime (minutes)'

mimsy=[[mimsy],[memory()]]

;restore,sporeout,/verb
;fits=float(fits)



tm=systime(/sec)
find_reg_newest, datacv=cvs[0],outfile='test_newer_rmask.sav'
print, (systime(/sec)-tm)/60., 'Find_reg runtime (minutes)'
stop
mimsy=[[mimsy],[memory()]]

tm=systime(/sec)
find_reg_new, m, n,datacv=cvs[1],dataf=sporeout,outfile=rmasks[1]
print, (systime(/sec)-tm)/60., 'Find_reg runtime (minutes)'

mimsy=[[mimsy],[memory()]]

tm=systime(/sec)
find_reg_new, m, n,datacv=cvs[2],dataf=sporeout,outfile=rmasks[2]
print, (systime(/sec)-tm)/60., 'Find_reg runtime (minutes)'

mimsy=[[mimsy],[memory()]]

tm=systime(/sec)
find_reg_new, m, n,datacv=cvs[3],dataf=sporeout,outfile=rmasks[3]
print, (systime(/sec)-tm)/60., 'Find_reg runtime (minutes)'

mimsy=[[mimsy],[memory()]]

tm=systime(/sec)
np=find_pigs_4(spout=spout,cvs=cvs,rmasks=rmasks,dataf=dataf,pig_sp=pig_sp)
print,(systime(/sec)-tm)/60.,'FindPigs_4 Runtime (minutes)'

mimsy=[[mimsy],[memory()]]


stop
end
