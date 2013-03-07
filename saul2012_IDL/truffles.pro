pro truffles,datacube=datacube,outstring=outstring,chin=chin,keepedges=keepedges


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;; Set Up Names of Files ;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

datastring=datacube

outstring=outstring

dataf='dataf'+outstring+'.sav'
datacoords='datacoords'+outstring+'.sav'
sporeout='spores'+outstring+'.sav'
spout='pig_sp'+outstring+'.sav'
pigfitout='pig_fit'+outstring+'.sav'
pigplotout='pig_plots'+outstring+'.ps'



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

top = '/hpc/30days/astro/users/drs2125/checkrun/'
datafdir= 'dataf/'
sporedir = 'spore/'
cv = 'cv/'
rmask = 'rmask/'
pigspout='pig_sp/'
pigfit= 'pigfit/'
pigplots= 'pigplots/'

save,dataf,datacoords,outstring,sporeout,spout,cvs,rmasks,wts,pigfitout,pigplotout,f='names'+outstring+'.sav'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; Save Cubes and Coords in Save File ;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

extract_coords,datastring,ra,dec,vels,fits
ra=float(ra)
dec=float(dec)
vels=float(vels)
save,ra,dec,vels,f=datafdir+datacoords

save,fits,f=datafdir+dataf

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;; Spore Does Galactic Subtraction ;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (file_exists(datastring) eq 1 and file_exists(top + sporedir + sporeout) eq 0) then begin
help,/mem
help
tm=systime(/sec)
spore,fits,chin=chin,ra=ra,dec=dec,vels=vels,outfile=sporedir+sporeout
print, (systime(/sec)-tm)/60., 'SPORE runtime (minutes)'

endif

fits=0b
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;; CV_2D Convolves the Cubes ;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (file_exists(top + sporedir + sporeout) eq 1 and file_exists(top + cv + cvs[0]) eq 0) then begin
    tm=systime(/sec)
    cv_2d,chin=chin,fwhm=7,v_fwhm=5,dataf=sporedir+sporeout,$
        outfile=cv+cvs[0],wtout=cv+wts[0],noisescale=0.25
    print, (systime(/sec)-tm)/60., 'CV_2D runtime (minutes)'
endif

if (file_exists(top + sporedir + sporeout) eq 1 and file_exists(top + cv + cvs[1]) eq 0) then begin
    tm=systime(/sec)
    cv_2d,chin=chin,fwhm=18,v_fwhm=5,dataf=sporedir+sporeout,$
        outfile=cv+cvs[1],wtout=cv+wts[1],noisescale=0.2
    print, (systime(/sec)-tm)/60., 'CV_2D runtime (minutes)'
endif

if (file_exists(top + sporedir + sporeout) eq 1 and file_exists(top + cv + cvs[2]) eq 0) then begin
    tm=systime(/sec)
    cv_2d,chin=chin,fwhm=7,v_fwhm=15,dataf=sporedir+sporeout,$
        outfile=cv+cvs[2],wtout=cv+wts[2],noisescale=0.2
    print, (systime(/sec)-tm)/60., 'CV_2D runtime (minutes)'
endif

if (file_exists(top + sporedir + sporeout) eq 1 and file_exists(top + cv + cvs[3]) eq 0) then begin
    tm=systime(/sec)
    cv_2d,chin=chin,fwhm=18,v_fwhm=15,dataf=sporedir+sporeout,$
        outfile=cv+cvs[3],wtout=cv+wts[3],noisescale=1./7.
    print, (systime(/sec)-tm)/60., 'CV_2D runtime (minutes)'
 endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; FIND_REG ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


if (file_exists(top + cv + cvs[0]) eq 1 and file_exists(top + rmask + rmasks[0]) eq 0) then begin
    tm=systime(/sec)
    find_reg, datacv=cv+cvs[0],outfile=rmask+rmasks[0]
    print, (systime(/sec)-tm)/60., 'Find_reg runtime (minutes)'
endif


if (file_exists(top + cv + cvs[1]) eq 1 and file_exists(top + rmask + rmasks[1]) eq 0) then begin
    tm=systime(/sec)
    find_reg, datacv=cv+cvs[1],outfile=rmask+rmasks[1]
    print, (systime(/sec)-tm)/60., 'Find_reg runtime (minutes)'
endif


if (file_exists(top + cv + cvs[2]) eq 1 and file_exists(top + rmask + rmasks[2]) eq 0) then begin
    tm=systime(/sec)
    find_reg, datacv=cv+cvs[2],outfile=rmask+rmasks[2]
    print, (systime(/sec)-tm)/60., 'Find_reg runtime (minutes)'
endif


if (file_exists(top + cv + cvs[3]) eq 1 and file_exists(top + rmask + rmasks[3]) eq 0) then begin
    tm=systime(/sec)
    find_reg, datacv=cv+cvs[3],outfile=rmask+rmasks[3]
    print, (systime(/sec)-tm)/60., 'Find_reg runtime (minutes)'
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;; FIND_PIGS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (file_exists(top + rmask + rmasks[0]) eq 1 and file_exists(top + pigspout + spout) eq 0) then begin
    tm=systime(/sec)
    np=find_pigs(spout=pigspout+spout,cvs=cv+cvs,rmasks=rmask+rmasks,$
        dataf=datafdir+dataf,datacoords=datafdir+datacoords)
    print,(systime(/sec)-tm)/60.,'Find_Pigs Runtime (minutes)'
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;; PIG_PROPS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (file_exists(top + pigspout +spout) eq 1 and file_exists(top + pigfit + pigfitout) eq 0) then begin
    pig_props,outstring,readdir='/hpc/30days/astro/users/drs2125/checkrun/',$
        pigspout=pigspout,rmask=rmask,sporedir=sporedir,datafdir=datafdir,writedir=pigfit
    print,'finished pig fits'
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;if (file_exists(top + pigfit +pigfitout) eq 1 and file_exists(top + pigplots +pigplotout) eq 0) then begin
;plot_pig_prosps_3D,outstring,0,pigplots,'/hpc/30days/astro/users/drs2125/checkrun/',pigspout=pigspout,rmask=rmask,sporedir=sporedir,datafdir=datafdir,pigfit=pigfit
;print,'finished pigplots'
;endif

print,'DONE!'

end
