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

save,dataf,datacoords,outstring,sporeout,spout,cvs,rmasks,wts,f='names'+outstring+'.sav'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; Save Cubes and Coords in Save File ;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;extract_coords,datastring,ra,dec,vels,fits
;ra=float(ra)
;dec=float(dec)
;vels=float(vels)
;save,ra,dec,vels,f=datacoords


;save,fits,f=dataf



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;; Spore Does Galactic Subtraction ;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

tm=systime(/sec)
;spore,fits,chin=chin,ra=ra,dec=dec,vels=vels
print, (systime(/sec)-tm)/60., 'SPORE runtime (minutes)'

;;;;;;;;;;;;;;;;;;;;;;;;; Slice Off Edges ;;;;;;;;;;;;;;;;;;;;;;;;;;;
;sz=size(fits)
;if ~keyword_set(keepedges) then begin
;fits=temporary(fits[30:sz[1]-30-1,30:sz[2]-30-1,*])
;endif


;save,fits,f=sporeout
;fits=0b



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;; CV_ALL_2D Convolves the Cubes ;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


tm=systime(/sec)
cv_all_2d,chin=chin,fwhm=7,v_fwhm=5,dataf=sporeout,outfile=cvs[0],wtout=wts[0]
print, (systime(/sec)-tm)/60., 'CV_ALL runtime (minutes)'


tm=systime(/sec)
cv_all_2d,chin=chin,fwhm=18,v_fwhm=5,dataf=sporeout,outfile=cvs[1],wtout=wts[1]
print, (systime(/sec)-tm)/60., 'CV_ALL runtime (minutes)'


tm=systime(/sec)
cv_all_2d,chin=chin,fwhm=7,v_fwhm=15,dataf=sporeout,outfile=cvs[2],wtout=wts[2]
print, (systime(/sec)-tm)/60., 'CV_ALL runtime (minutes)'


tm=systime(/sec)
cv_all_2d,chin=chin,fwhm=18,v_fwhm=15,dataf=sporeout,outfile=cvs[3],wtout=wts[3]
print, (systime(/sec)-tm)/60., 'CV_ALL runtime (minutes)'


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

restore,cvs[0]
cv/=4
save,cv,f=cvs[0]

restore,cvs[1]
cv/=5
save,cv,f=cvs[1]

restore,cvs[2]
cv/=5
save,cv,f=cvs[2]

restore,cvs[3]
cv/=7
save,cv,f=cvs[3]


end


tm=systime(/sec)
find_reg_newest, datacv=cvs[0],outfile=rmasks[0]
print, (systime(/sec)-tm)/60., 'Find_reg runtime (minutes)'


tm=systime(/sec)
find_reg_newest, datacv=cvs[1],outfile=rmasks[1]
print, (systime(/sec)-tm)/60., 'Find_reg runtime (minutes)'


tm=systime(/sec)
find_reg_newest, datacv=cvs[2],outfile=rmasks[2]
print, (systime(/sec)-tm)/60., 'Find_reg runtime (minutes)'


tm=systime(/sec)
find_reg_newest, datacv=cvs[3],outfile=rmasks[3]
print, (systime(/sec)-tm)/60., 'Find_reg runtime (minutes)'



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;; Slice off Edges of fits/ra/dec/vels ;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

restore,dataf
restore,datacoords
sz=size(fits)
if ~keyword_set(keepedges) then begin
fits=temporary(fits[30:sz[1]-30-1,30:sz[2]-30-1,*])
ra=ra[30:sz[1]-30-1,30:sz[2]-30-1]
dec=dec[30:sz[1]-30-1,30:sz[2]-30-1]
endif
sz=size(fits)
nx=41
fits=temporary(fits[nx/2:sz[1]-nx/2-1,nx/2:sz[2]-nx/2-1,*])
ra=ra[nx/2:sz[1]-nx/2-1,nx/2:sz[2]-nx/2-1]
dec=dec[nx/2:sz[1]-nx/2-1,nx/2:sz[2]-nx/2-1]


save,fits,f=dataf
fits=0b
save,ra,dec,vels,f=datacoords


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

tm=systime(/sec)
np=find_pigs_4(spout=spout,cvs=cvs,rmasks=rmasks,dataf=dataf,pig_sp=pig_sp,datacoords=datacoords)
print,(systime(/sec)-tm)/60.,'FindPigs_4 Runtime (minutes)'




end
