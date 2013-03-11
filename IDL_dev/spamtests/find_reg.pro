pro find_reg, p, f, fake=fake,dataf=dataf,fits=fits,cv=cv,rmask=rmask,minpeak=minpeak,water=water,datacv=datacv,outfile=outfile

;This is step 2 in the pig-searching code suite
;
;The point of this code is to take a 3D unsharp masked data cube of
;the 'PIGS' region and find possible pig candidates. The plan is to
;mask the data, then find a bright peak, grow a mask for the bright
;peak until it reaches, say, 4x the noise in the map, and find out how
;large that region is. Then repeat with the next brightest position,
;outside that region, repeating until all regions above some amplitude
;are covered
; settibng fk to '_inj' uses the pig injected cubes from cv_all


;Define Growth Measurement Array (Temp by Destry Jan. 23)
followgrowth=fltarr(1000,1000)


; for fake data set to '_inj', else = ''
if keyword_set(fake) then fk = '_inj' else fk = ''

; get the original data set
if keyword_set(fits) then begin
tlp=fits
endif else begin
tlp = readfits(dataf)
endelse

; the new, convolved, unsharpmasked data set
if ~keyword_set(cv) then begin
	if ~keyword_set(datacv) then begin
		restore, '~destrys/cv/cv'+fk+'_' + string(p, f='(I2.2)') +'_' + string(f, f='(I2.2)')+ '.sav', /ver
	endif else begin
		restore,datacv,/verb
	endelse
endif

szcv = size(cv)
dv = szcv[3]
dx=szcv[1]
dy=szcv[2]
; 2d mask
mask2d = fltarr(dx, dy)

; where there was no data
wh = where((tlp[*,*,dv/2] lt (-200)) or (finite(tlp[*,*,dv/2]) EQ 0) , ct)
if ct ne 0 then mask2d[wh] =.1
tlp = 0.
; grown to cover the areas where the convolution contaminated the map
bigmask = gromask(gromask(gromask(gromask(gromask(gromask(gromask(gromask(gromask(gromask(mask2d))))))))))

; 3d
allmask = rebin(reform(bigmask, dx, dy, 1), dx, dy, dv)

; now, 1 means is good, 0 means ignore
allmask = 1- allmask

;also get rid of velocity edges
allmask[*, *, 0:15] = 0.
allmask[*, *, dv-16:*] = 0.

; make a copy to house the region data
rmask = allmask

; a copy to mask out in loop
lmask = allmask

; coordinates
xxx = rebin(reform(findgen(dx), dx, 1, 1), dx, dy, dv)
yyy = rebin(reform(findgen(dy), 1, dy, 1), dx, dy, dv)
vvv = rebin(reform(findgen(dv), 1, 1, dv), dx, dy, dv)


; preset values for where the coverage can spread to, and which peaks
; to select
help,/mem
if ~keyword_set(minpeak) then minpeak = 5.
if ~keyword_set(water) then cutoff=3. else cutoff=water
mx = 100
q = 2
while mx gt minpeak do begin
    print, q
    pk = max(cv*lmask, xx,/nan)
    qmask = allmask*0.
    rmask[xx] = q
    qmask[xx] = 1.
    grow = 1
    tnm = 0.
    st = 0
    itter = 0.  ;counter for followgrowth

    wq = xx
    while grow eq 1 do begin
        ;wq = where(rmask eq q)
        tsub = systime(/sec)
        mmx = minmax(xxx[wq])
        mmy = minmax(yyy[wq])
        mmv = minmax(vvv[wq])
        ; subregion containing small area
        sr_cv =       cv[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)]
        sr_rmask = rmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)]                   
        sr_qmask = qmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)]
        newmask1 = shift(sr_qmask, 1, 0, 0) + sr_qmask + shift(sr_qmask, -1, 0, 0) 
        newmask2 = shift(newmask1, 0, 1, 0) + newmask1 + shift(newmask1, 0, -1, 0) 
        newmask3 = shift(newmask2, 0, 0, 1) + newmask2 + shift(newmask2, 0, 0, -1) 
        newmask3 = (newmask3 < 1)
        wh = where((sr_cv lt cutoff) or (finite(sr_cv) eq 0), ct)
        if ct ne 0 then newmask3(where((sr_cv lt cutoff) or (finite(sr_cv) Eq 0))) = 0.
 
        wh = where((sr_rmask ne q) and (sr_rmask ne 1), ct)
        if ct ne 0 then newmask3[wh] = 0.    
        rmn = rmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)]
        rmnwh = rmn[where(newmask3 eq 1)]
        
        rmn[where(newmask3 eq 1)] = q 

     
        if max(rmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)] gt rmn) eq 1 then stop
        rmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)] = rmn

        qmn = qmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)]
        qmn[where(newmask3 eq 1)] = 1.
  
          
        qmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)] = qmn
        wqmn = where(qmn eq 1)
        convert_wh, wqmn, (size(qmn))[1], (size(qmn))[2], (size(qmn))[3], (mmx[0]-1) > 0, (mmy[0]-1) > 0,(mmv[0]-1) > 0, szcv[1], szcv[2], szcv[3], wq

if (total(newmask3) - tnm) eq 0 then grow = 0
        followgrowth[itter,q]=total(newmask3)-tnm ;FOLOWGROWTH populate

        tnm = total(newmask3)
     endwhile
help,/mem
    lmask[where(rmask eq q)] = 0.
    mx = max(cv*lmask,/nan)
    q++
endwhile

if ~keyword_set(outfile) then begin
save, rmask, p, f, f='~/cv/rmask'+fk+'_' + string(p, f='(I2.2)') +'_' + string(f, f='(I2.2)')+ '.sav'
endif else begin
save,rmask,p,f,f=outfile
endelse


end
