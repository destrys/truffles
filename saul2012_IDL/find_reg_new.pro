pro find_reg_new, dataf=dataf,fits=fits,rmask=rmask,minpeak=minpeak,water=water,datacv=datacv,outfile=outfile

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


; get the original data set
if ~keyword_set(fits) then begin
restore,dataf,/verb
fits = float(fits)
endif

szcv = size(fits)
dv = szcv[3]
dx=szcv[1]
dy=szcv[2]
; 2d mask
mask2d = bytarr(dx, dy)

; where there was no data
wh = where((fits[*,*,dv/2] lt (-200)) or (finite(fits[*,*,dv/2]) EQ 0) , ct)
if ct ne 0 then mask2d[wh] =1

fits=0b


; the new, convolved, unsharpmasked data set
if ~keyword_set(datacv) then begin
	print,'Set DATACV so I have some data to work with'
endif else begin
	restore,datacv,/verb
endelse


; grown to cover the areas where the convolution contaminated the map
bigmask = gromask(gromask(gromask(gromask(gromask(gromask(gromask(gromask(gromask(gromask(mask2d))))))))))
bigmask = 1b-bigmask
; 3d
for i=0,szcv[3]-1 do begin
cv[*,*,i]=cv[*,*,i]*bigmask
endfor

;also get rid of velocity edges
cv[*, *, 0:15] = 0
cv[*, *, dv-16:*] = 0

; make a copy to house the region data
rmask = byte(cv*0)+1b

; preset values for where the coverage can spread to, and which peaks
; to select

if ~keyword_set(minpeak) then minpeak = 5.
if ~keyword_set(water) then cutoff=3. else cutoff=water
mx = 100b
q = 2b

help,/mem
while mx gt minpeak do begin
    print, q
    pk = max(cv, xx,/nan)
    qmask = byte(rmask*0b)
    rmask[xx] = q
    qmask[xx] = 1b
    grow = 1
    tnm = 0
    st = 0

    wq = xx
tm=systime(/sec)/60.
    while grow eq 1 do begin
	vinds=floor(wq/(dx*dy))
	yinds=(wq-(vinds*dx*dy))/dx
	xinds=(wq-(vinds*dx*dy)) mod dx
	mmx=minmax(xinds)
	mmy=minmax(yinds)
	mmv=minmax(vinds)
        
	; subregion containing small area
        sr_cv =       cv[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)]
        sr_rmask = rmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)]                   
        sr_qmask = qmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)]
        newmask1 = shift(sr_qmask, 1, 0, 0) + sr_qmask + shift(sr_qmask, -1, 0, 0) 
        newmask2 = shift(newmask1, 0, 1, 0) + newmask1 + shift(newmask1, 0, -1, 0) 
        newmask3 = shift(newmask2, 0, 0, 1) + newmask2 + shift(newmask2, 0, 0, -1) 
        newmask3 = (newmask3 < 1)
        wh = where((sr_cv lt cutoff) or (finite(sr_cv) eq 0), ct)
        if ct ne 0 then newmask3[wh] = 0
 
        wh = where((sr_rmask ne q) and (sr_rmask ne 1), ct)
        if ct ne 0 then newmask3[wh] = 0    
        rmn = rmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)]
	wh2=where(newmask3 EQ 1)
        rmnwh = rmn[wh2]
        rmn[wh2] = q 

     
        if max(rmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)] gt rmn) eq 1 then stop
        rmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)] = rmn

        qmn = qmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)]
        qmn[where(newmask3 eq 1)] = 1
  
          
        qmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)] = qmn
        wqmn = where(qmn eq 1)
        convert_wh, wqmn, (size(qmn))[1], (size(qmn))[2], (size(qmn))[3], (mmx[0]-1) > 0, (mmy[0]-1) > 0,(mmv[0]-1) > 0, szcv[1], szcv[2], szcv[3], wq

if (total(newmask3) - tnm) eq 0 then grow = 0

        tnm = total(newmask3)

     endwhile
help,/mem    
cv[where(rmask eq q)] = 0
    mx = max(cv,/nan)
    q++
    if q EQ 255 then begin
	rmask=fix(rmask)
	q=fix(q)
    endif

endwhile

if keyword_set(outfile) then begin
save,rmask,f=outfile
endif


end
