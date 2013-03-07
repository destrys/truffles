
; ORIGINAL CODE DESCRIPTION:
;The point of this code is to take a 3D unsharp masked data cube of
;the 'PIGS' region and find possible pig candidates. The plan is to
;mask the data, then find a bright peak, grow a mask for the bright
;peak until it reaches, say, 4x the noise in the map, and find out how
;large that region is. Then repeat with the next brightest position,
;outside that region, repeating until all regions above some amplitude
;are covered
; settibng fk to '_inj' uses the pig injected cubes from cv_all

; adapted from find_reg (above) to do a watertable analysis.

dx = 512
dy = 512
dv = 543

for f=3, 3 do begin
    for p=20, 28 do begin
; get the original data set
tlp = readfits('~goldston/Narrow/' + cname(p, f) + '_N.fits', nslice=1024)
;tlpall = fits[1500:2500, *, *]

; the new, convolved, unsharpmasked data set

; 2d mask
mask2d = fltarr(dx, dy) +1

; where there was no data
mask2d[where(tlp lt (-200))] = 0.

; 3d
allmask = rebin(reform(mask2d, dx, dy, 1), dx, dy, dv)

; coordinates
xxx = rebin(reform(findgen(dx), dx, 1, 1), dx, dy, dv)
yyy = rebin(reform(findgen(dy), 1, dy, 1), dx, dy, dv)
vvv = rebin(reform(findgen(dv), 1, 1, dv), dx, dy, dv)

mx = 100
q = 2
levels = findgen(10)
cutoff= 0.5

for i=0, n_levels do begin
    
    rmask = allmask
    qmask = rmask
    cv = allmask
    cv(where(tlp lt levels[i])) = 0.

    while mx ge 0 do begin
        pk = max(rmask, xx)
        rmask[xx] = q
        qmask[xx] = 1
        grow = 1
        tnm = 0.
        st = 0
    while grow eq 1 do begin
        wq = where(rmask eq q)
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

        wh = where(sr_cv lt cutoff, ct)
        if ct ne 0 then newmask3(where(sr_cv lt cutoff)) = 0.
 
        wh = where((sr_rmask ne q) and (sr_rmask ne 1), ct)
        if ct ne 0 then newmask3[wh] = 0.    
 
        rmn = rmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)]
        rmnwh = rmn[where(newmask3 eq 1)]
 ;       print, rmnwh(uniq(rmnwh, sort(rmnwh)))
        
        rmn[where(newmask3 eq 1)] = q 
            
        if max(rmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)] gt rmn) eq 1 then stop
        rmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)] = rmn
        
        qmn = qmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)]
        qmn[where(newmask3 eq 1)] = 1.

          
        qmask[(mmx[0]-1) > 0: (mmx[1]+1) < (dx-1), (mmy[0]-1) > 0: (mmy[1]+1) < (dy-1),(mmv[0]-1) > 0: (mmv[1]+1) < (dv-1)] = qmn
        if (total(newmask3) - tnm) eq 0 then grow = 0

        tnm = total(newmask3)
        

    endwhile
    lmask[where(rmask eq q)] = 0.
    mx = max(cv*lmask)
    q++
    print, q
;    fasthist, rmask(where(rmask gt 1))
endwhile

save, rmask, p, f, f='~goldston/cv/rmask'+fk+'_' + string(p, f='(I2.2)') +'_' + string(f, f='(I2.2)')+ '.sav'

endfor
endfor

end
