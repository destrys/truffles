; Go look for all the rmask regions that are small in ra/dec/vel space
; as save them into a structure called pig_sp
function find_pigs, m, n, fake=fake,dataf=dataf,ra=ra,dec=dec,fits=fits,vels=vels,cv=cv,pig_sp=pig_sp,rmask=rmask
;+
;NAME: 
;     FIND_PIGS
;PURPOSE:
;     Look through the rmask output of FIND_REG and select the 
;     regions that are 'small'. 
;     These Regions are saved into a structure called pig_sp
;INPUT:
;     m: Cube Index #
;     n: Cube Index #
;KEYWORDS:
;     /FAKE: Indicates that there are INnjected Guinea Pigs in the data.
;     DATAF: Filename of Data Cube if different from CNAME format.
;     RA   : RA Array
;     Dec  : Dec Array
;     Fits : Data Cube
;     Vels : Velocity Array
;     Rmask: Masked Region Cube
;     pig_sp: Clouds structure.
;HISTORY:
;     thepast:  Written by JEGP
;     Jan. 12, 2010: Documented (Destry)
;                  : Changed 'small' constraints
;     Jan. 14, 2010: Inserted DATAF keyword. (Destry)
;     Jan. 22, 2010: Added Output Keywords pig_sp
;                    Added Input Keywords ra,dec,fits,vels,rmask (Destry)  
;-
if keyword_set(fake) then fk = '_inj' else fk=''


dx = 40
dy = 40
dv = 200
voff = 100
; now find all pig candidates:

if ~keyword_set(cv) then begin
restore, '~/cv/rmask'+fk+'_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'
endif

sz = size(rmask)
xxx = rebin(reform(findgen(sz[1]), sz[1], 1, 1), sz[1], sz[2], sz[3])
yyy = rebin(reform(findgen(sz[2]), 1, sz[2], 1), sz[1], sz[2], sz[3])
vvv = rebin(reform(findgen(sz[3]), 1, 1, sz[3]), sz[1], sz[2], sz[3])

mx = max(rmask)
amiapig = fltarr(mx+1)
pigmask = rmask*0.

for j=2,mx do begin
    if j/10 eq j/10. then print, '[' + string(m, f='(I2.2)') + ',' + string(n, f='(I2.2)') + '], ' + string(j) + '/' +string(mx)
    mmx = minmax(xxx[where(rmask eq j)])
    mmy = minmax(yyy[where(rmask eq j)])
    mmv = minmax(vvv[where(rmask eq j)])

    szx = max(xxx[where(rmask eq j)]) - min(xxx[where(rmask eq j)])
    szy = max(yyy[where(rmask eq j)]) - min(yyy[where(rmask eq j)])
    szv = max(vvv[where(rmask eq j)]) - min(vvv[where(rmask eq j)])

    ; are we next to an edge?
    rsm = rmask[(mmx[0]-1) > 0: (mmx[1]+1) < (sz[1]-1), (mmy[0]-1) > 0: (mmy[1]+1) < (sz[2]-1),(mmv[0]-1) > 0: (mmv[1]+1) < (sz[3]-1)]

    newmask = rsm*0
    newmask(where(rsm eq j)) = 1.
    newmask1 = shift(newmask, 1, 0, 0) + newmask + shift(newmask, -1, 0, 0) 
    newmask2 = shift(newmask1, 0, 1, 0) + newmask1 + shift(newmask1, 0, -1, 0) 
    newmask3 = shift(newmask2, 0, 0, 1) + newmask2 + shift(newmask2, 0, 0, -1) 
    newmask3 = (newmask3 < 1)
    whh = where((newmask3 eq 1) and (rsm eq 0), cct)
    if ((szx le dx) and (szy le dy) and (szv le dv) and (cct eq 0)) then amiapig[j] = 1.
    pigmask([where(rmask eq j)]) = amiapig[j]
endfor

display, total(pigmask, 3)

;plot 'em
delv = 0.736122839/4.

;;; opening cube again.... ;;;;;
if ~keyword_set(fits) then begin
	if keyword_set(dataf) then begin
		extract_coords,dataf,ra,dec,vels,fits
	endif else begin
		extract_coords, cname(m, n) + '_W.fits', ra, dec, vels,fits
	endelse
endif
tlp=fits
if fk eq '_inj' then begin
    restore, '~/cv/cv'+fk+'_' + string(p, f='(I2.2)') +'_' + string(f, f='(I2.2)')+ '.sav', /ver
    cv =0

    injector, tlp, vels, ra, dec, n_elements(g_inj), v_rng, ra_rng, dec_rng, sz_rng, as_rng, a_rng, pa_rng, fwhm_rng, g_inj, useginj=1
endif
    
if total(amiapig) eq 0 then return, total(amiapig)

pig_sp = replicate({ra:0., dec:0., sp:fltarr(sz[3]), vels:fltarr(sz[3]), rm_num:0., spm:fltarr(sz[3]), osp:fltarr(sz[3])}, total(amiapig))
pgn = 0

smbound = 10.


for j=2, mx do begin
    if amiapig[j] eq 1 then begin
        xc = mean(xxx(where(rmask eq j)))
        yc = mean(yyy(where(rmask eq j)))
        xmx = (xc+smbound) < (sz[1] -1)
        xmn = (xc-smbound) > 0
        ymx = (yc+smbound) < (sz[2] -1)
        ymn = (yc-smbound) > 0
        msk = fltarr(xmx-xmn+1, ymx-ymn+1, sz[3])
        szm = size(msk)
        msk(where(rmask[xmn:xmx,ymn:ymx, *] eq j)) = 1.
        imgmask = total(msk, 3) < 1
        holemask = gromask(gromask(gromask(imgmask)))
        offmask = gromask(gromask(gromask(holemask))) - holemask
        wh = where((rmask[xmn:xmx, ymn:ymx, sz[3]/2]) eq 0, ct)
        if ct ne 0 then offmask[wh] = 0.
        spmask = total(total(msk, 1), 1) < 1
        spect = total(total(rebin(reform(imgmask, szm[1], szm[2], 1), szm[1], szm[2], szm[3])*tlp[xmn:xmx,ymn:ymx, *], 1), 1)/total(imgmask)
        offspect = total(total(rebin(reform(offmask, szm[1], szm[2], 1), szm[1], szm[2], szm[3])*tlp[xmn:xmx,ymn:ymx, *], 1), 1)/total(offmask)
        plot, vels, spect, title=string(m) + ', ' + string(n) + ', ' + string(j)
        oplot, vels(where(spmask eq 1)), spect(where(spmask eq 1)), color=100
        pig_sp[pgn].ra = mean((ra[xmn:xmx,ymn:ymx])(where(imgmask eq 1)))
        pig_sp[pgn].dec = mean((dec[xmn:xmx,ymn:ymx])(where(imgmask eq 1)))
        pig_sp[pgn].sp = spect
        pig_sp[pgn].vels = vels
        pig_sp[pgn].spm = spmask
        pig_sp[pgn].osp = offspect
        pig_sp[pgn].rm_num = j
        pgn++
    endif
endfor

save, pig_sp, f='~/cv/pig'+fk+'_sp_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'

return, total(amiapig)


end
