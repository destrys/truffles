pro plot_pig, i, j, npix=npix, ps=ps, exam=exam, points=points, nonum=nonum, gcheck=gcheck, hideg=hideg, old=old

if keyword_set(gcheck) then begin
    restore, '~/cv/glitch.sav'
endif
restore, '~/cv/pf_all_inj.sav';, /ver
wh = where((idall[*, 2] ne 0) and (gi_all.as lt 2))
whf = where((idall[*, 2] ne 0) and (alog10(idall[*, 0]) lt (-2.5)) and (gi_all.as lt 2) )
pf_all_inj = pf_all[idall[whf, 1]]
;restore, '~goldston/cv/pf.sav';, /ver
;kp = pf.kp
if keyword_set(old) then restore, '~/cv/cv_v2/pf_all.sav' else restore, '~goldston/cv/pf_all.sav';, /ver
if keyword_set(old) then restore, '~/cv/cv_v2/rngs_inj.sav' else restore, '~/cv/rngs_inj.sav'


as =  pf_all.twodfit[3]/pf_all.twodfit[4] > pf_all.twodfit[4]/pf_all.twodfit[3]
whr = where( (pf_all.gfit[1] gt -50) and  (pf_all.gfit[1] lt 50) and  (pf_all.gfit[2]*2*sqrt(2*alog(2)) gt 1) and (pf_all.gfit[2]*2*sqrt(2*alog(2)) lt 5) and (pf_all.twodfit[3]*2*sqrt(2*alog(2)) gt 3) and (pf_all.twodfit[2]*2*sqrt(2*alog(2)) lt 20) and (pf_all.twodfit[0] lt 3) and (as lt 2) )

if keyword_set(gcheck) then begin
    whgl = where(glitch[whr] eq 1, ctg)
    if ctg ne 0 then whg = whr(whgl)
endif

if keyword_set(hideg) then begin
    if ctg ne 0 then whr = whr(where(glitch(whr) eq 0))
    ctg = 0
endif


;twodfit:
; [0] amplitude, [1] x-position, [2] y-position, [3] x-sigma, [4] y-sigma,  [5] position angle, [6] offset, [7] dt/dx, [8] dt/dy

axes = [i, j]
titles = strarr(2)
units = strarr(2)
ranges = fltarr(2,2)

inj_vals = fltarr(2, n_elements(wh))
rvals = fltarr(2, n_elements(pf_all))
fit_inj_vals = fltarr(2, n_elements(whf))

; choose axes
for c = 0, 1 do begin
    case axes[c] of
        ; velocity
        1: begin
            titles[c] = 'V_LSR'
            units[c] = '[km/s]'
            ranges[c, *] = [-50, 50]; v_rng
            inj_vals[c, *] = gi_all[wh].v
            fit_inj_vals[c, *] = gi_all[whf].v
            
            rvals[c, *] = pf_all.gfit[1]
        end
        ; Size
        2: begin
            titles[c] = 'sqrt(a*b)'
            units[c] = '[FWHM, arcmin]'
            ranges[c, *] = [3, 9]; sz_rng
            inj_vals[c, *] = gi_all[wh].sz*sqrt(gi_all[wh].as)
            fit_inj_vals[c, *] = sqrt(pf_all_inj.twodfit[3]*pf_all_inj.twodfit[4])*2*sqrt(2*alog(2)); gi_all[whf].sz
            rvals[c, *] = sqrt(pf_all.twodfit[3]*pf_all.twodfit[4])*2*sqrt(2*alog(2))
        end
        ; velocity width
        3: begin
            titles[c] = 'FWHM'
            units[c] = '[km/s]'
            ranges[c, *] = fwhm_rng
            inj_vals[c, *] = gi_all[wh].fwhm
            fit_inj_vals[c, *] = gi_all[whf].fwhm
            rvals[c, *] = pf_all.gfit[2]*2*sqrt(2*alog(2))
        end
        ; Amplitude
        4: begin
            titles[c] = 'Amplitude'
            units[c] = '[K]'
            ranges[c, *] = [0, 2];a_rng
            inj_vals[c, *] = gi_all[wh].a
            fit_inj_vals[c, *] = pf_all_inj.twodfit[0]
            rvals[c, *] = pf_all.twodfit[0]
        end
       ; RA
        5: begin
            titles[c] = 'RA'
            units[c] = '[deg]'
            ranges[c, *] = [130., 250.]
            inj_vals[c, *] = gi_all[wh].ra
            fit_inj_vals[c, *] = gi_all[whf].ra
            rvals[c, *] = pf_all.ra
        end
        ; Dec
        6: begin
            titles[c] = 'Dec'
            units[c] = '[deg]'
            ranges[c, *] = [-2, 32.3]
            inj_vals[c, *] = gi_all[wh].dec
            fit_inj_vals[c, *] = gi_all[whf].dec
            rvals[c, *] = pf_all.dec
        end
        ; Aspect ratio
        7: begin
            titles[c] = 'Aspect Ratio'
            units[c] = ''
            ranges[c, *] = [1, 2]
            inj_vals[c, *] = gi_all[wh].as
            fit_inj_vals[c, *] = gi_all[whf].as
            rvals[c, *] = pf_all.twodfit[3]/pf_all.twodfit[4] > pf_all.twodfit[4]/pf_all.twodfit[3]
        end
        ; l
        8: begin
            titles[c] = 'Galactic l'
            units[c] = '[deg]'
            ranges[c, *] = [0, 360]
            glactc,  gi_all[wh].ra/15., gi_all[wh].dec, 2000, l, b, 1
            inj_vals[c, *] = l
            glactc,  gi_all[whf].ra/15., gi_all[whf].dec, 2000, lf, bf, 1
            fit_inj_vals[c, *] = lf
            glactc, pf_all.ra/15., pf_all.dec, 2000, lpf, bpf, 1
            rvals[c, *] = lpf
        end
        ; b
        9: begin
            titles[c] = 'Galactic b'
            units[c] = '[deg]'
            ranges[c, *] = [40., 90.]
            glactc,  gi_all[wh].ra/15., gi_all[wh].dec, 2000, l, b, 1
            inj_vals[c, *] = b
            glactc,  gi_all[whf].ra/15., gi_all[whf].dec, 2000, lf, bf, 1
            fit_inj_vals[c, *] = bf
            glactc, pf_all.ra/15., pf_all.dec, 2000, lpf, bpf, 1
            rvals[c, *] = bpf
        end

        ; position angle
        10: begin
            titles[c] = 'Position Angle'
            units[c] = '[deg]'
            ranges[c, *] = [0, 180]
            inj_vals[c, *] = gi_all[wh].pa;  *(gi_all[wh].as gt 1.5)
            fit_inj_vals[c, *] = gi_all[whf].pa;  *(gi_all[whf].pa gt 1.5)
            rvals[c, *] = ((pf_all.twodfit[5] + (pf_all.twodfit[3] gt pf_all.twodfit[4])*90. + 180) mod 180);  *( (pf_all.twodfit[3]/pf_all.twodfit[4] > pf_all.twodfit[4]/pf_all.twodfit[3]) gt 1.5)
        end


    endcase
endfor

if not(keyword_set(npix)) then npix = [4, 6]

xrng = v_rng
yrng = fwhm_rng

found =  hist_2d(  fit_inj_vals[0, *],fit_inj_vals[1, *], bin1=(ranges[0, 1]-ranges[0, 0])/npix[0], bin2= (ranges[1,1]-ranges[1,0])/npix[1], min1=ranges[0,0], min2=ranges[1, 0],max1=ranges[0, 1]-1d-6, max2=ranges[1,1]-1d-6)

all =  hist_2d(  inj_vals[0, *],inj_vals[1, *], bin1=(ranges[0, 1]-ranges[0, 0])/npix[0], bin2= (ranges[1,1]-ranges[1,0])/npix[1], min1=ranges[0,0], min2=ranges[1, 0], max1=ranges[0, 1]-1d-6, max2=ranges[1,1]-1d-6)

if keyword_set(ps) then psopen, titles[0]+ '_vs_' + titles[1] + '.eps', /encapsulated
loadct, 0, /sil
disp, float(found)/float(all), (findgen(npix[0]))*(ranges[0,1]-ranges[0,0])/npix[0] + ranges[0,0],  (findgen(npix[1]))*(ranges[1,1]-ranges[1, 0])/npix[1] + ranges[1,0], position=[0.1, 0.1, 0.8, 0.9], xtitle=titles[0] + ' ' + units[0], ytitle=titles[1] + ' ' + units[1], charsize=2

if keyword_set(points) then begin
    loadct, 13, /sil
    circle, /fill
    oplot,  inj_vals[0, *],inj_vals[1, *], color=50, psym=8, symsize=0.5
    oplot,  inj_vals[0, where(alog10(idall[wh, 0]) lt (-2.5))],inj_vals[1, where(alog10(idall[wh, 0]) lt (-2.5))], color=120, psym=8, symsize=0.5
endif
    
oplot, rvals[0, whr], rvals[1, whr], psym=1, thick=3, symsize=2, color=255
oplot, rvals[0, whr], rvals[1, whr], psym=1, thick=2, symsize=1, color=0

loadct, 13, /sil

if i eq 5 and j eq 6 then begin
    cover = readfits('/share/galfa/html/DR1-S/allT.fits', hdr)
    extast, hdr, astr
    xy2ad, findgen(21632), fltarr(21632), astr, xax, y0
    xy2ad, fltarr(2432), findgen(2432), astr, x0, yax
    mc = cover*0
    
    restore, '~goldston/cv/docubes.sav'
    restore, '/share/galfa/html/DR1-S/tilenames.sav'
    for q=0, n_elements(xw)-1 do begin
        mc(where((tilex eq xw[q]) and (tiley eq yw[q]))) = 1.
    endfor
    cvr = cover*mc gt 0.
    contour, cvr, xax, yax,/overplot
endif

    
if not (keyword_set(nonum)) then begin 
    for i=0, n_elements(whr)-1 do xyouts,  rvals[0, whr[i]], rvals[1, whr[i]], string(i, f='(I2.2)'),color=100, charsize=2, charthick=2
endif

plots, rvals[0, whr], rvals[1, whr], psym=1, thick=2, symsize=1;, color=kp[whr]*255

if keyword_set(gcheck) then  begin
    if ctg ne 0 then oplot, rvals[0, whg], rvals[1, whg], psym=1, thick=2, symsize=1, color=100
endif



loadct, 0, /sil

pbang = !p
xbang = !x
ybang = !y


fr = minmax(found/all)
;stop
disp, rebin(reform(findgen(100), 1, 100), 2, 100), findgen(2), (findgen(100)/99.)*(fr[1]-fr[0])+fr[0], /noerase, position=[0.85, 0.1, 0.9, 0.9], xticklen=1d-6, xtickname=replicate(' ' , 10)


if (keyword_set(exam) and (not (keyword_set(ps)))) then begin
    window, 1, xsi=200, ysi=200, xp=0, yp=1700
    wset, 0
    ch = 0
    !y = ybang
    !x = xbang
    !p = pbang
    

    !mouse.button = 1
    while !mouse.button ne 2 do begin
        cursor, x, y, /change
   
        mn = min( sqrt((x-rvals[0, whr])^2 + (y-rvals[1, whr])^2), ch1)
        if keyword_set(gcheck) then begin
            if !mouse.button eq 4 then begin
                glitch[whr[ch1]] = (glitch[whr[ch1]] + 1) mod 2
                loadct, 13, /sil
                plots, rvals[0, whr[ch1]], rvals[1, whr[ch1]], psym=1, thick=2, symsize=1, color=255- 155*glitch[whr[ch1]]
                loadct, 0, /sil
                !mouse.button = 1.
                wait, 1
            endif
        endif

     if ch1 ne ch then begin
            ch = ch1
            pbang = !p
            xbang = !x
            ybang = !y
            wset, 1
            display, pf_all[whr[ch]].postagestamp,pf_all[whr[ch]].ra - (findgen(31)-15)/60.,pf_all[whr[ch]].dec + (findgen(31)-15)/60., title=ch, xtitle='RA [deg]', ytitle='Dec [deg]'
            wset, 0
            !y = ybang
            !x = xbang
            !p = pbang
        endif
    endwhile
endif


if keyword_set(ps) then psclose

if keyword_set(gcheck) then begin
    save, glitch, f='~goldston/cv/glitch.sav'
endif


end
