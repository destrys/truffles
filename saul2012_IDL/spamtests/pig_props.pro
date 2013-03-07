pro pig_props, m, n, fake=fake,dataf=dataf

; step 4 in the pig-searching code. This determines the properties
; of the pigs, using the mpfit_2dgauss code and mpfit.
; again, fk = '_inj' uses the injected pigs code.

dp = 1
lp = 2
voff=100
delv = 0.736122839/4.
parinfo = replicate({fixed:0., limited:fltarr(2), limits:fltarr(2)}, 9)
parinfo[1].limited=[1, 1]
parinfo[2].limited=[1, 1]
parinfo[1].limits=[-2, 2]
parinfo[2].limits=[-2, 2]
if keyword_set(fake) then fk = '_inj' else fk = ''


;for n=3, 4 do begin
;    for m=19, 28 do begin
;        if (not( (m ge 27) and (n eq 4)) or (fk eq '_inj')) then begin 
        print, m, n

        restore, '~/cv/pig'+fk+'_sp_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'
        restore, '~/cv/cv'+fk+'_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'
        szcv = size(cv)
        szv = szcv[3]
if keyword_set(dataf) then begin
extract_coords,dataf,ra,dec,vels,fits
endif else begin 
       extract_coords, cname(m, n) + '_W.fits', ra, dec, vels,fits
endelse
        pssize = 15
       pig_fit = replicate({gfit:fltarr(6), sigs:fltarr(6), twodfit:fltarr(9), ra:0., dec:0., postagestamp:fltarr(pssize*2+1, pssize*2+1)}, n_elements(pig_sp))

       ; tlp = fits[*, *, 1024-voff/delv:1023+voff/delv]
tlp=fits
        fits = 0
        ;vels = vels[ 1024-voff/delv:1023+voff/delv]
        
        ; re inject the pigs
        if fk eq '_inj' then begin
            injector, tlp, vels, ra, dec, n_elements(g_inj), v_rng, ra_rng, dec_rng, sz_rng, as_rng, a_rng, pa_rng, fwhm_rng, g_inj, useginj=1
        endif

        for j=0, n_elements(pig_sp)-1 do begin
            ; finding a range over which to fit the on-off spectrum
            mn = (min(where(pig_sp[j].spm eq 1))-10./delv) > 0
            mx = (max(where(pig_sp[j].spm eq 1))+10./delv) < (szv-1)
            wvf = findgen(mx-mn+1)+mn

            gf = gaussfit( vels[wvf], (pig_sp[j].sp- pig_sp[j].osp)[wvf], fit, estimates=[1, mean(vels[wvf]), total(pig_sp[j].spm)*delv/2., 0, 0,0], sig=sig)
            pig_fit[j].gfit = fit
            pig_fit[j].sigs = sig
            
            img_on = total(tlp[*, *, where(pig_sp[j].spm eq 1)], 3)/total(pig_sp[j].spm) ; in K km/s
            d_off = 5
            off_plus = where(pig_sp[j].spm eq 1) + total(pig_sp[j].spm) + d_off
            off_minus = where(pig_sp[j].spm eq 1) - total(pig_sp[j].spm) - d_off
            nw = n_elements(where(pig_sp[j].spm eq 1))
            whbdmn = where(off_minus lt 0, ctm)
            if (ctm ne 0) and (ctm lt nw) then off_minus = off_minus[where(off_minus ge 0)]
            whbdpl = where(off_plus gt (szv-1), ctp)
            if (ctp ne 0) and (ctp lt nw) then off_plus = off_plus[where(off_plus le (szv-1))]
            if (ctp eq nw) then img_off =  total(tlp[*, *, off_minus], 3)/n_elements(off_minus)
            if (ctm eq nw) then img_off =  total(tlp[*, *, off_plus], 3)/n_elements(off_plus)
            if (ctm ne nw) and (ctp ne nw) then img_off =  total(tlp[*, *, [off_plus,off_minus]], 3)/n_elements([off_plus,off_minus])

            img = img_on-img_off

            xx = interpol(findgen(512), ra[*, 0],pig_sp[j].ra)
            yy = interpol(findgen(512), dec[0, *],pig_sp[j].dec)
            

            img_bord = fltarr(512+pssize*2, 512+pssize*2)
            img_bord[pssize:pssize+511, pssize:pssize+511] = img
            psimg = img_bord[xx :xx+pssize*2, yy:yy+pssize*2] 
            
            sir = 7
            subimg = img[(xx-sir) > 0 :(xx+sir) < 511, (yy-sir) > 0 :(yy+sir) < 511] 
            subra = ra[(xx-sir) > 0 :(xx+sir) < 511, (yy-sir) > 0 :(yy+sir) < 511] 
            subdec = dec[(xx-sir) > 0 :(xx+sir) < 511, (yy-sir) > 0 :(yy+sir) < 511] 
            szsub = size(subimg)


            fa = {x:subimg, err:(subimg*0.+1)}
            
            params = mpfit('mpfit_2dgauss', [1., 0., 0., 4., 4., 0., 0., 0., 0.], parinfo=parinfo, functargs=fa)
            
         
            pig_fit[j].twodfit = params
            pig_fit[j].ra = subra[params[1]+(szsub[1]-1.)/2., 0]
            pig_fit[j].dec = subdec[0, params[2]+(szsub[2]-1.)/2.]
            pig_fit[j].postagestamp = psimg
        endfor
        
        save, pig_fit, f='~/cv/pig'+fk+'_fit_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'

;    endfor
;endfor


end
