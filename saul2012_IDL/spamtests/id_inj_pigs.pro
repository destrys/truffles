pro id_inj_pigs, m, n

fk = '_inj'

delv = 0.736122839/4.

;for n=3, 4 do begin
;    for m=19, 28 do begin
 ;       if (not( (m ge 27) and (n eq 4)) or (fk eq '_inj')) then begin 
            restore, '~goldston/cv/pig_inj_fit_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'
            restore, '~goldston/cv/cv_inj_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'
            restore, '~goldston/cv/rmask_inj_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'
            cv =0.
            extract_coords, '~goldston/Narrow/' + cname(m, n) + '_N.fits', rr, dd, vv,fits
            tlp = fits[*, *, 1024-50/delv:1023+50/delv]
            vels =  vv[ 1024-50/delv:1023+50/delv]
            fits=0
     
            id = fltarr(n_elements(g_inj), 4)
            dra =  fltarr(n_elements(g_inj))
            ddec =  fltarr(n_elements(g_inj))
            dv =  fltarr(n_elements(g_inj))
            da =  fltarr(n_elements(g_inj))
            dpa =  fltarr(n_elements(g_inj))
            
            for i=0, n_elements(g_inj)-1 do begin
; use scaling such that distance is more or less in units of FWHM
                ; found that dv not a good discriminant.
                mn = min( sqrt( ((g_inj[i].ra - pig_fit.ra)^2 + (g_inj[i].dec - pig_fit.dec)^2)*(3./60.)^2.), xx); + 3.^2.*(g_inj[i].v - pig_fit.gfit[1])^2), xx)
                xp = (interpol(findgen(512), ra[*, 0], g_inj[i].ra)  <  511) > 0
                yp = (interpol(findgen(512), dec[0, *], g_inj[i].dec) < 511) > 0
                vp = (interpol(findgen(n_elements(vels)), vels, g_inj[i].v) < (n_elements(vels)-1) ) > 0
                
                id[i, 2] = rmask[xp, yp, vp]
                id[i, 3] = tlp[xp, yp, vp]
                id[i, 0] = mn
                id[i, 1] = xx
                dra[i] = (g_inj[i].ra - pig_fit.ra)[xx]
                ddec[i] =  (g_inj[i].dec - pig_fit.dec)[xx]
                dv[i] = (g_inj[i].v - pig_fit.gfit[1])[xx]
                da[i] = (g_inj[i].a - pig_fit.twodfit[0])[xx]
                dpa[i] = (g_inj[i].pa - pig_fit.twodfit[5])[xx]
   
            endfor
;            stop
            save, id, f='~goldston/cv/pig_id_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'
            
                                ;       endif
;    endfor
;endfor


end        
      
