pro plot_pig_prosps, m, n
delv = 0.736122839/4.

dp = 1
lp = 2

theta = findgen(50)/49.*2*!pi
;for n=3, 4 do begin
;    for m=19, 28 do begin
;        if not( (m ge 27) and (n eq 4)) then begin 
        print, m, n
        psopen, 'ppig_maps' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)'), /color

        !p.multi=[0, 4, 4]
        loadct, 0;, 41 , /sil, file='~/colors1.tbl'

        restore, '~goldston/cv/pig_fit_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'
        restore, '~goldston/cv/pig_sp_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'
        extract_coords, '~goldston/Narrow/' + cname(m, n) + '_N.fits', ra, dec, vels,fits
        
        tlp = fits[*, *, 1024-50/delv:1023+50/delv]
        fits = 0
        vels = vels[ 1024-50/delv:1023+50/delv]
        szps = (size(pig_fit[0].postagestamp))[1]
        for j=0, n_elements(pig_sp)-1 do begin
            aa = pig_sp[j].ra
            dd = pig_sp[j].dec
            glactc, aa/15., dd, 2000, l, b, 1
            right = aa le mean(ra)
            top = dd ge mean(dec)
            x0 = right*256.
            x1 = right*256+255
            y0 = top*256
            y1 = top*256+255
            img =  total(tlp[x0:x1, y0:y1, min(where(pig_sp[j].spm eq 1)):max(where(pig_sp[j].spm eq 1))], 3) > 0.
            mximg = max(img)
            display, mximg-img, reform(ra[x0:x1, 0]), reform(dec[0, y0:y1]), xtitle='RA [deg]', ytitle='dec [deg]', title='pig ' + string(m, f='(I2.2)') + ', ' + string(n, f='(I2.2)') + ', ' + string(pig_sp[j].rm_num, f='(I3.3)') + ' l='  + string(l,  f='(I3.3)') + ', b=' + string(b,  f='(I3.3)')
            oplot, [aa - dp, aa -lp], [dd, dd]
            oplot, [aa + dp, aa + lp], [dd, dd]
            oplot, [aa, aa], [dd -dp, dd -lp]
            oplot, [aa, aa], [dd +dp, dd +lp]


            xx = interpol(findgen(512), ra[*, 0],pig_sp[j].ra)
            yy = interpol(findgen(512), dec[0, *],pig_sp[j].dec)
            subra = ra[(xx-15) > 0 :(xx+15) < 511, (yy-15) > 0 :(yy+15) < 511] 
            subdec = dec[(xx-15) > 0 :(xx+15) < 511, (yy-15) > 0 :(yy+15) < 511] 
           
            parm = pig_fit[j].twodfit
            parm[6:8] = 0.
            rslt = mpfit_2dgauss(parm, x=subra, model=model)

;            contour, model, subra[*, 0], subdec[0, *], level=max(model)/2., /overplot
            
            ;oplot,pig_fit[j].twodfit[4] +  sin(pig_fit[j].twodfit[6])*cos(theta)*pig_fit[j].twodfit[2] + cos(pig_fit[j].twodfit[6])*sin(theta)*pig_fit[j].twodfit[3], pig_fit[j].twodfit[5] + cos(pig_fit[j].twodfit[6])*cos(theta)*pig_fit[j].twodfit[2] - sin(pig_fit[j].twodfit[6])*sin(theta)*pig_fit[j].twodfit[3]
            
            display, pig_fit[j].postagestamp, reverse(findgen(szps)- (szps-1)/2.)/60.+pig_sp[j].ra, (findgen(szps)-(szps-1)/2.)/60.+pig_sp[j].dec, xtitle='RA [deg]', ytitle='dec [deg]'
            rslt = mpfit_2dgauss(parm, x=pig_fit[j].postagestamp, model=model)
            
            contour, model,reverse(findgen(szps)- (szps-1)/2.)/60.+pig_sp[j].ra, (findgen(szps)-(szps-1)/2.)/60.+pig_sp[j].dec, level=max(model)/2., /overplot

            plot, pig_sp[j].vels, pig_sp[j].sp, xtitle = 'V_LSR [km/s]', ytitle='T_B [K]', title='ON'
            loadct, 13,/sil     ;, file='~/colors1.tbl'
            wh = where(pig_sp[j].spm eq 1)
            oplot,  pig_sp[j].vels[wh], pig_sp[j].sp[wh], color=255
            loadct, 0;41 , /sil, file='~/colors1.tbl'

            plot, pig_sp[j].vels, pig_sp[j].sp- pig_sp[j].osp, xtitle = 'V_LSR [km/s]', ytitle='T_B [K]', title='ON-OFF'
            off_cent =  pig_fit[j].gfit[3] +  pig_fit[j].gfit[4]*pig_fit[j].gfit[1] + pig_fit[j].gfit[5]*pig_fit[j].gfit[1]^2 
            oplot, [pig_fit[j].gfit[1]-pig_fit[j].gfit[2],pig_fit[j].gfit[1]+pig_fit[j].gfit[2]], off_cent+[pig_fit[j].gfit[0], pig_fit[j].gfit[0]]/2.
            loadct, 13, /sil, file='~/colors1.tbl'
            oplot,  pig_sp[j].vels[wh], (pig_sp[j].sp- pig_sp[j].osp)[wh], color=255
            loadct, 0;41 , /sil, file='~/colors1.tbl'
        endfor
  
        psclose
;    endfor
;endfor



end
