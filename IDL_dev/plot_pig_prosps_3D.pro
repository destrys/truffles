pro plot_pig_prosps_3D,  outstring, embed, path1, path2, ghosts=ghosts,datafdir=datafdir,pigspout=pigspout,sporedir=sporedir,pigfit=pigfit,rmask=rmask
delv = 0.736122839/4.


;dataf_fn='dataf'+outstring+'.sav'
datac_fn=datafdir+'datacoords'+outstring+'.sav'
pig_sp_fn=pigspout+'pig_sp'+outstring+'.sav'
truff_fn = sporedir+'spores'+outstring+'.sav'
pig_fit_fn =pigfit+'pig_fit'+outstring+'.sav'
if embed eq 0 then img_fn =path1+'pig_plots'+outstring else img_fn =path1+'pig_plots'+outstring+'_embed'

rmasks=rmask+['rmask'+outstring+'_7.5.sav',$
'rmask'+outstring+'_18.5.sav',$
'rmask'+outstring+'_7.15.sav',$
'rmask'+outstring+'_18.15.sav']

dp = 0.5
lp = 1

theta = findgen(50)/49.*2*!pi

;       psopen, '/home/janagrc/ppig_maps_new_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)'), /color
psopen, img_fn, /color

!p.multi=[0, 3, 4]
loadct, 0
restore, pig_fit_fn ;'/home/janagrc/cv/pig_fit_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'
restore, pig_sp_fn
;restore, dataf_fn               ;'~destrys/truffles/cubes/dataf_268+10_again.sav'
restore, truff_fn
sz=size(fits)
fits = fits[20:sz[1]-21, 20:sz[2]-21,*] 


restore, datac_fn               ;'~destrys/truffles/cubes/datacoords_268+10_again.sav'
sz = size(fits)
szps = (size(pig_fit[0].postagestamp))[1]
whp = findgen(n_elements(pig_sp))
if not keyword_set(ghosts) then begin
   whp = whp[where(pig_sp.primary eq 1)]
   pig_fit= pig_fit[where(pig_sp.primary eq 1)]
   pig_sp= pig_sp[where(pig_sp.primary eq 1)]
endif
if embed then begin
   whp = whp[where(pig_sp.embed eq 1)]
   pig_fit= pig_fit[where(pig_sp.embed eq 1)]
   pig_sp= pig_sp[where(pig_sp.embed eq 1)]
endif else begin
   whp = whp[where(pig_sp.embed eq 0)]
   pig_fit= pig_fit[where(pig_sp.embed eq 0)]
   pig_sp= pig_sp[where(pig_sp.embed eq 0)]
endelse

;; NB: making a subselection and sorting!!
restore, 'labmask.sav'

find_lab_val, ll, bb, v, lab, pig_sp.ra, pig_sp.dec, pig_sp.v, klab

online=0
if online then begin
restore, 'labmask.sav'

find_lab_val, ll, bb, v, lab, pig_sp.ra, pig_sp.dec, pig_sp.v, klab

sub_sel= where(klab/pig_sp.mxsig gt 0.1)
pig_sp = pig_sp(sub_sel)
pig_fit = pig_fit(sub_sel)
pgsrt = reverse(sort(klab/pig_sp.mxsig))
pig_sp = pig_sp[pgsrt]
pig_fit = pig_fit[pgsrt]
endif
noedge = 1.
if noedge then begin

sub_sel= where( (pig_sp.v lt 500) and (pig_sp.v gt (-500)) )
pig_sp = pig_sp(sub_sel)
pig_fit = pig_fit(sub_sel)
klab=klab[sub_sel]
endif
;;; ++++++++

for j=0, n_elements(pig_sp)-1 do begin
           
           aa = pig_sp[j].ra
           dd = pig_sp[j].dec
           glactc, aa/15., dd, 2000, l, b, 1
           right = aa le mean(ra)
           top = dd ge mean(dec)
           x0 = right*sz[1]/2
           x1 = right*sz[1]/2+sz[1]/2-1
           y0 = top*sz[2]/2
           y1 = top*sz[2]/2+sz[2]/2-1
                             
           if (n_elements(where(pig_sp[j].spm eq 1)) eq 1) then begin
              img = fits[x0:x1, y0:y1]
            endif else begin
               img =  total(fits[x0:x1, y0:y1, min(where(pig_sp[j].spm eq 1)):max(where(pig_sp[j].spm eq 1))], 3)
            endelse
            mximg = max(img)

            if pig_sp[j].embed eq 0 then emb = '' else emb = 'Embedded, '
            if pig_sp[j].primary eq 1 then prm = 'Primary' else prm = 'Ghost of ' + string(pig_sp[j].ghostof,f='(I3.2)')
            
           ; display, mximg-img, reform(ra[x0:x1,0]),reform(dec[0, y0:y1]), xtitle='RA [deg]', ytitle='dec [deg]', title='Pig='+ string(j,f='(I3.3)') + ' Emb=' + string(pig_sp[j].embed,f='(I1.1)')+' Prim='+string(pig_sp[j].primary,f='(I2.2)')+' Gst=' + string(pig_sp[j].ghostof,f='(I3.2)') + ' Ker=' + string(pig_sp[j].kern_num,f='(I1.1)')
            display, -1*img, reform(ra[x0:x1,0]),reform(dec[0, y0:y1]), xtitle='RA [deg]', ytitle='dec [deg]', title='Pig='+ string(whp[j],f='(I3.3)') +', ' + emb + prm + ', ' + 'Ker=' + string(pig_sp[j].kern_num,f='(I1.1)')
            oplot, [aa - dp, aa -lp], [dd, dd]
            oplot, [aa + dp, aa + lp], [dd, dd]
            oplot, [aa, aa], [dd -dp, dd -lp]
            oplot, [aa, aa], [dd +dp, dd +lp]

            ; psarea
            oplot, [aa - 15.5/60.,aa - 15.5/60.,aa + 15.5/60.,aa + 15.5/60.,aa - 15.5/60.],[dd + 15.5/60.,dd - 15.5/60.,dd - 15.5/60.,dd + 15.5/60.,dd + 15.5/60.]

            print,pig_sp[j].ra
            xx = interpol(findgen(sz[1]), ra[*, 0],pig_sp[j].ra)
            print,xx
            yy = interpol(findgen(sz[2]), dec[0, *],pig_sp[j].dec)
            ;print,(xx-15) > 0,(xx+15) < 511,(yy-15) < 0, (yy+15) < 511
            subra = ra[(xx-15) > 0 :(xx+15) < (sz[1]-1), (yy-15) > 0 :(yy+15) < (sz[2]-1)] 
            subdec = dec[(xx-15) > 0 :(xx+15) < (sz[1]-1), (yy-15) > 0 :(yy+15) < (sz[2]-1)] 
           
            parm = pig_fit[j].threedfit
            parm[8:*] = 0.
            data = fltarr(pig_fit[j].px[1]-pig_fit[j].px[0]+1, pig_fit[j].py[1]-pig_fit[j].py[0]+1, pig_fit[j].pv[1]-pig_fit[j].pv[0]+1) 
            rslt = mpfit_3DG(parm, x=data, model=model)
            ramod = ra[pig_fit[j].px[0]:pig_fit[j].px[1], pig_fit[j].py[0]:pig_fit[j].py[1]]
            decmod = dec[pig_fit[j].px[0]:pig_fit[j].px[1], pig_fit[j].py[0]:pig_fit[j].py[1]]
            model = total(model, 3)
mmax=max(-1*pig_fit[j].postagestamp*pig_fit[j].psrmask)
mmin=min(-1*pig_fit[j].postagestamp*pig_fit[j].psrmask)
mdiff=mmax-mmin

            display, -1*pig_fit[j].postagestamp, reverse(findgen(szps)- (szps-1)/2.)/60.+pig_sp[j].ra, (findgen(szps)-(szps-1)/2.)/60.+pig_sp[j].dec, xtitle='RA [deg]', ytitle='dec [deg]', title='sig= '+string(pig_sp[j].mxsig)+ ', LABT = ' + string(klab[j]),max=mmax+0.5*mdiff,min=mmin-0.5*mdiff
            if (pig_sp[j].primary eq 1) and (pig_sp[j].embed eq 0) then loadct, 13, /sil
;            contour, model, reform(ramod[*, 0]), reform(decmod[0, *]), level=max(model)/2., /overplot, c_thick=2
            contour, pig_fit[j].psrmask, reverse(findgen(szps)- (szps-1)/2.)/60.+pig_sp[j].ra, (findgen(szps)-(szps-1)/2.)/60.+pig_sp[j].dec, level=[0.5], /overplot, c_thick=2
            loadct, 0, /sil
                   
; new version that shows the spectrum
; inside the FWHM ellipse and fit
            scmod = fits[pig_fit[j].px[0]:pig_fit[j].px[1], pig_fit[j].py[0]:pig_fit[j].py[1], *]
            szsm = size(scmod)
            modmask = model*0.
            mmmod = minmax(model)
            modmask(where(model ge (mmmod[1]+mmmod[0])/2.)) = 1.
;average spectrum within ellipse
            onsp = total(total(scmod*rebin(reform(modmask, szsm[1], szsm[2], 1), szsm[1], szsm[2], szsm[3]), 1), 1)/total(modmask)
            parm = pig_fit[j].threedfit
            rslt = mpfit_3DG(parm, x=data, model=modeltilt)
            onspmod = total(total(modeltilt*rebin(reform(modmask, szsm[1], szsm[2], 1), szsm[1], szsm[2], pig_fit[j].pv[1]- pig_fit[j].pv[0]+1), 1), 1)/total(modmask)

            wh = where(pig_sp[j].spm eq 1)
            vrange = fltarr(2)
            vrange[0] = min((pig_sp[j].vels)[wh])-30
            vrange[1] = max((pig_sp[j].vels)[wh])+30            
            plot, pig_sp[j].vels, onsp, xtitle = 'V_LSR [km/s]', ytitle='T_B [K]', title='fit ellipse spectrum',/ystyle, xra=vrange, yrange=[(min(onspmod) - 3) > 0, max(onspmod)+1]
            loadct, 13,/sil     ;, file='~/colors1.tbl'
            oplot,vels[pig_fit[j].pv[0]:pig_fit[j].pv[1]], onspmod, color=255, thick=4
            loadct, 0, /sil

            if (min(pig_sp[j].osp) eq !VALUES.F_NAN) then begin
                plot, pig_sp[j].vels, pig_sp[j].sp- pig_sp[j].osp, xtitle = 'V_LSR [km/s]', ytitle='T_B [K]', title='ON-OFF'
                off_cent =  pig_fit[j].gfit[3] +  pig_fit[j].gfit[4]*pig_fit[j].gfit[1] + pig_fit[j].gfit[5]*pig_fit[j].gfit[1]^2 
                oplot, [pig_fit[j].gfit[1]-pig_fit[j].gfit[2],pig_fit[j].gfit[1]+pig_fit[j].gfit[2]], off_cent+[pig_fit[j].gfit[0], pig_fit[j].gfit[0]]/2.
                loadct, 13, /sil, file='~/colors1.tbl'
                oplot,  pig_sp[j].vels[wh], (pig_sp[j].sp- pig_sp[j].osp)[wh], color=255, thick=3
                loadct, 0;41 , /sil, file='~/colors1.tbl'

             endif
              
        endfor
  
        psclose
;    endfor
;endfor


end
