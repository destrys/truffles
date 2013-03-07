pro pig_props_3d, outstring, writedir=writedir,readdir=readdir,datafdir=datafdir,pigspout=pigspout,sporedir=sporedir,rmask=rmask

; step 4 in the pig-searching code. This determines the properties
; of the pigs, using the mpfit_2dgauss code and mpfit.
; again, fk = '_inj' uses the injected pigs code.

if ~keyword_set(readdir) then readdir=''

voff=100
delv = 0.736122839
fk=''
dataf_fn=datafdir+'dataf'+outstring+'.sav'
datac_fn=datafdir+'datacoords'+outstring+'.sav'
pig_sp_fn=pigspout+'pig_sp'+outstring+'.sav'
if keyword_set(writedir) then pig_fit_fn =writedir +'/pig_fit'+outstring+'.sav' else pig_fit_fn ='pig_fit'+outstring+'.sav'
truff_fn = sporedir+'spores'+outstring+'.sav'
rmasks=[rmask+'rmask'+outstring+'_7.5.sav',$
rmask+'rmask'+outstring+'_18.5.sav',$
rmask+'rmask'+outstring+'_7.15.sav',$
rmask+'rmask'+outstring+'_18.15.sav']

restore, pig_sp_fn 
szv = 2048


restore, truff_fn

restore, datac_fn

sz=size(ra)
if ~keyword_set(keepedges) then begin
ra=ra[30:sz[1]-30-1,30:sz[2]-30-1]
dec=dec[30:sz[1]-30-1,30:sz[2]-30-1]
endif
sz=size(ra)
nx=41
ra=ra[nx/2:sz[1]-nx/2-1,nx/2:sz[2]-nx/2-1]
dec=dec[nx/2:sz[1]-nx/2-1,nx/2:sz[2]-nx/2-1]
fits=fits[nx/2:sz[1]-nx/2-1,nx/2:sz[2]-nx/2-1,*]






sz = size(fits)

       pssize = 15

       pig_fit = replicate({mxsig:0.,embed:0,ghostof:0,primary:0,kern_num:0,rm_num:0, px:fltarr(2), py:fltarr(2), pv:fltarr(2),threedfit:fltarr(12),threedfiterr:fltarr(12), ra:0., dec:0., v:0., postagestamp:fltarr(pssize*2+1, pssize*2+1),stampmax:fltarr(pssize*2+1,pssize*2+1), specint:fltarr(sz[3]),specmax:fltarr(sz[3]),psrmask:fltarr(pssize*2+1, pssize*2+1), vel:0.,gflux:0.,sflux:0.,FWHM:0.,T_bmax:0.,N_HI:0.,M_HI:0., raROI:0., decROI:0., lROI:0., bROI:0., vLSRROI:0., vGSRROI:0., tpeakROI:0., NHpeakROI:0., VFWHMROI:0.,VFWHMROI_cent:0., PAROI:0., totalfluxROI:0., FWHMmajROI:0., FWHMminROI:0.}, n_elements(pig_sp))
       
        ; get set up to look through the kernals
        kern = -1
        for j=0, n_elements(pig_sp)-1 do begin
           loop_bar, j,  n_elements(pig_sp)
           kcur = pig_sp[j].kern_num
           if kern ne kcur then begin
              kern = kcur
              restore, rmasks[kern]
           endif
           ; find the cloud, masked by the fit
           padp = 3.
           padv=3.
           ; area to look for cloud, with a little padding
           px = [(pig_sp[j].px[0]-padp) > 0,(pig_sp[j].px[1]+padp) < (sz[1]-1)]
           py = [(pig_sp[j].py[0]-padp) > 0,(pig_sp[j].py[1]+padp) < (sz[2]-1)]
           pv = [(pig_sp[j].pv[0]-padv) > 0,(pig_sp[j].pv[1]+padv) < (sz[3]-1)]
           ; get the cloud, only rmask area has data
           cld = (rmask[px[0]:px[1],py[0]:py[1],pv[0]:pv[1]] eq  pig_sp[j].rm_num)*fits[px[0]:px[1],py[0]:py[1],pv[0]:pv[1]]
           whinfn = where(finite(cld) eq 0, ct)
           if ct ne 0 then cld[whinfn] = 0.
           ROI = (rmask[px[0]:px[1],py[0]:py[1],pv[0]:pv[1]] eq  pig_sp[j].rm_num)
           errmask = 1-(rmask[px[0]:px[1],py[0]:py[1],pv[0]:pv[1]] eq  pig_sp[j].rm_num)
           qmask = where(rmask[px[0]:px[1],py[0]:py[1],pv[0]:pv[1]] eq  pig_sp[j].rm_num)
           ; try destry's histogram trick
           dx = px[1]-px[0] + 1
           dy = py[1]-py[0]+1
           for q =0, 2 do begin
              h = Histogram([qmask,qmask+1,qmask-1,qmask+dx,qmask-dx,qmask+dx*dy,qmask-dx*dy], OMIN=omin)
              qmask = Where(h GT 0) + omin
           endfor
           errmask2 = errmask*0.
           errmask2[qmask] = 1.
           
           szcld= size(cld)
           ; axes over cloud area
           racld = ra[px[0]:px[1],py[0]:py[1]]
           deccld = dec[px[0]:px[1],py[0]:py[1]]
           vcld = vels[pv[0]:pv[1]]
           fa = {x:cld, err:(errmask*100.+1)}
           ; position of the peak, and value
           xm =  interpol(findgen(px[1]-px[0]+1), racld[*, 0],pig_sp[j].ra)
           ym =interpol(findgen(py[1]-py[0]+1), deccld[0, *],pig_sp[j].dec)
           vm = interpol(findgen(pv[1]-pv[0]+1), vcld,pig_sp[j].v)
           tm= cld[xm, ym, vm]

           parinfo = replicate({fixed:0., limited:fltarr(2), limits:fltarr(2)}, 12)
           ; no tip/tilt to start
           parinfo[9].fixed = 1.
           parinfo[10].fixed = 1.
           parinfo[11].fixed = 1.
           parinfo[0].limited=[1, 0]
           parinfo[1].limited=[1, 1]
           parinfo[2].limited=[1, 1]
           parinfo[3].limited=[1, 1]
           parinfo[8].limited=[1, 0]
           parinfo[8].limits=[0, 1000]

           parinfo[0].limits=[0, 100]
         ;  parinfo[1].limits=[(px[1]-px[0])*0.33, (px[1]-px[0])*0.66]
         ;  parinfo[2].limits=[(py[1]-py[0])*0.33, (py[1]-py[0])*0.66]
        ;   parinfo[3].limits=[(pv[1]-pv[0])*0.33, (pv[1]-pv[0])*0.66]
           parinfo[1].limits=[0, (px[1]-px[0])]
           parinfo[2].limits=[0, (py[1]-py[0])]
           parinfo[3].limits=[0, (pv[1]-pv[0])]

           params1 = mpfit('mpfit_3DG', [tm, xm, ym, vm,(px[1]-px[0])/4.,(py[1]-py[0])/4 , (pv[1]-pv[0])/4, 0., 0., 0., 0., 0.], parinfo=parinfo, functargs=fa, perror=pe, /quiet)
           
           ; now use output to fit the real data with no mask
           
           ; add tip/tilt 
           parinfo[9].fixed = 0.
           parinfo[10].fixed = 0.
           parinfo[11].fixed = 0.
           fa.x = fits[px[0]:px[1],py[0]:py[1],pv[0]:pv[1]]
           fa.err = (1-errmask2)*100.+1.
           params= mpfit('mpfit_3DG', params1, parinfo=parinfo, functargs=fa, perror=pe, /quiet)
           pargauss = params
           pargauss[8:*]=0

           rslt = mpfit_3DG(pargauss,x=cld, model=model)
           pig_fit[j].px = px
           pig_fit[j].py = py
           pig_fit[j].pv = pv           
           pig_fit[j].rm_num = pig_sp[j].rm_num
	   pig_fit[j].kern_num = pig_sp[j].kern_num
	   pig_fit[j].mxsig = pig_sp[j].mxsig
	   pig_fit[j].embed=pig_sp[j].embed
	   pig_fit[j].ghostof=pig_sp[j].ghostof
	   pig_fit[j].primary=pig_sp[j].primary
	   pig_fit[j].threedfit = params
;           pig_fit[j].threedfiterr = pe
           pig_fit[j].ra = racld[params[1], 0]
           pig_fit[j].dec = deccld[0, params[2]]
           pig_fit[j].v =vcld[params[3]]
           pig_fit[j].T_bmax = tm
	   pig_fit[j].specmax=fits[xm+px[0],ym+py[0],*]
           pig_fit[j].specint=pig_sp[j].sp
	   pig_fit[j].N_HI =1.82e18*max(total(model, 3))*delv
           ; ROI calcuations
           pig_fit[j].tpeakROI = max(cld)
           pig_fit[j].NHpeakROI = max(total(cld, 3))*1.82e18*delv
           pig_fit[j].raROI = pig_sp[j].ra
           pig_fit[j].decROI = pig_sp[j].dec
           pig_fit[j].vlsrROI  = pig_sp[j].v
           glactc, pig_fit[j].raROI/15., pig_fit[j].decROI, 2000, l, b, 1
           pig_fit[j].lROI = l
           pig_fit[j].bROI = b
           pig_fit[j].vgsrROI = pig_fit[j].vlsrROI + sin(l*!pi/180.)*cos(b*!pi/180.)*220. 

           pig_fit[j].totalfluxROI = total(cld)*delv ; arcmin^2 km/s K
           vvvcld = rebin(reform(vcld, 1, 1, pv[1]-pv[0]+1), px[1]-px[0]+1, py[1]-py[0]+1, pv[1]-pv[0]+1)
           mmv = minmax(vvvcld(where(cld ge pig_fit[j].tpeakROI/2.)))
           pig_fit[j].vfwhmROI = mmv[1]-mmv[0]
           pig_fit[j].vfwhmROI_cent = max(vcld(where(cld[xm, ym, *] ge max(cld[xm, ym, *])/2.))) - min(vcld(where(cld[xm, ym, *] ge  max(cld[xm, ym, *])/2.)))

           ; half-maximum mask
           hmmask = cld*0
           hmmask(where(cld ge max(cld, /nan)/2.)) = 1.
           hmmask2d = total(hmmask, 3) ge 1.
           rasroi = racld(where(hmmask2d eq 1))
           decsroi = deccld(where(hmmask2d eq 1))
           findaxes, rasroi, decsroi, PA, majax, minax
           pig_fit[j].PAROI= PA
           pig_fit[j].FWHMmajROI= majax
           pig_fit[j].FWHMminROI= minax

            ; Get the postage stamp
            tlpsz = size(fits[*,*,where(pig_sp[j].spm eq 1)],/n_dimensions)
            if (tlpsz eq 2) then begin
               img_on = fits[*,*,where(pig_sp[j].spm eq 1)]
            endif else begin
               img_on = total(fits[*, *, where(pig_sp[j].spm eq 1)], 3)/total(pig_sp[j].spm) ; in K km/s
            endelse
	    d_off = 5
            off_plus = where(pig_sp[j].spm eq 1) + total(pig_sp[j].spm) + d_off
            off_minus = where(pig_sp[j].spm eq 1) - total(pig_sp[j].spm) - d_off
            nw = n_elements(where(pig_sp[j].spm eq 1))
            whbdmn = where(off_minus lt 0, ctm)
            if (ctm ne 0) and (ctm lt nw) then off_minus = off_minus[where(off_minus ge 0)]
            whbdpl = where(off_plus gt (szv-1), ctp)
            if (ctp ne 0) and (ctp lt nw) then off_plus = off_plus[where(off_plus le (szv-1))]

            if (ctp eq nw) then begin
              if (n_elements(off_minus) gt 1) then begin
                img_off =  total(fits[*, *, off_minus], 3)/n_elements(off_minus)
              endif else begin
                img_off = fits[*, *, off_minus]/n_elements(off_minus)
              endelse
            endif
 
            if (ctm eq nw) then begin
              if (n_elements(off_plus) gt 1) then begin
                img_off =  total(fits[*, *, off_plus], 3)/n_elements(off_plus)
              endif else begin
                img_off = fits[*,*,off_plus]/n_elements(off_plus)
              endelse            
            endif

            if (ctm ne nw) and (ctp ne nw) then img_off =  total(fits[*, *, [off_plus,off_minus]], 3)/n_elements([off_plus,off_minus])

;            img = img_on-img_off
	    img=img_on
            xx = interpol(findgen(sz[1]), ra[*, 0],pig_sp[j].ra)
            if (xx lt 0) then begin
              xx=0
            endif
            yy = interpol(findgen(sz[2]), dec[0, *],pig_sp[j].dec)
            if (yy lt 0) then begin
              yy=0
            endif

            img_bord = fltarr(sz[1]+pssize*2, sz[2]+pssize*2)
            img_bord[pssize:pssize+sz[1]-1, pssize:pssize+sz[2]-1] = img
            psimg = img_bord[xx :xx+pssize*2, yy:yy+pssize*2] 
            pig_fit[j].postagestamp = psimg
 
	    img_bord*= 0
	    img_bord[pssize:pssize+sz[1]-1,pssize:pssize+sz[2]-1] = fits[*,*,vm+pv[0]]
	    pig_fit[j].stampmax=img_bord[xx:xx+pssize*2,yy:yy+pssize*2]

           msk_rm = fltarr(sz[1], sz[2])
            ; where the cloud is in the sub region
            mskcld = total(cld, 3) ne 0
            msk_rm[px[0]:px[1],py[0]:py[1]] = mskcld
            mask_bord = fltarr(sz[1]+pssize*2, sz[2]+pssize*2)
            mask_bord[pssize:pssize+sz[1]-1, pssize:pssize+sz[2]-1] = msk_rm
            maskimg = mask_bord[xx :xx+pssize*2, yy:yy+pssize*2] 
            pig_fit[j].psrmask =maskimg
        endfor
        
        save, pig_fit, f=pig_fit_fn

end
