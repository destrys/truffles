pro pig_match,clouds,cube=cube,pig_sp=pig_sp,match=match
;+
; NAME:
;       PIG_MATCH
;
; PURPOSE:
;       Match up the pigs and get efficiency stats
;
;
;
;
;
;
;
;
;
;
;
;
;-

;;;;; Restore Data ;;;;
restore,'~/cv/REAL_15_03.sav',/verb

extract_coords,cube,ra,dec,vels,fits

;;;;;; MAtch em' up, if requested ;;;;;;;
if keyword_set(match) then begin
found=matchy(g_inj,pig_sp,_extra=_extra)
endif else begin
restore,clouds,/verb
endelse

;;;;; Find backgroud values ;;;;;
decsz=n_elements(dec[0,*])
rasz=n_elements(ra[*,0])
g_injsz=n_elements(g_inj)
velsz=n_elements(vels)
foundsz=n_elements(found)

ndec=rebin(dec[0,*],g_injsz,decsz)
nra=rebin(reform(ra[*,0],1,rasz),g_injsz,rasz)
nvels=rebin(reform(vels,1,velsz),g_injsz,velsz)
jra=rebin(g_inj.ra,g_injsz,rasz)
jdec=rebin(g_inj.dec,g_injsz,decsz)
jvel=rebin(g_inj.v,g_injsz,velsz)
scrap=min(abs(ndec-jdec),dec_ind,dim=2,/nan)
scrap=min(abs(nra-jra),ra_ind,dim=2,/nan)
scrap=min(abs(nvels-jvel),vel_ind,dim=2,/nan)

dec_ind/=g_injsz
ra_ind/=g_injsz
vel_ind/=g_injsz

g_inj_bg=fits[ra_ind,dec_ind,vel_ind]

ndec=rebin(dec[0,*],foundsz,decsz)
nra=rebin(reform(ra[*,0],1,rasz),foundsz,rasz)
nvels=rebin(reform(vels,1,velsz),foundsz,velsz)
jra=rebin(found.ra,foundsz,rasz)
jdec=rebin(found.dec,foundsz,decsz)
jvel=rebin(found.v,foundsz,velsz)
scrap=min(abs(ndec-jdec),dec_ind,dim=2,/nan)
scrap=min(abs(nra-jra),ra_ind,dim=2,/nan)
scrap=min(abs(nvels-jvel),vel_ind,dim=2,/nan)

dec_ind/=foundsz
ra_ind/=foundsz
vel_ind/=foundsz

found_bg=fits[ra_ind,dec_ind,vel_ind]








stop
;;;;;;; Brightness ;;;;;
psopen,'eff_a',/encapsulate,/times,/inches,xsize=7,ysize=4,/bold
effhist,g_inj.a,found.a,binsize=0.05,xtit='Peak (K)',ytit='Recovered Fraction',yran=[0,1],min=0,xthick=5,ythick=5,thick=5,font=0,fillcolor=100
psclose

psopen,'eff_vlsr',/encapsulate,/times,/inches,xsize=7,ysize=4,/bold
effhist,g_inj.v,found.v,binsize=40.,xtit=textoidl('V_{LSR} (km s^{-1}))',font=0),ytit='Recovered Fraction',yran=[0,1],min=-400.,xthick=5,ythick=5,thick=5,font=0,fillcolor=100
psclose

psopen,'eff_vfwhm',/encapsulate,/times,/inches,xsize=7,ysize=4,/bold
effhist,g_inj.fwhm,found.fwhm,binsize=1.,xtit=textoidl('FWHM (km s^{-1})',font=0),ytit='Recovered Fraction',yran=[0,1],min=0.,xthick=5,ythick=5,thick=5,font=0,fillcolor=100
psclose

psopen,'eff_size',/encapsulate,/times,/inches,xsize=7,ysize=4,/bold
effhist,g_inj.sz,found.sz,binsize=1.,xtit='Size (arcsec)',ytit='Recovered Fraction',yran=[0,1],min=0.,xthick=5,ythick=5,thick=5,font=0,fillcolor=100
psclose

psopen,'eff_asp',/encapsulate,/times,/inches,xsize=7,ysize=4,/bold
effhist,g_inj.as,found.as,binsize=0.2,xtit='Aspect Ratio',ytit='Recovered Fraction',yran=[0,1],min=0.,xthick=5,ythick=5,thick=5,font=0,fillcolor=100
psclose

psopen,'back_tot',/encapsulate,/times,/inches,xsize=7,ysize=4,/bold
dshist,g_inj_bg,binsize=0.3,xtit='Bachgorund (K)',ytit='# of Spam clouds', yran=[0,600],min=-1.,xthick=5,ythick=5,thick=5,font=0,fillcolor=100
psclose

stop
psopen,'eff_back',/encapsulate,/times,/inches,xsize=7,ysize=4,/bold
effhist,g_inj_bg,found_bg,binsize=0.1,xtit='Background (K)',ytit='Recovered Fraction',yran=[0,1],min=-1.,xthick=5,ythick=5,thick=5,font=0,fillcolor=100
psclose

stop



END
