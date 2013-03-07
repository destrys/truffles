pro cv_all_2d, i, j, fake=fake,dataf=dataf,fwhm=fwhm,ra=ra,dec=dec,cv=cv
;+
; NAME: CV_ALL_2D
;
; PURPOSE:
; This is step one in the code to look for 'pigs' in the GALFA-HI
; data set, by convolving a reduced velocity range cube with 
; a matched filter function
; It also allows the injection of 'fake' pigs with the injector
; function
;
; This is the 2D version.
;
; INPUTS:
;        I: RA  (degrees) ?
;        J: Dec (degrees) ?
;
; KEYWORDS:
;        FAKE: Triggers Guinea Pig Injection
;        DATAF: Filename of Data Cube. If set, I & J are used for 
;               the output file only.
;        FWHM: FUll Width Half Max of Convolving Gaussian.
;              Default=3.4
;
; OUTPUTS:
;        Writes save file to : '~/cv/cv_i_j.sav'
;                         or : '~/cv/cv_inj_i_j.sav' if /fake is set
;        This save file includes:
;                         RA: 512 x 512 array (RA values)
;                        DEC: 512 x 512 array (DEC values)
;                         CV: 512 x 512 x 1087 array (Convolved Cube)
; CALLS:
;        GROMASK (gsr)
;        EXTRACT_COORDS (gsr)
;        CNAME (gsr)
;        INJECTOR (optional)
;
; CALLED BY:
;        BLANKETS.PRO
;
; HISTORY:
;        thepast:      Written (JEGP)
;        Jan. 9, 2010: Documented (Destry)
;        Jan.10, 2010: Negative points introduced by READFITS replaced
;                      with NANs (Destry)
;        Jan.11, 2010: Added DATAF Keyword (Destry)
;        Jan.16, 2010: Adpated from CV_ALL, changed to 2D convols. (Destry)
;        Jan.17, 2010: Added FWHM Keyword (Destry)
;
;-
if keyword_set(fake) then inj = 1 else inj = 0

;;;;;;;;;;;;;;;;;;;; Set Initial Parameters ;;;;;;;;;;;;;;;;;;;;;;;;;
dv = 0.736122839/4.     ; channel width km/s  ??
nv = 31                 ; ??
voff =100               ; Velocity range to search


;;;;;;;;;;;;;;;;;;;; Build Filtering Function ;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;; Need to Pull Parameters Out of Here ;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;; So They Can be Varied Easily ;;;;;;;;;;;;;;;;;;;
xx = rebin(reform(findgen(21)-10, 21, 1), 21, 21)
yy = rebin(reform(findgen(21)-10, 1, 21), 21, 21)
vv = rebin(reform((findgen(31)-15)*dv, 1, 1, 31), 21, 21)

if ~keyword_set(fwhm) then fwhm=3.4

posgau = 2*exp( (-0.5)*(xx/(fwhm/(2*sqrt(2*alog(2)))))^2)*exp( (-0.5)*(yy/(fwhm/(2*sqrt(2*alog(2)))))^2)*(vv/(2.0/(2*sqrt(2*alog(2)))))^2)

wid=1.5

neggau=  exp( (-0.5)*(xx/(wid*fwhm/(2*sqrt(2*alog(2)))))^2)*exp( (-0.5)*(yy/(wid*fwhm/(2*sqrt(2*alog(2)))))^2)*(vv/(wid*2.0/(2*sqrt(2*alog(2)))))^2)

neggau = neggau*total(posgau)/total(neggau)

totgau = posgau-neggau
totgau = totgau/max(totgau)


;;;; These are Inputs for Injector, Consider Moving Down ;;;;;;;;;;;;
v_rng = [-voff, voff]
sz_rng = [3, 6]
as_rng = [1,2]
a_rng = [0, 3.0]
pa_rng = [0, 180]
fwhm_rng = [1, 5]

;;;;;;;;;;;;;;;;; Read In FITS File ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;; Consider Separating this ;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;; Into another function ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        extract_coords, '~/Narrow/' + cname(i,j) + '_N.fits', ra, dec, vels, fits

IF keyword_set(dataf) then begin
extract_coords,dataf,ra,dec,vels,fits
endif else begin
        ;; Changed to open Wide Cube instead of Narrow..... (Jan 9 - Destry)
extract_coords,cname(i,j)+'_W.fits',ra,dec,vels,fits
endelse
;;;;;;;;;;;;;;; Fix Negative Values ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fits[where(fits LT 0)]=!values.f_nan

;;;;;;;;;;;;;;;;; Slice out Channels of interest ;;;;;;;;;;;;;;;;;;;;
;tlpall = fits[*, *, 1024-voff/dv:1023+voff/dv ]
;vels = vels[1024-voff/dv:1023+voff/dv ]
ra_rng = [ra[16, 0], ra[511-16, 0]]
dec_rng = [dec[0, 16], dec[0, 511-16]]

;removed slicing -Jan 14 (destry)
tlpall=fits
fits=0    ; ??

;;;;;;;;;;;;;;;;;;;;;;; Optional Injection ;;;;;;;;;;;;;;;;;;;;;;;;;;
if inj then injector, tlpall, vels, ra, dec, 1000, v_rng, ra_rng, dec_rng, sz_rng, as_rng, a_rng, pa_rng, fwhm_rng, g_inj


;;;;;;;;;;;;;;;;;;;;;;;;; Convolve It! ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;sz=size(tlpall)
;cv=tlpall*0
;for inx=0,sz[3]-1 do begin
;cv[*,*,inx] = convol(tlpall[*,*,inx], totgau,/nan)
;endfor

cv = convol(tlpall, totgau,/nan)
;;;;;;;;;;;;;;;;;;;;;;; Save Output ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;; Consider Removing this for speed ;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;; Maybe Make it a keyword option ;;;;;;;;;;;;;;;;;;;;;;
if inj then save, ra, dec, cv, v_rng, ra_rng, dec_rng, sz_rng, as_rng, a_rng, pa_rng, fwhm_rng, g_inj , f= '~/cv/cv_inj_' + string(i, f='(I2.2)') +'_' + string(j, f='(I2.2)')+ '.sav' else save, dec, ra, cv, f= '~/cv/cv_' + string(i, f='(I2.2)') +'_' + string(j, f='(I2.2)')+ '.sav'


tlpall=0.  ;??



end

@
