pro cv_2d,dataf=dataf,fwhm=fwhm,v_fwhm=v_fwhm,cv=cv,outfile=outfile,fits=fits,chin=chin,noisescale=noisescale,wtout=wtout,_extra=_extra
;+
; NAME: CV_2D
;
; PURPOSE:
; This is step one in the code to look for 'pigs' in the GALFA-HI
; data set, by convolving a reduced velocity range cube with 
; a matched filter function
;
; This is the 2D version.
;
; INPUTS:
;        I: RA  (degrees) ?
;        J: Dec (degrees) ?
;
; KEYWORDS:
;        DATAF: Filename of Data Cube. 
;            
;        FWHM: Spatial Full Width Half Max of Convolving Gaussian.
;              Default=7'
;
;        V_FWHM: Velocity Full Width Half Max of Convolving Gaussian.
;		Default=5 km/s
;
;
; OUTPUTS:
;        Writes save files to OUTFILE and WTOUT
;                       
;        OUTFILE includes:
;                         RA: 512 x 512 array (RA values)
;                        DEC: 512 x 512 array (DEC values)
;                         CV: 512 x 512 x 2048 array (Convolved Cube)
;        WTOUT includes:
;			  WT: Output from WTFIND.PRO
;
;
; CALLS:
;        WTFIND (truffles)
;
; CALLED BY:
;        TRUFFLES.PRO
;
; HISTORY:
;        thepast:      Written (JEGP)
;        Jan. 9, 2010: Documented (Destry)
;        Jan.10, 2010: Negative points introduced by READFITS replaced
;                      with NANs (Destry)
;        Jan.11, 2010: Added DATAF Keyword (Destry)
;        Jan.16, 2010: Adpated from CV_ALL, changed to 2D convols. (Destry)
;        Jan.17, 2010: Added FWHM Keyword (Destry)
;        Feb.13, 2010: General Housekeeping and inserted WTFIND (Destry)
;-
;;;;;;;;;;;;;;;;;;;; Set Initial Parameters ;;;;;;;;;;;;;;;;;;;;;;;;;

dv = 0.736122839     ; channel width km/s
wid=1.2                 ; Width of Negative Gaussian relative to Positive Gaussian

;;;;;;;;;;;;;;;;;;;;;;;;;; Set Scaling ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if ~keyword_set(fwhm) then fwhm=7.0
if ~keyword_set(v_fwhm) then v_fwhm=5.0
scale=1
if (fwhm EQ 7) and (v_fwhm EQ 5) then scale=0.120881
if (fwhm EQ 7) and (v_fwhm EQ 15) then scale=0.106187
if (fwhm EQ 18) and (v_fwhm EQ 5) then scale=0.116593
if (fwhm EQ 18) and (v_fwhm EQ 15) then scale=0.102864
if (fwhm NE 18) and (fwhm NE 7) then print,'INVALID Kernel Size, Scaling will be incorrect'
if (v_fwhm NE 5) and (v_fwhm NE 15) then print,'INVALID Kernel FWHM, Scaling will be incorrect'

;;;;;;;;;;;;;;;;;;;; Build Filtering Function ;;;;;;;;;;;;;;;;;;;;;;;

nx=41    ; Spatial Size of Convovling box
nnv=41   ; Velocity Size of Convolving box

xxx = rebin(reform(findgen(nx)-nx/2, nx, 1, 1), nx, nx, nnv)
yyy = rebin(reform(findgen(nx)-nx/2, 1, nx, 1), nx, nx, nnv)
vvv = rebin(reform((findgen(nnv)-nnv/2)*dv, 1, 1, nnv), nx, nx, nnv)

d2posgau = 2*exp( (-0.5)*(xxx[*,*,0]/(fwhm/(2*sqrt(2*alog(2)))))^2)*exp( (-0.5)*(yyy[*,*,0]/(fwhm/(2*sqrt(2*alog(2)))))^2)
v2posgau = 2*exp( (-0.5)*(vvv[0,0,*]/(v_fwhm/(2*sqrt(2*alog(2)))))^2)

d2neggau = exp( (-0.5)*(xxx[*,*,0]/(wid*fwhm/(2*sqrt(2*alog(2)))))^2)*exp( (-0.5)*(yyy[*,*,0]/(wid*fwhm/(2*sqrt(2*alog(2)))))^2)
v2neggau = exp( (-0.5)*(vvv[0,0,*]/(wid*v_fwhm/(2*sqrt(2*alog(2)))))^2)

d2neggau *= total(d2posgau)/total(d2neggau)
v2neggau *= total(v2posgau)/total(v2neggau)

dscale=max(d2posgau-d2neggau)
vscale=max(v2posgau-v2neggau)

d2posgau /=dscale
d2neggau /=dscale
v2posgau /=vscale
v2neggau /=vscale

;;;;;;;;;;;;;;;;; Build 3d kernel for wtfind.pro ;;;;;;;;;;;;;;;;;;;;;

posgau = 2*exp( (-0.5)*(xxx/(fwhm/(2*sqrt(2*alog(2)))))^2)*exp( (-0.5)*(yyy/(fwhm/(2*sqrt(2*alog(2)))))^2)*exp( (-0.5)*(vvv/(v_fwhm/(2*sqrt(2*alog(2)))))^2)

neggau=  exp( (-0.5)*(xxx/(wid*fwhm/(2*sqrt(2*alog(2)))))^2)*exp( (-0.5)*(yyy/(wid*fwhm/(2*sqrt(2*alog(2)))))^2)*exp( (-0.5)*(vvv/(wid*v_fwhm/(2*sqrt(2*alog(2)))))^2)

sccc=total(posgau)/total(neggau)

neggau = neggau*total(posgau)/total(neggau)

totgau = posgau-neggau
totgau = totgau/max(totgau)


posgau=0b
neggau=0b
xxx=0b
yyy=0b
vvv=0b

;;;;;;;;;;;;;;;;; Read In FITS File ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;; Consider Separating this ;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;; Into another function ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if ~keyword_set(fits) then begin
	IF keyword_set(dataf) then begin
		restore,dataf,/verb
	endif else print,'Need to set dataf keyword so I can find file!'
endif

;;;;;;;;;;;;;;; Fix Negative Values ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

bad_ind=where((fits LT -30) or (finite(fits) eq 0),bct)
if bct GT 0 then fits[bad_ind]=!values.f_nan


;;;;;;;;;;;;;;;; Find Watertable Array ;;;;;;;;;;;;;;;;;;;;;;;;;

wtfind,fits,totgau,wt2d_pos,chin=chin
wtfind,fits,totgau,wt2d_neg,chin=chin-1024
wt2d=wt2d_pos<wt2d_neg

save,wt2d,wt2d_neg,wt2d_pos,chin,file=wtout
wt2d_neg=0
wt2d_pos=0
;;;;;;;;;;;;;;;;;;;;;;;;; Convolve It! ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

sz=size(fits)

cv=fits

for i=0,sz[3]-1 do begin
	cv[*,*,i]=convol(fits[*,*,i],d2posgau,/nan)
	fits[*,*,i]=convol(fits[*,*,i],d2neggau,/nan)
endfor

cv *= scale
fits *= scale

for i=0,sz[1]-1 do begin
	for j=0,sz[2]-1 do begin
		cv[i,j,*]=convol(cv[i,j,*],v2posgau,/nan)
		fits[i,j,*]=convol(fits[i,j,*],v2neggau,/nan)
	endfor
endfor

cv-=fits
fits=0b

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;; Divide by the watertable to get s/N ;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
for i=0,sz[3]-1 do begin
cv[*,*,i]=cv[*,*,i]/wt2d
endfor

;;;;;;;;;;;;;;; Fix NaN Values ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if bct GT 0 then cv[bad_ind]=0
bad_ind=0b

;;;;;;;;;;;;;;;;;;;;;;; Save Output ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;; Consider Removing this for speed ;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;; Maybe Make it a keyword option ;;;;;;;;;;;;;;;;;;;;;;

cv=temporary(cv[nx/2:sz[1]-nx/2-1,nx/2:sz[2]-nx/2-1,*])
sz=size(cv)
cv[0,*,*]=0
cv[sz[1]-1,*,*]=0
cv[*,0,*]=0
cv[*,sz[2]-1,*]=0


if keyword_set(noisescale) then cv*=noisescale else noisescale=1



if keyword_set(outfile) then begin
save,cv,noisescale,f=outfile

endif

end
