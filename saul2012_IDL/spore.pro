pro spore,fits,chin=chin,ra=ra,dec=dec,vels=vels,dataf=dataf,outfile=outfile,nosmooth=nosmooth,keepedges=keepedges

;+
; NAME: SPORE
;
;
;
;
;
;-

;;;;; Load File ;;;;;;;

if keyword_set(dataf) then begin
	restore,dataf,/verb
	fits=float(fits)
endif

;;;;;;;;;;;;;;; Fix Negative Values ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

bad_ind=where(fits LT -30,ct)
if ct GT 0 then fits[bad_ind]=!values.f_nan
bad_ind=0b

;;;;;; Stamp Smooth ;;;;;;
help,/mem
help
if ~keyword_set(nosmooth) then begin
	fits=wax_off(fits,/verb)
endif


if ~keyword_set(chin) then begin

;;;;; Select Watertable channels ;;;;;
;; Get rid of nans
na=where(finite(fits) EQ 0,ct)
scrap=fits
if ct GT 0 then scrap[na]=0

print, 'choose a range over which to find the typical noise'
print, 'first channel'
cruise, scrap, ch1
wait, 0.5
print, 'last channel'
cruise, scrap, ch2

scrap=0

; if they chose backwards, switch 'em
if ch1 gt ch2 then begin
   chs = ch1
   ch1 = ch2
   ch2 = chs
endif

chin=[ch1,ch2]
endif


;;;;;;; Unless Told otherwise, Slice off edges ;;;;;;

sz=size(fits)
if ~keyword_set(keepedges) then begin
fits=temporary(fits[30:sz[1]-30-1,30:sz[2]-30-1,*])
endif



;;;;;; Optional Save ;;;;;;;


if keyword_set(outfile) then save,fits,chin,f=outfile

end
