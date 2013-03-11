pro spore,fits,coinsize=coinsize,outfile=outfile,keepedges=keepedges,verbose=verbose

;+
; NAME: SPORE
;
; PURPOSE: Median smooths the 1st and 2nd dimensions of a 
;          datacube.
;
; SYNTAX: spore,fits,[coinsize=60,outfile='filename.sav',/keepedges,/verbose]
;
; INPUTS: fits - datacube
;         coinsize - size (in pixels) of square median smooth
;         outfile - 
;         /keepedges - 
;         /verbose - 
; CALLS: WAX_OFF
;
; CALLED BY: TRUFFLES
;
; NOTES: Wants NANs where there is no data.
;
; HISTORY: Written in 2010 by Destry
;          Simplified and sped up - 3.11.13 - Destry
;-

; Initialize
if ~keyword_set(coinSize) then coinSize=60
halfCoin=fix(coinSize)/2
coinSize=halfcoin*2

if keyword_set(verbose) then begin
print,'BEGINNING SPORES'
help,/mem
help
endif

;;;;;; Stamp Smooth ;;;;;;
fits=wax_off(fits,coinsize,/verb)

;;;;;;; Unless Told otherwise, Slice off edges ;;;;;;

sz=size(fits)
if ~keyword_set(keepedges) then begin
fits=temporary(fits[halfCoin:sz[1]-halfCoin-1,halfCoin:sz[2]-halfcoin-1,*])
endif

;;;;;; Optional Save ;;;;;;;

if keyword_set(outfile) then save,fits,f=outfile

end
