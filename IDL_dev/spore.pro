pro spore,fits,coinsize=coinsize,outfile=outfile,keepedges=keepedges,_ref_extra=ex

;+
; NAME: SPORE
;
; PURPOSE: Median smooths the 1st and 2nd dimensions of a 
;          datacube.
;
; SYNTAX: spore,fits,[coinsize=60,outfile='filename.sav',/keepedges,/verbose]
;
; INPUTS: fits - datacube
;         coinsize - size (in pixels) of square median smooth.
;                    default = 60
;         outfile - if present, saves smoothed cube to 'filename.sav'
;         /keepedges - if present, doesnot trim edges, otherwise
;                      removes half the smoothing area size from each
;                      edge of the cube
;         /verbose - if present, prints extra info
;
; CALLS: WAX_OFF
;
; CALLED BY: TRUFFLES
;
; NOTES: Wants NANs where there is no data.
;   For a 512 x 512 x ? array, each slice takes
;   ~7 seconds on a Macbook Pro 2.66 GHz i7
;   For a GALFA cube (512x512x2048) that's 4hr. 
;
; HISTORY: Written in 2010 by Destry
;          Simplified and sped up - 3.11.13 - Destry
;-

; Initialize
if ~keyword_set(coinSize) then coinSize=60
halfCoin=fix(coinSize)/2

if keyword_set(verbose) then begin
print,'BEGINNING SPORES'
help,/mem
help
endif

;;;;;; 2D Median Filter ;;;;;;
fits=wax_off(fits,coinsize,_extra='verbose')

;;;;;;; Unless Told otherwise, Slice off edges ;;;;;;
sz=size(fits)
if ~keyword_set(keepedges) then begin
fits=temporary(fits[halfCoin:sz[1]-halfCoin-1,halfCoin:sz[2]-halfcoin-1,*])
endif

;;;; Update Header? ;;;;

;;;;;; Optional Save ;;;;;;;

if keyword_set(outfile) then save,fits,f=outfile

end
