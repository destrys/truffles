function wax_off,cv,coinsize,verbose=verbose
;+
; FUNCTION: WAX_OFF
;
; PURPOSE: Apply median smooth to a 3D datacube
;
; SYNTAX: cube = wax_off(cube,coinsize[,/verbose])
;
; CALLED BY: SPORE
;
; NOTES:
;   For a 512 x 512 x ? array, each slice takes
;   ~7 seconds on a Macbook Pro 2.66 GHz i7
;   For a GALFA cube (512x512x2048) that's 4hr.
;
; HISTORY:
;    Created by Destry in 2010
;    Cleaned and sped up - 3.11.2013 - Destry
;-

;Init
sz=size(cv)
halfCoin=coinsize/2


if keyword_set(verbose) then begin
print,'WAX_OFF Starting'
print,systime()
lasttm=systime(/sec)
help,/mem
help
endif

for i=0,sz[3]-1 do begin
    if keyword_set(verbose) then begin
        print,'slice ',i
        tm=systime(/sec)
        print,'lastloop '+strtrim(tm-lasttm,1)+' seconds'
        lasttm=tm
    endif

    cv[*,*,i]-=median(cv[*,*,i],coinsize)

    cv[0:halfCoin-1,*,i]=rebin(cv[halfCoin,*,i],halfcoin,sz[2])
    cv[sz[1]-halfCoin:sz[1]-1,*,i]=rebin(cv[sz[1]-halfCoin-1,*,i],halfCoin,sz[2])

    cv[*,0:halfCoin-1,i]=rebin(cv[*,halfCoin,i],sz[1],coinsize/2)
    cv[*,sz[2]-halfCoin:sz[2]-1,i]=rebin(cv[*,sz[2]-halfCoin-1,i],sz[1],coinsize/2)
endfor

return,cv

end
