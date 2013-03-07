pro find_reg,rmask=rmask,minpeak=minpeak,water=water,datacv=datacv,outfile=outfile

;This is step 2 in the pig-searching code suite
;
;The point of this code is to take a 3D unsharp masked data cube of
;the 'PIGS' region and find possible pig candidates. The plan is to
;mask the data, then find a bright peak, grow a mask for the bright
;peak until it reaches, say, 4x the noise in the map, and find out how
;large that region is. Then repeat with the next brightest position,
;outside that region, repeating until all regions above some amplitude
;are covered
; settibng fk to '_inj' uses the pig injected cubes from cv_all


; the new, convolved, unsharpmasked data set
if ~keyword_set(datacv) then begin
	print,'Set DATACV so I have some data to work with'
endif else begin
	restore,datacv,/verb
endelse

szcv = size(cv)
dv = szcv[3]
dx=szcv[1]
dy=szcv[2]

; make a copy to house the region data
rmask = byte(cv*0)

; preset values for where the coverage can spread to, and which peaks
; to select

if ~keyword_set(minpeak) then minpeak = 6.
if ~keyword_set(water) then cutoff=4. else cutoff=water
mx = max(cv,xx,/nan)
q = 2b

while mx gt minpeak do begin
    oldcount=0
    print, q
    qmask = long([xx])
    grow = 1
while grow eq 1 do begin

if n_elements(qmask) lt 3.e6 then begin
	h = Histogram([xx,qmask+1,qmask-1,qmask+dx,qmask-dx,qmask+dx*dy,qmask-dx*dy], OMIN=omin)
	qmask = Where(h GT 0) + omin
endif else begin
        h = Histogram([xx,qmask+1,qmask-1], OMIN=omin)
	tempqmask= Where(h GT 0)+ omin
        h = Histogram([tempqmask,qmask+dx*dy,qmask-dx*dy], OMIN=omin)
	tempmask = Where(h GT 0)+ omin
	h = Histogram([tempqmask,qmask+dx,qmask-dx], OMIN=omin)
        tempqmask=0b
	qmask = Where(h GT 0) + omin
	print,'used double hists'
endelse
	wh=where(cv[qmask] ge cutoff,count)
	qmask=qmask[wh]
	if count gt oldcount then begin
		oldcount=count
	endif else begin
		grow=0
	endelse
endwhile
	cv[qmask] = 0
	rmask[qmask]=q
	mx = max(cv,xx,/nan)
	q++
    if q EQ 255 then begin
	rmask=fix(rmask)
	q=fix(q)
    endif

endwhile

if keyword_set(outfile) then begin
save,rmask,f=outfile
endif


end
