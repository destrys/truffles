pro find_reg_newer,rmask=rmask,minpeak=minpeak,water=water,datacv=datacv,outfile=outfile

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
	crap=where(finite(cv) eq 0,ct)
if ct gt 0 then cv[crap]=0
endelse

szcv = size(cv)
dv = szcv[3]
dx=szcv[1]
dy=szcv[2]

cv[*, *, 0:15] = 0
cv[*, *, dv-16:*] = 0


; make a copy to house the region data
rmask = byte(cv*0)

; preset values for where the coverage can spread to, and which peaks
; to select

if ~keyword_set(minpeak) then minpeak = 5.
if ~keyword_set(water) then cutoff=3. else cutoff=water
mx = max(cv,xx,/nan)
q = 2b


while mx gt minpeak do begin
oldcount=0
	qloop=0
    print, q
    qmask = ulong([xx])
    grow = 1
tm=systime(/seco)/60.
    while grow eq 1 do begin
;print,qloop
;qloop++
newq=qmask+1
   h = Histogram([newq, qmask], OMIN=omin)
   qmask = Where(h GT 0) + omin

    ;mina = Min(newq, Max=maxa)
    ;uniqq = Where((Histogram(newq, Min=mina, Max=maxa) NE 0) AND $
    ;         ( Histogram(qmask, Min=mina, Max=maxa) EQ 0), count)
;if count GT 0 then newinds=newq[uniqq]

newq=qmask-1
   h = Histogram([newq, qmask], OMIN=omin)
   qmask = Where(h GT 0) + omin
;    mina = Min(newq, Max=maxa)
;    uniqq = Where((Histogram(newq, Min=mina, Max=maxa) NE 0) AND $
;             ( Histogram(qmask, Min=mina, Max=maxa) EQ 0), count)
;if count GT 0 then begin


newq=qmask+dx
   h = Histogram([newq, qmask], OMIN=omin)
   qmask = Where(h GT 0) + omin
;    mina = Min(newq, Max=maxa)
;    uniqq = Where((Histogram(newq, Min=mina, Max=maxa) NE 0) AND $
;             ( Histogram(qmask, Min=mina, Max=maxa) EQ 0), count)
;if count GT 0 then newinds=[newinds,newq[uniqq]]

newq=qmask-dx
   h = Histogram([newq, qmask], OMIN=omin)
   qmask = Where(h GT 0) + omin
;    mina = Min(newq, Max=maxa)
;    uniqq = Where((Histogram(newq, Min=mina, Max=maxa) NE 0) AND $
;             ( Histogram(qmask, Min=mina, Max=maxa) EQ 0), count)
;if count GT 0 then newinds=[newinds,newq[uniqq]]

newq=qmask+dx*dy
   h = Histogram([newq, qmask], OMIN=omin)
   qmask = Where(h GT 0) + omin
;    mina = Min(newq, Max=maxa)
;    uniqq = Where((Histogram(newq, Min=mina, Max=maxa) NE 0) AND $
;             ( Histogram(qmask, Min=mina, Max=maxa) EQ 0), count)
;if count GT 0 then newinds=[newinds,newq[uniqq]]

newq=qmask-dx*dy
   h = Histogram([newq, qmask], OMIN=omin)
   qmask = Where(h GT 0) + omin
;    mina = Min(newq, Max=maxa)
;    uniqq = Where((Histogram(newq, Min=mina, Max=maxa) NE 0) AND $
;             ( Histogram(qmask, Min=mina, Max=maxa) EQ 0), count)
;if count GT 0 then newinds=[newinds,newq[uniqq]]

;newinds=newinds[uniq(newinds,sort(newinds))]

wh=where(cv[qmask] ge cutoff,count)
if count gt oldcount then begin
qmask=qmask[wh]
oldcount=count
endif else begin
grow=0
endelse

     endwhile
print,systime(/seco)/60.-tm,'Time for loop'
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
