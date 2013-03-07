; Go look for all the rmask regions that are small in ra/dec/vel space
; as save them into a structure called pig_sp
function find_pigs_4, datacoords=datacoords,spout=spout,dataf=dataf,cvs=cvs,rmasks=rmasks,keepedges=keepedges
;+
;NAME: 
;     FIND_PIGS
;PURPOSE:
;     Look through the rmask output of FIND_REG and select the 
;     regions that are 'small'. 
;     These Regions are saved into a structure called pig_sp
;INPUT:
;     m: Cube Index #
;     n: Cube Index #
;KEYWORDS:
;     /FAKE: Indicates that there are INnjected Guinea Pigs in the data.
;     DATAF: Filename of Data Cube if different from CNAME format.
;     RA   : RA Array
;     Dec  : Dec Array
;     Fits : Data Cube
;     Vels : Velocity Array
;     Rmask: Masked Region Cube
;     pig_sp: Clouds structure.
;HISTORY:
;     thepast:  Written by JEGP
;     Jan. 12, 2010: Documented (Destry)
;                  : Changed 'small' constraints
;     Jan. 14, 2010: Inserted DATAF keyword. (Destry)
;     Jan. 22, 2010: Added Output Keywords pig_sp
;                    Added Input Keywords ra,dec,fits,vels,rmask (Destry)  
;-


;;;;; Set Sizes for Extended Regions ;;;;
dx = 38
dy = 38
dv = 200
voff = 100    ; number of channels on edges to ignore


restore,rmasks[0],/verb

;;; opening ra/dec/vels  ;;;;;
restore,datacoords
sz=size(ra)
if ~keyword_set(keepedges) then begin
ra=ra[30:sz[1]-30-1,30:sz[2]-30-1]
dec=dec[30:sz[1]-30-1,30:sz[2]-30-1]
endif
sz=size(ra)
nx=41
ra=ra[nx/2:sz[1]-nx/2-1,nx/2:sz[2]-nx/2-1]
dec=dec[nx/2:sz[1]-nx/2-1,nx/2:sz[2]-nx/2-1]

sz = size(rmask)
extendo=byte(rmask*0)
q=0
npig=0
for i=0,3 do begin
if (i eq 0) or (i eq 2) then begin
dx=36
dy=36
endif else begin
dx=38
dy=38
endelse

if i NE 0 then restore,rmasks[i],/verb
restore,cvs[i],/verb

mx = max(rmask)
amiapig = fltarr(mx+1)
nelcl = n_elements(where(rmask gt 1))
if n_elements(whpigs) eq 0 then begin
   whpigs = lonarr(nelcl)
   pignum = whpigs
   pignumr = whpigs
   pigsig = fltarr(nelcl)
endif else begin
   whpigs = [whpigs,lonarr(nelcl)]
   pignum = [pignum, lonarr(nelcl)]
   pignumr = [pignumr, lonarr(nelcl)]
   pigsig = [pigsig,fltarr(nelcl)]
endelse

for j=2,mx do begin
    if j/10 eq j/10. then print, string(j) + '/' +string(fix(mx))

    wwrmask= where(rmask eq j)
        vinds=wwrmask/(sz[1]*sz[2])
        yinds=(wwrmask-(vinds*sz[1]*sz[2]))/sz[1]
        xinds=(wwrmask-(vinds*sz[1]*sz[2])) mod sz[1]
        mmx=minmax(xinds)
        mmy=minmax(yinds)
        mmv=minmax(vinds)
    szx = mmx[1] - mmx[0]
    szy = mmy[1] - mmy[0]
    szv = mmv[1] - mmv[0]

    ; are we next to an edge?
    ; FIX THIS TO DETECT edge of observed region

    if (szx le dx) and (szy le dy) and (szv le dv) then amiapig[j] = 1.
    if (mmx[0] le 1) or (mmx[1] ge (sz[1]-2)) or (mmy[0] le 1) or (mmy[1] ge (sz[2]-2)) then amiapig[j] = 0


;if the object is too 'big'
if amiapig[j] eq 0 then begin
; add to the region in which things are extended
extendo[wwrmask] += byte(2^i)
endif

;record positions and numbers of pigs
if amiapig[j] EQ 1 then begin
   nelpig = n_elements(wwrmask)
   whpigs[q:q+nelpig-1] = wwrmask
   pignum[q:q+nelpig-1] = npig
   pignumr[q:q+nelpig-1] = j
   pigsig[q:q+nelpig-1] = cv[wwrmask]
   
   q = q+nelpig
   npig++

endif

endfor

whpigs = whpigs[0:q-1]
pignum = pignum[0:q-1]
pignumr = pignumr[0:q-1]
pigsig = pigsig[0:q-1]
if i eq 0 then kerpig = byte(pignum*0) else kerpig = [kerpig, bytarr(n_elements(pignum)-n_elements(kerpig))+i]
endfor
;plot 'em
delv = 0.736122839

pig_sp = replicate({ra:0., dec:0., v:0., sp:fltarr(sz[3]), vels:fltarr(sz[3]), rm_num:0., spm:fltarr(sz[3]), osp:fltarr(sz[3]),embed:0,primary:0,ghostof:0L,kern_num:0,mxsig:0.,px:intarr(2),py:intarr(2),pv:intarr(2)}, npig)

;populate
pig_sp[pignum].kern_num=kerpig

; label embedded pigs as such.
embedwh = where(extendo[whpigs] gt 0, ct)
if ct ne 0 then begin
      pig_sp[pignum[embedwh]].embed =1
endif
pig_sp.primary = -1
; figure out primacy
mxsig = fltarr(npig)
for f = 0, npig-1 do begin
   region=where(pignum eq f)
   mxsig[f] = max(pigsig[region],place)
   if (n_elements(place) gt 1) then place=place[0]
   spot=whpigs[region[place]]
   placev=floor(spot/(sz[1]*sz[2]))
   placey=(floor(spot-sz[1]*sz[2]*placev))/sz[1]
   placex=floor(spot-sz[1]*sz[2]*placev) mod sz[1]
   pig_sp[f].ra=ra[placex,placey]
   pig_sp[f].dec=dec[placex,placey]
   pig_sp[f].v=vels[placev]
   pig_sp[f].mxsig=mxsig[f]
   wwrmask=whpigs[region]
        vinds=wwrmask/(sz[1]*sz[2])
        yinds=(wwrmask-(vinds*sz[1]*sz[2]))/sz[1]
        xinds=(wwrmask-(vinds*sz[1]*sz[2])) mod sz[1]
   pig_sp[f].px=minmax(xinds)
   pig_sp[f].py=minmax(yinds)
   pig_sp[f].pv=minmax(vinds)
endfor

primcube = byte(extendo*0.)
whundefarea = bytarr(q)+1b
while min(pig_sp.primary) eq (-1) do begin
   print, n_elements(where(pig_sp.primary eq -1)), n_elements(pig_sp)
; where do we not know if a clouds is primary?
   mx = max(mxsig*(pig_sp.primary EQ -1), pcloud)
   ; the primary cloud
   pig_sp[pcloud].primary =1
   pig_sp[pcloud].ghostof =-1
   ; the area in 3-space where the primary cloud lives
   pc_area = whpigs[where(pignum eq pcloud)]
 
;;;;; attempting without primcube;;;;;;
;;; need to compare pc_area to whpigs
;pc_min=min(pc_area)
;pc_max=max(pc_area)
;hh=where((histogram(pc_area,min=pc_min,max=pc_max) ne 0) and (histogram(whpigs,min=pc_min,max=pc_max,rev=revrev) ne 0),whct)
;if whct gt 0 then begin
;h_ghosts=uniq(sort(pignum[revrev[hh]]

;endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  whundefarea[where(pignum EQ pcloud)] =0b
   primcube[pc_area] = 1b
   whg = where(primcube[whpigs]*whundefarea eq 1b, ct)
   if ct ne 0 then begin
      ghosts = pignum[whg]
      pig_sp[ghosts].ghostof = pcloud
      pig_sp[ghosts].primary = 0
      ughost = ghosts[uniq(ghosts, sort(ghosts))]
      for i=0, n_elements(ughost)-1 do begin
         whundefarea[where(pignum eq ughost[i])] =0b
      endfor
   endif
   primcube[pc_area] = 0b

endwhile


smbound = 20.

rmask=0b
cv=0b
primcube=0b

;;; restoring original data cube ;;;
restore,dataf
sz=size(fits)
if ~keyword_set(keepedges) then begin
fits=temporary(fits[30:sz[1]-30-1,30:sz[2]-30-1,*])
endif
nx=41
sz=size(fits)
fits=temporary(fits[nx/2:sz[1]-nx/2-1,nx/2:sz[2]-nx/2-1,*])
sz=size(fits)

for i=0, 3 do begin

pigsr_i = pignumr[where(kerpig eq i)]
whpigs_i= whpigs[where(kerpig eq i)]
pigs_i=pignum[where(kerpig eq i)]
pgi = pigs_i(uniq(pigs_i, sort(pigs_i)))
mx = max(pgi)
amiapig = fltarr(mx+1)
amiapig[pgi] = 1.
ra=float(ra)
dec=float(dec)
for jjj=0, n_elements(pgi)-1 do begin
j=pgi[jjj]
jr=max(pignumr[where(pignum eq j)])
whpigs_i_r=whpigs[where(pignum eq j)]
        loc=where((ra eq pig_sp[j].ra) and (dec eq pig_sp[j].dec))
        yc=loc/sz[1]
	xc=loc mod sz[1]
        xmx = (xc+smbound) < (sz[1] -1)
        xmn = (xc-smbound) > 0
        ymx = (yc+smbound) < (sz[2] -1)
        ymn = (yc-smbound) > 0
        vinds=whpigs_i_r/(sz[1]*sz[2])
        yinds=((whpigs_i_r-(vinds*sz[1]*sz[2]))/sz[1])-ymn[0]
        xinds=((whpigs_i_r-(vinds*sz[1]*sz[2])) mod sz[1])-xmn[0]
        msk = fltarr(xmx-xmn+1, ymx-ymn+1, sz[3])
        szm = size(msk)
	whpigs_i_r=vinds*szm[1]*szm[2]+yinds*szm[1]+xinds
        msk[whpigs_i_r] = 1.
        imgmask = total(msk, 3) < 1
        holemask = gromask(gromask(gromask(imgmask)))
        offmask = gromask(gromask(gromask(holemask))) - holemask
        wh = where((msk[*, *, sz[3]/2]) eq 0, ct)
        if ct ne 0 then offmask[wh] = 0.
        spmask = total(total(msk, 1), 1) < 1
        spect = total(total(rebin(reform(imgmask, szm[1], szm[2], 1), szm[1], szm[2], szm[3])*fits[xmn:xmx,ymn:ymx, *], 1), 1)/total(imgmask)
        offspect = total(total(rebin(reform(offmask, szm[1], szm[2], 1), szm[1], szm[2], szm[3])*fits[xmn:xmx,ymn:ymx, *], 1), 1)/total(offmask)
;        plot, vels, spect, title=string(j)
;        oplot, vels(where(spmask eq 1)), spect(where(spmask eq 1)), color=100
        pig_sp[j].sp = spect
        pig_sp[j].vels = vels
        pig_sp[j].spm = spmask
        pig_sp[j].osp = offspect
        pig_sp[j].rm_num = jr
;help,/mem
endfor

endfor

save, pig_sp,extendo, f=spout

return, npig

end
