pro cruise_c3, cube, v, mask=mask, vrng=vrng, cont1=cont1, cont2=cont2,level=level, ra0=ra0, dec0=dec0, xdata=xdata, ydata=ydata, ct1=ct1,ct2=ct2, rmasks=rmasks,REF_EXTRA=_extra

bad=where((finite(cube) EQ 0) or (cube LT -30),count)
if count GT 0 then cube[bad]=0

sz = size(cube)
if not keyword_set(mask) then mask = fltarr(sz[1], sz[2]) + 1
if not keyword_set(vrng) then vrng = findgen(sz[3])
if not keyword_set(ra0) then ra0 = findgen(sz[1])
if not keyword_set(dec0) then dec0 = findgen(sz[2])

rmaskstr=replicate({rmask:bytarr(sz[1],sz[2],sz[3])},4)
for i=0,3 do begin
restore,rmasks[i],/verb
rmaskstr[i].rmask=rmask
endfor
device,dec=0

velit=0
kernit=0
stopit=0
userinput=''
v=sz[3]/2
kern=0
kern2=1

loadct,0,/sil
display, cube[*, *, sz[3]/2]*mask, ra0, dec0, aspect=1.
      cimg = reform(rmaskstr[kern].rmask[*, *, v])
        cimg2=reform(rmaskstr[kern2].rmask[*,*,v])
       levels = cimg(uniq(cimg, sort(cimg))) - 0.5
        levels=[2]
        loadct,1,/sil
       contour, cimg, level=levels, /overplot, c_color=220
        loadct,6,/sil
        contour,cimg2, level=levels, /overplot, c_color=120
        loadct,0,/sil





while stopit EQ 0 do begin
read,userinput,prompt='What do you want to do? Options are: Vel, Kern, or Stop: '
if userinput EQ 'Vel' then velit=1 else begin
if userinput EQ 'Kern' then kernit = 1 else begin
if userinput EQ 'Stop' then stopit =1 else begin
print,'I dont understand, try again'
endelse
endelse
endelse




while velit EQ 1 do begin
print,'Move mouse over image. Drag Left Mouse button to change displayed Velocity. Right Click to exit.'
!mouse.button=0.
while (!mouse.button eq 0) or (!mouse.button eq 1) do begin
    cursor, x, y, /change, /norm
if !mouse.button eq 1 then  v = (x*sz[3] < sz[3]) > 0.
    display, cube[*, *, v]*mask, ra0, dec0, aspect=1.
    if keyword_set(cont1) then begin
       cimg = reform(rmaskstr[kern].rmask[*, *, v])
	cimg2= reform(rmaskstr[kern2].rmask[*,*,v])
       levels = cimg(uniq(cimg, sort(cimg))) - 0.5
       levels=[2]
       loadct, 1, /sil
       contour, cimg, level=levels, /overplot, c_color=220
	loadct,6,/sil
	contour,cimg2, level=levels, /overplot, c_color=120
       loadct, 0, /sil
     endif
;if (!mouse.button eq 1) then g=(g+1) mod 4
;;    if keyword_set(cont2) then begin
;       cimg = reform(cont2[*, *, v])
;       levels = cimg(uniq(cimg, sort(cimg))) - 0.4
;       if keyword_set(ct2) then loadct, ct2, /sil
;       contour, cimg, level=levels, /overplot, c_color=findgen(n_elements(levels))/n_elements(levels)*255. 
;       if keyword_set(ct2) then loadct, 0, /sil
;    endif


    if keyword_set(xdata) then oplot, xdata, ydata, _EXTRA=_extra
    xyouts, 0.1, 0.1, vrng(v), /norm, charsize=1.5, charthick=1.5
    xyouts, 0.1,0.15, 'Kernel'+string(kern),/norm,charsize=1.5,charthick=1.5
endwhile
velit=0
endwhile


while kernit GT 0 do begin
read,kernwhich,prompt='Which Kernel Contour do you want to change? 0 for Blue, 1 for Green, -1 to Exit: '
if kernwhich EQ -1 then kernit=0
while kernwhich eq 0 do begin
read,kern3,prompt='Which Kernel would you like to overplot as the blue contour(0-3) Type -1 to stop: '
if (kern3 GE 0) and (kern3 LT 4) then kern=kern3 else begin
if kern3 eq -1 then kernwhich=-2 else begin
print,'I do not understand, try again.'
endelse
endelse
    display, cube[*, *, v]*mask, ra0, dec0, aspect=1.
      cimg = reform(rmaskstr[kern].rmask[*, *, v])
	cimg2=reform(rmaskstr[kern2].rmask[*,*,v])
       levels = cimg(uniq(cimg, sort(cimg))) - 0.5
        levels=[2]
	loadct,1,/sil
       contour, cimg, level=levels, /overplot, c_color=220
	loadct,6,/sil
	contour,cimg2, level=levels, /overplot, c_color=120
	loadct,0,/sil
endwhile
while kernwhich eq 1 do begin
read,kern3,prompt='Which Kernel would you like to overplot as the green contour(0-3) Type -1 to stop: '
if (kern3 GE 0) and (kern3 LT 4) then kern2=kern3 else begin
if kern3 eq -1 then kernwhich=-2 else begin
print,'I do not understand, try again.'
endelse
endelse
    display, cube[*, *, v]*mask, ra0, dec0, aspect=1.
      cimg = reform(rmaskstr[kern].rmask[*, *, v])
	cimg2=reform(rmaskstr[kern2].rmask[*,*,v])
       levels = cimg(uniq(cimg, sort(cimg))) - 0.5
        levels=[2]
	loadct,1,/sil
       contour, cimg, level=levels, /overplot, c_color=220
	loadct,6,/sil
	contour,cimg2, level=levels, /overplot, c_color=120
	loadct,0,/sil

endwhile

endwhile


endwhile

end
