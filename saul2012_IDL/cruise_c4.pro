pro cruise_c4, dataf=dataf, datacoords=datacoords,mask=mask, vrng=vrng, cont1=cont1, cont2=cont2,level=level, ra0=ra0, dec0=dec0, xdata=xdata, ydata=ydata, rmasks=rmasks,spfile=spfile,REF_EXTRA=_extra,cruisedata=cruisedata,cruiseout=cruiseout

if keyword_set(cruisedata) then begin
restore,cruisedata,/verb
endif else begin

restore,dataf,/verb
restore,datacoords,/verb
cube=temporary(fits)
;extract_coords,dataf,ra,dec,vels,cube

bad=where((finite(cube) EQ 0) or (cube LT -30),count)
if count GT 0 then cube[bad]=0

sz = size(cube)
if not keyword_set(mask) then mask = fltarr(sz[1], sz[2]) + 1
if not keyword_set(vrng) then vrng = findgen(sz[3])
if not keyword_set(ra0) then ra0 = findgen(sz[1])
if not keyword_set(dec0) then dec0 = findgen(sz[2])
restore,spfile,/verb
rmaskstr=replicate({rmask:intarr(sz[1],sz[2],sz[3])},4)
;rmaskstr_em=replicate({rmask:intarr(sz[1],sz[2],sz[3])},4)
for i=0,3 do begin
restore,rmasks[i],/verb
whkern=where(pig_sp.kern_num EQ i,np)
for j=0,np-1 do begin
whhh=where(rmask eq pig_sp[whkern[j]].rm_num)
if pig_sp[whkern[j]].embed eq 1 then begin
	if pig_sp[whkern[j]].primary eq 1 then begin
;		rmaskstr_em[i].rmask[whhh]=j+1001 
		rmaskstr[i].rmask[whhh]=j+1001
	endif else begin
;		rmaskstr_em[i].rmask[whhh]=j+1
		rmaskstr[i].rmask[whhh]=j+1
	endelse
endif else begin
	if pig_sp[whkern[j]].primary eq 1 then begin
;		rmaskstr[i].rmask[whhh]=j+1001
		rmaskstr[i].rmask[whhh]=j+3001
	endif else begin
;		rmaskstr[i].rmask[whhh]=j+1
		rmaskstr[i].rmask[whhh]=j+2001
	endelse
endelse
endfor
endfor

if keyword_set(cruiseout) then begin
;save,cube,ra,dec,vels,extendo,pig_sp,rmaskstr,rmaskstr_em,f=cruiseout
save,cube,ra,dec,vels,extendo,pig_sp,rmaskstr,f=cruiseout
endif
endelse

sz = size(cube)
if not keyword_set(mask) then mask = fltarr(sz[1], sz[2]) + 1
if not keyword_set(vrng) then vrng = findgen(sz[3])
if not keyword_set(ra0) then ra0 = findgen(sz[1])
if not keyword_set(dec0) then dec0 = findgen(sz[2])





device,dec=0

levels=[0.5,999,1998,1999,2999]
lines=[1,0,0,1,0]
thick=[1,1,1,3,3]
color1=[90,90,0,90,90]
color2=[40,40,0,40,40]
red=200
ctable=12
velit=0
kernit=0
stopit=0
infoit=0
userinput=''
v=sz[3]/2
kern=0
kern2=1
setovel=1


loadct,0,/sil
display, cube[*, *, sz[3]/2]*mask, ra0, dec0, aspect=1.
        loadct,ctable,/sil
	contour, rmaskstr[kern].rmask[*,*,v], level=levels, c_linesty=lines,/overplot, c_color=color1,c_thick=thick
        contour,rmaskstr[kern2].rmask[*,*,v], level=levels,c_linesty=lines, /overplot, c_color=color2,c_thick=thick
	contour, extendo[*,*,v], level=[0.5],/overplot, c_color=red
	loadct,0,/sil




while stopit EQ 0 do begin
read,userinput,prompt='What do you want to do? Options are: Vel, Kern, Info, or Stop: '
if userinput EQ 'Vel' then velit=1 else begin
if userinput EQ 'Kern' then kernit = 1 else begin
if userinput EQ 'Stop' then stopit =1 else begin
if userinput EQ 'Info' then infoit =1 else begin
print,'I dont understand, try again'
endelse
endelse
endelse
endelse

while infoit EQ 1 do begin
print,'Left Click on a location for information'
cursor,x,y,/down,/data
x=floor(x)
y=floor(y)
print,'x=',x,' y=',y,' v=',v
whkern=where(pig_sp.kern_num EQ kern,np)
pigindex=(rmaskstr[kern].rmask[x,y,v] mod 1000)-1
if pigindex GE 0 then begin
	print,'Blue Contour is pig #',whkern[pigindex]
help,pig_sp[whkern[pigindex]],/str
endif
whkern=where(pig_sp.kern_num EQ kern2,np)
pigindex=(rmaskstr[kern2].rmask[x,y,v] mod 1000)-1
if pigindex GE 0 then begin
	print,'Green Contour is pig #',whkern[pigindex]
help,pig_sp[whkern[pigindex]],/str
endif

infoit=0
endwhile




while velit EQ 1 do begin
print,'Move mouse over image. Drag Left Mouse button to change displayed Velocity. Right Click to exit.'
!mouse.button=0.
while (!mouse.button eq 0) or (!mouse.button eq 1) do begin
    cursor, x, y, /change, /norm
if !mouse.button eq 1 then  begin
	if setovel eq 1 then begin
		ovel=v
		oxel=x
		setovel=0
	endif
	v = ((ovel + ((x-oxel)*sz[3])/4) < sz[3]) > 0.
endif else begin
	setovel=1
endelse
    display, cube[*, *, v]*mask, ra0, dec0, aspect=1.

        loadct,ctable,/sil
        contour, rmaskstr[kern].rmask[*,*,v], level=levels, c_linesty=lines,/overplot, c_color=color1,c_thick=thick
        contour,rmaskstr[kern2].rmask[*,*,v], level=levels,c_linesty=lines, /overplot, c_color=color2,c_thick=thick
        contour, extendo[*,*,v], level=[0.5],/overplot, c_color=red
        loadct,0,/sil





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
        loadct,ctable,/sil
        contour, rmaskstr[kern].rmask[*,*,v], level=levels, c_linesty=lines,/overplot, c_color=color1,c_thick=thick
        contour,rmaskstr[kern2].rmask[*,*,v], level=levels,c_linesty=lines, /overplot, c_color=color2,c_thick=thick
        contour, extendo[*,*,v], level=[0.5],/overplot, c_color=red
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
        loadct,ctable,/sil
        contour, rmaskstr[kern].rmask[*,*,v], level=levels, c_linesty=lines,/overplot, c_color=color1,c_thick=thick
        contour,rmaskstr[kern2].rmask[*,*,v], level=levels,c_linesty=lines, /overplot, c_color=color2,c_thick=thick
        contour, extendo[*,*,v], level=[0.5],/overplot, c_color=red
        loadct,0,/sil





endwhile

endwhile


endwhile

end
