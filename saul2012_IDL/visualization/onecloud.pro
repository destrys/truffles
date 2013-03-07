pro onecloud,suffix=suffix

restore,'names'+suffix+'.sav'

restore,'dataf'+suffix+'.sav'
restore,'datacoords'+suffix+'.sav'
restore,spout
sz=size(fits)
stopit=0
kern=0
setovel=1
red=200
ctable=12
color1=90
device,dec=0

;restore,cvs[kern]
restore,rmasks[kern]
;fits=temporary(cv)
while stopit EQ 0 do begin
read,userinput,prompt='Which Cloud do you want to look at? (-1 to Stop): '

if userinput gt 0 then begin
cloud=userinput

if pig_sp[cloud].kern_num ne kern then begin
kern=pig_sp[cloud].kern_num
;restore,cvs[kern]
restore,rmasks[kern]
;fits=temporary(cv)
endif
sz=size(fits)
velit=1

px0=(pig_sp[cloud].px[0]-10) > 0
px1=(pig_sp[cloud].px[1]+10) < sz[1]-1
py0=(pig_sp[cloud].py[0]-10) > 0
py1=(pig_sp[cloud].py[1]+10) < sz[2]-1
pv0=(pig_sp[cloud].pv[0]-10) > 0
pv1=(pig_sp[cloud].pv[1]+10) < sz[3]-1

cube=fits[px0:px1,py0:py1,pv0:pv1]
rcube=(rmask[px0:px1,py0:py1,pv0:pv1] eq pig_sp[cloud].rm_num)*1
nrcube=((rmask[px0:px1,py0:py1,pv0:pv1] ne pig_sp[cloud].rm_num) and (rmask[px0:px1,py0:py1,pv0:pv1] ne 0))*1

excube=(extendo[px0:px1,py0:py1,pv0:pv1] gt 0)
ra0=ra[px0:px1,py0:py1]
dec0=dec[px0:px1,py0:py1]

sz=size(cube)
v=sz[3]/2
;ra0=ra0[0:sz[1]-1]
;dec0=reform(dec0[0,0:sz[2]-1])
ra0=findgen(sz[1])
dec0=findgen(sz[2])

while velit EQ 1 do begin
print,'Move mouse over image. Drag Left Mouse button to change displayed Velocity. Right Click to exit.'
!mouse.button=0.
d_min=min(cube[*, *, v])
d_max=max(cube[*,*,v])
while (!mouse.button eq 0) or (!mouse.button eq 1) do begin
    cursor, x, y, /change, /norm
if !mouse.button eq 1 then  begin
        if setovel eq 1 then begin
                ovel=v
                oxel=x
                setovel=0
        endif
        v = ((ovel + ((x-oxel)*sz[3])) < (sz[3]-1)) > 0.
endif else begin
        setovel=1
endelse
	loadct,0,/sil

    display, cube[*, *, v], ra0, dec0,min=d_min,max=d_max, aspect=1.
        loadct,ctable,/sil
       contour, nrcube[*,*,v], level=[1],/overplot, c_color=red,c_thick=3
        contour, rcube[*,*,v], level=[1],/overplot, c_color=color1,c_thick=3
        contour, excube[*,*,v], level=[1],/overplot, c_color=red
        loadct,0,/sil






endwhile
velit=0
endwhile

endif else stopit=1
endwhile

end
