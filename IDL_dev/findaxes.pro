pro findaxes, ra, dec, PA, majax, minax
; makes the assumption that the area is small
; ra and dec are in degrees
; ditching the cos dec for simplicity for now...
x = ra;*cos(mean(dec)*!pi/180.)
y = dec
nel = n_elements(ra)
p = fltarr(2, nel)
p[0, *] = x
p[1, *] = y
thetas = findgen(180)*!pi/180.
height = fltarr(180)
width = height
for i=0, 179 do begin
   RMAT = [[cos(thetas[i]),(-1)*sin(thetas[i])],[sin(thetas[i]),cos(thetas[i])]]
   p1 = RMAT#p
   width[i] = max(p1[0, *]) -min(p1[0, *])
   height[i] = max(p1[1, *]) -min(p1[1, *])
endfor
majax = (max(width, xx))*60. ; arcmin
minax = (height(xx))*60. ; arcmin
PA = 180-(findgen(180))[xx] ; how should this really be defined?

end
