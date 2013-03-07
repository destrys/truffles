readcol, 'Coordinate_Data.txt', name, num, ra1, ra2, ra3, dec1, dec2, dec3, v, f=['A, F, F, F,F, F, F, F, F']

nick = replicate({name:'', num:0., ra:0., dec:0, v:0}, 91)

nick.num = num
nick.name=name
nick.v = v
for i=0, 90 do nick[i].ra = ten([ra1[i], ra2[i], ra3[i]])*15.
for i=0, 90 do nick[i].dec = ten([dec1[i], dec2[i], dec3[i]])

restore, 'pf_all.sav'
restore, 'glitch.sav', /ver
pf = pf_all[where(glitch eq 0)]

md = fltarr(91)

glactc, nick.ra/15, nick.dec, 2000, l, b, 1
whpf = fltarr(91)

for i=0, 90 do begin
    md[i] = min(sqrt((60*(nick[i].ra - pf.ra))^2 + (60*(nick[i].dec - pf.dec))^2 + (1.*(nick[i].v - pf.gfit[1]))^2), xx)
    whpf[i] = xx
endfor

wh = where((b gt 45) and (abs(nick.v) lt 45))

end
