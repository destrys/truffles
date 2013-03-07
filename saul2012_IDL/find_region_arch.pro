pro find_region_arch, ras, decs, projs
zero=0
if n_elements(ras) eq 0 then begin
    zero = 1
    ras = [0]
    decs = [0]
    projs = ['']
endif
; on a plot in degrees ra, degrees, dec
cursor, x, y
dx = 4.
SP = {RA0:[(x/15.-dx/60./15) > 0.,(x/15.+dx/60./15) < 24], DEC0:[y-dx/60., y+dx/60.]}
RF=['PROJNUM']

rt = searcharchive(SP, RF)

if tag_exist(rt, 'projnum') then begin
    names = rt.projnum
    uname = names(uniq(names, sort(names)))
    nun = fltarr(n_elements(uname))
    for i=0, n_elements(uname)-1 do begin
        nun[i] = n_elements(where(names eq uname[i]))
    endfor
    mx = max(nun, xx)
    projs = [projs, uname[xx]]
    ras = [ras, x]
    decs = [decs, y]
endif

if zero then begin
    ras=ras[1:*]
    decs = decs[1:*]
    projs = projs[1:*]
endif

end
