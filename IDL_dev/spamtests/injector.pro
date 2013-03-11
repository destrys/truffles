; add guinea pigs over a range of velocities, positions, widths,
; aspect ratios, amplitudes, position angles, velocity widths

pro injector, data, vels, ras, decs, np, v_rng, ra_rng, dec_rng, sz_rng, as_rng, a_rng, pa_rng, fwhm_rng, g_inj, useginj=useginj

if not (keyword_set(useginj)) then g_inj = replicate({v:0., ra:0., dec:0., sz:0., as:0., a:0., pa:0., fwhm:0.}, np)
for i=0, np-1 do begin
    if not (keyword_set(useginj)) then begin
        g_inj[i].v = randomu(seed)*(max(v_rng) -min(v_rng)) + min(v_rng)
        g_inj[i].ra = randomu(seed)*(max(ra_rng) -min(ra_rng)) + min(ra_rng)
        g_inj[i].dec = randomu(seed)*(max(dec_rng) -min(dec_rng)) + min(dec_rng)
        g_inj[i].sz = randomu(seed)*(max(sz_rng) -min(sz_rng)) + min(sz_rng)
        g_inj[i].as = randomu(seed)*(max(as_rng) -min(as_rng)) + min(as_rng)
        g_inj[i].a = randomu(seed)*(max(a_rng) -min(a_rng)) + min(a_rng)
        g_inj[i].pa = randomu(seed)*(max(pa_rng) -min(pa_rng)) + min(pa_rng)
        g_inj[i].fwhm = randomu(seed)*(max(fwhm_rng) -min(fwhm_rng)) + min(fwhm_rng)
    endif

   ; positions in the data cube
    x = interpol(findgen(512), ras[*, 0], g_inj[i].ra)
    y = interpol(findgen(512), decs[0, *], g_inj[i].dec)
    v = interpol(findgen(n_elements(vels)), vels, g_inj[i].v)
    cube = fltarr( 1-floor((x-g_inj[i].sz*6) > 0) +floor((x+g_inj[i].sz*6) < 511), 1-floor((y-g_inj[i].sz*6) > 0) +floor((y+g_inj[i].sz*6) < 511),1-floor((v-g_inj[i].fwhm*6) > 0)+floor((v+g_inj[i].fwhm*6) < (n_elements(vels)-1)))
    szc = size(cube)

    rac = ras[(x-g_inj[i].sz*6) > 0:(x+g_inj[i].sz*6) < 511, (y-g_inj[i].sz*6) > 0:(y+g_inj[i].sz*6) < 511]
    decc = decs[(x-g_inj[i].sz*6) > 0:(x+g_inj[i].sz*6) < 511, (y-g_inj[i].sz*6) > 0:(y+g_inj[i].sz*6) < 511]
    vc = vels[(v-g_inj[i].fwhm*6) > 0: (v+g_inj[i].fwhm*6) < (n_elements(vels)-1)]
    

    p0 = [g_inj[i].a, (g_inj[i].ra - mean(rac))*60*(-1),(g_inj[i].dec - mean(decc))*60., g_inj[i].sz*g_inj[i].as/(2*sqrt(2*alog(2))),g_inj[i].sz/(2*sqrt(2*alog(2))), g_inj[i].pa, 0, 0, 0]
    pr = mpfit_2dgauss(p0, x=total(cube, 3), model=model)
    vgauss = exp((-1)*((vc-g_inj[i].v)/(g_inj[i].fwhm/(2*sqrt(2*alog(2)))))^2)
    for j=0, n_elements(vc)-1 do cube[*, *, j] = model*vgauss[j] 

    data[(x-g_inj[i].sz*6) > 0:(x+g_inj[i].sz*6) < 511, (y-g_inj[i].sz*6) > 0:(y+g_inj[i].sz*6) < 511,(v-g_inj[i].fwhm*6) > 0: (v+g_inj[i].fwhm*6) < (n_elements(vels)-1)] = data[(x-g_inj[i].sz*6) > 0:(x+g_inj[i].sz*6) < 511, (y-g_inj[i].sz*6) > 0:(y+g_inj[i].sz*6) < 511,(v-g_inj[i].fwhm*6) > 0: (v+g_inj[i].fwhm*6) < (n_elements(vels)-1)] + cube
endfor



end
