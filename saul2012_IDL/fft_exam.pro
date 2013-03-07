pro fft_exam, j, srng, img, medpwr

i1 = fltarr(512, 512)

for i=srng[0], srng[1] do begin
    img = readfits('~/Narrow/' + cname(j,3) + '_N.fits', hdr, nslice=i, /sil)
    i1 = i1+img
endfor

i2 = fltarr(512, 512)

for i=srng[0], srng[1] do begin
    img = readfits('~/Narrow/' + cname(j,4) + '_N.fits', hdr, nslice=i, /sil )
    i2 = i2+img
endfor

y0 = 100

img = fltarr(512, 512)

img[*,0:511-y0-16] = i1[*,y0:511-16]
img[*,512-y0-16:*] = i2[*, 16:y0+32-1]

whh = where(img lt -100, ct)

if ct ne 0 then img(whh) = 0.

pwr = (float(fft(img)*conj(fft(img))))[0:255, 0:255]

xx = rebin(reform(findgen(256), 256, 1), 256, 256)
yy = rebin(reform(findgen(256), 1, 256), 256, 256)
rr = sqrt(xx^2 + yy^2)
wh = where(rr le 255)
medpwr = fltarr(256)
pwrwh = pwr[wh]
rrwh=  rr[wh]

for i=0, 255 do medpwr[i] = median(pwrwh[where((rrwh ge i) and (rrwh lt (i+1)))])

end

pro fftloop, j, meds, amp, vs

nloop = 64
dv = 8
v0 = 1024-dv*nloop/2.

vs = (findgen(nloop)-nloop/2+0.5)*dv*0.184

meds = fltarr(256, nloop)
amp = fltarr(nloop)
for i=0, nloop-1 do begin
    loop_bar, i, nloop
    fft_exam, j, [v0+i*dv, v0+(i+1)*dv-1], img, medpwr
    meds[*, i] = medpwr
    amp[i] = mean(img)
endfor

end
