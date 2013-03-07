pro makenoiseconv

fits=readfits('GALFA_HI_RA+DEC_NOISE_2.fits',hdr,/noscale)

dv = 0.736122839
nx= 41
nv = 41
sz=size(fits)

fwhm=3.4
nx=41    ; Spatial Size of Convovling box
nnv=41   ; Velocity Size of Convolving box

xxx = rebin(reform(findgen(nx)-nx/2, nx, 1, 1), nx, nx, nnv)
yyy = rebin(reform(findgen(nx)-nx/2, 1, nx, 1), nx, nx, nnv)

d2posgau = 2*exp( (-0.5)*(xxx[*,*,0]/(fwhm/(2*sqrt(2*alog(2)))))^2)*exp( (-0.5)*(yyy[*,*,0]/(fwhm/(2*sqrt(2*alog(2)))))^2)


for i=0,sz[3]-1 do begin
        fits[*,*,i]=convol(fits[*,*,i],d2posgau,/nan)
endfor

stop

writefits,'GALFA_HI_RA+DEC_NOISE_3.4.fits',fits,hdr


end
