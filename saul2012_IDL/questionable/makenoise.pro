pro makenoise

fits=readfits('GALFA_HI_RA+DEC_124.00+26.35_W.fits',hdr,/noscale)

dv = 0.736122839
nx= 41
nv = 41
sz=size(fits)

fits=randomn(seed,sz[1],sz[2],sz[3])*.1
fits[0:sz[1]/2,*,*]*=5
fits/=0.0066227
fits=fix(fits)

writefits,'GALFA_HI_RA+DEC_NOISE_2.fits',fits,hdr


end
