pro addacloud

fits=readfits('GALFA_HI_RA+DEC_132.00+26.35_W_big2.fits',hdr,/noscale)

dv = 0.736122839
nx= 41
nv = 41
xxx = rebin(reform(dindgen(nx)-nx/2, nx, 1, 1), nx, nx, nv)
yyy = rebin(reform(dindgen(nx)-nx/2, 1, nx, 1), nx, nx, nv)
vvv = rebin(reform((dindgen(nv)-nv/2)*dv, 1, 1, nv), nx, nx, nv)




fwhm=20.
v_fwhm=10.
posgau = 2*exp( (-0.5)*(xxx/(fwhm/(2*sqrt(2*alog(2)))))^2)*exp( (-0.5)*(yyy/(fwhm/(2*sqrt(2*alog(2)))))^2)*exp( (-0.5)*(vvv/(v_fwhm/(2*sqrt(2*alog(2)))))^2)

posgau/=max(posgau)
posgau/=0.0066227
posgau=fix(posgau)
fits[619-70-20:619-70+20,400-20:400+20,800-nv/2:800+nv/2]+=posgau

writefits,'GALFA_HI_RA+DEC_132.00+26.35_W_big2_plus.fits',fits,hdr


fits=readfits('GALFA_HI_RA+DEC_124.00+26.35_W_big2.fits',hdr,/noscale)

fits[69-20:69+20,400-20:400+20,800-nv/2:800+nv/2]+=posgau


writefits,'GALFA_HI_RA+DEC_124.00+26.35_W_big2_plus.fits',fits,hdr


end
