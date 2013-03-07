pro find_raw

restore,'~/cv/rmask_15_03.sav'
extract_coords,'~/truffles/cubes/GALFA_HI_RA+DEC_124.00+26.35_W_SPAM_0.fits',ra,dec,vels,fits

n_c=max(rmask)
ras=fltarr(2,n_c)
decs=ras
vs=ras

ramask=rebin(ra,512,512,2048)
decmask=rebin(dec,512,512,2048)
vmask=rebin(reform(vels,1,1,2048),512,512,2048)


for i = 2,max(rmask) do begin


;Find range
interest=where(rmask eq i)
    ras[*,i-1] = minmax(ramask[interest])
    decs[*,i-1] = minmax(decmask[interest])
    vs[*,i-1] = minmax(vmask[interest])

endfor

stop
end
