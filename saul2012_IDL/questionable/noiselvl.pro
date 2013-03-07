pro noiselvl

!p.multi=[0,2,3]
restore,'pig_sp_NOISE.sav',/verb
g=where(pig_sp.primary eq 0,ct_np)
print,'secondarys', ct_np
g=where(pig_sp.kern_num eq 0,ct)
pig_sp_0=pig_sp[g]
print,'kernel 0',ct
g=where(pig_sp.kern_num eq 1,ct)
pig_sp_1=pig_sp[g]
print,'kernel 1',ct
g=where(pig_sp.kern_num eq 2,ct)
pig_sp_2=pig_sp[g]
print,'kernel 2',ct
g=where(pig_sp.kern_num eq 3,ct)
pig_sp_3=pig_sp[g]
print,'kernel 3',ct

stop

hh=histogram(pig_sp_0.mxsig,binsize=0.1,min=4)
bins=findgen(n_elements(hh))*0.1+4+0.05
plot,bins,hh,psym=10,xtit='Max S/N in Convolved Cube',ytit='Number of Clouds',font=0,thick=3,xthick=3,ythick=3

hh=histogram(pig_sp_1.mxsig,binsize=0.1,min=4)
bins=findgen(n_elements(hh))*0.1+4+0.05
plot,bins,hh,psym=10,xtit='Max S/N in Convolved Cube',ytit='Number of Clouds',font=0,thick=3,xthick=3,ythick=3

hh=histogram(pig_sp_2.mxsig,binsize=0.1,min=4)
bins=findgen(n_elements(hh))*0.1+4+0.05
plot,bins,hh,psym=10,xtit='Max S/N in Convolved Cube',ytit='Number of Clouds',font=0,thick=3,xthick=3,ythick=3

hh=histogram(pig_sp_3.mxsig,binsize=0.1,min=4)
bins=findgen(n_elements(hh))*0.1+4+0.05
plot,bins,hh,psym=10,xtit='Max S/N in Convolved Cube',ytit='Number of Clouds',font=0,thick=3,xthick=3,ythick=3



stop

restore,'wt_124+26_noise_18.5.sav',/verb
extract_coords,'GALFA_HI_RA+DEC_124.00+26.35_W.fits',ra,dec,vels
sz=size(ra)
ra=ra[30:sz[1]-30-1,30:sz[2]-30-1]
dec=dec[30:sz[1]-30-1,30:sz[2]-30-1] 
tosave=where(dec[0,0:451] gt 24)
ra=ra[*,min(tosave):max(tosave)]
dec=dec[*,min(tosave):max(tosave)]
wt2d=wt2d[*,min(tosave):max(tosave)] 
vin=where(abs(vels) lt 600)
vels=vels[vin]
vran=minmax(vels)
raran=minmax(ra)
decran=minmax(dec)

minnn=min(wt2d[where(wt2d ne 0)])

display,wt2d,ra[*,0],dec[0,*],min=minnn,xtit='RA',ytit='Dec',font=0,thick=3
oplot,pig_sp.ra,pig_sp.dec,psym=4,color=255,thick=3


plot,bins,hh,psym=10,xtit='Max S/N in Convolved Cube',ytit='Number of Clouds',font=0,thick=3,xthick=3,ythick=3

plot,pig_sp.v,pig_sp.dec,xran=vran,yran=decran,/xsty,/ysty,psym=4,$
	xtit='Velocity (km/s)',ytit='Dec',font=0,thick=3,xthick=3,ythick=3

plot,pig_sp.v,pig_sp.ra,xran=vran,yran=raran,/xsty,/ysty,psym=4,$
	xtit='Velocity (km/s)',ytit='Dec',font=0,thick=3,xthick=3,ythick=3

plot,pig_sp.kern_num,psym=4.,/xsty,/ysty,ytit='Cloud Kernel',$
	xtit='Cloud Number',font=0,thick=3,xthick=3,ythick=3
plot,[0,100],[0,100],/nodata,font=0
xyouts,2,80,'Total Negative Primary Clouds'+string(n_elements(pig_sp)),font=0
xyouts,2,70,'Total Non-Primary Clouds'+string(ct_np),font=0



stop
end
