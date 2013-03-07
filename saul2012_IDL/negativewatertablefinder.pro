pro negativewatertablefinder


restore,'cv_124+26_spam_0_18.15.sav',/verb

cv=cv[100:300,100:300,1600:1800]
;reg1=cv[100:250,100:400,200:600]
;reg2=cv[100:250,100:400,1448:1848]
hh=histogram(cv,binsize=0.1,min=-40)
bins=findgen(n_elements(hh))*0.1-40+0.05
plot,bins,hh,psym=10,/ylog,xran=[-40,40],yran=[1,10.^8],xtit='convolved voxel value',$
	ytit='Number of Voxels',tit='18arcmin 15km/s 124+26 spam0'
;hh=histogram(reg2,binsize=0.1,min=-40)
;bins=findgen(n_elements(hh))*0.1-40+0.05
;oplot,bins,hh,psym=10

print,stddev(cv)

fitbins=bins[where(abs(bins) lt 10 and abs(bins) gt .3)]
fithh=hh[where(abs(bins) lt 10 and abs(bins) gt .3)]

scrap=gaussfit(fitbins,fithh,ghgh,chisq=chi,nterms=3)

oplot,fitbins,scrap,thick=2

print,ghgh
;print,stddev(reg2)
stop

end
