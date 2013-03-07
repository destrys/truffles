function matchy,spam,pig_sp,savef=savef
;+
;
;
;-

sz=size(spam)
n_spam=sz[1]

pmatch=intarr(n_elements(pig_sp))
smatch=intarr(n_spam)
p_vels=pig_sp.vels*pig_sp.spm
p_vels[where(p_vels EQ 0)]=!values.f_nan
vmin=min(p_vels,dim=1,/nan)
vmax=max(p_vels,dim=1,/nan)
for i=0,n_spam-1 do begin
space_dist=sqrt((pig_sp.ra-spam[i].ra)^2+(pig_sp.dec-spam[i].dec)^2)*60.
scrap=(space_dist LT spam[i].sz) and (vmax GT spam[i].v) and (vmin LT spam[i].v)

smatch[i]=total(scrap)
pmatch+=scrap
endfor

found=spam[where(smatch GT 0)]

stop
return,found
end
