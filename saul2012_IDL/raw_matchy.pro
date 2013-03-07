pro raw_matchy

; investigating raw regions matching


restore,'rawregions.save',/verb
restore,'REAL_15_03.sav',/verb

found_ind=intarr(n_elements(g_inj))
nmatches=intarr(n_elements(ras[0,*]))

ras[*,2]=0
decs[*,2]=0

for i=0,n_elements(g_inj)-1 do begin

inds=where(((g_inj[i].ra LT ras[1,*]) and (g_inj[i].ra GT ras[0,*]) and (g_inj[i].dec LT decs[1,*]) and (g_inj[i].dec GT decs[0,*]) and (g_inj[i].v LT vs[1,*]) and (g_inj[i].v GT vs[0,*])),ct)

found_ind[i]=ct
if ct GT 0 then nmatches[inds]+=1

endfor
found=g_inj[where(found_ind GT 0)]
stop

end

