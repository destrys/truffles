function wax_off,cv,verb=verb

sz=size(cv)
coinsize=60
help,/mem
help
wax=cv*0
for i=0,sz[3]-1 do begin
	if keyword_set(verb) then print,'slice ',i
	wax[*,*,i]=median(cv[*,*,i],coinsize)
left=rebin(wax[coinsize/2,*,i],coinsize/2,sz[2])
right=rebin(wax[sz[1]-coinsize/2-1,*,i],coinsize/2,sz[2])
wax[0:coinsize/2-1,*,i]=left
wax[sz[1]-coinsize/2:sz[1]-1,*,i]=right
bottom=rebin(wax[*,coinsize/2,i],sz[1],coinsize/2)
top=rebin(wax[*,sz[2]-coinsize/2-1,i],sz[1],coinsize/2)
wax[*,0:coinsize/2-1,i]=bottom
wax[*,sz[2]-coinsize/2:sz[2]-1,i]=top
endfor
wax=cv-temporary(wax)

return,wax

end
