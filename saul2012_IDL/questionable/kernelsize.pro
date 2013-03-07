pro kernelsize,fwhm=fwhm,v_fwhm=v_fwhm 
;;;;;;;;;;;;;;;;;;;; Set Initial Parameters ;;;;;;;;;;;;;;;;;;;;;;;;;

dv = 0.736122839     ; channel width km/s
wid=1.2                 ; Width of Negative Gaussian relative to Positive Gaussian

;;;;;;;;;;;;;;;;;;;;;;;;;; Set Scaling ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if ~keyword_set(fwhm) then fwhm=7.0
if ~keyword_set(v_fwhm) then v_fwhm=5.0
scale=1
if (fwhm EQ 7) and (v_fwhm EQ 5) then scale=0.120881
if (fwhm EQ 7) and (v_fwhm EQ 15) then scale=0.106187
if (fwhm EQ 18) and (v_fwhm EQ 5) then scale=0.116593
if (fwhm EQ 18) and (v_fwhm EQ 15) then scale=0.102864
if (fwhm NE 18) and (fwhm NE 7) then print,'INVALID Kernel Size, Scaling will be incorrect'
if (v_fwhm NE 5) and (v_fwhm NE 15) then print,'INVALID Kernel FWHM, Scaling will be incorrect'

;;;;;;;;;;;;;;;;;;;; Build Filtering Function ;;;;;;;;;;;;;;;;;;;;;;;

nx=41    ; Spatial Size of Convovling box
nnv=41   ; Velocity Size of Convolving box

xxx = rebin(reform(findgen(nx)-nx/2, nx, 1, 1), nx, nx, nnv)
yyy = rebin(reform(findgen(nx)-nx/2, 1, nx, 1), nx, nx, nnv)
vvv = rebin(reform((findgen(nnv)-nnv/2)*dv, 1, 1, nnv), nx, nx, nnv)

d2posgau = 2*exp( (-0.5)*(xxx[*,*,0]/(fwhm/(2*sqrt(2*alog(2)))))^2)*exp( (-0.5)*(yyy[*,*,0]/(fwhm/(2*sqrt(2*alog(2)))))^2)
v2posgau = 2*exp( (-0.5)*(vvv[0,0,*]/(v_fwhm/(2*sqrt(2*alog(2)))))^2)

d2neggau = exp( (-0.5)*(xxx[*,*,0]/(wid*fwhm/(2*sqrt(2*alog(2)))))^2)*exp( (-0.5)*(yyy[*,*,0]/(wid*fwhm/(2*sqrt(2*alog(2)))))^2)
v2neggau = exp( (-0.5)*(vvv[0,0,*]/(wid*v_fwhm/(2*sqrt(2*alog(2)))))^2)

d2neggau *= total(d2posgau)/total(d2neggau)
v2neggau *= total(v2posgau)/total(v2neggau)

dscale=max(d2posgau-d2neggau)
vscale=max(v2posgau-v2neggau)

d2posgau /=dscale
d2neggau /=dscale
v2posgau /=vscale
v2neggau /=vscale

;;;;;;;;;;;;;;;;; Build 3d kerel for wtfind.pro ;;;;;;;;;;;;;;;;;;;;;

posgau = 2*exp( (-0.5)*(xxx/(fwhm/(2*sqrt(2*alog(2)))))^2)*exp( (-0.5)*(yyy/(fwhm/(2*sqrt(2*alog(2)))))^2)*exp( (-0.5)*(vvv/(v_fwhm/(2*sqrt(2*alog(2)))))^2)

neggau=  exp( (-0.5)*(xxx/(wid*fwhm/(2*sqrt(2*alog(2)))))^2)*exp( (-0.5)*(yyy/(wid*fwhm/(2*sqrt(2*alog(2)))))^2)*exp( (-0.5)*(vvv/(wid*v_fwhm/(2*sqrt(2*alog(2)))))^2)

sccc=total(posgau)/total(neggau)

neggau = neggau*total(posgau)/total(neggau)

totgau = posgau-neggau
totgau = totgau/max(totgau)

print,'total = ',total(abs(totgau))
f=where(totgau gt 0.1,ct)
print,'number of voxels about 0.01 = ',ct
end
