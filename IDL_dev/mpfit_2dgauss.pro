; fit 2 3D gaussians
; use data coordinates for simplicity

function mpfit_2dgauss, p, x=x, err=err, model=model

if not keyword_set(err) then err = 1.
sz = size(x)

xx = rebin(reform(findgen(sz[1])-(sz[1]-1.)/2., sz[1], 1, 1), sz[1], sz[2])
yy = rebin(reform(findgen(sz[2])-(sz[2]-1.)/2., 1, sz[2], 1), sz[1], sz[2])

; parameters: 
; [0] amplitude, [1] x-position, [2] y-position, [3] x-sigma, [4] y-sigma,  [5] position angle, [6] offset, [7] dt/dx, [8] dt/dy
ellipr = fltarr(sz[1], sz[2])
ellipr =  sqrt((( (cos(p[5]*!pi/180.)*(xx-p[1]) + sin(p[5]*!pi/180.)*(yy-p[2])))/(p[3]))^2+((( (-1)*sin(p[5]*!pi/180.)*(xx-p[1]) + cos(p[5]*!pi/180.)*(yy-p[2])))/(p[4]))^2)


model = p[0]*exp((-0.5)*ellipr^2) +  p[6] + p[7]*xx + p[8]*yy



;p = reform(p, 8*2)
return, reform((x-model), sz[1]*sz[2])/err
end
