; fit a 3D gaussian plus piston, tip, tilt
; use data coordinates for simplicity

function mpfit_3DG, p, x=x, err=err, model=model
q = p

if not keyword_set(err) then err = 1.
if not keyword_set(err) then print, 'nerr'
sz = size(x)

xx = rebin(reform(findgen(sz[1]), sz[1], 1, 1), sz[1], sz[2], sz[3])
yy = rebin(reform(findgen(sz[2]), 1, sz[2], 1), sz[1], sz[2], sz[3])
vv = rebin(reform(findgen(sz[3]), 1, 1, sz[3]), sz[1], sz[2], sz[3])

; parameters: 
; [0] amplitude, [1] x-position, [2] y-position, [3] v-position, [4]
; x-sigma, [5] y-sigma, [6] v-sigma, [7] position angle, [8] piston,
; [9] slope with x, [10] slope with y, [11], slope with v
ellipr = fltarr(sz[1], sz[2], sz[3])
ellipr =  sqrt((( (cos(q[7]*!pi/180.)*(xx-q[1]) + sin(q[7]*!pi/180.)*(yy-q[2])))/(q[4]))^2+((( (-1)*sin(q[7]*!pi/180.)*(xx-q[1]) + cos(q[7]*!pi/180.)*(yy-q[2])))/(q[5]))^2 + ((vv-q[3])/q[6])^2)

model = fltarr(sz[1], sz[2], sz[3])
model = q[0]*exp((-0.5)*ellipr^2)+q[8]+ q[9]*(xx-sz[1]/2.) +  q[10]*(yy-sz[2]/2.) +  q[11]*(vv-sz[3]/2.)

return, reform((x-model), sz[1]*sz[2]*sz[3])/err
end
