pro convert_wh, whold, xo, yo, zo, xp, yp, zp, xn, yn, zn, whnew 
;make integers
xo = floor(xo)
yo = floor(yo)
zo = floor(zo)
xn = floor(xn)
yn = floor(yn)
zn = floor(zn)
xp = floor(xp)
yp = floor(yp)
zp = floor(zp)

whnx = whold mod xo
whny = whold/xo mod yo
whnz = whold/(xo*yo)

whnew = (whnx + xp) + (whny + yp)*xn + (whnz+zp)*(xn*yn)
end
