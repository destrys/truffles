pro region_wrapper

restore,'~/cv_15_03_1_5.fits',/verb

cv=wax_off(cv)

tm=systime(/sec)
find_reg, i, j, fake=fake,_extra=_extra,cv=cv,rmask=rmask,fits=fits
print,(systime(/sec)-tm)/60.,'Find_reg runtime (minutes)'


restore,'~/cv_15_03_1_15.fits',/verb



stop
end

