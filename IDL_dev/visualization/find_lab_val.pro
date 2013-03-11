pro find_lab_val, ll, bb, v, lab, ora, odec, ov, klab

  l = (reform(ll, 721, 361))[*, 0]
  b = reform((reform(bb, 721, 361))[0,*])
  glactc, ora/15., odec, 2000, ol, ob, 1
  olw = ((ol +180) mod 360)-180  
  klab = interpolate(lab, interpol(findgen(721), l, olw), interpol(findgen(361), b, ob),  interpol(findgen(891), v, ov))
  
end
