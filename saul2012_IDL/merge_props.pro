fks = ['', '_inj']

for k = 0, 0 do begin
    if k eq 0 then begin
        xw = [replicate(17, 3), replicate(18, 5), replicate(19, 5), replicate(20, 5),  replicate(21, 5), replicate(22, 5), replicate(23, 4),replicate(24, 4), replicate(25, 4),  replicate(26, 4),  replicate(27, 3),  replicate(28, 4),  replicate(29, 5),  replicate(30, 2)] 
          yw = [findgen(3)+2, findgen(5), findgen(5), findgen(5),  findgen(5), findgen(5), [0, 1, 3, 4] ,[0,1,3,4], [0,1,3,4],  [0,1,3,4], [0,1,3], [0,1,3,4],  findgen(5),  [3, 4]]
      endif
      if k eq 1 then restore, 'docubes.sav'
    
    fk = fks[k]
    q=0
    qi = 0
    
    for f=0, n_elements(xw)-1 do begin
;    for n=3, 4 do begin
;        for m=19, 28 do begin
        m = xw[f]
        n = yw[f]
                                ;       if (not( (m ge 27) and (n eq 4)) or (fk eq '_inj')) then begin 
        restore, '~/cv/pig'+fk+'_fit_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'
        psz = (size(pig_fit.postagestamp))[1]
        q = q+ n_elements(pig_fit)
        if k eq 1 then  begin
            restore, '~/cv/pig_id_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'
            qi = qi + (size(id))[1]
            
        endif
                                ;      endif
                                ;  endfor
                                ; endfor
    endfor

    pf_all = replicate({gfit:fltarr(6), sigs:fltarr(6), twodfit:fltarr(9), ra:0., dec:0., postagestamp:fltarr(psz, psz)},q)
    if k eq 1 then begin
        idall = fltarr(qi, 4)
        gi_all = replicate({v:0., ra:0., dec:0., sz:0., as:0., a:0., pa:0., fwhm:0.}, qi)
    endif
    p = 0
    pi = 0 
for f=0, n_elements(xw)-1 do begin
;    for n=3, 4 do begin
;        for m=19, 28 do begin
;            if (not( (m ge 27) and (n eq 4)) or (fk eq '_inj')) then begin 
    m = xw[f]
    n = yw[f]
    
    
    restore, '~/cv/pig'+fk+'_fit_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'
    pf_all[p:p+n_elements(pig_fit)-1] = pig_fit
    p = p+ n_elements(pig_fit)
    if k eq 1 then begin
        restore, '~/cv/pig_id_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'
        idall[pi:pi+(size(id))[1]-1, *] = id
        restore, '~/cv/cv_inj_' + string(m, f='(I2.2)') +'_' + string(n, f='(I2.2)')+ '.sav'
        gi_all[pi:pi+(size(id))[1]-1] = g_inj
        pi = pi + (size(id))[1]
    endif
    ;        endif
;endfor
;endfor
endfor

if k eq 0 then save, pf_all, f='~/cv/pf_all.sav' else save, pf_all,idall,gi_all, f='~/cv/pf_all_inj.sav'

endfor
end
