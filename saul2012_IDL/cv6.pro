pro cv4

;ra=['172.00','180.00','188.00','196.00','204.00','212.00','220.00']
;dec=['02.35','10.35','18.35','26.35','34.35']

ra =['044.00']
dec=['10.35']
;for i=0,6 do begin
;  for i=0,6 do begin
;   for j=0,2 do begin
i=0
j=0
    m=ra[i]
    n=dec[j]

    filein='/hpc/30days/astro/users/jmg2223/Wide/GALFA_HI_RA+DEC_'+ m + '+' + n + '_W_EXP.fits'
    ;fileout='GALFA_HI_RA+DEC_'+ m + '+' + n + '_W_EXP'

    truffles,datacube=filein,outstring='_' + m + '+' + n,chin=[1250,1800]
    print,m,n,'Complete'
;  endfor
;endfor


print,'Finished'

end
