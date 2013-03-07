pro progress

;progress = replicate({exp:0,spore:0,cv:0,},44,5)

progress = dblarr(45,5)

;loadct,13
;plot,[0,.25],[0,.25],psym=2,xrange=[0,45],yrange=[0,5],/xstyle

;for i=0,43 do begin
;for j=0,4 do begin

;print,'/hpc/30days/astro/users/jmg2223/Wide/GALFA_HI_RA+DEC_' + cname2(i,j) + '_W_EXP.fits'
;progress[i,j] = file_test('/hpc/30days/astro/users/jmg2223/Wide/GALFA_HI_RA+DEC_' + cname2(i,j) + '_W_EXP.fits') + file_test('/hpc/30days/astro/users/jmg2223/Run1/spores_' + cname2(i,j) + '.sav') + file_test('/hpc/30days/astro/users/jmg2223/Run1_cv/cv_' + cname2(i,j) + '_7.5.sav') + file_test('/hpc/30days/astro/users/jmg2223/Run1_cv/wt_' + cname2(i,j) + '_7.5.sav')+ file_test('/hpc/30days/astro/users/jmg2223/Run1_cv/rmask_'+ cname2(i,j) + '_7.5.sav')

;oplot,i,j,psym=2,color=((progress[i,j]/4)*255)

;endfor
;endfor
;print,progress
;set_plot,'ps'
;device,filename='/hpc/30days/astro/users/jmg2223/progress.ps'
;device,DECOMPOSED=0
;image = congrid(progress,880,100)
;imSize=SIZE(image)
;WINDOW,XSIZE=imSize[1],YSIZE=imSize[2]
;tvscl,image
;device,/close
;[where(progress eq 0)]
;,psym=2,xrange=[0,45],yrange=[0,5],/xstyle,/ystyle,color=10
;tv,progress[where(progress eq 1)],psym=2,color=50
;oplot,progress[where(progress eq 2)],psym=2,color=100
;oplot,progress[where(progress eq 3)],psym=2,color=150
;oplot,progress[where(progress eq 4)],psym-2,color=200
openw,2,'t_log.dat'

for i=0,44 do begin
for j=0,4 do begin
;i=2
;j=5

m=string(i, format='(I2.2)')
n=string(j, format='(I2.2)')
outstring = cname2(i,j)
openw,lun,'t_' + m + '_' + n + '.submit',/get_lun
printf,lun,("universe = vanilla")
printf,lun,("notification = Always")
;printf,lun,("notify_user = janagrc@gmail.com")
printf,lun,("+AccountingGroup = " + """group_astro""")
printf,lun,("executable= /usr/local/bin/idl")
printf,lun,("environment =" + """IDL_STARTUP=/hpc/scratch/astro/users/jmg2223/gsr2.6/start_ao.idl.pro GALFAPATH=/hpc/scratch/astro/users/jmg2223/gsr2.6/ CARLPATH=/hpc/scratch/astro/users/jmg2223/gsr2.6/carlpath/ GSRPATH=/hpc/scratch/astro/users/jmg2223/gsr2.6/ TREETOP=/hpc/scratch/astro/users/jmg2223/ PHILTOP=/hpc/scratch/astro/users/jmg2223/ PHILPATH=/share/galfa/gsr2.6/phil/""")
printf,lun,("initialdir = /hpc/30days/astro/users/drs2125/checkrun/")
printf,lun,("input=" + 't_' + m + '_' + n + ".idl")
printf,lun,("output =t_"+ m + '_' + n + ".out")
printf,lun,("error =t_"+ m + '_' + n + ".err")
printf,lun,("log =t_" + m + '_' + n + ".log")
printf,lun,("requirements = BigMem == TRUE")
printf,lun,("queue")

free_lun,lun

openw,lun,'t_'+ m + '_' + n + '.idl'
printf,lun,'i=',i
printf,lun,'j=',j
printf,lun,'.com mpfit_3DG'
printf,lun,'.com plot_pig_prosps_3D'
printf,lun,("filein = '/hpc/30days/astro/users/jmg2223/Wide/GALFA_HI_RA+DEC_" + outstring + "_W_EXP.fits'")
printf,lun,("truffles_check,datacube=filein,outstring='_" + outstring + "',chin=[1250,1800]")
free_lun,lun

;openw,lun,'t_'+ m + '_' + n + '.idl',/get_lun
;printf,lun,'t_'+ m + '_' + n
;free_lun,lun

printf,2,("condor_submit " + 't_' + m + '_' + n + ".submit")

endfor
endfor

end
