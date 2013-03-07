restore, '~goldston/cv/pf.sav'

psopen, 'pigstats'
!p.multi=[0, 2, 2]

gs =1
if gs eq 1 then begin
x0= 151.60294    
x1=   250.75624   
y0=    23.616753   
y1=    32.074200
lds_vel =  findgen(151)*1.0305-75*1.0305
ra = x0 + rebin(reform(findgen(x1-x0), x1-x0, 1), x1-x0, y1-y0)
dec = y0 + rebin(reform(findgen(y1-y0), 1, y1-y0), x1-x0, y1-y0)   
ldspec = ldsregion(ra, dec,  getenv('GSRPATH') +'savfiles/lds.sav')         
endif

plot, ((abs(pf.twodfit[3]) < 20.) > 0) > ((abs(pf.twodfit[4]) < 20.) > 0)*(2*sqrt(2*alog(2))), (pf.gfit[2]*2*sqrt(2*alog(2)) < 5) > 0, psym=1, xtitle='size [arcmin, FWHM]', ytitle='width [km/s FWHM]'                     
plots, ((abs(pf.twodfit[3]) < 20.) > 0) > ((abs(pf.twodfit[4]) < 20.) > 0)*(2*sqrt(2*alog(2))), (pf.gfit[2]*2*sqrt(2*alog(2)) < 5) > 0, psym=1, color=(-150)*pf.kp + 150

plot, (pf.gfit[1] < 50) > (-50), (pf.gfit[2]*2*sqrt(2*alog(2)) < 5) > 0, psym=1, xtitle='velocity [km/s]', ytitle='width [km/s FWHM]', xra=[-50, 50]
plots, (pf.gfit[1] < 50) > (-50), (pf.gfit[2]*2*sqrt(2*alog(2)) < 5) > 0, psym=1, color=(-150)*pf.kp + 150

plot, ((abs(pf.twodfit[3]) < 20) > 0) > ((abs(pf.twodfit[4]) < 20) > 0)*(2*sqrt(2*alog(2))), ((pf.gfit[0] >  0) <  100), psym=1, xtitle='size [arcmin, FWHM]', ytitle='Amplitude [km/s FWHM]'                     
plots, ((abs(pf.twodfit[3]) < 20.) > 0) > ((abs(pf.twodfit[4]) < 20..l) > 0)*(2*sqrt(2*alog(2))), ((pf.gfit[0] >  0) <  100), psym=1, color=(-150)*pf.kp + 150

plot, lds_vel, ldspec, title='average spectrum', xtitle='v_LSR, [km/s]', ytitle='T_B [K]', xra=[-60, 60]

psclose

end

