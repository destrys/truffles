pro cv_comb

;
; CV_COMB
;      combines cv cubes of different convoled sizes
;
;

; Scaling Factors
fifteen=19.2923
ten=8.6476
five=2.16263

; Buidl cubes
best_cv=fltarr(512,512,2048)

restore,'~/cv/cv_15_03_0_3.4.sav'

best_cv=cv


restore,'~/cv_15_03_0_5.sav'

best_cv=best_cv>cv

restore,'~/cv_15_03_0_10.sav'

best_cv=best_cv>cv


restore,'~/cv_15_03_0_15.sav'


cv=best_cv>cv

stop

end


