;+
; what do you wrap pigs in? blankets, of course.
;
; NAME: BLANKETS_2D
;
; PURPOSE: Wrapper for Pigs Routines (to be renamed truffles in thefuture)
;
; INPUTS:
;      I:  RA  (degrees) ?
;      J:  Dec (degrees) ?
;
; KEYWORDS:
;      FAKE: Triggers Guinea Pig Injection
;
; OUTPUTS:
;      Writes Convolved Data to : '~/cv/cv(_inj)_i_j.sav'
;      Prints amount of time it took to run.
;
; CALLS:
;      CV_ALL
;      FIND_REG
;      FIND_PIGS
;      PIG_PROPS
;      ID_INJ_PIGS (optional)
;
; HISTORY:
;      thepast      : Written in by JEGP
;      Jan. 9, 2010 : Documented by Destry
;-

pro blankets_2d, i, j, dataf, fake=fake,_extra=_extra

tm = systime(/sec)

;;;;;;;;;;;; Convolve Cube with Matched Filter Function ;;;;;;;;;;;;;

multi_cv,i,j,dataf,cv=cv,ra=ra,dec=dec,/nosave

save, dec, ra, cv, f= '~/cv/cv_' + string(i, f='(I2.2)') +'_' + string(j, f='(I2.2)')+ '.sav'
print, (systime(/sec)-tm)/60., 'MULTI_CV runtime (minutes)'
tm=systime(/sec)
stop
;;;;;;;;;;;; Find Candidate Regions in Convolved Data ;;;;;;;;;;;;;;;

find_reg, i, j, fake=fake
print,(systime(/sec)-tm)/60.,'Find_reg runtime (minutes)'
stop
;;;;;;;;;;;; Identify Small Candidate Regions ;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;; Insert Those Regions into Structure ;;;;;;;;;;;;;;;;;;;;
tm=systime(/sec)
np = find_pigs(i, j, fake=fake,_extra=_extra)
print,(systime(/sec)-tm)/60.,'Find_pigs runtime (minutes)'
stop

    
    ;;;;;;;;;;;;;;; Find Pig Properties Using MPFIT ;;;;;;;;;;;;;;;;;

if np ne 0 then begin   
tm=systime(/sec)
    pig_props, i, j, fake=fake,_extra=_extra
;    plot_pig_prosps, i, j
print,(systime(/sec)-tm)/60.,'Pig_props runtime (minutes)'
endif


;;;;;;;;;;;; If Using Injected Pigs, Do Something ;;;;;;;;;;;;;;;;;;;

if keyword_set(fake) then id_inj_pigs, i, j


end



