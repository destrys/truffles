;+
; what do you wrap pigs in? blankets, of course.
;
; NAME: BLANKETS
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

pro blankets, i, j, fake=fake,_extra=_extra

tm = systime(/sec)

;;;;;;;;;;;; Convolve Cube with Matched Filter Function ;;;;;;;;;;;;;

cv_all, i, j, fake=fake,_extra=_extra,ra=ra,dec=dec,vels=vels,cv=cv,fits=fits

print, (systime(/sec)-tm)/60., 'CV_ALL runtime (minutes)'
tm=systime(/sec)
stop
;;;;;;;;;;;; Find Candidate Regions in Convolved Data ;;;;;;;;;;;;;;;

find_reg, i, j, fake=fake,_extra=_extra,cv=cv,rmask=rmask,fits=fits
print,(systime(/sec)-tm)/60.,'Find_reg runtime (minutes)'
stop
;;;;;;;;;;;; Identify Small Candidate Regions ;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;; Insert Those Regions into Structure ;;;;;;;;;;;;;;;;;;;;
tm=systime(/sec)
np = find_pigs(i, j, fake=fake,ra=ra,dec=dec,vels=vels,fits=fits,rmask=rmask,_extra=_extra)
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



