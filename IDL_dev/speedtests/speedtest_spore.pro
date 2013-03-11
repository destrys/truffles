pro speedtest_spore
;+
; FUNCTION: SPEEDTEST_SPORE
;
; PURPOSE: Test the speed of the spore routine.
;
; HISTORY: Written on March 11 by Destry.
;-

;Initialize array
; GALFA-sized cube: 512 x 512 x 2048

testcube = randomn(123,512,512,100)

tm=systime(/sec)
spore,testcube,coinsize=59,/verbose
print, (systime(/sec)-tm), 'SPORE runtime'

end
