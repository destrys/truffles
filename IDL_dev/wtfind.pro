pro wtfind, data, kernel, wt2d,chin=chin,chout=chout,outfile=outfile

szd = size(data)
szk = size(kernel)
; first determine a range over which to fit for noise

if keyword_set(chin) then begin
ch1=chin[0]
ch2=chin[1]
endif else begin
print, 'choose a range over which to find the typical noise'
print, 'first channel'
cruise, data, ch1
wait, 0.5
print, 'last channel'
cruise, data, ch2
endelse

; if they chose backwards, switch 'em
if ch1 gt ch2 then begin
   chs = ch1
   ch1 = ch2
   ch2 = chs
endif

chout=[ch1,ch2]

; find the typical channel-to-channel noise per pixel
noise = fltarr(szd[1], szd[2])
for i=0, szd[1]-1 do begin
   for j=0, szd[2]-1 do begin
      noise[i, j] = stddev(data[i, j, ch1:ch2])
   endfor
endfor

; find the effect of the beam in 3d
kernel2D = total(kernel^2, 3)

; find the water table level, i.e. the 1-sigma level to divide the
; convolved cube by
wt2d=sqrt(convol(noise^2., kernel2d, /nan))

if keyword_set(outfile) then begin
save,chout,noise,wt2d,f=outfile
endif

end
