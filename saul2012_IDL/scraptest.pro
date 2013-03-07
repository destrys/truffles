pro scraptest,config=config

if keyword_set(config) then begin
plot,[3,4]
!mouse.button=0
print,'Position your mouse over the display window and click your left mouse button'
cursor,x,y,/up
button1=!mouse.button
!mouse.button=0
print,'Position your mouse over the display window and click your right mouse button'
cursor,x,y,/up
button2=!mouse.button
endif





!mouse.button=0.
while !mouse.button ne button2 do begin
while !mouse.button eq 0 do begin
	cursor, x, y, /change, /norm
print,'nobutton'
wait,0.2
endwhile
while !mouse.button eq button1 do begin
    cursor, x, y, /change, /norm
print,'button1'
wait,0.2
endwhile
endwhile

end
