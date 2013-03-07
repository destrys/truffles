pro effhist,data1,data2,binsize=binsize,fillcolor=fillcolor,oplot=oplot,xrange=xrange,_extra=extra

;;; Histogram Data ;;;;;

if ~keyword_set(binsize) then binsize=0

hist1data=histogram(data1,binsize=binsize,_extra=extra,omin=omin,omax=omax)
if ~keyword_set(binsize) then binsize=(omax-omin)/n_elements(hist1data)
hist2data=histogram(data2,binsize=binsize,nbins=n_elements(hist1data),min=omin,_extra=_extra)
zerr=where(hist1data EQ 0,ct)
if ct GT 0 then hist1data[zerr] =1
fhistdata=float(hist2data)/float(hist1data)



IF KeyWord_Set(xrange) then xind=xrange else xind=[omin,omax]
IF ~KeyWord_Set(oplot) Then Begin
	plot,xind,[0,max(fhistdata)*1.1],/nodata,_extra=extra
ENDIF

start = omin
endpt = start + binsize

FOR j=0,N_Elements(fhistdata)-1 DO BEGIN
        x = [start, start, endpt, endpt, start]
        y = [0, fhistdata[j], fhistdata[j], 0, 0]
	if keyword_set(fillcolor) then begin
		polyfill,x,y,color=fillcolor,/data
	ENDIF
	OPLOT, x, y, _extra=extra
        start = start + binsize
        endpt = start + binsize
ENDFOR
stop
end

