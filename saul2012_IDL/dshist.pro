pro dshist,data,binsize=binsize,fillcolor=fillcolor,oplot=oplot,xrange=xrange,_extra=extra

;;; Histogram Data ;;;;;
if ~keyword_set(binsize) then binsize=0

histdata=histogram(data,binsize=binsize,_extra=extra,omin=omin,omax=omax)

if ~keyword_set(binsize) then binsize=(omax-omin)/n_elements(histdata)


IF KeyWord_Set(xrange) then xind=xrange else xind=[omin,omax]
IF ~KeyWord_Set(oplot) Then Begin
	plot,xind,[0,max(histdata)*1.1],/nodata,_extra=extra
ENDIF

start = omin
endpt = start + binsize

FOR j=0,N_Elements(histdata)-1 DO BEGIN
        x = [start, start, endpt, endpt, start]
        y = [0, histdata[j], histdata[j], 0, 0]
	if keyword_set(fillcolor) then begin
		polyfill,x,y,color=fillcolor,/data
	ENDIF
	OPLOT, x, y, _extra=extra
        start = start + binsize
        endpt = start + binsize
ENDFOR

end

