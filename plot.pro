
pro lc, ps=ps
	; Plots the lightcurve
	
	if keyword_set(ps) then open_ps,'lc.ps'
	
	readcol, 'out/prof', t, F, format=('F,F,X')

	plot, t, F, /xlog, /ylog, xtitle=textoidl('Time (s)'), $
			ytitle=textoidl('Luminosity (erg s^{-1})'),xrange=[30,3e4],$
			yrange=[3e36,2e38], xstyle=1,ystyle=1, linestyle=0, charsize=1.5,/nodata

	oplot, t, F
	
	if keyword_set(ps) then close_ps
	
end


pro prof2, delay=delay, png=png
	; shows a movie of the temperature profile evolving over time

;	window,0,retain=2, xsize=600, ysize=600

	!p.multi=[0,1,2,0,0]

	readcol, 'out/prof', tt,FF, format=('F,F,X,X,X')
	
	openr, lun, 'out/out', /get_lun

	; first find out how many grid points
	ngrid=0
	readf, lun, ngrid, format='(I0)'
	print, 'Number of grid points=',ngrid

	count=0

	while (not eof(lun)) do begin

		; get the time
		readf, lun, time
		print,"t=",time
		

		data=dblarr(13,ngrid)
		readf, lun, data
		y=data(0,*)
		T=data(1,*)
		F=data(2,*)
		beta=data(12,*)

		erase
		plot, y, T, /xlog, /ylog,charsize=1.2, ytitle=textoidl('T (K)'),$
				xtitle=textoidl('y (g cm^{-2})'), yrange=[1e8,1e10],ystyle=1, $
				xrange=[3d8,3d14], xstyle=1
		xyouts, 5d12, 5d9, 't='+string(time)+' s', charsize=1.4

		if keyword_set(delay) then begin
			for i=1L,delay do begin
			endfor
		endif

		tt2=tt[where(tt le time)]
		ff2=ff[where(tt le time)]
		plot, tt,FF, linestyle=1, /xlog, /ylog, xrange=[10.0,1e5], xstyle=1,$
					ytitle=textoidl('Flux (erg cm^{-2} s^{-1})'), $
					xtitle=textoidl('Time (s)'), charsize=1.2,$
					yrange=[2d36,2d38], ystyle=1
		if (n_elements(tt2) gt 1) then begin
				oplot,tt2,FF2,linestyle=0
		endif

		count++

		if keyword_set(png) then begin
			image = TVRD(0,0,!D.X_Size,!D.Y_Size,True=1, Order=order)
			filename=string(format='("png/",I03,".png")',count)
			Write_PNG,filename,image
		endif

	endwhile

	free_lun,lun

	!p.multi=[0,1,1,0,0]

end






pro open_ps, name
	!p.font=0
	set_plot, 'ps'
	device,filename=name,/color,/times,/encapsul
	!p.multi=[0,1,1,0,0]
	!p.charsize=1.2
	!p.thick=3
	!x.thick=3
	!y.thick=3
end

pro close_ps
	device,/close
    set_plot,'x'
	!p.font=-1
	!p.thick=1
	!x.thick=1
	!y.thick=1	
end



