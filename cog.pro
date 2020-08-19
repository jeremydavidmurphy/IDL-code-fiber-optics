PRO cog, list

; This routine is just for plotting purposes. It accepts a list of text files
; (the fibername_FRD.txt files created in frd.pro) and does some plotting and
; analysis on them.)

;*************************************************************************
plotps = 'on'
;*************************************************************************

readcol,list,f='a',files,silent=1
n0 = n_elements(files)

colors = intarr(n0)
cs = floor(255.0/n0)
for j=0,n0-1 do colors[j] =  (j*cs)       

set_plot,'x'
for j=0,n0-1 do begin
   filename = files[j]
   readcol,filename,f='f,f',ee,fnum,silent=1
   n1 = n_elements(ee)
   fibername = strsplit(filename,'_',/extract)
   fibername = fibername[0]

   if (j eq 0) then begin
      window,0,retain=2
      device,decomposed=0
      loadct,27,/silent
      plot,fnum,ee,xrange=[10,3],xtitle='Output F/#',ytitle='Encircled Energy',$
           /nodata,xstyle=1,ystyle=1,yrange=[10,100]
      loadct,4,/silent
      oplot,[3.65,3.65],[10,100],color=110,thick=2
      oplot,[3.47,3.47],[10,100],color=150,thick=2
      oplot,[8.26,8.26],[10,100],color=180,thick=2
      loadct,33,/silent
   endif
   oplot,fnum,ee,color=colors[j],psym=-1
   wait,1.0
   xyouts,0.12,(0.90-(j*0.03)),fibername,color=colors[j],/normal,charsize=1.5
endfor

;now the ps file is generated
if (plotps eq 'on') then begin
   set_plot,'ps'
   for j=0,n0-1 do begin
      filename = files[j]
      readcol,filename,f='f,f',ee,fnum,silent=1
      n1 = n_elements(ee)
      fibername = strsplit(filename,'_',/extract)
      fibername = fibername[0]
      if (j eq 0) then begin
         device,file='compareFnum.ps',/color
         loadct,0,/silent
         plot,fnum,ee,xrange=[10,3],xtitle='Output F/#',ytitle='Encircled Energy',$
              /nodata,xstyle=1,ystyle=1,charthick=2,xthick=2,ythick=2,yrange=[10,100],$
              title='Comparison of Output FRD'
         loadct,33,/silent
      endif
      oplot,fnum,ee,color=colors[j],psym=-1,thick=1.0
      xyouts,0.17,(0.89-(j*0.03)),fibername,color=colors[j],/normal,charthick=2
      if (j eq n0-1) then begin
         loadct,4,/silent
         oplot,[3.65,3.65],[10,100],color=110,thick=2
         oplot,[3.47,3.47],[10,100],color=150,thick=2
         oplot,[8.26,8.26],[10,100],color=180,thick=2
         device,/close_file
      endif
   endfor
endif

print,'The next ENTER deletes the figure!'
pause
while !d.window ne -1 do wdelete, !d.window

stop
END
