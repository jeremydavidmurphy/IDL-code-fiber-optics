PRO Pfib_trans, list

; This routine is used to plot the output of transmission.pro, which
; calculates the transmission values for optical fibers.

; LIST: This is a list of all the fiber_trans.txt files to be plotted

readcol,list,format='a',flist
n0 = n_elements(flist)

colors = intarr(n0)
for j=0,n0-1 do colors[j] = floor(195./n0)*j+40

data = fltarr(9,n0,3)
bg = fltarr(9,n0,3)
names = strarr(n0)
for j=0,n0-1 do begin
    t = strsplit(flist[j],'_',/extract)
    names[j] = t[0]
endfor
for j=0,n0-1 do begin
    readcol,flist[j],format='f',a ,b,c,d,e,f,g,h,i
    if (j eq 0) then wave = [a[0],b[0],c[0],d[0],e[0],f[0],g[0],h[0],i[0]]
    data[*,j,0] = [a[5],b[5],c[5],d[5],e[5],f[5],g[5],h[5],i[5]]
    data[*,j,1] = [a[6],b[6],c[6],d[6],e[6],f[6],g[6],h[6],i[6]]
    data[*,j,2] = [a[7],b[7],c[7],d[7],e[7],f[7],g[7],h[7],i[7]]
    bg[*,j,0] = [a[2],b[2],c[2],d[2],e[2],f[2],g[2],h[2],i[2]]
    bg[*,j,1] = [a[4],b[4],c[4],d[4],e[4],f[4],g[4],h[4],i[4]]
    bg[*,j,2] = bg[*,j,0] / bg[*,j,1]
endfor
stop
set_plot,'x'
window,0,retain=2,xsize=420,ysize=350,xpos=0,ypos=55
device,decomposed=0
loadct,0
plot,wave,data[*,0,0],xrange=[min(wave)-20,max(wave)+20],$
  xstyle=1,xtitle='wavelength (A)',ytitle='Transmission',$
  title='Unsubtracted Transmission Values',psym=-1,$
  yrange=[min(data[*,*,0])-0.05,max(data[*,*,0])+0.05],ystyle=1
loadct,27
for j=1,n0-1 do begin
    oplot,wave,data[*,j,0],psym=-1,color=colors[j-1]
endfor

window,1,retain=2,xsize=420,ysize=350,xpos=420,ypos=55
device,decomposed=0
loadct,0
plot,wave,data[*,0,1],xrange=[min(wave)-20,max(wave)+20],$
  xstyle=1,xtitle='wavelength (A)',ytitle='Transmission',$
  title='Subtracted Transmission Values',psym=-1,$
  yrange=[min(data[*,*,1])-0.05,max(data[*,*,1])+0.05],ystyle=1
loadct,27
for j=1,n0-1 do begin
    oplot,wave,data[*,j,1],psym=-1,color=colors[j-1]
endfor

window,2,retain=2,xsize=420,ysize=350,xpos=840,ypos=55
device,decomposed=0
loadct,0
plot,wave,data[*,0,2],xrange=[min(wave)-20,max(wave)+20],$
  xstyle=1,xtitle='wavelength (A)',ytitle='Transmission',$
  title='Subtracted Transmission Values',psym=-1,$
  yrange=[min(data[*,*,2])-0.05,max(data[*,*,2])+0.05],ystyle=1
loadct,27
for j=1,n0-1 do begin
    oplot,wave,data[*,j,2],psym=-1,color=colors[j-1]
endfor

wset,0
loadct,0
xyouts,0.8,0.5,names[0],/normal
loadct,27
for j=1,n0-1 do xyouts,0.8,0.5-(j*0.03),names[j],color=colors[j],$
  /normal

wset,1
loadct,0
xyouts,0.8,0.5,names[0],/normal
loadct,27
for j=1,n0-1 do xyouts,0.8,0.5-(j*0.03),names[j],color=colors[j],$
  /normal

wset,2
loadct,0
xyouts,0.8,0.5,names[0],/normal
loadct,27
for j=1,n0-1 do xyouts,0.8,0.5-(j*0.03),names[j],color=colors[j],$
  /normal

pause
wdelete,0,1,2

ans = ''
print,'Make some ps files? (y/n)'
read,ans
if (ans eq 'y') then begin
    print,'Which one? (1,2,3)'
    read,ans
    if (ans eq '1') then titl='Unsubtracted Transmission Values'
    if (ans eq '2') then titl='Subtracted Transmission Values'
    if (ans eq '3') then titl='Subtracted Transmission Values'
    ans = uint(ans)-1
    
    dataout = data[*,*,ans]
    set_plot,'ps'
    device,file='trans.ps',/color
    loadct,0
    plot,wave,dataout[*,0],psym=-1,title=titl,$
      xrange=[min(wave)-20,max(wave)+20],xstyle=1,$
      xtitle='wavelength (A)',ytitle='Transmission',$
      yrange=[min(dataout)-0.05,max(dataout)+0.05],ystyle=1,$
      xthick=3,ythick=3,charthick=3,thick=2
    xyouts,0.8,0.5,names[0],/normal,charthick=2
    loadct,4
    for j=1,n0-1 do begin
        oplot,wave,dataout[*,j],psym=-1,color=colors[j-1],thick=2
        xyouts,0.8,0.5-(j*0.03),names[j],color=colors[j],/normal,charthick=2
    endfor
    device,/close_file
endif

stop
END
