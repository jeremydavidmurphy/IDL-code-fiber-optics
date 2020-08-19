PRO Pfrd,file

; This routine is used to plot a set of FRD values...

nfibers = 6

n1 = strarr(1)
n2 = n1 & n3 = n1 & n4 = n1 & n5 = n1 & n6 = n1
d1 = fltarr(30)
d2 = d1 & d3 = d1 & d4 = d1 & d5 = d1 & d6 = d1
f1 = fltarr(30)
f2 = f1 & f3 = f1 & f4 = f1 & f5 = f1 & f6 = f1


data = fltarr(30,2,nfibers)
names = strarr(nfibers)
colors = intarr(nfibers)
for j=0,nfibers-1 do colors[j] = 60 + (j * floor(165.0/nfibers))

openr,5,file
readf,5,n1,d1,f1,n2,d2,f2,n3,d3,f3,n4,d4,f4,n5,d5,f5,n6,d6,f6
free_lun,5
names = [n1,n2,n3,n4,n5,n6]
data[*,0,0] = d1
data[*,1,0] = f1
data[*,0,1] = d2
data[*,1,1] = f2
data[*,0,2] = d3
data[*,1,2] = f3
data[*,0,3] = d4
data[*,1,3] = f4
data[*,0,4] = d5
data[*,1,4] = f5
data[*,0,5] = d6
data[*,1,5] = f6

yup = 1.0
ydn = 0.65
xup = max(data[*,1,*])+0.1
xdn = min(data[*,1,*])-0.1

pnames = ['No Bend','15 cm','10 cm','5 cm','3 cm','1.5 cm']
set_plot,'ps'
device,file='Pfrd.ps',/color
loadct,0,/silent

plot,data[*,1,0],data[*,0,0],psym=-1,yrange=[ydn,yup],xrange=[xup,xdn],/xstyle,/ystyle,$
  xtitle='Output F-ratio',ytitle='Encircled Energy',/nodata,xthick=3,ythick=3,charthick=4,$
  charsize=1.3,thick=3,title='Bend Radius Tests for F&T 266 micron Fiber'
loadct,4,/silent
for j=0,nfibers-1 do begin
    oplot,data[*,1,j],data[*,0,j],psym=-1,color=colors[j],thick=3
    ttt = strsplit(names[j],'.',/extract)
    xyouts,0.2,0.85-(j*0.04),/normal,pnames[j],charthick=3,charsize=1.0,color=colors[j]
endfor
device,/close_file
set_plot,'x'



stop
END
