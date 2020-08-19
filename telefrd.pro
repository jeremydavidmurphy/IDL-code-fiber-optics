pro telefrd, list

p1 = 480
p2 = 500

;the pixel index for the 'bottom' and 'top' range
xb1 = 0
xb2 = 199
xt1 = 1855
xt2 = 2047

readcol, list ,format='a',files

n3 = n_elements(files)
cstep = floor(190/n3)

temp = readfits(files[0])
n1 = n_elements(temp[*,0])
n2 = n_elements(temp[0,*])

trans = fltarr(n2,n3)

data = fltarr(n1,n2,n3)
for j=0,n3-1 do begin
    print,'Reading in frame '+files[j]
    data[*,*,j] = readfits(files[j],/silent)
endfor

bot = fltarr(xb2-xb1+1,n3) & nbot = bot
top = fltarr(xt2-xt1+1,n3) & ntop = top
xbot = findgen(xb2-xb1+1) + 1
xtop = findgen(xt2-xt1+1) + xt1+1

for j=0,n3-1 do begin
    frame = data[*,*,j]
    strip = frame[p1:p2,*]
    trans[*,j] = median(strip,dimension=1)
;    trans[*,j] = total(strip,1)
    bot[*,j] = trans[xb1:xb2,j]
    top[*,j] = trans[xt1:xt2,j]
    nbot[*,j] = bot[*,j]/max(bot[*,j])
    ntop[*,j] = top[*,j]/max(top[*,j])
endfor

window,0,retain=2
device,decomposed=0
loadct,0
plot,xbot,nbot[*,0],title='Fibers # 1 to 22',xtitle='Pixel',/nodata,$
  yrange=[0,1.05],ystyle=1,xrange=[xb1-20,xb2],xstyle=1

loadct,27
for j=0,n3-1 do begin
    oplot,xbot,nbot[*,j],color=60+j*cstep,thick=1.5
    xyouts,0.12,0.90-(j*0.04),files[j],color=60+j*cstep,charsize=1.2,/normal
    pause
endfor

window,2,retain=2
device,decomposed=0
loadct,0
plot,xtop,ntop[*,0],title='Fibers # 226 to 246',xtitle='Pixel',/nodata,$
  yrange=[0,1.05],ystyle=1,xstyle=1,xrange=[xt1,xt2+20]

loadct,27
for j=0,n3-1 do begin
    oplot,xtop,ntop[*,j],color=60+j*cstep,thick=1.5
    xyouts,0.85,0.90-(j*0.04),files[j],color=60+j*cstep,charsize=1.2,/normal
    pause
endfor
print,'Next ENTER deletes the plots...'
pause
wdelete,0,2

stop
end
