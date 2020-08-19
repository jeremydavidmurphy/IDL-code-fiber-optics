
pro telefrdQ, list

;this is a Quantitative version of the frd tests on telescope. 

;p1 = 150
;p2 = 170
;p1 = 480
;p2 = 500
;p1 = 900
;p2 = 920

;the pixel index for the 'bottom' and 'top' range
xb1 = 0
xb2 = 2000
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

ans1=''
ans2=''
print,'Plot top or bottom? (t/b)'
read,ans1

print,'Normalized locally? (y/n)'
read,ans2

for k=0,2 do begin
if(k eq 0) then p1=45
if(k eq 0) then p2=70
if(k eq 1) then p1=480
if(k eq 1) then p2=500
if(k eq 2) then p1=930
if(k eq 2) then p2=950

if(ans2 eq 'n') then begin
for j=0,n3-1 do begin
    frame = data[*,*,j]
    strip = frame[p1:p2,*]
    trans[*,j] = median(strip,dimension=1)
    bot[*,j] = trans[xb1:xb2,j]
    top[*,j] = trans[xt1:xt2,j]
;    nbot[*,j] = bot[*,j]/median(trans[*,j])
;    ntop[*,j] = top[*,j]/median(trans[*,j])
    nbot[*,j] = bot[*,j]/max(trans[*,j])
    ntop[*,j] = top[*,j]/max(trans[*,j])
endfor
endif
if(ans2 eq 'y') then begin
for j=0,n3-1 do begin
    frame = data[*,*,j]
    strip = frame[p1:p2,*]
    trans[*,j] = median(strip,dimension=1)
    bot[*,j] = trans[xb1:xb2,j]
    top[*,j] = trans[xt1:xt2,j]
    nbot[*,j] = bot[*,j]/max(bot[*,j])
    ntop[*,j] = top[*,j]/max(top[*,j])
endfor
endif


if (k eq 0) then window,0,retain=2
if (k eq 1) then window,1,retain=2
if (k eq 2) then window,2,retain=2
if (k eq 0) then ttl='BLUE'
if (k eq 1) then ttl='GREEN'
if (k eq 2) then ttl='RED'

;window,0,retain=2
device,decomposed=0
loadct,0

if (ans1 eq 'b') then begin
    plot,xbot,nbot[*,0],title=ttl+' Fibers # 246 to 224',xtitle='Pixel',/nodata,$
      yrange=[0,1.05],ystyle=1,xrange=[xb1-20,xb2],xstyle=1
    
    loadct,27
    for j=0,n3-1 do begin
        oplot,xbot,nbot[*,j],color=60+j*cstep,thick=1.5
        xyouts,0.12,0.90-(j*0.04),files[j],color=60+j*cstep,charsize=1.2,/normal
        pause
    endfor
endif

if (ans1 eq 't') then begin
    window,0,retain=2
    device,decomposed=0
    loadct,0
    plot,xtop,ntop[*,0],title=ttl+' Fibers # 22 to 1',xtitle='Pixel',/nodata,$
      yrange=[0,1.05],ystyle=1,xstyle=1,xrange=[xt1,xt2+20]
    
    loadct,27
    for j=0,n3-1 do begin
        oplot,xtop,ntop[*,j],color=60+j*cstep,thick=1.5
        xyouts,0.85,0.90-(j*0.04),files[j],color=60+j*cstep,charsize=1.2,/normal
        pause
    endfor
endif
endfor
;print,'Enter a first guess for fiber #1 peak position:'
;read,ans


print,'Next ENTER deletes the plots...'
pause
wdelete,0,1,2

stop
end
