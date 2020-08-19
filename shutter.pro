pro shutter,list

readcol,list,format='a',files

test = readfits(files[0])
n1 = n_elements(test[*,0])
n2 = n_elements(test[0,*])
n3 = n_elements(files)

data = dblarr(n1,n2,n3)
pxa = dblarr(n1,n3)
pya = pxa
pxm = pxa
pym = pxa
std_prox = dblarr(n3)
std_proy = dblarr(n3)

for j=0,n3-1 do begin
    data[*,*,j] = readfits(files[j])
endfor

for j=0,n1-1 do begin
    for k=0,n3-1 do begin
        pxm[j,k] = median(data[j,*,k])
        pym[j,k] = median(data[*,j,k])
        pxa[j,k] = mean(data[j,*,k])
        pya[j,k] = mean(data[*,j,k])
    endfor
endfor

indxa = bsort(pxa[30,*])
indya = bsort(pya[30,*])
indxm = bsort(pxm[30,*])
indym = bsort(pym[30,*])
std_proxa = files[indxa]
std_proya = files[indya]
std_proxm = files[indxm]
std_proym = files[indym]


;set_plot,'x'

;window,0,retain=2
;device,decomposed=0

set_plot,'ps'
device,file='xa.ps',/color
loadct,0
plot,pxa[*,0],title='X-PROFILE (average)',xrange=[5,n1-5],xstyle=1,yrange=[1420,1460],ystyle=1,$
  xthick=2,ythick=2,charthick=2,charsize=1.3,/nodata
loadct,27
for j=0,n3-1 do oplot,pxa[*,j],color=j
device,/close_file

;window,1,retain=2
;device,decomposed=0
device,file='ya.ps',/color
loadct,0
plot,pya[*,0],title='Y-PROFILE (average)',xrange=[5,n2-5],xstyle=1,yrange=[1420,1460],ystyle=1,$
  xthick=2,ythick=2,charthick=2,charsize=1.3,/nodata
loadct,27
for j=0,n3-1 do oplot,pya[*,j],color=j
device,/close_file

;window,2,retain=2
;device,decomposed=0
device,file='xm.ps',/color
loadct,0
plot,pxm[*,0],title='X-PROFILE (median)',xrange=[5,n1-5],xstyle=1,yrange=[1420,1460],ystyle=1,$
  xthick=2,ythick=2,charthick=2,charsize=1.3,/nodata
loadct,27
for j=0,n3-1 do oplot,pxm[*,j],color=j
device,/close_file

;window,3,retain=2
;device,decomposed=0
device,file='ym.ps',/color
loadct,0
plot,pym[*,0],title='Y-PROFILE (median)',xrange=[5,n2-5],xstyle=1,yrange=[1420,1460],ystyle=1,$
  xthick=2,ythick=2,charthick=2,charsize=1.3,/nodata
loadct,27
for j=0,n3-1 do oplot,pym[*,j],color=j
device,/close_file

print,'The top five suspects are:'
for j=0,4 do begin
    print,std_proxa[j],' ',std_proya[j],'  ',std_proxm[j],'  ',std_proym[j]
endfor

;print,'The next ENTER deletes the plots..'
;pause
;wdelete,0,1,2,3

stop

end
