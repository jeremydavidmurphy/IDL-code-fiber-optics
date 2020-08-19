PRO vp27stability, list

; This routine is used to compare flat field frames from various
; nights.

readcol,list,f='a',files
n0 = n_elements(files)
temp = fltarr(n0)
array = fltarr(2048,n0)
subbed = fltarr(2048,n0-1)

for j=0,n0-1 do begin
   frame = readfits(files[j],h,/silent)
   n1 = n_elements(frame[*,100])
   if n1 lt 1200 then size = 's' else size = 'l'
   temp[j] = sxpar(h,'INSTRTEM',count=c)
   if size eq 's' then slice = frame[450:550,*]
   if size eq 'l' then slice = frame[900:1000,*]
   mslice = median(slice,dim=1)
   nslice = mslice/max(mslice)
   array[*,j] = nslice
   if j gt 0 then subbed[*,j-1] = 100*(array[*,j]-array[*,j-1])
endfor

n1 = 0
n2 = 100
for j=0,n0-1 do begin
   set_plot,'ps'
   name = 'VP27_flats-'+strn(j+1)+'.ps'
   device,file=name,/color
   loadct,0
   plot,array[n1:n2,0],ystyle=9,yrange=[0,1],ytitle='Normalized Transmission',$
        position=[0.1,0.12,0.9,0.95],/nodata
   oplot,array[n1:n2,0],color=60
;   oplot,array[n1:n2,1],color=150
   if j gt 0 then begin
      loadct,4
      axis,yaxis=1,ytitle='Difference (%)',yrange=[-5,5],/ys,/save,color=150
      oplot,[0,100],[0,0],color=150
      oplot,sub[n1:n2,j],color=150,thick=3
;      device,/close
;   set_plot,'x'
;   n1 = n2
;   n2 = n2 + 100
;pause
endfor
device,/close
;window,0,retain=2,xsize=700,ysize=500
;loadct,0
;plot,array[0:100,0]
;loadct,33
;for j=0,n0-1 do begin
;   print,files[j],' ',strn(temp[j])
;   oplot,array[0:100,j],color=j*8
;   pause
;endfor
stop
END

