; This routine plots profiles of fiber spots. It uses the
; locatecenterF routine to determine the center.

;COMPILE locatecenterF
PRO profileP, list

overplot = 'y' ;set to something else if you do not want to overplot all the spots.
readcol,list,f='a',silent=1,files

n0 = n_elements(files)
center = fltarr(4,n0)
colors = intarr(n0)
for j=0,n0-1 do colors[j] = 60 + (j * floor(165.0/n0))

window,0,retain=2
device,decomposed=0
window,2,retain=2
device,decomposed=0
for j=0,n0-1 do begin
    frame = readfits(files[j],/silent)
    nx = n_elements(frame[*,0])
    ny = n_elements(frame[0,*])
    if (j eq 0) then begin
        xslice = fltarr(nx,n0)
        yslice = fltarr(ny,n0)
    endif
    center[*,j] = locatecenterF(frame)
    xslice[*,j] = median(frame[*,center[1,j]-10:center[1,j]+10],dim=2)
    yslice[*,j] = median(frame[center[0,j]-10:center[0,j]+10,*],dim=1)
    wset,2
    if j eq 0 then begin
        loadct,0,/silent
        plot,xslice[*,j],title='X AXIS',xrange=[0,nx],/xstyle,xtitle='X AXIS'
        xyouts,0.15,0.9,files[j],/normal,charsize=1.2
    endif
    if j gt 0 and overplot eq 'y' then begin
        loadct,4,/silent
        oplot,xslice[*,j],color=colors[j]
        xyouts,0.15,0.9-(j*0.05),files[j],/normal,color=colors[j],charsize=1.2
    endif else begin
        loadct,0,/silent
        plot,xslice[*,j],title='X AXIS',xrange=[0,nx],/xstyle
        xyouts,0.15,0.9,files[j],/normal,charsize=1.2
    endelse
    wset,0
    if j eq 0 then begin
        loadct,0,/silent
        plot,yslice[*,j],title='Y AXIS',xrange=[0,ny],/xstyle,xtitle='Y AXIS'
        xyouts,0.15,0.9,files[j],/normal,charsize=1.2
    endif
    if j gt 0 and overplot eq 'y' then begin
        loadct,4,/silent
        oplot,yslice[*,j],color=colors[j] 
        xyouts,0.15,0.9-(j*0.05),files[j],/normal,color=colors[j],charsize=1.2
        endif else begin
        loadct,0,/silent
        plot,yslice[*,j],title='Y AXIS',xrange=[0,ny],/xstyle
        xyouts,0.15,0.9,files[j],/normal,charsize=1.2
    endelse
endfor

print,'The next ENTER deletes the figures...'
pause
while !d.window ne -1 do wdelete, !d.window

ans = ''
print,'Plot to a PS file? (y/n):'
read,ans

if (ans eq 'y') then begin
    set_plot,'ps'
    device,file='profile.ps',/color
    !p.multi = [0,2,1,0,1]
    loadct,0,/silent
    plot,xslice[*,0],yrange=[0,max(xslice)],xtitle='X AXIS',position=[0.08,0.1,0.48,0.92],$
      charsize=0.6,xrange=[0,nx],/xstyle,/nodata
    loadct,4,/silent
    for j=0,n0-1 do begin
        oplot,xslice[*,j],color=colors[j]
        xyouts,0.1,0.87-(j*0.025),files[j],color=colors[j],/normal,charsize=0.6
    endfor
    loadct,0,/silent
    plot,yslice[*,0],yrange=[0,max(yslice)],xtitle='Y AXIS',position=[0.55,0.1,0.95,0.92],$
      charsize=0.6,xrange=[0,ny],/xstyle,/nodata
    loadct,4,/silent
    for j=0,n0-1 do begin
        oplot,yslice[*,j],color=colors[j]
    endfor
    device,/close_file
    !p.multi = [0,0,0,0,0]
    set_plot,'x'
    !p.multi = [0,0,0,0,0]
endif else print,'Very Well...'


stop
END
