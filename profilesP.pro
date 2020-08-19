; This routine is designed to read in a list of spots, optionally
; subtract a dark, and overplot a fiber profile.

;COMPILE: locatecenterF

;DATALIST: The list of fits files you want to work on
;TEMPLIST: The temperature output.

pro profilesP,datalist, templist, DARK=dark

readcol,templist,f='x,x,a,x,f',ttime,ttemp
n2 = n_elements(ttime)
tttime = fltarr(n2)
for j=0,n2-1 do begin
    h = strsplit(ttime[j],':',/extract)
    hr = float(h[0])
    min = float(h[1])
    sec = float(h[2])
    tttime[j] = hr + (min/60.0) + (sec/3600.0)
endfor

readcol,datalist,f='a',files,silent=1
n0 = n_elements(files)
;
files = reverse(files)
ftemp1 = fltarr(n0)
ftemp2 = fltarr(n0)

if (n_elements(dark) ne 0) then begin
    darkframe = readfits(dark,/silent)
    darkS = 'on'
endif else darkS = 'off'

for j=0,n0-1 do begin
    print,'Working on frame '+files[j]
    frame = readfits(files[j],h,/silent)
    if (j eq 0) then begin
        n1 = n_elements(frame[*,0])
        profiles = fltarr(n1,n0)
        times = fltarr(n0)
    endif
    if (darkS eq 'on') then frame = frame - darkframe
    cent = locatecenterF(frame)
    profiles[*,j] = median(frame[*,cent[1]-10:cent[1]+10],dim=2)

    cntr = 0
    out = 'lost'
    repeat begin
        hh = strsplit(h[cntr],' ',/extract)
        if (hh[0] eq 'DATE-OBS=') then out = 'found' else cntr = cntr + 1
    endrep until (out eq 'found')
    hhh = strsplit(hh[1],'T',/extract)
    hhh = hhh[1]
    hhhh = strsplit(hhh,':',/extract)
    hr = float(hhhh[0])
    min = float(hhhh[1])
    sec = float(hhhh[2])
    times[j] = hr - 5 + (min/60.0) + (sec/3600.0)

    for k=0,n2-1 do begin
        thistime = times[j]-12.0
        i1 = where(tttime le thistime)
        i2 = where(tttime ge thistime)
        t1 = mean([ttemp[max(i1)],ttemp[min(i2)]])
        t2 = mean([tttime[max(i1)],tttime[min(i2)]])
    endfor
    ftemp1[j] = t1
    print,files[j]+' '+strn(thistime)+' '+strn(t1)+' '+strn(t2)
endfor

ftemp2 = interpol(ttemp,tttime,times-12)

colors = intarr(n0)
for j=0,n0-1 do colors[j] = 60 + (j * floor(165.0/n0))

set_plot,'x'
window,0,retain=2
device,decomposed=0
loadct ,0
plot,profiles[*,0],xtitle='Pixels',ytitle='CCD Counts',$
  yrange=[min(profiles),max(profiles)+1000],ystyle=1,/nodata,$
  xrange=[100,1400],xstyle=1
loadct,4
;for j=0,n0-1 do begin
;    oplot,profiles[*,j],color=colors[j]
;    if (j gt 0 and j lt 24) then xyouts,0.12,0.93-(j*0.03),files[j],color=colors[j],/normal
;    if (j gt 25 and j lt 49) then xyouts,0.22,0.9-((j-26)*0.03),files[j],color=colors[j],/normal
;    if (j gt 50 and j lt 80) then xyouts,0.85,0.9-((j-51)*0.03),files[j],color=colors[j],/normal
;endfor
for j=0,n0-1 do begin
    oplot,profiles[*,j],color=ftemp1[j]*3
    if (j lt 24) then xyouts,0.12,0.93-(j*0.03),files[j],color=ftemp1[j]*3,/normal
    if (j gt 25 and j lt 49) then xyouts,0.22,0.9-((j-26)*0.03),files[j],color=ftemp1[j]*3,/normal
    if (j gt 50 and j lt 80) then xyouts,0.85,0.9-((j-51)*0.03),files[j],color=ftemp1[j]*3,/normal
;    wait,1
endfor

window,2,retain=2
device,decomposed=0
loadct ,0
plot,profiles[*,0],xtitle='Pixels',ytitle='CCD Counts',$
  yrange=[0,10000],ystyle=1,/nodata,$
  xrange=[1000,1200],xstyle=1
loadct,4
for j=0,n0-1 do begin
    oplot,profiles[*,j],color=ftemp1[j]*3
    if (j lt 24) then xyouts,0.12,0.93-(j*0.03),files[j],color=ftemp1[j]*3,/normal
    if (j gt 25 and j lt 49) then xyouts,0.22,0.9-((j-26)*0.03),files[j],color=ftemp1[j]*3,/normal
    if (j gt 50 and j lt 80) then xyouts,0.85,0.9-((j-51)*0.03),files[j],color=ftemp1[j]*3,/normal
;    wait,1
endfor

;pause
while !d.window ne -1 do wdelete, !d.window

set_plot,'ps'
device,file='profiles1.ps',/color
loadct ,0
plot,profiles[*,0],xtitle='Pixels',ytitle='CCD Counts',$
  yrange=[min(profiles),max(profiles)+1000],ystyle=1,/nodata,$
  xthick=3,ythick=3,charthick=3,xrange=[0,1300],xstyle=1,$
  position=[0.15,0.13,0.95,0.91],title='VP3 Temperature Tests'
xyouts,0.18,0.83+0.035,'Frame        Temp',/normal
loadct,4
for j=0,n0-1 do begin
    oplot,profiles[*,j],color=ftemp1[j]*3
    xyouts,0.18,0.83-(j*0.035),files[j],color=ftemp1[j]*3,/normal
    xyouts,0.33,0.83-(j*0.035),strn(ftemp2[j]),color=ftemp1[j]*3,/normal
endfor
device,/close_file

device,file='profiles2.ps',/color
loadct ,0
plot,profiles[*,0],xtitle='Pixels',ytitle='CCD Counts',$
  yrange=[0,10000],ystyle=1,/nodata,xrange=[1000,1200],xstyle=1,$
  xthick=3,ythick=3,charthick=3,position=[0.13,0.13,0.95,0.91],$
  title='VP3 Temperature Tests'
xyouts,0.18,0.83+0.035,'Frame        Temp',/normal
loadct,4
for j=0,n0-1 do begin
    oplot,profiles[*,j],color=ftemp1[j]*3
    xyouts,0.18,0.83-(j*0.035),files[j],color=ftemp1[j]*3,/normal
    xyouts,0.33,0.83-(j*0.035),strn(ftemp2[j]),color=ftemp1[j]*3,/normal
endfor
device,/close_file

set_plot,'x'
stop
end
