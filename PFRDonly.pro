; This routine is used to plot the output of the trans_FRD routine.

PRO PFRDonly, FILE=file, TITLE=title1

; FILE: The trans_FRD_all.dat file (NOT the trimmed one...). If not
; used, the code will look for 'trans_FRD_all.dat'

; This code plots the output of trans_FRD.pro. It generates 3 ps
; files. A transmission figure, An FRD COG figure showing all
; wavelengths, and a FRD COG figure for just 480 nm (to cut down on
; the glare). To set this up, you must have a list of fibers called 

; THIS CODE NEEDS A LIST OF THE FIBERS (CALLED FIBERDIR.LIST) AND A
; LIST OF THE WAVELENGTHS (CALLED WAVE.LIST) IN THE CALLING DIRECTORY.
offset = 0.2 ;this is a number that should be ADDED to the f/ratio values 
               ;from the BEFORE tests to correct for a systematic offset of the input f/#
;offset = 0.0
;***********************************************************************
stepsize = 364.0
if (n_elements(file) eq 0) then file = 'trans_FRD_all.dat'

readcol,'fiberdir.list',format='a',fibers,/silent
n0 = n_elements(fibers)

readcol,'wave.list',format='f',wavelength,/silent
n1 = n_elements(wavelength)

iEE = findgen(90)*3.0
iFnum = iEE + 1
ichi = iEE + 2
EE = dblarr(n1,90,n0)
Fnum = EE
chi2 = EE

for j=0,n0-1 do begin ;a loop over each fiber tested
    fiber = fibers[j]
    skp = j * stepsize
    if (n1 eq 1) then begin
        readcol,file,comment='w',format='x,d',$
          w350,skipline=skp,/silent
        data = [[w350]]
    endif
    if (n1 eq 2) then begin
        readcol,file,comment='w',format='x,d,d',$
          w350,w365,skipline=skp,/silent
        data = [[w350],[w365]]
    endif
    if (n1 eq 3) then begin
        readcol,file,comment='w',format='x,d,d,d',$
          w350,w365,w380,skipline=skp,/silent
        data = [[w350],[w365],[w380]]
    endif
    if (n1 eq 4) then begin
        readcol,file,comment='w',format='x,d,d,d,d',$
          w350,w365,w380,w400,skipline=skp,/silent
        data = [[w350],[w365],[w380],[w400]]
    endif
    if (n1 eq 5) then begin
        readcol,file,comment='w',format='x,d,d,d,d,d',$
          w350,w365,w380,w400,w440,skipline=skp,/silent
        data = [[w350],[w365],[w380],[w400],[w440]]
    endif
    if (n1 eq 6) then begin
        readcol,file,comment='w',format='x,d,d,d,d,d,d',$
          w350,w365,w380,w400,w440,w480,skipline=skp,/silent
        data = [[w350],[w365],[w380],[w400],[w440],[w480]]
    endif
    if (n1 eq 7) then begin
        readcol,file,comment='w',format='x,d,d,d,d,d,d,d',$
          w350,w365,w380,w400,w440,w480,w520,skipline=skp,/silent
        data = [[w350],[w365],[w380],[w400],[w440],[w480],[w520]]
    endif
    if (n1 eq 8) then begin
        readcol,file,comment='w',format='x,d,d,d,d,d,d,d,d',$
          w350,w365,w380,w400,w440,w480,w520,w560,skipline=skp,/silent
        data = [[w350],[w365],[w380],[w400],[w440],[w480],[w520],[w560]]
    endif
    if (n1 eq 9) then begin
        readcol,file,comment='w',format='x,d,d,d,d,d,d,d,d,d',$
          w350,w365,w380,w400,w440,w480,w520,w560,w600,skipline=skp,/silent
        data = [[w350],[w365],[w380],[w400],[w440],[w480],[w520],[w560],[w600]]
    endif

    n2 = n_elements(w350)

    data = transpose(data)
    for k=0,n1-1 do EE[k,*,j] = data[k,iEE] * 100.0
    for k=0,n1-1 do Fnum[k,*,j] = data[k,iFnum] + offset
    for k=0,n1-1 do chi2[k,*,j] = data[k,ichi]
endfor

colors = findgen(n0) * floor(255.0/n0)
if n1 gt 1 then colors2 = findgen(n1) * floor(255.0/n1) else colors2 = 150

for j=0,n0-1 do begin
    if (j eq 0) then begin
        loadct,0
        plot,Fnum[0,*,j],EE[0,*,j],xrange=[10,2],/xstyle,yrange=[0,100],/ystyle,/nodata,$
          xtitle='Output F-ratio',ytitle='% Encircled Energy'
        oplot,[3.3,3.3],[0,100],thick=2
        oplot,[3.65,3.65],[0,100],thick=2
        loadct,33
    endif
    for k=0,n1-1 do oplot,Fnum[k,*,j],EE[k,*,j],color=colors2[k]
    wait,0.2
endfor

set_plot,'ps'
device,file='FRDall.ps',/color
for j=0,n0-1 do begin
    if (j eq 0) then begin
        loadct,0
        plot,Fnum[0,*,j],EE[0,*,j],xrange=[10,2],/xstyle,yrange=[0,100],/ystyle,/nodata,$
          xtitle='Output F-ratio',ytitle='% Encircled Energy',xthick=2,ythick=2,$
          charthick=2,thick=2,title=title1
        xyouts,3.69,10,'F/#: 3.65',orientation=90,charsize=0.8,charthick=2
        xyouts,3.34,10,'F/#: 3.30',orientation=90,charsize=0.8,charthick=2
        loadct,4
        oplot,[3.3,3.3],[0,100],thick=2,color=150
        oplot,[3.65,3.65],[0,100],thick=2,color=110
        loadct,33
    endif
    for k=0,n1-1 do xyouts,0.17,0.9-(k*0.03),strn(wavelength[k])+' nm',/normal,color=colors2[k],charthick=2,charsize=0.8
    for k=0,n1-1 do oplot,Fnum[k,*,j],EE[k,*,j],color=colors2[k],thick=1
endfor
device,/close_file

device,file='FRD480.ps',/color
for j=0,n0-1 do begin
    if (j eq 0) then begin
        loadct,0
        plot,Fnum[0,*,j],EE[0,*,j],xrange=[10,2],/xstyle,yrange=[0,100],/ystyle,/nodata,$
          xtitle='Output F-ratio',ytitle='% Encircled Energy',xthick=2,ythick=2,$
          charthick=2,thick=2,title=title1
        xyouts,3.69,10,'F/#: 3.65',orientation=90,charsize=0.8,charthick=2
        xyouts,3.34,10,'F/#: 3.30',orientation=90,charsize=0.8,charthick=2
        loadct,4
        oplot,[3.3,3.3],[0,100],thick=2,color=150
        oplot,[3.65,3.65],[0,100],thick=2,color=110
        loadct,33
    endif
    oplot,Fnum[0,*,j],EE[0,*,j],thick=1,color=colors[j]
    xyouts,0.17,0.9-(j*0.03),fibers[j]+':',/normal,color=colors[j],charthick=2,charsize=0.8
    fprint = sigfig(Fnum[0,89,j],4)
    xyouts,0.22,0.9-(j*0.03),'F/'+fprint+' (95%)',/normal,color=colors[j],charthick=2,charsize=0.8
endfor
device,/close_file
set_plot,'x'

pause
while !d.window ne -1 do wdelete, !d.window

stop
END
