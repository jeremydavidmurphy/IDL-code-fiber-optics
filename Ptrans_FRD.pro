; This routine is used to plot the output of the trans_FRD routine.

PRO Ptrans_FRD, FILE=file, TITLE=title1

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

bundle = 'vp6t'

;***********************************************************************
;turn me on to plot the raw VIRUS fiber transmission curves
plotraw = 'off'
comp = 'grad78'
if (plotraw eq 'on') then begin
    if (comp eq 'danu') then begin
        readcol,'/home/danu/murphy/fibers/VIRUS_fiber_transmission_curves/list',f='a',Tfiles
        readcol,'/home/danu/murphy/fibers/VIRUS_fiber_transmission_curves/'+Tfiles[0],skipline=1,silent=1,f='i,x,f',w1,t1
        readcol,'/home/danu/murphy/fibers/VIRUS_fiber_transmission_curves/'+Tfiles[1],skipline=1,silent=1,f='i,x,f',w2,t2
        readcol,'/home/danu/murphy/fibers/VIRUS_fiber_transmission_curves/'+Tfiles[2],skipline=1,silent=1,f='i,x,f',w3,t3
        readcol,'/home/danu/murphy/fibers/VIRUS_fiber_transmission_curves/'+Tfiles[3],skipline=1,silent=1,f='i,x,f',w4,t4
    endif
    if (comp eq 'grad78') then begin
        readcol,'/scratch/grad78/murphy/fibers/VIRUS_fiber_transmission_curves/list',f='a',Tfiles
        readcol,'/scratch/grad78/murphy/fibers/VIRUS_fiber_transmission_curves/'+Tfiles[0],skipline=1,silent=1,f='i,x,f',w1,t1
        readcol,'/scratch/grad78/murphy/fibers/VIRUS_fiber_transmission_curves/'+Tfiles[1],skipline=1,silent=1,f='i,x,f',w2,t2
        readcol,'/scratch/grad78/murphy/fibers/VIRUS_fiber_transmission_curves/'+Tfiles[2],skipline=1,silent=1,f='i,x,f',w3,t3
        readcol,'/scratch/grad78/murphy/fibers/VIRUS_fiber_transmission_curves/'+Tfiles[3],skipline=1,silent=1,f='i,x,f',w4,t4
    endif
endif
;***********************************************************************
readcol,'fiberdir.list',format='a',fibers,/silent
n0 = n_elements(fibers)

readcol,'wave.list',format='f',wavelength,/silent
n1 = n_elements(wavelength)

;readcol,'vp6t_trans_ave.txt',f='f,f,f',wave,vp6t1,vp6t2

if (bundle eq 'vp6t') then begin
    stepsize = 378.0
    iEE1 = findgen(8)*5.0
    iEE2 = findgen(82)*3.0 + 38.0
    iEE = [iEE1,iEE2]
    iFnum = iEE + 1
    ichi = iEE + 2
    itrans = findgen(7)*5.0 + 4.0
    ikept = itrans - 1.0
    trans = dblarr(n1,7,n0)
    keptflux = trans
    EE = dblarr(n1,90,n0)
    Fnum = EE
    chi2 = EE
endif
if (bundle eq 'vp5t') then begin
    stepsize = 138.0
    iEE1 = findgen(8)*5.0
    iEE2 = findgen(21)*3.0 + 38.0
    iEE = [iEE1,iEE2]
    iFnum = iEE + 1
    ichi = iEE + 2
    itrans = findgen(7)*5.0 + 4.0
    ikept = itrans - 1.0
    trans = dblarr(n1,7,n0)
    keptflux = trans
    EE = dblarr(n1,29,n0)
    Fnum = EE
    chi2 = EE
endif

if (n_elements(file) eq 0) then file = 'trans_FRD_all.dat'

for j=0,n0-1 do begin ;a loop over each fiber tested
    fiber = fibers[j]
    skp = j * stepsize
    readcol,file,comment='w',format='x,d,d,d,d,d,d,d,d,d',$
      w350,w365,w380,w400,w440,w480,w520,w560,w600,skipline=skp,/silent
    n2 = n_elements(w350)

    data = [[w350],[w365],[w380],[w400],[w440],[w480],[w520],[w560],[w600]]
    data = transpose(data)
    for k=0,n1-1 do trans[k,*,j] = data[k,itrans] * 100.0
    for k=0,n1-1 do keptflux[k,*,j] = data[k,ikept]
    for k=0,n1-1 do EE[k,*,j] = data[k,iEE] * 100.0
    for k=0,n1-1 do Fnum[k,*,j] = data[k,iFnum] + offset
    for k=0,n1-1 do chi2[k,*,j] = data[k,ichi]
endfor

; The F-out values are now smoothed...
for j=0,n0-1 do begin
    for k=0,n1-1 do Fnum[k,*,j] = smooth(Fnum[k,*,j],7)
endfor

colors = findgen(n0) * floor(255.0/n0)
window,0,retain=2
device,decomposed=0

loadct,0
plot,wavelength,median(trans[*,*,0],dim=2),psym=-1,xrange=[325,700],/xstyle,$
  yrange=[50,100],/nodata,xtitle='Wavelength (nm)',ytitle='% Fiber Transmission',$
  symsize=1.5

loadct,33
for j=0,n0-1 do begin
    oplot,wavelength,median(trans[*,*,j],dim=2),psym=-1,symsize=1.5,color=colors[j]
    for k=0,6 do oplot,wavelength,trans[*,k,j],psym=1,color=colors[j]
    print,fibers[j]
    xyouts,0.85,0.9-(j*0.03),fibers[j],/normal,color=colors[j]
    xyouts,0.90,0.9-(j*0.03),strn(median(trans[n1-1,*,j],/even)),/normal,color=colors[j]
    wait,0.2
endfor

if (n_elements(title1) ne 1) then title1='No Title'

;Now a ps file is made of the same thing...
set_plot,'ps'
device,file='TRANS.ps',/color
loadct,0
plot,wavelength,median(trans[*,*,0],dim=2),psym=-1,xrange=[325,750],/xstyle,$
  yrange=[50,100],/nodata,xtitle='Wavelength (nm)',ytitle='% Fiber Transmission',$
  symsize=1.5,xthick=2,ythick=2,charthick=2,thick=2,title=title1
if (plotraw eq 'on') then begin
    loadct,4
    oplot,w1,t1,thick=2,color=110
    oplot,w2,t2,thick=2,color=110
;    oplot,w3,t3,thick=2,color=110
    oplot,w4,t4,thick=2,color=110
endif
loadct,33
for j=0,n0-1 do begin
    oplot,wavelength,median(trans[*,*,j],dim=2),psym=-1,symsize=1.5,color=colors[j],$
      thick=2
    for k=0,6 do oplot,wavelength,trans[*,k,j],psym=1,color=colors[j],thick=2
    xyouts,0.72,0.9-(j*0.03),fibers[j],/normal,color=colors[j],charthick=2,charsize=0.8
    p = sigfig(median(trans[n1-1,*,j]),4)
    xyouts,0.77,0.9-(j*0.03),strn(p)+'% at 600 nm',/normal,charthick=2,color=colors[j],charsize=0.8
endfor
device,/close_file

set_plot,'x'
window,2,retain=2
device,decomposed=0
colors2 = findgen(n1) * floor(255.0/n1)

if (bundle eq 'vp6t') then begin
    yr=[0,100]
    xr=[10,3]
endif
if (bundle eq 'vp5t') then begin
    yr=[65,100]
    xr=[4.5,3]
endif

for j=0,n0-1 do begin
    if (j eq 0) then begin
        loadct,0
        plot,Fnum[0,*,j],EE[0,*,j],xrange=[xr[0],xr[1]],/xstyle,yrange=[yr[0],yr[1]],/ystyle,/nodata,$
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
        plot,Fnum[0,*,j],EE[0,*,j],xrange=[xr[0],xr[1]],/xstyle,yrange=[yr[0],yr[1]],/ystyle,/nodata,$
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

if (n1 gt 1) then f = n1-1 else f = 0
device,file='FRD.ps',/color
for j=0,n0-1 do begin
    if (j eq 0) then begin
        loadct,0
        plot,Fnum[0,*,j],EE[0,*,j],xrange=[xr[0],xr[1]],/xstyle,yrange=[yr[0],yr[1]],/ystyle,/nodata,$
          xtitle='Output F-ratio',ytitle='% Encircled Energy',xthick=2,ythick=2,$
          charthick=2,thick=2,title=title1
        xyouts,3.69,10,'F/#: 3.65',orientation=90,charsize=0.8,charthick=2
        xyouts,3.34,10,'F/#: 3.30',orientation=90,charsize=0.8,charthick=2
        loadct,4
        oplot,[3.3,3.3],[0,100],thick=2,color=150
        oplot,[3.65,3.65],[0,100],thick=2,color=110
        loadct,33
    endif
    oplot,Fnum[f,*,j],EE[f,*,j],thick=1,color=colors[j]
    xyouts,0.17,0.9-(j*0.03),fibers[j]+':',/normal,color=colors[j],charthick=2,charsize=0.8
    fprint = sigfig(Fnum[f,89,j],4)
    xyouts,0.22,0.9-(j*0.03),'F/'+fprint+' (95%)',/normal,color=colors[j],charthick=2,charsize=0.8
endfor
device,/close_file

if (n1 gt 1) then begin
device,file='FRDblue_v_red.ps',/color
for j=0,n0-1 do begin
    if (j eq 0) then begin
        loadct,0
        plot,Fnum[0,*,j],EE[0,*,j],xrange=[xr[0],xr[1]],/xstyle,yrange=[yr[0],yr[1]],/ystyle,/nodata,$
          xtitle='Output F-ratio',ytitle='% Encircled Energy',xthick=2,ythick=2,$
          charthick=2,thick=2,title=title1
        xyouts,3.69,10,'F/#: 3.65',orientation=90,charsize=0.8,charthick=2
        xyouts,3.34,10,'F/#: 3.30',orientation=90,charsize=0.8,charthick=2
        loadct,4
        oplot,[3.3,3.3],[0,100],thick=2,color=150
        oplot,[3.65,3.65],[0,100],thick=2,color=110
    endif
    oplot,Fnum[0,*,j],EE[0,*,j],thick=1,color=60
    oplot,Fnum[n1-1,*,j],EE[n1-1,*,j],thick=1,color=150
    xyouts,0.17,0.9-(j*0.03),fibers[j]+':',/normal,color=colors[j],charthick=2,charsize=0.8
    fprint = sigfig(Fnum[5,89,j],4)
    xyouts,0.22,0.9-(j*0.03),'F/'+fprint+' (95%)',/normal,color=colors[j],charthick=2,charsize=0.8
endfor
device,/close_file
endif
set_plot,'x'

pause
while !d.window ne -1 do wdelete, !d.window

stop
END
