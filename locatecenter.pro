function locatecenter, data
on_error,2

cutoffp = 0.20 ;the percentage of the peak values to consider the radius for centering.

;******************************************************************
plottingF = 'on'
;plottingF = 'off'
;******************************************************************


nn0 = n_elements(data[*,0])
nn1 = n_elements(data[0,*])
nn3 = nn0 * nn1


x_axis = indgen(nn0,/long)
ellipse = fltarr(2,120)

igood  = 0
oneimage = data
stdimage = oneimage[bsort(oneimage)]
cut = cutoffp * mean(stdimage[nn3-1000:nn3-10])
print,'The cutoff radius value is '+strn(cut)
index = where(oneimage GE cut)
t2 = total(oneimage(index))
t3 = total(oneimage)
percent = t2/t3
ellipse[*,*] = Fit_Ellipse(index, XSize=nn0, YSize=nn1, CENTER=c1, AXES=c2)
center = [round(c1),round(c2)]
if (plottingF eq 'on') then begin
    loadct,37
    window,0,retain=2,xsize=nn0/winsize,ysize=nn1/winsize
    TVImage, BytScl(oneimage, Top=!D.Table_Size-3)
    plots,ellipse[*,*]/winsize,/Device,Color=FSC_Color('yellow'),thick=2
endif

return,center

stop
end
