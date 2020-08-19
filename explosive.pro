PRO explosive, filelist

readcol,filelist,format='a',list

n0 = n_elements(list)

test = readfits(list[0],/silent)
n1 = float(n_elements(test[*,0]))
n2 = float(n_elements(test[0,*]))

values = fltarr(2,n0-1)

window,0,retain=2,xsize=n1/5.0,ysize=n2/5.0
device,decomposed = 0
loadct, 37

basem = readfits('base.fit',/silent)
for j=0,n0-1 do begin
    one = readfits(list[j],/silent)
;    two = readfits(list[j+1],/silent)
    three = float(one) / float(basem)
    tvimage,bytscl(three,top=!d.table_size-3)
;    tvscale,three
;    std = three[bsort(three)]
;    i = where(std lt 0.95 or std gt 1.05,count)

;    if (count gt 0) then begin
;        diff = total(abs(std[i] - 1.0))
        diff = total(abs(three - 1.0))
        values[0,j] = stddev(three)
        values[1,j] = diff
;        plot,std[i],psym=3
;        xyouts,0.2,0.2,strn(diff),/normal,charsize=2.0,charthick=2
        print,strn(median(three))+' '+strn(stddev(three))+' '+strn(diff)
;        wset,0
;    endif
    pause
endfor
window,2,retain=2
device,decomposed=0
loadct,0

plot,values[1,*],values[0,*],psym=2,/nodata
loadct,4
for j=0,3 do plots,values[1,j],values[0,j],psym=2,color=150
for j=4,8 do plots,values[1,j],values[0,j],psym=2,color=255
    

stop
END
