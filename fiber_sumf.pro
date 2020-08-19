function fiber_sumf, data, maskname, aperture

;This function just returns both the total and the median value for a
;fiber value in a VP frame. The output is an array of size
;2,n_fibers. The first column is TOTAL and the second is MEDIAN

step = (aperture-1)/2.0

n1 = n_elements(data[*,0])

readcol,maskname,format='i,x',y
n3 = n_elements(y)
indexarr = y - 1.0

dataout = fltarr(2,n3)

for jj=0,n3-1 do begin ;a loop through each fiber
    index = uint(indexarr[jj])
    rowD = data[*,index-step:index+step]
    if (median(rowD) eq -666) then begin
        dataout[0,jj] = 0.0
        dataout[1,jj] = 0.0
    endif else begin
        dataout[0,jj] = total(rowD)
        dataout[1,jj] = median(rowD)
    endelse
endfor

return,dataout

stop
end
