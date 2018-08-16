pro getdata,nlat,nt,dir,job,field,u,swend=swend

dummy=fltarr(nlat)
u=fltarr(nlat,nt)

datadir='/home/flare-s2/scott/bobdata'
datadir='~/bobdata'

if keyword_set(swend) then begin
openr,1,datadir+'/'+dir+'/'+job+'/zavg.'+field,/swap_endian,/f77_unformatted
endif else begin
openr,1,datadir+'/'+dir+'/'+job+'/zavg.'+field,/f77_unformatted
endelse

for it=0,nt-1 do begin

readu,1,dummy
u(*,it)=reform(dummy,nlat)

endfor

close,1

end
