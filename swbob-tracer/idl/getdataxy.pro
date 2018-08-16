pro getdataxy,nlon,dir,job,fld,tstr,u,swend=swend

nlat=nlon/2

dummy=fltarr(long(nlon)*nlat)
u=fltarr(nlon,nlat)

datadir='~/dataloc/bb'
;datadir='/scratch/rks/bobdata'

if keyword_set(swend) then begin
openr,1,datadir+'/'+dir+'/'+job+'/'+fld+'.'+tstr,/swap_endian,/f77_unformatted
endif else begin
openr,1,datadir+'/'+dir+'/'+job+'/'+fld+'.'+tstr   ;,/f77_unformatted
endelse

readu,1,dummy
u(*,*)=reform(dummy,nlon,nlat)

close,1

print,'got data: '+fld+'.'+tstr

end
