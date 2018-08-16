pro getdataxyres,nlon,job,field,it,u

nlat=nlon/2

dummy=fltarr(long(nlon)*nlat)
u=fltarr(nlon,nlat)

openr,1,'../jobs/'+job+'/'+job+'.'+field+'.'+string(it,format='(i4.4)'),/f77_unformatted

readu,1,dummy
u(*,*)=reform(dummy,nlon,nlat)

close,1

print,'got data: '+field

end

; update to following format:

; pro getdataxyz,nlon,nlevp,dir,job,field,tstr,u

; nlat=nlon/2

; nlev=nlevp
; if field eq 'w' then nlev=nlevp+1
; if field eq 'p' then nlev=nlevp-1

; dummy=fltarr(long(nlon)*nlat*nlev)
; u=fltarr(nlon,nlat,nlev)

; ;openr,1,'/home/flare-s2/scott/bobdata/'+dir+'/'+job+'/'+job+'.'+field+'.'+tstr,/f77_unformatted,/swap_endian
; openr,1,'/home/flare-s2/scott/bobdata/'+dir+'/'+job+'/'+job+'.'+field+'.'+tstr,/f77_unformatted

; readu,1,dummy
; u(*,*,*)=reform(dummy,nlon,nlat,nlev)

; close,1

; print,'got data: '+field+'.'+tstr

; end
