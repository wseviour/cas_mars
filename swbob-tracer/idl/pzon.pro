dir='swjprad'
res='T170'

job='taurad1_m0_a00010_ld1000_r0.zz.'+res
job='j3_tau1_a00010_ld1000_r0.zzj.'+res
;job='j3_tau.1_a00001_ld1000_r0.zzj.'+res
psfile='test.ps'

nt=1000

umax=0.05

dimn=1

iplan=1
day=[9.9,10.7,18]*3600
rad=[7.1,6.0,2.5]*1e7
ufac=rad(iplan)/day(iplan)


tsubfac=10
ysubfac=2

if res eq 'T341' then nlon=1024
if res eq 'T170' then nlon=512
if res eq 'T85' then nlon=256
if res eq 'T42' then nlon=128
nlat=nlon/2

xs=findgen(nlon)/nlon*360
;ys=(1-findgen(nlat)/(nlat-1)*2)*90
tmp=fltarr(3,nlat/2)
openr,1,'../grids/GRID.'+res & readf,1,tmp & close,1
lattmp=tmp(1,*)
;lats=[rotate(lattmp,1),-rotate(lattmp,3)]     ; corrected 4/20/06, lattmp=colat
lats=[!pi/2-rotate(lattmp,3),-!pi/2+rotate(lattmp,1)]
ys=lats*180/!pi
ts=findgen(nt)
fcor=4*!pi*sin(lats)

;--------------------------------------------------------------------
set_plot,'ps'
!p.font=0
device,filename=psfile
device,set_font='Times-Italic'   ; assigns times-italic to use font !20
device,xsize=18,ysize=10,xoffset=1.5,yoffset=3,/times,/color
!p.multi=0 & !p.noerase=1
!x.style=1 & !x.minor=-1 & !y.style=1 & !y.minor=-1 & !x.title='' & !y.title=''
!p.charsize=1 & !x.charsize=.6 & !y.charsize=.6
!p.thick=1 & !x.thick=1 & !y.thick=1
erase
np=1 & ncol=1. & nrow=np/ncol & gap=0.07
;--------------------------------------------------------------------

rgb=fltarr(3,256) & tmp=fltarr(3)
openr,1,'~/cs/graphics/mycolmap'
for i=0,255 do begin
    readf,1,tmp & rgb(*,i)=tmp
endfor
close,1
tvlct,rgb(0,*),rgb(1,*),rgb(2,*)
tvlct,0,0,0
nc=40
cols = (255*indgen(nc+1))/nc

hue=0.5
red=fltarr(256) & grn=red & blu=red
ramp0=2*findgen(128)
ramp1=2*(127-findgen(128))
grn(0:127)=ramp0 & grn(128:255)=ramp1
blu(0:127)=ramp0+hue*ramp1 & blu(128:255)=ramp1
red(0:127)=ramp0 & red(128:255)=ramp1+hue*ramp0
tvlct,red,grn,blu
;tvlct,0,0,0
nc=40
cols = (255*indgen(nc+1))/nc

!p.title=''

letter=['(a) ','(b) ','(c) ','(d) ','(e) ','(f) ']

;  u  ---------------------------------------------------

getdata,nlat,nt,dir,job,'u',zon

tse=findgen(nt+1)
zone=fltarr(nlat,nt+1) 
zone(*,0)=0 & zone(*,1:nt)=zon
zonty=transpose(zone)

amax=max([abs(min(zon)),abs(max(zon))]) & clim=amax
zlevs=clim*(2*findgen(nc+1)-nc)/nc & zlevs(0)=-amax & zlevs(nc)=amax

!x.title='!20t!7'  & !x.range=[0,nt] & !x.ticks=1
!y.title='!9f!7'   & !y.range=[-90,90] & !y.ticks=6 & !y.style=1


!p.position=[0.05,0.05,0.70,0.95] + [0,1,0,-1]*gap

contour,zonty,tse,ys,levels=zlevs,c_color=cols,/cell_fill
xyouts,nt/2,95,alignment=0.5,'!20u!7(!9f!7,!20t!7)',size=0.9

;  ue  ---------------------------------------------------

!y.title=''  & !y.style=5
!x.title=''  &

if dimn eq 0 then !x.range=[-1,1]*umax  & !x.ticks=2
if dimn eq 1 then !x.range=[-1,1]*umax*ufac  & !x.ticks=2

!p.position=[0.78,0.05,0.98,.95] + [0,1,0,-1]*gap

if dimn eq 0 then plot,zon(*,nt-1),ys
if dimn eq 1 then plot,zon(*,nt-1)*ufac,ys,/noclip



;label='!20L!DD!N!7 = '+lstrs(il)
;xpos=0 & ypos=.98
;xyouts,xpos,ypos,label,size=1.0,/normal

device,/close

end
