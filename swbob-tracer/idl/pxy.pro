; test different spherical projections

dir='swtracer'
res='T170'

job='zz-tracer.c0000.'+res

it=20

fld='q'

tstr=string(it,format='(i5.5)')
psstem=fld+tstr+'_'+res
psfile=psstem+'.ps'

if res eq 'T2730' then nlon=8192
if res eq 'T1365' then nlon=4096
if res eq 'T682' then nlon=2048
if res eq 'T341' then nlon=1024
if res eq 'T170' then nlon=512
if res eq 'T85' then nlon=256
if res eq 'T42' then nlon=128
nlat=nlon/2

xs=findgen(nlon)/nlon*360
tmp=fltarr(3,nlat/2)
openr,1,'../grids/GRID.'+res & readf,1,tmp & close,1
lattmp=tmp(1,*)
lats=[!pi/2-rotate(lattmp,3),-!pi/2+rotate(lattmp,1)]
ys=lats*180/!pi

omega=7.292e-5
fcor=2*omega*sin(lats)

;--------------------------------------------------------------------
set_plot,'ps'
!p.font=0
device,filename=psfile
device,set_font='Times-Italic'   ; assigns times-italic to use font !20
device,xsize=12,ysize=10,xoffset=1,yoffset=3,/times,/color
!p.multi=0 & !p.noerase=1
!x.style=1 & !x.minor=-1 & !y.style=1 & !y.minor=-1 & !x.title='' & !y.title=''
!p.charsize=1 & !x.charsize=.8 & !y.charsize=.8
!p.thick=1 & !x.thick=1 & !y.thick=1
erase
np=1 & ncol=1. & nrow=np/ncol & sep=0.05
;--------------------------------------------------------------------

!x.range=[0,360] & !y.range=[0,90]
!y.ticks=0 ;& !y.minor=4 & !y.tickv=[10,30,50,70]

;!x.title='longitude' & !y.title='latitude'
!p.title=''


letter=['(a) ','(b) ','(c) ','(d) ','(e) ','(f) ']

xs=findgen(nlon+1)/nlon*360

if fld eq 'h' then begin
    getdataxy,nlon,dir,job,fld,tstr,qxy
    print,'min/max '+fld+' = ',min(qxy),max(qxy)
    qxy=shift(qxy,[nlon/2,0])
    xs=xs-xs(nlon/2)
    !x.range=[-50,50] & !y.range=[20,70]
    qxy0=qxy(0,*)
    for i=0,nlat-1 do qxy(*,i)=qxy(*,i)-qxy0(i)
;    gravity=9.80616
;    qxy=qxy/gravity     ;;; hxy already output in meters (see bldsf.F)
    zlevs=(findgen(10)+1)*100
endif


getdataxy,nlon,dir,job,fld,tstr,qxy

!x.range=[0,360] & !y.range=[-90,90]
!x.ticks=6 & !y.ticks=4

if fld eq 'q' then begin
   loadct,3
   nc=64
   cols = 255*(nc-abs(nc-2*indgen(nc+1)))/nc
   qlim=2.0
   qlevs=qlim*(2*findgen(nc+1)-nc)/nc
endif

if fld eq 's' then begin
   loadct,1
   nc=64
   cols = 255*(nc-indgen(nc+1))/nc
   qlim=1.0
   qlevs=qlim*(2*findgen(nc+1)-nc/4)/nc
endif

iq=where(qxy le -qlim)
if min(iq) ge 0 then qxy(iq)=-qlim

qxye=fltarr(nlon+1,nlat)
qxye(0:nlon-1,*)=qxy & qxye(nlon,*)=qxy(0,*)
qxy=qxye

ip=0

col=ip mod ncol & row=ip/fix(ncol)
!p.position=[col/ncol,(nrow-1-row)/nrow,(col+1)/ncol,(nrow-row)/nrow] $
		 +[1.,ncol/nrow,-1.,-ncol/nrow]*sep

contour,qxy,xs,ys,levels=qlevs,c_labels=0,c_color=cols,/cell_fill

device,/close


spawn, 'pstopnm -portrait -xborder=0 -yborder=0 '+psfile
spawn, 'pnmtopng '+psstem+'001.ppm > '+psstem+'.png'
spawn, 'rm '+psfile
spawn, 'rm '+psstem+'001.ppm'
;spawn, 'display '+psstem+'.png'+' &'
;spawn,'gv '+psfile+' &'

end
