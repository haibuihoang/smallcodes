'reinit'

'open uv.ctl'

'set gxout stream'
'set ccolor 15'
'd u;v'

*determin plotting area
'q gxinfo'
_x1=subwrd(result,14)
_x2=subwrd(result,16)
_y1=subwrd(result,20)
_y2=subwrd(result,22)
_w=_x2-_x1
_h=_y2-_y1
_rxmax=5000*100  ;*meters
_rymax=5000*100  ;*meters

res=read('adv.txt')
numrec=subwrd(res,2)
say "num:" numrec
rec=1
c=1
while(rec<=numrec)
  res=read('adv.txt')
  px=subwrd(res,2);py=subwrd(res,3)
  ;*say rec': 'px','py
  ;*calcule coordinate
  x = _x1+_w*px/_rxmax
  y = _y1+_h*py/_rymax
  ;*say x':'y
  '!sleep 0.1'
  'set line 'c
  'draw mark 2 'x' 'y' 0.1'
  c=c+1
  if (c>16)
    c=1
  endif
  rec=rec+1
endwhile


