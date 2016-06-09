file='evol.ctl'
endT=1138
dt=10
R=1000
Dr=100

'reinit'
'open 'file

'set display color white'
'set mproj off'
t=1
i=1
while (t<endT)
   'set t ' t
   'c'
   'set xaxis 0 1000 100'
   'set rbrange -2  22';'set cint 2';'set cthick 5'
   'd v'
   h=int((t-1)/60)
   m=t-h*60-1
   
 
   'draw title Tangential wind at ' h ':'m
 
 
   if (i<10)
   'printim  img/v00' i '.png'
   else 
      if (i<100)
      'printim  img/v0' i '.png'
      else
      'printim  img/v' i '.png'
      endif
   endif       
   
   
   t=t+dt
   i=i+1
   say 't=' t ', i=' i
   
endwhile


'close 1'

function int(stuff)
  res = ''
  i = 1
  c = substr(stuff,i,1)
  while (c!='' & ('x'%c)!='x.') 
    res = res%c
    i = i + 1
    c = substr(stuff,i,1)
  endwhile
return res
  