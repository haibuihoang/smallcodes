;Function for smooth 2d fields using cosine weighting function 

undef("smoothie")
function smoothie(f[*][*]:float,dx:float,r:float)
begin
  pi=3.14159265358979
  pio2=pi/2
  dims=dimsizes(f)
  ny=dims(0)
  nx=dims(1)
  print(dims)
  print("Smoothing field")
  fnew=f
  dij = tointeger(ceil(r/dx))
  twodij=2*dij

;"caching the weight matrix"
  weird = new((/2*dij+1,2*dij+1/),"float")
  do i=0,2*dij
    do j=0,2*dij
      rx=(i-dij)*dx
      ry=(j-dij)*dx
      rad=sqrt(rx^2+ry^2)/r
      if (rad.le.1) then
         weird(j,i) = cos(rad*pio2)^2
      else
         weird(j,i) = 0.
      end if
    end do
  end do

  do j=0,ny-1
    jj1=0     ;This is the weight matrix indices
    jj2=twodij
    j1=j-dij  ;This is the matrix in the fields
    j2=j+dij
    if (j1.lt.0) then
      jj1=-j1
      j1=0
    end if
    if (j2.gt.ny-1) then
       jj2=twodij+ny-1-j2
       j2=ny-1
    end if
    print(j+"/"+ny)
    do i=0,nx-1i
       ii1=0
       ii2=twodij
       i1=i-dij
       i2=i+dij
       if (i1.lt.0) then
         ii1=-i1
         i1=0
       end if
       if (i2.gt.nx-1) then
         ii2=twodij+nx-1-i2
         i2=nx-1
       end if

       w=sum(weird(jj1:jj2,ii1:ii2))
       ff=sum(weird(jj1:jj2,ii1:ii2)*f(j1:j2,i1:i2))
       fnew(j,i)=ff/w
  end do
  end do

  return fnew
end


