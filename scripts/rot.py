from paraview.simple import *
from math import pi,sin,cos,sqrt
import sys

Render()
SaveScreenshot('/work2/data/neuro/ims/neur_%05d.png'%(0))

rv = GetActiveView()
cam = rv.GetActiveCamera()

t = 0.
n = 1000
t0 = 0.
t1 = 2.*pi
dt = (t1 - t0)/(n - 1.)
cdt = cos(dt)
sdt = sin(dt)
i = 0
while i < n:
  sys.stderr.write('.')
  if ((i + 1) % 80) == 0:
    sys.stderr.write('\n')
  pos = list(cam.GetPosition())
  fp = rv.CenterOfRotation
  x0 = pos[0] - fp[0]
  z0 = pos[2] - fp[2]
  x1 = x0*cdt - z0*sdt
  z1 = x0*sdt + z0*cdt
  x1 += fp[0]
  z1 += fp[2]
  pos[0] = x1
  pos[2] = z1
  cam.SetPosition(pos)
  cam.SetFocalPoint(fp)
  t += dt
  #sys.stderr.write('t=%g\n'%(t*180./pi))
  ct = cos(t)
  st = sin(t)
  cam.SetViewUp([-st,0.,ct])
  Render()
  SaveScreenshot('/work2/data/neuro/ims/neur_%05d.png'%(i+1))
  i += 1


sys.stderr.write('done!')
