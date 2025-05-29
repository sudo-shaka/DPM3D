import DPM3D
import PlotDPM3D
import matplotlib.pyplot as plt
from progressbar import progressbar

C1 = DPM3D.Cell(x1=5,y1=4,z1=1.5,calA0=1.05,VertexRecursion=2,r0=1.5,Kv=0.0,Ka=0.0,Kb=0)
C2 = DPM3D.Cell(x1=2.5,y1=4,z1=1.5,calA0=1.05,VertexRecursion=2,r0=1.5,Kv=0.0,Ka=0.0,Kb=0)

T = DPM3D.Tissue([C1,C2],0.35)
for ci in range(T.NCELLS):
  T.Cells[ci].Ks = 0.0
  T.Cells[ci].l0 *= 3
T.Kre = 1.0
T.Kat = 2.0;

#T.AttractionMethod = "JunctionCatch"
#T.AttractionMethod = "JunctionSlip"
T.AttractionMethod = "General"
nsteps = 100; nout = 25; dt = 0.01
for i in progressbar(range(nout)):
  PlotDPM3D.plottissue2D(T)
  plt.savefig("/tmp/"+str(i)+"test_out.png")
  plt.close()
  for _ in range(nsteps):
    T.UpdateShapeForces()
    T.InteractingUpdate()
    for ci in range(T.NCELLS):
      T.Cells[ci].StickToSurface(0.25,1.0)
    T.EulerUpdate(dt)
