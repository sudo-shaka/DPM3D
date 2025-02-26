import DPM3D
from PlotDPM3D import plottissue2D
import matplotlib.pyplot as plt
from progressbar import progressbar

C1 = DPM3D.Cell(x1=8,y1=8,z1=1,calA0=1.05,VertexRecursion=2,r0=1.3,Kv=0.7,Ka=2.0,Kb=0)
C2 = DPM3D.Cell(x1=0,y1=0,z1=1,calA0=1.05,VertexRecursion=2,r0=1.5,Kv=0.7,Ka=2.0,Kb=0)

T = DPM3D.Tissue([C1,C2]*12,0.35)
T.MonolayerDisperse()
for ci in range(T.NCELLS):
  T.Cells[ci].Ks = 5.0
T.Kre = 100.0
T.Kat = 2;
T.AttractionMethod = "JunctionSlip"
nsteps = 10; nout = 25; dt = 0.01
for i in progressbar(range(nout)):
  for _ in range(nsteps):
    T.UpdateShapeForces()
    T.InteractingUpdate()
    for ci in range(T.NCELLS):
      T.Cells[ci].StickToSurface(0.25,1.0)
    T.EulerUpdate(dt)
  if i == nout - 1:
    plottissue2D(T)
    plt.savefig("test_out.png")
    plt.close()
