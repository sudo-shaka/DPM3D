import DPM3D
import PlotDPM3D
import matplotlib.pyplot as plt
from progressbar import progressbar

C1 = DPM3D.Cell(x1=8,y1=8,z1=1,calA0=1.05,VertexRecursion=3,r0=1.3,Kv=0.7,Ka=1.5,Kb=0)
C2 = DPM3D.Cell(x1=0,y1=0,z1=1,calA0=1.05,VertexRecursion=3,r0=1.5,Kv=0.7,Ka=1.5,Kb=0)

T = DPM3D.Tissue([C1,C2]*12,0.35)
T.MonolayerDisperse()
for ci in range(T.NCELLS):
  T.Cells[ci].Ks = 7.0
T.Kre = 50
T.Kat = 5;
T.AttractionMethod = "JunctionSlip"
nsteps = 50; nout = 25; dt = 0.001
for i in progressbar(range(nout)):
  for _ in range(nsteps):
    T.UpdateShapeForces()
    T.InteractingUpdate()
    for ci in range(T.NCELLS):
      T.Cells[ci].StickToSurface(0.25,1.0)
    T.EulerUpdate(dt)
    #PlotDPM3D.plotVessel3D(T)
  if i == nout - 1 or i == 1:
    PlotDPM3D.plottissue2D(T)
    plt.savefig("/tmp/"+str(i)+"test_out.png")
    plt.close()
