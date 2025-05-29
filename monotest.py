#/usr/bin/python3

import DPM3D
import PlotDPM3D
import matplotlib.pyplot as plt
from progressbar import progressbar

#contruct cells
C1 = DPM3D.Cell(x1=8,y1=8,z1=1.3,calA0=1.05,VertexRecursion=2,r0=1.2,Kv=2.0,Ka=0.5,Kb=0)
C2 = DPM3D.Cell(x1=0,y1=0,z1=1.5,calA0=1.05,VertexRecursion=2,r0=1.5,Kv=2.0,Ka=0.5,Kb=0)

#construct tissue
T = DPM3D.Tissue([C1,C2]*16,0.3) #array of cells and packing fraction
T.MonolayerDisperse() #disperse cells randomly in 2D (X/Y)
for ci in range(T.NCELLS): 
  T.Cells[ci].Ks = 2.0 #Set the adherance to surface stiffness (surface as z=0)
T.Kre = 10.0 #stiffness to prevent cell overlapping

#integration variables
nsteps = 500; nout = 5; dt = 0.001

print("Images will be saved to /tmp/")

#updating and saving plots
for i in progressbar(range(nout)):
  for _ in range(nsteps):
    T.UpdateShapeForces()
    T.InteractingUpdate()
    for ci in range(T.NCELLS):
      T.Cells[ci].StickToSurface(0.25,1.0)
    T.EulerUpdate(dt)
  PlotDPM3D.plottissue2D(T)
  plt.savefig("/tmp/"+str(i)+"test_out.png")
  plt.close()
