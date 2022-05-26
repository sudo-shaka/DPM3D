import DPM3D
from PlotDPM3D import plottissue
from matplotlib import pyplot as plt
from progressbar import progressbar
import imageio

Cell = DPM3D.Cell(x1=0,y1=0,z1=0,calA0=1.0, VertexRecursion=2,r0=1,Kv=1.0,Ka=1.0, Kb = 0.05)

CellList = [Cell]*5;

T = DPM3D.Tissue(Cells=CellList, phi0=0.7)
T.TissueDisperse()
T.Kc = 1.0;

with imageio.get_writer('/tmp/out.gif',mode='I') as writer:
  for i in progressbar(range(90)):
    T.EulerUpdate(nsteps=10,dt=0.005)
    plottissue(T)
    plt.savefig('/tmp/'+str(i)+'.png', bbox_inches='tight')
    plt.close()
