import DPM3D
import matplotlib.pyplot as plt
import imageio
from progressbar import progressbar
from IO import OutCellToObj as saveObj

Cell = DPM3D.Cell(x1=0,y1=0.0,z1=10,calA0=1.3,VertexRecursion=2,r0=10,Kv = 1.0, Ka = 1.0, Kb = 0.0)
Cell.Ks = 1.0;
nout = 50; dt = 0.1; nsteps = 1000;

for i in progressbar(range(nout)):
    for j in range(nsteps):
        Cell.EulerUpdate(1,dt);
        Cell.StickToSurface(0,10)
    saveObj("out.obj",Cell)
