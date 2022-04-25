import DPM3D
from PlotDPM3D import plotcell
import matplotlib.pyplot as plt
import imageio
from progressbar import progressbar
import os
#os.system('find /tmp/ -maxdepth 1 -iname "*.png" -delete')

Cell = DPM3D.Cell(x1=0,y1=0.0,z1=0.0,calA0=1.05,VertexRecursion=3,Kv = 1.0, Ka = 0.05, Kb = 0.05)

#Cell.VertexFIRE(1.0,1e-4,10000,0.01)

positions = Cell.Positions;
Z = [positions[i].z for i in range(len(positions))]

minz = min(Z);
minz *= 0.8;
Cell.nSurfacePoints = 500;
Cell.SurfaceInit(minz);

nout = 180; dt = 0.01; nsteps = 1000;

for i in range(nout):
    for j in range(nsteps):
        Cell.EulerUpdate(1,dt);
        Cell.UpdateShapeForces();
        if j % 100 ==0:
            Cell.StickToSurface(minz,0.1)
            Cell.SurfaceInit(minz)
        else:
            Cell.Crawling();
        #Cell.StickToSurface(minz,0.1)
        #Cell.Crawling();
    #print(Cell.GetCalA0(),Cell.GetVolume(),Cell.GetSA())
    #Cell.SurfaceInit(minz)
    print(Cell.GetVolume(),Cell.GetSA())
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(projection='3d')
    ax.set_ylim(-1.5,1.5)
    ax.set_zlim(minz,minz+2.25)
    ax.set_xlim(-1.5,1.5)
    ax.view_init(10, -120+(i))
    #Xs = [Cell.SurfacePositions[n].x for n in range(len(Cell.SurfacePositions))]
    #Ys = [Cell.SurfacePositions[n].y for n in range(len(Cell.SurfacePositions))]
    #ax.scatter(Xs,Ys,[minz]*len(Cell.SurfacePositions),color='grey',alpha = 0.2)
    for tri in Cell.TriangleIndex:
        x = [Cell.Positions[j].x for j in tri];
        y = [Cell.Positions[j].y for j in tri];
        z = [Cell.Positions[j].z for j in tri];
        x.append(x[0]); y.append(y[0]); z.append(z[0]);
        ax.plot(x,y,z,color='black');
        x.pop(); y.pop(); z.pop();
        F = [abs(Cell.Forces[i].x) + abs(Cell.Forces[i].y) + abs(Cell.Forces[i].z) for i in tri];
        ax.scatter(x,y,z,c=F,cmap='coolwarm')
    plt.xlabel('x'); plt.ylabel('y')
    plt.savefig('/tmp/'+str(i)+'.png')
    plt.close()

with imageio.get_writer('/tmp/out.gif',mode='I') as writer:
    for i in range(nout):
        filename = '/tmp/'+str(i)+'.png'
        image = imageio.imread(filename)
        writer.append_data(image)
    #return the monolayer
    print('Image is saved to "/tmp/out.gif"')