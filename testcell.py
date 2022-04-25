import DPM3D
import matplotlib.pyplot as plt
import imageio
from progressbar import progressbar
import os

Cell = DPM3D.Cell(x1=0,y1=0.0,z1=0.0,calA0=1.05,VertexRecursion=3,Kv = 1.0, Ka = 0.05, Kb = 0.05)

nout = 180; dt = 0.01; nsteps = 100;

for i in range(nout):
    for j in range(nsteps):
        Cell.EulerUpdate(1,dt);
        Cell.UpdateShapeForces();

    print(Cell.GetVolume(),Cell.GetSA())
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(projection='3d')
    ax.view_init(10, -120+(i))
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