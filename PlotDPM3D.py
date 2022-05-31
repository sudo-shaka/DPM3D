import matplotlib.pyplot as plt
from numpy import random, ceil

def plotcell(Cell):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for tri in Cell.TriangleIndex:
        x = [Cell.Positions[i].x for i in tri];
        y = [Cell.Positions[i].y for i in tri];
        z = [Cell.Positions[i].z for i in tri];
        ax.plot(x,y,z,color='black');
        F = [abs(Cell.Forces[i].x) + abs(Cell.Forces[i].y) + abs(Cell.Forces[i].z) for i in tri];
        ax.scatter(x,y,z,c=F,cmap='coolwarm')

    #plt.show()

def plottissue(Tissue):
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(projection='3d')
    random.seed(10)
    r1 = random.rand(Tissue.NCELLS)
    r2 = random.rand(Tissue.NCELLS)
    r3 = random.rand(Tissue.NCELLS)
    for ci in range(Tissue.NCELLS):
        Cell = Tissue.Cells[ci];
        for tri in Cell.TriangleIndex:
            x = [Cell.Positions[i].x for i in tri]
            y = [Cell.Positions[i].y for i in tri]
            z = [Cell.Positions[i].z for i in tri]
            plot = True
            for i in range(3):
                if x[i] < 0:
                    x[i] += Tissue.L*(round(abs(x[i]/Tissue.L))+1)
                    plot = False
                if y[i] < 0:
                    y[i] += Tissue.L*(round(abs(y[i]/Tissue.L))+1)
                    plot = False
                if z[i] < 0:
                    z[i] += Tissue.L
                    plot = False
                if x[i] > Tissue.L:
                    x[i] -= Tissue.L * ceil((x[i]-Tissue.L)/Tissue.L)
                    plot = False
                if y[i] > Tissue.L:
                    y[i] -= Tissue.L * ceil((y[i]-Tissue.L)/Tissue.L)
                    plot = False
                if z[i] > Tissue.L:
                    z[i] -= Tissue.L * ceil((z[i]-Tissue.L)/Tissue.L)
                    plot = False

            xt = sorted(x);yt = sorted(y);zt = sorted(z)
            if plot == False:
                if (abs(xt[0] - xt[-1])**2 + abs(yt[0] - yt[-1])**2 + abs(zt[0] - zt[-1])**2) < Tissue.L:
                    plot = True
            ax.set_xlim(0,Tissue.L);
            ax.set_ylim(0,Tissue.L);
            ax.set_zlim(0,Tissue.L);
            if plot:
                ax.plot(x,y,z,color=(r1[ci],r2[ci],r3[ci]))
def plottissue2D(Tissue):
    fig = plt.figure(figsize=(10,10))
    random.seed(10)
    r1 = random.rand(Tissue.NCELLS)
    r2 = random.rand(Tissue.NCELLS)
    r3 = random.rand(Tissue.NCELLS)
    for ci in range(Tissue.NCELLS):
        Cell = Tissue.Cells[ci];
        for tri in Cell.TriangleIndex:
            x = [Cell.Positions[i].x for i in tri]
            y = [Cell.Positions[i].y for i in tri]
            plot = True
            for i in range(3):
                if x[i] < 0:
                    x[i] += Tissue.L*(round(abs(x[i]/Tissue.L))+1)
                    plot = False
                if y[i] < 0:
                    y[i] += Tissue.L*(round(abs(y[i]/Tissue.L))+1)
                    plot = False
                    plot = False
                if x[i] > Tissue.L:
                    x[i] -= Tissue.L * ceil((x[i]-Tissue.L)/Tissue.L)
                    plot = False
                if y[i] > Tissue.L:
                    y[i] -= Tissue.L * ceil((y[i]-Tissue.L)/Tissue.L)
                    plot = False

            xt = sorted(x);yt = sorted(y)
            if plot == False:
                if (abs(xt[0] - xt[-1])**2 + abs(yt[0] - yt[-1])**2) < Tissue.L:
                    plot = True
            if plot:
                plt.plot(x,y,color=(r1[ci],r2[ci],r3[ci]))
            plt.scatter(x,y,color=(r1[ci],r2[ci],r3[ci]))
            plt.axis('equal')
