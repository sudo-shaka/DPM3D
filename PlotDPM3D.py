import matplotlib.pyplot as plt
import numpy as np
from numpy import random, ceil
import gc

def PlotECM(ECM,fig):
  ax = fig.add_subplot(projection='3d')
  X = [ECM.Positions[i].x for i in range(ECM.NP)]
  Y = [ECM.Positions[i].y for i in range(ECM.NP)]
  Z = [ECM.Positions[i].z for i in range(ECM.NP)]
  ax.scatter(X,Y,Z,color='black')

def plotcell(Cell):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for tri in Cell.TriangleIndex:
        x = [Cell.Positions[i].x for i in tri];
        y = [Cell.Positions[i].y for i in tri];
        z = [Cell.Positions[i].z for i in tri];
        ax.plot(x,y,z,color='black');
        F = [abs(Cell.Forces[i].x) + abs(Cell.Forces[i].y) + abs(Cell.Forces[i].z) for i in tri];
        ax.scatter(x,y,z,c=F,cmap='coolwarm',animated=True)

    #plt.show()

def plottissue(Tissue):
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(projection='3d')
    random.seed(1)
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
                ax.plot(x,y,z,color=(r1[ci],r2[ci],r3[ci]),animated=True)


def  MapToCylinder(point, radius, scale):
    theta = point.x * scale;
    x = (radius-point.z) * np.cos(theta)
    z = (radius-point.z) * np.sin(theta)
    y = point.y

    return x,y,z


def plotVessel3D(Tissue):
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(projection='3d')
    ax.set_xlim(-Tissue.L/2,Tissue.L/2);
    ax.set_ylim(-Tissue.L/2,Tissue.L/2);
    ax.set_zlim(0,Tissue.L);
    r1 = random.rand(Tissue.NCELLS)
    r2 = random.rand(Tissue.NCELLS)
    r3 = random.rand(Tissue.NCELLS)
    for ci in range(Tissue.NCELLS):
        [x,y,z] = Tissue.GetVesselPositions(ci)
        ax.scatter(x,y,z,color=(r1[ci],r2[ci],r3[ci]))


def plotVessel(Tissue):
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(projection='3d')
    random.seed(1)
    ax.set_xlim(-Tissue.L/2,Tissue.L/2);
    ax.set_ylim(-Tissue.L/2,Tissue.L/2);
    ax.set_zlim(0,Tissue.L);
    r1 = random.rand(Tissue.NCELLS)
    r2 = random.rand(Tissue.NCELLS)
    r3 = random.rand(Tissue.NCELLS)
    radius = Tissue.L/(2*np.pi)
    scale = (2*np.pi)/Tissue.L
    for ci in range(Tissue.NCELLS):
        Cell = Tissue.Cells[ci];
        for tri in Cell.TriangleIndex:
            points = [Cell.Positions[i] for i in tri];
            for p in points:
                x,y,z = MapToCylinder(p,radius,scale)
                ax.scatter(x,y,z, color=(r1[ci],r2[ci],r3[ci]),animated=True)


def plottissue2D(Tissue):
    plt.figure(figsize=(10,10))
    random.seed(1)
    r1 = random.rand(Tissue.NCELLS)
    r2 = random.rand(Tissue.NCELLS)
    r3 = random.rand(Tissue.NCELLS)
    for ci in range(Tissue.NCELLS):
        Cell = Tissue.Cells[ci];
        for tri in Cell.TriangleIndex:
            x = [Cell.Positions[i].x for i in tri]
            y = [Cell.Positions[i].y for i in tri]
            j = [Cell.isJunction[i] for i in tri]
            a = [Cell.isFocalAdhesion[i] for i in tri]
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
                if j[i]:
                    plt.scatter(x[i],y[i], color='black',animated=True)
                elif a[i]:
                    plt.scatter(x[i],y[i], color='red',animated=True)
                else:
                    plt.scatter(x[i],y[i], color=(r1[ci],r2[ci],r2[ci]),animated=True)

            xt = sorted(x);yt = sorted(y)
            if plot == False:
                if (abs(xt[0] - xt[-1])**2 + abs(yt[0] - yt[-1])**2) < Tissue.L:
                    plot = True
            if plot:
                plt.plot(x,y,color=(r1[ci],r2[ci],r3[ci]),animated=True)
            plt.xlim([0,Tissue.L])
            plt.ylim([0,Tissue.L])
    del x,y,plot,xt,yt
    gc.collect()
