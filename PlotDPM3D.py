import matplotlib.pyplot as plt
import numpy as np
from numpy import random

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
            x = np.mod([Cell.Positions[i].x for i in tri],Tissue.L)
            y = np.mod([Cell.Positions[i].y for i in tri],Tissue.L)
            z = np.mod([Cell.Positions[i].z for i in tri],Tissue.L)
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

def plotjunction2D(Tissue):
    plt.figure(figsize=(10,10))
    for ci in range(Tissue.NCELLS):
        cell = Tissue.Cells[ci]
        for tri in cell.TriangleIndex:
            isJunction = [cell.isJunction[i] for i in tri]
            for i in range(3):
                if isJunction[i]:
                    f = abs(cell.Forces[tri[i]].x) + abs(cell.Forces[tri[i]].y) + abs(cell.Forces[tri[i]].z)
                    plt.scatter(np.mod(cell.Positions[tri[i]].x,Tissue.L),np.mod(cell.Positions[tri[i]].y,Tissue.L),c=f,cmap = 'coolwarm')


def plottissue2D(Tissue):
    plt.figure(figsize=(10,10))
    random.seed(1)
    r1 = random.rand(Tissue.NCELLS)
    r2 = random.rand(Tissue.NCELLS)
    r3 = random.rand(Tissue.NCELLS)
    for ci in range(Tissue.NCELLS):
        pos = Tissue.Cells[ci].GetPositions()
        if np.isnan(pos).any():
            print("NaN values found!")
        x,y = pos[0], pos[1]
        if(Tissue.PBC):
            x,y = np.mod(pos[0],Tissue.L), np.mod(pos[1],Tissue.L)
        for tri in Tissue.Cells[ci].TriangleIndex:
            vx,vy = [x[i] for i in tri], [y[i] for i in tri]
            if max(vx) - min(vx) <= Tissue.L/2 and max(vy) - min(vy) <= Tissue.L/2:
                plt.plot(vx,vy, color = (r1[ci],r2[ci],r3[ci]))
            plt.scatter(x,y,s=3, color = (r1[ci],r2[ci],r3[ci]))
    if(Tissue.PBC):
        plt.xlim([0,Tissue.L])
        plt.ylim([0,Tissue.L])
    else:
        plt.axis('equal')
