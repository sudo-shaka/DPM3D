import matplotlib.pyplot as plt

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
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for ci in range(Tissue.NCELLS):
        Cell = Tissue.Cells[ci];
        for tri in Cell.TriangleIndex:
            x = [Cell.Positions[i].x for i in tri];
            y = [Cell.Positions[i].y for i in tri];
            z = [Cell.Positions[i].z for i in tri];
            ax.plot(x,y,z,color='black');
            F = [abs(Cell.Forces[i].x) + abs(Cell.Forces[i].y) + abs(Cell.Forces[i].z) for i in tri];
            ax.scatter(x,y,z,c=F,cmap='coolwarm')