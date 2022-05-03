import matplotlib.pyplot as plt

def plotcell(Cell):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for tri in Cell.TriangleIndex:
        x = [Cell.Positions[i].x for i in tri];
        y = [Cell.Positions[i].y for i in tri];
        z = [Cell.Positions[i].z for i in tri];
        ax.plot(x,y,z);
        F = [abs(Cell.Forces[i].x) + abs(Cell.Forces[i].y) + abs(Cell.Forces[i].z) for i in tri];
        ax.scatter(x,y,z,c=F,cmap='coolwarm')

    plt.show()