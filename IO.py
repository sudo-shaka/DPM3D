def OutCellToObj(filename,Cell):
  X = [Cell.Positions[i].x for i in range(Cell.NV)]
  Y = [Cell.Positions[i].y for i in range(Cell.NV)]
  Z = [Cell.Positions[i].z for i in range(Cell.NV)]
  string = ["v {:.4f} {:.4f} {:.4f}\n".format(X[i], Y[i], Z[i]) for i in range(Cell.NV)]
  with open(filename,'w') as out:
    out.write("# OBJ file\n")
    out.write("o Cell\n")
    for i in range(Cell.NV):
      out.write(string[i])
    for i in range(len(Cell.TriangleIndex)):
      face = Cell.TriangleIndex[i]
      face = [f+1 for f in face]
      fstring = "f "+str(face[0])+" "+str(face[1])+" "+str(face[2]) + "\n"
      out.write(fstring)

  out.close()


def OutCellToFile(filename,Cell):
  X = [Cell.Positions[i].x for i in range(Cell.NV)]
  Y = [Cell.Positions[i].y for i in range(Cell.NV)]
  Z = [Cell.Positions[i].z for i in range(Cell.NV)]
  Fv = Cell.GetVolumeForces()
  Fa = Cell.GetAreaForces()
  Fb = Cell.GetBendingForces()
  Fv = [abs(Fv[i].x)+abs(Fv[i].y)+abs(Fv[i].z) for i in range(len(Fv))]
  Fa = [abs(Fa[i].x)+abs(Fa[i].y)+abs(Fa[i].z) for i in range(len(Fa))]
  Fb = [abs(Fb[i].x)+abs(Fb[i].y)+abs(Fb[i].z) for i in range(len(Fv))]
  with open(filename,'a') as out:
    out.write(str(X)[1:-1]+"\n")
    out.write(str(Y)[1:-1]+"\n")
    out.write(str(Z)[1:-1]+"\n")
    out.write(str(Fv)[1:-1]+"\n")
    out.write(str(Fa)[1:-1]+"\n")
    out.write(str(Fb)[1:-1]+"\n")
  out.close()
  with open(filename+'_shape','a') as f:
    f.write(str(Cell.GetSA())+','+str(Cell.GetVolume())+'\n')
  f.close()

def OutTissueToFile(filename,Tissue):
  with open(filename,'w') as out:
    out.write('\n')
  out.close()
  for ci in range(Tissue.NCELLS):
    X = [Tissue.Cells[ci].Positions[i].x for i in range(Tissue.Cells[ci].NV)]
    Y = [Tissue.Cells[ci].Positions[i].y for i in range(Tissue.Cells[ci].NV)]
    Z = [Tissue.Cells[ci].Positions[i].z for i in range(Tissue.Cells[ci].NV)]
    with open(filename,'a') as out:
      out.write(str(X)[1:-1]+"\n")
      out.write(str(Y)[1:-1]+"\n")
      out.write(str(Z)[1:-1]+"\n")
    out.close()
  with open(filename,'a') as out:
    out.write('\n')
  out.close()
