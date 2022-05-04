def OutCellToFile(filename,Cell):
  X = [Cell.Positions[i].x for i in range(Cell.NV)]
  Y = [Cell.Positions[i].y for i in range(Cell.NV)]
  Z = [Cell.Positions[i].z for i in range(Cell.NV)]
  with open(filename,'a') as out:
    out.write(str(X)[1:-1]+"\n")
    out.write(str(Y)[1:-1]+"\n")
    out.write(str(Z)[1:-1]+"\n")

