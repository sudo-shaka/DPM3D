def OutCellToFile(filename,Cell):
  X = [Cell.Positions[i].x for i in range(Cell.NV)]
  Y = [Cell.Positions[i].y for i in range(Cell.NV)]
  Z = [Cell.Positions[i].z for i in range(Cell.NV)]
  with open(filename,'a') as out:
    out.write(str(X)[1:-1]+"\n")
    out.write(str(Y)[1:-1]+"\n")
    out.write(str(Z)[1:-1]+"\n")

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
