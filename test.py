import DPM3D
def celltest():
  Kv = 0.7; Ka = 2.0; Kb = 0.0
  Cell =DPM3D.Cell(x1=0,y1=0,z1=5.0,calA0=1.05,VertexRecursion=2,r0=1.0,Kv=Kv,Ka=Ka,Kb=Kb)
  Z = [Cell.Positions[i].z for i in range(len(Cell.Positions))]
  minz = min(Z);

  for i in range(100000):
    Cell.StickToSurface(minz,0.4)
    Cell.EulerUpdate(1,0.001)

def ttest():
  Kv = 0.7; Ka = 2.0; Kb = 0.0
  Cell =DPM3D.Cell(x1=0,y1=0,z1=5.0,calA0=1.05,VertexRecursion=2,r0=1.0,Kv=Kv,Ka=Ka,Kb=Kb)
  Z = [Cell.Positions[i].z for i in range(len(Cell.Positions))]
  minz = min(Z);
  T = DPM3D.Tissue([Cell]*10,0.35)
  T.MonolayerDisperse()
  for i in range(1000):
    T.StickToSurface(minz,0.4)
    T.EulerUpdate(1,0.001)


ttest()
