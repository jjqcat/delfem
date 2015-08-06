### How to get total number of nodes in mesh ? ###

Msh::CMesher2D have a function to get a array of node coordinates:

```
const std::vector<Com::CVector2D>& Msh::CMesher2D::GetVectorAry() const
```

and please get number of points as

```
CMesher2D mesh;
unsigned int nnode = mesh.GetVectorAry().size(); // number of nodes
```


### How to get valus from 3D deformation field? ###

get value of each element in element array (ID=id\_ea) in field (ID=id\_field\_val) will be:

```
const CField& field_val = world.GetField(id_field_val);
const CElemAry::CElemSeg& es_c = field_val.GetElemSeg(id_ea,CORNER,true,world);
const CNodeAry::CNodeSeg& ns_v = field_val.GetNodeSeg(CORNER,true,world,VALUE);
const CNodeAry::CNodeSeg& ns_c = field_val.GetNodeSeg(CORNER,false,world,VALUE);
for(unsigned int ielem=0;ielem<ea.Size();ielem++){
  unsigned int noes[4]; es_c.GetNodes(ielem,noes); // get node index for each element
  double coords[4][3];
  double disp[4][3];
  for(unsigned int ino=0;ino<4;ino++){
    ns_c.GetValue(noes[ino],coords[ino]); // get undeformed coord for each element nodes
    ns_v.GetValue(noes[ino],disp[ino]); // get displacement for each element nodes
  }
}
```