# ２次元メッシュを切る #

## 要素長さを指定してメッシュを切る ##
```
Cad::CCadObj2D cad_2d;
Msh::CMesher2D(cad_2d,0.1) // cut mesh with elem length=0.1
```


## 出来るだけ少ない要素の数でメッシュを切る ##
```
Msh::CMesher2D mesh(cad_2d);
```

## メッシュサイズを大まかに指定してメッシュを切る ##



## 指定した領域にメッシュを切る ##



# FEMとの接続 #

```
unsigned int id_base = world.AddMesh(mesh);
```

# 切り込みをいれたメッシュの生成 #