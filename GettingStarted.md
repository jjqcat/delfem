# 解析プログラム全体の流れ #

有限要素法解析の手続きは以下のようになっています
  * 形状を作る
  * 形状からメッシュを切る
  * FEM離散場にメッシュを設定する
  * 方程式クラスにFEM離散場を設定
  * 方程式に物性値を設定
  * 方程式に境界条件を設定
  * 方程式を解く
  * 可視化のための描画クラスを設定

<img src='http://delfem.googlecode.com/files/fem_procedure.jpg' width='400' height='300'>



以下それぞれの段階について，詳しく説明させていただきます．<br>
<br>
<h1>形状を作る</h1>

ここではFEMのための形状をDelFEMを使って，どのように作るかについて説明します．有限要素法解析を行うためには，メッシュが切られてなければなりません．DelFEMでは内部の２次元のメッシュ生成クラスとして(Msh::Mesher2D)が用意されています．このクラスを使ってメッシュを切るためには，抽象２次元形状クラス(Cad::ICad2D)から形状クラスを派生させる必要があります．ICad2Dの派生クラスを自分で定義しても良いのですが，DelFEM内に既に用意されている(Cad::CCadObj2D)を使うと便利です．<br>
<br>
２次元形状クラス(CCadObj2D)は単純な操作：例えば，多角形を作る，点を辺やループの内部に作る，点と点を繋ぎ合わせる，点や辺を削除する．．．などで形状を作れるクラスです．ここでは代表的な操作を紹介したいと思います．<br>
<br>
<h3>多角形の追加</h3>
次の例は正方形を作る例です．２次元座標の列を(Cad::CadObj2D::AddPolygon)を使って指定すると，その座標の頂点を持つ多角形のループが生成されIDが返ります．頂点の列は時計回りでも反時計回りでも構いませんが，交差をしてはループが発生されません．ループの発生に失敗した場合は０が返ります．<br>
<pre><code>// 正方形を作る<br>
Cad::CCadObj2D cad_2d;<br>
std::vector&lt;Com::CVector2D&gt; vec_ary;<br>
vec_ary.push_back( Com::CVector2D(-0.0,-1.0) );<br>
vec_ary.push_back( Com::CVector2D( 1.0,-1.0) );<br>
vec_ary.push_back( Com::CVector2D( 1.0, 0.0) );<br>
vec_ary.push_back( Com::CVector2D(-0.0, 0.0) );<br>
unsigned int id_l = cad_2d.AddPolygon( vec_ary ); // id_lはループのID<br>
</code></pre>


<h1>メッシュを作る</h1>

(Cad::ICad2D)を継承するクラスなら(Msh::CMesher2D)を使うと簡単に２次元のメッシュを生成することができます．<br>
<br>
形状全体に一様な要素の長さでメッシュを切る場合は，コンストラクタでそれを行うことができます．<br>
<br>
<pre><code>Cad::CCadObj2D cad_2d;<br>
//<br>
// edit cad_2d here<br>
//<br>
Msh::CMesher2D mesh2d(cad_2d,0.04); // create mesh with element length=0.04<br>
</code></pre>


<h1>FEM離散場にメッシュを設定する</h1>

作成されたメッシュを有限要素法に設定します．下の例はメッシュクラスCMesher2Dのインスタンスmesh2dを有限要素法に設定している例です．CFieldWorld::AddMesh()の引数にメッシュを代入すると，のこのメッシュをFEM内で表すunsigned int型のハンドルが返ります．失敗した場合は０が返ります．<br>
<br>
<pre><code>Fem::Field::CFemField world;<br>
const unsigned int id_base = world.AddMesh( mesh2d ); // set mesh into fem field<br>
</code></pre>

FEMに設定した時に得られるハンドルからCFieldWorld::GetIDConverter()という関数を呼ぶとCad上の要素IDと要素配列IDを関連づけてくれるオブジェクトCIDConvEAMshCadが返ります．<br>
<br>
<pre><code> const CIDConvEAMshCad conv = world.GetIDConverter(id_base); // get ID converter<br>
</code></pre>

<h1>方程式クラスにFEM離散場を設定</h1>

<h1>方程式に物性値を設定</h1>

<h1>方程式に境界条件を設定</h1>

方程式に境界条件を設定します．偏微分方程式の境界条件としては固定境界条件と応力境界条件があります．ここでは固定境界条件の設定について解説します．<br>
<br>
方程式クラスからに要素配列を固定境界条件として設定します．戻り値としてunsigned int型の固定境界条件上の部分場を表すハンドルが返ります．<br>
<pre><code>unsigned int CEqnSystem::AddFixElemAry(unsigned int,CFieldWorld&amp;);<br>
unsigned int CEqnSystem::AddFixElemAry(std::vector&lt;unsigned int&gt;&amp;,CFieldWorld&amp;);<br>
</code></pre>

固定境界条件を設定したデフォルトの状態では値が０の状態で固定されています．場管理クラス(CFieldWorld)からこの部分場を表す場のクラス(CField)を取得して値を設定することで，固定された場に値を設定することができます．例えばCadの４番のEdgeに相当する要素配列に'sin(t*PI*2*0.1)'で表される単振動するy変位を与えるのは次のようなコードになります．<br>
<br>
<pre><code>unsigned int id_ea_bc1 = conv.GetIdEA_fromCad(4,Cad::EDGE); // element array at boundary<br>
unsigned int id_field_bc1 = solid.AddFixElemAry(id_bc_ea1,world);<br>
{ // set Value to boundary condition<br>
 CField&amp; bc1_field = world.GetField(id_field_bc0);<br>
 bc1_field.SetValue("sin(t*PI*2*0.1)", 1,Fem::Field::VALUE, world,true);<br>
}<br>
</code></pre>

方程式に要素配列のIDを設定する以外にも，境界条件を設定したい場の部分的な場のハンドルを取得して方程式クラスにそれを設定する方法もあります．<br>
<br>
<pre><code>unsigned int id_ea_bc = conv.GetIdEA_fromCad(1,Cad::VERTEX);<br>
unsigned int id_field_bc = world.GetPartialField(id_field,id_ea_bc); // get handle of part of field(ID:id_field)<br>
fluid.AddFixField(id_field_bc,world); // set equation the fixed boundary condition<br>
</code></pre>

Fem::Field::CField::SetValueの第５変数のboolは，この設定された値を保存するかどうかです．この値がfalseの場合はこの式が場にこの関数を呼んだ時一度だけ値が設定されます．trueの場合は設定された値は保存され，Fem::Field::CFieldWorld::FieldValueExec(double)関数が呼ばれる度に，引数のdoubleの時間がtに代入されて場の値に設定されます．<br>
<br>
<pre><code>world.FieldValueExec(cur_time); // update the value of boundary condition with current time <br>
</code></pre>

<h1>方程式を解く</h1>

全ての連立方程式クラス(Fem::Eqn::CEqnSystem...)は抽象連立方程式クラス(Fem::Eqn::CEqnSystem)を継承しています.このクラスはSolve()というメンバ関数があるので，方程式を解くためにはこの関数を呼ぶだけで構いません．<br>
<br>
<pre><code>solid.Solve(world)<br>
</code></pre>

動的なクラスでは，方程式を解くことは時間発展を計算することに他なりません．それもこのSolve()を呼ぶことで，予め次のFem::Eqn::CEqnSystem::SetTimeIntegrationParameter()という関数を呼ぶことで，方程式に定められた時間刻みだけ時間発展をさせることができます．<br>
<br>
<pre><code>solid.SetTimeIntegrationParameter(dt); // Set time step<br>
</code></pre>

<h1>可視化のための描画を設定</h1>

可視化クラスは抽象クラス，Com::CDrawerを継承しています．Com::CDrawerクラスは描画関数Draw()の他，バウンディングボックスなどを取得する関数GetBoundingBox()などが宣言されています．<br>
<br>
場を可視化するクラスはCom::CDrawerを継承する抽象クラスCom::CDrawer_Fieldを継承しています．これらはUpdate(const Fem::Field::CFieldWorld&)という関数が宣言されており，場の値を更新します．<br>
<br>
<table><thead><th>Fem::Field::View::CDrawerFaceクラス </th><th> 場の領域を描画するクラス </th></thead><tbody>
<tr><td>Fem::Field::View::CDrawerEdgeクラス </td><td> 場のメッシュの辺を描画するクラス </td></tr>
<tr><td>Fem::Field::View::CDrawerVectorクラス </td><td> ベクトル場を矢印を引いて描画するクラス </td></tr>