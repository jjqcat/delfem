# Table of contents #



---

# Overview #

DelFEM is a project aiming to provide handy finite element analysis environment for non-expert users. This project include 2-dimensional modeling, meshing and finite element analysis(FEA) components. The FEA supports various type of analysis including fluid, solid ,thermal and their coupling. Users can run FEA simulation through simple object oriented C++ problem description.


---

# Wiki #

  * Getting Started : [GettingStarted](GettingStarted.md)
  * Frequently Asked Question : [FAQ](FAQ.md)
  * Architecture of DelFEM : [Architecture](Architecture.md)


---

# Features #

  * Seamless integration of modeling, meshing and FEM
  * fast computation
  * flexible FEM that supports many PDEs
  * portability

|Modeler| 2D multi-region non-manifold|
|:------|:----------------------------|
|Mesher | 2D Triangle element, 3D Tetrahedral element(extruded from Triangle)|
|FEM (Scalar)|Poisson,Diffusion,Advection-Diffusion|
|FEM (Solid)|Linear Solid,St.Venant-Kirchhoff Material,incompressive hyperelastic material|
|FEM (Fluid)|Stokes Fluid,Navier-Stokes Fluid|
|FEM (Complex)|Helmholtz                    |
|Rigid Body|3D rigid body,６ types of constraint model, coupling analysis of rigid and elastic material|
|linear solver|CG,ILU(0) preconditioned CG,BiCGSTAB,ILU(0) preconditioned BiCGSTAB|
|Visualizatin|OpenGL                       |
|GUI    | gui with Qt is included in demo|



---

# Gallerly #

<img src='http://delfem.googlecode.com/files/solid3d.jpg' width='200' height='150'>  <img src='http://delfem.googlecode.com/files/fluid2d.jpg' width='300' height='150'>

<img src='http://delfem.googlecode.com/files/fluid2d_thinwall.jpg' width='260' height='150'>  <img src='http://delfem.googlecode.com/files/solid2d_crack.jpg' width='300' height='150'>

<img src='http://delfem.googlecode.com/files/helmholtz2d.jpg' width='150' height='150'>  <img src='http://delfem.googlecode.com/files/dkt_shell.jpg' width='200' height='150'>
<img src='http://delfem.googlecode.com/files/stvk-rigid.jpg' width='200' height='150'>

<img src='http://delfem.googlecode.com/files/qt_demo_fluid2d.png ' width='400' height='200'>

<hr />
<h1>Videos</h1>

<a href='http://www.youtube.com/watch?feature=player_embedded&v=KEou8W4z72g' target='_blank'><img src='http://img.youtube.com/vi/KEou8W4z72g/0.jpg' width='240' height=180 /></a>  <a href='http://www.youtube.com/watch?feature=player_embedded&v=HCiAH62gYfE' target='_blank'><img src='http://img.youtube.com/vi/HCiAH62gYfE/0.jpg' width='240' height=180 /></a><br>
<br>
<a href='http://www.youtube.com/watch?feature=player_embedded&v=iIPqyikaqJQ' target='_blank'><img src='http://img.youtube.com/vi/iIPqyikaqJQ/0.jpg' width='240' height=180 /></a>  <a href='http://www.youtube.com/watch?feature=player_embedded&v=WIfYP90L7jM' target='_blank'><img src='http://img.youtube.com/vi/WIfYP90L7jM/0.jpg' width='240' height=180 /></a><br>
<br>
more movies are in <a href='http://www.youtube.com/user/umet55#g/u'>http://www.youtube.com/user/umet55#g/u</a>

<hr />
<h1>Coding Example(Non static analysis of St.Venant-Kirchhoff material)</h1>

Header<br>
<pre><code>#include "delfem/cad_obj2d.h"<br>
#include "delfem/mesher2d.h"<br>
<br>
#include "delfem/field.h"<br>
#include "delfem/field_world.h"<br>
#include "delfem/drawer_field_face.h"<br>
#include "delfem/drawer_field_edge.h"<br>
<br>
#include "delfem/eqnsys_solid.h"<br>
</code></pre>

Problem Setting<br>
<pre><code> Cad::CCadObj2D cad_2d;<br>
 { // set model shape<br>
  std::vector&lt;Com::CVector2D&gt; vec_ary;<br>
  vec_ary.push_back( Com::CVector2D(0.0,0.0) );<br>
  vec_ary.push_back( Com::CVector2D(5.0,0.0) );<br>
  vec_ary.push_back( Com::CVector2D(5.0,1.0) );<br>
  vec_ary.push_back( Com::CVector2D(0.0,1.0) );<br>
  cad_2d.AddPolygon( vec_ary );<br>
 }<br>
 Msh::CMesher2D mesh2d(cad_2d,0.1); // build mesh with elem lengh = 0.1 from shape<br>
 Fem::Field::CFemField world;<br>
 const unsigned int id_base = world.AddMesh( mesh2d ); // set mesh into fem field<br>
 const CIDConvEAMshCad conv = world.GetIDConverter(id_base); // get ID converter<br>
<br>
 // set field to equation<br>
 Fem::Eqn::CEqnSystem_Solid2d solid;<br>
 solid.UpdateDomain_Field(id_base, world);<br>
 solid.SetSaveStiffMat(false);	<br>
 solid.SetStationary(false);<br>
 // set material parameter<br>
 solid.SetYoungPoisson(10.0,0.3,true);<br>
 solid.SetGeometricalNonlinear(true);<br>
 solid.SetGravitation(0.0,0.0);<br>
 solid.SetTimeIntegrationParameter(dt,0.7);<br>
<br>
 unsigned int id_field_bc0 = solid.AddFixElemAry(conv.GetIdEA_fromCad(2,Cad::EDGE),world);<br>
 unsigned int id_field_bc1 = solid.AddFixElemAry(conv.GetIdEA_fromCad(4,Cad::EDGE),world);<br>
 FieldValueSetter fvs(id_field_bc0,world);<br>
 world.SetMathExp("sin(t*PI*2*0.1)", 1,Fem::Field::VALUE, world);<br>
<br>
 // set drawer array<br>
 std::vector&lt;CDrawer*&gt; drawer_ary;<br>
 id_field_disp = solid.GetIdField_Disp();<br>
 drawer_ary.PushBack( new View::CDrawerFace(id_field_disp,false,world) );<br>
 drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,false,world) );<br>
 drawer_ary.PushBack( new View::CDrawerEdge(id_field_disp,true ,world) );<br>
</code></pre>

OpenGL drawing function<br>
<pre><code> cur_time += dt;　// update current time<br>
 fvs.Execute(cur_time,world); // update the value of boundary condition<br>
 solid.Solve(world);　// Solve FEM<br>
 drawer_ary.Update(world);　// update drawer<br>
 drawer_ary.Draw();　// OpenGL drawing<br>
</code></pre>


<hr />
<h1>ライセンス形態</h1>

なるべく沢山の人に有限要素法を"理解して"使ってもらることを目指して，オープンソフトウェアでの開発を行っています．<br>
<br>
<b>"DelFEMはLGPLライセンスversion3に基づいて配布しています"</b>

<h2>できること</h2>

<ul><li>配布せずに，個人や会社内で利用することには何の制約もありません．<br>
</li><li>動的リンクを用いる場合ならLGPLソフトウェアから以外でも，本体部分のソースを公開する必要なく配布することができます．</li></ul>


<h2>しなければならないこと</h2>

<ul><li>DelFEMを使用していることを明記しなければなりません。<br>
</li><li>配布を求められた際にDelFEMの配布を行う義務があります。<br>
</li><li>DelFEM自体に手を加えた場合には、手を加えた部分のソースコードの公開が必要です。<br>
</li><li>スタティックリンクした際には、独自開発部分のソースコードを配布するか、独自開発部分のオブジェクトを配布する必要があります。<br>
</li><li>作成したアプリケーションについて、リバースエンジニアリングを禁止してはなりません。</li></ul>

その他，オープンソースのメリットについては以下に詳しくまとまっています．<br>
<br>
何故、オープンソースなのか？ : <a href='http://www.jp.redhat.com/opensource/os_jp.html'>http://www.jp.redhat.com/opensource/os_jp.html</a>

LGPLについて，詳しくは以下を参考にしてください．<br>
<br>
GNU 劣等一般公衆利用許諾契約書 : <a href='http://www.opensource.jp/lesser/lgpl.ja.html'>http://www.opensource.jp/lesser/lgpl.ja.html</a>

<hr />
<h1>How to build?</h1>

<h2>build with Visual Studio</h2>

<del>ライブラリ単体でビルドするのではなく，プログラムと一緒にビルドするような仕様になっています．インクルードファイルのディレクトリをVisualStudioに設定してください．例えばVsualC++6.0では[ビルド]メニューバー　>>　[オプション]メニュー >> [ディレクトリ]タブ >> [インクルードファイル]コンボボックス　でDelFEM/includeの場所を指定して下さい．開発プロジェクトに 静的ライブラリ(DelFEM/lib/stlib)をインポートして，依存関係を設定して下さい．詳しいやり方はサンプルプログラム(DelFEM/test_glut/)にならって下さい</del>

開発者がMacに乗り換えたために，Visual Studioのプロジェクトファイルは用意していません．ご了承下さい．<br>
<br>
<h2>build with MinGW</h2>

<ol><li>MinGWのインストールされたディレクトリの中のMinGW/include/の下にdelfemという名前のフォルダを作る．<br>
</li><li>作成したMinGW/include/delfemのフォルダの下にDelFEM/include/delfemの内容をコピーする．<br>
</li><li>DelFEM/MakefileをMakeする<br>
</li><li>makeして生成されるDelFEM/lib/libdfm.aをMinGW/lib/以下にコピーする．</li></ol>

以上の作業は，コマンドプロンプトでコマンドを打つ替りにバッチファイル  DelFEM/make_with_mingw.batを実行しても実現されます．<br>
<br>
<h2>Build with Mac XCode</h2>

XCodeを用いたデモプログラムを纏めたプロジェクトファイルが/test_glut/tests.xcodeprojにあります．<br>
<br>
XCodeを用いたDelFEMの静的ライブラリは/lib/xcodelib_coreにあります．これをリンクして使用してください．<br>
<br>
<h2>Build Demo</h2>

DelFEM/test_glut/にデモがあります．<br>
<br>
デモをコンパイルするにはFreeGLUTが必要です．<br>
以下を参考にインストールして下さい．<br>
<a href='http://www.transmissionzero.co.uk/software/freeglut-devel/'>http://www.transmissionzero.co.uk/software/freeglut-devel/</a>


<hr />
<h1>Updates</h1>

ver 1.2.3 :二次元 CADでの形状のインタラクティブな編集デモと，初歩的な三次元CADのデモを公開．二次元CADに初歩的な３次Bezier関数の機能を追加<br>
<br>
ver 1.2.2 : CFieldValueSetterクラスの追加，CCadObj2Dへの形状処理結果の取得クラスの追加など様々なリファクタリングを行いより使いやすくしました．<br>
<br>
ver 1.2.1 : Qtを用いたGUIつきのデモを追加しました．<br>
<br>
ver 1.2.0 : 剛体解析，剛体と弾性体の連成機能，非圧縮性超弾性体の解析，MacのXCodeでの開発環境を追加しました<br>
<br>
<br>
<br>
<br>
<br>
<hr />
<h1>Acknowledgements</h1>

<h2>支援</h2>

IPA未踏ユース２００８年度上半期の支援をうけて作られています．<br>
<br>
「インタラクティブなＵＩ備えた統合型設計解析ソフトウェアの開発」： <a href='http://www.ipa.go.jp/jinzai/mitou/2008/2008_1/youth/gaiyou/t-02.html'>http://www.ipa.go.jp/jinzai/mitou/2008/2008_1/youth/gaiyou/t-02.html</a>

<h2>使わせていただいたソース</h2>

以下のソースを使わせていただいています．<br>
<br>
"Uglyfont" Soji Yamakawaさん作 : <a href='http://homepage3.nifty.com/ysflight/uglyfont/uglyfontj.html'>http://homepage3.nifty.com/ysflight/uglyfont/uglyfontj.html</a>

どうもありがとうございます！<br>
<br>
<h2>貢献してくれた方々</h2>

PENGUINITIS様： <a href='http://www.geocities.jp/penguinitis2002/'>http://www.geocities.jp/penguinitis2002/</a>

各種環境でのビルド，開発環境の構築方法，インタラクティブな３D梁の変形デモなどなどについてWebページを作って教えて頂きました．<br>
<br>
malibu-bulldog様： <a href='http://d.hatena.ne.jp/malibu-bulldog/'>http://d.hatena.ne.jp/malibu-bulldog/</a>

Ubuntuでのビルドについて教えて頂きました．<br>
<br>
<br>
その他，色々な人からアドバイスを頂いたり助けてもらいました．どうもありがとうございます！