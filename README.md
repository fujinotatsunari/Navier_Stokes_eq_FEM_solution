2次元長方形領域における無次元化Navier-Stokes方程式の有限要素法(FEM)によるシミュレーションプログラムです.オブジェクト指向と数値計算の勉強のため作成しました.

説明:
Mesh_generate:(C++)
計算領域と計算モデル,初期条件を指定し, 作成します.
モデルはcavity流れ(cavity), backstep流れ(backstep),角柱周り流れ(pillar)の3つがチュートリアルとして提供されます.(pillarがうまく動きません. 現在修正中)

calculation:(C++)
作成したmeshデータからNavier-Stokes方程式を有限要素法で数値的に解きます.
スキームはHSMAC法, 直接完全陰解法(未実装), 線形化直接完全陰解法です.

viewer:(Python)
calculationで出力されたデータ(U,V,Pnode,magnitude)をコンター図(Pnode, magnitude), ベクトル場図(U,V), 流線図(U,V)にして可視化します.


使い方:
visual studioで上のプログラムを起動します. ターミナルに表示されるパラメータを入力し, 計算を行います.
計算に伴うディレクトリは基本自動で作られます.(C:以下に作られる)
ただし, viewerのみ手動で出力したいデータをディレクトリに入れないといけません.(C:/Result/2d_Navier_Stokes_eq/input/内)



動作環境:
Visual Studio 2022
