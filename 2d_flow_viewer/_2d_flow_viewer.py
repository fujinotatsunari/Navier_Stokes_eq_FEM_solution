

import numpy as np
import matplotlib.pyplot as plt
import csv

from numpy import sin, cos, pi, sqrt

#スカラー場
class Scalar:
    def __init__(self, value):
        self.value = value #節点に対応する値

    #ファイル名からデータをインプット
    def inputfile(self, filename):
        data = Inputcsv(filename)
        self.value = data.array2df()

#ベクター場
class Vector:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    #ファイル名からデータをインプット
    def inputfile(self, Ufile, Vfile):
        Udata = Inputcsv(Ufile)
        Vdata = Inputcsv(Vfile)
        
        self.x = Udata.array2df()
        self.y = Vdata.array2df()

        #ファイル名からデータをインプット(1時限配列版)
    def inputfile1d(self, Ufile, Vfile):
        Udata = Inputcsv(Ufile)
        Vdata = Inputcsv(Vfile)
        
        self.x = Udata.array1df()
        self.y = Vdata.array1df()
#座標データ(正方形格子)
class SquareGrid2d:
    def __init__(self, xb, xt, yb, yt, xnode, ynode):
        self.xb = xb
        self.xt = xt
        self.yb = yb
        self.yt = yt
        self.xnode = xnode
        self.ynode = ynode

    #x軸の提供
    def x(self):
        x = np.linspace(self.xb, self.xt, self. xnode) #x軸を確保
        return x

    #y軸の提供
    def y(self):
        y = np.linspace(self.yb, self.yt, self. ynode) #y軸を確保
        return y

    def mesh(self):
        mesh = SquareGrid2d(self.xb, self.xt, self.yb, self.yt, self.xnode, self.ynode)
        x = mesh.x()
        y = mesh.y()
        X,Y = np.meshgrid(x,y)
        return X,Y


#csv読み込みデータ
class Inputcsv:
    def __init__(self, name):
        self.name = name #ファイル名
        
    #csvの行r, 列cを指定してその要素を返す(float)
    def elemf(self, r, c):
        openfile= 'C:/Result/2d_Navier_Stokes_eq/input/' + self.name
        array = np.loadtxt(openfile, delimiter=',')
        return array[r][c]


    #csvデータを一次元配列として返す
    def array1df(self):
        openfile= 'C:/Result/2d_Navier_Stokes_eq/input/' + self.name
        array = np.loadtxt(openfile, delimiter=',')
        array.ravel()
        return array

    #csvデータを２次元配列で返す
    def array2df(self):
        openfile= 'C:/Result/2d_Navier_Stokes_eq/input/' + self.name
        array = np.loadtxt(openfile, delimiter=',')
        return array

    #列を取得
    def colarray1d(self, col):
        openfile= 'C:/Result/2d_Navier_Stokes_eq/input/' + self.name
        with open(openfile) as f:
            reader = csv.reader(f)
            arr = [row for row in reader]
            arr_T = [list(x) for x in zip(*arr)]
            
            arr_T[col] = list(filter(lambda a:a != '',arr_T[col]))#空白を削除
            float_list = map(float, arr_T[col])
            return list(float_list)

#長方形グリッドの生成関数
def mesh_gen():
    print("x領域左端 xb->")
    xb = float(input())
    print("x領域右端 xt->")
    xt = float(input())
    print("y領域左端 yb->")
    yb = float(input())
    print("y領域右端 yt->")
    yt = float(input())
    print("x方向節点数 xnode->")
    xnode = int(input())
    print("y方向節点数 ynode->")
    ynode = int(input())

    mesh = SquareGrid2d(xb, xt, yb, yt, xnode, ynode)
    return mesh

#コンター図の描画
def plot_contor():

    grid = mesh_gen()
    v = Scalar(None)
    print("input file name->")
    filename = input()
    v.inputfile(filename)
    X,Y = np.meshgrid(grid.x(), grid.y())
   
    plt.contourf(X, Y, v.value, cmap = 'jet')
    plt.xlim(grid.xb, grid.xt)
    plt.ylim(grid.yb, grid.yt)
    plt.axis('equal')
    plt.xlabel('X')
    plt.ylabel('Y');
    plt.colorbar()
    plt.show()



#ベクトル場の描画
def plot_vector():
     grid = mesh_gen()
     V = Vector(None, None)
     print("x-component input file name->")
     Ufile = input()
     print("y-component input file name->")
     Vfile = input()
     V.inputfile(Ufile, Vfile)
     X,Y = np.meshgrid(grid.x(), grid.y())
     u = V.x
     v = V.y
     plt.quiver(X, Y, u, v, np.sqrt(u**2 + v**2), cmap = 'jet')
     plt.xlim(grid.xb, grid.xt)
     plt.ylim(grid.yb, grid.yt)
     plt.axis('equal')
     plt.xlabel('X')
     plt.ylabel('Y')
     plt.show()


#流線の描画
def plot_streamline():
    grid = mesh_gen()
    V = Vector(None, None)
    print("x-component input file name->")
    Ufile = input()
    print("y-component input file name->")
    Vfile = input()
    V.inputfile(Ufile, Vfile)
    dx = (grid.xt - grid.xb)/(grid.xnode - 1)
    dy = (grid.yt - grid.yb)/(grid.ynode - 1)

    X,Y = np.meshgrid(grid.x(), grid.y())
    u = V.x
    v = V.y
    plt.streamplot(X, Y, u, v, color=(np.sqrt(u**2 + v**2)), cmap ='jet')
    plt.xlim(grid.xb, grid.xt)
    plt.ylim(grid.yb, grid.yt)
    plt.axis('equal')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()



print("2次元流れ場の可視化")
print("0:圧力場, 1:流れ場(スカラー), 2:流れ場(ベクトル) ->")
flag = int(input())
if flag == 0:
    plot_contor()

elif flag == 1:
    plot_contor()

elif flag == 2:
    print("可視化の方法")
    print("0:ベクトル図, 1:流線図 ->")
    fflag = int(input())
    if fflag == 0:
        plot_vector()
    elif fflag == 1:
        plot_streamline()

