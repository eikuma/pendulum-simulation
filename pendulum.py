# 課題4，5

import matplotlib.pyplot as plt
import numpy as np
import math

L = 1.0+13/100  # 振り子の長さ
M = 1.0  # 質点の質量
# M = 1000
G = 9.8  # 重力加速度
B = 0.8  # 粘性摩擦係数[kg*m^2/s]
H = 0.01  # 刻み幅
TIME = 20  # 総シミュレーション時間[s]


def func_dTheta1(Theta1, Theta2):
    return Theta2


def func_dTheta2(Theta1, Theta2):
    return -(G/L)*math.sin(Theta1)-(B/(M*L**2))*Theta2


# 状態変数の定義&初期化
Theta1 = 1.0  # 質点の位置　Theta1
Theta2 = 0.0  # 質点の速度　Theta2
t = 0  # 経過時間に相当

# 一時記憶の変数
k1_Theta1 = 0.0
k2_Theta1 = 0.0
k3_Theta1 = 0.0
k4_Theta1 = 0.0
k1_Theta2 = 0.0
k2_Theta2 = 0.0
k3_Theta2 = 0.0
k4_Theta2 = 0.0

# 積分法の結果を格納するリスト
Theta1_list = [0]*(int(TIME/H)+1)
Theta2_list = [0]*(int(TIME/H)+1)

# シミュレーション
for t in range(0, int(TIME/H)+1):

    k1_Theta1 = func_dTheta1(Theta1, Theta2)
    k1_Theta2 = func_dTheta2(Theta1, Theta2)

    k2_Theta1 = func_dTheta1(Theta1+(k1_Theta1*H/2), Theta2+(k1_Theta2*H/2))
    k2_Theta2 = func_dTheta2(Theta1+(k1_Theta1*H/2), Theta2+(k1_Theta2*H/2))

    k3_Theta1 = func_dTheta1(Theta1+(k2_Theta1*H/2), Theta2+(k2_Theta2*H/2))
    k3_Theta2 = func_dTheta2(Theta1+(k2_Theta1*H/2), Theta2+(k2_Theta2*H/2))

    k4_Theta1 = func_dTheta1(Theta1+k3_Theta1*H, Theta2+k3_Theta2*H)
    k4_Theta2 = func_dTheta2(Theta1+(k2_Theta1*H/2), Theta2+(k2_Theta2*H/2))

    # Theta1, Theta2 を更新する
    Theta1 = Theta1+(k1_Theta1+2*k2_Theta1+2*k3_Theta1+k4_Theta1)*H/6
    Theta2 = Theta2+(k1_Theta2+2*k2_Theta2+2*k3_Theta2+k4_Theta2)*H/6

    # リストに格納
    Theta1_list[t] = Theta1
    Theta2_list[t] = Theta2


# グラフの表示
plt.plot(np.arange(0, TIME+H/2, H), Theta1_list, label="angular")
plt.plot(np.arange(0, TIME+H/2, H), Theta2_list, label="angular velocity")
plt.xlabel("time[s]")
plt.ylabel("angular[rad],angular velocity[rad/s]")
plt.legend()
plt.show()
