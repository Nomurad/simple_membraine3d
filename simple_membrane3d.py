import numpy as np
import copy
import csv 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import math
import time
from numpy import linalg as LA
# プログラムパスは~/lpthw/simple_spring.py
def prepare(): # 空リストに各粒子の番号と初期座標を突っ込む
    input_loaded = np.loadtxt('input_mempos.csv', delimiter=',', skiprows=1, dtype=float)
    for _id, (x,y,z) in enumerate(input_loaded):
        particles.append([_id, np.array([x,y,z])])
    print(particles)

    input_unload = np.loadtxt('input_mempos.csv', delimiter=',', skiprows=1, dtype=float)
    for __id, (x,y,z) in enumerate(input_unload):
        particles_0.append([_id, np.array([x,y,z])])
    print(particles_0)

    input_force = np.loadtxt('input_dist_force.csv', delimiter=',', skiprows=6, dtype=float)
    for ___id, prs in enumerate(input_force):
        forces.append(prs[2])

    input_velo = np.loadtxt('input_velo.csv', delimiter=',', dtype=float)
    for ____id, (v_x,v_y,v_z) in enumerate(input_velo):
        velo.append(np.array([v_x,v_y,v_z]))
    

def init_calc(): # 各粒子の初期加速度と初期速度を求めたい。
    for i, p0 in enumerate(forces): # 初期加速度の設定
        # print(i)
        if i <= 0 or i >= len(forces)-1:
            acce.append([i, np.array([0.0, 0.0, 0.0])])
        else: #初期外力（加速度）をここに打ち込む
            j = np.array([0.0, 1.0, 0.0]) * forces[i] # 本チャンの３次元計算の時はこことNをいじる。
            acce.append([i, j]) #ここまでで各粒子についての加速度ベクトルが求まる。
    print(acce)

            



def calculation():
    def f_v(_a, _vs, _Ps, _Ps0): # _aはスカラ, _vsはベクトル, _Ps, _Ps0は３行2列の行列
        """ 
            [point positions]
                0 - 1 - 2
                |   |   |
                3 - 4 - 5
                |   |   |
                6 - 7 - 8

        """
        #中心点から各点へのベクトル
        x40 = (_Ps[0]-_Ps[4]) 
        x41 = (_Ps[1]-_Ps[4]) 
        x42 = (_Ps[2]-_Ps[4]) 
        x43 = (_Ps[3]-_Ps[4]) 
        x45 = (_Ps[5]-_Ps[4]) 
        x46 = (_Ps[6]-_Ps[4]) 
        x47 = (_Ps[7]-_Ps[4]) 
        x48 = (_Ps[8]-_Ps[4])
        print('x41:', x41)
        print('x43:', x43)
        print('x47:', x47)
        print('x45:', x45)
        x040 = (_Ps0[0]-_Ps0[4]) 
        x041 = (_Ps0[1]-_Ps0[4]) 
        x042 = (_Ps0[2]-_Ps0[4]) 
        x043 = (_Ps0[3]-_Ps0[4]) 
        x045 = (_Ps0[5]-_Ps0[4]) 
        x046 = (_Ps0[6]-_Ps0[4]) 
        x047 = (_Ps0[7]-_Ps0[4]) 
        x048 = (_Ps0[8]-_Ps0[4]) 
        #中心点周りの面の面積
        s4103 = LA.norm(np.cross(x41,x40))/2\
                +LA.norm(np.cross(x40,x43))/2
        s4367 = LA.norm(np.cross(x43,x46))/2\
                +LA.norm(np.cross(x46,x47))/2
        s4785 = LA.norm(np.cross(x47,x48))/2\
                +LA.norm(np.cross(x48,x45))/2
        s4521 = LA.norm(np.cross(x45,x42))/2\
                +LA.norm(np.cross(x42,x41))/2
        s04103 = LA.norm(np.cross(x041,x040))/2\
                +LA.norm(np.cross(x040,x043))/2
        s04367 = LA.norm(np.cross(x043,x046))/2\
                +LA.norm(np.cross(x046,x047))/2
        s04785 = LA.norm(np.cross(x047,x048))/2\
                +LA.norm(np.cross(x048,x045))/2
        s04521 = LA.norm(np.cross(x045,x042))/2\
                +LA.norm(np.cross(x042,x041))/2
        
        #各方向への平均面積(ここだけ反時計回り順で設定してる)
        S_iminus = (s4103 + s4367) / 2 #43方向
        S_Jminus = (s4103 + s4521) / 2 #41方向
        S_iplus = (s4785 + s4521) / 2 #45方向
        S_Jplus = (s4785 + s4367) / 2 #47方向
        S_iminus0 = (s04103 + s04367) / 2 #43方向
        S_Jminus0 = (s04103 + s04521) / 2 #41方向
        S_iplus0 = (s04785 + s04521) / 2 #45方向
        S_Jplus0 = (s04785 + s04367) / 2 #47方向
        # 各方向への厚み
        h_iminus = h_0 / ((poisson/(1-poisson) * (S_iminus - S_iminus0) / S_iminus0) + 1) #43方向
        h_Jminus = h_0 / ((poisson/(1-poisson) * (S_Jminus - S_Jminus0) / S_Jminus0) + 1) #41方向
        h_iplus = h_0 / ((poisson/(1-poisson) * (S_iplus - S_iplus0) / S_iplus0) + 1) #45方向
        h_Jplus = h_0 / ((poisson/(1-poisson) * (S_Jplus - S_Jplus0) / S_Jplus0) + 1) #47方向
        # 各断片の重心
        g401 = (_Ps[4] + _Ps[0] + _Ps[1]) / 3
        g430 = (_Ps[4] + _Ps[3] + _Ps[0]) / 3
        g436 = (_Ps[4] + _Ps[3] + _Ps[6]) / 3
        g467 = (_Ps[4] + _Ps[6] + _Ps[7]) / 3
        g478 = (_Ps[4] + _Ps[7] + _Ps[8]) / 3
        g485 = (_Ps[4] + _Ps[8] + _Ps[5]) / 3
        g452 = (_Ps[4] + _Ps[5] + _Ps[2]) / 3
        g421 = (_Ps[4] + _Ps[2] + _Ps[1]) / 3
        g0401 = (_Ps0[4] + _Ps0[0] + _Ps0[1]) / 3
        g0430 = (_Ps0[4] + _Ps0[3] + _Ps0[0]) / 3
        g0436 = (_Ps0[4] + _Ps0[3] + _Ps0[6]) / 3
        g0467 = (_Ps0[4] + _Ps0[6] + _Ps0[7]) / 3
        g0478 = (_Ps0[4] + _Ps0[7] + _Ps0[8]) / 3
        g0485 = (_Ps0[4] + _Ps0[8] + _Ps0[5]) / 3
        g0452 = (_Ps0[4] + _Ps0[5] + _Ps0[2]) / 3
        g0421 = (_Ps0[4] + _Ps0[2] + _Ps0[1]) / 3
        # 各断片面積
        s410 = LA.norm(np.cross(x41,x40))/2
        s403 = LA.norm(np.cross(x40,x43))/2
        s436 = LA.norm(np.cross(x43,x46))/2
        s467 = LA.norm(np.cross(x46,x47))/2
        s478 = LA.norm(np.cross(x47,x48))/2
        s485 = LA.norm(np.cross(x48,x45))/2
        s452 = LA.norm(np.cross(x45,x42))/2
        s421 = LA.norm(np.cross(x42,x41))/2
        s0410 = LA.norm(np.cross(x041,x040))/2
        s0403 = LA.norm(np.cross(x040,x043))/2
        s0436 = LA.norm(np.cross(x043,x046))/2
        s0467 = LA.norm(np.cross(x046,x047))/2
        s0478 = LA.norm(np.cross(x047,x048))/2
        s0485 = LA.norm(np.cross(x048,x045))/2
        s0452 = LA.norm(np.cross(x045,x042))/2
        s0421 = LA.norm(np.cross(x042,x041))/2
        # 四角の重心
        g4103 = (s410*g401 + s403*g430) / (s410 + s403)
        g4367 = (s436*g436 + s467*g467) / (s436 + s467)
        g4785 = (s478*g478 + s485*g485) / (s478 + s485)
        g4521 = (s452*g452 + s421*g421) / (s452 + s421)
        g04103 = (s0410*g0401 + s0403*g0430) / (s0410 + s0403)
        g04367 = (s0436*g0436 + s0467*g0467) / (s0436 + s0467)
        g04785 = (s0478*g0478 + s0485*g0485) / (s0478 + s0485)
        g04521 = (s0452*g0452 + s0421*g0421) / (s0452 + s0421)
        # 各重心間の距離
        Lj20 = LA.norm(g4521 - g4103)
        Lj06 = LA.norm(g4103 - g4367)
        Lj68 = LA.norm(g4367 - g4785)
        Lj82 = LA.norm(g4785 - g4521)
        
        # ひずみ
        eps_i41 = (LA.norm(x41) - LA.norm(x041)) / LA.norm(x041)
        eps_J41 = (LA.norm(g4521 - g4103) - LA.norm(g04521 - g04103)) / LA.norm(g04521 - g04103)
        eps_i43 = (LA.norm(x43) - LA.norm(x043)) / LA.norm(x043)
        eps_J43 = (LA.norm(g4103 - g4367) - LA.norm(g04103 - g04367)) / LA.norm(g04103 - g04367)
        eps_i47 = (LA.norm(x41) - LA.norm(x041)) / LA.norm(x041)
        eps_J47 = (LA.norm(g4367 - g4785) - LA.norm(g04367 - g04785)) / LA.norm(g04367 - g04785)
        eps_i45 = (LA.norm(x41) - LA.norm(x041)) / LA.norm(x041)
        eps_J45 = (LA.norm(g4785 - g4521) - LA.norm(g04785 - g04521)) / LA.norm(g04785 - g04521)
        # 張力
        F_T41 = (young_modulus * h_Jminus * Lj20 * (eps_i41 + poisson * eps_J41) / (1 - poisson**2))*x41/LA.norm(x41)
        F_T43 = (young_modulus * h_iminus * Lj06 * (eps_i43 + poisson * eps_J43) / (1 - poisson**2))*x43/LA.norm(x43)
        F_T47 = (young_modulus * h_Jplus * Lj68 * (eps_i47 + poisson * eps_J47) / (1 - poisson**2))*x47/LA.norm(x47)
        F_T45 = (young_modulus * h_iplus * Lj82 * (eps_i45 + poisson * eps_J45) / (1 - poisson**2))*x45/LA.norm(x45)
        # せん断ひずみ
        gamma513 = (math.acos((np.dot(x45,x41))/(LA.norm(x45)*LA.norm(x41))) - math.acos((np.dot(x045,x041))/(LA.norm(x045)*LA.norm(x041)))\
            + math.acos((np.dot(x43,x41))/(LA.norm(x43)*LA.norm(x41))) - math.acos((np.dot(x043,x041))/(LA.norm(x043)*LA.norm(x041))))/2
        gamma137 = (math.acos((np.dot(x41,x43))/(LA.norm(x41)*LA.norm(x43))) - math.acos((np.dot(x041,x043))/(LA.norm(x041)*LA.norm(x043)))\
            + math.acos((np.dot(x43,x47))/(LA.norm(x43)*LA.norm(x47))) - math.acos((np.dot(x043,x047))/(LA.norm(x043)*LA.norm(x047))))/2
        gamma375 = (math.acos((np.dot(x47,x43))/(LA.norm(x47)*LA.norm(x43))) - math.acos((np.dot(x047,x043))/(LA.norm(x047)*LA.norm(x043)))\
            + math.acos((np.dot(x45,x47))/(LA.norm(x45)*LA.norm(x47))) - math.acos((np.dot(x045,x047))/(LA.norm(x045)*LA.norm(x047))))/2
        gamma751 = (math.acos((np.dot(x47,x45))/(LA.norm(x47)*LA.norm(x45))) - math.acos((np.dot(x047,x045))/(LA.norm(x047)*LA.norm(x045)))\
            + math.acos((np.dot(x45,x41))/(LA.norm(x45)*LA.norm(x41))) - math.acos((np.dot(x045,x041))/(LA.norm(x045)*LA.norm(x041))))/2
        # せん断力
        F_S41 = ((young_modulus * h_Jminus * LA.norm(x41) * gamma513)/(2 * (1 + poisson)))*x41/LA.norm(x41)
        F_S43 = ((young_modulus * h_Jminus * LA.norm(x43) * gamma137)/(2 * (1 + poisson)))*x43/LA.norm(x43)
        F_S47 = ((young_modulus * h_Jminus * LA.norm(x47) * gamma375)/(2 * (1 + poisson)))*x47/LA.norm(x47)
        F_S45 = ((young_modulus * h_Jminus * LA.norm(x45) * gamma751)/(2 * (1 + poisson)))*x45/LA.norm(x45)
        #[m/s] Re : , Temp : 20 [℃], Prs : 0.1013 [MPa], Chord : 1 [m]
        # C_v = 1
        # 曲げ力の方向
        # dx_J = _Ps[7][0] - _Ps[1][0]
        # dy_J = _Ps[7][1] - _Ps[1][1]
        # thetaJ = np.arctan2(dy_J, dx_J) - np.pi / 2
        # e_J = np.array([np.cos(thetaJ), np.sin(thetaJ)])
        # dx_J0 = _Ps0[7][0] - _Ps0[1][0]
        # dy_J0 = _Ps0[7][1] - _Ps0[1][1]
        # thetaJ0 = np.arctan2(dy_J0, dx_J0) - np.pi / 2
        # e_J0 = np.array([np.cos(thetaJ0), np.sin(thetaJ0)])
        # J方向の曲げ力
        n_J = np.cross(x47,x41)/LA.norm(np.cross(x41,x47))
        l_Jalfa = LA.norm(_Ps[1] - _Ps[7])
        cos_Jalfa = (LA.norm(x41)**2 + LA.norm(x47)**2 - l_Jalfa**2) / (2 * LA.norm(x41) * LA.norm(x47))
        if cos_Jalfa > 1.0:
            cos_Jalfa = 1.0
        elif cos_Jalfa < -1.0:
            cos_Jalfa = -1.0
        sin_Jalfa = math.sqrt(1 - cos_Jalfa**2)
        CJa2 = math.sqrt((cos_Jalfa + 1)/2)
        SJa2 = math.sqrt((1 - cos_Jalfa)/2)
        zJC = (_Ps[7][2]-_Ps[1][2])/(_Ps[7][0]-_Ps[1][0]) * (_Ps[4][0]-_Ps[1][0]) + _Ps[1][2] #曲げ力の方向の場合わけに必要
        if _Ps[4][2] > zJC:
            e_j = np.dot(np.array([[CJa2 + (n_J[0]**2) * (1 - CJa2), n_J[0] * n_J[1] * (1 - CJa2) + n_J[2] * SJa2, n_J[0] * n_J[2] * (1 - CJa2) - n_J[1] * SJa2],\
                                    [n_J[1] * n_J[0] * (1 - CJa2) - n_J[2] * SJa2, CJa2 + (n_J[1]**2) * (1 - CJa2), n_J[1] * n_J[2] * (1 - CJa2) + n_J[0] * SJa2],\
                                    [n_J[2] * n_J[0] * (1 - CJa2) + n_J[1] * SJa2, n_J[2] * n_J[1] * (1 - CJa2) - n_J[0] * SJa2, CJa2 + (n_J[2]**2) * (1 - CJa2)]]), (_Ps[7] - _Ps[4])/LA.norm(_Ps[7] - _Ps[4]))
        else:
            e_j = np.dot(np.array([[CJa2 + (n_J[0]**2) * (1 - CJa2), n_J[0] * n_J[1] * (1 - CJa2) - n_J[2] * SJa2, n_J[0] * n_J[2] * (1 - CJa2) + n_J[1] * SJa2],\
                                    [n_J[1] * n_J[0] * (1 - CJa2) + n_J[2] * SJa2, CJa2 + (n_J[1]**2) * (1 - CJa2), n_J[1] * n_J[2] * (1 - CJa2) - n_J[0] * SJa2],\
                                    [n_J[2] * n_J[0] * (1 - CJa2) - n_J[1] * SJa2, n_J[2] * n_J[1] * (1 - CJa2) + n_J[0] * SJa2, CJa2 + (n_J[2]**2) * (1 - CJa2)]]), (_Ps[7] - _Ps[4])/LA.norm(_Ps[7] - _Ps[4]))
        d_etha_J = (2 * sin_Jalfa / l_Jalfa) - (2 * math.sqrt(1 - np.dot(x041,x047)**2/(LA.norm(x041)*LA.norm(x047))**2)/(LA.norm(x041 - x047)))

        n_i = np.cross(x45,x43)/LA.norm(np.cross(x43,x45))  
        cos_ialfa = np.dot(x43,x45) / (LA.norm(x43) * LA.norm(x45))
        sin_ialfa = math.sqrt(1 - cos_ialfa**2)
        Cia2 = math.sqrt((cos_ialfa + 1)/2)
        Sia2 = math.sqrt((1 - cos_ialfa)/2)
        ziC = (_Ps[5][2]-_Ps[3][2])/(_Ps[5][0]-_Ps[3][0]) * (_Ps[4][0]-_Ps[3][0]) + _Ps[3][2]
        if _Ps[4][2] > ziC:
            e_i = np.dot(np.array([[Cia2 + (n_i[0]**2) * (1 - Cia2), n_i[0] * n_i[1] * (1 - Cia2) + n_i[2] * Sia2, n_i[0] * n_i[2] * (1 - Cia2) - n_i[1] * Sia2],\
                                    [n_i[1] * n_i[0] * (1 - Cia2) - n_i[2] * Sia2, Cia2 + (n_i[1]**2) * (1 - Cia2), n_i[1] * n_i[2] * (1 - Cia2) + n_i[0] * Sia2],\
                                    [n_i[2] * n_i[0] * (1 - Cia2) + n_i[1] * Sia2, n_i[2] * n_i[1] * (1 - Cia2) - n_i[0] * Sia2, Cia2 + (n_i[2]**2) * (1 - Cia2)]]), (_Ps[7] - _Ps[4])/LA.norm(_Ps[7] - _Ps[4]))
        else:
            e_i = np.dot(np.array([[Cia2 + (n_i[0]**2) * (1 - Cia2), n_i[0] * n_i[1] * (1 - Cia2) - n_i[2] * Sia2, n_i[0] * n_i[2] * (1 - Cia2) + n_i[1] * Sia2],\
                                    [n_i[1] * n_i[0] * (1 - Cia2) + n_i[2] * Sia2, Cia2 + (n_i[1]**2) * (1 - Cia2), n_i[1] * n_i[2] * (1 - Cia2) - n_i[0] * Sia2],\
                                    [n_i[2] * n_i[0] * (1 - Cia2) - n_i[1] * Sia2, n_i[2] * n_i[1] * (1 - Cia2) + n_i[0] * Sia2, Cia2 + (n_i[2]**2) * (1 - Cia2)]]), (_Ps[5] - _Ps[4])/LA.norm(_Ps[5] - _Ps[4]))
        d_etha_i = (2 * sin_ialfa / LA.norm(x45 - x43)) - (2 * math.sqrt(1 - np.dot(x043,x045)**2/(LA.norm(x043)*LA.norm(x045))**2)/(LA.norm(x043 - x045)))


        l_J = (Lj20 + Lj06 + Lj68 + Lj82) / 4
        h = (h_iminus + h_iplus + h_Jminus + h_Jplus) / 4
        I = (l_J * h**3) / 12
        M_i = (young_modulus * I * (d_etha_i + poisson * d_etha_J)/(1 - poisson**2))
        M_J = (young_modulus * I * (d_etha_J + poisson * d_etha_i)/(1 - poisson**2))
        #曲げ力
        F_Bi = M_i / LA.norm(x43) + M_i / LA.norm(x45) * e_i
        F_BJ = M_J / LA.norm(x41) + M_J / LA.norm(x47) * e_j
        #空気力
        # S = (S_iminus + S_iplus + S_Jminus + S_Jplus) / 4
        # F_A = p * S
        F_A = np.array([0.0, 0.0, -0.1]) * _a

        # 運動方程式（支配方程式）
        S_0 = (S_iminus0 + S_iplus0 + S_Jminus0 + S_Jplus0) / 4
        F_T = F_T41 + F_T43 + F_T45 + F_T47
        F_S = F_S41 + F_S43 + F_S45 + F_S47
        F_B = F_Bi + F_BJ
        return (F_T + F_S + F_B + F_A) / (rho * h_0 * S_0) - c * _vs

    for i, (x,y) in enumerate(particles):
        
        if i%11 == 0 or i%11 == 10: #変位を与えない点の条件を与える。
            pass
        elif i >= 1 and i <= 9:
            pass
        elif i >= 111 and i <= 119:
            pass
        else:
            acc = acce[i][1] # 粒子２の加速度
            v2 = velo[i] # 粒子2の速度
            p1 = particles[i-12][1] # 粒子1(右隣)の座標
            p2 = particles[i-11][1] # 粒子2(真ん中)の座標
            p3 = particles[i-10][1] # 粒子3(左隣)の座標
            p4 = particles[i-1][1]
            p5 = particles[i][1]
            p6 = particles[i+1][1]
            p7 = particles[i+10][1]
            p8 = particles[i+11][1]
            p9 = particles[i+12][1]
            Ps = np.array([p1, p2, p3, p4, p5, p6, p7, p8, p9]) # 対象粒子、隣り合う粒子座標のセットリスト
            p1_0 = particles_0[i-12][1] # 粒子1の初期座標
            p2_0 = particles_0[i-11][1] # 粒子2の初期座標
            p3_0 = particles_0[i-10][1] # 粒子3の初期座標
            p4_0 = particles_0[i-1][1] 
            p5_0 = particles_0[i][1] 
            p6_0 = particles_0[i+1][1]
            p7_0 = particles_0[i+10][1] 
            p8_0 = particles_0[i+11][1] 
            p9_0 = particles_0[i+12][1] 
            Ps0 = np.array([p1_0, p2_0, p3_0, p4_0, p5_0, p6_0, p7_0, p8_0, p9_0]) # 初期座標
        

            if i <= 0 or i >= len(particles)-1:
                v2[0] = 0.0
                v2[1] = 0.0
                v2[2] = 0.0
            else:
                # print("cal2")
                # print("Vs: ", Vs)
                # kr1 = dt * Vs
                # kv1 = dt * np.array([np.array([0.0,0.0]), f_v(acc, v2, Ps, Ps0), np.array([0.0,0.0])])
                # kr2 = dt * (Vs + kv1 / 2)
                # kv2 = dt * np.array([np.array([0.0,0.0]), f_v(acc, v2 + kv1[1] / 2, Ps + kr1 / 2, Ps0), np.array([0.0,0.0])])
                # kr3 = dt * (Vs + kv2 / 2)
                # kv3 = dt * np.array([np.array([0.0,0.0]), f_v(acc, v2 + kv2[1] / 2, Ps + kr2 / 2, Ps0), np.array([0.0,0.0])])
                # kr4 = dt * (Vs + kv3)
                # kv4 = dt * np.array([np.array([0.0,0.0]), f_v(acc, v2 + kv3[1], Ps + kr3, Ps0), np.array([0.0,0.0])])
                kr1 = dt * v2
                kv1 = dt * f_v(acc, v2, Ps, Ps0)
                kr2 = dt * (v2 + kv1/2)
                kv2 = dt * f_v(acc, v2+kv1/2, Ps, Ps0)
                kr3 = dt * (v2 + kv2/2)
                kv3 = dt * f_v(acc, v2+kv2/2, Ps, Ps0)
                kr4 = dt * (v2 + kv3)
                kv4 = dt * f_v(acc, v2+kv3, Ps, Ps0)


                particles[i][1] += (kr1 + 2 * kr2 + 2 * kr3 + kr4) / 6
                velo[i] += (kv1 + 2 * kv2 + 2 * kv3 + kv4) / 6
                # print("acc: ", acc)
                # print("velo: ", velo[1])
                

if __name__ == "__main__":
    time = 0.0
    dt = 0.001
    rho = 910 # ρ
    h_0 = 0.0004 #厚み 
    poisson = 0.5
    young_modulus = 2 * 0.84 * 10**6 * (1 + poisson)
    c = 1000
    velo = []
    acce = []
    history = []
    particles = []
    particles_0 = []
    forces = []

    fig = plt.figure()
    ims = []
    ax = fig.add_subplot(111, projection='3d')


    prepare()
    init_calc()
    # while time <= 0.01:

    #     X = []
    #     Y = []
    #     Z = []
    #     if time%10.0 <=0.01:
    #         history.append([time, ',', particles[1][1][0], ',',particles[1][1][1]])
            
    #     calculation()
    #     time += dt
    # for i, P in enumerate(particles):

    #     X.append(P[1][0])
    #     Y.append(P[1][1])
    #     Z.append(P[1][2])
    # ims = []
    # im = ax.scatter(X,Y,Z,c='c')
    # np.savetxt("time_history.csv", history, fmt="%s")
    # plt.show