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
    # init list
    particles = []
    particles_0 = []
    forces = []
    velo = []

    input_loaded = np.loadtxt('input_mempos.csv', delimiter=',', skiprows=1, dtype=float)
    for _id, (x,y,z) in enumerate(input_loaded):
        particles.append([_id, np.array([x,y,z])])
    # print(particles)

    input_unload = np.loadtxt('input_mempos.csv', delimiter=',', skiprows=1, dtype=float)
    for _id, (x,y,z) in enumerate(input_unload):
        particles_0.append([_id, np.array([x,y,z])])
    # print(particles_0)

    input_force = np.loadtxt('input_dist_force.csv', delimiter=',', skiprows=6, dtype=float)
    for prs in (input_force):
        forces.append(prs[2])

    input_velo = np.loadtxt('input_velo.csv', delimiter=',', dtype=float)
    for (v_x,v_y,v_z) in (input_velo):
        velo.append(np.array([v_x,v_y,v_z]))

    return particles, particles_0, forces, velo

def init_acceleration(): # 各粒子の初期加速度と初期速度を求めたい。
    acce = []
    for i, p0 in enumerate(forces): # 初期加速度の設定
        # print(i)
        if i <= 0 or i >= len(forces)-1:
            acce.append([i, np.array([0.0, 0.0, 0.0])])
        else: #初期外力（加速度）をここに打ち込む
            j = np.array([0.0, 1.0, 0.0]) * forces[i] # 本チャンの３次元計算の時はこことNをいじる。
            acce.append([i, j]) #ここまでで各粒子についての加速度ベクトルが求まる。
    print("acc = ", acce)
    return acce

class Index_iterator():
    def __init__(self, start, fin):
        self.index = start
        self.index_max = fin
    
    def reset(self):
        self.index = 1

    def get_indexes(self, start_idx, num):
        self.index = start_idx

        for i in range(num):
            yield self.index
            self.index += 1
            if self.index > 8:
                self.reset()

    def get_indexes_reverse(self, start_idx, num):
        self.index = start_idx

        for i in range(num):
            yield self.index
            self.index -= 1
            if self.index < 1:
                self.index = 8
    
    def __iter__(self):
        return self
    
    def __next__(self):
        if self.index > 8:
            self.index = 1
        return self.index


def calculation():
    def f_v(_a, _vs, _Ps, _Ps0): # _aはスカラ, _vsはベクトル, _Ps, _Ps0は３行2列の行列
        """ 
            [point positions]
                8 - 1 - 2
                |(4)|(1)|
                7 - 0 - 3
                |(3)|(2)|
                6 - 5 - 4

        """
        center_pos = _Ps[0]
        center_pos_0 = _Ps0[0]
        idx_iter = Index_iterator(1, 8)
        #中心点から各点へのベクトル
        x = []
        x0 = []
        for p in (_Ps):
            x.append(p - center_pos)
        for p in _Ps(_Ps0):
            x0.append(p - center_pos_0)

        x01 = (_Ps[1]-center_pos) 
        x02 = (_Ps[2]-center_pos) 
        x03 = (_Ps[3]-center_pos) 
        x04 = (_Ps[4]-center_pos) 
        x05 = (_Ps[5]-center_pos) 
        x06 = (_Ps[6]-center_pos) 
        x07 = (_Ps[7]-center_pos) 
        x08 = (_Ps[8]-center_pos)
        print('p_id', center_pos, end='\t')
        print('x01:', x01, end="\t")
        print('x03:', x03, end="\t")
        print('x05:', x05, end="\t")
        print('x07:', x07)
        x001 = (_Ps0[1]-_Ps0[0]) 
        x002 = (_Ps0[2]-_Ps0[0]) 
        x003 = (_Ps0[3]-_Ps0[0]) 
        x004 = (_Ps0[4]-_Ps0[0]) 
        x005 = (_Ps0[5]-_Ps0[0]) 
        x006 = (_Ps0[6]-_Ps0[0]) 
        x007 = (_Ps0[7]-_Ps0[0]) 
        x008 = (_Ps0[8]-_Ps0[0]) 
        
        #中心点周りの面の面積
        def calc_area(j,k,l):
            s = LA.norm(np.cross(x[j],x[k]))/2 \
                + LA.norm(np.cross(x[k],x[l]))/2
            return s

        s = []
        s0 = []
        hen = [1,3,5,7]
        for i in range(4):
            j,k,l = [n for n in idx_iter.get_indexes(start_idx=hen[i], 3)]
            s[i] = calc_area(j,k,l)
            s0[i] = calc_area(j,k,l)

        # s0123 = LA.norm(np.cross(x[1],x[2]))/2\
        #         +LA.norm(np.cross(x[2],x[3]))/2
        # s4367 = LA.norm(np.cross(x[3],x[4]))/2\
        #         +LA.norm(np.cross(x[4],x[5]))/2
        # s4785 = LA.norm(np.cross(x[5],x[6]))/2\
        #         +LA.norm(np.cross(x[6],x[7]))/2
        # s4521 = LA.norm(np.cross(x[7],x[8]))/2\
        #         +LA.norm(np.cross(x[8],x[1]))/2
        # s04103 = LA.norm(np.cross(x0[1],x0[2]))/2\
        #         +LA.norm(np.cross(x0[2],x0[3]))/2
        # s04367 = LA.norm(np.cross(x0[3],x0[4]))/2\
        #         +LA.norm(np.cross(x0[4],x0[7]))/2
        # s04785 = LA.norm(np.cross(x0[7],x0[8]))/2\
        #         +LA.norm(np.cross(x0[8],x0[5]))/2
        # s04521 = LA.norm(np.cross(x0[5],x0[2]))/2\
        #         +LA.norm(np.cross(x0[2],x0[1]))/2
        
        #各方向への平均面積(ここだけ反時計回り順で設定してる)
        S_iminus = (s[1] + s[2]) / 2 #43方向
        S_Jminus = (s[1] + s[4]) / 2 #41方向
        S_iplus = (s[3] + s[4]) / 2 #45方向
        S_Jplus = (s[3] + s[2]) / 2 #47方向
        S_iminus0 = (s0[1] + s0[2]) / 2 #43方向
        S_Jminus0 = (s0[1] + s0[4]) / 2 #41方向
        S_iplus0 = (s0[3] + s0[4]) / 2 #45方向
        S_Jplus0 = (s0[3] + s0[2]) / 2 #47方向
        # 各方向への厚み
        h_iminus = h_0 / ((poisson/(1-poisson) * (S_iminus - S_iminus0) / S_iminus0) + 1) #43方向
        h_Jminus = h_0 / ((poisson/(1-poisson) * (S_Jminus - S_Jminus0) / S_Jminus0) + 1) #41方向
        h_iplus = h_0 / ((poisson/(1-poisson) * (S_iplus - S_iplus0) / S_iplus0) + 1) #45方向
        h_Jplus = h_0 / ((poisson/(1-poisson) * (S_Jplus - S_Jplus0) / S_Jplus0) + 1) #47方向
        # 各断片の重心
        g = []
        kado = [2,4,6,8]
        hen = [1,3,5,7]
        for i in range(len(kado)):
            _kado = kado[i]
            _hen1, _ = [idx for idx in idx_iter.get_indexes_reverse(_kado, 2)]
            _hen2, _ = [idx for idx in idx_iter.get_indexes(_kado, 2)]
            _hen = [_hen1, _hen2]
            _g1 = (center_pos + _Ps[_kado] + _Ps[_hen1])/3
            _g2 = (center_pos + _Ps[_kado] + _Ps[_hen2])/3
            g.append([_g1, _g2])

        g401 = (center_pos + _Ps[0] + _Ps[1]) / 3
        g430 = (center_pos + _Ps[3] + _Ps[0]) / 3
        g436 = (center_pos + _Ps[3] + _Ps[6]) / 3
        g467 = (center_pos + _Ps[6] + _Ps[7]) / 3
        g478 = (center_pos + _Ps[7] + _Ps[8]) / 3
        g485 = (center_pos + _Ps[8] + _Ps[5]) / 3
        g452 = (center_pos + _Ps[5] + _Ps[2]) / 3
        g421 = (center_pos + _Ps[2] + _Ps[1]) / 3
        g0401 = (_Ps0[4] + _Ps0[0] + _Ps0[1]) / 3
        g0430 = (_Ps0[4] + _Ps0[3] + _Ps0[0]) / 3
        g0436 = (_Ps0[4] + _Ps0[3] + _Ps0[6]) / 3
        g0467 = (_Ps0[4] + _Ps0[6] + _Ps0[7]) / 3
        g0478 = (_Ps0[4] + _Ps0[7] + _Ps0[8]) / 3
        g0485 = (_Ps0[4] + _Ps0[8] + _Ps0[5]) / 3
        g0452 = (_Ps0[4] + _Ps0[5] + _Ps0[2]) / 3
        g0421 = (_Ps0[4] + _Ps0[2] + _Ps0[1]) / 3
        
        # 各断片面積
        triangle_area = []
        kado = [2,4,6,8]
        for i in range(len(kado)):
            j, k = [idx for idx in idx_iter.get_indexes_reverse(kado[i], 1)]
            _s1 = LA.norm(np.cross(x[j],x[k]))/2
            j, k = [idx for idx in idx_iter.get_indexes(kado[i], 1)]
            _s2 = LA.norm(np.cross(x[j],x[k]))/2
            triangle_area.append([_s1, _s2])

        s410 = LA.norm(np.cross(x[1],x[2]))/2
        s403 = LA.norm(np.cross(x[2],x[3]))/2
        s436 = LA.norm(np.cross(x[3],x[4]))/2
        s467 = LA.norm(np.cross(x[4],x[5]))/2
        s478 = LA.norm(np.cross(x[5],x[6]))/2
        s485 = LA.norm(np.cross(x[6],x[7]))/2
        s452 = LA.norm(np.cross(x[7],x[8]))/2
        s421 = LA.norm(np.cross(x[8],x[1]))/2
        s0410 = LA.norm(np.cross(x0[1],x0[2]))/2
        s0403 = LA.norm(np.cross(x0[2],x0[3]))/2
        s0436 = LA.norm(np.cross(x0[3],x0[4]))/2
        s0467 = LA.norm(np.cross(x0[4],x0[5]))/2
        s0478 = LA.norm(np.cross(x0[5],x0[6]))/2
        s0485 = LA.norm(np.cross(x0[6],x0[7]))/2
        s0452 = LA.norm(np.cross(x0[7],x0[8]))/2
        s0421 = LA.norm(np.cross(x0[8],x0[1]))/2
        # 四角の重心

        center_g_square = []
        for i in range(len(g)):
            _g = (triangle_area[i][0]*g[i][0] + triangle_area[i][1]*g[i][1])/(triangle_area[i][0] + triangle_area[i][1])
            center_g.append(_g)
        g4103 = (s410*g401 + s403*g430) / (s410 + s403)
        g4367 = (s436*g436 + s467*g467) / (s436 + s467)
        g4785 = (s478*g478 + s485*g485) / (s478 + s485)
        g4521 = (s452*g452 + s421*g421) / (s452 + s421)
        g04103 = (s0410*g0401 + s0403*g0430) / (s0410 + s0403)
        g04367 = (s0436*g0436 + s0467*g0467) / (s0436 + s0467)
        g04785 = (s0478*g0478 + s0485*g0485) / (s0478 + s0485)
        g04521 = (s0452*g0452 + s0421*g0421) / (s0452 + s0421)
        # 各重心間の距離
        Lj82 = LA.norm(g4521 - g4103)
        Lj24 = LA.norm(g4103 - g4367)
        Lj46 = LA.norm(g4367 - g4785)
        Lj68 = LA.norm(g4785 - g4521)
        
        # ひずみ
        eps_i41 = (LA.norm(x01) - LA.norm(x041)) / LA.norm(x041)
        eps_J41 = (LA.norm(g4521 - g4103) - LA.norm(g04521 - g04103)) / LA.norm(g04521 - g04103)
        eps_i43 = (LA.norm(x03) - LA.norm(x043)) / LA.norm(x043)
        eps_J43 = (LA.norm(g4103 - g4367) - LA.norm(g04103 - g04367)) / LA.norm(g04103 - g04367)
        eps_i47 = (LA.norm(x01) - LA.norm(x041)) / LA.norm(x041)
        eps_J47 = (LA.norm(g4367 - g4785) - LA.norm(g04367 - g04785)) / LA.norm(g04367 - g04785)
        eps_i45 = (LA.norm(x01) - LA.norm(x041)) / LA.norm(x041)
        eps_J45 = (LA.norm(g4785 - g4521) - LA.norm(g04785 - g04521)) / LA.norm(g04785 - g04521)
        # 張力
        F_T1 = (young_modulus * h_Jminus * Lj82 * (eps_i41 + poisson * eps_J41) / (1 - poisson**2))*x01/LA.norm(x01)
        F_T3 = (young_modulus * h_iminus * Lj24 * (eps_i43 + poisson * eps_J43) / (1 - poisson**2))*x03/LA.norm(x03)
        F_T5 = (young_modulus * h_Jplus * Lj46 * (eps_i47 + poisson * eps_J47) / (1 - poisson**2))*x05/LA.norm(x05)
        F_T7 = (young_modulus * h_iplus * Lj68 * (eps_i45 + poisson * eps_J45) / (1 - poisson**2))*x07/LA.norm(x07)
        # せん断ひずみ
        gamma513 = (math.acos((np.dot(x07,x01))/(LA.norm(x07)*LA.norm(x01))) - math.acos((np.dot(x045,x041))/(LA.norm(x045)*LA.norm(x041)))\
            + math.acos((np.dot(x03,x01))/(LA.norm(x03)*LA.norm(x01))) - math.acos((np.dot(x043,x041))/(LA.norm(x043)*LA.norm(x041))))/2
        gamma137 = (math.acos((np.dot(x01,x03))/(LA.norm(x01)*LA.norm(x03))) - math.acos((np.dot(x041,x043))/(LA.norm(x041)*LA.norm(x043)))\
            + math.acos((np.dot(x03,x05))/(LA.norm(x03)*LA.norm(x05))) - math.acos((np.dot(x043,x047))/(LA.norm(x043)*LA.norm(x047))))/2
        gamma375 = (math.acos((np.dot(x05,x03))/(LA.norm(x05)*LA.norm(x03))) - math.acos((np.dot(x047,x043))/(LA.norm(x047)*LA.norm(x043)))\
            + math.acos((np.dot(x07,x05))/(LA.norm(x07)*LA.norm(x05))) - math.acos((np.dot(x045,x047))/(LA.norm(x045)*LA.norm(x047))))/2
        gamma751 = (math.acos((np.dot(x05,x07))/(LA.norm(x05)*LA.norm(x07))) - math.acos((np.dot(x047,x045))/(LA.norm(x047)*LA.norm(x045)))\
            + math.acos((np.dot(x07,x01))/(LA.norm(x07)*LA.norm(x01))) - math.acos((np.dot(x045,x041))/(LA.norm(x045)*LA.norm(x041))))/2
        # せん断力
        F_S41 = ((young_modulus * h_Jminus * LA.norm(x01) * gamma513)/(2 * (1 + poisson)))*x01/LA.norm(x01)
        F_S43 = ((young_modulus * h_Jminus * LA.norm(x03) * gamma137)/(2 * (1 + poisson)))*x03/LA.norm(x03)
        F_S47 = ((young_modulus * h_Jminus * LA.norm(x05) * gamma375)/(2 * (1 + poisson)))*x05/LA.norm(x05)
        F_S45 = ((young_modulus * h_Jminus * LA.norm(x07) * gamma751)/(2 * (1 + poisson)))*x07/LA.norm(x07)
        
        # J方向の曲げ力
        n_j_cross = np.cross(x05, x01)
        if any(n_j_cross):
            n_J = n_j_cross/LA.norm(n_j_cross)
        else: 

        l_Jalfa = LA.norm(_Ps[1] - _Ps[7])
        cos_Jalfa = (LA.norm(x01)**2 + LA.norm(x05)**2 - l_Jalfa**2) / (2 * LA.norm(x01) * LA.norm(x05))
        if cos_Jalfa > 1.0:
            cos_Jalfa = 1.0
        elif cos_Jalfa < -1.0:
            cos_Jalfa = -1.0
        sin_Jalfa = math.sqrt(1 - cos_Jalfa**2)
        CJa2 = math.sqrt((cos_Jalfa + 1)/2)
        SJa2 = math.sqrt((1 - cos_Jalfa)/2)
        zJC = (_Ps[7][2]-_Ps[1][2])/(_Ps[7][0]-_Ps[1][0]) * (center_pos[0]-_Ps[1][0]) + _Ps[1][2] #曲げ力の方向の場合わけに必要
        if center_pos[2] > zJC:
            e_j = np.dot(np.array([[CJa2 + (n_J[0]**2) * (1 - CJa2), n_J[0] * n_J[1] * (1 - CJa2) + n_J[2] * SJa2, n_J[0] * n_J[2] * (1 - CJa2) - n_J[1] * SJa2],\
                                    [n_J[1] * n_J[0] * (1 - CJa2) - n_J[2] * SJa2, CJa2 + (n_J[1]**2) * (1 - CJa2), n_J[1] * n_J[2] * (1 - CJa2) + n_J[0] * SJa2],\
                                    [n_J[2] * n_J[0] * (1 - CJa2) + n_J[1] * SJa2, n_J[2] * n_J[1] * (1 - CJa2) - n_J[0] * SJa2, CJa2 + (n_J[2]**2) * (1 - CJa2)]]), (_Ps[7] - center_pos)/LA.norm(_Ps[7] - center_pos))
        else:
            e_j = np.dot(np.array([[CJa2 + (n_J[0]**2) * (1 - CJa2), n_J[0] * n_J[1] * (1 - CJa2) - n_J[2] * SJa2, n_J[0] * n_J[2] * (1 - CJa2) + n_J[1] * SJa2],\
                                    [n_J[1] * n_J[0] * (1 - CJa2) + n_J[2] * SJa2, CJa2 + (n_J[1]**2) * (1 - CJa2), n_J[1] * n_J[2] * (1 - CJa2) - n_J[0] * SJa2],\
                                    [n_J[2] * n_J[0] * (1 - CJa2) - n_J[1] * SJa2, n_J[2] * n_J[1] * (1 - CJa2) + n_J[0] * SJa2, CJa2 + (n_J[2]**2) * (1 - CJa2)]]), (_Ps[7] - center_pos)/LA.norm(_Ps[7] - center_pos))
        d_etha_J = (2 * sin_Jalfa / l_Jalfa) - (2 * math.sqrt(1 - np.dot(x041,x047)**2/(LA.norm(x041)*LA.norm(x047))**2)/(LA.norm(x041 - x047)))

        n_i = np.cross(x07,x03)/LA.norm(np.cross(x03,x07))  
        cos_ialfa = np.dot(x03,x07) / (LA.norm(x03) * LA.norm(x07))
        sin_ialfa = math.sqrt(1 - cos_ialfa**2)
        Cia2 = math.sqrt((cos_ialfa + 1)/2)
        Sia2 = math.sqrt((1 - cos_ialfa)/2)
        ziC = (_Ps[5][2]-_Ps[3][2])/(_Ps[5][0]-_Ps[3][0]) * (center_pos[0]-_Ps[3][0]) + _Ps[3][2]
        if center_pos[2] > ziC:
            e_i = np.dot(np.array([[Cia2 + (n_i[0]**2) * (1 - Cia2), n_i[0] * n_i[1] * (1 - Cia2) + n_i[2] * Sia2, n_i[0] * n_i[2] * (1 - Cia2) - n_i[1] * Sia2],\
                                    [n_i[1] * n_i[0] * (1 - Cia2) - n_i[2] * Sia2, Cia2 + (n_i[1]**2) * (1 - Cia2), n_i[1] * n_i[2] * (1 - Cia2) + n_i[0] * Sia2],\
                                    [n_i[2] * n_i[0] * (1 - Cia2) + n_i[1] * Sia2, n_i[2] * n_i[1] * (1 - Cia2) - n_i[0] * Sia2, Cia2 + (n_i[2]**2) * (1 - Cia2)]]), (_Ps[7] - center_pos)/LA.norm(_Ps[7] - center_pos))
        else:
            e_i = np.dot(np.array([[Cia2 + (n_i[0]**2) * (1 - Cia2), n_i[0] * n_i[1] * (1 - Cia2) - n_i[2] * Sia2, n_i[0] * n_i[2] * (1 - Cia2) + n_i[1] * Sia2],\
                                    [n_i[1] * n_i[0] * (1 - Cia2) + n_i[2] * Sia2, Cia2 + (n_i[1]**2) * (1 - Cia2), n_i[1] * n_i[2] * (1 - Cia2) - n_i[0] * Sia2],\
                                    [n_i[2] * n_i[0] * (1 - Cia2) - n_i[1] * Sia2, n_i[2] * n_i[1] * (1 - Cia2) + n_i[0] * Sia2, Cia2 + (n_i[2]**2) * (1 - Cia2)]]), (_Ps[5] - center_pos)/LA.norm(_Ps[5] - center_pos))
        d_etha_i = (2 * sin_ialfa / LA.norm(x07 - x03)) - (2 * math.sqrt(1 - np.dot(x043,x045)**2/(LA.norm(x043)*LA.norm(x045))**2)/(LA.norm(x043 - x045)))


        l_J = (Lj20 + Lj06 + Lj68 + Lj82) / 4
        h = (h_iminus + h_iplus + h_Jminus + h_Jplus) / 4
        I = (l_J * h**3) / 12
        M_i = (young_modulus * I * (d_etha_i + poisson * d_etha_J)/(1 - poisson**2))
        M_J = (young_modulus * I * (d_etha_J + poisson * d_etha_i)/(1 - poisson**2))
        #曲げ力
        F_Bi = M_i / LA.norm(x03) + M_i / LA.norm(x07) * e_i
        F_BJ = M_J / LA.norm(x01) + M_J / LA.norm(x05) * e_j
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
            """ 
            [point positions]
                8 - 1 - 2
                |   |   |
                7 - 0 - 3
                |   |   |
                6 - 5 - 4
            """
            acc = acce[i][1] # 粒子0の加速度
            v2 = velo[i] # 粒子0の速度
            p0 = particles[i][1] # 粒子0(中心
            p1 = particles[i-11][1] # 粒子1（12時方向
            p2 = particles[i-12][1] # 粒子2(1時半方向
            p3 = particles[i-1][1] # 粒子3（3時方向
            p4 = particles[i+10][1] # 粒子4(4時半
            p5 = particles[i+11][1]　# 粒子5(6時
            p6 = particles[i+12][1] # 粒子6(7時半
            p7 = particles[i+1][1] # 粒子7(9時
            p8 = particles[i-10][1] # 粒子8(10時半
            Ps = np.array([p0, p1, p2, p3, p4, p5, p6, p7, p8]) # 対象粒子、隣り合う粒子座標のセットリスト
            p0_0 = particles_0[i][1] # 粒子0(中心
            p1_0 = particles_0[i-11][1] # 粒子1（12時方向
            p2_0 = particles_0[i-12][1] # 粒子2(1時半方向
            p3_0 = particles_0[i-1][1] # 粒子3（3時方向
            p4_0 = particles_0[i+10][1] # 粒子4(4時半
            p5_0 = particles_0[i+11][1]　# 粒子5(6時
            p6_0 = particles_0[i+12][1] # 粒子6(7時半
            p7_0 = particles_0[i+1][1] # 粒子7(9時
            p8_0 = particles_0[i-10][1] # 粒子8(10時半 
            Ps0 = np.array([p0_0, p1_0, p2_0, p3_0, p4_0, p5_0, p6_0, p7_0, p8_0]) # 初期座標
        

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
    history = []

    fig = plt.figure()
    ims = []
    ax = fig.add_subplot(111, projection='3d')

    particles, particles_0, forces, velo = prepare()
    acce = init_acceleration()
    while time <= 0.01:

        X = []
        Y = []
        Z = []
        if time%10.0 <=0.01:
            history.append([time, particles[1][1][0], particles[1][1][1]])
            
        calculation()
        time += dt
    for i, P in enumerate(particles):

        X.append(P[1][0])
        Y.append(P[1][1])
        Z.append(P[1][2])
    ims = []
    im = ax.scatter(X,Y,Z,c='c')
    np.savetxt("time_history.csv", history, delimiter=",")
    plt.show()