import math
import random
import numpy as np
import matplotlib.pyplot as plt
import copy
beta = 0.4
J = 1
num_size = 128



def de_sigma(num_size):
    # 确定初始spin序列
    sigma = np.random.rand(num_size)
    for i in range(num_size):

        if sigma[i] > 0.5:
            sigma[i] = 1
        else:
            sigma[i] = -1

    sigma = np.array(sigma)
    sigma[num_size - 1] = sigma[0]  # periodic boundary condition
    print(sum(sigma))

    return sigma


sigma = de_sigma(num_size)


# print(sigma)


def cal_energy(sequence):
    # 计算总能量
    energy = 0
    for i in range(num_size):

        bond_energy = -J * sequence[i] * sequence[(i + 1) % num_size]
        energy = energy + bond_energy
    return energy





def energy_difference(sequence, site):
    difference = 2 * J * sequence[site] * (
            sequence[(site + 1) % num_size] + sequence[(site- 1) % num_size])  # deltaE=2j*(onsite)*(left +right)
    return difference


#
def metro(initial):
    site=np.random.randint(num_size)
    #print(site)
    delta_E = energy_difference(initial, site)
    p0 = math.exp(-beta * delta_E)
    prob = min(1.0, p0)
    u = np.random.random()  # u[0]为随机数（0，1），x/(0,1)为一个随机变量，p(x<u[0])=u[0]
    if prob > u:  # 有prob的概率发生某件事
        initial[site] = -initial[site]  # flip
       # print('yeas','prob:',prob,' u',u)
    else:
        initial[site] = initial[site]
    return initial
E_aver_thermal = []


steps=10000000
storage = []

sigma_used = sigma
i = 0
bin=200
while i < steps:
    # print('sig', sigma_used)

    newpat = metro(sigma_used)
   # print('new', newpat)
    new_pat_use=copy.copy(newpat)

    storage.append(new_pat_use)  #python 对象引用容易出bug

   # print('storage',storage)
    i = i + 1
    sigma_used = newpat

storage_sample=storage[::200]
energy_list=[cal_energy(i) for i in storage_sample]
mean_energy_list=[]
step_list=[]
print(len(energy_list))
for i in range(len(energy_list)):
    interval=400
    steps_0=i
    print(steps_0)
    # energy_dots=energy_list[int(steps_0):int(steps_0)+1]   #取某个step数目附近1000个样本点
    # mean_energy=np.mean(energy_dots)               #求这1000个样本点的平均值
    mean_energy_list.append(np.mean(energy_list[:i]))
    step_list.append(steps_0)

# Etot = 0
# #接下来看能量随steps变化
# #bin
# bin=1000
# E_list_bin=[sigma]
# #print('storage',storage)
# step_list=[]
# #print(storage[:100])
#
# print('start')
# for i in range(400):
#     steps_0=i*80000
#
#     E_list_bin=storage[:steps_0+1:bin]
#     energy_list=[cal_energy(i) for i in E_list_bin]
#
#     mean_energy=np.mean(energy_list)
#     mean_energy=copy.copy(mean_energy)
#     E_aver_thermal.append(mean_energy)  #某个step下算出来的平均能量
#     x=copy.copy(steps_0)          #将该次步数放入x轴列表中
#     step_list.append(x)
#

plt.plot(step_list[100:], mean_energy_list[100:])
plt.xlabel('steps')
plt.ylabel('energy')
plt.show()


