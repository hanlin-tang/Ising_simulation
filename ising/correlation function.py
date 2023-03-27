import math
import random
import numpy as np
import matplotlib.pyplot as plt
import copy

J = 1
num_size =128
steps=10000000
bin = 500
beta=0.8
def de_sigma(num_size):
    # 确定初始spin序列
    sigma = np.random.rand(num_size)
    for i in range(num_size):

        if sigma[i] > 0.5:
            sigma[i] = 1
        else:
            sigma[i] = -1

    sigma = np.array(sigma)

    return sigma


sigma = de_sigma(num_size)


#print(sigma)


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


def metro(initial,beta):

    site=np.random.randint(num_size)
    #print(site)


    delta_E = energy_difference(initial, site)
    p0 = math.exp(-beta * delta_E)
    prob = min(1.0, p0)
    # #print(site+1)
    # print(initial)
    u = np.random.random()  # u[0]为随机数（0，1），x/(0,1)为一个随机变量，p(x<u[0])=u[0]

    if prob > u:  # 有prob的概率发生某件事
        initial[site] = -initial[site]  # flip
       # print('yeas','prob:',prob,' u',u)
    else:
        initial[site] = initial[site]
    #print(initial,'prob:',prob)
    #print()



    return initial

def calculate_correlation(sequence,n):
    corre=sequence[0]*sequence[n]
    return corre


storage = []
sigma_used = sigma
i = 0

while i < steps:
    # print('sig', sigma_used)

    newpat = metro(sigma_used,beta)
    # print('new', newpat)
    new_pat_use = copy.copy(newpat)

    storage.append(new_pat_use)  # python 对象引用容易出bug

    # print('storage',storage)
    i = i + 1
    sigma_used = newpat

storage_sample = storage[::bin]

storage_sample_minus_first=storage_sample[1500:]  #去除前面1000个数据点


def get_correlation(n):

   # print(storage_sample_minus_first)
    result_corre_list = [calculate_correlation(i,n) for i in storage_sample_minus_first]  #计算每个构型的关联
   # print('storage_sample_minus_first[1:]',storage_sample_minus_first[1:])
    #print('energy_list', energy_list)
    results_cor_mean=np.mean(result_corre_list)




    return results_cor_mean   #
# i=0
# alist=[]
# while i<10:
#     #print(get_energy(0.8)
#     i=i+1
#     print('i',i)
#     alist.append(get_energy(0.8))
# mean=np.mean(alist)
# print('mean',mean)
n_list=[n for n in range(40)]

cor_n_list=[get_correlation(n) for n in n_list]
theory_list=[math.tanh(beta)**i for i in n_list]
plt.plot(n_list,cor_n_list,'r+',label='simulation')
plt.plot(n_list,theory_list,'b-',label='theory')
plt.xlabel('beta')
plt.ylabel('correlation function')
plt.title('steps:{}'.format(steps))
plt.legend()          # show label
plt.show()


