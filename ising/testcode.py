import math
import random
import copy
J=1
a=[-1.,  1.,  1., -1., -1., -1., -1.,  1.,  1.,  1.,  1.,  1.,  1.,
       -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1.,
       -1., -1., -1., -1., -1., -1.,  1.,  1.,  1.,  1.,  1.,  1., -1.,
       -1., -1.,  1.,  1.,  1., -1., -1., -1., -1., -1., -1., -1., -1.,
       -1., -1., -1.,  1.,  1.,  1., -1., -1., -1., -1., -1.,  1.,  1.,
        1.,  1., -1., -1., -1.,  1.,  1.,  1.,  1.,  1., -1., -1., -1.,
       -1., -1., -1., -1., -1., -1., -1., -1., -1., -1.,  1.,  1.,  1.,
        1.,  1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1.,
        1.,  1.,  1.,  1.,  1., -1., -1.,  1.,  1.,  1.,  1., -1., -1.,
        1.,  1., -1., -1., -1.,  1.,  1.,  1.,  1., -1., -1.]
def cal_energy(sequence):
    # 计算总能量
    energy = 0
    for i in range(len(a)):

        bond_energy = -J * sequence[i] * sequence[(i + 1) % len(a)]
        #print(bond_energy)
        energy = energy + bond_energy
    return energy
e=cal_energy(a)
num_size=10
import numpy as np
def metro(initial,beta):

    site=np.random.randint(num_size)
    #print(site)
    p0 = beta
    # prob = min(1.0, p0)
    # #print(site+1)
    # print(initial)
    print()
    print('initial',initial)
    u = np.random.random()  # u[0]为随机数（0，1），x/(0,1)为一个随机变量，p(x<u[0])=u[0]
    final=copy.copy(initial)
    print('initial01', initial)
    final[site] = -initial[site]   #python赋值后的列表中的元素在内存中还是引用的原来的变量名，改了这个，initial会跟着变
    print('initial02', initial)
    if p0 > u:  # 有prob的概率发生某件事
        final=final  # flip
       # print('yeas','prob:',prob,' u',u)
        print('yes')
    else:
        final=initial
        print('no')
    print('final',final,'initial',initial)
    #print(initial,'prob:',prob)
    #print()



    return final
alist=[]
for i in range(1000):
    alist.append(np.mean(metro([1,1,1,1,1,1,1,1,1,1],0.5)))
print(np.mean(alist))
print(np.mean([1,1,1,1,1,1,1,1,1,-1]))

print()
a=[1,2,3,4]
b=[6,7,8,9]
a[0]=-b[0]
print('a',a,'b',b)
x=[1,2,3,4]
y=[2,3,4,5]
import matplotlib.pyplot as plt
plt.plot(x,y,'r+',label='sdfsdfsdf',)
plt.plot(x,y,'b-',label='sdfddddddsdfsdf')
plt.legend()
plt.show()

print(math.asinh(4))
print(1/math.sinh(4))