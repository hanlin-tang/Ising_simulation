# -*-coding:utf-8-*-
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import copy



def de_sigma(num_size):
    # 确定初始spin序列, 用i, j, 标记二维 ising
    sigma_list=[]     #创健new list
    for i in range(num_size):
        sigma = np.random.rand(num_size)



        for j in range(num_size):

            if sigma[j] > 0.5:
                sigma[j] = 1
            else:
                sigma[j] = -1
        print('sig',sigma)
        sigma_list.append(sigma)




    return sigma_list




#apply cluster algri
def choo_neigh_(sig1, sig2):
    #P0=max(sig1*sig2,0)*0  #无论如何翻转改了都很小，高温极限，不相关
    P0=max(1-math.exp(-2*beta * J*sig1*sig2),0)
    #print('p0',P0,'SIGMA',sig1,'sig2',sig2)
    compare = np.random.random()

    if P0 > compare:
        a=0
    else:
        a=1
    return a    #a=0 加入，反之不加入

def sub_list(a,b):   #减去列表操作
    final=[i for i in a if i not in b ]
    return final

def find_four_bonds(i,j):    #找到周围的四个格点
    up = (i - 1) % num_size  # up, j
    down = (i + 1) % num_size  # down, j
    right = (j - 1) % num_size  # i right
    left = (j + 1) % num_size  # i left
    cluster_ready = [[up, j], [down, j], [i, right], [i, left]]
    #print('clusterreadu',cluster_ready)
    return cluster_ready

def get_tot_energy(sigma_list):
    E_lsit = []
    for i in range(num_size):
        for j in range(num_size):
            site = find_four_bonds(i,j)
            energy_site = sigma_list[i][j]*sigma_list[site[0][0]][site[0][1]]+\
                          sigma_list[i][j]*sigma_list[site[1][0]][site[1][1]]+\
                          sigma_list[i][j]*sigma_list[site[2][0]][site[2][1]]+\
                          sigma_list[i][j]*sigma_list[site[3][0]][site[3][1]]
            E_lsit.append(energy_site)
            tot_energy = np.sum(E_lsit)/2
    return tot_energy


def form_small_four_cluster(i,j,sigma_list, counted_list):  #输入格点，找四个bond，减去考虑过的, 再判断是否加入
    cluster_ready= find_four_bonds(i,j)    #找到这个格点附近的四个bond, 并在下一步减去之前考虑过的bond
    #print('slusready', cluster_ready, 'counted',counted_list)
    clus_with_no_dupli=sub_list(cluster_ready, counted_list)
    #print('clus_with_no_dupli',clus_with_no_dupli)

    useful_clus=[k for k in clus_with_no_dupli if choo_neigh_(sigma_list[k[0]][k[1]], sigma_list[i][j])==0]  #将符合要求的指标装入备选列表
    return useful_clus

def form_cluster(sigma_list,i0,j0):    #return cluster index list


   # print('ini', i0, j0)
    cluster_list=[[i0,j0]] #候选列表

    singletimecluster=[[i0,j0]]   # 初始选择的任意格点
    while 1:    #找到cluster
        singletime=[]
        #counter=[]  #初始化counter
        for i in singletimecluster:     # singletimecluster: 上一次加入的格点

            this_time_newlist = form_small_four_cluster(i[0], i[1], sigma_list, cluster_list) #新的可以加入的格点

            #this = copy.copy(this_time_newlist)    #重新分配地址
            singletime = singletime+copy.copy(this_time_newlist)  #加入这个格点四周允许的sites
        #counter1=copy.copy(counter)   # 这次考虑过的bond,无论是否有加入
       #counted_list=counted_list+copy.copy(singletime) #counter1   #所有考虑过的bond
        singletimecluster = np.unique(np.array(copy.copy(singletime)) ,axis = 0).tolist()
        #print('每次加入：',singletimecluster)
        cluster_list =np.unique(np.array(cluster_list + singletimecluster) ,axis = 0).tolist()#加入到cluster中

        #print('sig',singletimecluster)
        if len(singletimecluster) == 0:
             break
    #print('len of cluster',len(cluster_list))
    return cluster_list
# print('x',sigma_list,'i0',i0,'j0',j0)
# print([i for i in cluster_list])

def cal_ma_ideal(beta):
    a = np.sinh(2*beta)
    m = 1-(a)**(-4)
    m = m**(1/8)
    return m



beta_list = [0.03+0.02*i for i in range(30)]
xx_list = [0.001*i for i in range(600)]

yy_list = []
for i in xx_list:
    if i < 0.442:
        yy_list.append(0)
    else:
        yy_list.append(cal_ma_ideal(i))



num_list=[7,8,17]
J = 1
for kk in num_list:
    num_size = kk
    sigma_list=de_sigma(num_size)
    ulist=[]
    maglist = []
    cv_list = []
    #sus_list = []
    mean_energy_list = []
    for beta in beta_list:
        steps= 10000
        #print('beta', beta)
        mlist=[]
        EE_list = []

        for j in range(steps):           #开始集团更新
            #print(j)

            i0 = random.choice(range(num_size))
            j0 = random.choice(range(num_size))

            clusterlist = copy.copy(form_cluster(sigma_list,i0,j0))
            value=copy.copy(sigma_list[i0][j0])

            for i in clusterlist:   #反转cluster中的每一个site
                #print('start',sigma_list[i[0]][i[1]])
                sigma_list[i[0]][i[1]]=-value
            #for i in sigma_list:
                #print('fliped',i)
               # print('end',sigma_list[i[0]][i[1]])
            E_tot = get_tot_energy(sigma_list)
            m=sum(np.array(sigma_list).flatten())
            EE_list.append(E_tot)
            mlist.append(abs(m))
        mlist=mlist[5000:]
        EE_list = EE_list[5000:]
        #print('end')
        #print(mlist)
        mag = np.mean(mlist)
        mean_energy = np.mean(EE_list)
        energy_suq_mean = np.mean([i*i for i in EE_list])


        maglist.append(mag/kk**2)
        #print(maglist)
        mean_energy_list.append(mean_energy)
        cv_list.append(beta * beta * (energy_suq_mean - mean_energy ** 2))
        #sus_list.append((energy_suq_mean - mean_energy ** 2)*beta)


        #plt.hist(mlist,bins=20)
        maggg = [i / (num_size ** 2) for i in mlist]
        magcube=[i**4 for i in maggg]
        magsqure=[i*i for i in maggg]

        u=1-(np.mean(magcube))/(3*(np.mean(magsqure))**2)
        ulist.append(u)

    plt.plot(beta_list, ulist,'*', label = num_size)
    plt.legend()
    plt.xlabel('beta')
    plt.ylabel('u')
plt.show()