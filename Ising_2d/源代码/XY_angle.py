# -*-coding:utf-8-*-
# -*-coding:utf-8-*-
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import copy
from copy import deepcopy



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

def gen_xy(num):   #create original list
    s_list = []
    for i in range(num):
        sigma = np.random.rand(num)

        s_list.append(sigma*2*3.14)
    return s_list



#apply cluster algri
def choo_neigh_(sig1, sig2,ref_angle1):
    #P0=max(sig1*sig2,0)*0  #无论如何翻转改了都很小，高温极限，不相关
    #P0=max(1-math.exp(-2*beta * J*sig1*sig2),0)
    P0 = max(1-math.exp(-2*beta*J*np.cos(sig1-ref_angle1)*np.cos(sig2-ref_angle1)),0)
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
            energy_site = np.cos(sigma_list[i][j]-sigma_list[site[0][0]][site[0][1]])+\
                          np.cos(sigma_list[i][j]-sigma_list[site[1][0]][site[1][1]])+\
                          np.cos(sigma_list[i][j]-sigma_list[site[2][0]][site[2][1]])+\
                          np.cos(sigma_list[i][j]-sigma_list[site[3][0]][site[3][1]])
            E_lsit.append(energy_site)
            tot_energy = np.sum(E_lsit)/2
    return tot_energy

def get_energy_config(sigma_list):
    E_site = deepcopy(sigma_list)
    for i in range(num_size):
        for j in range(num_size):
            site = find_four_bonds(i, j)
            energy_site = np.cos(sigma_list[i][j] - sigma_list[site[0][0]][site[0][1]]) + \
                          np.cos(sigma_list[i][j] - sigma_list[site[1][0]][site[1][1]]) + \
                          np.cos(sigma_list[i][j] - sigma_list[site[2][0]][site[2][1]]) + \
                          np.cos(sigma_list[i][j] - sigma_list[site[3][0]][site[3][1]])
            E_site[i][j] = energy_site

    return E_site

def form_small_four_cluster(i,j,sigma_list, counted_list,ref):  #输入格点，找四个bond，减去考虑过的, 再判断是否加入
    cluster_ready= find_four_bonds(i,j)    #找到这个格点附近的四个bond, 并在下一步减去之前考虑过的bond
    #print('slusready', cluster_ready, 'counted',counted_list)
    clus_with_no_dupli=sub_list(cluster_ready, counted_list)
    #print('clus_with_no_dupli',clus_with_no_dupli)

    useful_clus=[k for k in clus_with_no_dupli if choo_neigh_(sigma_list[k[0]][k[1]], sigma_list[i][j],ref)==0]  #将符合要求的指标装入备选列表
    return useful_clus

def form_cluster(sigma_list,i0,j0,ref0):    #return cluster index list


   # print('ini', i0, j0)
    cluster_list=[[i0,j0]] #候选列表

    singletimecluster=[[i0,j0]]   # 初始选择的任意格点
    while 1:    #找到cluster
        singletime=[]
        #counter=[]  #初始化counter
        for i in singletimecluster:     # singletimecluster: 上一次加入的格点

            this_time_newlist = form_small_four_cluster(i[0], i[1], sigma_list, cluster_list,ref0) #新的可以加入的格点

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


def angle_conv(a):
    b = (a%6.28)/3.14*180
    return b
beta_list = [111.89]





num_list=[50]
x = np.arange(0,num_list[0],1)
y = np.arange(0,num_list[0],1)

xx, yy = np.meshgrid(x, y, sparse=True)
J = 1
for kk in num_list:
    num_size = kk
    sigma_list=gen_xy(num_size)     #generate new config
    original = deepcopy(sigma_list)

    ulist=[]
    maglist = []
    cv_list = []
    #sus_list = []
    mean_energy_list = []
    for beta in beta_list:
        steps= 50
        #print('beta', beta)
        mlist=[]
        EE_list = []
        step_list =[]
        for j in range(steps):           #开始集团更新

            step_list.append(j)

            i0 = random.choice(range(num_size))
            j0 = random.choice(range(num_size))

            ref_angle = 2*3.14159*random.random()     #random generating ref angle


            clusterlist = copy.copy(form_cluster(sigma_list,i0,j0,ref_angle))
            #print(len(clusterlist))
            #value=copy.copy(sigma_list[i0][j0])
            a = random.random()
            if a>0:

                for i in clusterlist:   #反转cluster中的每一个site
                    value = copy.copy(sigma_list[i[0]][i[1]])
                    #print('start',sigma_list[i[0]][i[1]])
                    #sigma_list[i[0]][i[1]]=-value           # flip
                    #print('ref',angle_conv(ref_angle),angle_conv(sigma_list[i[0]][i[1]]))
                    sigma_list[i[0]][i[1]] = 3.14159-value+2*ref_angle
                # print('ref', angle_conv(ref_angle),angle_conv(sigma_list[i[0]][i[1]]))
                # print()
                #for i in sigma_list:
                #print('fliped',i)
               # print('end',sigma_list[i[0]][i[1]])
            # E_tot = get_tot_energy(sigma_list)
            #
            # m=sum(np.cos(np.array(sigma_list).flatten()))
            # EE_list.append(E_tot/kk**2)
            # mlist.append(abs(m))

        #plt.plot(np.array(step_list),mlist, label = 'steps ')
        # mlist=mlist[800:]
        # EE_list = EE_list[800:]    #energy ensemble at certain temperature

        #print('end')
        #print(mlist)
        # mag = np.mean(mlist)
        # mean_energy = np.mean(EE_list)
        # energy_suq_mean = np.mean([i*i for i in EE_list])


        #maglist.append(mag/kk**2)
        #print(maglist)
        #mean_energy_list.append(mean_energy)
        #cv_list.append(beta * beta * (energy_suq_mean - mean_energy ** 2))
        #sus_list.append((energy_suq_mean - mean_energy ** 2)*beta)


        #plt.hist(mlist,bins=20)
        # maggg = [i / (num_size ** 2) for i in mlist]
        # magcube=[i**4 for i in maggg]
        # magsqure=[i*i for i in maggg]

        #u=1-(np.mean(magcube))/(3*(np.mean(magsqure))**2)
        # ulist.append(u)

    # plt.plot(beta_list, cv_list,'*', label = 'lattice = {}'.format(num_size))
    # plt.xlabel('beta')
    # plt.ylabel('Cv')
    # plt.legend()
    fig, ax = plt.subplots(figsize=(4, 4))
    angle_list = np.array(sigma_list)
    plt.imshow((angle_list % 6.28) / 6.28 * 360, cmap='twilight_shifted', interpolation='nearest', origin='upper')
    cbar = plt.colorbar()

    # 设置colorbar范围
    cbar.mappable.set_clim(0, 360)
    fig2, ax2 = plt.subplots(figsize=(4, 4))
    xx, yy = np.meshgrid(x, y, sparse=True)
    u = np.cos(sigma_list)
    v = np.sin(sigma_list)
    #print('sig',u)
    u1 = np.cos(original)
    v1 = np.sin(original)
    energy_config = get_energy_config(sigma_list)


    plt.imshow((np.array(original) % 6.28) / 6.28 * 360, cmap='twilight_shifted', interpolation='nearest',
               origin='upper')
    cbar = plt.colorbar()

    #设置colorbar范围
    cbar.mappable.set_clim(0,360)


    ax2.quiver(xx, yy, u1, v1, color = 'blue')
    ax.quiver(xx, yy, u, v)
plt.show()