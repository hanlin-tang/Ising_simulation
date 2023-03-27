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

def form_small_four_cluster(i,j,sigma_list):  #输入格点，找四个bond，减去考虑过的, 再判断是否加入
    cluster_ready= find_four_bonds(i,j)    #找到这个格点附近的四个bond, 并在下一步减去之前考虑过的bond
    #print('slusready', cluster_ready, 'counted',counted_list)
 #   clus_with_no_dupli=sub_list(cluster_ready, counted_list)
    #print('clus_with_no_dupli',clus_with_no_dupli)

    useful_clus=[k for k in cluster_ready if choo_neigh_(sigma_list[k[0]][k[1]], sigma_list[i][j])==0]  #将符合要求的指标装入备选列表
    return useful_clus

def form_cluster(sigma_list,i0,j0):    #return cluster index list

    counted_list=[] #所有考虑过的bond
   # print('ini', i0, j0)
    cluster_list=[[i0,j0]] #候选列表

    singletimecluster=[[i0,j0]]   # 初始选择的任意格点
    while 1:    #找到cluster
        singletime=[]
        counter=[]  #初始化counter
        for i in singletimecluster:     # singletimecluster: 上一次加入的格点

            this_time_newlist = form_small_four_cluster(i[0], i[1], sigma_list) # 新的可以加入的格点
            #print('this',this_time_newlist)
            this_time_newlist = sub_list(this_time_newlist, counted_list)
            counter=counter+find_four_bonds(i[0],i[1])   # 这个格点对应的所有bond
            this = copy.copy(this_time_newlist)     # 重新分配地址
            singletime = singletime+this   # 加入这个格点四周允许的sites
        counter1=copy.copy(counter)   # 这次考虑过的bond,无论是否有加入
        counted_list=counted_list+counter1   #所有考虑过的bond
        singletimecluster = copy.copy(singletime)
        #print('每次加入：',singletimecluster)
        cluster_list = cluster_list + singletimecluster #加入到cluster中
        #print('sig',singletimecluster)
        if len(singletimecluster) == 0:
             break
    #print('len of cluster',len(cluster_list))
    return cluster_list
# print('x',sigma_list,'i0',i0,'j0',j0)
# print([i for i in cluster_list])



beta_list = [0.27+i*0.5/30 for i in range(30)]   # 构建温度序列

num_list=[5,6,8]
J = 1
for kk in num_list:
    num_size = kk



    sigma_list=de_sigma(num_size)

    y_list = []
    x_list = []

    for beta in beta_list:    # 遍历温度list
        print(beta)
        T = 1 / beta
        t = (T/2.09-1)  # 约化温度
        if beta<0.44:     #高低温热化时间不同，cluster能不能长起来的问题
            steps=30000
        else:
            steps=800

        # if beta < 0.4:
        #
        #     steps= 40000
        # elif 0.4<beta<0.5:
        #     steps= 20000
        # else:
        #     steps= 5000

        print('beta', beta)
        mlist = []


        for j in range(steps):

            #print(j)

            i0 = random.choice(range(num_size))
            j0 = random.choice(range(num_size))

            clusterlist = copy.copy(form_cluster(sigma_list,i0,j0))
            value=copy.copy(sigma_list[i0][j0])

            for i in clusterlist:   #反转cluster中的每一个site

                sigma_list[i[0]][i[1]]=-value

            m=sum(np.array(sigma_list).flatten())

            mlist.append(abs(m)/num_size**2)     #Z2对称性，记得加绝对值！！！

        if beta<0.44:
            mlist = mlist[8000:]  #int(steps/10)  去掉没热化的部分
        else:
            mlist = mlist[50:]
        average_m = sum(mlist)/len(mlist)
        y = abs(average_m) * num_size ** (1 / 8)
        y_list.append(y) # 通过磁化计算y

        x = t * num_size  # 通过调节温度来调x
        print(x)
        x_list.append(x)
        # print('x', x_list)
        # print('y', y_list)
    tlist=[1/i for i in beta_list]
    if num_size==8:
        plt.plot( x_list, y_list,'+',label='L = {}'.format(num_size))
    else:
        plt.plot(x_list, y_list, '-', label='L = {}'.format(num_size))
    plt.legend()
    plt.xlabel('X')
    plt.ylabel('Y')
plt.show()


















