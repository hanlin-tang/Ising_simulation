import math
import random
import numpy as np
import matplotlib.pyplot as plt
import copy



def de_sigma(num_size):
    # confirm initial spin lattice, label 2d ising by i,j
    sigma_list=[]     # create new list
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

    P0=max(1-math.exp(-2*beta * J*sig1*sig2),0)
    #print('p0',P0,'SIGMA',sig1,'sig2',sig2)
    compare = np.random.random()

    if P0 > compare:
        a=0
    else:
        a=1
    return a    #a=0 add vice versa

def sub_list(a,b):   #list substraction
    final=[i for i in a if i not in b ]
    return final

def find_four_bonds(i,j):    #find nearest bonds
    up = (i - 1) % num_size  # up, j
    down = (i + 1) % num_size  # down, j
    right = (j - 1) % num_size  # i right
    left = (j + 1) % num_size  # i left
    cluster_ready = [[up, j], [down, j], [i, right], [i, left]]
    #print('clusterreadu',cluster_ready)
    return cluster_ready

def form_small_four_cluster(i,j,sigma_list, counted_list):
    cluster_ready= find_four_bonds(i,j)
    #print('slusready', cluster_ready, 'counted',counted_list)
    clus_with_no_dupli=sub_list(cluster_ready, counted_list)
    #print('clus_with_no_dupli',clus_with_no_dupli)

    useful_clus=[k for k in clus_with_no_dupli if choo_neigh_(sigma_list[k[0]][k[1]], sigma_list[i][j])==0]
    return useful_clus


def form_cluster(sigma_list,i0,j0):    #return cluster index list


   # print('ini', i0, j0)
    cluster_list=[[i0,j0]] #candidate list

    singletimecluster=[[i0,j0]]   # initialize first lattice
    while 1:    # find cluster
        singletime=[]

        for i in singletimecluster:     # singletimecluster: lattice last time join

            this_time_newlist = form_small_four_cluster(i[0], i[1], sigma_list, cluster_list) # new qualified lattice


            singletime = singletime+copy.copy(this_time_newlist)  # add 4 bonds

        singletimecluster = np.unique(np.array(copy.copy(singletime)) ,axis = 0).tolist()

        cluster_list =np.unique(np.array(cluster_list + singletimecluster) ,axis = 0).tolist()# add into cluster

        #print('sig',singletimecluster)
        if len(singletimecluster) == 0:
             break
    #print('len of cluster',len(cluster_list))
    return cluster_list



beta_list = [0.1+ i/40 for i in range(25)]
num_list=[4,6,8,9]
J = 1
for kk in num_list:
    num_size = kk
    sigma_list=de_sigma(num_size)
    ulist=[]
    mag_temp = []
    for beta in beta_list:
        steps = 2000
        print('beta', beta)
        mlist = []
        EE_list = []

        for j in range(steps):  # 开始集团更新
            # print(j)

            i0 = random.choice(range(num_size))
            j0 = random.choice(range(num_size))

            clusterlist = copy.copy(form_cluster(sigma_list, i0, j0))
            value = copy.copy(sigma_list[i0][j0])

            for i in clusterlist:  # 反转cluster中的每一个site
                # print('start',sigma_list[i[0]][i[1]])
                sigma_list[i[0]][i[1]] = -value
            # for i in sigma_list:
            # print('fliped',i)
            # print('end',sigma_list[i[0]][i[1]])
            E_tot = get_tot_energy(sigma_list)
            m = sum(np.array(sigma_list).flatten())
            EE_list.append(E_tot)
            mlist.append(abs(m))
        mlist=mlist[1000:]
        #print('end')
        magnetism_mean = np.mean(mlist)
        mag_temp.append(magnetism_mean)


       # plt.hist(mlist,bins=20)
        magcube=[(i/kk**2)**4 for i in mlist]
        magsqure=[(i/kk**2)**2 for i in mlist]
        u=1-(sum(magcube)/steps)/3*(sum(magsqure)/steps)**2
        ulist.append(u)

    plt.plot(beta_list, ulist,'*')
    plt.xlabel('beta')
    plt.ylabel('u')
plt.show()


















