#!/usr/bin/env python
# coding: utf-8

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.cm
import math
import networkx as nx


def get_peak(data_dir,contour_start,threshold):  
    dic, data = ng.bruker.read_pdata(data_dir, scale_data=False)
    Col_num=data.shape[1]
    Row_num=data.shape[0]
    # plot parameters
    cmap = matplotlib.cm.Blues_r    # contour map (colors to use for contours)
    contour_num     = 20        # number of contour levels
    contour_factor  = 1.20      # scaling factor between contour levels
    # calculate contour levels
    cl = contour_start * contour_factor ** np.arange(contour_num) 
    # detect all peaks with a threshold
    peaks = ng.peakpick.pick(data, pthres=threshold, algorithm="connected", msep=[1, 1])

    x = peaks["X_AXIS"]
    y =peaks["Y_AXIS"]
    a=plt.scatter(x, y, marker=".", s=11, color="r",zorder=1)
    b=plt.contour(data, cl, cmap=cmap, extent=(0, data.shape[1] - 1, 0, data.shape[0] - 1),alpha=0.5,zorder=0)
    # add markers for peak positions
    plt.xlim(0,Col_num)
    plt.ylim(0,Row_num)
    plt.show()
    return(x,y)

def get_half_data(x,y,data_dir):
    middle_x=[]
    middle_y=[]
    New_x=[]
    New_y=[]
    for i in range(len(x)):
        if y[i]-x[i]>=2:
            New_x.append(x[i])
            New_y.append(y[i])
        elif abs(y[i]-x[i])<=2:
            middle_x.append(x[i])
            middle_y.append(y[i])
        else:
            continue
    plt.scatter(New_x, New_y, marker=".", s=11, color="r",zorder=1)
    plt.scatter(middle_x, middle_y, marker=".", s=11, color="black",zorder=1)
    dic, data = ng.bruker.read_pdata(data_dir, scale_data=False)
    Col_num=data.shape[1]
    Row_num=data.shape[0]
    plt.xlim(0,Col_num)
    plt.ylim(0,Row_num)
    plt.show()
    return(New_x,New_y)


#单位转化
def unit_conversion(x_pts,y_pts,data_dir):
    dic, data = ng.bruker.read_pdata(data_dir, scale_data=False)
    udic = ng.bruker.guess_udic(dic, data)
    x_ppm_total=[]
    x_ppm=[]
    y_ppm_total=[]
    y_ppm=[]
    for a in x_pts:
        for i in a:
            delta_x=-udic[1]["sw"]/(udic[1]["size"]*udic[1]["obs"])
            first_x=udic[1]["car"]/udic[1]["obs"]-delta_x*udic[1]["size"]/2
            chemical_shift_x=(i * delta_x) + first_x
            x_ppm.append( round(chemical_shift_x, 2))
        x_ppm_total.append(x_ppm)
        x_ppm=[]
    for b in y_pts:
        for j in b:
            delta_y=-udic[0]["sw"]/(udic[0]["size"]*udic[0]["obs"])
            first_y=udic[0]["car"]/udic[0]["obs"]-delta_y*udic[0]["size"]/2
            chemical_shift_y=(j * delta_y) + first_y
            y_ppm.append(round(chemical_shift_y, 2))
        y_ppm_total.append(y_ppm)
        y_ppm=[]
    return (x_ppm_total,y_ppm_total)


def group_compare(group1,group2,threshold):
    hit=[]
    for a in group1:
        for b in group2:
            if abs(a-b)<=threshold:
                hit.append(a)
            else:
                continue
    if len(hit)==0:
        return("N")
    else:
        return("Y")

def group_find(x,y,threshold):
    #根据碳的距离来分组
    x=np.array(x)
    y=np.array(y)
    groups=[]
    group=[]
    C_list=y.tolist()
    #每次总是选出列表的第一个值
    while len(C_list)>0:
        a=C_list[0]
        group.append(a)
        C_list.remove(a)
        remove_index=[]
        #筛选出和前面选出的值距离相差1以内的值，组成一个组
        for i,b in enumerate(C_list):
            if abs(a-b)<=threshold:
                group.append(b)
                remove_index.append(i)
            else:
                continue
        #从列表中去掉已经成功分组的值
        indexes = remove_index
        for index in sorted(indexes, reverse=True):
            del C_list[index]
        #把一个组汇总在总列表中
        groups.append(group)
        #清空组，进行下一次循环
        group=[]

    #根据碳信号的分组来对氢信号进行分组
    proton_list=x.tolist()
    groups_H=[]
    count=0
    for i in range(len(groups)):
        group_H=proton_list[count:count+len(groups[i])]
        count=count+len(groups[i])
        groups_H.append(group_H)
        group_H=[]

    #根据氢信号来进一步合并分组：
    #如果组与组之间有相同的氢信号，说明这两组来自于同一个化合物
    result_row=[]
    result_table=[]
    group_New=groups_H
    group_New_C=groups
    result_table_C=[]
    index=[]
    for x,a in enumerate(group_New):
        if len(result_table)==0:
            result_table.append(a)
            result_table_C.append(group_New_C[x])
        else:
            for order,b in enumerate(result_table):
                result=group_compare(a,b,threshold)
                if result=="Y":
                    index.append(order)
                else:
                    continue
            #比较结束后，把a放在列表的最后
            result_table.append(a)
            result_table_C.append(group_New_C[x])
            #相同的组合并在a组中去
            for i in index:
                result_table[-1]=result_table[-1]+result_table[i]
                result_table_C[-1]=result_table_C[-1]+result_table_C[i]
            #删除被移动的组，要倒着删除索引值才对
            index.reverse()
            for j in index:
                del result_table[j]
                del result_table_C[j]
            index=[]
    return (result_table,result_table_C)

def connection_network_HC(x_ppm,y_ppm,x_ppm1,y_ppm1,cor_BC):
    #根据分组绘制网络图
    fontsize=20
    nodesize=1800
    fig = plt.figure(figsize=(10, 10))
    len_group=len(x_ppm)*3
    for i in range(len(x_ppm)):
        plt.subplot(math.ceil(len_group/3),3,3*i+1)
        G=nx.Graph()
        New_data=x_ppm[i]+y_ppm[i]
        #每一个组画一个网络图，删除重复的元素，首先确定网络图的节点
        lst = []
        for el in New_data:
            if lst.count(el) < 1:
                lst.append(el)
        point =lst
        G.add_nodes_from(point)
        #然后确定节点的连接关系
        edglist=[]
        n=len(x_ppm[i])
        for j in range(n):
            edglist.append((x_ppm[i][j],y_ppm[i][j]))
        G.add_edges_from(edglist)
        #接下来定义网络图，画图
        position = nx.circular_layout(G)
        nx.draw_networkx_nodes(G,position, nodelist=point, node_size= nodesize, node_color="yellow",edgecolors="yellow")
        nx.draw_networkx_edges(G,position,width=3.0)
        nx.draw_networkx_labels(G,position,font_size=fontsize)
        #设置边缘
        x_values, y_values = zip(*position.values())
        x_max = max(x_values)
        x_min = min(x_values)
        y_max = max(y_values)
        y_min = min(y_values)
        x_margin = (x_max - x_min) * 0.25
        y_margin = (y_max - y_min) * 0.25
        plt.xlim(x_min - x_margin, x_max + x_margin)
        plt.ylim(y_min - y_margin, y_max + y_margin)
        plt.title("Group"+str(i+1)+" "+"COSY")
    for i in range(len(x_ppm1)):
        plt.subplot(math.ceil(len_group/3),3,3*i+2)
        G=nx.Graph()
        New_data=x_ppm1[i]+y_ppm1[i]
        #每一个组画一个网络图，删除重复的元素，首先确定网络图的节点
        lst = []
        for el in New_data:
            if lst.count(el) < 1:
                lst.append(el)
        point =lst
        G.add_nodes_from(point)
        #然后确定节点的连接关系
        edglist=[]
        n=len(x_ppm1[i])
        for j in range(n):
            edglist.append((x_ppm1[i][j],y_ppm1[i][j]))
        G.add_edges_from(edglist)
        #接下来定义网络图，画图
        position = nx.circular_layout(G)
        nx.draw_networkx_nodes(G,position, nodelist=point, node_size= nodesize, node_color="red",edgecolors="yellow")
        nx.draw_networkx_edges(G,position,width=3.0)
        nx.draw_networkx_labels(G,position,font_size=fontsize)
        #设置边缘
        x_values, y_values = zip(*position.values())
        x_max = max(x_values)
        x_min = min(x_values)
        y_max = max(y_values)
        y_min = min(y_values)
        x_margin = (x_max - x_min) * 0.25
        y_margin = (y_max - y_min) * 0.25
        plt.xlim(x_min - x_margin, x_max + x_margin)
        plt.ylim(y_min - y_margin, y_max + y_margin)
        plt.title("Group"+str(i+1)+" "+"HSQC-COSY")
    for i in range(len(cor_BC)):
        plt.subplot(math.ceil(len_group/3),3,3*i+3)
        G=nx.Graph()
        New_data=cor_BC[i]
        #每一个组画一个网络图，删除重复的元素，首先确定网络图的节点
        lst = [i+1]
        for el in New_data:
            if lst.count(el) < 1:
                lst.append(el)
        point =lst
        G.add_nodes_from(point)
        #然后确定节点的连接关系
        edglist=[]
        n=len(cor_BC[i])
        for j in range(n):
            edglist.append((i+1,cor_BC[i][j]))
        G.add_edges_from(edglist)
        #接下来定义网络图，画图
        position = nx.circular_layout(G)
        nx.draw_networkx_nodes(G,position, nodelist=point,  node_size= nodesize, node_color="green")
        nx.draw_networkx_edges(G,position,style="dashed")
        nx.draw_networkx_labels(G,position,font_size=fontsize)
        #设置边缘
        x_values, y_values = zip(*position.values())
        x_max = max(x_values)
        x_min = min(x_values)
        y_max = max(y_values)
        y_min = min(y_values)
        x_margin = (x_max - x_min) * 0.25
        y_margin = (y_max - y_min) * 0.25
        plt.xlim(x_min - x_margin, x_max + x_margin)
        plt.ylim(y_min - y_margin, y_max + y_margin)
        plt.title("HMBC correlations of Group"+str(i+1))
    plt.tight_layout(pad=1)
    plt.savefig("netwrok.png")
    plt.show()


# #根据分组绘制网络图
# def connection_network_BC(cor_total):
#     fig = plt.figure(figsize=(16, 16))
#     len_group=len(cor_total)
#     for i in range(len_group):
#         plt.subplot(math.ceil(len_group/3),3,i+1)
#         G=nx.Graph()
#         New_data=cor_total[i]
#         #每一个组画一个网络图，删除重复的元素，首先确定网络图的节点
#         lst = [i+1]
#         for el in New_data:
#             if lst.count(el) < 1:
#                 lst.append(el)
#         point =lst
#         G.add_nodes_from(point)
#         #然后确定节点的连接关系
#         edglist=[]
#         n=len(cor_total[i])
#         for j in range(n):
#             edglist.append((i+1,cor_total[i][j]))
#         G.add_edges_from(edglist)
#         #接下来定义网络图，画图
#         position = nx.circular_layout(G)
#         nx.draw_networkx_nodes(G,position, nodelist=point, node_size=900, node_color="green")
#         nx.draw_networkx_edges(G,position)
#         nx.draw_networkx_labels(G,position)
#     plt.tight_layout(pad=1)
#     plt.show()


#结合QC图谱将相连的氢转化为相连的碳
def COSY_conversion1(ppm_COSY,x_ppm_QC,y_ppm_QC,threshold):
    C_cor=[]
    C_cor_total=[]
    count=0
    for i in range(len(ppm_COSY)):
        for j in range(len(ppm_COSY[i])):
            delta_list=[]
            for a in range(len(x_ppm_QC[0])):
                delta=abs(ppm_COSY[i][j]-x_ppm_QC[0][a])
                delta_list.append(delta)
            if min(delta_list)<=threshold:
                index=delta_list.index(min(delta_list))
                C_cor.append(y_ppm_QC[0][index])
                count=1
            if count==0:
                C_cor.append("NAN")
            count=0
        C_cor_total.append(C_cor)
        C_cor=[]
    return(C_cor_total)

#结合QC图谱将相连的氢转化为相连的碳
def COSY_conversion(ppm_COSY,x_ppm_QC,y_ppm_QC,threshold):
    C_cor=[]
    C_cor_total=[]
    count=0
    for i in range(len(ppm_COSY)):
        for j in range(len(ppm_COSY[i])):
            for a in range(len(x_ppm_QC[0])):
                if abs(ppm_COSY[i][j]-x_ppm_QC[0][a])<=threshold:
                    C_cor.append(y_ppm_QC[0][a])
                    count=1
                    break
                else:
                    continue
            if count==0:
                C_cor.append("NAN")
            count=0
        C_cor_total.append(C_cor)
        C_cor=[]
    return(C_cor_total)

def delete_NAN(Carbon_x,Carbon_y,x_ppm_COSY,y_ppm_COSY):
    index1=[]
    index2=[]
    #删除每一行中值为“NAN”的元素
    for i in range(len(Carbon_x)):
        for j in range(len(Carbon_x[i])):
            if Carbon_x[i][j]=="NAN":
                index1.append(i)
                index2.append(j)
    index1.reverse()
    index2.reverse()
    for a in range(len(index1)):
        del Carbon_x[index1[a]][index2[a]] 
        del Carbon_y[index1[a]][index2[a]]
        del x_ppm_COSY[index1[a]][index2[a]]
        del y_ppm_COSY[index1[a]][index2[a]]
    #删除空的行
    index3=[]
    for i in range(len(Carbon_x)):
        if len(Carbon_x[i])==0:
            index3.append(i)
    index3.reverse()
    for b in range(len(index3)):
        del Carbon_x[index3[b]]
        del Carbon_y[index3[b]]
        del x_ppm_COSY[index3[b]]
        del y_ppm_COSY[index3[b]]
    return (Carbon_x,Carbon_y,x_ppm_COSY,y_ppm_COSY)

def find_bc_cor (Proton_sys,Carbon_sys,x_ppm_BC,y_ppm_BC,H_threshold,C_threshold):
    cor_total_H2C=[]
    cor1=[]
    #寻找和本自旋系统的碳有相关的氢，并把其系统编号存在对应的数组中
    for i in range(len(Carbon_sys)):
        for j in range(len(Carbon_sys[i])):
            for a in range(len(y_ppm_BC[0])):
                #如果BC相关表中有，则看对应的氢，判断这个氢是否是本自旋系统的
                if abs(Carbon_sys[i][j]-y_ppm_BC[0][a])<=C_threshold:
                    for m in range(len(Proton_sys)):
                        for n in range(len(Proton_sys[m])):
                            if abs(x_ppm_BC[0][a]-Proton_sys[m][n])<=H_threshold:
                                cor1.append(m+1)
        cor_total_H2C.append(cor1)
        cor1=[]

    cor_total_C2H=[]
    cor2=[]
    #寻找和本自旋系统的氢有相关的碳，并把其系统编号存在对应的数组中
    for i in range(len(Proton_sys)):
        for j in range(len(Proton_sys[i])):
            for a in range(len(x_ppm_BC[0])):
                #如果BC相关表中有，则看对应的碳，判断这个碳是否是本自旋系统的
                if abs(Proton_sys[i][j]-x_ppm_BC[0][a])<=H_threshold:
                    for m in range(len(Carbon_sys)):
                        for n in range(len(Carbon_sys[m])):
                            if abs(y_ppm_BC[0][a]-Carbon_sys[m][n])<=C_threshold:
                                cor2.append(m+1)
        cor_total_C2H.append(cor2)
        cor2=[]

    #合并两部分的相关
    cor_total=[]
    cor_total_1=[]
    for i in range(len(cor_total_H2C)):
        cor_total_1=cor_total_H2C[i]+cor_total_C2H[i]
        cor_total.append(cor_total_1)
        cor_total_1=[]
    return(cor_total)

def get_csv_data(CarbonShiftDir,QC_H,QC_C,COSY_H1,COSY_H2,COSY_Cx,COSY_Cy,BC_H,BC_C,threshold_H=0.05,threshold_C=0.5):
    BC_H2C=[]
    BC_C2C=[]
    #将BC中碳氢相关转化为碳碳相关（这个地方匹配就可能会出现错误，导致后面翻译出来的BC相关不准确，改为最小距离匹配可能好些）
    #在这一步中QC中没有找到H的BC杂质都被去除了
    for i in range(len(BC_H)):
        DeltaList=[]
        for j in range(len(QC_H)):
            Delta=abs(BC_H[i]-QC_H[j])
            DeltaList.append(Delta)
        if min(DeltaList)<=threshold_H:
            #可能存在QC_C并不在目标化合物化学位移值里面，所以转化的数据也需要进一步清洗
            BC_H2C.append(QC_C[DeltaList.index(min(DeltaList))])
            BC_C2C.append(BC_C[i])
    
    #接下来进行BC中的碳数据的清洗，进一步去除提取到的杂质信号
    #BC_H2C是根据QC转化的，而QC的碳来自13C NMR，所以BC_H2C的碳不许要清洗
    #接下来是清洗BC-C2C的信号
    a=pd.read_csv(CarbonShiftDir)
    CarbonShift=a.iloc[:,0].to_list()
    BC_H2C_1=[]
    BC_C2C_1=[]
    for i in range(len(BC_C2C)):
        DeltaList=[]
        for j in range(len(CarbonShift)):
            Delta=abs(BC_C2C[i]-CarbonShift[j])
            DeltaList.append(Delta)
        if min(DeltaList)<=threshold_C:
            BC_H2C_1.append(BC_H2C[i])
            BC_C2C_1.append(CarbonShift[DeltaList.index(min(DeltaList))])
    
    #将BC-H2C按照降序进行排序
    a=[BC_H2C_1,BC_C2C_1]
    b=pd.DataFrame(a)
    df=pd.DataFrame(b.values.T,columns=["BC_H2C","BC_C2C"])
    df = df.sort_values(by=['BC_H2C'], ascending=[False])
    BC_H2C_1=df.iloc[:,0].to_list()
    BC_C2C_1=df.iloc[:,1].to_list()
    
    #将根据COSY相关得到的二维列表转化为1维列表
    COSY_H1_1=[]
    COSY_H2_1=[]
    for i in range(len(COSY_H1)):
        COSY_H1_1=COSY_H1_1+COSY_H1[i]
        COSY_H2_1= COSY_H2_1+COSY_H2[i]
    
    #将根据COSY相关得到的二维列表转化为1维列表
    COSY_C1=[]
    COSY_C2=[]
    for i in range(len(COSY_Cx)):
        COSY_C1=COSY_C1+COSY_Cx[i]
        COSY_C2=COSY_C2+COSY_Cy[i]
        
        
    data=[]
    data.append(QC_H)
    data.append(QC_C)
    data.append(COSY_H1_1)
    data.append(COSY_H2_1)
    data.append(BC_H)
    data.append(BC_C)
    data.append(COSY_C1)
    data.append(COSY_C2)
    data.append(BC_H2C_1)
    data.append(BC_C2C_1)
    b=pd.DataFrame(data)
    index_colums = ['QC_H','QC_C','COSY_H1','COSY_H2','BC_H','BC_C','COSY_C1',
                    'COSY_C2','BC_H2C','BC_C2C']
    b=pd.DataFrame(b.values.T,columns=index_colums)
    return(b)

def figure_show(result_table,result_table_C):
    #重新读取1D NMR数据
    data_dir =dir_H
    dicH, dataH = ng.bruker.read_pdata(data_dir, scale_data=False)

    udicH = ng.bruker.guess_udic(dicH, dataH)
    # 实现格式转换成NMRPipe文件，这个新文件的参数文件夹没有改变
    H = ng.convert.converter()
    H.from_bruker(dicH, dataH, udicH)
    ng.pipe.write("1H_pipe.fid", *H.to_pipe(), overwrite=True)

    #重新读取新文件,获取单位转化对象
    dicH, dataH = ng.pipe.read("1H_pipe.fid")
    ucH = ng.pipe.make_uc(dicH, dataH)
    dic, data = ng.bruker.read_pdata(dir_2D, scale_data=False)
    Col_num=data.shape[1]
    Row_num=data.shape[0]
    F1P_1=dic.get('procs').get('F1P')
    F2P_1=dic.get('procs').get('F2P')
    F1P_2=dic.get('proc2s').get('F1P')
    F2P_2=dic.get('proc2s').get('F2P')
    #绘图
    #绘制氢谱
    fig=plt.figure(figsize=(12,8))
    plt.subplot2grid((5,5),(0,1),colspan=5)
    front_size=16
    plt.plot(ucH.ppm_scale(), dataH)
    plt.xlim(10,0)
    plt.ylim(0,2*1e8)
    plt.yticks([])
    ax=plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    #绘制二维相关谱
    plt.subplot2grid((5,5),(1,1),colspan=5,rowspan=5)
    color=["red","blue","green","brown","cyan","olive","purple","pink","orange","statebule","lightgreen"]
    for i in range(len(result_table)):
        for j in range(len(result_table[i])):
            x=(1024-result_table[i][j])/1024*( 9.59+0.42)-0.29
            y=(1024-result_table_C[i][j])/1024*( 9.59+0.42)-0.29
            plt.scatter(x,y)  
            plt.xlim(10,0)
            plt.ylim(10,0)
#     plt.yticks([])
#     plt.xticks([])
    plt.rcParams['xtick.direction']="in"
    plt.rcParams['ytick.direction']="in"

    #绘制纵向氢谱      
    plt.subplot2grid((5,5),(1,0),rowspan=5)
    plt.plot(dataH,ucH.ppm_scale())
    plt.ylim(10,0)
    plt.xlim(2*1e8,0)
    plt.xticks([])
    ax=plt.gca()
    ax.yaxis.set_ticks_position('right')
    ax=plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    plt.tight_layout(pad=0)

#如果QC峰中的碳化学位移没有在1D NMR中出现，则删除，这样有利提高后面匹配的准确性
#还是用最近邻匹配加阈值
#把QC中的碳替换成了1D NMR中的碳，有利于后面校对数据
def QCMove (CarbonShiftDir,x_ppm_QC,y_ppm_QC,therohold):
    """Remove the noesy signal in HSQC peaks"""
    a=pd.read_csv(CarbonShiftDir)
    CarbonShift=a.iloc[:,0].to_list()
    x_ppm_QC_1=[]
    y_ppm_QC_1=[]
    for i in range(len(y_ppm_QC)):
        DeltaList=[]
        for j in range(len(CarbonShift)):
            Delta=abs(y_ppm_QC[i]-CarbonShift[j])
            DeltaList.append(Delta)
        if min(DeltaList)<=therohold:
            x_ppm_QC_1.append(x_ppm_QC[i])
            y_ppm_QC_1.append(CarbonShift[DeltaList.index(min(DeltaList))])
    return (x_ppm_QC_1, y_ppm_QC_1)
def threshold_input(s,b):
    if s=="" and b=="HSQC":
        contour_start=5000000
        threshold=5000000
    elif s=="" and b=="COSY":
        contour_start=5000000
        threshold=5000000
    elif s=="" and b=="HMBC":
        contour_start=5000000
        threshold=5000000
    else:
        s = s.split(' ')
        a=list(filter(lambda x:x!='',s))
        a=[float(i) for i in a]
        a=[int(i) for i in a]
        contour_start=a[0]
        threshold=a[1]
    print("contour_start:"+str(contour_start)+"threshold:"+str(threshold))
    return(contour_start,threshold) 
#误差可能来源于两方面，第一方面是峰提取的时候带来的误差。第二是QC信号转化BC和QC的时候可能出现的误差
#第二类误差可能是QC中两个碳信号隔得太近，导致COSY和BC在将H转化为C时，发生了错配，这时候COSY和BC可能是距离很近的碳之间的错误转换
#第二类误差也可能是QC中两个氢信号隔得太近，导致COSY和BC在将H转化为C时，发生了错配，这时候COSY和BC可能是距离很远的碳之间的错误转换
def poss_error(thresohold_H,thresohold_C):
    cor_data=pd.read_csv("cor_data.csv")
    QC_H=cor_data["QC_H"].tolist()
    QC_C=cor_data["QC_C"].tolist()
    QC_H= [a_ for a_ in QC_H if a_ == a_]
    QC_C=[a_ for a_ in QC_C if a_ == a_]
    poss_error1=[]
    for i in range(len(QC_H)):
        for j in range(len(QC_H)):
            if 0<abs(QC_H[i]-QC_H[j])<=thresohold_H:
                a=(QC_C[i],QC_C[j],"H")
                poss_error1.append(a)
    poss_error2=[]
    QC_C=list(set(QC_C))
    for i in range(len(QC_C)):
        for j in range(len(QC_C)):
            if 0<abs(QC_C[i]-QC_C[j])<=thresohold_C:
                b=(QC_C[i],QC_C[j],"C")
                poss_error2.append(b)
    return(poss_error1,poss_error2)
    
def NMR_Data_2D_Process(NMR_Data_2D):
    peakdata=pd.read_csv(NMR_Data_2D)
    x_ppm_QC=peakdata["HSQC_H"].tolist()
    y_ppm_QC=peakdata["HSQC_C"].tolist()
    x_ppm_COSY1=peakdata["COSY_H1"].tolist()
    y_ppm_COSY1=peakdata["COSY_H2"].tolist()
    x_ppm_BC=peakdata["HMBC_H"].tolist()
    y_ppm_BC=peakdata["HMBC_C"].tolist()
    list1=[x_ppm_QC,y_ppm_QC,x_ppm_COSY1,y_ppm_COSY1, x_ppm_BC,y_ppm_BC]
    list2=[]
    for a in list1:
        list3=a[0].split(",")
        list3=[float(i) for i in list3]
        list2.append(list3)
    return (list2)


class analysis_2DNMR (object): 
    def __init__(self,CarbonListFile):
        self.CarbonListFile=CarbonListFile
    def autoanalysis(self,dir_QC,dir_COSY,dir_BC):
        self.dir_QC=dir_QC
        self.dir_COSY=dir_COSY
        self.dir_BC=dir_BC
        #QC峰提取，并绘图，根据结果调整阈值，设置峰提取的阈值,contour_start设置得越低，原图中点越多
        a=0
        while a==0:
            QC_parameter=input("Please enter the value of contour_start(Default 5000000) and threshold (Default 5000000) for HSQC")
            contour_start,threshold=threshold_input(QC_parameter,b="HSQC")
            x_pts_QC,y_pts_QC=get_peak(self.dir_QC,contour_start=contour_start,threshold =threshold)
            #QC单位转换：pts转化为ppm
            x_ppm_QC,y_ppm_QC=unit_conversion([x_pts_QC],[y_pts_QC], dir_QC)
            #去除掉多余的QC相关信息
            x_ppm_QC,y_ppm_QC=QCMove(self.CarbonListFile,x_ppm_QC[0],y_ppm_QC[0],0.5)
            #将QC相关信息储存在csv文件中
            HSQC_Peaks=[]
            HSQC_Peaks.append(x_ppm_QC)
            HSQC_Peaks.append(y_ppm_QC)
            HSQC_Peaks=pd.DataFrame(HSQC_Peaks)
            HSQC_Peaks=pd.DataFrame(HSQC_Peaks.values.T,columns=["QC_H","QC_C"])
            HSQC_Peaks.to_csv("HSQC_Peaks.csv",index=False)
            #因为QC实在太重要了，可以在这一步检查一下QC是不是标正确了，如果不正确可以人为核对一下，有无缺失，有无错配
            task=input(" enter “C”to continue or enter any key to Re-extract HSQC")
            if task=='C':
                a=1
            else:
                a=0
        #重新读取核对后的QC信号
        df=pd.read_csv("HSQC_Peaks.csv")
        x_ppm_QC=df.iloc[:,0].tolist()
        y_ppm_QC=df.iloc[:,1].tolist()
        x_ppm_QC=[x_ppm_QC]
        y_ppm_QC=[y_ppm_QC]
        
        e=0
        while e==0:
            COSY_parameter=input("Please enter the value of contour_start(Default 5000000) and threshold (Default 5000000) for 1H 1HCOSY")
            contour_start,threshold=threshold_input(COSY_parameter,b="COSY")
            #COSY峰提取
            x_pts_COSY,y_pts_COSY=get_peak(self.dir_COSY,contour_start =contour_start,threshold =threshold)
            #删除掉一半的COSY峰,并且重新绘图
            x_pts1,y_pts1=get_half_data(x_pts_COSY,y_pts_COSY,dir_COSY)
            task=input(" enter “C”to continue or enter any key to Re-extract 1H 1H COSY")
            if task=='C':
                e=1
            else:
                e=0
        #对COSY相关信号进行分组,
        result_table_H1, result_table_H2=group_find(x_pts1,y_pts1,threshold =1)
        #COSY单位转换
        x_ppm_COSY,y_ppm_COSY=unit_conversion(result_table_H1,result_table_H2, dir_COSY)
        #绘制氢氢相关网络连接图，但是这里的x_ppm_COSY,y_ppm_COSY可能有些是噪音或者卫星峰的信号，最好是删除噪音以后再画图
        #connection_network(x_ppm_COSY,y_ppm_COSY)
        
        #结合QC图谱将相连的氢转化为相连的碳
        Carbon_x=COSY_conversion(x_ppm_COSY,x_ppm_QC,y_ppm_QC,threshold=0.02)
        Carbon_y=COSY_conversion(y_ppm_COSY,x_ppm_QC,y_ppm_QC,threshold=0.02)
        #删除匹配不上的位置，因为有横竖两组信号，所以要删两次
        Carbon_x1,Carbon_y1,x_ppm_COSY1,y_ppm_COSY1=delete_NAN(Carbon_x,Carbon_y,x_ppm_COSY,y_ppm_COSY)
        Carbon_y2,Carbon_x2,x_ppm_COSY2,y_ppm_COSY2=delete_NAN(Carbon_y1,Carbon_x1,x_ppm_COSY1,y_ppm_COSY1)
        #HMBC分析：判断自旋系统的邻接关系
        #首先合并同一个自旋系统里的碳和氢
        Carbon_sys=[]
        Proton_sys=[]
        for i in range(len(Carbon_x2)):
            Carbon_sys.append(Carbon_x2[i]+Carbon_y2[i])
            Proton_sys.append(x_ppm_COSY2[i]+y_ppm_COSY2[i])
            
        #提取HMBC的相关信号,并且将单位转化为ppm
        f=0
        while f==0:
            BC_parameter=input("Please enter the value of contour_start(Default 5000000) and threshold (Default 5000000) for HMBC")
            contour_start,threshold=threshold_input(BC_parameter,b="HMBC")
            x_pts_BC,y_pts_BC=get_peak(self.dir_BC,contour_start =contour_start,threshold =threshold)
            task=input(" enter “C”to continue or enter any key to Re-extract HMBC")
            if task=='C':
                f=1
            else:
                f=0
        x_ppm_BC,y_ppm_BC=unit_conversion([x_pts_BC], [y_pts_BC],dir_BC)
        cor_data=get_csv_data(self.CarbonListFile,x_ppm_QC[0],y_ppm_QC[0],x_ppm_COSY,y_ppm_COSY,
                              Carbon_x2,Carbon_y2,x_ppm_BC[0],y_ppm_BC[0],threshold_H=0.02)
        #分析自旋系统之间的BC相关关系，若自旋系统和其它自旋系统有相关，则在列表中记录下其他自旋系统的编号
        cor_total=find_bc_cor (Proton_sys,Carbon_sys,x_ppm_BC,y_ppm_BC,H_threshold=0.02,C_threshold=0.2)
        #可视化所有分析结果
        connection_network_HC(x_ppm_COSY2,y_ppm_COSY2,Carbon_x2,Carbon_y2,cor_total)
        self.cor_data=cor_data
        cor_data.to_csv("cor_data.csv")
        self.poss_error1,self.poss_error2=poss_error(0.03,0.5)
    def autoanalysisGUI(self,dir_QC,dir_COSY,dir_BC,contour_start_QC,contour_start_COSY,contour_start_BC,threshold_QC,threshold_COSY,threshold_BC):
		    self.dir_QC=dir_QC
		    self.dir_COSY=dir_COSY
		    self.dir_BC=dir_BC
		    #QC峰提取，并绘图，根据结果调整阈值，设置峰提取的阈值,contour_start设置得越低，原图中点越多
		    x_pts_QC,y_pts_QC=get_peak(self.dir_QC,contour_start_QC,threshold_QC)
		    #QC单位转换：pts转化为ppm
		    x_ppm_QC,y_ppm_QC=unit_conversion([x_pts_QC],[y_pts_QC], dir_QC)
		    #去除掉多余的QC相关信息
		    x_ppm_QC,y_ppm_QC=QCMove(self.CarbonListFile,x_ppm_QC[0],y_ppm_QC[0],0.5)
		    x_ppm_QC=[x_ppm_QC]
		    y_ppm_QC=[y_ppm_QC]
		    #COSY峰提取
		    x_pts_COSY,y_pts_COSY=get_peak(self.dir_COSY,contour_start =contour_start_COSY,threshold =threshold_COSY)
		    #删除掉一半的COSY峰,并且重新绘图
		    x_pts1,y_pts1=get_half_data(x_pts_COSY,y_pts_COSY,dir_COSY)
		
		    #对COSY相关信号进行分组,
		    result_table_H1, result_table_H2=group_find(x_pts1,y_pts1,threshold =1)
		    #COSY单位转换
		    x_ppm_COSY,y_ppm_COSY=unit_conversion(result_table_H1,result_table_H2, dir_COSY)
		    #绘制氢氢相关网络连接图，但是这里的x_ppm_COSY,y_ppm_COSY可能有些是噪音或者卫星峰的信号，最好是删除噪音以后再画图
		    #connection_network(x_ppm_COSY,y_ppm_COSY)
		
		    #结合QC图谱将相连的氢转化为相连的碳
		    Carbon_x=COSY_conversion(x_ppm_COSY,x_ppm_QC,y_ppm_QC,threshold=0.02)
		    Carbon_y=COSY_conversion(y_ppm_COSY,x_ppm_QC,y_ppm_QC,threshold=0.02)
		    #删除匹配不上的位置，因为有横竖两组信号，所以要删两次
		    Carbon_x1,Carbon_y1,x_ppm_COSY1,y_ppm_COSY1=delete_NAN(Carbon_x,Carbon_y,x_ppm_COSY,y_ppm_COSY)
		    Carbon_y2,Carbon_x2,x_ppm_COSY2,y_ppm_COSY2=delete_NAN(Carbon_y1,Carbon_x1,x_ppm_COSY1,y_ppm_COSY1)
		    #HMBC分析：判断自旋系统的邻接关系
		    #首先合并同一个自旋系统里的碳和氢
		    Carbon_sys=[]
		    Proton_sys=[]
		    for i in range(len(Carbon_x2)):
		        Carbon_sys.append(Carbon_x2[i]+Carbon_y2[i])
		        Proton_sys.append(x_ppm_COSY2[i]+y_ppm_COSY2[i])
		    #提取HMBC的相关信号,并且将单位转化为ppm
		    x_pts_BC,y_pts_BC=get_peak(self.dir_BC,contour_start =contour_start_BC,threshold =threshold_BC)
		    x_ppm_BC,y_ppm_BC=unit_conversion([x_pts_BC], [y_pts_BC],dir_BC)
		    cor_data=get_csv_data(self.CarbonListFile,x_ppm_QC[0],y_ppm_QC[0],x_ppm_COSY,y_ppm_COSY,
		                          Carbon_x2,Carbon_y2,x_ppm_BC[0],y_ppm_BC[0],threshold_H=0.02)
		    #分析自旋系统之间的BC相关关系，若自旋系统和其它自旋系统有相关，则在列表中记录下其他自旋系统的编号
		    cor_total=find_bc_cor (Proton_sys,Carbon_sys,x_ppm_BC,y_ppm_BC,H_threshold=0.02,C_threshold=0.2)
		    #可视化所有分析结果
		    connection_network_HC(x_ppm_COSY2,y_ppm_COSY2,Carbon_x2,Carbon_y2,cor_total)
		    self.cor_data=cor_data
		    cor_data.to_csv("cor_data.csv")
		    self.poss_error1,self.poss_error2=poss_error(0.03,0.5)
    def semi_autoanalysis(self,peakdatadir,threshold_H=0.015,threshold_C=0.2):
        peakdata=NMR_Data_2D_Process(peakdatadir)
        x_ppm_QC=peakdata[0]
        y_ppm_QC=peakdata[1]
        x_ppm_COSY1=peakdata[2]
        y_ppm_COSY1=peakdata[3]
        x_ppm_BC=peakdata[4]
        y_ppm_BC=peakdata[5]
        #去除掉多余的QC相关信息
        x_ppm_QC,y_ppm_QC=QCMove(self.CarbonListFile,x_ppm_QC,y_ppm_QC,0.5)
        #对COSY相关信号进行分组,
        x_ppm_COSY,y_ppm_COSY=group_find(x_ppm_COSY1,y_ppm_COSY1,threshold =threshold_H)
        print(x_ppm_COSY)
        print(y_ppm_COSY)
        #结合QC图谱将相连的氢转化为相连的碳
        Carbon_x=COSY_conversion1(x_ppm_COSY,[x_ppm_QC],[y_ppm_QC],threshold=threshold_H)
        Carbon_y=COSY_conversion1(y_ppm_COSY,[x_ppm_QC],[y_ppm_QC],threshold=threshold_H)
        #删除匹配不上的位置，因为有横竖两组信号，所以要删两次
        Carbon_x1,Carbon_y1,x_ppm_COSY,y_ppm_COSY=delete_NAN(Carbon_x,Carbon_y,x_ppm_COSY,y_ppm_COSY)
        Carbon_y2,Carbon_x2,x_ppm_COSY2,y_ppm_COSY2=delete_NAN(Carbon_y1,Carbon_x1,x_ppm_COSY,y_ppm_COSY)
        #将所有相关信息储存
        cor_data=get_csv_data(self.CarbonListFile,x_ppm_QC,y_ppm_QC,[x_ppm_COSY1],[y_ppm_COSY1],
                              Carbon_x2,Carbon_y2,x_ppm_BC,y_ppm_BC,threshold_H=0.02,threshold_C=0.2)
        #HMBC分析：判断自旋系统的邻接关系
        #首先合并同一个自旋系统里的碳和氢
        
        Carbon_sys=[]
        Proton_sys=[]
        for i in range(len(Carbon_x2)):
            Carbon_sys.append(Carbon_x2[i]+Carbon_y2[i])
            Proton_sys.append(x_ppm_COSY2[i]+y_ppm_COSY2[i])
        #分析自旋系统之间的BC相关关系，若自旋系统和其它自旋系统有相关，则在列表中记录下其他自旋系统的编号
        cor_total=find_bc_cor (Proton_sys,Carbon_sys,[x_ppm_BC],[y_ppm_BC],H_threshold=threshold_H,C_threshold=threshold_C)
        #可视化所有分析结果
        connection_network_HC(x_ppm_COSY,y_ppm_COSY,Carbon_x2,Carbon_y2,cor_total)
        self.cor_data=cor_data
        cor_data.to_csv("cor_data.csv")
        self.poss_error1,self.poss_error2=poss_error(0.03,0.5)
