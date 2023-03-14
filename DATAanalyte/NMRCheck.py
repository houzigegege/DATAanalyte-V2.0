#!/usr/bin/env python
# coding: utf-8

from rdkit import Chem
from rdkit.Chem import GetAdjacencyMatrix
from rdkit.Chem import  Draw
from rdkit.Chem.Draw import IPythonConsole 
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions 
import pandas as pd
import numpy as np
import os
import PySimpleGUI as sg
from itertools import permutations
import matplotlib.pyplot as plt
import heapq
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
from io import BytesIO
from PIL import Image
from cairosvg import svg2png

#将化合物用图表示，并且获得化合物的编号信息
#comp_dir：化合物mol文件
#Carbon_label：化合物化学位移所匹配的序号
#结构检查的思路：找到BC相关信号的碳化学位移，找到碳的编号，
#根据编号查找距离矩阵看看这两个碳的距离。如果距离大于2则说明：结构有问题或者归属有问题
class NMRcheck(object):
    def __init__(self,comp_dir,CarbonListFile,CorrelationFile):
        self.comp_dir=comp_dir
        self.CarbonListFile=CarbonListFile
        self.CorrelationFile=CorrelationFile
    def get_comp_figure(self):
        m= Chem.MolFromMolFile(self.comp_dir)
        opts =  DrawingOptions()
        opts.atomLabelFontSize=200
        opts.includeAtomNumbers=True
        opts.bondLineWidth=2
        opts.dblBondLengthFrac=0.8
        opts.dotsPerAngstrom=1000
        self.fig=Draw.MolToImage(m,options=opts,size=(200,200))
    def InputFileCreate(self):
        df=pd.read_csv(self.CarbonListFile)
        df1=df.iloc[:,0]
        df1.to_csv("InputFile.csv",index=False)
    def C_check(self):
        """检查不同类型的碳原子数量是否正确"""
        mol= Chem.MolFromMolFile(self.comp_dir)
        CH3_list=[]
        CH2_list=[]
        CH_list=[]
        C_list=[]
        Atoms=mol.GetAtoms()
        for i,atom in enumerate(Atoms):
            #首先获得原子元素符号，判断其是否为碳
            a=atom.GetSymbol()
            #如果是碳进一步判断碳的类型
            if a=="C":
                if atom.GetTotalNumHs()==3:
                    CH3_list.append(i)
                elif atom.GetTotalNumHs()==2:
                    CH2_list.append(i)
                elif atom.GetTotalNumHs()==1:
                    CH_list.append(i)
                elif atom.GetTotalNumHs()==0:
                    C_list.append(i)
        #获得化合物的碳信息
        df_shift=pd.read_csv(self.CarbonListFile)
        Type_list=df_shift.iloc[:,1].tolist()
        CH3_real=Type_list.count("q")
        CH2_real=Type_list.count("t")
        CH_real=Type_list.count("d")
        C_real=Type_list.count("s")
        #比较实验数据和推测的结构中的数据
        C_check_result=[]
        C_check_result.append([CH3_real,len(CH3_list),CH3_real==len(CH3_list)])
        C_check_result.append([CH2_real,len(CH2_list),CH2_real==len(CH2_list)])
        C_check_result.append([CH_real,len(CH_list),CH_real==len(CH_list)])
        C_check_result.append([C_real,len(C_list),C_real==len(C_list)])
        self.C_check_result=C_check_result
    def BC_check(self,Carbon_label_BC,threshold):
        """检查结构和HMBC信号是否符合
          Carbon_label_BC:与化学位移列表相对应的碳的编号，这个编号是通过RDkit打开结构以后可以看到
          threshold:碳化学位移匹配的阈值
          """
        self.Carbon_label=Carbon_label_BC
        self.threshold=threshold
        #获得化合物的距离矩阵
        comp= Chem.MolFromMolFile(self.comp_dir)
        Distance = np.array(Chem.GetDistanceMatrix(comp), "d")
        DistanceMatrix =pd.DataFrame(Distance)
        #获得化合物的化学位移信息
        df_shift=pd.read_csv(self.CarbonListFile)
        chemical_shift=df_shift.iloc[:,0].tolist()
        #获得化合物的BC相关信息
        df1=pd.read_csv(self.CorrelationFile)
        BC_H2C=df1["BC_H2C"].tolist()
        BC_C2C=df1["BC_C2C"].tolist()
        BC_H2C= [a_ for a_ in BC_H2C if a_ == a_]
        BC_C2C= [a_ for a_ in BC_C2C if a_ == a_]
        #首先找到每一组BC和与化学位移的差值，看看其最小值是不是都小于阈值
        Normal=[]
        Abnormal=[]
        for i in range(len(BC_H2C)):
            delta_H2C=[]
            delta_C2C=[]
            #首先计算每组和化学位移的差值
            for j in range(len(chemical_shift)):
                a=abs(BC_H2C[i]-chemical_shift[j])
                b=abs(BC_C2C[i]-chemical_shift[j])
                delta_H2C.append(a)
                delta_C2C.append(b)
            #print(i,min(delta_H2C),min(delta_C2C))
            #获取最小值的索引所对应的碳原子编号
            min_H2C=Carbon_label_BC[delta_H2C.index(min(delta_H2C))]
            min_C2C=Carbon_label_BC[delta_C2C.index(min(delta_C2C))]
            #获取最小值的对应的化学位移
            shift_H2C=chemical_shift[delta_H2C.index(min(delta_H2C))]
            shift_C2C=chemical_shift[delta_C2C.index(min(delta_C2C))]
            #print(i,min_H2C, min_C2C)
            #判断最小差值小于阈值，则去距离矩阵中寻找距离值
            if min(delta_H2C)<=threshold and min(delta_C2C)<=self.threshold:
                C_distance1=DistanceMatrix.iloc[min_H2C,min_C2C]
                #print(i,C_distance1)
                if C_distance1<=2:
                    Normal.append([min_H2C,min_C2C,C_distance1,shift_H2C,shift_C2C])
                else:
                    Abnormal.append([min_H2C,min_C2C,C_distance1,shift_H2C,shift_C2C])
        
        self.Normal_BC=Normal
        self.Abnormal_BC=Abnormal
        Abnormal=pd.DataFrame(Abnormal)
        Normal=pd.DataFrame(Normal)
        if len(Abnormal)>0:
            #获取异常的碳原子序号
            Abnormal_label1=Abnormal.iloc[:,0].tolist()
            Abnormal_label2=Abnormal.iloc[:,1].tolist()
            Abnormal_label=Abnormal_label1+Abnormal_label2
            ls=set(Abnormal_label)
            Abnormal_label=list(ls)
            self.Abnormal_BClabel=Abnormal_label
        else:
            self.Abnormal_BClabel=[]
    def COSY_check (self,Carbon_label,threshold):
        self.Carbon_label=Carbon_label
        self.threshold=threshold
        comp= Chem.MolFromMolFile(self.comp_dir)
        #获得化合物的距离矩阵
        Distance = np.array(Chem.GetDistanceMatrix(comp), "d")
        DistanceMatrix =pd.DataFrame(Distance)
        #获得化合物的化学位移信息
        df_shift=pd.read_csv(self.CarbonListFile)
        chemical_shift=df_shift.iloc[:,0].tolist()
        #获得化合物的COSY相关信息
        df1=pd.read_csv(self.CorrelationFile)
        COSY_C1=df1["COSY_C1"].tolist()
        COSY_C2=df1["COSY_C2"].tolist()
        COSY_C1= [a_ for a_ in COSY_C1 if a_ == a_]
        COSY_C2= [a_ for a_ in COSY_C2 if a_ == a_]
        #首先找到每一组BC和与化学位移的差值，看看其最小值是不是都小于阈值
        Normal=[]
        Abnormal=[]
        for i in range(len(COSY_C1)):
            delta_COSY1=[]
            delta_COSY2=[]
            #首先计算每组和化学位移的差值
            for j in range(len(chemical_shift)):
                a=abs(COSY_C1[i]-chemical_shift[j])
                b=abs(COSY_C2[i]-chemical_shift[j])
                delta_COSY1.append(a)
                delta_COSY2.append(b)
            #print(i,min(delta_H2C),min(delta_C2C))
            #获取最小值的索引所对应的碳原子编号
            min_COSY1=Carbon_label[delta_COSY1.index(min(delta_COSY1))]
            min_COSY2=Carbon_label[delta_COSY2.index(min(delta_COSY2))]
            #获取最小值的对应的化学位移
            shift_COSY1=chemical_shift[delta_COSY1.index(min(delta_COSY1))]
            shift_COSY2=chemical_shift[delta_COSY2.index(min(delta_COSY2))]
            #print(i,min_H2C, min_C2C)
            #判断最小差值小于阈值，则去距离矩阵中寻找距离值
            if min(delta_COSY1)<=threshold and min(delta_COSY2)<=self.threshold:
                C_distance1=DistanceMatrix.iloc[min_COSY1,min_COSY2]
                #print(i,C_distance1)
                if C_distance1<=1:
                    Normal.append([min_COSY1,min_COSY2,C_distance1,shift_COSY1,shift_COSY2])
                else:
                    Abnormal.append([min_COSY1,min_COSY2,C_distance1,shift_COSY1,shift_COSY2])
        self.Normal_COSY=Normal
        self.Abnormal_COSY=Abnormal
        #print(Normal,Abnormal)
        Abnormal=pd.DataFrame(Abnormal)
        Normal=pd.DataFrame(Normal)
        #获取异常的碳原子序号
        if len(Abnormal)>0:
            Abnormal_label1=Abnormal.iloc[:,0].tolist()
            Abnormal_label2=Abnormal.iloc[:,1].tolist()
            Abnormal_label=Abnormal_label1+Abnormal_label2
            ls=set(Abnormal_label)
            Abnormal_label=list(ls)
            self.Abnormal_COSYlabel=Abnormal_label
        else: 
            self.Abnormal_COSYlabel=[]
    def CorrectP(self):
#         df=pd.DataFrame(self.C_check_result)
#         a=df.iloc[:,3].tolist()
#         if "False" in a:
#             CorrectP=0
#         else:
        self.CorrectP_COSY=len(self.Normal_COSY)/(len(self.Normal_COSY)+len(self.Abnormal_COSY))
        self.CorrectP_BC=len(self.Normal_BC)/(len(self.Normal_BC)+len(self.Abnormal_BC))
        self.CorrectP_Total=(len(self.Normal_COSY)+len(self.Normal_BC))/(len(self.Normal_COSY)+len(self.Abnormal_COSY)+len(self.Normal_BC)+len(self.Abnormal_BC))
    def NMRassign(self,InputFilePath,threshold):
        a=pd.read_csv(InputFilePath).iloc[:,1].to_list()
        Carbon_label=[]
        for i in a:
            if i==i:
                Carbon_label.append(int(i))
            else:
                Carbon_label.append("nan")
        df_shift=pd.read_csv(self.CarbonListFile)
        shift_list=df_shift.iloc[:,0].tolist()
        type_list=df_shift.iloc[:,1].tolist()
        mol= Chem.MolFromMolFile(self.comp_dir)
        Atoms=mol.GetAtoms()
        carbon_label_all=[]
        for i,atom in enumerate(Atoms):
            #首先获得原子元素符号，判断其是否为碳
            a=atom.GetSymbol()
            #如果是碳进一步判断碳的类型
            if a=="C":
                carbon_label_all.append(i)
        #获得没有匹配的碳编号
        carbon_label_un= [x for x in carbon_label_all if x not in Carbon_label]
        Carbon_type=[]
        for label in carbon_label_un:
            atom=mol.GetAtomWithIdx(label)
            proton_nums=atom.GetTotalNumHs()
            if  proton_nums==3:
                Carbon_type.append("q")
            if  proton_nums==2:
                Carbon_type.append("t")
            if  proton_nums==1:
                Carbon_type.append("d")
            if  proton_nums==0:
                Carbon_type.append("s")
        #然后还需要知道剩下的化学位移：
        shift_un=[]#没有匹配的化学位移
        label_un=[]#没有匹配的化学位移对应的碳类型
        order_un=[]#没有匹配的化学位移的顺序
        for i,label in enumerate(Carbon_label):
            if label=="nan":
                shift_un.append(shift_list[i])
                label_un.append(type_list[i])
                order_un.append(i)
        #现在需要穷举出所有的情况，对每一种情况进行BC和COSY的检查
        #首先穷举出所有情况，即穷举出有多少种Carbon_label的排列情况，固定化学位移不变。
        #以上问题转化为：第一步：n个label_for_CH2只能放在n带有t标签的空中，和m个label_for_CH只能放在m个d标签的空中，一共有多少种情况
        #第二步，把每种情况放进Carbon_label列表中，并进行BC和COSY的检查
        label_for_CH3=[]
        label_for_CH2=[]
        label_for_CH=[]
        label_for_C=[]
        for i,Ctype in enumerate(Carbon_type):
            if Ctype=="q":
                label_for_CH3.append(carbon_label_un[i])
            if Ctype=="t":
                label_for_CH2.append(carbon_label_un[i])
            if Ctype=="d":
                label_for_CH.append(carbon_label_un[i])
            if Ctype=="s":
                label_for_C.append(carbon_label_un[i])
        #由于Carbon_label的顺序是按照CH3，CH2，CH，C的顺序来排列的，下面的顺序也要对应
        #对应之后直接就可以把生成的情况对应放进去了，很方便。
        all_possible=[]
        for i in permutations(label_for_CH3):
            for j in permutations(label_for_CH2):
                for m in permutations(label_for_CH):
                    for n in permutations(label_for_C):
                        single_possible=i+j+m+n
                        all_possible.append(list(single_possible))
        #分别把每种情况放于Carbon_label中,生成所有情况的Carbon_label的顺序Carbon_label_all
        Carbon_label_all=[]
        self.CorrectP_Total_list=[]
        self.CorrectP_BC_list=[]
        self.CorrectP_COSY_list=[]
        
        self.Normal_BC_list=[]
        self.Normal_COSY_list=[]
        self.Normal_BClabel_list=[]
        self.Normal_COSYlabel_list=[]
        
        self.Abnormal_BC_list=[]
        self.Abnormal_COSY_list=[]
        self.Abnormal_BClabel_list=[]
        self.Abnormal_COSYlabel_list=[]
        for i in range(len(all_possible)):
            for j in range(len(order_un)):
                Carbon_label[order_un[j]]=all_possible[i][j]
            #将每种情况都存在列表中
            Carbon_label_all.append(Carbon_label[:])
            self.Carbon_label_all=Carbon_label_all
            #对每种情况进行BC和COSY的检查
            self.BC_check(Carbon_label,threshold)
            self.COSY_check(Carbon_label,threshold)
            
            self.Normal_BC_list.append(self.Normal_BC[:])
            self.Normal_COSY_list.append(self.Normal_COSY[:])

            
            self.Abnormal_BC_list.append(self.Abnormal_BC[:])
            self.Abnormal_COSY_list.append(self.Abnormal_COSY[:])
            self.Abnormal_BClabel_list.append(self.Abnormal_BClabel[:])
            self.Abnormal_COSYlabel_list.append(self.Abnormal_COSYlabel[:])
            #计算每种情况的CorrectP
            self.CorrectP()
            self.CorrectP_Total_list.append(round(self.CorrectP_Total,2))
            self.CorrectP_BC_list.append(round(self.CorrectP_BC,2))
            self.CorrectP_COSY_list.append(round(self.CorrectP_COSY))
    def visualization_check(self):
        #获取分数CorrectP_Total_list前10，并且可视化
        testList =self.CorrectP_Total_list
        tmp = zip(range(len(testList)), testList)
        large10 = heapq.nlargest(10, tmp, key=lambda x:x[1])
        #获得得分前10项
        CorrectP_Total_TOP10=[]
        CorrectP_BC_TOP10=[]
        CorrectP_COSY_TOP10=[]
        labels__TOP10=[]
        for i in range(len(large10)):
            labels__TOP10.append(large10[i][0])
            CorrectP_Total_TOP10.append(large10[i][1])
            CorrectP_BC_TOP10.append(self.CorrectP_BC_list[large10[i][0]])
            CorrectP_COSY_TOP10.append(self.CorrectP_COSY_list[large10[i][0]])
        #开始绘图
        plt.figure(figsize=(4,3))
        x = np.arange(len(labels__TOP10))  # x轴刻度标签位置
        width = 0.25  # 柱子的宽度
        plt.bar(x - width, CorrectP_Total_TOP10, width, label='CorrectP_Total',color="red")
        plt.bar(x, CorrectP_BC_TOP10, width, label='CorrectP_BC',color="green")
        plt.bar(x + width, CorrectP_COSY_TOP10, width, label='CorrectP_COSY',color="blue")
        plt.xticks(x,rotation=45,fontsize=15,labels=labels__TOP10)
        plt.yticks(fontsize=15)
        plt.ylabel('Scores',fontsize=15)
        plt.savefig("visualization_check.png")
        plt.show()
        self.CorrectP_Total_TOP10=CorrectP_Total_TOP10
        self.CorrectP_BC_TOP10=CorrectP_BC_TOP10
        self.CorrectP_COSY_TOP10=CorrectP_COSY_TOP10
        self.labels__TOP10=labels__TOP10
    def abnormal_figure(self,num):
        m= Chem.MolFromMolFile(self.comp_dir)
        view = rdMolDraw2D.MolDraw2DSVG(400, 300)
        tm = rdMolDraw2D.PrepareMolForDrawing(m)
        view.DrawMolecule(tm, highlightAtoms=list(self.Abnormal_COSYlabel_list[num]+self.Abnormal_BClabel_list[num]))
        view.FinishDrawing()
        svg = view.GetDrawingText()
        svg2png(bytestring=svg, write_to="./abnormal_figure.png")
        self.img = Image.open("./abnormal_figure.png")
    def abnormal_figure_check(self):
		    m= Chem.MolFromMolFile(self.comp_dir)
		    view = rdMolDraw2D.MolDraw2DSVG(400, 300)
		    tm = rdMolDraw2D.PrepareMolForDrawing(m)
		    view.DrawMolecule(tm, highlightAtoms=list(self.Abnormal_COSYlabel+self.Abnormal_BClabel))
		    view.FinishDrawing()
		    svg = view.GetDrawingText()
		    svg2png(bytestring=svg, write_to="./abnormal_figure.png")
		    self.img = Image.open("./abnormal_figure.png")
