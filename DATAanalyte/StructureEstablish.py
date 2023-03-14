#!/usr/bin/env python
# coding: utf-8



#!/usr/bin/python3
# -*- coding: utf-8 -*-
import string
import itertools
import copy
import os
import sys
import configparser
import numpy as np
import networkx as nx
from rdkit import rdBase, Chem
from rdkit.Chem import AllChem, Draw,CombineMols
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
from matplotlib.colors import ColorConverter
import networkx as nx
import argparse
import multiprocessing
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions 
import pandas as pd
import heapq
import matplotlib.pyplot as plt
import math
from io import BytesIO
from PIL import Image
from cairosvg import svg2png

class GraphNode:
    def __init__(self, name: string, num: int):
        self.m_conNodeIndex: list[int] = list()
        self.m_name: string = name
        self.m_num: int = num

    def __str__(self):
        return '<' + self.m_name + ',' + str(self.m_num) + '>'
def printNodesComb(snodes: list):
    snodeslen = len(snodes)
    lvertex = list()
    for snode in snodes:
        vertex = [0] * len(snodes)
        for nodeindex in snode.m_conNodeIndex:
            vertex[nodeindex] = 1
        lvertex.append(vertex)
#     print('-' * 3 * snodeslen)
#     for vertex in lvertex:
#         print(vertex)
#     print('-' * 3 * snodeslen)
    return(lvertex)
def searchNodesComb(snodes: list, nodeindex: int = 0):
    # 计算当前节点剩余连接数量
    num = snodes[nodeindex].m_num - len(snodes[nodeindex].m_conNodeIndex)
    # 计算未连接完成的节点
    lnodes = list()
    for i, snode in enumerate(snodes):
        if nodeindex == i:
            continue
        snodenum = snode.m_num - len(snode.m_conNodeIndex)
        if snodenum > 0:
            lnodes.append(snode)
    if  len(lnodes) > 0:
        # 求出剩余节点组合列表
        rltnodes = list(itertools.combinations(lnodes, num))
        for rtnodes in rltnodes:
            snodesex = copy.deepcopy(snodes)
            for node in list(rtnodes):
                for i in range(len(snodes)):
                    if snodesex[i].m_name == node.m_name:
                        snodesex[nodeindex].m_conNodeIndex.append(i)
                        snodesex[i].m_conNodeIndex.append(nodeindex)

            if nodeindex + 1 < len(snodesex):
                searchNodesComb(snodesex, nodeindex + 1)
    else:
        if len(lnodes) == 0:
            if num == 0:
                result=printNodesComb(snodes)
                AdjacencyMatrixResults.append(result)
                
def inputData(strnodename,strnodenum):
    #定义一个全局变量来获取最后的结果
    global  AdjacencyMatrixResults
    AdjacencyMatrixResults=[]
    # 获取顶点名及度
    strnodename = strnodename.replace(' ', '')
    strnodenum = strnodenum.replace(' ', '')
    nodename = strnodename.split(',')
    nodenum = strnodenum.split(',')
    bdigit = True
    for num in nodenum:
        if not num.isdigit():
            bdigit = False
            break
    snodes = list()
    for i in range(len(nodename)):
        node = GraphNode(nodename[i], int(nodenum[i]))
        snodes.append(node)
    # 查找所有可能的组合
    searchNodesComb(snodes)
    
    
#将邻接矩阵转化为距离矩阵，同时把邻接矩阵转化为化学结构
def MolFromGraphs(node_list, adjacency_matrix):
    # create empty editable mol object
    mol = Chem.RWMol()

    # add atoms to mol and keep track of index
    node_to_idx = {
            }
    for i in range(len(node_list)):
        a = Chem.Atom(node_list[i])
        molIdx = mol.AddAtom(a)
        node_to_idx[i] = molIdx

    # add bonds between adjacent atoms
    for ix, row in enumerate(adjacency_matrix):
        for iy, bond in enumerate(row):

            # only traverse half the matrix
            if iy <= ix:
                continue

            # add relevant bond type (there are many more of these)
            if bond == 0:
                continue
            elif bond == 1:
                bond_type = Chem.rdchem.BondType.SINGLE
                mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
    # Convert RWMol to Mol object
    mol = mol.GetMol()
    smi = Chem.MolToSmiles(mol)
#     opts = DrawingOptions()
#     opts.includeAtomNumbers=True
#     fig=Draw.MolToImage(mol,options=opts,size=(300,300))
    distance = np.array(Chem.GetDistanceMatrix(mol), "d")
    return(distance,smi)

def BCsimplify(chemical_shift,threshold):
    df1=pd.read_csv("cor_data.csv")
    BC_H2C=df1["BC_H2C"].tolist()
    BC_C2C=df1["BC_C2C"].tolist()
    BC_H2C= [a_ for a_ in BC_H2C if a_ == a_]
    BC_C2C= [a_ for a_ in BC_C2C if a_ == a_]
    BC_H2C_simp=[]
    BC_C2C_simp=[]
    for i in range(len(BC_H2C)):
        delta_H2C=[]
        delta_C2C=[]
        #首先计算每组和化学位移的差值
        for j in range(len(chemical_shift)):
            a=abs(BC_H2C[i]-chemical_shift[j])
            b=abs(BC_C2C[i]-chemical_shift[j])
            delta_H2C.append(a)
            delta_C2C.append(b)
        if min(delta_H2C)<=threshold and min(delta_C2C)<=threshold:
            BC_H2C_simp.append(BC_H2C[i])
            BC_C2C_simp.append(BC_C2C[i])
    return(BC_H2C_simp,BC_C2C_simp)


def BC_check(chemical_shift,BC_H2C_simp,BC_C2C_simp,DistanceMatrix,threshold):
    """chemical_shift,列表,化学位移值
       BC_H2C_simp/BC_C2C_simp,简化后的二维相关信息
       DistanceMatrix 化合物的距离矩阵
    """
    #获得化合物的距离矩阵
    DistanceMatrix =pd.DataFrame(DistanceMatrix)
    #获得化合物的BC相关信息
    BC_H2C=BC_H2C_simp
    BC_C2C=BC_C2C_simp
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
        #获取最小值的索引
        min_H2C=delta_H2C.index(min(delta_H2C))
        min_C2C=delta_C2C.index(min(delta_C2C))
        #获取最小值对应的化学位移
        shift_H2C=chemical_shift[min_H2C]
        shift_C2C=chemical_shift[min_C2C]
        #判断最小差值小于阈值,表示这组BC信号和一维的碳信号匹配上了，则去距离矩阵中寻找距离值
        if min(delta_H2C)<=threshold and min(delta_C2C)<=threshold:
            C_distance1=DistanceMatrix.iloc[min_H2C,min_C2C]
            #print(i,C_distance1)
            if C_distance1<=2:
                Normal.append([min_H2C,min_C2C,C_distance1,shift_H2C,shift_C2C])
            else:
                Abnormal.append([min_H2C,min_C2C,C_distance1,shift_H2C,shift_C2C])
    Abnormal=pd.DataFrame(Abnormal)
    Normal=pd.DataFrame(Normal)
    Pcorr=len(Normal)/(len(Normal)+len(Abnormal))
    return(Pcorr)
def MarkImpossible(ShiftList,ComfirmList,DistanceMatrix):
    DistanceMatrix =pd.DataFrame(DistanceMatrix)
    for i in range(len(ComfirmList)):
        a=ShiftList.index(ComfirmList[i][0])
        b=ShiftList.index(ComfirmList[i][1])
        distance=DistanceMatrix.iloc[a,b]
        if distance!=1:
            label="False"
            break
        else:
            label="True"
    return(label)      
#把结果以网络图的形式展示
def DrawNetwork(AdjacencyMatrixResults,num):
    a=np.array(AdjacencyMatrixResults)
    G=nx.Graph()
    G=nx.from_numpy_matrix(np.array(a[num]))
    position = nx.circular_layout(G)
    nx.draw_networkx_nodes(G,position, node_color="yellow",edgecolors="yellow")
    nx.draw_networkx_edges(G,position,width=3.0)
    nx.draw_networkx_labels(G,position)
    
#将邻接矩阵转化为距离矩阵，同时把邻接矩阵转化为化学结构
def StructureLabel(node_list, adjacency_matrix):
    # create empty editable mol object
    mol = Chem.RWMol()
    # add atoms to mol and keep track of index
    node_to_idx = {
            }
    for i in range(len(node_list)):
        a = Chem.Atom(node_list[i])
        molIdx = mol.AddAtom(a)
        node_to_idx[i] = molIdx

    # add bonds between adjacent atoms
    for ix, row in enumerate(adjacency_matrix):
        for iy, bond in enumerate(row):

            # only traverse half the matrix
            if iy <= ix:
                continue

            # add relevant bond type (there are many more of these)
            if bond == 0:
                continue
            elif bond == 1:
                bond_type = Chem.rdchem.BondType.SINGLE
                mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
    # Convert RWMol to Mol object
    mol = mol.GetMol() 
    opts = DrawingOptions()
    opts.atomLabelFontSize=500
    opts.includeAtomNumbers=True
    opts.bondLineWidth=2
    opts.dblBondLengthFrac=0.8
    opts.dotsPerAngstrom=1000
    fig=Draw.MolToImage(mol,options=opts,size=(300,300))
    return(fig,mol)

def add_atom(mol, idx,AtomType=9):
    '''给定分子与原子索引，在该原子上加一个原子
    '''
    mol = Chem.RWMol(mol)
    for idx in idx:
        mol.AddAtom(Chem.Atom(AtomType))
        mol.AddBond(mol.GetNumAtoms()-1,idx,Chem.BondType.SINGLE)
    return Chem.RemoveHs(mol)
def add_atom_single(mol, idx,AtomType=9):
    '''给定分子与原子索引，在该原子上加一个原子
    '''
    mol = Chem.RWMol(mol)
    mol.AddAtom(Chem.Atom(AtomType))
    mol.AddBond(mol.GetNumAtoms()-1,idx,Chem.BondType.SINGLE)
    return Chem.RemoveHs(mol)
def add_atom_double(mol, idx,AtomType=9):
    '''给定分子与原子索引，在该原子上加一个原子
    '''
    mol = Chem.RWMol(mol)
    mol.AddAtom(Chem.Atom(AtomType))
    mol.AddBond(mol.GetNumAtoms()-1,idx,Chem.BondType.DOUBLE)
    return Chem.RemoveHs(mol)
#函数一，获取marker邻居原子的index, 注意marker只能是一个单键连接核上的原子，否则邻居会多于一个
def get_neiid_bysymbol(mol,marker):
    try:
        for atom in mol.GetAtoms():
            if atom.GetSymbol()==marker:
                neighbors=atom.GetNeighbors()
                if len(neighbors)>1:
                    print ('Cannot process more than one neighbor, will only return one of them')
                atom_nb=neighbors[0]
                return atom_nb.GetIdx()
    except Exception as e:
        print (e)
        return None
#函数二，获取marker原子的index
def get_id_bysymbol(mol,marker):
    for atom in mol.GetAtoms():
        if atom.GetSymbol()==marker:
            return atom.GetIdx()

def combine2frags(mol_a,mol_b,maker_a,maker_b):
    """Combine two fragments"""
    #将两个待连接分子置于同一个对象中
    merged_mol = Chem.CombineMols(mol_a,mol_b)
    bind_pos_a=get_neiid_bysymbol(merged_mol,maker_a)
    bind_pos_b=get_neiid_bysymbol(merged_mol,maker_b)
    #转换成可编辑分子，在两个待连接位点之间加入单键连接，特殊情形需要其他键类型的情况较少，需要时再修改
    ed_merged_mol= Chem.EditableMol(merged_mol)
    ed_merged_mol.AddBond(bind_pos_a,bind_pos_b,order=Chem.rdchem.BondType.SINGLE)
    #将图中多余的marker原子逐个移除，先移除marker a
    marker_a_idx=get_id_bysymbol(merged_mol,maker_a)
    ed_merged_mol.RemoveAtom(marker_a_idx)
    #marker a移除后原子序号变化了，所以又转换为普通分子后再次编辑，移除marker b
    temp_mol = ed_merged_mol.GetMol()
    marker_b_idx=get_id_bysymbol(temp_mol,maker_b)
    ed_merged_mol=Chem.EditableMol(temp_mol)
    ed_merged_mol.RemoveAtom(marker_b_idx)
    final_mol = ed_merged_mol.GetMol()
    return final_mol
def combine1frags(merged_mol,maker_a,maker_b):
    """Combine one fragments"""
    #将两个待连接分子置于同一个对象中
    bind_pos_a=get_neiid_bysymbol(merged_mol,maker_a)
    bind_pos_b=get_neiid_bysymbol(merged_mol,maker_b)
    #转换成可编辑分子，在两个待连接位点之间加入单键连接，特殊情形需要其他键类型的情况较少，需要时再修改
    ed_merged_mol= Chem.EditableMol(merged_mol)
    ed_merged_mol.AddBond(bind_pos_a,bind_pos_b,order=Chem.rdchem.BondType.SINGLE)
    #将图中多余的marker原子逐个移除，先移除marker a
    marker_a_idx=get_id_bysymbol(merged_mol,maker_a)
    ed_merged_mol.RemoveAtom(marker_a_idx)
    #marker a移除后原子序号变化了，所以又转换为普通分子后再次编辑，移除marker b
    temp_mol = ed_merged_mol.GetMol()
    marker_b_idx=get_id_bysymbol(temp_mol,maker_b)
    ed_merged_mol=Chem.EditableMol(temp_mol)
    ed_merged_mol.RemoveAtom(marker_b_idx)
    final_mol = ed_merged_mol.GetMol()
    return final_mol


# In[14]:


class StructureEstablish():
    """strnodename: The chimical shift list of the carbon;
       strnodenum: the degree list;
       NMR2DFlie: the files containing BC correlation information;
       ComfirmList: the carbon you can sure linked to eath other, write as [(121.1,23.3),(122.1,44.5)]"""
    def __init__(self,strnodename,strnodenum,NMR2DFlie,ComfirmList):
        self.strnodename=strnodename
        self.strnodenum=strnodenum
        self.NMR2DFlie=NMR2DFlie
        self.ComfirmList=ComfirmList
    def AllStructure(self,threshold):
        strnodename=self.strnodename
        strnodenum=self.strnodenum
        ComfirmList=self.ComfirmList
        #第一步：穷举出所有的邻接矩阵
        inputData(strnodename,strnodenum)
        if len(AdjacencyMatrixResults)==0:
            print("根据输入的信息无法拼成一个化合物碎片，请重新输入")
        else:
            #准备：简化BC相关的信号，为了减少计算量。只筛选出含有输入化学位移的BC信号
            ShiftList=strnodename.split(',')
            ShiftList=[float(i) for i in ShiftList]
            BC_H2C_simp,BC_C2C_simp=BCsimplify(ShiftList,threshold)
            
            #计算node_list表示有多少个碳
            node_list=[]
            for i in range(len(AdjacencyMatrixResults[0])):
                node_list.append(6)
            #第二步：分别把每个邻接矩阵转化为距离矩阵,首先去除掉肯定不可能的情况
            AdjacencyMatrixRemain=[]
            DistanceMatrixRemain=[]
            SmileRemain=[]
            for i in range(len(AdjacencyMatrixResults)):
                #2.1把邻接矩阵转化为距离矩阵
                DistanceMatrix,smi=MolFromGraphs(node_list,AdjacencyMatrixResults[i])
                #2.3 去除不可能的情况（COSY和经验）ComfirmList：二维数组
                if len(ComfirmList)==0:
                    label="True"
                else:
                    label=MarkImpossible(ShiftList,ComfirmList,DistanceMatrix)
                if label=="True":
                    AdjacencyMatrixRemain.append(AdjacencyMatrixResults[i])
                    DistanceMatrixRemain.append(DistanceMatrix)
                    SmileRemain.append(smi)
            #第三步：根据BC相关给剩下的情况打分
            CheckResult=[]
            for i in range(len(DistanceMatrixRemain)):
                Pcorr=BC_check(ShiftList,BC_H2C_simp,BC_C2C_simp,DistanceMatrixRemain[i],threshold)
                CheckResult.append(Pcorr)
            self.BC_H2C_simp=BC_H2C_simp
            self.BC_C2C_simp=BC_C2C_simp
            self.CheckResult=CheckResult
            self.ShiftList=ShiftList
            self.AdjacencyMatrixRemain=AdjacencyMatrixRemain
            self.SmileRemain=SmileRemain
    #第三步：展示出前num个的结构，连接网络图，及Pccor
    def ResultDisplay(self,num):
        """num:display the result of Top(num)"""
        testList =self.CheckResult
        tmp = zip(range(len(testList)), testList)
        TopNum= heapq.nlargest(num, tmp, key=lambda x:x[1])
        plt.figure(figsize=(6,4))
        for i in range(len(TopNum)):
            x=str(TopNum[i][0])
            y=TopNum[i][1]
            plt.bar("Situ"+str(x),y,color="black")
            #plt.text(x,y+0.5,y,ha='center',color="red")
            plt.xticks(rotation=45)
        plt.savefig("establishfig1.png")
        plt.tick_params(labelsize=20)
    def CheckFragmentOne(self,num,label=False):
        """num:the index of the one to be show"""
        plt.subplot(1,2,1)
        DrawNetwork(self.AdjacencyMatrixRemain,num)
        plt.subplot(1,2,2)
        mol=Chem.MolFromSmiles(self.SmileRemain[num])
        print("Chemical shift:")
        print(self.ShiftList)
        print("smiles of the fragment:")
        print(self.SmileRemain[num])
        Draw.MolToFile(mol,"structure.png",size=(200, 200))
        img=plt.imread("structure.png")
        plt.imshow(img)
        plt.axis("off")
        plt.show()
    def ShowFragment(self,num):
        print("Chemical shift:")
        print(self.ShiftList)
        print("smiles of the fragment:")
        print(self.SmileRemain[num])
        node_list=[]
        for i in range(len(self.AdjacencyMatrixRemain[0])):
            node_list.append(6)
        self.fig,mol=StructureLabel(node_list, self.AdjacencyMatrixRemain[num])
        self.fig.save("molecule1.png")

        
    def CombineFrag(self,num,MolForKnown,SiteNew,SiteKnown,AtomNumbers=False,
                   LabelFontSize=500,bondLineWidth=2):
        """combine the two fragments
        num：the index of unknown fragment;
        MolForKnown：the Mol file of known fragment;
        SiteNew: the positions of the sites to be combined, written as (1,2)"""
        #第一步，在计算得到的片段上添加标记
        node_list=[]
        for i in range(len(self.AdjacencyMatrixRemain[0])):
            node_list.append(6)
        self.fig,mol=StructureLabel(node_list, self.AdjacencyMatrixRemain[num])
        MolForNew=add_atom(mol, idx=(SiteNew[0],),AtomType=9)
        MolForNew=add_atom(MolForNew, idx=(SiteNew[1],),AtomType=19)
        #第二步，在已知片段上添加标记,两个标记一定不能一样，因为这两个片段合并在一起以后没法再索引再区分了，所以不能一样
        MolForKnown=add_atom(MolForKnown, idx=(SiteKnown[0],),AtomType=11)
        MolForKnown=add_atom(MolForKnown, idx=(SiteKnown[1],),AtomType=37)
        #mol = mol.GetMol() 
        #第三步，拼接两个片段,得到一个大片段
        FinalMol=combine2frags(MolForNew,MolForKnown,maker_a='F',maker_b='Na')
        #第四步，缝合大片段
        FinalMol=combine1frags(FinalMol,maker_a="K",maker_b="Rb")
        #绘图
        opts = DrawingOptions()
        opts.atomLabelFontSize= LabelFontSize
        opts.includeAtomNumbers=AtomNumbers
        opts.bondLineWidth=bondLineWidth
        opts.dblBondLengthFrac=0.8
        opts.dotsPerAngstrom=1000
        self.fig1=Draw.MolToImage(FinalMol,options=opts,size=(300,300))
        self.CombineMol=FinalMol
        self.FinalSmiles=Chem.MolToSmiles(FinalMol)
        print("smiles of the combined structure:", self.FinalSmiles)
    def AddOther(self,mol,ListAtom,DoubleBondList):
        """mol:the mol file of compound;
        ListAtom:[(a,b，c),] a:the site b:the atom type; c：the bond type"""
        #根据ListAtom注释信息添加原子
        if len(ListAtom)>0:
            for i in range(len(ListAtom)):
                if ListAtom[i][2]=="single":
                    mol=add_atom_single(mol, idx=ListAtom[i][0],AtomType=ListAtom[i][1])
                elif ListAtom[i][2]=="double":
                    mol=add_atom_double(mol, idx=ListAtom[i][0],AtomType=ListAtom[i][1])
            #根据DoubleBondList添加双键
#         mol = Chem.RWMol(mol)
#         if len(DoubleBondList)>0:
#             for i in range(len(DoubleBondList)):
#                 mol.AddBond(DoubleBondList[i][0], DoubleBondList[i][1], Chem.BondType.SINGLE)
        #绘图
        opts = DrawingOptions()
        opts.atomLabelFontSize=500
        opts.includeAtomNumbers=True
        opts.bondLineWidth=2
        opts.dblBondLengthFrac=0.8
        opts.dotsPerAngstrom=1000
        self.fig2=Draw.MolToImage(mol,options=opts,size=(300,300))
        FinalSmiles=Chem.MolToSmiles(mol)
        print("smiles of the combined structure:", FinalSmiles)
        self.AddOtherAtomMol=mol
    def ShowTopNFragment(self,num):
        testList =self.CheckResult
        tmp = zip(range(len(testList)), testList)
        TopNum= heapq.nlargest(num, tmp, key=lambda x:x[1])
        #plt.figure(figsize=(6,4),dpi=300)
        node_list=[]
        for i in range(len(task.AdjacencyMatrixRemain[0])):
            node_list.append(6)
        for i in range(len(TopNum)):
            plt.subplot(math.ceil(num/3),3,i+1)
            fig=StructureLabel(node_list, task.AdjacencyMatrixRemain[TopNum[i][0]])
            fig.save("molecule.png")
            img=plt.imread("molecule.png")
            plt.imshow(img)
            plt.axis("off")
