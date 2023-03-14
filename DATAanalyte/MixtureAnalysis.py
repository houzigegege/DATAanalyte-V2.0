#!/usr/bin/env python
# coding: utf-8

import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import heapq
from adjustText import adjust_text

class mixture_analyte(object):
    def __init__(self,CarbonFileName,CompareResult):
        self.explt_data=pd.read_csv(CarbonFileName)
        self.hit_data=pd.read_csv(CompareResult)
    def mix_analyte(self,num,width=0.5,fontsize=10,x_left=200,x_right=0):
        """Information on known compounds that may be present in the mixture to be analyzed
        num：The number of the compound in the database
        width：The width of the spectral signal 
        """
        #根据数据库数据绘图
        plt.figure(figsize=(12.5, 3))
        height=10
        database=pd.read_excel("Inhouse_database_1.xls",sheet_name="C_DEPT_NMR")
        lib_shift_data=database.iloc[:,num*2-2:num*2-1].iloc[:,0].tolist()
        lib_type_data=database.iloc[:,num*2-1:num*2].iloc[:,0].tolist()
        for n,peak in enumerate(lib_shift_data):
            if str(lib_type_data[n])=="q": 
                plt.bar(peak, height=height,width=width,color="r")

            if str(lib_type_data[n])=="t": 
                plt.bar(peak, height=-height,width=width,color="green")

            if str(lib_type_data[n])=="d": 
                plt.bar(peak, height=height,width=width,color="blue")

            if str(lib_type_data[n])=="s": 
                plt.bar(peak, height=-height,width=width,color="orange")      
        plt.xlim(200,0)
        plt.ylim(-15,15)
        plt.xticks(fontsize =15)
        plt.yticks(fontsize =15)
        plt.hlines(0,xmin=0,xmax=200,linestyle="-", color="black")
        plt.xlabel('ppm', fontsize =15)
        plt.savefig("MixFigure1.png")
        plt.show()
        
        
        #获得待分析化合物的化学位移数据、碳类型数据和pts数据，这些数据通过combine_data生成后，存在"test_data1.csv"中
        explt_shift_data=self.explt_data.iloc[:,0]
        explt_type_data=self.explt_data.iloc[:,1]
        explt_height_data=self.explt_data.iloc[:,2]
        #获得要拆分的数据,即和数据库匹配上的实验数据
        hit_data_one=self.hit_data.iloc[:,num*2-1:num*2].iloc[:,0].tolist()
        #拆分数据
        explt_shift_remain=[]
        explt_type_remain=[]
        explt_height_remain=[]
        explt_shift_del=[]
        explt_type_del=[]
        explt_height_del=[]
        for i,peak in enumerate(explt_shift_data):
            if peak in hit_data_one:
                explt_shift_del.append(peak)
                explt_type_del.append(explt_type_data[i])
                explt_height_del.append(explt_height_data[i])
                hit_data_one.remove(peak)
            else:
                explt_shift_remain.append(peak)
                explt_type_remain.append(explt_type_data[i])
                explt_height_remain.append(explt_height_data[i])
       
    
        #重新绘图
        plt.figure(figsize=(12.5, 3))
        #peaks = ng.peakpick.pick(data, pthres=threshold, algorithm="downward")
        
        #绘制母图
        for i,Ctype in enumerate(explt_type_remain):
            height =explt_height_remain[i]
            ppm=explt_shift_remain[i]
            if Ctype=="q":
                plt.bar(ppm, height,width=width,color="gray")
                #plt.text(ppm,height,round(ppm,1), ha="center", va="center",rotation=90,fontsize=fontsize,color="black")
            if Ctype=="t":
                plt.bar(ppm, -height,width=width,color="gray")
                #plt.text(ppm,-height,round(ppm,1), ha="center", va="center",rotation=90,fontsize=fontsize,color="black")
            if Ctype=="d":
                plt.bar(ppm, height,width=width,color="gray")
                #plt.text(ppm,height,round(ppm,1), ha="center", va="center",rotation=90,fontsize=fontsize,color="black")
            if Ctype=="s":
                plt.bar(ppm, -height,width=width,color="gray")
                #plt.text(ppm,-height,round(ppm,1), ha="center", va="center",rotation=90,fontsize=fontsize,color="black")
        for i,Ctype in enumerate(explt_type_del):
            height =explt_height_del[i]
            ppm=explt_shift_del[i]
            if Ctype=="q":
                plt.bar(ppm, height,width=width,color="r")
                #plt.text(ppm,height,round(ppm,1), ha="center", va="center",rotation=90,fontsize=fontsize,color="black")
            if Ctype=="t":
                plt.bar(ppm, -height,width=width,color="g")
                #plt.text(ppm,-height,round(ppm,1), ha="center", va="center",rotation=90,fontsize=fontsize,color="black")
            if Ctype=="d":
                plt.bar(ppm, height,width=width,color="b")
                #plt.text(ppm,height,round(ppm,1), ha="center", va="center",rotation=90,fontsize=fontsize,color="black")
            if Ctype=="s":
                plt.bar(ppm, -height,width=width,color="orange")
                #plt.text(ppm,-height,round(ppm,1), ha="center", va="center",rotation=90,fontsize=fontsize,color="black")
        plt.hlines(0,0,200,linestyle="-", color="black")
        plt.xlim(x_left, x_right)
        #plt.ylim(y_bottom,y_top)
        plt.xlabel('ppm', fontsize =15)
        #plt.ylabel('Hight', fontsize =15)
        plt.tick_params(labelsize=15)
        plt.savefig("MixFigure2.png")
        plt.show()

        #绘制子图
        plt.figure(figsize=(12.5, 3))
        for i,Ctype in enumerate(explt_type_del):
            height =explt_height_del[i]
            ppm=explt_shift_del[i]
            if Ctype=="q":
                plt.bar(ppm, height,width=width,color="red")
                plt.text(ppm,height,round(ppm,1), ha="center", va="center",rotation=90,fontsize=fontsize,color="black")
            if Ctype=="t":
                plt.bar(ppm, -height,width=width,color="green")
                plt.text(ppm,-height,round(ppm,1), ha="center", va="center",rotation=90,fontsize=fontsize,color="black")
            if Ctype=="d":
                plt.bar(ppm, height,width=width,color="blue")
                plt.text(ppm,height,round(ppm,1), ha="center", va="center",rotation=90,fontsize=fontsize,color="black")
            if Ctype=="s":
                plt.bar(ppm, -height,width=width,color="orange")
                plt.text(ppm,-height,round(ppm,1), ha="center", va="center",rotation=90,fontsize=fontsize,color="black")
        plt.hlines(0,0,200,linestyle="-", color="black")
        plt.xlim(x_left, x_right)
        plt.xlabel('ppm', fontsize =15)
        #plt.ylabel('Hight', fontsize =15)
        plt.tick_params(labelsize=15)
        plt.savefig("MixFigure3.png")
        plt.show()
    
    def mark_known_gray(self,com_list,width=0.5,fontsize=10,x_left=200,x_right=0,shift="False",adjust="False"):
        #获得待分析化合物的化学位移数据、碳类型数据和pts数据，这些数据通过combine_data生成后，存在"test_data1.csv"中
        explt_shift_data=self.explt_data.iloc[:,0]
        explt_type_data=self.explt_data.iloc[:,1]
        explt_height_data=self.explt_data.iloc[:,2]
        #获得要表为灰色的已知化合物的数据
        hit_data_all=[]
        for num in com_list:
            hit_data_one=self.hit_data.iloc[:,num*2-1:num*2].iloc[:,0].tolist()
            hit_data_all= hit_data_all+hit_data_one
        #print(hit_data_all)
        #拆分数据
        explt_shift_remain=[]
        explt_type_remain=[]
        explt_height_remain=[]
        explt_shift_gray=[]
        explt_type_gray=[]
        explt_height_gray=[]
        for i,peak in enumerate(explt_shift_data):
            if peak in hit_data_all:
                explt_shift_gray.append(peak)
                explt_type_gray.append(explt_type_data[i])
                explt_height_gray.append(explt_height_data[i])
                hit_data_all.remove(peak)
            else:
                explt_shift_remain.append(peak)
                explt_type_remain.append(explt_type_data[i])
                explt_height_remain.append(explt_height_data[i]) 
        #绘图
        plt.figure(figsize=(10, 3), dpi=300)
        text_height1=[]
        for i,Ctype in enumerate(explt_type_remain):
            height =explt_height_remain[i]
            ppm=explt_shift_remain[i]
            if Ctype=="q":
                plt.bar(ppm, height,width=width,color="red")
                text_height1.append(height)
            if Ctype=="t":
                plt.bar(ppm, -height,width=width,color="red")
                text_height1.append(-height)
            if Ctype=="d":
                plt.bar(ppm, height,width=width,color="red")
                text_height1.append(height)
            if Ctype=="s":
                plt.bar(ppm, -height,width=width,color="red")
                text_height1.append(-height)
        if shift=="True":
            texts=[plt.text(explt_shift_remain[i],text_height1[i],round(explt_shift_remain[i],1),ha="center", va="center",rotation=90) for i in range(len(explt_shift_remain))]
            if adjust=="True":
                adjust_text(texts, 
                    arrowprops=dict(arrowstyle='->',lw= 1,color='gray'),
                              only_move={'text': 'x'},
                             expand_text=(3, 1),
                            #expand_objects=(3, 1.2),
                            force_text=(0.75, 0), 
                            force_objects=(1, 0),)
        text_height2=[]
        for i,Ctype in enumerate(explt_type_gray):
            height =explt_height_gray[i]
            ppm=explt_shift_gray[i]
            if Ctype=="q":
                plt.bar(ppm, height,width=width,color="gray")
                text_height2.append(height)
            if Ctype=="t":
                plt.bar(ppm, -height,width=width,color="gray")
                text_height2.append(-height)
            if Ctype=="d":
                plt.bar(ppm, height,width=width,color="gray")
                text_height2.append(height)
            if Ctype=="s":
                plt.bar(ppm, -height,width=width,color="gray")
                text_height2.append(-height)
        if shift=="True":
            texts=[plt.text(explt_shift_gray[i],text_height2[i],round(explt_shift_gray[i],1),ha="center", va="center",rotation=90) for i in range(len(explt_shift_gray))]
            if adjust=="True":
                adjust_text(texts, 
                    arrowprops=dict(arrowstyle='->',lw= 1,color='gray'),
                              only_move={'text': 'x'},
                             expand_text=(3, 1),
                            #expand_objects=(3, 1.2),
                            force_text=(0.75, 0), 
                            force_objects=(1, 0),)
        plt.hlines(0,0,200,linestyle="-", color="black")
        plt.xlim(x_left, x_right)
        #plt.ylim(y_bottom,y_top)
        plt.xlabel('ppm', fontsize =15)
        #plt.ylabel('Hight', fontsize =15)
        plt.tick_params(labelsize=15)
        plt.show()
    def move_known(self,com_list,width=0.5,fontsize=10,x_left=200,x_right=0,shift="False",adjust="False"):
        #获得待分析化合物的化学位移数据、碳类型数据和pts数据，这些数据通过combine_data生成后，存在"test_data1.csv"中
        explt_shift_data=self.explt_data.iloc[:,0]
        explt_type_data=self.explt_data.iloc[:,1]
        explt_height_data=self.explt_data.iloc[:,2]
        #获得要标为灰色的已知化合物的数据
        hit_data_all=[]
        for num in com_list:
            hit_data_one=self.hit_data.iloc[:,num*2-1:num*2].iloc[:,0].tolist()
            hit_data_all= hit_data_all+hit_data_one
        #拆分数据
        explt_shift_remain=[]
        explt_type_remain=[]
        explt_height_remain=[]
        explt_shift_gray=[]
        explt_type_gray=[]
        explt_height_gray=[]
        for i,peak in enumerate(explt_shift_data):
            if peak in hit_data_all:
                explt_shift_gray.append(peak)
                explt_type_gray.append(explt_type_data[i])
                explt_height_gray.append(explt_height_data[i])
                hit_data_all.remove(peak)
            else:
                explt_shift_remain.append(peak)
                explt_type_remain.append(explt_type_data[i])
                explt_height_remain.append(explt_height_data[i]) 
        #绘图
        plt.figure(figsize=(12.5, 3))
        text_height=[]
        for i,Ctype in enumerate(explt_type_remain):
            height =explt_height_remain[i]
            ppm=explt_shift_remain[i]
            if Ctype=="q":
                plt.bar(ppm, height,width=width,color="red")
                text_height.append(height)
            if Ctype=="t":
                plt.bar(ppm, -height,width=width,color="g")
                text_height.append(-height)
            if Ctype=="d":
                plt.bar(ppm, height,width=width,color="b")
                text_height.append(height)
            if Ctype=="s":
                plt.bar(ppm, -height,width=width,color="orange")
                text_height.append(-height)
        if shift=="True":
            texts=[plt.text(explt_shift_remain[i],text_height[i],round(explt_shift_remain[i],1),ha="center", va="center",rotation=90) for i in range(len(explt_shift_remain))]
            if adjust=="True":
                adjust_text(texts, 
                    arrowprops=dict(arrowstyle='->',lw= 1,color='gray'),
                              only_move={'text': 'x'},
                             expand_text=(3, 1),
                            #expand_objects=(3, 1.2),
                            force_text=(0.75, 0), 
                            force_objects=(1, 0),)
        plt.hlines(0,0,200,linestyle="-", color="black")
        plt.xlim(x_left, x_right)
        #plt.ylim(y_bottom,y_top)
        plt.xlabel('ppm', fontsize =15)
        #plt.ylabel('Hight', fontsize =15)
        plt.tick_params(labelsize=15)
        plt.savefig("MixFigure4.png")
        plt.show()
    #def split_by_indensity(self,threshold=):
    def high_content_split(self,width=0.5,fontsize=10,x_left=200,x_right=0,shift="False",adjust="False"):
        #获得待分析化合物的化学位移数据、碳类型数据和pts数据，这些数据通过combine_data生成后，存在"test_data1.csv"中
        explt_shift_data=self.explt_data.iloc[:,0].tolist()
        explt_type_data=self.explt_data.iloc[:,1].tolist()
        explt_height_data=self.explt_data.iloc[:,2].tolist()
       
        #获得CH，CH2，CH3的参考高度，选择最高的一个信号
        refer1=max(explt_height_data)
        #获得C的参考高度：
        height_of_C=[]
        for i,type_C in enumerate(explt_type_data):
            if type_C=="s":
                height_of_C.append(explt_height_data[i])
        refer2=max(height_of_C)
        #根据参考高度拆分图谱
        #先拆分CH，CH2，CH3
        explt_shift_remain1=[]
        explt_type_remain1=[]
        explt_height_remain1=[]
        explt_shift_split=[]
        explt_type_split=[]
        explt_height_split=[]
        for i,height in enumerate(explt_height_data):
            if height>0.6*refer1:
                explt_shift_split.append(explt_shift_data[i])
                explt_type_split.append(explt_type_data[i])
                explt_height_split.append(explt_height_data[i])
            else:
                explt_shift_remain1.append(explt_shift_data[i])
                explt_type_remain1.append(explt_type_data[i])
                explt_height_remain1.append(explt_height_data[i])
        #然后拆分C
        C_ppm=[]
        C_type=[]
        C_height=[]
        explt_shift_remain2=[]
        explt_type_remain2=[]
        explt_height_remain2=[]
        for i,type_C in enumerate(explt_type_remain1):
            if type_C=="s" and explt_height_remain1[i]>0.6*refer2:
                explt_shift_split.append(explt_shift_remain1[i])
                explt_type_split.append(explt_type_remain1[i])
                explt_height_split.append(explt_height_remain1[i])
            else:
                explt_shift_remain2.append(explt_shift_remain1[i])
                explt_type_remain2.append(explt_type_remain1[i])
                explt_height_remain2.append(explt_height_remain1[i])
        #绘图
        plt.figure(figsize=(10, 3), dpi=300)
        text_height=[]
        for i,Ctype in enumerate(explt_type_remain2):
            height =explt_height_remain2[i]
            ppm=explt_shift_remain2[i]
            
            if Ctype=="q":
                plt.bar(ppm, height,width=width,color="red")
                text_height.append(height)
            if Ctype=="t":
                plt.bar(ppm, -height,width=width,color="g")
                text_height.append(-height)
            if Ctype=="d":
                plt.bar(ppm, height,width=width,color="b")
                text_height.append(height)
            if Ctype=="s":
                plt.bar(ppm, -height,width=width,color="orange")
                text_height.append(-height)
        if shift=="True":
            texts=[plt.text(explt_shift_remain2[i],text_height[i],round(explt_shift_remain2[i],1),ha="center", va="center",rotation=90) for i in range(len(explt_shift_remain2))]
            if adjust=="True":
                adjust_text(texts, 
                    arrowprops=dict(arrowstyle='->',lw= 1,color='gray'),
                              only_move={'text': 'x'},
                             expand_text=(3, 1),
                            #expand_objects=(3, 1.2),
                            force_text=(0.75, 0), 
                            force_objects=(1, 0),
                    )
        plt.hlines(0,0,200,linestyle="-", color="black")
        plt.xlim(x_left, x_right)
        #plt.ylim(y_bottom,y_top)
        plt.xlabel('ppm', fontsize =15)
        #plt.ylabel('Hight', fontsize =15)
        plt.tick_params(labelsize=15)
        plt.show()
        
        plt.figure(figsize=(10, 3), dpi=300)
        text_height1=[]
        for i,Ctype in enumerate(explt_type_split):
            height =explt_height_split[i]
            ppm=explt_shift_split[i]
            if Ctype=="q":
                plt.bar(ppm, height,width=width,color="red")
                text_height1.append(height)
            if Ctype=="t":
                plt.bar(ppm, -height,width=width,color="g")
                text_height1.append(-height)
            if Ctype=="d":
                plt.bar(ppm, height,width=width,color="b")
                text_height1.append(height)
            if Ctype=="s":
                plt.bar(ppm, -height,width=width,color="orange")
                text_height1.append(-height)
        if shift=="True":
            texts=[plt.text(explt_shift_split[i],text_height1[i],round(explt_shift_split[i],1),ha="center", va="center",rotation=90) for i in range(len(explt_shift_remain2))]
            if adjust=="True":
                adjust_text(texts, 
                    arrowprops=dict(arrowstyle='->',lw= 1,color='gray'),
                              only_move={'text': 'x'},
                             expand_text=(3, 1),
                            #expand_objects=(3, 1.2),
                            force_text=(0.75, 0), 
                            force_objects=(1, 0),)
        plt.hlines(0,0,200,linestyle="-", color="black")
        plt.xlim(x_left, x_right)
        #plt.ylim(y_bottom,y_top)
        plt.xlabel('ppm', fontsize =15)
        #plt.ylabel('Hight', fontsize =15)
        plt.tick_params(labelsize=15)
        plt.show()
        
