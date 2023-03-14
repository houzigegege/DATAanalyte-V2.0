#!/usr/bin/env python
# coding: utf-8

# In[2]:


#!/usr/bin/env python
# coding: utf-8

# -*- coding: utf-8 -*-
"""
Created on Wed May 11 09:14:34 2022
@author: Administrator
"""
import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import heapq
import PySimpleGUI as sg
import math
import matplotlib.cm
import networkx as nx
from PIL import Image
import os


os.chdir("D:\DATAanalyte2.0")

import DATAanalyte
from DATAanalyte.NMR1D import *
from DATAanalyte.MixtureAnalysis import *
from DATAanalyte.NMR2D import *
from DATAanalyte.NMRCheck import *
from DATAanalyte.StructureEstablish import *

def combine_cor(list1,list2):
    list_combine=[]
    for i in range(len(list1)):
        a=str(list1[i])+"/"+str(list2[i])
        list_combine.append(a)
    return(list_combine)

#-------------------------------------页面布局————————————————————————
menu_def=[
    ['File','open file::openkey'],
    ['Edit','edit1'],
    ['Help','about']
     ]
layout11=[
          [sg.T("STEP1:Please select/enter the path of the NMR data:")],
          [sg.FolderBrowse(button_text="13C_NMR"),sg.In("",size=(30,1))],
          [sg.FolderBrowse(button_text="DEPT90"),sg.In("",size=(30,1))],
          [sg.FolderBrowse(button_text="DEPT135"),sg.In("",size=(30,1))],
          #step2
          [sg.T("STEP2:get peak")],
          [sg.T("threshold_C:"),sg.In(3000000.0,size=(20,1),key="-threshold_C-")],
          [sg.T("threshold_C90:"),sg.In(90000000.0,size=(20,1),key="-threshold_C90-")],
          [sg.T("threshold_C135pos:"),sg.In(50000000.0,size=(20,1),key="-threshold_C135pos-")],
          [sg.T("threshold_C135neg:"),sg.In(-50000000.0,size=(20,1),key="-threshold_C135neg-")],
          [sg.B("Submit",key="-submit1-"),sg.B("Combine",key="-combine-")],
          #step3
          [sg.T("STEP3:solvent removal:")],
          #[sg.T("sovent:"),sg.In("Chloroform",key="-solvent-",size=(10,1)),sg.T("type:"),sg.In("C",size=(10,1),key="-type-")],
          [sg.T("sovent:"),sg.Combo(["Chloroform","Methanol"],key="-solvent-",size=(10,1)),
           sg.T("type:"),sg.Combo(["CH3","CH2","CH","C"],size=(10,1),key="-type-"),
           sg.B("Submit",key=("-submit2-"))],
          #step4
          [sg.T("STEP4:impurity removal:")],
          [sg.T("type"),sg.Combo(["CH3","CH2","CH","C"],size=(5,1),key="-typeimpurity-"),
           sg.T("threshold:"),sg.In("",size=(5,1),key="-thresholdimpurity1-"),
           sg.In("1e7",size=(5,1),key="-thresholdimpurity2-"),
           sg.B("Submit",key=("-submit3-"))],
          #step5
          [sg.T("STEP5:impurity removal by label:")],
          [sg.T("type"),sg.Combo(["CH3","CH2","CH","C"],size=(10,1),key="-typeimpurity-"),
           sg.T("labels"),sg.In("",size=(10,1),key="-thresholdimpurity-"),
           sg.B("Submit",key="-submit4-")],
          #step6
          [sg.T("STEP6:data analysis:")],
          [sg.B("13C NMR",key="-CNMRcompare-"),sg.B("13C DEPT NMR",key="-CDEPTNMRcompare-")],
          #step7
          [sg.T("STEP7:select compare:"),sg.In(size=(10,1),key="-compname-"),sg.B("Submit",key="-submit5-")],
          #mixture 
          [sg.T("____________________________________________________________")],
          [sg.T("Mixture:highlight compound"),sg.In(size=(10,1),key="-compnum-"),sg.B("Submit",key="-submit6-")],
          [sg.T("Mixture:move known"),sg.In(size=(10,1),key="-complist-"),sg.B("Submit",key="-submit7-")]
         ]
          
layout12=[[sg.T("The analysis result of 13C DEPT NMR",font=("Arial",16),background_color="white",text_color="black",pad=((250,10),(20,0)))],
          [sg.Image("./fig1.png",tooltip=("13C NMR, DEPT90 and DEPT135 peaks of the sample"),key="-fig1-",pad=((60,0),(0,0),),
          background_color='white',right_click_menu=["1",['save as::-fig1-']]),
          sg.Image("./fig3.png",tooltip=("Top 10 compare result of the sample "),background_color='white',key="-fig3-",pad=((0,10),(0,0)),right_click_menu=["1",['save as::-fig3-']])],
          #[sg.T("From:"),sg.In(200,size=(10,1)),sg.T("To:"),sg.In(0,size=(10,1)),sg.B("Set")],
          [sg.Image("./fig2.png",tooltip=("The combined 13C DEPT NMR.red:CH3; green:CH2; blue:CH; orange:C"),key="-fig2-",pad=((0,0),(0,0)),background_color='white',right_click_menu=["1",['save as::-fig2-',"properties::-fig3-"]])],
          #[sg.T("From:"),sg.In(200,size=(10,1)),sg.T("To:"),sg.In(0,size=(10,1)),sg.B("Set")],
          [sg.Image("./fig4.png",tooltip=("Compare result of between sample and the compound from database"),key="-fig4-",pad=((0,0),(0,0)),background_color='white',right_click_menu=["1",['save as::-fig4-']])],
          [sg.T("---------------------------------------------------------------------------------------------Mixture Analysis Result-------------------------------------------------------------------------")],
          [sg.Image("./MixFigure1.png",tooltip=("The selected compound from database"),key="-MixFigure1-",pad=((0,0),(0,0)),background_color='white',right_click_menu=["1",['save as::-MixFigure1-']])],
          [sg.Image("./MixFigure2.png",tooltip=("The 13C DEPT NMR of the mixture"),key="-MixFigure2-",pad=((0,0),(0,0)),background_color='white',right_click_menu=["1",['save as::-MixFigure2-']])],
          [sg.Image("./MixFigure3.png",tooltip=("The hitted peaks of the mixture"),key="-MixFigure3-",pad=((0,0),(0,0)),background_color='white',right_click_menu=["1",['save as::-MixFigure3-']])],
          [sg.Image("./MixFigure4.png",tooltip=("The remain peaks after moving known compound"),key="-MixFigure4-",pad=((0,0),(0,0)),background_color='white',right_click_menu=["1",['save as::-MixFigure4-']])]
         ]
          

layout21=[
          [sg.T("1.START WITH RAW DATA:",text_color="white",font="Arial")],
          #[sg.T("Select file path and corresponding peak extraction parameters：")],
          [sg.FileBrowse(button_text="13C_data"),sg.In("",size=(30,1))],
          [sg.FolderBrowse(button_text="dir_COSY"),sg.In("",size=(30,1))],
          [sg.T("contour start COSY"),sg.In("5000000",size=(10,0),key="-startCOSY-"),sg.T("threshold COSY"),sg.In("20000000",size=(10,0),key="-thresholdCOSY-")],
          [sg.FolderBrowse(button_text="dir_HSQC"),sg.In("",size=(30,1))],
          [sg.T("contour start QC"),sg.In("5000000",size=(10,0),key="-startQC-"),sg.T("threshold QC"),sg.In("20000000",size=(10,0),key="-thresholdQC-")],
          [sg.FolderBrowse(button_text="dir_HMBC"),sg.In("",size=(30,1))],
          [sg.T("contour start HMBC"),sg.In("5000000",size=(10,0),key="-startBC-"),sg.T("threshold HMBC"),sg.In("20000000",size=(10,0),key="-thresholdBC-")],
          [sg.B("Submit",key="-submit21-")],
          [sg.T("Peak extraction results:"),sg.B("HSQC",pad=((0,0),(0,0)),key="-figQC-"),sg.B("COSY",pad=((0,0),(0,0)),key="-figCOSY-"),sg.B("HMBC",pad=((0,0),(0,0)),key="-figBC-")],
          [sg.Image("./test_2d1.png",key="-fig2D-",background_color="cyan")],
          [sg.T("2.START WITH PEAK DATA:",text_color="white",font="Arial")],
          #[sg.T("Select the file containing peak information：")],
          [sg.FileBrowse(button_text="1D_data"),sg.In("",size=(30,1))],
          [sg.FileBrowse(button_text="2D_data",target="-IN-"),sg.In("",size=(30,1),key="-IN-"),sg.B("Submit",key="-submit22")]
          ]

listQC=[]
listCOSY=[]
listBC=[]
listCOSY2C=[]
listBC2C=[]

layout22=[
          [sg.T("2D NMR DATA:",text_color="white",font="Arial")],
          [sg.T("1H 1H COSY:",pad=((0,0),(0,0))),sg.T("HMBC:",pad=((30,60),(0,0))),sg.T("HSQC:")],
          [sg.LB(listCOSY,size=(15,15),key="-LBCOSY-",background_color="lightcyan"),sg.LB(listBC,size=(15,15),key="-LBBC-",background_color="lightcyan"),sg.LB(listQC,size=(15,15),key="-LBQC-",background_color="lightcyan")],
          [sg.T("1H 1H COSY conversion: "),sg.T("HMBC conversion:"),sg.T("Poss_errors")],
          [sg.LB(listCOSY2C,size=(15,16),key="-LBCOSYtrans-",background_color="khaki"),
           sg.LB(listBC2C,size=(15,16),key="-LBBCtrans-",background_color="khaki"),
           sg.LB(listBC2C,size=(15,16),key="-LBPossErrors-",background_color="khaki")]
           ]
          
layout23=[
          [sg.Image("./test_2d2.png",key="-network-",pad=((0,0),(0,0)),background_color='white',right_click_menu=["1",['save as::-net-']])]
          ]

layout31=[
          [sg.T("STEP1:Enter the compound to be checked:")],
          [sg.FileBrowse(button_text="1D_data_check"),sg.In("",size=(30,1))],
          [sg.FileBrowse(button_text="2D_data_check"),sg.In("",size=(30,1))],
          [sg.FileBrowse(button_text="structure",target="-INcheck-"),sg.In("",size=(30,1),key="-INcheck-"),sg.B("submit",key="-submit31-")],
          [sg.Image("structure.png",key="-structure-")],
          [sg.T("STEP2:Give the label of the carbon(right list):"),sg.B("submit",key="-submit32-")],
          [sg.FileBrowse(button_text="InputFileCheck",target="-INassign2-"),sg.In("",size=(20,1),key="-INassign2-"),sg.B("submit",key="-submit33-")],
          [sg.Image("./abnormal_figure.png",key="-abnormalfigure-")]
          ]
compdata=pd.read_csv("NMR-1D.csv")
list_shift=compdata.iloc[:,0].tolist()
list_type=compdata.iloc[:,1].tolist()
layout32=[
          [sg.T(str(list_shift[i]),text_color="black",background_color="gray"),
           sg.T(list_type[i],text_color="black",background_color="gray"),
           sg.In("",size=(10,1),key="-"+"CLabel"+str(i)+"-")] for i in range(len(list_shift))                                                                  
         ]
Normal_BC=[]
Abnornal_BC=[]
Abnormal_BClabel=[]
Normal_COSY=[]
Abnornal_COSY=[]
Abnormal_COSYlabel=[]
layout33=[[sg.T("Normal BC:"),sg.T("Abnornal BC:",pad=((120,120),(0,0))),sg.T("Abnormal BC label:")],
       [sg.LB(Normal_BC,size=(25,15),key="-NormalBC-"),sg.LB(Abnornal_BC,size=(25,15),key="-AbnormalBC-"),
       sg.LB(Abnormal_BClabel,size=(25,15),key="-AbnormalBClabel-")],
       [sg.T("Normal COSY:"),sg.T("Abnornal COSY:",pad=((105,90),(0,0))),sg.T("Abnormal COSY label:")],
       [sg.LB(Normal_COSY,size=(25,15),key="-NormalCOSY-"),sg.LB(Abnornal_COSY,size=(25,15),key="-AbnormalCOSY-"),
       sg.LB(Abnormal_COSYlabel,size=(25,15),key="-AbnormalCOSYlabel-")]
        ]

layout41=[
          [sg.FileBrowse(button_text="1D_data_assigned"),sg.In("",size=(30,1))],
          [sg.FileBrowse(button_text="2D_data_assigned"),sg.In("",size=(30,1))],
          [sg.T("STEP1:Enter the compound to be assigned:")],
          [sg.FileBrowse(button_text="Structure",target="-INassign-"),sg.In("",size=(30,1),key="-INassign-"),sg.B("submit",key="-submit41-")],
          [sg.Image("structureAssign.png",key="-structureAssign-")],
          [sg.T("STEP2:Give the label of the carbon(right list):"),sg.B("submit",key="-submit42-")],
          [sg.FileBrowse(button_text="InputFile",target="-INassign1-"),sg.In("",size=(30,1),key="-INassign1-"),sg.B("submit",key="-submit43-")],
          [sg.Image("visualization_check.png",key="-visualizationcheck-",background_color="white")]
          ]
layout42=[
          [sg.T(str(list_shift[i]),text_color="black",background_color="gray"),
           sg.T(list_type[i],text_color="black",background_color="gray"),
           sg.In("",size=(10,1),key="-"+"CLabela"+str(i)+"-")] for i in range(len(list_shift))                                                                  
         ]
blank=[]
layout43=[[sg.T("CorrectP_NO./Total/BC/COSY:"),sg.T("Abnormal BC:",pad=((40,60),(0,0))),sg.T("Abnormal BC label:")],
          [sg.LB(blank,key="-CorrectP-",size=(25,20)),sg.LB(blank,key="-AbnormalBC1-",size=(25,20)),sg.LB(blank,key="-AbnormalBClabel1-",size=(10,20))],
          [sg.T("NO."),sg.In("",key="-singlesitu-",size=(15,1)),sg.B("submit",key="-submit44-")],
          [sg.T("shift&label"),sg.T("Abnormal COSY:",pad=((150,50),(0,0))),sg.T("Abnormal COSY label:",pad=((0,0),(0,0)))],
          [sg.LB(blank,key="-onesituation-",size=(25,20)),sg.LB(blank,key="-AbnormalCOSY1-",size=(25,20)),sg.LB(blank,key="-AbnormalCOSYlabel1-",size=(10,20))]
         ]
layout51=[[sg.T("INTEGRATE BY INTERVAL")],
          [sg.T("STEP1:1H NMR data input:")],
          [sg.FileBrowse(button_text="1HNMR data",target="-INcheck2-"),sg.In("",size=(30,1),key="-INcheck2-")],
          [sg.T("STEP2:Parameters forintegration:")],
          [sg.T("Integration interval:"),sg.In("0.5",size=(20,1),key="-interval-")],
          [sg.T("Minimum chemical shift:"),sg.In("0",size=(20,1),key="-minppm-")],
          [sg.T("Maximum chemical shift:"),sg.In("10",size=(20,1),key="-maxppm-")],
          [sg.B("submit",key="-submit51-")],
          [sg.T("INTEGRATE BY COMPOUND")],
          [sg.T("STEP1:1H NMR and compound data input:")],
          [sg.FileBrowse(button_text="1HNMR data",target="-INcheck3-"),sg.In("",size=(30,1),key="-INcheck3-")],
          [sg.FileBrowse(button_text="compound information",target="-INcheck4-"),sg.In("",size=(30,1),key="-INcheck4-")],
          [sg.B("submit",key="-submit52-")],
          [sg.T("STEP2:Parameters for quantitative:")],
          [sg.T("Concentration of internal standard (mol/L)"),sg.In("",key="-C0-",size=(10,1))],
          [sg.T("Proton numbers of internal standard"),sg.In("",key="-N0-",size=(10,1))],
          [sg.T("Volume of deuterated reagent(L)"),sg.In("",key="-volume-",size=(10,1))],
          [sg.T("Weight of the sample (g)"),sg.In("",key="-weight-",size=(10,1))],
          [sg.B("submit",key="-submit53-")]
         ]
layout52=[[sg.ML("abc",size=(70,20),key="-ingrateresult-")]
          
         ]
layout61=[
          [sg.T("STEP1:2D NMR data input:")],
          [sg.FileBrowse(button_text="2DNMR data",target="-INestablish1-"),sg.In("",size=(30,1),key="-INestablish1-")],
          [sg.T("STEP2:Input chemical shift/degree/thereshold:")],
          [sg.T("Chemical shift:"),sg.In("",size=(40,1),key="-EstablishShift-")],
          [sg.T("Degree:"),sg.In("",size=(40,1),key="-EstablishDegree-")],
          [sg.T("thereshold:"),sg.In("0.8",size=(20,1),key="-Establishtheres-")],
          [sg.T("ComfirmList"),sg.In("",size=(40,1),key="-ComfirmList-")],
          [sg.B("submit",key="-submit61-")],
          [sg.T("The analysis result:")],
          [sg.Image("establishfig1.png",key="-establishfig1-")],
          [sg.T("STEP3:Check fragment:"),sg.In("",size=(20,1),key="-EstablishIn1-"),sg.B("submit",key="-submit62-")],
          [sg.T("The smiles of the select fragment:")],
          [sg.In("",size=(50,1),key="-SmileFragment-")],
#           [sg.T("STEP1:1H NMR and compound data input:")],
#           [sg.FileBrowse(button_text="1HNMR data",target="-INcheck3-"),sg.In("",size=(30,1),key="-INcheck3-")],
#           [sg.FileBrowse(button_text="compound information",target="-INcheck4-"),sg.In("",size=(30,1),key="-INcheck4-")],
#           [sg.B("submit",key="-submit52-")],
#           [sg.T("STEP2:Parameters for quantitative:")],
#           [sg.T("Concentration of internal standard (mol/L)"),sg.In("",key="-C0-",size=(10,1))],
#           [sg.T("Proton numbers of internal standard"),sg.In("",key="-N0-",size=(10,1))],
#           [sg.T("Volume of deuterated reagent(L)"),sg.In("",key="-volume-",size=(10,1))],
#           [sg.T("Weight of the sample (g)"),sg.In("",key="-weight-",size=(10,1))],
#           [sg.B("submit",key="-submit53-")]
         ]
layout62=[
          [sg.T("The select fragment:")],
          [sg.Image("establishfig2.png",key="-establishfig2-")],
          [sg.T("STEP4 Input the smiles of the known fragment:")],
          [sg.In("",size=(40,1),key="-EstablishIn2-"),sg.B("submit",key="-submit63-")],
          [sg.Image("establishfig3.png",key="-establishfig3-")],
          [sg.T("STEP5 Link two fragments:")],
          [sg.T("site-A1"),sg.In("",size=(10,1),key="-siteA1-"),sg.T("site-A2"),sg.In("",size=(10,1),key="-siteA2-")],
          [sg.T("site-B1"),sg.In("",size=(10,1),key="-siteB1-"),sg.T("site-B2"),sg.In("",size=(10,1),key="-siteB2-")],
          [sg.B("submit",key="-submit64-")],
          [sg.T("The combined structure:")],
          [sg.Image("establishfig4.png",key="-establishfig4-")],
          [sg.T("smiles:"),sg.In("",size=(40,1),key="-combinedsmiles-")],
          ]
#——————————-----------整体布局----------------------------------------------------------------------
layout1=[[sg.Menu(menu_def)],
         [sg.B("1DNMR",button_color="red",pad=((0,0),(0,0)),key="-func1-"),
          sg.B("2DNMR",button_color="green",pad=((0,0),(0,0)),key="-func22-"),
          sg.B("NMRassign",button_color="blue",pad=((0,0),(0,0)),key="-func4-"),
          sg.B("NMRcheck",button_color="orange",pad=((0,0),(0,0)),key="-func3-"),
          sg.B("Establish",button_color="cyan",pad=((0,0),(0,0)),key="-func6-")
          ],
         [sg.Frame("",layout11,visible=True,size=(400,750),key="-element11-"),sg.Col(layout12,key="-element12-", element_justification="left",vertical_alignment="left",size=(900,750), vertical_scroll_only=True,background_color='white',visible=True,scrollable=True),
         sg.Frame("",layout21,size=(400,750),visible=False,key="-element21-",),sg.Frame("",layout22,visible=False,size=(400,750),key="-element22-"),sg.Col(layout23,key="-element23-",visible=False,size=(660,750),background_color='white',scrollable=True),
         sg.Frame("",layout31,size=(400,750),visible=False,key="-element31-",background_color="white"),sg.Col(layout32,key="-element32-", element_justification="left",vertical_alignment="left",size=(200,750), vertical_scroll_only=True,background_color='gray',visible=False,scrollable=True),
         sg.Frame("",layout33,size=(550,750),visible=False,key="-element33-"),
         sg.Frame("",layout41,size=(400,750),visible=False,key="-element41-",background_color="white"),
         sg.Col(layout42,key="-element42-", element_justification="left",vertical_alignment="left",size=(200,750), vertical_scroll_only=True,background_color='gray',visible=False,scrollable=True),
         sg.Frame("",layout43,size=(550,750),visible=False,key="-element43-"),
         sg.Frame("",layout61,size=(550,750),visible=False,key="-element61-"),
         sg.Col(layout62,size=(550,750),element_justification="left",vertical_alignment="left",visible=False,vertical_scroll_only=True,background_color='gray',scrollable=True,key="-element62-")
         ]]
#--------------------------------------------------------------------------------------------------------------
window=sg.Window('DATAanalyte',layout1,keep_on_top=False,element_justification='left',resizable=True)
win2_active = False
while True:
    event,values=window.read()
    print(event)
    print(values)
    if event==None:
        break
#——————————————————————窗口切换————————————————————————————————
    if event=="-func1-":
        window["-element11-"].update(visible=True)
        window["-element12-"].update(visible=True)
        window["-element21-"].update(visible=False)
        window["-element22-"].update(visible=False)
        window["-element23-"].update(visible=False)
        window["-element31-"].update(visible=False)
        window["-element32-"].update(visible=False)
        window["-element33-"].update(visible=False)
        window["-element41-"].update(visible=False)
        window["-element42-"].update(visible=False)
        window["-element43-"].update(visible=False)
        window["-element61-"].update(visible=False)
        window["-element62-"].update(visible=False)
    if event=="-func22-":
        window["-element21-"].update(visible=True)
        window["-element22-"].update(visible=True)
        window["-element23-"].update(visible=True)
        window["-element11-"].update(visible=False)
        window["-element12-"].update(visible=False)
        window["-element31-"].update(visible=False)
        window["-element32-"].update(visible=False)
        window["-element33-"].update(visible=False)
        window["-element41-"].update(visible=False)
        window["-element42-"].update(visible=False)
        window["-element43-"].update(visible=False)
        window["-element61-"].update(visible=False)
        window["-element62-"].update(visible=False)
    if event=="-func3-":
        window["-element11-"].update(visible=False)
        window["-element12-"].update(visible=False)
        window["-element21-"].update(visible=False)
        window["-element22-"].update(visible=False)
        window["-element23-"].update(visible=False)
        window["-element31-"].update(visible=True)
        window["-element32-"].update(visible=True)
        window["-element33-"].update(visible=True)
        window["-element41-"].update(visible=False)
        window["-element42-"].update(visible=False)
        window["-element43-"].update(visible=False)
        window["-element61-"].update(visible=False)
        window["-element62-"].update(visible=False)
    if event=="-func4-":
        window["-element11-"].update(visible=False)
        window["-element12-"].update(visible=False)
        window["-element21-"].update(visible=False)
        window["-element22-"].update(visible=False)
        window["-element23-"].update(visible=False)
        window["-element31-"].update(visible=False)
        window["-element32-"].update(visible=False)
        window["-element33-"].update(visible=False)
        window["-element41-"].update(visible=True)
        window["-element42-"].update(visible=True)
        window["-element43-"].update(visible=True)
        window["-element61-"].update(visible=False)
        window["-element62-"].update(visible=False)
    if event=="-func6-":
        window["-element11-"].update(visible=False)
        window["-element12-"].update(visible=False)
        window["-element21-"].update(visible=False)
        window["-element22-"].update(visible=False)
        window["-element23-"].update(visible=False)
        window["-element31-"].update(visible=False)
        window["-element32-"].update(visible=False)
        window["-element33-"].update(visible=False)
        window["-element41-"].update(visible=False)
        window["-element42-"].update(visible=False)
        window["-element43-"].update(visible=False)
        window["-element61-"].update(visible=True)
        window["-element62-"].update(visible=True)

        
#——————————————————————Page1——————————————————————————————————————
    #增加第二个窗口用于图片的属性设置
    if not win2_active and event == "properties::-fig3-":
        win2_active = True
        list1=["On","Off"]
        layout_win2=[[sg.T("X from"),sg.In("",size=(10,0),key="-xleft-"),sg.T("to"),sg.In("",size=(10,0),key="-xright-")],
             [sg.T("Y from"),sg.In("",size=(10,0),key="-yneg-"),sg.T("to"),sg.In("",size=(10,0),key="-ypos-")],
                    [sg.T("Label:"),sg.Combo(["False","True"],size=(10,1),key="-label-"),sg.T("Chemical shift:"),sg.Combo(["False","True"],size=(10,1),key="-shift-")],
                    [sg.T("Color reset:")],
                    [sg.ColorChooserButton(button_text="CH3",target="-CH3color-"),sg.In("#ff0000",size=(10,0),key="-CH3color-"),
                     sg.ColorChooserButton(button_text="CH2",target="-CH2color-"),sg.In("#008000",size=(10,0),key="-CH2color-")],
                    [sg.ColorChooserButton(button_text="CH",target="-CHcolor-"),sg.In("#0000ff",size=(10,0),key="-CHcolor-"),
                     sg.ColorChooserButton(button_text="C",target="-Ccolor-"),sg.In("#ff8000",size=(10,0),key="-Ccolor-")],
                    [sg.B("Set",key="-setproperty-")]
                    ]
        win2 = sg.Window('properties', layout_win2)
        
    if win2_active:
        event2,values2 = win2.read()
        if event2 == sg.WIN_CLOSED or event2 == 'Exit':
            win2_active  = False
            win2.close()
    if event=="-submit1-":
        try:
            data_dir=values['13C_NMR']
            data_dir135=values['DEPT135']
            data_dir90=values["DEPT90"]
            test=CNMR_process(data_dir,data_dir135,data_dir90)
            test.get_13CDEPT_peak(threshold_C=float(values['-threshold_C-']),threshold_C90=float(values['-threshold_C90-']),
                                  threshold_C135pos=float(values['-threshold_C135pos-']),threshold_C135neg=float(values['-threshold_C135neg-']),DrawWin="True")
            window['-fig1-'].update("./fig1.png")
        except Exception as e:
             sg.popup(e)
    
    if event=="-combine-":
        try:
            test.sort_Ctype_mindelta()
            test.CNMR_reconstract(width=0.6,fontsize='medium',chemical_shift="True",label="False",summary="True",x_left=200,x_right=-2,y_top=0.6*1e8,y_bottom=-0.6*1e8,DrawWin="True")
            window['-fig2-'].update("./fig2.png")
        except Exception as e:
             sg.popup(e)
    if event=="-submit2-":
        try:
            if values["-type-"]=="CH3":
                test.solvent_remove(test.hit_CH3,type="CH3", solvent=values["-solvent-"])
                test.CNMR_reconstract(width=0.6,fontsize='medium',chemical_shift="True",label="False",summary="True",x_left=200,x_right=-2,y_top=0.6*1e8,y_bottom=-0.6*1e8,DrawWin="True")
                window['-fig2-'].update("./fig2.png")
            if values["-type-"]=="CH2":
                test.solvent_remove(test.hit_CH2,type="CH2", solvent=values["-solvent-"]) 
                test.CNMR_reconstract(width=0.6,fontsize='medium',chemical_shift="True",label="False",summary="True",x_left=200,x_right=-2,y_top=0.6*1e8,y_bottom=-0.6*1e8,DrawWin="True")
                window['-fig2-'].update("./fig2.png")
            if values["-type-"]=="CH":
                test.solvent_remove(test.hit_CH,type="CH", solvent=values["-solvent-"])
                test.CNMR_reconstract(width=0.6,fontsize='medium',chemical_shift="True",label="False",summary="True",x_left=200,x_right=-2,y_top=0.6*1e8,y_bottom=-0.6*1e8,DrawWin="True")
                window['-fig2-'].update("./fig2.png")
            if values["-type-"]=="C":
                test.solvent_remove(test.Hit_C_or_unhited,type="C", solvent=values["-solvent-"])
                test.CNMR_reconstract(width=0.6,fontsize='medium',chemical_shift="True",label="False",summary="True",x_left=200,x_right=-2,y_top=0.6*1e8,y_bottom=-0.6*1e8,DrawWin="True")
                window['-fig2-'].update("./fig2.png")
        except Exception as e:
             sg.popup(e)
    if event=="-submit3-":
        try:
            if values["-typeimpurity-"]=="CH3":
                test.impurity_removal(test.hit_CH3,type="CH3",rate=0.5,Standsrd_num=0,method="absolute",threshold=float(values["-thresholdimpurity1-"])*float(values["-thresholdimpurity2-"]))
                test.CNMR_reconstract(width=0.6,fontsize='medium',chemical_shift="True",label="False",summary="True",x_left=200,x_right=-2,y_top=0.6*1e8,y_bottom=-0.6*1e8)
                window['-fig2-'].update("./fig2.png")
            if values["-typeimpurity-"]=="CH2":
                test.impurity_removal(test.hit_CH2,type="CH2",rate=0.5,Standsrd_num=0,method="absolute",threshold=float(values["-thresholdimpurity1-"])*float(values["-thresholdimpurity2-"]))
                test.CNMR_reconstract(width=0.6,fontsize='medium',chemical_shift="True",label="False",summary="True",x_left=200,x_right=-2,y_top=0.6*1e8,y_bottom=-0.6*1e8)
                window['-fig2-'].update("./fig2.png")
            if values["-typeimpurity-"]=="CH":
                test.impurity_removal(test.hit_CH,type="CH",rate=0.5,Standsrd_num=0,method="absolute",threshold=float(values["-thresholdimpurity1-"])*float(values["-thresholdimpurity2-"]))
                test.CNMR_reconstract(width=0.6,fontsize='medium',chemical_shift="True",label="False",summary="True",x_left=200,x_right=-2,y_top=0.6*1e8,y_bottom=-0.6*1e8)
                window['-fig2-'].update("./fig2.png")
            if values["-typeimpurity-"]=="C":
                test.impurity_removal(test.Hit_C_or_unhited,type="C",rate=0.5,Standsrd_num=0,method="absolute",threshold=float(values["-thresholdimpurity1-"])*float(values["-thresholdimpurity2-"]))
                test.CNMR_reconstract(width=0.6,fontsize='medium',chemical_shift="True",label="False",summary="True",x_left=200,x_right=-2,y_top=0.6*1e8,y_bottom=-0.6*1e8)
                window['-fig2-'].update("./fig2.png")
        except Exception as e:
             sg.popup(e)
    if event=="-CNMRcompare-":
        try:
            C_DEPT_data=pd.read_excel("Inhouse_database_1.xls",sheet_name="C_DEPT_NMR")
            test.combine_data(CarbonFileName="NMR-1D.csv")
            test.carbon_compare_new(C_DEPT_data,CompareResult="CompareResult.csv",therohold=0.8)
            window['-fig3-'].update("./fig3.png")
        except Exception as e:
             sg.popup(e)
    if event=="-CDEPTNMRcompare-":
        try:
            C_DEPT_data=pd.read_excel("Inhouse_database_1.xls",sheet_name="C_DEPT_NMR")
            test.combine_data(CarbonFileName="NMR-1D.csv")
            test.carbon_compare_new(C_DEPT_data,CompareResult="CompareResult.csv",therohold=0.8)
            window['-fig3-'].update("./fig3.png")
        except Exception as e:
             sg.popup(e)
    if event=="-submit5-":
        try:
            database=C_DEPT_data
            test.compare_select(int(values["-compname-"]),database)
            window['-fig4-'].update("./fig4.png")
        except Exception as e:
             sg.popup(e)
    if event=="-submit6-":
        try:
            CarbonFileName="NMR-1D.csv"
            CompareResult="CompareResult.csv"
            mixtask=mixture_analyte(CarbonFileName,CompareResult)
            num=int(values["-compnum-"])
            mixtask.mix_analyte(num)
            window['-MixFigure1-'].update("./MixFigure1.png")
            window['-MixFigure2-'].update("./MixFigure2.png")
            window['-MixFigure3-'].update("./MixFigure3.png")
        except Exception as e:
            sg.popup(e)
    if event=="-submit7-":
        try:
            num=[int(values["-compnum-"])]
            mixtask.move_known(num)
            window['-MixFigure4-'].update("./MixFigure4.png")
        except Exception as e:
            sg.popup(e)
#______________________________________page2______________________________________________
    if event=="-submit21-":
        try:
            CarbonListFile=values["13C_data"]
            dir_QC=values["dir_HSQC"]
            dir_COSY=values["dir_COSY"]
            dir_BC=values["dir_HMBC"]
            contour_start_QC=int(values["-startQC-"])
            contour_start_COSY=int(values["-startCOSY-"])
            contour_start_BC=int(values["-startBC-"])
            threshold_QC=int(values["-thresholdQC-"])
            threshold_COSY=int(values["-thresholdCOSY-"])
            threshold_BC=int(values["-thresholdBC-"])
           
            test2D=analysis_2DNMR(CarbonListFile)
            test2D.autoanalysisGUI(dir_QC,dir_COSY,dir_BC,contour_start_QC,contour_start_COSY,contour_start_BC,threshold_QC,threshold_COSY,threshold_BC)
            print("123")
            window["-fig2D-"].update("./QC.png")
            #重新读取相关数据，并把数据放在二维相关的列表中
            cor_data=pd.read_csv("cor_data.csv")
            QC_H=cor_data["QC_H"].tolist()
            QC_C=cor_data["QC_C"].tolist()
            COSY_H1=cor_data["COSY_H1"].tolist()
            COSY_H2=cor_data["COSY_H2"].tolist()
            BC_H=cor_data["BC_H"].tolist()
            BC_C=cor_data["BC_C"].tolist()
            COSY_C1=cor_data["COSY_C1"].tolist()
            COSY_C2=cor_data["COSY_C2"].tolist()
            BC_H2C=cor_data["BC_H2C"].tolist()
            BC_C2C=cor_data["BC_C2C"].tolist()
            #去除Nan
            QC_H= [a_ for a_ in QC_H if a_ == a_]
            QC_C=[a_ for a_ in QC_C if a_ == a_]
            COSY_H1= [a_ for a_ in COSY_H1 if a_ == a_]
            COSY_H2=[a_ for a_ in COSY_H2 if a_ == a_]
            BC_H= [a_ for a_ in BC_H if a_ == a_]
            BC_C=[a_ for a_ in BC_C if a_ == a_]
            COSY_C1= [a_ for a_ in COSY_C1 if a_ == a_]
            COSY_C2=[a_ for a_ in COSY_C2 if a_ == a_]
            BC_H2C= [a_ for a_ in BC_H2C if a_ == a_]
            BC_C2C=[a_ for a_ in BC_C2C if a_ == a_]
            #合并相关数据
            QC=combine_cor(QC_H,QC_C)
            COSY=combine_cor(COSY_H1,COSY_H2)
            BC=combine_cor(BC_H,BC_C)
            COSY_transfer=combine_cor(COSY_C1,COSY_C2)
            BC_transfer=combine_cor(BC_H2C,BC_C2C)
            #更新列表
            window["-LBCOSY-"].update(COSY)
            window["-LBBC-"].update(BC)
            window["-LBQC-"].update(QC)
            window["-LBCOSYtrans-"].update(COSY_transfer)
            window["-LBBCtrans-"].update(BC_transfer)
            #更新网络图
            window["-network-"].update("netwrok.png")
        except Exception as e:
             sg.popup(e)
    if event=="-submit22":
        try:
            CarbonListFile=values["1D_data"]
            peakdatadir=values["2D_data"]
            test2D=analysis_2DNMR(CarbonListFile)
            test2D.semi_autoanalysis(peakdatadir)
            #重新读取相关数据，并把数据放在二维相关的列表中
            cor_data=pd.read_csv("cor_data.csv")
            QC_H=cor_data["QC_H"].tolist()
            QC_C=cor_data["QC_C"].tolist()
            COSY_H1=cor_data["COSY_H1"].tolist()
            COSY_H2=cor_data["COSY_H2"].tolist()
            BC_H=cor_data["BC_H"].tolist()
            BC_C=cor_data["BC_C"].tolist()
            COSY_C1=cor_data["COSY_C1"].tolist()
            COSY_C2=cor_data["COSY_C2"].tolist()
            BC_H2C=cor_data["BC_H2C"].tolist()
            BC_C2C=cor_data["BC_C2C"].tolist()
            #去除Nan
            QC_H= [a_ for a_ in QC_H if a_ == a_]
            QC_C=[a_ for a_ in QC_C if a_ == a_]
            COSY_H1= [a_ for a_ in COSY_H1 if a_ == a_]
            COSY_H2=[a_ for a_ in COSY_H2 if a_ == a_]
            BC_H= [a_ for a_ in BC_H if a_ == a_]
            BC_C=[a_ for a_ in BC_C if a_ == a_]
            COSY_C1= [a_ for a_ in COSY_C1 if a_ == a_]
            COSY_C2=[a_ for a_ in COSY_C2 if a_ == a_]
            BC_H2C= [a_ for a_ in BC_H2C if a_ == a_]
            BC_C2C=[a_ for a_ in BC_C2C if a_ == a_]
            #合并相关数据
            QC=combine_cor(QC_H,QC_C)
            COSY=combine_cor(COSY_H1,COSY_H2)
            BC=combine_cor(BC_H,BC_C)
            COSY_transfer=combine_cor(COSY_C1,COSY_C2)
            BC_transfer=combine_cor(BC_H2C,BC_C2C)
            print(test2D.poss_error1)
            print(test2D.poss_error2)
            poss_errors=test2D.poss_error1+test2D.poss_error2
            #更新列表
            window["-LBCOSY-"].update(COSY)
            window["-LBBC-"].update(BC)
            window["-LBQC-"].update(QC)
            window["-LBCOSYtrans-"].update(COSY_transfer)
            window["-LBBCtrans-"].update(BC_transfer)
            window["-LBPossErrors-"].update(poss_errors)
            #更新网络图
            window["-network-"].update("netwrok.png")
            
        except Exception as e:
             sg.popup(e)
    if event=="-figQC-":
        window["-fig2D-"].update("./QC.png")
    if event=="-figCOSY-":
        window["-fig2D-"].update("./COSY.png")
    if event=="-figBC-":
        window["-fig2D-"].update("./BC.png")

    #保存图片
    if event=="save as::-fig1-":
        #result=sg.PopupGetFolder("Please select the path:")
        try:
            path = sg.popup_get_file('Please select the path:',save_as=True,default_extension='png')
            figure= Image.open("fig1.png") 
            figure.save(path)
        except Exception as e:
            sg.popup(e)
    if event=="save as::-fig2-":
        try:
            path = sg.popup_get_file('Please select the path:',save_as=True,default_extension='png')
            figure=Image.open("fig2.png") 
            figure.save(path)
        except Exception as e:
            sg.popup(e)
    if event=="save as::-fig3-":
        try:
            path = sg.popup_get_file('Please select the path:',save_as=True,default_extension='png')
            figure= Image.open("fig3.png") 
            figure.save(path)
        except Exception as e:
            sg.popup(e)
    if event=="save as::-fig4-":
        try:
            path = sg.popup_get_file('Please select the path:',save_as=True,default_extension='png')
            figure= Image.open("fig4.png") 
            figure.save(path)
        except Exception as e:
            sg.popup(e)
    if event=="save as::-net-":
        try:
            path = sg.popup_get_file('Please select the path:',save_as=True,default_extension='png')
            figure= Image.open("network.png") 
            figure.save(path)
        except Exception as e:
            sg.popup(e)
#     if event=="properties::-fig3-":
#         try:
#             window2=sg.Window('properties',layout_win2,keep_on_top=False,element_justification='left',resizable=True)
#             while True:
#                 event1,values1=window.read()
#                 print(event1)
#                 print(values1)
#             window2.close()
#         except Exception as e:
#             sg.popup(e)
#-----------------------------------Page3-----------------------------------------------------
    if event=="-submit31-":
        try:
            comp_dir=values["structure"]
            CarbonListFile=values["1D_data_check"]
            CorrelationFile=values["2D_data_check"]
            check_task=NMRcheck(comp_dir,CarbonListFile,CorrelationFile)
            check_task.get_comp_figure()
            check_task.fig.save('structure.png')
            window["-structure-"].update('structure.png')
        except Exception as e:
            sg.popup(e)
    if event=="-submit32-":
        try:
            Carbon_label=[]
            for i in range(len(list_shift)):
                Carbon_label.append(int(values["-CLabel"+str(i)+"-"]))
            a=[list_shift, Carbon_label]
            df=pd.DataFrame(a)
            df=pd.DataFrame(df.values.T)
            df.to_csv("InputFile.csv",index=None)
            check_task.BC_check(Carbon_label,0.5)
            check_task.COSY_check(Carbon_label,0.5)
            print(check_task.Normal_COSY)
            check_task.abnormal_figure_check()
            window["-NormalBC-"].update(check_task.Normal_BC)
            window["-AbnormalBC-"].update(check_task.Abnormal_BC)
            window["-AbnormalBClabel-"].update(check_task.Abnormal_BClabel)
            window["-NormalCOSY-"].update(check_task.Normal_COSY)
            window["-AbnormalCOSY-"].update(check_task.Abnormal_COSY)
            window["-AbnormalCOSYlabel-"].update(check_task.Abnormal_COSYlabel)
            window["-abnormalfigure-"].update('abnormal_figure.png')
        except Exception as e:
            sg.popup(e)
    if event=="-submit33-":
        try:
            InputFilePath1=values["InputFileCheck"]
            df=pd.read_csv(InputFilePath1)
            Carbon_label=df.iloc[:,1].tolist()
            check_task.BC_check(Carbon_label,0.5)
            check_task.COSY_check(Carbon_label,0.5)
            print(check_task.Normal_COSY)
            check_task.abnormal_figure_check()
            window["-NormalBC-"].update(check_task.Normal_BC)
            window["-AbnormalBC-"].update(check_task.Abnormal_BC)
            window["-AbnormalBClabel-"].update(check_task.Abnormal_BClabel)
            window["-NormalCOSY-"].update(check_task.Normal_COSY)
            window["-AbnormalCOSY-"].update(check_task.Abnormal_COSY)
            window["-AbnormalCOSYlabel-"].update(check_task.Abnormal_COSYlabel)
            window["-abnormalfigure-"].update('abnormal_figure.png')
        except Exception as e:
            sg.popup(e)
#-----------------------------------Page4------------------------------------------------------
    if event=="-submit41-":
        try:
            comp_dir=values["Structure"]
            CarbonListFile=values["1D_data_assigned"]
            CorrelationFile=values["2D_data_assigned"]
            assign_task=NMRcheck(comp_dir,CarbonListFile,CorrelationFile)
            assign_task.get_comp_figure()
            assign_task.fig.save('structureAssign.png')
            window["-structureAssign-"].update('structureAssign.png')
        except Exception as e:
            sg.popup(e)
    if event=="-submit42-":
        try:
            #获取Carbon_label
            Carbon_label=[]
            for i in range(len(list_shift)):
                if values["-CLabela"+str(i)+"-"]=="":
                    Carbon_label.append("nan")
                else:
                    Carbon_label.append(int(values["-CLabela"+str(i)+"-"]))
            #穷举每种情况，计算
            a=[list_shift, Carbon_label]
            df=pd.DataFrame(a)
            df=pd.DataFrame(df.values.T)
            df.to_csv("InputFile.csv",index=None)
            InputFilePath="InputFile.csv"
            assign_task.NMRassign(InputFilePath,0.5)
            #可视化
            assign_task.visualization_check()
            #把每种情况的CorrectP的汇总结果放在同一个表中
            text_LB1=[]
            for i in range(len(assign_task.CorrectP_Total_list)):
                per=[i,assign_task.CorrectP_Total_list[i],assign_task.CorrectP_BC_list[i], assign_task.CorrectP_COSY_list[i]]
                text_LB1.append(per)
                per=[]
            #更新列表和TOP10的柱状图
            window["-visualizationcheck-"].update('visualization_check.png')
            window["-CorrectP-"].update(text_LB1)  
        except Exception as e:
            sg.popup(e)
            
    if event=="-submit43-":
        try:  
            InputFilePath=values["InputFile"]
            assign_task.NMRassign(InputFilePath,0.5)
            #可视化
            assign_task.visualization_check()
            #把每种情况的CorrectP的汇总结果放在同一个表中
            text_LB1=[]
            for i in range(len(assign_task.CorrectP_Total_list)):
                per=[i,assign_task.CorrectP_Total_list[i],assign_task.CorrectP_BC_list[i], assign_task.CorrectP_COSY_list[i]]
                text_LB1.append(per)
                per=[]
            #更新列表和TOP10的柱状图
            window["-visualizationcheck-"].update('visualization_check.png')
            window["-CorrectP-"].update(text_LB1)  
        except Exception as e:
            sg.popup(e)
        
    if event=="-submit44-":

        try:
            label=int(values["-singlesitu-"])
            window["-AbnormalBC1-"].update(assign_task.Abnormal_BC_list[label])
            window["-AbnormalCOSY1-"].update(assign_task.Abnormal_COSY_list[label])
            window["-AbnormalBClabel1-"].update(assign_task.Abnormal_BClabel_list[label])
            window["-AbnormalCOSYlabel1-"].update(assign_task.Abnormal_COSYlabel_list[label])
            #输出该情况的排序
            df_shift=pd.read_csv(CarbonListFile)
            shift_list=df_shift.iloc[:,0].tolist()
            type_list=df_shift.iloc[:,1].tolist()
            Carbon_label_for_one=assign_task.Carbon_label_all[label]
            shift_label=[]
            for i in range(len(shift_list)):
                a=[shift_list[i],type_list[i],Carbon_label_for_one[i]]
                shift_label.append(a)
                a=[]
            window["-onesituation-"].update(shift_label)
        except Exception as e:
            sg.popup(e)
            
    if event=="-submit61-":

        try:
            strnodename=values["-EstablishShift-"]
            strnodenum=values["-EstablishDegree-"]
            NMR2DFlie=values["-INestablish1-"]
            thereshold=float(values["-Establishtheres-"])
            a=values["-ComfirmList-"]
            b=a.split(';')
            ComfirmList=[]
            if a!="":
                for i,link in enumerate(b):
                    c=link.split(',')
                    d=(float(c[0]),float(c[1]))
                    ComfirmList.append(d)
            task=StructureEstablish(strnodename,strnodenum,NMR2DFlie,ComfirmList)
            task.AllStructure(thereshold)
            task.ResultDisplay(10)
            window["-establishfig1-"].update('establishfig1.png')
            
            
        except Exception as e:
            sg.popup(e)
            
    if event=="-submit62-":
        try:
            num=int(values["-EstablishIn1-"])
            task.ShowFragment(num)
            window["-SmileFragment-"].update(task.SmileRemain[num])
            window["-establishfig2-"].update("molecule1.png")
            
        except Exception as e:
            sg.popup(e)
    if event=="-submit63-":
        try:
            SmileKnown=values["-EstablishIn2-"]
            mol_b= Chem.MolFromSmiles(SmileKnown)
            opts = DrawingOptions()
            opts.atomLabelFontSize=500
            opts.includeAtomNumbers=True
            opts.bondLineWidth=2
            opts.dblBondLengthFrac=0.8
            opts.dotsPerAngstrom=1000
            figb=Draw.MolToImage(mol_b,options=opts,size=(300,300))
            figb.save("KnownStructure.png")
            window["-establishfig3-"].update("KnownStructure.png")
            
        except Exception as e:
            sg.popup(e)
    if event=="-submit64-":
        try:
            siteA1=int(values["-siteA1-"])
            siteA2=int(values["-siteA2-"])
            siteB1=int(values["-siteB1-"])
            siteB2=int(values["-siteB2-"])
            num2=int(values["-EstablishIn1-"])
            SiteNew=(siteA1,siteA2)
            SiteKnown=(siteB1,siteB2)
            task.CombineFrag(num2,mol_b,SiteNew=SiteNew,SiteKnown=SiteKnown,AtomNumbers=True,
                   LabelFontSize=500,bondLineWidth=2)
            task.fig1.save("Combined.png")
            window["-establishfig4-"].update("Combined.png")
            window["-combinedsmiles-"].update(task.FinalSmiles)
        except Exception as e:
            sg.popup(e)


# In[ ]:




