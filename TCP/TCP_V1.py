# Two-color pyrometer
"""
Created on 11 July 2020 by Chih-Ting Chen
Allow user to choose between the KL method and the ratio method
Step0: Prompt information about the input
Step1: Calibration chart
Step2: Apparant Temperature 
Step3: True Temprature
Step3: Abel Inversion (if required)
"""
import matplotlib.pyplot as plt
import tkinter as tk
import numpy as np
import png
import math
from numpy.polynomial import polynomial as P
## Step: Input information
def save_entry_fields():
    print("Dark current file: %s\n Max. Temp.: %s\n Min. Temp.: %s\n Temp. interval: %s\n Image format: %s\n bit depth: %s\n Flame file: %s\n Method: %s\n Split from center: %s" % (e1.get(),e2.get(),e3.get(),e4.get(),e5.get(),e6.get(),e7.get(),e8.get(),e9.get()))
    global InputInfo
    InputInfo = [e1.get(),e2.get(),e3.get(),e4.get(),e5.get(),e6.get(),e7.get(),e8.get(),e9.get(),e10.get(),e11.get()]
master = tk.Tk()
tk.Label(master, text="File name of the dark current image:").grid(row=0)
tk.Label(master, text="Maximum temperature of calibration:").grid(row=1)
tk.Label(master, text="Minimum temperature of calibration:").grid(row=2)
tk.Label(master, text="Temperature interval:").grid(row=3)
tk.Label(master, text="Image format:").grid(row=4)
tk.Label(master, text="Image bit depth:").grid(row=5)
tk.Label(master, text="File name of the flame image:").grid(row=6)
tk.Label(master, text="Ratio or KL method:").grid(row=7)
tk.Label(master, text="Image split along the center? Yes/No").grid(row=8)
tk.Label(master, text="Filter 1 wavelength in nm:").grid(row=9)
tk.Label(master, text="Filter 2 wavelength in nm:").grid(row=10)

e1=tk.Entry(master)
e2=tk.Entry(master)
e3=tk.Entry(master)
e4=tk.Entry(master)
e5=tk.Entry(master)
e6=tk.Entry(master)
e7=tk.Entry(master)
e8=tk.Entry(master)
e9=tk.Entry(master)
e10=tk.Entry(master)
e11=tk.Entry(master)

e1.insert(10,"darkCurrent")
e2.insert(10,"1000")
e3.insert(10,"840")
e4.insert(10,"10")
e5.insert(10,".png")
e6.insert(10,"16")
e7.insert(10,"f1")
e8.insert(10,"Ratio")
e9.insert(10,"Yes")
e10.insert(10,"690")
e11.insert(10,"750")

e1.grid(row=0,column=1)
e2.grid(row=1,column=1)
e3.grid(row=2,column=1)
e4.grid(row=3,column=1)
e5.grid(row=4,column=1)
e6.grid(row=5,column=1)
e7.grid(row=6,column=1)
e8.grid(row=7,column=1)
e9.grid(row=8,column=1)
e10.grid(row=9,column=1)
e11.grid(row=10,column=1)
tk.Button(master,text='Quit',command=master.quit).grid(row=11,column=0,sticky=tk.W,pady=4)
tk.Button(master,text='Save & Show',command=save_entry_fields).grid(row=11,column=1,sticky=tk.W,pady=4)
master.mainloop()
master.destroy()

## Step: load calibration images and subtract the dark current
NumImg =int((int(InputInfo[1])-int(InputInfo[2]))/int(InputInfo[3])+1)
temp = []
in_temp = []
caliImg =[]
for x in range(int(InputInfo[1]),int(InputInfo[2])-int(InputInfo[3]),-int(InputInfo[3])):
    temp.append(x+273)
    in_temp.append(1/float(x+273))
    CalibFileName = str(x)+InputInfo[4]
    GrayScl = 2**int(InputInfo[5])-1
    reader1 = png.Reader(filename='images/%s' %(CalibFileName))
    reader2 = png.Reader(filename='images/%s' %(InputInfo[0]+InputInfo[4]))
    pngdata1 = reader1.read()
    pngdata2 = reader2.read()
    px_array1 = np.vstack(map(np.uint16,pngdata1[2]))
    px_array2 = np.vstack(map(np.uint16, pngdata2[2]))
    px_arrayR = np.subtract(px_array1/GrayScl,px_array2/GrayScl)
    caliImg.append(px_arrayR) #The first image has the brightest field
## Step: Select analysis region for calibration
plt.figure()
CImP1 = plt.imshow(caliImg[0],cmap='summer')
plt.title('T1')
plt.colorbar(CImP1)
# select reference points for two image fields
def tellme(s):
    plt.title(s,fontsize=12,weight='bold')
    print(s)
CREF = []
CImb = []
tellme('Select one point in each field as reference.')
CREF = np.asarray(plt.ginput(2, timeout =-1))
tellme('Select three points to define image bounds.')
CImb = np.asarray(plt.ginput(3, timeout =-1))
# draw the defined rectangular region
plt.plot([int(CImb[0,0]),int(CImb[1,0])],[int(CImb[0,1]),int(CImb[0,1])],linestyle='--',linewidth=1.5,color='w')
plt.plot([int(CImb[1,0]),int(CImb[1,0])],[int(CImb[0,1]),int(CImb[2,1])],linestyle='--',linewidth=1.5,color='w')
plt.plot([int(CImb[1,0]),int(CImb[0,0])],[int(CImb[2,1]),int(CImb[2,1])],linestyle='--',linewidth=1.5,color='w')
plt.plot([int(CImb[0,0]),int(CImb[0,0])],[int(CImb[0,1]),int(CImb[2,1])],linestyle='--',linewidth=1.5,color='w')
# Calculate and drawthe analysis region bounds for the right image field
[Cdelx,Cdely]=int(CREF[0,0])-int(CREF[1,0]),int(CREF[0,1])-int(CREF[1,1])
CImb_r = np.asarray([[int(CImb[0,0])-Cdelx,int(CImb[0,1])-Cdely],[int(CImb[1,0])-Cdelx,int(CImb[1,1])-Cdely],[int(CImb[2,0])-Cdelx,int(CImb[2,1])-Cdely]])
plt.plot([int(CImb_r[0,0]),int(CImb_r[1,0])],[int(CImb_r[0,1]),int(CImb_r[0,1])],linestyle='--',linewidth=1.5,color='w')
plt.plot([int(CImb_r[1,0]),int(CImb_r[1,0])],[int(CImb_r[0,1]),int(CImb_r[2,1])],linestyle='--',linewidth=1.5,color='w')
plt.plot([int(CImb_r[1,0]),int(CImb_r[0,0])],[int(CImb_r[2,1]),int(CImb_r[2,1])],linestyle='--',linewidth=1.5,color='w')
plt.plot([int(CImb_r[0,0]),int(CImb_r[0,0])],[int(CImb_r[0,1]),int(CImb_r[2,1])],linestyle='--',linewidth=1.5,color='w')
plt.plot()
plt.show()
# Calculate the number of pixels are included
[NpixX,NpixY] = abs(int(CImb[0,0])-int(CImb[1,0])),abs(int(CImb[1,1])-int(CImb[2,1]))
NumPixs = NpixX*NpixY
# Extract data and average from the specified regions
lnAvgCali1 = []   # filter 1 (left field)
lnAvgCali2 = []   # filter 2 (right field)
for filNum in range(0,2,1):
    for T in range(0,NumImg,1):
        AvgN = 0    # initiate the index variable
        for x in range(0,NpixX,1):
            for y in range(0,NpixY,1):
                if filNum == 0:
                    AvgN = float(AvgN)+float(caliImg[T][y+int(CImb[0,1]),x+int(CImb[0,0])])/NumPixs
                else:
                    AvgN = float(AvgN)+float(caliImg[T][y+int(CImb_r[0,1]),x+int(CImb_r[0,0])])/NumPixs
        if filNum == 0:
            lnAvgCali1.append(math.log(AvgN))
        else:
            lnAvgCali2.append(math.log(AvgN))
## Step: ln(I) chart creation
plt.figure(figsize=[6.6, 4.8])
plt.plot(in_temp,lnAvgCali1,'r^',label='%s nm filter' % InputInfo[9])
plt.plot(in_temp,lnAvgCali2,'go',label='%s nm filter' % InputInfo[10])
# linear fit
c1, stats1 = P.polyfit(in_temp,lnAvgCali1,1,full=True)
c2, stats2 = P.polyfit(in_temp,lnAvgCali2,1,full=True)
plt.plot([in_temp[0],in_temp[NumImg-1]],[float(in_temp[0])*c1[1]+c1[0],float(in_temp[NumImg-1])*c1[1]+c1[0]],label='Linear fit of %s nm' % InputInfo[9])
plt.plot([in_temp[0],in_temp[NumImg-1]],[float(in_temp[0])*c2[1]+c2[0],float(in_temp[NumImg-1])*c2[1]+c2[0]],label='Linear fit of %s nm' % InputInfo[10])
# calculate the equivalent wavelength
c2_planck = (1.98644*10**-25)/(1.3806*10**-23)  # c2=hc/k
Lambda_eqv1 = (-c2_planck/c1[1])*10**9  # in nm
Lambda_eqv2 = (-c2_planck/c2[1])*10**9  # in nm
plt.text(0.00079,-3.7,'Equivalent \u03BB for %s filter = %d nm' %(InputInfo[9], Lambda_eqv1),fontsize=12)
plt.text(0.00079,-3.9,'Equivalent \u03BB for %s filter = %d nm' %(InputInfo[10], Lambda_eqv2),fontsize=12)
plt.legend(loc="upper right")
plt.title('Instrument Calibration')
plt.xlabel('1/T, 1/K')
plt.ylabel('log(intensity), a.u.')
plt.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
plt.show()
