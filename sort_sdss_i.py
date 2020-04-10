# !/usr/bin/env python3
import shutil
import csv
import pandas as pd
import fnmatch
import numpy as np
import matplotlib.pyplot as plt
from   astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from  astropy.io import fits
import os
import glob
import random
import gc

#open catalogue file for reading

fp=open('catalogue.csv','r')
reader=csv.reader(fp)
rows=list(reader)
num_lines: int=len(rows)

#print("This is the catalogue file of Ultra Diffuse Galaxies:")
#print("The number of lines in the file are:", num_lines)
#print("The file table format is:")
#print (rows[0])
#myrow=rows[1]
#myrow1=rows[num_lines-1]
#print("The first data set line is :")
#print(myrow)
#print("The last data set line is :")
#print(myrow1)
#print("This data file has data from SEXID:", myrow[1], "to SEXID:", myrow1[1] ,"\n")
#curdir=os.getcwd()
#print(curdir)

#make directory : where to put the sorted lists
os.makedirs('Sorted_lists_SDSS_i', exist_ok=True)


#file_header="KIDSFIELD,SEXID,RA,DEC,MSRB,Reff,MAG_R,SERSI"
#WARNING: do not open catalogue.csv in any other program/window, it alters the file stream.

#df.round({'Column_name':decimal places, etc})

#Surface Brightness Sorting
def zero():
    df = pd.read_csv('catalogue.csv', sep=',', skipinitialspace=True)
    #column_order = ['KIDSFIELD', 'SEXID', 'RA', 'DEC', 'MSRB', 'Reff', 'MAG_R', 'SERSI']
    df=df[['KIDSFIELD', 'SEXID', 'RA', 'DEC', 'MSRB', 'Reff', 'MAG_R', 'SERSI']] #Load by column
    df.sort_values('MSRB',inplace=True)#InPlace for writing the updated data frame. IMPORTANT!!!
    os.makedirs('Sorted_lists_SDSS_i/msrb',exist_ok=True)
    df.round({'RA': 7, 'DEC': 7, 'MSRB': 7, 'Reff': 7, 'MAG_R': 7, 'SERSI': 7}).to_csv('Sorted_lists_SDSS_i/msrb/sorted_msrb_catg.csv',index = False)

'''
def one():
    df = pd.read_csv('catalogue.csv', sep=',', skipinitialspace=True )
    column_order = ['KIDSFIELD', 'SEXID', 'RA', 'DEC', 'MSRB', 'Reff', 'MAG_R', 'SERSI']
    df = df[['KIDSFIELD', 'SEXID', 'RA', 'DEC', 'MSRB', 'Reff', 'MAG_R', 'SERSI']]
    df.sort_values('Reff', inplace=True)
    os.makedirs('Sorted_lists_SDSS_i/reff',exist_ok=True)
    df[column_order].to_csv('Sorted_lists_SDSS_i/reff/sorted_reff_catg.csv',index = False)

def two():
    df = pd.read_csv('catalogue.csv', sep=',', skipinitialspace=True )
    column_order = ['KIDSFIELD', 'SEXID', 'RA', 'DEC', 'MSRB', 'Reff', 'MAG_R', 'SERSI']
    df = df[['KIDSFIELD', 'SEXID', 'RA', 'DEC', 'MSRB', 'Reff', 'MAG_R', 'SERSI']]
    df.sort_values('MAG_R', inplace=True)
    os.makedirs('Sorted_lists_SDSS_i/magtot',exist_ok=True)
    df[column_order].to_csv('Sorted_lists_SDSS_i/magtot/sorted_magtot_catg.csv',index = False)
'''
#Sersic Index sorting
def three():
    df = pd.read_csv('catalogue.csv', sep=',', skipinitialspace=True)
    #column_order = ['KIDSFIELD', 'SEXID', 'RA', 'DEC', 'MSRB', 'Reff', 'MAG_R', 'SERSI']
    df = df[['KIDSFIELD', 'SEXID', 'RA', 'DEC', 'MSRB', 'Reff', 'MAG_R', 'SERSI']]#Load by column
    df.sort_values('SERSI', inplace=True) #InPlace for writing the updated data frame. IMPORTANT!!!
    os.makedirs('Sorted_lists_SDSS_i/seri',exist_ok=True)
    df.round({'RA': 7, 'DEC': 7, 'MSRB': 7, 'Reff': 7, 'MAG_R': 7, 'SERSI': 7}).to_csv('Sorted_lists_SDSS_i/seri/sorted_seri_catg.csv',index = False)

loopcycle=False
os.chdir('/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/')
file_list = [f for f in os.listdir('.') if f.endswith('.fits')]
os.chdir('/home/soumyananda/Dropbox/Goswami_REVISED/CATALOGRADEC_SPLICED/')

while loopcycle==False:
        print("PROGRAM STARTING \n")
        inpval = 0 #int(input("The number you have entered is: \n "))
        if inpval==0:
            #execute sort by criteria 0
            zero()
            loopcycle=True
            os.chdir('/home/soumyananda/Dropbox/Goswami_REVISED/CATALOGRADEC_SPLICED/Sorted_lists_SDSS_i/msrb/')
            f2= open('sorted_msrb_catg.csv', 'r')
            f2r = pd.read_csv(f2, sep=',', skipinitialspace=True)
            num_start = f2r.iloc[1][4]
            #print(num_start)
            t= num_lines-2
            num_end = f2r.iloc[t][4]
            #print(num_end)
            val = (num_end - num_start) / 12
            #print(val)
            os.makedirs('/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/STACKED_IMAGES_MSRB/',exist_ok=True)
            os.chdir("/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/STACKED_IMAGES_MSRB/")
            f=open("countsinfosdssi.txt","w+")
            f.write("x,y,yerr,num_files\n")
            os.chdir('/home/soumyananda/Dropbox/Goswami_REVISED/CATALOGRADEC_SPLICED/Sorted_lists_SDSS_i/msrb/')
            for i in range(1,10):
                resl=num_start + val*(i-1) #incremental dynamic setting for splitting files
                resu = num_start + val * i
                #print(resl,resu)
                #code for binned files for a range of MSRB/SER etc
                os.chdir('/home/soumyananda/Dropbox/Goswami_REVISED/CATALOGRADEC_SPLICED/Sorted_lists_SDSS_i/msrb/')
                filename = "msrb_part_" + str(i) + ".csv"
                mrc = open(filename, 'w+')
                mrc.write('KIDSFIELD, SEXID, RA, DEC,MSRB, Reff, MAG_R, SERSI')
                mrc.write('\n')
                for j in range(1, num_lines-2):
                    if  resl <= f2r.iloc[j][4] <= resu:
                        spc=str(f2r.iloc[j][0])+' , '+str(f2r.iloc[j][1])+' , '+str(round(f2r.iloc[j][2],6))+' , '+str(round(f2r.iloc[j][3],6))+' , '+str(f2r.iloc[j][4])+' , '+str(f2r.iloc[j][5])+' , '+str(f2r.iloc[j][6])+' , '+str(round(f2r.iloc[j][7],4))
                        mrc.write(spc)
                        mrc.write('\n')
                mrc.close()
                os.chdir('/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/')
                z = "Sorted_by_msrb_" + str(i)
                os.makedirs(z, exist_ok=True)
                os.chdir('/home/soumyananda/Dropbox/Goswami_REVISED/CATALOGRADEC_SPLICED/Sorted_lists_SDSS_i/msrb/')
                drc = open(filename, 'r')
                dfz = pd.read_csv(drc, delimiter=',')
                #print(dfz)
                n = len(dfz)
                for k in range(1, n):
                    btemp1= float(dfz.iloc[k][2])
                    ctemp1= float(dfz.iloc[k][3])
                    os.chdir('/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/')
                    for h in range(0,len(file_list)):
                        pt=file_list[h]
                        s=pt.split("_")
                        s1= round(float(s[0]),6)
                        s2= round(float(s[1]),6)
                        if np.isclose(btemp1, s1, rtol=1e-4, atol=1e-5, equal_nan=False) and np.isclose(ctemp1, s2, rtol=1e-4, atol=1e-5, equal_nan=False):
                            #print(btemp1-s1)
                            #print(ctemp1-s2)
                            temp = "/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/" + z + "/"
                            shutil.copy2(str(pt), temp)
                image_concat = []
                temp1 = "/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/" + z + "/"
                os.chdir(temp1)
                image_list = glob.glob('*.fits')
                number_files = len(image_list)
                print('Actual no of files copied=',number_files)
                print("number of objects in file= ", n)
                err_arr = []
                #Comment Beginning
                for counter in range(0, 10):
                    image_concat1 = []
                    os.chdir(temp1)
                    m1 = int(np.floor(number_files / 10))
                    t1 = random.randint(0, number_files)
                    if (number_files - t1 <= (number_files / 10)):
                        t1 = random.randint(0, int(np.floor(9 * number_files / 10)))
                        print('Difference < num_files/10 is true, adjusting upper lim..')
                    t2 = t1 + m1
                    # print("Value of t1=",t1)
                    # print("Value of t2=",t2)
                    for imgno in range(t1, t2):
                        hdul1 = fits.open(str(image_list[imgno]), ignore_missing_end=True,memmap=True)
                        # print(str(image_list[imgno]))
                        image_concat1.append(hdul1[0].data.copy())
                        hdul1.close()
                        del hdul1
                        gc.collect()
                    final_image_temp = np.zeros(shape=image_concat1[0].shape)
                    for image1 in image_concat1:
                        final_image_temp += image1
                    hdu1 = fits.PrimaryHDU(final_image_temp)
                    s111 = "stack_sdssi6as_temp"+ str(counter) + ".fits"
                    os.makedirs(
                        '/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/STACKED_IMAGES_MSRB/',
                        exist_ok=True)
                    os.chdir(
                        "/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/STACKED_IMAGES_MSRB/")
                    hdu1.writeto(s111, overwrite=True)
                    hdul_list_temp = fits.open(s111)
                    imgd_temp = hdul_list_temp[0].data
                    x1 = np.matrix(imgd_temp)
                    x2 = x1.sum()
                    #print("The signal count is=",x2)
                    cts=(x2 /m1) * number_files
                    #print("Total count for all members projected=",cts)
                    err_arr.append(cts)
                    os.remove(s111)
                    #print('Exiting Counter=',counter)
                os.chdir(temp1)
                #Comment ending
                for p in range(0, number_files):
                    hdul = fits.open(str(image_list[p]), ignore_missing_end=True,memmap=True)
                    #print(str(image_list[p]))
                    #print("-------------hdul-----------")
                    #print(log10(np.abs(err_arr[counter1])))
                    #print(type(log10(np.abs(err_arr[counter1]))))
                    image_concat.append(hdul[0].data.copy())
                    hdul.close()
                    del hdul
                    gc.collect()
                final_image = np.zeros(shape=image_concat[0].shape)
                for image in image_concat:
                    final_image += image
                hdu = fits.PrimaryHDU(final_image)
                s11 = "stack_sdssi6as_" + str(i) + ".fits"
                s12 = "stack_sdssi6as_" + str(i) + ".png"
                os.chdir("/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/STACKED_IMAGES_MSRB/")
                shutil.rmtree(temp1)
                hdu.writeto(s11, overwrite=True)
                image_data = fits.getdata(s11, ext=0)
                #print(image_data.shape)
                hdul_list = fits.open(s11)
                imgd = hdul_list[0].data
                x = np.matrix(imgd)
                print("The summed value for the image is:",x.sum())
                t= 30.0-2.5*np.log10(x.sum())
                t_err_temp = 0
                for counter1 in range(0, 10):
                    print("The error counter value is:", err_arr[counter1])
                    t_err_temp = t_err_temp + ((-2.5*np.log10(np.abs(err_arr[counter1]))+30.0) -t)**2.0
                t_err = np.sqrt(t_err_temp / 9.0)
                mid=0.5*(resl+resu)
                d=str(mid)+","+str(t)+","+str(t_err)+","+str(number_files)+"\n"
                f.write(d)
                plt.figure(num=i, figsize=(8, 6), dpi=600, facecolor='w', edgecolor='k')
                plt.imshow(image_data, cmap='gray')
                tempstr = "BY_MSRB" + str(resl) + "_" + str(resu)
                plt.title(tempstr)
                plt.savefig(s12)
                mrc.close()
            f.close()
        if inpval==0:
            os.chdir('/home/soumyananda/Dropbox/Goswami_REVISED/CATALOGRADEC_SPLICED/')
            #execute sort by criteria 3
            three()
            loopcycle=True
            os.chdir('/home/soumyananda/Dropbox/Goswami_REVISED/CATALOGRADEC_SPLICED/Sorted_lists_SDSS_i/seri/')
            f2= open('sorted_seri_catg.csv', 'r')
            f2r = pd.read_csv(f2, sep=',', skipinitialspace=True)
            num_start = f2r.iloc[1][7]
            #print(num_start)
            t= num_lines-2
            num_end = f2r.iloc[t][7]
            #print(num_end)
            val = (num_end - num_start) / 12
            #print(val)
            os.makedirs('/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/STACKED_IMAGES_SERI/',exist_ok=True)
            os.chdir("/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/STACKED_IMAGES_SERI/")
            f=open("countsinfosdssi1.txt","w+")
            f.write("x,y,yerr,num_files\n")
            os.chdir('/home/soumyananda/Dropbox/Goswami_REVISED/CATALOGRADEC_SPLICED/Sorted_lists_SDSS_i/seri/')
            for i in range(1,10):
                resl=num_start + val*(i-1) #incremental dynamic setting for splitting files
                resu = num_start + val * i
                #print(resl,resu)
                #code for binned files for a range of MSRB/SER etc
                os.chdir('/home/soumyananda/Dropbox/Goswami_REVISED/CATALOGRADEC_SPLICED/Sorted_lists_SDSS_i/seri/')
                filename = "seri_part_" + str(i) + ".csv"
                mrc = open(filename, 'w+')
                mrc.write('KIDSFIELD, SEXID, RA, DEC,MSRB, Reff, MAG_R, SERSI')
                mrc.write('\n')
                for j in range(1, num_lines-2):
                    if  resl <= f2r.iloc[j][7] <= resu:
                        spc=str(f2r.iloc[j][0])+' , '+str(f2r.iloc[j][1])+' , '+str(round(f2r.iloc[j][2],6))+' , '+str(round(f2r.iloc[j][3],6))+' , '+str(f2r.iloc[j][4])+' , '+str(f2r.iloc[j][5])+' , '+str(f2r.iloc[j][6])+' , '+str(round(f2r.iloc[j][7],4))
                        mrc.write(spc)
                        mrc.write('\n')
                mrc.close()
                os.chdir('/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/')
                z = "Sorted_by_seri_" + str(i)
                os.makedirs(z, exist_ok=True)
                os.chdir('/home/soumyananda/Dropbox/Goswami_REVISED/CATALOGRADEC_SPLICED/Sorted_lists_SDSS_i/seri/')
                drc = open(filename, 'r')
                dfz = pd.read_csv(drc, delimiter=',')
                #print(dfz)
                n = len(dfz)
                for k in range(1, n):
                    btemp1= float(dfz.iloc[k][2])
                    ctemp1= float(dfz.iloc[k][3])
                    os.chdir('/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/')
                    for h in range(0,len(file_list)):
                        pt=file_list[h]
                        s=pt.split("_")
                        s1= round(float(s[0]),6)
                        s2= round(float(s[1]),6)
                        if np.isclose(btemp1, s1, rtol=1e-4, atol=1e-5, equal_nan=False) and np.isclose(ctemp1, s2, rtol=1e-4, atol=1e-5, equal_nan=False):
                            #print(btemp1-s1)
                            #print(ctemp1-s2)
                            temp = "/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/" + z + "/"
                            shutil.copy2(str(pt), temp)
                image_concat = []
                temp1 = "/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/" + z + "/"
                os.chdir(temp1)
                image_list = glob.glob('*.fits')
                number_files = len(image_list)
                print('Actual no of files copied=',number_files)
                print("number of objects in file= ", n)
                err_arr = []
                #Comment Beginning
                for counter in range(0, 10):
                    image_concat1 = []
                    os.chdir(temp1)
                    m1 = int(np.floor(number_files / 10))
                    t1 = random.randint(0, number_files)
                    if (number_files - t1 <= (number_files / 10)):
                        t1 = random.randint(0, int(np.floor(9 * number_files / 10)))
                        print('Difference < num_files/10 is true, adjusting upper lim..')
                    t2 = t1 + m1
                    # print("Value of t1=",t1)
                    # print("Value of t2=",t2)
                    for imgno in range(t1, t2):
                        hdul1 = fits.open(str(image_list[imgno]), ignore_missing_end=True,memmap=True)
                        # print(str(image_list[imgno]))
                        image_concat1.append(hdul1[0].data.copy())
                        hdul1.close()
                        del hdul1
                        gc.collect()
                    final_image_temp = np.zeros(shape=image_concat1[0].shape)
                    for image1 in image_concat1:
                        final_image_temp += image1
                    hdu1 = fits.PrimaryHDU(final_image_temp)
                    s111 = "stack_sdssi6as_temp"+ str(counter) + ".fits"
                    os.makedirs(
                        '/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/STACKED_IMAGES_SERI/',
                        exist_ok=True)
                    os.chdir(
                        "/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/STACKED_IMAGES_SERI/")
                    hdu1.writeto(s111, overwrite=True)
                    hdul_list_temp = fits.open(s111)
                    imgd_temp = hdul_list_temp[0].data
                    x1 = np.matrix(imgd_temp)
                    x2 = x1.sum()
                    #print("The signal count is=",x2)
                    cts=(x2 / m1) * number_files
                    #print("Total count for all members projected=",cts)
                    err_arr.append(cts)
                    os.remove(s111)
                    #print('Exiting Counter=',counter)
                os.chdir(temp1)
                #Comment ending
                for p in range(0, number_files):
                    hdul = fits.open(str(image_list[p]), ignore_missing_end=True,memmap=True)
                    #print(str(image_list[p]))
                    #print("-------------hdul-----------")
                    #print(log10(np.abs(err_arr[counter1])))
                    #print(type(log10(np.abs(err_arr[counter1]))))
                    image_concat.append(hdul[0].data.copy())
                    hdul.close()
                    del hdul
                    gc.collect()
                final_image = np.zeros(shape=image_concat[0].shape)
                for image in image_concat:
                    final_image += image
                hdu = fits.PrimaryHDU(final_image)
                s11 = "stack_sdssi6as_" + str(i) + ".fits"
                s12 = "stack_sdssi6as_" + str(i) + ".png"
                os.chdir("/home/soumyananda/Dropbox/Goswami_REVISED/DATA_FITS_KIDS_GAMA_VST/SDSS_i_FITS_6AS/STACKED_IMAGES_SERI/")
                shutil.rmtree(temp1)
                hdu.writeto(s11, overwrite=True)
                image_data = fits.getdata(s11, ext=0)
                #print(image_data.shape)
                hdul_list = fits.open(s11)
                imgd = hdul_list[0].data
                x = np.matrix(imgd)
                print("The summed value for the image is:",x.sum())
                t= 30.0-2.5*np.log10(x.sum())
                t_err_temp = 0
                for counter1 in range(0, 10):
                    print("The error counter value is:", err_arr[counter1])
                    t_err_temp = t_err_temp + ((-2.5*np.log10(np.abs(err_arr[counter1]))+30.0) -t)**2.0
                t_err = np.sqrt(t_err_temp / 9.0)
                mid=0.5*(resl+resu)
                d=str(mid)+","+str(t)+","+str(t_err)+","+str(number_files)+"\n"
                f.write(d)
                plt.figure(num=i + 9, figsize=(8, 6), dpi=600, facecolor='w', edgecolor='k')
                plt.imshow(image_data, cmap='gray')
                tempstr = "BY_SERI" + str(resl) + "_" + str(resu)
                plt.title(tempstr)
                plt.savefig(s12)
                mrc.close()
            f.close()

fp.close()