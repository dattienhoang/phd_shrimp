import math
import time
import scipy.optimize as opt
import scipy.ndimage
from scipy import signal, misc
import numpy as np
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib
import matplotlib.gridspec as gridspec
#import matplotlib.patches as mpatches


#import pandas as pd
import trackpy as tp
import pims
import os
import sys

import wx
from sys import exit

#from pandas import DataFrame, Series  # for convenience


wildcard = "Tiff file (*.tif; *.tiff)|*.tif;*.tiff|" \
         "All files (*.*)|*.*"
         
wildcardDataFile = "txt file (*.txt)|*.txt|" \
         "All files (*.*)|*.*"
         


StartingPanel = 1           # 1: AFShrimp, 2: DataFileAnalysis, 3: Image&IntensityPlots


 
#######################################################
# feature find threshold
ThrethholdIntensity = 40000         #Typically 50000



ColorScaleMin = 1100
ColorScaleMax = 2200


SingleFeature_x = 385
SingleFeature_y = 81



additionalFramesN = 1 


FilePrefixInput = '3F_'

#######################################################



#######################################################

magnification = 92.6 # nm / pixel

### Adjacent Shrimp ###
AdjacentShrimp_Start_FrameN = 8
AdjacentShrimp_End_FrameN = 112

NadjacentFrameNs = 6


#######################################################














def MessageBox(self, msg, title):
    """ Display Info Dialog Message """
    font = wx.Font(14, wx.MODERN, wx.NORMAL, wx.NORMAL)
    style = wx.OK | wx.ICON_INFORMATION | wx.STAY_ON_TOP
    dialog = wx.MessageDialog(self, msg, title, style)
    dialog.CenterOnParent()
    dialog.SetFont(font)
    result = dialog.ShowModal()
    if result == wx.ID_OK:
        print dialog.GetFont().GetFaceName()
    dialog.Destroy()
    return












def restart_program():
    """Restarts the current program.
    Note: this function does not return. Any cleanup action (like
    saving data) must be done before calling this function."""
    python = sys.executable
    os.execl(python, python, * sys.argv)




def TowD_PointDensityCalculation(xdata, ydata, oneDsubLength, oneDtotalLength):
    d = oneDsubLength
    l = oneDtotalLength
    
    N = int(l/d)
    
    arr = np.zeros((N,N))
    
    xave = np.mean(xdata)
    xstart = xave - d*(int(N/2))
    yave = np.mean(ydata)
    ystart = yave - d*(int(N/2))
    
    
    for nx in range(N):
        for ny in range(N):
            count = 0
            for k in range(len(xdata)):
                if (xstart+nx*d <= xdata[k] < xstart+(nx+1)*d) and (ystart+ny*d <= ydata[k] < ystart+(ny+1)*d):
                    count += 1
            arr[nx][ny] = count
                
    return arr
    

     
#TowD_PointDensityCalculation([1,1,1,2,3,4,5,5,5,5,1,1], [1,1,1,2,3,4,5,5,5,5,2,2],1,10)       

    
    
    

def Weighted_Ave_Err(values, errors):
    if len(values) != len(errors):
        print 'Error: different lengths'
        return 
        
    values = np.array(values, dtype = float)
    errors = np.array(errors, dtype = float)
    weights = 1/(errors**2)
    
    weightedAve = np.sum(values * weights)/np.sum(weights)
    
    weightedAveErr = np.sqrt(1/np.sum(weights) )
    
    return weightedAve, weightedAveErr
    
    


def AutoCorrFunc(arraydataTemp):
        # Calculating AutoCorrelationFunction
                
        arrayData = np.array(arraydataTemp)
        #yunbiased = arrayData-np.mean(arrayData)
        yunbiased = arrayData
        ynorm = np.sum(yunbiased**2)
        acf = np.correlate(yunbiased, yunbiased, "full")/ynorm
        acf = acf[len(acf)/2:] # use only second half
        #acf2 = np.append(acf[1:], [0])
        return acf


'''
def AutoCorrFunc(xTemp):
    x = np.array(xTemp)
    result = ( np.correlate(x, x, mode='full') )/(np.correlate(x,x))
    return result[result.size/2:]
'''



def RMS_distance_XY_Error(xtemp,ytemp, error):
    x = np.array(xtemp)
    y = np.array(ytemp)
    err = np.array(error)
    xave = np.mean(x)
    yave = np.mean(y)
    distanceSquare = (x - xave)**2 + (y - yave)**2
    distanceSquareErr = 1.414 * np.sqrt(distanceSquare) * err
    
    distanceSquareAve = np.mean(distanceSquare)
    distanceSquareAveErr = np.sqrt(np.sum(distanceSquareErr**2))
    
    RMS = np.sqrt(distanceSquareAve)
    RMS_err = distanceSquareAveErr/(1.414*np.sqrt(distanceSquareAve))
    
    
    #print 'distanceSquare ',distanceSquare
    return RMS, RMS_err
    



def RMS_distance_XY(xtemp,ytemp):
    x = np.array(xtemp)
    y = np.array(ytemp)
    xave = np.mean(x)
    yave = np.mean(y)
    distanceSquare = (x - xave)**2 + (y - yave)**2
    #print 'distanceSquare ',distanceSquare
    return np.sqrt(np.mean(distanceSquare))
    



#print RMS_distance_XY([1,2,3], [1,2,3])


         
def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()
    
    
def SingleExp(x, A, tau, y0):
    return(A*np.exp(-x/tau) + y0)    
    
    
    
    
    
    
    
    
def SingleTwoFrameShrimp_Function(xpos, ypos, frame1, frame2, figFileInfo_text, image1Frames_str, image2Frames_str, colormap, noiseAmp):   
    
    
    
    

    try:
        plt.close('data1_fit_3D')
        plt.close('data2_fit_3D')
        plt.close('data3_fit_3D')
        
    except:
        pass
    
    
    
    
    
    
    
    
    print '\nframe1 \n', frame1
    
    print '\ntype(frame1) ', type(frame1)
    
    print 'len: frame1', len(frame1)
    print 'len: frame1[0]', len(frame1[0])
    
    
    #N_xPixel = len(frame1)
    #N_yPixel = len(frame1[0])
    
    
    ################################################
    x_ImageJ, y_ImageJ = xpos, ypos
    
    #xc = y_ImageJ
    #yc = N_xPixel - x_ImageJ +1
            
    xc = x_ImageJ - 1 # for trackpy 3.0
    yc = y_ImageJ # for trackpy 3.0
        
    ################################################

    

    #print 'filename = ', filename
    #frames = tp.TiffStack(filename)


    '''
    x2 = np.linspace(0, 511, 512)
    y2 = np.linspace(0, 511, 512)
    x2, y2 = np.meshgrid(x2, y2)
    '''
    
    
    NCenterSubPixel = 10
    
    
    if noiseAmp != 0:
        noise1 = np.random.normal(loc=0, scale=noiseAmp, size=(21,21))
        noise2 = np.random.normal(loc=0, scale=noiseAmp, size=(21,21))    
    else:
        noise1 = 0
        noise2 = 0


    frameTemp1 = frame1[yc-NCenterSubPixel:yc+NCenterSubPixel+1, xc-NCenterSubPixel:xc+NCenterSubPixel+1] + noise1
    frameTemp2 = frame2[yc-NCenterSubPixel:yc+NCenterSubPixel+1, xc-NCenterSubPixel:xc+NCenterSubPixel+1] + noise2
    frameTempDiff = np.array(frameTemp1, dtype=np.int32) - np.array(frameTemp2, dtype=np.int32)


    
    print 'frameTemp1 \n', frameTemp1
    
    print '\nframeTempDiff \n', frameTempDiff
    print 'np.sum(frameTempDiff) = ', np.sum(frameTempDiff)
    if np.sum(frameTempDiff) < 0:
        frameTempDiff *= (-1)
    
    
    
    
    x3 = np.linspace(0, 2*NCenterSubPixel, 2*NCenterSubPixel+1)
    y3 = np.linspace(0, 2*NCenterSubPixel, 2*NCenterSubPixel+1)
    x3, y3 = np.meshgrid(x3, y3)
    
    
    
    x3fit = np.linspace(0, 20, 210)
    y3fit = np.linspace(0, 20, 210)
    x3fit, y3fit = np.meshgrid(x3fit, y3fit)
    






    frameTemp1_raveled = frameTemp1.ravel()  
    frameTemp2_raveled = frameTemp2.ravel()  
    frameTempDiff_raveled = frameTempDiff.ravel()                


    # fit the data  for all frames
 
 
    k = 1
   
    Afit = 0.0
    x0fit = 0.0
    y0fit = 0.0
    z0fit = 0.0
    sigma_x = 0.0
    sigma_y = 0.0
    theta = 0.0
    
    
    
    
    
    #try:
        
    initial_guess = (6000, NCenterSubPixel+1,  NCenterSubPixel+1, 3, 3, 10, 2000)   
    popt, pcov = opt.curve_fit(twoD_Gaussian, (x3, y3), frameTemp1_raveled, p0 = initial_guess)
    data_fitted_frameTemp1 = twoD_Gaussian((x3fit, y3fit), *popt)

    Afit = popt[0]
    x0fit = popt[1]
    y0fit = popt[2]
    
    x0fitImageJ1 = xpos + (NCenterSubPixel - y0fit)
    y0fitImageJ1 = ypos - (NCenterSubPixel - x0fit)
    
    z0fit = popt[6]

    sigma_x = popt[3]
    sigma_y = popt[4]
    theta = popt[5]
    # AF gshrimp AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaErr
    AfitErr = np.sqrt(pcov[0,0])                    
    x0fitImageJErr = np.sqrt(pcov[1,1])
    y0fitImageJErr = np.sqrt(pcov[2,2]) 
    z0fitErr = np.sqrt(pcov[6,6]) 
    sigma_xErr = np.sqrt(pcov[3,3]) 
    sigma_yErr = np.sqrt(pcov[4,4]) 
    thetaErr = np.sqrt(pcov[5,5])
       
       
    ########### eccentricicity calculation  ############
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2) 
    
    #print 'a, b, c =', a,b,c

    A = a
    B = 2*b
    C = c
            
    eccentricityTemp1 = np.sqrt( 2*(np.sqrt((A-C)**2 + B**2)) / (A+C + np.sqrt(((A-C)**2 + B**2)))  )
    #print 'eccentricityTemp = ', eccentricityTemp    
    
    xye1_str = ( '('+ str(np.around(x0fitImageJ1,3)) +', ' + str(np.around(y0fitImageJ1,3)) + 
    '),  Err: (' + str(np.around(x0fitImageJErr,3))+', '+str(np.around(y0fitImageJErr,3)) + 
    ')\n    e: ' + str(np.around(eccentricityTemp1,3)) 
    )    
    
    xy0fitImageJ1Err = np.around( (x0fitImageJErr + y0fitImageJErr)/2.0 ,3)  
    
    
    BG1fit = z0fit
    A1fit = Afit
    
    print 'z0fit, Afit', z0fit, Afit
    




    
    popt, pcov = opt.curve_fit(twoD_Gaussian, (x3, y3), frameTemp2_raveled, p0 = initial_guess)
    data_fitted_frameTemp2 = twoD_Gaussian((x3fit, y3fit), *popt)

    Afit = popt[0]
    x0fit = popt[1]
    y0fit = popt[2]
    
    x0fitImageJ2 = xpos + (NCenterSubPixel - y0fit)
    y0fitImageJ2 = ypos - (NCenterSubPixel - x0fit)
    
    z0fit = popt[6]

    sigma_x = popt[3]
    sigma_y = popt[4]
    theta = popt[5]
        # AF gshrimp AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaErr
    AfitErr = np.sqrt(pcov[0,0])                    
    x0fitImageJErr = np.sqrt(pcov[1,1])
    y0fitImageJErr = np.sqrt(pcov[2,2]) 
    z0fitErr = np.sqrt(pcov[6,6]) 
    sigma_xErr = np.sqrt(pcov[3,3]) 
    sigma_yErr = np.sqrt(pcov[4,4]) 
    thetaErr = np.sqrt(pcov[5,5])
       
       
    ########### eccentricicity calculation  ############
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2) 
    
    #print 'a, b, c =', a,b,c

    A = a
    B = 2*b
    C = c
            
    eccentricityTemp2 = np.sqrt(2*(np.sqrt((A-C)**2 + B**2)) / (A+C + np.sqrt(((A-C)**2 + B**2)))  )
    #print 'eccentricityTemp = ', eccentricityTemp    
    
    xye2_str = ( '('+ str(np.around(x0fitImageJ2,3)) +', ' + str(np.around(y0fitImageJ2,3))+
    '),  Err: (' + str(np.around(x0fitImageJErr,3))+', '+str(np.around(y0fitImageJErr,3)) + 
     ')\n    e: ' + str(np.around(eccentricityTemp2,3)) 
    )
    
    xy0fitImageJ2Err = np.around( (x0fitImageJErr + y0fitImageJErr)/2.0 ,3)  
    
    
    
    
    
    
    BG2fit = z0fit
    A2fit = Afit
    
    print 'z0fit, Afit', z0fit, Afit
    
    
   


    
    
    
        
    initial_guess = (3000, NCenterSubPixel+1,  NCenterSubPixel+1, 3, 3, 10, 200)   
        
    popt, pcov = opt.curve_fit(twoD_Gaussian, (x3, y3), frameTempDiff_raveled, p0 = initial_guess)
    data_fitted_frameTempDiff = twoD_Gaussian((x3fit, y3fit), *popt)
  
    
    Afit = popt[0]
    x0fit = popt[1]
    y0fit = popt[2]
    
    x0fitImageJdiff = xpos + (NCenterSubPixel - y0fit)
    y0fitImageJdiff = ypos - (NCenterSubPixel - x0fit)
    
    z0fit = popt[6]

    sigma_x = popt[3]
    sigma_y = popt[4]
    theta = popt[5]


    # AF gshrimp AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaErr
    AfitErr = np.sqrt(pcov[0,0])                    
    x0fitImageJErr = np.sqrt(pcov[1,1])
    y0fitImageJErr = np.sqrt(pcov[2,2]) 
    z0fitErr = np.sqrt(pcov[6,6]) 
    sigma_xErr = np.sqrt(pcov[3,3]) 
    sigma_yErr = np.sqrt(pcov[4,4]) 
    thetaErr = np.sqrt(pcov[5,5])
                    
                    
    # AF gshrimp Percentage of the errors: AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaErr
    AfitErrP = np.around(100.0 * AfitErr / Afit, 0)                    
    x0fitImageJErrPx = np.around(np.sqrt(pcov[1,1]),2)
    y0fitImageJErrPx = np.around(np.sqrt(pcov[2,2]),2) 
    z0fitErrP = np.around(100.0 * z0fitErr / z0fit, 0)      
    sigma_xErrP = np.around(100.0 * sigma_xErr / sigma_x, 0)   
    sigma_yErrP = np.around(100.0 * sigma_yErr / sigma_y, 0)  
    thetaErrDeg = np.around(np.sqrt(pcov[5,5]),2)                       
    

                     
                    
    ########### eccentricicity calculation  ############
   
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2) 
    
    #print 'a, b, c =', a,b,c

    A = a
    B = 2*b
    C = c
            
    eccentricityTempdiff = np.sqrt(2*(np.sqrt((A-C)**2 + B**2)) / (A+C + np.sqrt(((A-C)**2 + B**2)))  )
    #print 'eccentricityTemp = ', eccentricityTemp                        

    xyediff_str = ( '('+ str(np.around(x0fitImageJdiff,3)) +', ' + str(np.around(y0fitImageJdiff,3)) +
    '),  Err: (' + str(np.around(x0fitImageJErr,3))+', '+str(np.around(y0fitImageJErr,3)) + 
    ')\n   e: ' + str(np.around(eccentricityTempdiff,3)) + ' '
    )
    
    xy0fitImageJdiffErr = np.around( (x0fitImageJErr + y0fitImageJErr)/2.0 ,3)  
    
    
    
    
    BG3fit = z0fit
    A3fit = Afit
    
    print 'z0fit, Afit', z0fit, Afit
    
    
   
   


    

    
    
    
    
    
    print 'figFileInfo_text ', figFileInfo_text
    
    
    fig_TwoFrameShrimp = plt.figure('TwoFrameShrimp')
    fig_TwoFrameShrimp.clear()
    
    fig_TwoFrameShrimp.subplots_adjust(top = 0.86, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.2, hspace = 0.2)
    plt.figtext(0.5, 0.95, figFileInfo_text + '\n' + 'Noise: ' + str(noiseAmp)
            ,ha='center', va='bottom', color='black', weight='normal', size=12)
    fig_TwoFrameShrimp.text(0.05, 0.96, 'x, y : (' + str(xpos)+', '+str(ypos)+')', ha="left", va="bottom", size="medium",color="red")
            
    
    mngr = plt.get_current_fig_manager()
    mngr.window.setGeometry(105,30,1500, 1000) # left two for position of the window, right two window size
    
    plt.subplot(231, aspect = 1.0)
    plt.title(image1Frames_str + '\n' + xye1_str)
    plt.imshow(frameTemp1, interpolation = 'None', cmap=colormap, origin = 'bottom')
    plt.colorbar()
    plt.contour(x3fit, y3fit, data_fitted_frameTemp1.reshape(210, 210), 8, colors='w')
    
    
    plt.subplot(232, aspect = 1.0)
    plt.title(image2Frames_str + '\n' + xye2_str)
    plt.imshow(frameTemp2, interpolation = 'None', cmap=colormap, origin = 'bottom')
    plt.colorbar()
    plt.contour(x3fit, y3fit, data_fitted_frameTemp2.reshape(210, 210), 8, colors='w')
    
    plt.subplot(233, aspect = 1.0)
    plt.title('SHRImP\n' + xyediff_str)
    plt.imshow(frameTempDiff, interpolation = 'None', cmap=colormap, origin = 'bottom')
    plt.colorbar()
    plt.contour(x3fit, y3fit, data_fitted_frameTempDiff.reshape(210, 210), 8, colors='w')
    
    
    
    
    
    
    
    
    
    plt.subplot(236, aspect = 1.0)
    plt.title('FIONA & SHRImP')

    markerColorAFS = ['r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k','r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k','r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k','r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k']
    
    plt.scatter(x0fitImageJ1, y0fitImageJ1, marker = '.',s= 100, edgecolors = markerColorAFS[0], c = markerColorAFS[0]
    , label = 'FIONA #1' )
    EmissionSite_Cir1 = plt.Circle((x0fitImageJ1, y0fitImageJ1), radius = xy0fitImageJ1Err, color = markerColorAFS[0], fill=False)
    fig_TwoFrameShrimp.gca().add_artist(EmissionSite_Cir1)

    plt.scatter(x0fitImageJ2, y0fitImageJ2, marker = '.',s= 100, edgecolors = markerColorAFS[1], c = markerColorAFS[1]
    , label = 'FIONA #2' )
    EmissionSite_Cir2 = plt.Circle((x0fitImageJ2, y0fitImageJ2), radius = xy0fitImageJ2Err, color = markerColorAFS[1], fill=False)
    fig_TwoFrameShrimp.gca().add_artist(EmissionSite_Cir2)

    
    plt.scatter(x0fitImageJdiff, y0fitImageJdiff, marker = 'x',s= 100, edgecolors = markerColorAFS[2], c = markerColorAFS[2]
    , label = 'SHRImP' )
    EmissionSite_Cir3 = plt.Circle((x0fitImageJdiff, y0fitImageJdiff), radius = xy0fitImageJdiffErr, color = markerColorAFS[2], fill=False)
    fig_TwoFrameShrimp.gca().add_artist(EmissionSite_Cir3)

    


    plt.legend(loc='upper right', prop={'size':12})
    plt.xlim( xpos - 2, xpos + 4 ) 
    plt.ylim( ypos - 2, ypos + 4 )
       
    
    
    
    
    
    plt.show()
    
    

    fig1_3d = plt.figure('data1_fit_3D')
    mngr = plt.get_current_fig_manager()
    mngr.window.setGeometry(20,30,600, 600) # left two for position of the window, right two window size
    ax1 = fig1_3d.gca(projection='3d')
    ax1.plot_surface(x3, y3, frameTemp1, rstride=1, cstride=1, cmap=colormap)
    ax1.plot_wireframe(x3fit, y3fit, data_fitted_frameTemp1.reshape(210, 210), rstride=10, cstride=10, colors='Lime', linewidths=2)
    ax1.set_zlim(BG1fit*0.9, BG1fit + A1fit*1.1)
    ax1.zaxis.set_major_locator(LinearLocator(10))
    ax1.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    
    
    
    fig2_3d = plt.figure('data2_fit_3D')
    mngr = plt.get_current_fig_manager()
    mngr.window.setGeometry(620,30,600, 600) # left two for position of the window, right two window size
    ax2 = fig2_3d.gca(projection='3d')
    ax2.plot_surface(x3, y3, frameTemp2, rstride=1, cstride=1, cmap=colormap)
    ax2.plot_wireframe(x3fit, y3fit, data_fitted_frameTemp2.reshape(210, 210), rstride=10, cstride=10, colors='Lime', linewidths=2)
    ax2.set_zlim(BG1fit*0.9, BG1fit + A1fit*1.1)
    ax2.zaxis.set_major_locator(LinearLocator(10))
    ax2.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
        
    

    

    fig2_3d = plt.figure('data3_fit_3D')
    mngr = plt.get_current_fig_manager()
    mngr.window.setGeometry(1200,30,600, 600) # left two for position of the window, right two window size
    ax3 = fig2_3d.gca(projection='3d')
    ax3.plot_surface(x3, y3, frameTempDiff, rstride=1, cstride=1, cmap=colormap)
    ax3.plot_wireframe(x3fit, y3fit, data_fitted_frameTempDiff.reshape(210, 210), rstride=10, cstride=10, colors='Lime', linewidths=2)
    ax3.set_zlim(BG1fit*0.0, BG1fit*0 + A1fit*1.1)
    ax3.zaxis.set_major_locator(LinearLocator(10))
    ax3.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    

    
    
    '''
    
    
    
    
    
    
    # fit the data 
    initial_guess = (6000,10,10,3,3,10,200)
    
    frameTemp0_raveled = frameTemp0.ravel()
    
    popt, pcov = opt.curve_fit(twoD_Gaussian, (x3, y3), frameTemp0_raveled, p0 = initial_guess)
    
    
    x3fit = np.linspace(0, 20, 210)
    y3fit = np.linspace(0, 20, 210)
    x3fit, y3fit = np.meshgrid(x3fit, y3fit)
    
    data_fitted = twoD_Gaussian((x3fit, y3fit), *popt)
    
    #fig, ax = plt.subplots(1, 1)
    
    plt.figure('data_fit_2D')
    mngr = plt.get_current_fig_manager()
    mngr.window.setGeometry(5,30,800, 800) # left two for position of the window, right two window size
    plt.hold(True)
    plt.imshow(frameTemp0.reshape(21, 21), interpolation = 'None', cmap=plt.cm.jet, origin='bottom',
        extent=(x3.min(), x3.max(), y3.min(), y3.max()))
    plt.contour(x3fit, y3fit, data_fitted.reshape(210, 210), 8, colors='w')
    
    
    
    #data2d = data_fitted.reshape(210,210)
    
    
    fig5 = plt.figure('data_fit_3D')
    mngr = plt.get_current_fig_manager()
    mngr.window.setGeometry(800,30,800, 800) # left two for position of the window, right two window size
    ax5 = fig5.gca(projection='3d')
    surf5 = ax5.plot_surface(x3, y3, frameTemp0, rstride=1, cstride=1, cmap=cm.jet)
    p = ax5.plot_wireframe(x3fit, y3fit, data_fitted.reshape(210, 210), rstride=10, cstride=10, colors='Lime')
    ax5.set_zlim(1000, 16000)
    ax5.zaxis.set_major_locator(LinearLocator(10))
    ax5.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    
    
    plt.show()
    
    
    
    
    # fit the data  for all frames
    initial_guess = (6000,10,10,3,3,10,200)
    
    
    
    matplotlib.rcParams.update({'font.size': 11})
    plt.figure('data_fit_2D_frame_0--', figsize=(19,9))
    k = 1
    m = 0
    A = 0
    x0fit = 0
    y0fit = 0
    z0fit = 0
    
    xyfit = '(' + str(np.around(x0fit,2)) + ', ' + str(np.around(y0fit,2)) + ')'

    eccentricityAll = []
    
    for n in range(29):
        print 'frame = ', n
        
        frameTemp0 = frames[n][yc-10:yc+11, xc-10:xc+11]
        
        frameTemp0_raveled = frameTemp0.ravel()
        
        try:
            popt, pcov = opt.curve_fit(twoD_Gaussian, (x3, y3), frameTemp0_raveled, p0 = initial_guess)
            print popt
            A = popt[0]
            x0fit = popt[1]
            y0fit = popt[2]
            
            x0fitImageJ = xpos - (10.0 - y0fit)
            y0fitImageJ = ypos - (10.0 - x0fit)
            
            z0fit = popt[6]
            xyfit = '(' + str(np.around(x0fitImageJ,2)) + ', ' + str(np.around(y0fitImageJ,2)) + ')'


            ########### eccentricicity calculation  ############
            sigma_x = popt[3]
            sigma_y = popt[4]
            theta = popt[5]
            z0fit = popt[6]
            
            a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
            b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
            c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2) 
            
            print 'a, b, c =', a,b,c
            
            
            A = a
            B = 2*b
            C = c
                    
            eccentiricityTemp = np.sqrt(2*(np.sqrt((A-C)**2 + B**2)) / (A+C + np.sqrt(((A-C)**2 + B**2)))  )
            print 'eccentiricityTemp = ', eccentiricityTemp
            
            eccentricityAll.append(eccentiricityTemp)
            
            ####################################################
            




  
        except:
            print 'fitting failed'
            eccentiricityTemp = None
            
        x3fit = np.linspace(0, 20, 210)
        y3fit = np.linspace(0, 20, 210)
        x3fit, y3fit = np.meshgrid(x3fit, y3fit)
        
        data_fitted = twoD_Gaussian((x3fit, y3fit), *popt)
        
        #fig, ax = plt.subplots(1, 1)
        
        
        plt.subplot(2,5,k)
        plt.subplots_adjust(top = 0.95, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.1, hspace = 0.2)
        plt.title('frame # ' + str(n) + '   e = ' + str(np.around(eccentiricityTemp, 2)) + '\nAfit = '+ str(np.around(A,-1)) + '    B.G. ='+str(np.around(z0fit)) 
        +'\nfit: x, y = ' + xyfit )
        plt.hold(True)
        plt.imshow(frameTemp0.reshape(21, 21), interpolation = 'None', cmap=plt.cm.jet, origin='bottom',
            extent=(x3.min(), x3.max(), y3.min(), y3.max()))
        plt.contour(x3fit, y3fit, data_fitted.reshape(210, 210), 8, colors='w')
    
        
        
        if k%10 == 0:
            m += 10
            plt.figure('data_fit_2D_frame_' + str(m + 1) + '--', figsize=(19,9) )
            k = 0
        
        k += 1
        
        
    plt.show()
    
    '''
    
     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
def ManyTwoDFits(filename, xpos, ypos, Nframe):   
    colorScaleMin = 1500
    colorScaleMax = 14000
    AbsoluteColorScaleYN = 'n' # y or n

    colorMap = 'jet' # default = jet or spectral

    print 'file name = ', os.path.basename(filename)    
    
    DataFrameAll = tp.TiffStack(filename)
    
    #N_xPixel = DataFrameAll.frame_shape[0]
    #N_yPixel = DataFrameAll.frame_shape[1]
    
    #print N_xPixel
    
    DataFrameOne = DataFrameAll[Nframe]
    
    #a = DataFrameOne[yc-10:yc+11, xc-10:xc+11]
    x3 = np.linspace(0, 20, 21)
    y3 = np.linspace(0, 20, 21)
    x3, y3 = np.meshgrid(x3, y3)
    
    
    
    # fit the data  for all frames
    initial_guess = (6000,10,10,3,3,10,1000)
    
    
    
    matplotlib.rcParams.update({'font.size': 11})
    plt.figure('data_fit_2D_Feature_N_0--', figsize=(19,9))
    plt.figtext(0.5, 0.97, str(os.path.basename(filename)),ha='center', color='black', weight='bold', size='medium')
    k = 1
    m = 0
    A = 0
    x0fit = 0
    y0fit = 0
    z0fit = 0
    
    xyfit = '(' + str(np.around(x0fit,2)) + ', ' + str(np.around(y0fit,2)) + ')'
    
    Nfeatures = len(xpos)
    AFitAll = []
    x0FitAll = []
    y0FitAll = []
    xSigmaFitAll = []
    ySigmaFitAll = []
    thetaFitAll = []
    z0FitAll = []
    eccentricityAll = []
    for n in range(Nfeatures):
        print 'Feature # = ', n
        
        ############################
        #xc = ypos[n]
        #yc = N_xPixel - xpos[n]
            
        xc = xpos - 1 # for trackpy 3.0
        yc = ypos # for trackpy 3.0
        ####################################

        
        print 'xpos, ypos = ', xpos[n], ypos[n]
        print 'xc, yc = ', xc, yc
        featureImg_temp = DataFrameOne[yc-10:yc+11, xc-10:xc+11]
        print 'featureImg_temp.size : ', featureImg_temp.size
        '''
        plt.figure('temp')
        plt.imshow(featureImg_temp)
        plt.show()
        exit()
        '''
        
        featureImg_temp_raveled = featureImg_temp.ravel()
        
        try:
            popt, pcov = opt.curve_fit(twoD_Gaussian, (x3, y3), featureImg_temp_raveled, p0 = initial_guess)
            print 'ManyTwoDFits'
            print popt
            Afit = popt[0]
            x0fit = popt[1]
            y0fit = popt[2]
            
            x0fitImageJ = xpos[n] - (10.0 - y0fit) 
            y0fitImageJ = ypos[n] - (10.0 - x0fit)
            #centerImg_temp = featureImg_temp[y0fit-1:y0fit+2, x0fit-1:x0fit+1]
            '''
            plt.imshow(centerImg_temp)
            plt.show()
            '''
            
            
            sigma_x = popt[3]
            sigma_y = popt[4]
            theta = popt[5]
            z0fit = popt[6]
            
            AFitAll.append(popt[0])
            x0FitAll.append(popt[1])
            y0FitAll.append(popt[2])
            xSigmaFitAll.append(popt[3])
            ySigmaFitAll.append(popt[4])
            thetaFitAll.append(popt[5])
            z0FitAll.append(popt[6])
                        
            
            a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
            b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
            c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2) 
            
            print 'a, b, c =', a,b,c
            
            
            A = a
            B = 2*b
            C = c
                    
            eccentricityTemp = np.sqrt(2*(np.sqrt((A-C)**2 + B**2)) / (A+C + np.sqrt(((A-C)**2 + B**2)))  )
            print 'eccentricityTemp = ', eccentricityTemp
            
            eccentricityAll.append(eccentricityTemp)
            
            xyfit = '(' + str(np.around(x0fitImageJ,2)) + ', ' + str(np.around(y0fitImageJ,2)) + ')'
            x3fit = np.linspace(0, 20, 201)
            y3fit = np.linspace(0, 20, 201)
            x3fit, y3fit = np.meshgrid(x3fit, y3fit)
            
            data_fitted = twoD_Gaussian((x3fit, y3fit), *popt)
            
            if n < 50 :
                plt.subplot(2,5,k)
                plt.subplots_adjust(top = 0.9, bottom = 0.02, left = 0.05, right = 0.95, wspace = 0.1, hspace = 0.2)
                plt.title('feature # ' + str(n) + '      e = '+ str(np.around(eccentricityTemp,2)) +
                    '\nAfit = '+ str(np.around(Afit,-1)) + '    B.G. ='+str(np.around(z0fit)) +
                    '\nfit: x, y = ' + xyfit )
                plt.hold(True)
    
        except:
            print 'fitting failed'
            if n < 50 :
                plt.subplot(2,5,k)
                plt.subplots_adjust(top = 0.9, bottom = 0.02, left = 0.05, right = 0.95, wspace = 0.1, hspace = 0.2)
                plt.title('feature # ' + str(n) + '\nno fit\n')
                plt.hold(True)
        #fig, ax = plt.subplots(1, 1)
        
        
        
        
        try:
            if n < 50 :
                if AbsoluteColorScaleYN == 'y':
                    plt.imshow(np.rot90(featureImg_temp_raveled.reshape(21, 21) ,3), 
                               interpolation = 'None', cmap=colorMap, origin='upper', vmin=colorScaleMin, vmax=colorScaleMax,                               
                               extent=(x3.min(), x3.max()+1, y3.min(), y3.max()+1))
                else:
                    plt.imshow(np.rot90(featureImg_temp_raveled.reshape(21, 21) ,3), 
                               interpolation = 'None', cmap=colorMap, origin='upper',
                               extent=(x3.min(), x3.max()+1, y3.min(), y3.max()+1))
                    
                plt.contour(x3fit, y3fit, data_fitted.reshape(201, 201), 8, colors='w')
    
        except:
            if n < 50 :
                x3 = np.linspace(0, 5, 6)
                y3 = np.linspace(0, 5, 6)
                x3, y3 = np.meshgrid(x3, y3)
                featureImg_temp = DataFrameOne[yc-2:yc+3, xc-2:xc+3]
                featureImg_temp_raveled = featureImg_temp.ravel()
                
                if AbsoluteColorScaleYN == 'y':
                    plt.imshow(np.rot90(featureImg_temp_raveled.reshape(5, 5) ,3), 
                               interpolation = 'None', cmap=colorMap, origin='upper', vmin=colorScaleMin, vmax=colorScaleMax,
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))               
                else:
                    plt.imshow(np.rot90(featureImg_temp_raveled.reshape(5, 5) ,3), 
                               interpolation = 'None', cmap=colorMap, origin='upper',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                    
                    
                print 'no fitting image'
                
                x3 = np.linspace(0, 20, 21)
                y3 = np.linspace(0, 20, 21)
                x3, y3 = np.meshgrid(x3, y3)        
        
        if n < 49 :
            if k%10 == 0:
                m += 10
                plt.figure('data_fit_2D_frame_' + str(m + 1) + '--', figsize=(19,9) )
                plt.figtext(0.5, 0.97, str(os.path.basename(filename)),ha='center', color='black', weight='bold', size='medium')
                k = 0
            
            k += 1
    
    #print 'eccentricityAll \n', eccentricityAll    

    plt.figure('eccentricity')
    plt.subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.1, hspace = 0.2)
    plt.figtext(0.5, 0.97, str(os.path.basename(filename)),ha='center', color='black', weight='bold', size='medium')
    mngr = plt.get_current_fig_manager()
    mngr.window.setGeometry(50,30,1200, 800)
    
    plt.subplot(2,2,1)
    plt.title('eccentricity plot for all features')
    plt.plot(eccentricityAll) 
    
    plt.subplot(2,2,3)
    plt.title('eccentricity histogram for all features')    
    plt.hist(eccentricityAll)
    
    
    
    
    # collecting only SBR > 2
    eccentricityForSBRgreaterOne = []
    for n in range(len(eccentricityAll)):
        if AFitAll[n] / float(z0FitAll[n]) > 1.5:
            eccentricityForSBRgreaterOne.append(eccentricityAll[n])
        else:
            pass
    
      
    plt.subplot(2,2,4)
    plt.title('eccentricity_Histogram_for_SBR_greater_than_ 1.5')    
    if len(eccentricityForSBRgreaterOne) != 0:
        plt.hist(eccentricityForSBRgreaterOne)
    
    
    plt.ion()
    plt.show()
    







def FastIntensityHistogramFunc(frame, xpos, ypos):

    print '\n Running Fast Intensity Histogram Func' 
    
    
    

    ########################################################
    # mismatching correction between python numpy and ImageJ

    
    #N_xPixel = frames.frame_shape[0]
    #N_yPixel = frames.frame_shape[1]

    #N_xPixel = len(frame[0])
    #N_yPixel = len(frame[0][0])
    
    
    x_ImageJ, y_ImageJ = np.array(xpos), np.array(ypos)
    
    #xc = y_ImageJ
    #yc = N_xPixel - x_ImageJ +1
    
        
    xc = x_ImageJ - 1 # for trackpy 3.0
    yc = y_ImageJ # for trackpy 3.0

    ########################################################
    
    

    

    centerIntensityMax5Ave = [] 
    

    frametemp = frame
    for n in range(len(xc)):
        
        frameCenterTemp = frametemp[yc[n]-4:yc[n]+5, xc[n]-4:xc[n]+5]
        print 'n, ', n
        print 'frameCenterTemp\n ', frameCenterTemp
        
        try:
            FiveMaxAveTemp = int(np.mean( np.sort(frameCenterTemp.ravel())[-5:]) )
            centerIntensityMax5Ave.append(FiveMaxAveTemp)
        except:
            pass
        
        

        
            

    print '\n Endding  Fast Intensity Histogram Func \n\n'       
    
    return centerIntensityMax5Ave
    
    '########## Fast Intensity Histogram Func ####################'
    '#####################################################################'










def FastIntensityTrajectoryImagePlotsFunc(loadedframes, filename, xpos, ypos, NframeStart, NtotalImages):

    print '\n Running Fast IntensityTrajectoryImagePlotsFunc' 
    print 'filename = ', filename
    frames = loadedframes
    
    frameStart = 0
    frameEnd = len(frames) - 1
    
    print 'len: frames ', len(frames)    
    
    
    
    ########################################################
    # mismatching correction between python numpy and ImageJ

    #N_xPixel = frames.frame_shape[0]
    #N_yPixel = frames.frame_shape[1]
    
    x_ImageJ, y_ImageJ = np.array(xpos), np.array(ypos)
    
    #xc = y_ImageJ
    #yc = N_xPixel - x_ImageJ +1
        
    xc = x_ImageJ - 1 # for trackpy 3.0
    yc = y_ImageJ # for trackpy 3.0

    ########################################################

    
    
    featureImages = [[] for _ in xrange(len(xc))] 
    featureImagesFrameN = [[] for _ in xrange(len(xc))] 
    
    
    centerIntensityIntegratedTrajectory = [[] for _ in xrange(len(xc))] 
    
    
    for k in range(frameStart,frameEnd+1):
        print '\n frame # ', k
        frametemp = frames[k]
        for n in range(len(xc)):
            
            frameLargeTemp = frametemp[yc[n]-7:yc[n]+8, xc[n]-7:xc[n]+8]
            frameCenterTemp = frametemp[yc[n]-6:yc[n]+7, xc[n]-6:xc[n]+7]
            
            
            BGave = ( np.sum(frameLargeTemp.ravel()) - np.sum(frameCenterTemp.ravel()) ) / 56
            
            
            IntegratedIntesityTemp = np.sum(frameCenterTemp.ravel()) - BGave*169
            
            centerIntensityIntegratedTrajectory[n].append(IntegratedIntesityTemp)
            
            if NframeStart <= k < NframeStart + NtotalImages:
                featureImages[n].append(frameLargeTemp)
                featureImagesFrameN[n].append(k)
                

    print '\n Endding  IntensityTrajectoryImagePlotsFunc\n\n'       
    
    return centerIntensityIntegratedTrajectory, featureImages, featureImagesFrameN
    '########## FastIntensityTrajectoryImagePlotsFunc ####################'
    '#####################################################################'
















'# do not use'
def IntensityTrajectoryPlotsFunc(filename, xpos, ypos):
    print '\n Runnging IntensityTrajectoryPlotsFunc'
    
    print 'filename = ', filename
    frames = tp.TiffStack(filename)
    
    
    frameStart = 0
    frameEnd = len(frames) - 1
    
    print 'len: frames ', len(frames)    
    
    
    
    ########################################################
    # mismatching correction between python numpy and ImageJ

    #N_xPixel = frames.frame_shape[0]
    #N_yPixel = frames.frame_shape[1]
    
    x_ImageJ, y_ImageJ = np.array(xpos), np.array(ypos)
    
    #xc = y_ImageJ
    #yc = N_xPixel - x_ImageJ +1
        
    xc = x_ImageJ - 1 # for trackpy 3.0
    yc = y_ImageJ # for trackpy 3.0

    ########################################################
    
    
    #print 'frameStart ', frameStart
    #print 'frameEnd ', frameEnd
    
    

    


    
    '''
    # calculating each feature for all frames => slow
    centerIntensityAveTrajectory = []
    centerIntensityMaxTrajectory = []
    for kTemp in range(len(xc)):
        print 'Intensity Trajectory feature # ', kTemp
        xctemp = xc[kTemp]
        yctemp = yc[kTemp]      
        
        centerIntensityAveTemp = []
        centerIntensityMaxTemp = []
        
        for nTemp in np.arange(frameStart, frameEnd+1):
            frameCenterTemp = frames[nTemp][yctemp-2:yctemp+3, xctemp-2:xctemp+3]
            centerIntensityAveTemp.append(np.sum(frameCenterTemp))
            centerIntensityMaxTemp.append(np.amax(frameCenterTemp))
            
        centerIntensityAveTrajectory.append(centerIntensityAveTemp)   
        centerIntensityMaxTrajectory.append(centerIntensityMaxTemp)    
    '''
    
    # calculating each frame for all features => 
    centerIntensityAveTrajectory = [[] for _ in xrange(len(xc))] 
    centerIntensityMaxTrajectory = [[] for _ in xrange(len(xc))]
    centerIntensityMax5AveTrajectory = [[] for _ in xrange(len(xc))]
    for k in range(frameStart,frameEnd+1):
        print 'frame # ', k
        frametemp = frames[k]
        for n in range(len(xc)):
            frameCenterTemp = frametemp[yc[n]-3:yc[n]+4, xc[n]-3:xc[n]+4]
            FiveMaxAveTemp = int(np.mean( np.sort(frameCenterTemp.ravel())[-5:]) )
            
            centerIntensityAveTrajectory[n].append(np.mean(frameCenterTemp))
            centerIntensityMaxTrajectory[n].append(np.amax(frameCenterTemp)) 
            centerIntensityMax5AveTrajectory[n].append(FiveMaxAveTemp)
            
           
    return centerIntensityAveTrajectory, centerIntensityMaxTrajectory, centerIntensityMax5AveTrajectory
    




























def FastFIONA(LoadedFrames, filepath, xpos, ypos, frameStart, frameEnd, SvsB_Tolerance, maxEcc, maxPer, maxPx, additionalFramesN):   
    # using all features for all frames by analyzing frame by frame

    
    #FastFIONA(self.filepath, xpos, ypos, int(self.Nframe.GetValue()), int(self.Nframe.GetValue()) + NframesFiona -1 , 1.1, 0.6, 10, 0.2, 0)
        
         


    
    framesTemp = LoadedFrames
    
    if additionalFramesN >= 1:
            
        frames = []
        frames.append( (framesTemp[0] + framesTemp[1])/2.0)
        for n in range(1, frameEnd+1):
            frames.append((framesTemp[n-1] + framesTemp[n] + framesTemp[n+1])/3.0)
            
    else:
        frames = framesTemp
        
        
    


    ########################################################
    # mismatching correction between python numpy and ImageJ

    
    #N_xPixel = frames.frame_shape[0]
    #N_yPixel = frames.frame_shape[1]

    #N_xPixel = len(frames[0])
    #N_yPixel = len(frames[0][0])
    
    
    x_ImageJ, y_ImageJ = np.array(xpos), np.array(ypos)
    
    #xc = y_ImageJ
    #yc = N_xPixel - x_ImageJ +1 
    
    xc = x_ImageJ - 1 # for trackpy 3.0
    yc = y_ImageJ # for trackpy 3.0
            
    
    ########################################################
    
    
    

    print '\n\n## fast FIONA only ##'
    print 'frameStart ', frameStart
    print 'frameEnd ', frameEnd
    print 'max Ecc = ', maxEcc
    print 'max Pencentage err = ', maxPer
    print 'max Pixel err = ', maxPx
    print 'additionalFramesN = ', additionalFramesN
 
    print 'data file = ', filepath
    




    
    
    
    
    TotalNframes = len(frames)
    print 'TotalNframes: ', TotalNframes
    
    if frameEnd >= TotalNframes:
        frameEnd = TotalNframes - 1
        print '\nframeEnd changed to ', TotalNframes - 1, '\n'
    
    



    NCenterSubPixel = 7  # number of pixels from the center
    
    
    
    #frameTemp0 = frame0[yc-NCenterSubPixel:yc+NCenterSubPixel + 1, xc-NCenterSubPixel:xc+NCenterSubPixel+1]
    x3 = np.linspace(0, 2 * NCenterSubPixel, 2* NCenterSubPixel + 1)
    y3 = np.linspace(0, 2 * NCenterSubPixel, 2*NCenterSubPixel+1)
    x3, y3 = np.meshgrid(x3, y3)
    
    



    centerImage = [[] for _ in xrange(len(xc))] # centerImage[Molecule#][Frame#]
    
    for frameTemp in np.arange(frameStart, frameStart+4): # from the 0th frame to the end frame.
        print 'calculating center image data for frame # ', frameTemp
        for k in range(len(xpos)):            

            frameCenterImageTemp = frames[frameTemp][yc[k]-7:yc[k]+8, xc[k]-7:xc[k]+8]
            centerImage[k].append(frameCenterImageTemp) # it's not matching to the actual frame #, only first a few frames are saved
        
 
          
          
          
            
    xposfit_FIONA_0 = [[] for _ in xrange(len(xc))] # with S/B
    yposfit_FIONA_0 = [[] for _ in xrange(len(xc))] # with S/B
    
    centerIntensityMax5Ave = [[] for _ in xrange(len(xc))] 
    
    FIONA_frameN = []


    
    for f1n in range(frameStart, frameEnd + 1):
        
        #################################################################################
        ###  FIONA only part                              ############

        FIONA_frameN.append(f1n)
        
        for MN_fiona in range(len(xpos)): #for each feature
        
        
            print 'frame # ', f1n, ',   FIONA for M# ', MN_fiona 
            
            
            
  
            frameFIONA = frames[f1n][yc[MN_fiona]-NCenterSubPixel:yc[MN_fiona]+NCenterSubPixel+1, xc[MN_fiona]-NCenterSubPixel:xc[MN_fiona]+NCenterSubPixel+1]


            frameFIONA_raveled = frameFIONA.ravel()
            
            
            tempMax5Ave = int(np.mean( np.sort(frameFIONA_raveled)[-5:])) 
            tempMin5Ave = float(np.mean( np.sort(frameFIONA_raveled)[:5])) 
            
            #tempSBratio = tempMax5Ave/tempMin5Ave

            centerIntensityMax5Ave[MN_fiona].append(tempMax5Ave)
            
            #print 'tempMax5Ave, tempMin5Ave, tempSBratio = ', tempMax5Ave, tempMin5Ave, tempSBratio
        
            
            #if tempSBratio < 2.0:
            #    print 'skipped\n'
            #    continue
            
                    
            # FIONA fit the data  for all frames
            initial_guess = (tempMax5Ave, NCenterSubPixel+1,  NCenterSubPixel+1, 3, 3, 10, tempMin5Ave)
 


            try:
                popt, pcov = opt.curve_fit(twoD_Gaussian, (x3, y3), frameFIONA_raveled, p0 = initial_guess)
                #print '\nFIONA popt: ', popt
                #print '\nFIONA pcov: \n', pcov

                x0fit = popt[1]
                y0fit = popt[2]
                
                
                
                x0fitImageJ = xpos[MN_fiona] + (NCenterSubPixel - y0fit)  # to correct for imageJ x position
                y0fitImageJ = ypos[MN_fiona] - (NCenterSubPixel - x0fit)
                
                
                xposfit_FIONA_0[MN_fiona].append(np.around(x0fitImageJ,2))
                yposfit_FIONA_0[MN_fiona].append(np.around(y0fitImageJ,2))
                


             
                ####################################################
            
            except:
                   
                xposfit_FIONA_0[MN_fiona].append(np.nan)
                yposfit_FIONA_0[MN_fiona].append(np.nan)
                
                print 'FIONA fitting failed'

                #print 'M#, datatemp FIONA: ' , MN_fiona , '   ' ,  datatempFIONA

    return (centerImage, xposfit_FIONA_0, yposfit_FIONA_0, FIONA_frameN, centerIntensityMax5Ave)


















  


    
def FIONAonlyByEachFrame(filepath, xpos, ypos, frameStart, frameEnd, SvsB_Tolerance, maxEcc, maxPer, maxPx, additionalFramesN):   
    # using all features for all frames by analyzing frame by frame


    
    print '\n\n## FIONA only All Frames By Each Frame ##'
    print 'frameStart ', frameStart
    print 'frameEnd ', frameEnd
    print 'max Ecc = ', maxEcc
    print 'max Pencentage err = ', maxPer
    print 'max Pixel err = ', maxPx
    print 'additionalFramesN = ', additionalFramesN
 
    print 'data file = ', filepath
    



    
    framesTemp = tp.TiffStack(filepath)
    
    if additionalFramesN >= 1:
            
        frames = []
        frames.append( (framesTemp[0] + framesTemp[1])/2.0)
        for n in range(1, frameEnd+1):
            frames.append((framesTemp[n-1] + framesTemp[n] + framesTemp[n+1])/3.0)
            
    else:
        frames = framesTemp
        
        
    
    ########################################################
    # mismatching correction between python numpy and ImageJ

    
    #N_xPixel = frames.frame_shape[0]
    #N_yPixel = frames.frame_shape[1]

    #N_xPixel = len(frames[0])
    #N_yPixel = len(frames[0][0])
    

    x_ImageJ, y_ImageJ = np.array(xpos), np.array(ypos)
    
    #xc = y_ImageJ
    #yc = N_xPixel - x_ImageJ +1 
        
    xc = x_ImageJ - 1 # for trackpy 3.0
    yc = y_ImageJ # for trackpy 3.0

    ########################################################


    
    
    
    
    
    
    TotalNframes = len(frames)
    print 'TotalNframes: ', TotalNframes
    
    if frameEnd >= TotalNframes:
        frameEnd = TotalNframes - 1
        print '\nframeEnd changed to ', TotalNframes - 1, '\n'
    
    
    #NStartFrame = 0
    #frame0 = frames[NStartFrame]
    
    x2 = np.linspace(0, 511, 512)
    y2 = np.linspace(0, 511, 512)
    x2, y2 = np.meshgrid(x2, y2)
    
    
    NCenterSubPixel = 7  # number of pixels from the center
    
    
    
    #frameTemp0 = frame0[yc-NCenterSubPixel:yc+NCenterSubPixel + 1, xc-NCenterSubPixel:xc+NCenterSubPixel+1]
    x3 = np.linspace(0, 2 * NCenterSubPixel, 2* NCenterSubPixel + 1)
    y3 = np.linspace(0, 2 * NCenterSubPixel, 2*NCenterSubPixel+1)
    x3, y3 = np.meshgrid(x3, y3)
    
    
    BG_IntensityAve = [[] for _ in xrange(len(xc))] 
    centerIntensityAve = [[] for _ in xrange(len(xc))] 
    centerIntensityMax = [[] for _ in xrange(len(xc))] 
    centerImage = [[] for _ in xrange(len(xc))] # centerImage[Molecule#][Frame#]
    
    for frameTemp in np.arange(0, frameEnd+1): # from the 0th frame to the end frame.
        print 'calculating center intensity data for frame # ', frameTemp
        for k in range(len(xpos)):            
            
            frameLargeTemp = frames[frameTemp][yc[k]-4:yc[k]+5, xc[k]-4:xc[k]+5]
            frameCenterTemp = frames[frameTemp][yc[k]-2:yc[k]+3, xc[k]-2:xc[k]+3]
            #FiveMaxAveCenterTemp = int(np.mean( np.sort(frameCenterTemp.ravel())[-5:]) )
            
            if (frameTemp >= frameStart) and (frameTemp <= frameStart + 4):
                frameCenterImageTemp = frames[frameTemp][yc[k]-7:yc[k]+8, xc[k]-7:xc[k]+8]
                centerImage[k].append(frameCenterImageTemp) # it's not matching to the actual frame #, only first a few frames are saved
            
            centerIntensityAve[k].append(np.sum(frameCenterTemp)/25)
            centerIntensityMax[k].append(np.amax(frameCenterTemp))
            BG_IntensityAve[k].append( (np.sum(frameLargeTemp) - np.sum(frameCenterTemp))/56 )
            
    SBratioAve = np.float16(np.array(centerIntensityAve)) / np.array(BG_IntensityAve)    
    SBratioMax = np.float16(np.array(centerIntensityMax)) / np.array(BG_IntensityAve)    
    print 'SBratioAve = ', SBratioAve   
    print 'SBratioMax = ', SBratioMax
            
    
            
    
    #print 'centerIntensityAve = ', centerIntensityAve
    #print 'centerIntensityAve len = ', len(centerIntensityAve)
    
    BG_AveTotal = np.average(BG_IntensityAve)
    BG_AveEachFrame = np.around(np.average(BG_IntensityAve, axis = 0), 0) 
    print 'BG_AveTotal = ', BG_AveTotal
    print 'BG_AveEachFrame = ', BG_AveEachFrame
    #print 'frameCenterTemp ', frameCenterTemp.ravel()


    
    xposfit_FIONA_0 = [[] for _ in xrange(len(xc))] # with S/B
    yposfit_FIONA_0 = [[] for _ in xrange(len(xc))] # with S/B
    
    xposfit_FIONA_1 = [[] for _ in xrange(len(xc))] # with S/B and max e
    yposfit_FIONA_1 = [[] for _ in xrange(len(xc))] # with S/B and max e
    
    xposfit_FIONA_2 = [[] for _ in xrange(len(xc))] # with S/B, max e, max AfitErr/Afit %, max px
    yposfit_FIONA_2 = [[] for _ in xrange(len(xc))] # with S/B, max e, max AfitErr/Afit %, max px
        
    
    
    eccentricityFor_FIONA = [[] for _ in xrange(len(xc))] 
        




    
    data_FIONA = [[] for _ in xrange(len(xc))] 
    
    
    for n in range(len(xc)):
        data_FIONA[n].append('M# ' + str(n) + ' (x,y): ' + str(x_ImageJ[n]) + ' ' + str(y_ImageJ[n])  + ' #AveFrames: ' + str(additionalFramesN)      )
        data_FIONA[n].append('frame#,  x0fitFionaImageJ, y0fitFionaImageJ, Afit, z0fit, sigma_x, sigma_y, theta, eccentricity, BG_IntensityAve, centerIntensityAve, centerIntensityMax, SBratioAve, SBratioMax ::: AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaErr') 
    
    
    
    for f1n in range(0, frameEnd + 1):
        
        #################################################################################
        ###  FIONA only part                              ############
        print '\n FIONA start for frame # ', f1n
        for MN_fiona in range(len(xpos)): #for each feature
        
            if SBratioAve[MN_fiona][f1n] <= SvsB_Tolerance :
                print 'Low SB, FIONA passSed for M# ', MN_fiona, '  frame # ', f1n
                
                pass

            else:           
                print '\nFIONA for M# ', MN_fiona, '  frame # ', f1n
                frameFIONA = frames[f1n][yc[MN_fiona]-NCenterSubPixel:yc[MN_fiona]+NCenterSubPixel+1, xc[MN_fiona]-NCenterSubPixel:xc[MN_fiona]+NCenterSubPixel+1]


                frameFIONA_raveled = frameFIONA.ravel()
            
            
                # FIONA fit the data  for all frames
                initial_guess = (6000, NCenterSubPixel+1,  NCenterSubPixel+1, 3, 3, 10, 2000)
 
                k = 1
               
                Afit = 0
                x0fit = 0
                y0fit = 0
                z0fit = 0
                
                try:
                    popt, pcov = opt.curve_fit(twoD_Gaussian, (x3, y3), frameFIONA_raveled, p0 = initial_guess)
                    #print '\nFIONA popt: ', popt
                    #print '\nFIONA pcov: \n', pcov
                    Afit = popt[0]
                    x0fit = popt[1]
                    y0fit = popt[2]
                    
                    print 'x0fit = ',x0fit
                    
                    x0fitImageJ = xpos[MN_fiona] + (NCenterSubPixel - y0fit)  # to correct for imageJ x position
                    y0fitImageJ = ypos[MN_fiona] - (NCenterSubPixel - x0fit)
                    
                    z0fit = popt[6]
                    
                    sigma_x = popt[3]
                    sigma_y = popt[4]
                    theta = popt[5]

                    
                    #FIONA AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaErr
                    AfitErr = np.sqrt(pcov[0,0])                    
                    x0fitImageJErr = np.sqrt(pcov[1,1])
                    y0fitImageJErr = np.sqrt(pcov[2,2]) 
                    z0fitErr = np.sqrt(pcov[6,6]) 
                    sigma_xErr = np.sqrt(pcov[3,3]) 
                    sigma_yErr = np.sqrt(pcov[4,4]) 
                    thetaErr = np.sqrt(pcov[5,5])
                                         

                    
                    # FIONA Percentage of the errors: AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaErr
                    AfitErrP = np.around(100.0 * AfitErr / Afit, 0)                    
                    x0fitImageJErrPx = np.around(np.sqrt(pcov[1,1]),2)
                    y0fitImageJErrPx = np.around(np.sqrt(pcov[2,2]),2) 
                    #z0fitErrP = np.around(100.0 * z0fitErr / z0fit, 0)      
                    #sigma_xErrP = np.around(100.0 * sigma_xErr/sigma_x, 0)   
                    #sigma_yErrP = np.around(100.0 * sigma_yErr/sigma_y, 0)  
                    #thetaErrDeg = np.around(np.sqrt(pcov[5,5]),2)
                    
                    
                    
                    
                    ########### eccentricicity calculation  ############
                    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
                    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
                    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2) 
                    
                    #print 'a, b, c =', a,b,c
                    
                    
                    A = a
                    B = 2*b
                    C = c
                            
                    eccentricityTemp = np.sqrt(2*(np.sqrt((A-C)**2 + B**2)) / (A+C + np.sqrt(((A-C)**2 + B**2)))  )
                    eccentricityFor_FIONA[MN_fiona].append(eccentricityTemp)
                    #print 'eccentricityTemp = ', eccentricityTemp
                    
                    #'frame#,  x0fitFionaImageJ, y0fitFionaImageJ, Afit, z0fit, sigma_x, sigma_y, theta, eccentricity
                    #, BG_IntensityAve, centerIntensityAve, centerIntensityMax, SBratioAve, SBratioMax)'
                    datatempFIONA = (str(f1n) + ' ' + str(x0fitImageJ) + ' ' + str(y0fitImageJ) + ' ' 
                            + str(Afit) +  ' ' + str(z0fit) + ' ' + str(sigma_x) + ' ' + str(sigma_y) + ' ' +  str(theta) + ' ' + str(eccentricityTemp) + ' '
                            + str(BG_IntensityAve[MN_fiona][f1n]) + ' ' + str(centerIntensityAve[MN_fiona][f1n]) + ' ' 
                            + str(centerIntensityMax[MN_fiona][f1n]) + ' ' + str(SBratioAve[MN_fiona][f1n]) + ' ' + str(SBratioMax[MN_fiona][f1n]) + ' Err: '
                            + str(AfitErr) +  ' ' + str(x0fitImageJErr) + ' ' + str(y0fitImageJErr) + ' ' + str(z0fitErr) + ' ' 
                            + str(sigma_xErr) + ' ' + str(sigma_yErr) + ' ' +  str(thetaErr)
                            )
                    
                    
                    
                    data_FIONA[MN_fiona].append(datatempFIONA)

                    
                    #print 'M#, datatemp FIONA: ' , MN_fiona , '   ' ,  datatempFIONA
                    
                    xposfit_FIONA_0[MN_fiona].append(np.around(x0fitImageJ,2))
                    yposfit_FIONA_0[MN_fiona].append(np.around(y0fitImageJ,2))
                    
                    if (eccentricityTemp <= maxEcc):
                        xposfit_FIONA_1[MN_fiona].append(np.around(x0fitImageJ,2))
                        yposfit_FIONA_1[MN_fiona].append(np.around(y0fitImageJ,2))
                        print '\nFIONA_1 OK'
                        

                        if (AfitErrP <= maxPer) & (x0fitImageJErrPx <= maxPx) & (y0fitImageJErrPx <= maxPx):
                            xposfit_FIONA_2[MN_fiona].append(np.around(x0fitImageJ,2))
                            yposfit_FIONA_2[MN_fiona].append(np.around(y0fitImageJ,2))
                            print '\nFIONA_2 OK'
                            
                            
                        else:
                            print 'FINOA_2 failed'
                            print 'AfitErrP ', AfitErrP
                            print 'x0fitImageJErrPx ', x0fitImageJErrPx
                            print 'y0fitImageJErrPx ', y0fitImageJErrPx
                            
                    else:
                        print 'FIONA_1 failed'
                        print 'eccentricityTemp ', eccentricityTemp
                        

                

                 
                    ####################################################
                
                except:
                    print 'FIONA fitting failed'
                    eccentricityTemp = 1.2
                    datatempFIONA = (str(f1n) + ' ' + str(999) + ' ' + str(999) + ' ' 
                            + str(999) +  ' ' + str(999) + ' ' + str(999) + ' ' + str(999) + ' ' +  str(999) + ' ' + str(eccentricityTemp) + ' '
                            + str(BG_IntensityAve[MN_fiona][f1n]) + ' ' + str(centerIntensityAve[MN_fiona][f1n]) + ' ' 
                            + str(centerIntensityMax[MN_fiona][f1n]) + ' ' + str(SBratioAve[MN_fiona][f1n]) + ' ' + str(SBratioMax[MN_fiona][f1n]) + ' Err: '       
                            + str(999) +  ' ' + str(999) + ' ' + str(999) + ' ' + str(999) + ' ' 
                            + str(999) + ' ' + str(999) + ' ' +  str(999)
                            )
                     
                    data_FIONA[MN_fiona].append(datatempFIONA)
                    eccentricityFor_FIONA[MN_fiona].append(eccentricityTemp)

                    
                    #print 'M#, datatemp FIONA: ' , MN_fiona , '   ' ,  datatempFIONA
                    
        
        print 'FIONA end\n'
        ######################  FIONA only part                              ############
        #################################################################################
      
        
                        
    

    


    return (centerImage, xposfit_FIONA_0, yposfit_FIONA_0, data_FIONA, xposfit_FIONA_1, yposfit_FIONA_1, xposfit_FIONA_2, yposfit_FIONA_2
            , eccentricityFor_FIONA, centerIntensityMax)











  


'''
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
'''



    
def FIONA_drift(filepath, xpos, ypos, frameStart, frameEnd, SvsB_Tolerance, maxEcc, maxPer, maxPx, additionalFramesN):   
    # using all features for all frames by analyzing frame by frame


    
    framesTemp = tp.TiffStack(filepath)
    
    if additionalFramesN == 1:
            
        frames = []
        frames.append( (framesTemp[0] + framesTemp[1])/2.0)
        for n in range(1, frameEnd+1):
            print 'summing frames n-1, n, n+1: ', n
            frames.append((framesTemp[n-1] + framesTemp[n] + framesTemp[n+1])/3.0)
            
    elif additionalFramesN == 0: 
        frames = framesTemp
        
    else:
        print 'error in additional adjacent frame numbers '
        return
        
        
    ########################################################
    # mismatching correction between python numpy and ImageJ


    #N_xPixel = len(frames[0])
    
    x_ImageJ, y_ImageJ = np.array(xpos), np.array(ypos)
    
    #xc = y_ImageJ
    #yc = N_xPixel - x_ImageJ +1 
        
    xc = x_ImageJ - 1 # for trackpy 3.0
    yc = y_ImageJ # for trackpy 3.0

    ########################################################
    

    
    
    
    
    
    TotalNframes = len(frames)
    print 'TotalNframes: ', TotalNframes
    
    if frameEnd >= TotalNframes:
        frameEnd = TotalNframes - 1
        print '\nframeEnd changed to ', TotalNframes - 1, '\n'
    
    
    #NStartFrame = 0
    #frame0 = frames[NStartFrame]
    
    x2 = np.linspace(0, 511, 512)
    y2 = np.linspace(0, 511, 512)
    x2, y2 = np.meshgrid(x2, y2)
    
    
    NCenterSubPixel = 7  # number of pixels from the center
    
    
    
    #frameTemp0 = frame0[yc-NCenterSubPixel:yc+NCenterSubPixel + 1, xc-NCenterSubPixel:xc+NCenterSubPixel+1]
    x3 = np.linspace(0, 2 * NCenterSubPixel, 2* NCenterSubPixel + 1)
    y3 = np.linspace(0, 2 * NCenterSubPixel, 2*NCenterSubPixel+1)
    x3, y3 = np.meshgrid(x3, y3)
    
    
    BG_IntensityAve = [[] for _ in xrange(len(xc))] 
    centerIntensityAve = [[] for _ in xrange(len(xc))] 
    centerIntensityMax = [[] for _ in xrange(len(xc))] 
    centerImage = [[] for _ in xrange(len(xc))] # centerImage[Molecule#][Frame#]
    
    for frameTemp in np.arange(0, frameEnd+1): # from the 0th frame to the end frame.
        print 'calculating center intensity data for frame # ', frameTemp
        for k in range(len(xpos)):            
            
            frameLargeTemp = frames[frameTemp][yc[k]-4:yc[k]+5, xc[k]-4:xc[k]+5]
            frameCenterTemp = frames[frameTemp][yc[k]-2:yc[k]+3, xc[k]-2:xc[k]+3]
            #FiveMaxAveCenterTemp = int(np.mean( np.sort(frameCenterTemp.ravel())[-5:]) )
            
            if (frameTemp >= frameStart) and (frameTemp <= frameStart + 4):
                frameCenterImageTemp = frames[frameTemp][yc[k]-7:yc[k]+8, xc[k]-7:xc[k]+8]
                centerImage[k].append(frameCenterImageTemp) # it's not matching to the actual frame #, only first a few frames are saved
            
            centerIntensityAve[k].append(np.sum(frameCenterTemp)/25)
            centerIntensityMax[k].append(np.amax(frameCenterTemp))
            BG_IntensityAve[k].append( (np.sum(frameLargeTemp) - np.sum(frameCenterTemp))/56 )
            
    SBratioAve = np.float16(np.array(centerIntensityAve)) / np.array(BG_IntensityAve)    
    SBratioMax = np.float16(np.array(centerIntensityMax)) / np.array(BG_IntensityAve)    
    print 'SBratioAve = ', SBratioAve   
    print 'SBratioMax = ', SBratioMax
            
    
            
    
    #print 'centerIntensityAve = ', centerIntensityAve
    #print 'centerIntensityAve len = ', len(centerIntensityAve)
    
    BG_AveTotal = np.average(BG_IntensityAve)
    BG_AveEachFrame = np.around(np.average(BG_IntensityAve, axis = 0), 0) 
    print 'BG_AveTotal = ', BG_AveTotal
    print 'BG_AveEachFrame = ', BG_AveEachFrame
    #print 'frameCenterTemp ', frameCenterTemp.ravel()


    
    xposfit_FIONA_0 = [[] for _ in xrange(len(xc))] # with S/B
    yposfit_FIONA_0 = [[] for _ in xrange(len(xc))] # with S/B
    
    xposfit_FIONA_1 = [[] for _ in xrange(len(xc))] # with S/B and max e
    yposfit_FIONA_1 = [[] for _ in xrange(len(xc))] # with S/B and max e
    
    xposfit_FIONA_2 = [[] for _ in xrange(len(xc))] # with S/B, max e, max AfitErr/Afit %, max px
    yposfit_FIONA_2 = [[] for _ in xrange(len(xc))] # with S/B, max e, max AfitErr/Afit %, max px
        

    xposfit_dx_FIONA_2 = [[] for _ in xrange(len(xc))] # with S/B, max e, max AfitErr/Afit %, max px
    yposfit_dy_FIONA_2 = [[] for _ in xrange(len(xc))] # with S/B, max e, max AfitErr/Afit %, max px
        


    FrameN_FIONA_2 = [[] for _ in xrange(len(xc))] # 
    
    FrameN_FIONA_2_drift = []
    fx2_dx_ave = []
    fy2_dy_ave = []
    
    eccentricityFor_FIONA = [[] for _ in xrange(len(xc))] 
        




    
    data_FIONA = [[] for _ in xrange(len(xc))] 
    
    for n in range(len(xc)):
        data_FIONA[n].append('M# ' + str(n) + ' (x,y): ' + str(x_ImageJ[n]) + ' ' + str(y_ImageJ[n])  + ' AdjacentFrame#s: ' + str(additionalFramesN)      )
        data_FIONA[n].append('frame#,  x0fitFionaImageJ, y0fitFionaImageJ, Afit, z0fit, sigma_x, sigma_y, theta, eccentricity, BG_IntensityAve, centerIntensityAve, centerIntensityMax, SBratioAve, SBratioMax ::: AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaErr') 
    
    
    
    for f1n in range(0, frameEnd + 1):
        
        #################################################################################
        ###  FIONA only part                              ############
        print '\n FIONA start for frame # ', f1n
        dx_for_ave = []
        dy_for_ave = []
        for MN_fiona in range(len(xpos)): #for each feature
        
            if SBratioAve[MN_fiona][f1n] <= SvsB_Tolerance :
                print 'Low SB, FIONA passSed for M# ', MN_fiona, '  frame # ', f1n
                
                pass

            else:           
                print '\nFIONA for M# ', MN_fiona, '  frame # ', f1n
                frameFIONA = frames[f1n][yc[MN_fiona]-NCenterSubPixel:yc[MN_fiona]+NCenterSubPixel+1, xc[MN_fiona]-NCenterSubPixel:xc[MN_fiona]+NCenterSubPixel+1]


                frameFIONA_raveled = frameFIONA.ravel()
                
                                
                tempMax5Ave = int(np.mean( np.sort(frameFIONA_raveled)[-5:])) 
                tempMin5Ave = float(np.mean( np.sort(frameFIONA_raveled)[:5])) 
                
            
            
                # FIONA fit the data  for all frames
            
            
                #initial_guess = (6000, NCenterSubPixel+1,  NCenterSubPixel+1, 3, 3, 10, 2000)
                
                initial_guess = (tempMax5Ave, NCenterSubPixel+1,  NCenterSubPixel+1, 3, 3, 10, tempMin5Ave)
 
                k = 1
               
                Afit = 0
                x0fit = 0
                y0fit = 0
                z0fit = 0
                
                try:
                    popt, pcov = opt.curve_fit(twoD_Gaussian, (x3, y3), frameFIONA_raveled, p0 = initial_guess)
                    #print '\nFIONA popt: ', popt
                    #print '\nFIONA pcov: \n', pcov
                    Afit = popt[0] 
                    x0fit = popt[1] 
                    y0fit = popt[2]
                    
                    print 'x0fit = ',x0fit
                    
                    x0fitImageJ = xpos[MN_fiona] + (NCenterSubPixel - y0fit)  # to correct for imageJ x position
                    y0fitImageJ = ypos[MN_fiona] - (NCenterSubPixel - x0fit)
                    
                    z0fit = popt[6]
                    
                    sigma_x = popt[3]
                    sigma_y = popt[4]
                    theta = popt[5]

                    
                    #FIONA AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaErr
                    AfitErr = np.sqrt(pcov[0,0])                    
                    x0fitImageJErr = np.sqrt(pcov[1,1])
                    y0fitImageJErr = np.sqrt(pcov[2,2]) 
                    z0fitErr = np.sqrt(pcov[6,6]) 
                    sigma_xErr = np.sqrt(pcov[3,3]) 
                    sigma_yErr = np.sqrt(pcov[4,4]) 
                    thetaErr = np.sqrt(pcov[5,5])
                                         

                    
                    # FIONA Percentage of the errors: AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaErr
                    AfitErrP = np.around(100.0 * AfitErr / Afit, 0)                    
                    x0fitImageJErrPx = np.around(np.sqrt(pcov[1,1]),2)
                    y0fitImageJErrPx = np.around(np.sqrt(pcov[2,2]),2) 
                    #z0fitErrP = np.around(100.0 * z0fitErr / z0fit, 0)      
                    #sigma_xErrP = np.around(100.0 * sigma_xErr/sigma_x, 0)   
                    #sigma_yErrP = np.around(100.0 * sigma_yErr/sigma_y, 0)  
                    #thetaErrDeg = np.around(np.sqrt(pcov[5,5]),2)
                    
                    
                    
                    
                    ########### eccentricicity calculation  ############
                    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
                    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
                    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2) 
                    
                    #print 'a, b, c =', a,b,c
                    
                    
                    A = a
                    B = 2*b
                    C = c
                            
                    eccentricityTemp = np.sqrt(2*(np.sqrt((A-C)**2 + B**2)) / (A+C + np.sqrt(((A-C)**2 + B**2)))  )
                    eccentricityFor_FIONA[MN_fiona].append(eccentricityTemp)
                    #print 'eccentricityTemp = ', eccentricityTemp
                    
                    #'frame#,  x0fitFionaImageJ, y0fitFionaImageJ, Afit, z0fit, sigma_x, sigma_y, theta, eccentricity
                    #, BG_IntensityAve, centerIntensityAve, centerIntensityMax, SBratioAve, SBratioMax)'
                    datatempFIONA = (str(f1n) + ' ' + str(x0fitImageJ) + ' ' + str(y0fitImageJ) + ' ' 
                            + str(Afit) +  ' ' + str(z0fit) + ' ' + str(sigma_x) + ' ' + str(sigma_y) + ' ' +  str(theta) + ' ' + str(eccentricityTemp) + ' '
                            + str(BG_IntensityAve[MN_fiona][f1n]) + ' ' + str(centerIntensityAve[MN_fiona][f1n]) + ' ' 
                            + str(centerIntensityMax[MN_fiona][f1n]) + ' ' + str(SBratioAve[MN_fiona][f1n]) + ' ' + str(SBratioMax[MN_fiona][f1n]) + ' Err: '
                            + str(AfitErr) +  ' ' + str(x0fitImageJErr) + ' ' + str(y0fitImageJErr) + ' ' + str(z0fitErr) + ' ' 
                            + str(sigma_xErr) + ' ' + str(sigma_yErr) + ' ' +  str(thetaErr)
                            )
                    
                    
                    
                    data_FIONA[MN_fiona].append(datatempFIONA)

                    
                    #print 'M#, datatemp FIONA: ' , MN_fiona , '   ' ,  datatempFIONA
                    
                    xposfit_FIONA_0[MN_fiona].append(np.around(x0fitImageJ,2))
                    yposfit_FIONA_0[MN_fiona].append(np.around(y0fitImageJ,2))
                    
                    if (eccentricityTemp <= maxEcc):
                        xposfit_FIONA_1[MN_fiona].append(np.around(x0fitImageJ,2))
                        yposfit_FIONA_1[MN_fiona].append(np.around(y0fitImageJ,2))
                        print '\nFIONA_1 OK'
                        

                        if (AfitErrP <= maxPer) & (x0fitImageJErrPx <= maxPx) & (y0fitImageJErrPx <= maxPx):
                            xposfit_FIONA_2[MN_fiona].append(np.around(x0fitImageJ,2))
                            yposfit_FIONA_2[MN_fiona].append(np.around(y0fitImageJ,2))
                            FrameN_FIONA_2[MN_fiona].append(f1n)


                            if len(FrameN_FIONA_2[MN_fiona]) == 0:
                                dxtemp = 0
                                dytemp = 0
                                
                            else:
                                dxtemp = x0fitImageJ - xposfit_FIONA_2[MN_fiona][0]
                                dytemp = y0fitImageJ - yposfit_FIONA_2[MN_fiona][0]
                                
                                dx_for_ave.append(dxtemp)
                                dy_for_ave.append(dytemp)




                            xposfit_dx_FIONA_2[MN_fiona].append(np.around(dxtemp,3))
                            yposfit_dy_FIONA_2[MN_fiona].append(np.around(dytemp,3))

                            
                            print '\nFIONA_2 OK'
                            
                            
                        else:
                            print 'FINOA_2 failed'
                            print 'AfitErrP ', AfitErrP
                            print 'x0fitImageJErrPx ', x0fitImageJErrPx
                            print 'y0fitImageJErrPx ', y0fitImageJErrPx
                            
                    else:
                        print 'FIONA_1 failed'
                        print 'eccentricityTemp ', eccentricityTemp
                        

                

                 
                    ####################################################
                
                except:
                    print 'FIONA fitting failed'
                    eccentricityTemp = 1.2
                    datatempFIONA = (str(f1n) + ' ' + str(999) + ' ' + str(999) + ' ' 
                            + str(999) +  ' ' + str(999) + ' ' + str(999) + ' ' + str(999) + ' ' +  str(999) + ' ' + str(eccentricityTemp) + ' '
                            + str(BG_IntensityAve[MN_fiona][f1n]) + ' ' + str(centerIntensityAve[MN_fiona][f1n]) + ' ' 
                            + str(centerIntensityMax[MN_fiona][f1n]) + ' ' + str(SBratioAve[MN_fiona][f1n]) + ' ' + str(SBratioMax[MN_fiona][f1n]) + ' Err: '       
                            + str(999) +  ' ' + str(999) + ' ' + str(999) + ' ' + str(999) + ' ' 
                            + str(999) + ' ' + str(999) + ' ' +  str(999)
                            )
                     
                    data_FIONA[MN_fiona].append(datatempFIONA)
                    eccentricityFor_FIONA[MN_fiona].append(eccentricityTemp)

                    
                    #print 'M#, datatemp FIONA: ' , MN_fiona , '   ' ,  datatempFIONA
                    
        FrameN_FIONA_2_drift.append(f1n)
        fx2_dx_ave.append(np.mean(dx_for_ave))
        fy2_dy_ave.append(np.mean(dy_for_ave))
        
        print 'FIONA end\n'
        ######################  FIONA only part                              ############
        #################################################################################
      
        
        
    return (FrameN_FIONA_2, xposfit_FIONA_2, yposfit_FIONA_2, data_FIONA, eccentricityFor_FIONA, xposfit_dx_FIONA_2, yposfit_dy_FIONA_2, FrameN_FIONA_2_drift, fx2_dx_ave, fy2_dy_ave)
            
        


    
def gshrimpAllFramesByEachFrame(filepath, xpos, ypos, frameStart, frameEnd, SvsB_Tolerance, maxEcc, maxPer, maxPx, additionalFramesN):   
    # using all features for all frames by analyzing frame by frame

    print '\n\n## gshrimp All Frames By Each Frame ##'
    print 'Processing gshrimpAllFramesByEachFrame'
    print 'frameStart ', frameStart
    print 'frameEnd ', frameEnd
    print 'max Ecc = ', maxEcc
    print 'max Pencentage err = ', maxPer
    print 'max Pixel err = ', maxPx
    print 'additionalFramesN = ', additionalFramesN
 
    print 'data file = ', filepath
    



    
    framesTemp = tp.TiffStack(filepath)
    
    if additionalFramesN == 1:
            
        frames = []
        frames.append( (framesTemp[0] + framesTemp[1])/2.0)
        for n in range(1, frameEnd+1):
            frames.append((framesTemp[n-1] + framesTemp[n] + framesTemp[n+1])/3.0)
            
    elif additionalFramesN == 0: 
        frames = framesTemp
        
    else:
        print 'error in additional adjacent frame numbers '
        return
        
        
        
    ########################################################
    # mismatching correction between python numpy and ImageJ

    
    #N_xPixel = frames.frame_shape[0]
    #N_yPixel = frames.frame_shape[1]

    #N_xPixel = len(frames[0])
    #N_yPixel = len(frames[0][0])
    
    
    #print 'len(frames[0]) = ', len(frames[0])
    #print 'frames[0] = ', frames[0]
    #print 'len(frames[0][0]) = ', len(frames[0][0])
    #print 'frames[0][0] = ', frames[0][0]
        
    
    x_ImageJ, y_ImageJ = np.array(xpos), np.array(ypos)
    
    #xc = y_ImageJ
    #yc = N_xPixel - x_ImageJ + 1
    
    xc = x_ImageJ - 1 # for trackpy 3.0
    yc = y_ImageJ # for trackpy 3.0
            
    ########################################################
    

    
    
    
    
    
    TotalNframes = len(frames)
    print 'TotalNframes: ', TotalNframes
    
    if frameEnd >= TotalNframes:
        frameEnd = TotalNframes - 1
        print '\nframeEnd changed to ', TotalNframes - 1, '\n'
    
    
    #NStartFrame = 0
    #frame0 = frames[NStartFrame]
    
    x2 = np.linspace(0, 511, 512)
    y2 = np.linspace(0, 511, 512)
    x2, y2 = np.meshgrid(x2, y2)
    
    
    NCenterSubPixel = 7  # number of pixels from the center
    
    
    
    #frameTemp0 = frame0[yc-NCenterSubPixel:yc+NCenterSubPixel + 1, xc-NCenterSubPixel:xc+NCenterSubPixel+1]
    x3 = np.linspace(0, 2 * NCenterSubPixel, 2* NCenterSubPixel + 1)
    y3 = np.linspace(0, 2 * NCenterSubPixel, 2*NCenterSubPixel+1)
    x3, y3 = np.meshgrid(x3, y3)
    
    
    BG_IntensityAve = [[] for _ in xrange(len(xc))] 
    centerIntensityAve = [[] for _ in xrange(len(xc))] 
    centerIntensityMax = [[] for _ in xrange(len(xc))] 
    centerImage = [[] for _ in xrange(len(xc))] # centerImage[Molecule#][Frame#]
    
    for frameTemp in np.arange(0, frameEnd+1): # from the 0th frame to the end frame.
        print 'calculating center intensity data for frame # ', frameTemp
        for k in range(len(xpos)):            
            
            frameLargeTemp = frames[frameTemp][yc[k]-4:yc[k]+5, xc[k]-4:xc[k]+5]
            frameCenterTemp = frames[frameTemp][yc[k]-2:yc[k]+3, xc[k]-2:xc[k]+3]
            #FiveMaxAveCenterTemp = int(np.mean( np.sort(frameCenterTemp.ravel())[-5:]) )
            
            if (frameTemp >= frameStart) and (frameTemp <= frameStart + 4):
                frameCenterImageTemp = frames[frameTemp][yc[k]-7:yc[k]+8, xc[k]-7:xc[k]+8]
                centerImage[k].append(frameCenterImageTemp) # it's not matching to the actual frame #, only first a few frames are saved
            
            centerIntensityAve[k].append(np.sum(frameCenterTemp)/25)
            centerIntensityMax[k].append(np.amax(frameCenterTemp))
            BG_IntensityAve[k].append( (np.sum(frameLargeTemp) - np.sum(frameCenterTemp))/56 )
            
    SBratioAve = np.float16(np.array(centerIntensityAve)) / np.array(BG_IntensityAve)    
    SBratioMax = np.float16(np.array(centerIntensityMax)) / np.array(BG_IntensityAve)    
    print 'SBratioAve = ', SBratioAve   
    print 'SBratioMax = ', SBratioMax
            
    
            
    
    #print 'centerIntensityAve = ', centerIntensityAve
    #print 'centerIntensityAve len = ', len(centerIntensityAve)
    
    BG_AveTotal = np.average(BG_IntensityAve)
    BG_AveEachFrame = np.around(np.average(BG_IntensityAve, axis = 0), 0) 
    print 'BG_AveTotal = ', BG_AveTotal
    print 'BG_AveEachFrame = ', BG_AveEachFrame
    #print 'frameCenterTemp ', frameCenterTemp.ravel()


    
    xposfit_FIONA_0 = [[] for _ in xrange(len(xc))] # with S/B
    yposfit_FIONA_0 = [[] for _ in xrange(len(xc))] # with S/B
    
    xposfit_FIONA_1 = [[] for _ in xrange(len(xc))] # with S/B and max e
    yposfit_FIONA_1 = [[] for _ in xrange(len(xc))] # with S/B and max e
    
    xposfit_FIONA_2 = [[] for _ in xrange(len(xc))] # with S/B, max e, max AfitErr/Afit %, max px
    yposfit_FIONA_2 = [[] for _ in xrange(len(xc))] # with S/B, max e, max AfitErr/Afit %, max px
        
    
    eccentricityFor_FIONA = [[] for _ in xrange(len(xc))] 
    eccentricityFor_gshrimp = [[] for _ in xrange(len(xc))] 
    
    
    xposfit_gshrimp0 = [[] for _ in xrange(len(xc))] 
    yposfit_gshrimp0 = [[] for _ in xrange(len(xc))] 
    
    xposfit_gshrimp1 = [[] for _ in xrange(len(xc))] # with S/B and max e
    yposfit_gshrimp1 = [[] for _ in xrange(len(xc))] # with S/B and max e
    
    xposfit_gshrimp2 = [[] for _ in xrange(len(xc))] # with S/B, max e, max AfitErr/Afit %, max px
    yposfit_gshrimp2 = [[] for _ in xrange(len(xc))] # with S/B, max e, max AfitErr/Afit %, max px
    
        




    
    data_FIONA = [[] for _ in xrange(len(xc))] 
    data_gshrimp = [[] for _ in xrange(len(xc))] 
    
    for n in range(len(xc)):
        data_gshrimp[n].append('M# ' + str(n) + ' (x,y): ' + str(x_ImageJ[n]) + ' ' + str(y_ImageJ[n]) + ' #AveFrames: ' + str(additionalFramesN)      )
        data_gshrimp[n].append('f1n f2n centerIntensityMax[f1n] centerIntensityMax[f2n] Afit x0fitImageJ y0fitImageJ z0fit sigma_x sigma_y theta eccentricity ::: AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaErr') 
        
        data_FIONA[n].append('M# ' + str(n) + ' (x,y): ' + str(x_ImageJ[n]) + ' ' + str(y_ImageJ[n])  + ' AdjacentFrame#s: ' + str(additionalFramesN)      )
        data_FIONA[n].append('frame#,  x0fitFionaImageJ, y0fitFionaImageJ, Afit, z0fit, sigma_x, sigma_y, theta, eccentricity, BG_IntensityAve, centerIntensityAve, centerIntensityMax, SBratioAve, SBratioMax ::: AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaErr') 
    
    
    
    for f1n in range(0, frameEnd + 1):
        
        #################################################################################
        ###  FIONA only part                              ############
        print '\n FIONA start for frame # ', f1n
        for MN_fiona in range(len(xpos)): #for each feature
        
            if SBratioAve[MN_fiona][f1n] <= SvsB_Tolerance :
                print 'Low SB, FIONA passSed for M# ', MN_fiona, '  frame # ', f1n
                
                pass

            else:           
                print '\nFIONA for M# ', MN_fiona, '  frame # ', f1n
                frameFIONA = frames[f1n][yc[MN_fiona]-NCenterSubPixel:yc[MN_fiona]+NCenterSubPixel+1, xc[MN_fiona]-NCenterSubPixel:xc[MN_fiona]+NCenterSubPixel+1]


                frameFIONA_raveled = frameFIONA.ravel()
                
                                
                tempMax5Ave = int(np.mean( np.sort(frameFIONA_raveled)[-5:])) 
                tempMin5Ave = float(np.mean( np.sort(frameFIONA_raveled)[:5])) 
                
            
            
                # FIONA fit the data  for all frames
            
            
                #initial_guess = (6000, NCenterSubPixel+1,  NCenterSubPixel+1, 3, 3, 10, 2000)
                
                initial_guess = (tempMax5Ave, NCenterSubPixel+1,  NCenterSubPixel+1, 3, 3, 10, tempMin5Ave)
 
                k = 1
               
                Afit = 0
                x0fit = 0
                y0fit = 0
                z0fit = 0
                
                try:
                    popt, pcov = opt.curve_fit(twoD_Gaussian, (x3, y3), frameFIONA_raveled, p0 = initial_guess)
                    #print '\nFIONA popt: ', popt
                    #print '\nFIONA pcov: \n', pcov
                    Afit = popt[0] 
                    x0fit = popt[1] 
                    y0fit = popt[2]
                    
                    print 'x0fit = ',x0fit
                    
                    x0fitImageJ = xpos[MN_fiona] + (NCenterSubPixel - y0fit)  # to correct for imageJ x position
                    y0fitImageJ = ypos[MN_fiona] - (NCenterSubPixel - x0fit)
                    
                    z0fit = popt[6]
                    
                    sigma_x = popt[3]
                    sigma_y = popt[4]
                    theta = popt[5]

                    
                    #FIONA AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaErr
                    AfitErr = np.sqrt(pcov[0,0])                    
                    x0fitImageJErr = np.sqrt(pcov[1,1])
                    y0fitImageJErr = np.sqrt(pcov[2,2]) 
                    z0fitErr = np.sqrt(pcov[6,6]) 
                    sigma_xErr = np.sqrt(pcov[3,3]) 
                    sigma_yErr = np.sqrt(pcov[4,4]) 
                    thetaErr = np.sqrt(pcov[5,5])
                                         

                    
                    # FIONA Percentage of the errors: AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaErr
                    AfitErrP = np.around(100.0 * AfitErr / Afit, 0)                    
                    x0fitImageJErrPx = np.around(np.sqrt(pcov[1,1]),2)
                    y0fitImageJErrPx = np.around(np.sqrt(pcov[2,2]),2) 
                    #z0fitErrP = np.around(100.0 * z0fitErr / z0fit, 0)      
                    #sigma_xErrP = np.around(100.0 * sigma_xErr/sigma_x, 0)   
                    #sigma_yErrP = np.around(100.0 * sigma_yErr/sigma_y, 0)  
                    #thetaErrDeg = np.around(np.sqrt(pcov[5,5]),2)
                    
                    
                    
                    
                    ########### eccentricicity calculation  ############
                    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
                    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
                    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2) 
                    
                    #print 'a, b, c =', a,b,c
                    
                    
                    A = a
                    B = 2*b
                    C = c
                            
                    eccentricityTemp = np.sqrt(2*(np.sqrt((A-C)**2 + B**2)) / (A+C + np.sqrt(((A-C)**2 + B**2)))  )
                    eccentricityFor_FIONA[MN_fiona].append(eccentricityTemp)
                    #print 'eccentricityTemp = ', eccentricityTemp
                    
                    #'frame#,  x0fitFionaImageJ, y0fitFionaImageJ, Afit, z0fit, sigma_x, sigma_y, theta, eccentricity
                    #, BG_IntensityAve, centerIntensityAve, centerIntensityMax, SBratioAve, SBratioMax)'
                    datatempFIONA = (str(f1n) + ' ' + str(x0fitImageJ) + ' ' + str(y0fitImageJ) + ' ' 
                            + str(Afit) +  ' ' + str(z0fit) + ' ' + str(sigma_x) + ' ' + str(sigma_y) + ' ' +  str(theta) + ' ' + str(eccentricityTemp) + ' '
                            + str(BG_IntensityAve[MN_fiona][f1n]) + ' ' + str(centerIntensityAve[MN_fiona][f1n]) + ' ' 
                            + str(centerIntensityMax[MN_fiona][f1n]) + ' ' + str(SBratioAve[MN_fiona][f1n]) + ' ' + str(SBratioMax[MN_fiona][f1n]) + ' Err: '
                            + str(AfitErr) +  ' ' + str(x0fitImageJErr) + ' ' + str(y0fitImageJErr) + ' ' + str(z0fitErr) + ' ' 
                            + str(sigma_xErr) + ' ' + str(sigma_yErr) + ' ' +  str(thetaErr)
                            )
                    
                    
                    
                    data_FIONA[MN_fiona].append(datatempFIONA)

                    
                    #print 'M#, datatemp FIONA: ' , MN_fiona , '   ' ,  datatempFIONA
                    
                    xposfit_FIONA_0[MN_fiona].append(np.around(x0fitImageJ,2))
                    yposfit_FIONA_0[MN_fiona].append(np.around(y0fitImageJ,2))
                    
                    if (eccentricityTemp <= maxEcc):
                        xposfit_FIONA_1[MN_fiona].append(np.around(x0fitImageJ,2))
                        yposfit_FIONA_1[MN_fiona].append(np.around(y0fitImageJ,2))
                        print '\nFIONA_1 OK'
                        

                        if (AfitErrP <= maxPer) & (x0fitImageJErrPx <= maxPx) & (y0fitImageJErrPx <= maxPx):
                            xposfit_FIONA_2[MN_fiona].append(np.around(x0fitImageJ,2))
                            yposfit_FIONA_2[MN_fiona].append(np.around(y0fitImageJ,2))
                            print '\nFIONA_2 OK'
                            
                            
                        else:
                            print 'FINOA_2 failed'
                            print 'AfitErrP ', AfitErrP
                            print 'x0fitImageJErrPx ', x0fitImageJErrPx
                            print 'y0fitImageJErrPx ', y0fitImageJErrPx
                            
                    else:
                        print 'FIONA_1 failed'
                        print 'eccentricityTemp ', eccentricityTemp
                        

                

                 
                    ####################################################
                
                except:
                    print 'FIONA fitting failed'
                    eccentricityTemp = 1.2
                    datatempFIONA = (str(f1n) + ' ' + str(999) + ' ' + str(999) + ' ' 
                            + str(999) +  ' ' + str(999) + ' ' + str(999) + ' ' + str(999) + ' ' +  str(999) + ' ' + str(eccentricityTemp) + ' '
                            + str(BG_IntensityAve[MN_fiona][f1n]) + ' ' + str(centerIntensityAve[MN_fiona][f1n]) + ' ' 
                            + str(centerIntensityMax[MN_fiona][f1n]) + ' ' + str(SBratioAve[MN_fiona][f1n]) + ' ' + str(SBratioMax[MN_fiona][f1n]) + ' Err: '       
                            + str(999) +  ' ' + str(999) + ' ' + str(999) + ' ' + str(999) + ' ' 
                            + str(999) + ' ' + str(999) + ' ' +  str(999)
                            )
                     
                    data_FIONA[MN_fiona].append(datatempFIONA)
                    eccentricityFor_FIONA[MN_fiona].append(eccentricityTemp)

                    
                    #print 'M#, datatemp FIONA: ' , MN_fiona , '   ' ,  datatempFIONA
                    
        
        print 'FIONA end\n'
        ######################  FIONA only part                              ############
        #################################################################################
      
        
        
        
        
        
        
        
        #################################################################################
        #########   All frame gsrimp calculation    #####################################
        print '\nAF gsrimp frameStartInter f1n: ', f1n
        
        if f1n < frameStart:
            print 'This frame skipped f1n: ', f1n
            continue
        

        for f2n in range(f1n +1 , frameEnd + 1):            
            
            for fN in range(len(xpos)): #for each feature
                if SBratioAve[fN][f1n] <= SvsB_Tolerance :
                    print '\npass, SBratioAve[fN][f1n]', fN, f1n, SBratioAve[fN][f1n]
                    pass
                elif SBratioAve[fN][f2n] <= SvsB_Tolerance :
                    print '\npass, SBratioAve[fN][f2n]', fN, f2n, SBratioAve[fN][f2n]
                    pass
                else:
                        
                    print '\ng1: f1n = ', f1n, '    f2n = ', f2n, '   feature #', fN
                    a1 = frames[f1n][yc[fN]-NCenterSubPixel:yc[fN]+NCenterSubPixel+1, xc[fN]-NCenterSubPixel:xc[fN]+NCenterSubPixel+1]
                    a2 = frames[f2n][yc[fN]-NCenterSubPixel:yc[fN]+NCenterSubPixel+1, xc[fN]-NCenterSubPixel:xc[fN]+NCenterSubPixel+1]
                    
                    #print 'centerIntensityMax[f1n] - centerIntensityMax[f2n] : ', centerIntensityMax[fN][f1n] - centerIntensityMax[fN][f2n]
                    
                    if int(centerIntensityMax[fN][f1n]) - int(centerIntensityMax[fN][f2n]) >= 0:
                        cframe = np.array(a1, dtype=np.int32) - np.array(a2, dtype=np.int32)
                        #print 'frame #', f2n, ': subtraction is +'
                    else:
                        cframe = np.array(a2, dtype=np.int32) - np.array(a1, dtype=np.int32)
                        #print 'frame #', f2n, ': subtraction is -'                       
                
    
                    frameTemp0_raveled = cframe.ravel()                
                
                    # fit the data  for all frames
                    initial_guess = (2000, NCenterSubPixel+1,  NCenterSubPixel+1, 3, 3, 10, 200)
     
                    k = 1
                   
                    Afit = 0.0
                    x0fit = 0.0
                    y0fit = 0.0
                    z0fit = 0.0
                    sigma_x = 0.0
                    sigma_y = 0.0
                    theta = 0.0
                    
                    try:
                        popt, pcov = opt.curve_fit(twoD_Gaussian, (x3, y3), frameTemp0_raveled, p0 = initial_guess)
                        #print '\nAF shrimp popt: ', popt
                        #print '\nAF shrimp pcov: \n', pcov
                        
                        Afit = popt[0]
                        x0fit = popt[1]
                        y0fit = popt[2]
                        
                        x0fitImageJ = xpos[fN] + (NCenterSubPixel - y0fit)  # to correct for imageJ x position
                        y0fitImageJ = ypos[fN] - (NCenterSubPixel - x0fit)
                        
                        z0fit = popt[6]
    
                        sigma_x = popt[3]
                        sigma_y = popt[4]
                        theta = popt[5]


                        # AF gshrimp AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaErr
                        AfitErr = np.sqrt(pcov[0,0])                    
                        x0fitImageJErr = np.sqrt(pcov[1,1])
                        y0fitImageJErr = np.sqrt(pcov[2,2]) 
                        z0fitErr = np.sqrt(pcov[6,6]) 
                        sigma_xErr = np.sqrt(pcov[3,3]) 
                        sigma_yErr = np.sqrt(pcov[4,4]) 
                        thetaErr = np.sqrt(pcov[5,5])
                                        
                                        
                        # AF gshrimp Percentage of the errors: AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaErr
                        AfitErrP = np.around(100.0 * AfitErr / Afit, 0)                    
                        x0fitImageJErrPx = np.around(np.sqrt(pcov[1,1]),2)
                        y0fitImageJErrPx = np.around(np.sqrt(pcov[2,2]),2) 
                        #z0fitErrP = np.around(100.0 * z0fitErr / z0fit, 0)      
                        #sigma_xErrP = np.around(100.0 * sigma_xErr / sigma_x, 0)   
                        #sigma_yErrP = np.around(100.0 * sigma_yErr / sigma_y, 0)  
                        #thetaErrDeg = np.around(np.sqrt(pcov[5,5]),2)                       
                        

                        '''
                        print 'AfitErrP ', AfitErrP
                        print 'x0fitImageJErrPx ', x0fitImageJErrPx
                        print 'y0fitImageJErrPx ', y0fitImageJErrPx
                        print 'z0fitErrP ', z0fitErrP
                        print 'sigma_xErrP ', sigma_xErrP
                        print 'sigma_yErrP ', sigma_yErrP
                        print 'thetaErrDeg ', thetaErrDeg
                        '''                        
                                        
                        ########### eccentricicity calculation  ############
                       
                        a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
                        b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
                        c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2) 
                        
                        #print 'a, b, c =', a,b,c
            
                        A = a
                        B = 2*b
                        C = c
                                
                        eccentricityTemp = np.sqrt(2*(np.sqrt((A-C)**2 + B**2)) / (A+C + np.sqrt(((A-C)**2 + B**2)))  )
                        #print 'eccentricityTemp = ', eccentricityTemp                        
                    
                        datatemp = (str(f1n) + ' ' + str(f2n) + ' ' + str(centerIntensityMax[fN][f1n]) + ' ' + str(centerIntensityMax[fN][f2n]) + ' ' 
                                + str(Afit) + ' ' + str(x0fitImageJ) + ' ' + str(y0fitImageJ) + ' ' + str(z0fit) + ' ' 
                                + str(sigma_x) + ' ' + str(sigma_y) + ' ' +  str(theta) + ' ' + str(eccentricityTemp) + ' Err: '  
                                + str(AfitErr) +  ' ' + str(x0fitImageJErr) + ' ' + str(y0fitImageJErr) + ' ' + str(z0fitErr) + ' ' 
                                + str(sigma_xErr) + ' ' + str(sigma_yErr) + ' ' +  str(thetaErr)
                                )
                         
                        data_gshrimp[fN].append(datatemp)
                        eccentricityFor_gshrimp[fN].append(eccentricityTemp)
                        
                        #print 'datatemp: ', 'M# ' , fN , ' ' ,  datatemp
                        
                        #eccentricityAll[fN].append(eccentricityTemp)

                        xposfit_gshrimp0[fN].append(np.around(x0fitImageJ,2))
                        yposfit_gshrimp0[fN].append(np.around(y0fitImageJ,2))
                        
                        if (eccentricityTemp <= maxEcc) & (Afit > 100):
                            
                            xposfit_gshrimp1[fN].append(np.around(x0fitImageJ,2))
                            yposfit_gshrimp1[fN].append(np.around(y0fitImageJ,2))
                            print '\ng1 OK'                           
                            
                            
                            if (AfitErrP <= maxPer) & (x0fitImageJErrPx <= maxPx) & (y0fitImageJErrPx <= maxPx):
                                xposfit_gshrimp2[fN].append(np.around(x0fitImageJ,2))
                                yposfit_gshrimp2[fN].append(np.around(y0fitImageJ,2))
                                print '\ng2 OK'
                            else:
                                print '\ng2 failed'
                                print 'AfitErrP ', AfitErrP
                                print 'x0fitImageJErrPx ', x0fitImageJErrPx
                                print 'y0fitImageJErrPx ', y0fitImageJErrPx
                                
                        else:
                            print '\ng1 failed'
                            print 'eccentricityTemp ', eccentricityTemp
                            print 'Afit ', Afit
                                
                     
                        ####################################################
                    
                    except:
                        print '\nAF shrimp fitting failed\n'
                        eccentricityTemp = 1.2
                        datatemp = (str(f1n) + ' ' + str(f2n) + ' ' + str(centerIntensityMax[fN][f1n]) + ' ' + str(centerIntensityMax[fN][f2n]) + ' ' 
                                + str(999) + ' ' + str(999) + ' ' + str(999) + ' ' + str(999) + ' ' 
                                + str(999) + ' ' + str(999) + ' ' +  str(999) + ' ' + str(eccentricityTemp)   + ' Err: '
                                + str(999) +  ' ' + str(999) + ' ' + str(999) + ' ' + str(999) + ' ' 
                                + str(999) + ' ' + str(999) + ' ' +  str(999)
                                )
                         
                        data_gshrimp[fN].append(datatemp)
                        eccentricityFor_gshrimp[fN].append(eccentricityTemp)
                        
                        #print 'datatemp: ', 'M# ' , fN , ' ' ,  datatemp
                        
    

    
    return (xposfit_gshrimp1, yposfit_gshrimp1, centerImage, data_gshrimp, eccentricityFor_gshrimp, xposfit_gshrimp0, yposfit_gshrimp0
        , xposfit_FIONA_0, yposfit_FIONA_0, data_FIONA, xposfit_gshrimp2, yposfit_gshrimp2, xposfit_FIONA_1, yposfit_FIONA_1, xposfit_FIONA_2, yposfit_FIONA_2
        , eccentricityFor_FIONA)








'''
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
'''


  
class MainTracking(wx.Frame):
    #ImageFilePath = ''
 
    #----------------------------------------------------------------------
    def __init__(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, "afSHRImP_GUI", size=(1000, 1000), pos=(0,0))
                                  
        self.MainPanel = wx.Panel(self, wx.ID_ANY)
        
        

        #self.btn_SwitchToAFShrimp = wx.Button(self.MainPanel, pos=(300,10), label="Switch To afShrimp")
        #self.btn_SwitchToAFShrimp.Bind(wx.EVT_BUTTON, self.SwitchToAFShrimp)
        self.btn_SwitchToDataFileAnalysis = wx.Button(self.MainPanel, pos=(500,10), label="Switch To Data File Analysis")
        self.btn_SwitchToDataFileAnalysis.Bind(wx.EVT_BUTTON, self.SwitchToDataFileAnalysis)
        self.btn_SwitchToIntensityTrajectoryPlots = wx.Button(self.MainPanel, pos=(700,10), label="Switch To Intensity Trajectory Plots with Images")
        self.btn_SwitchToIntensityTrajectoryPlots.Bind(wx.EVT_BUTTON, self.SwitchToIntensityTrajectoryPlots)


        
        
        #self.xposImageJ = ''
        #self.yposImageJ = ''
        
        '''
        wx.StaticText(self.MainPanel, 1, "Feature position (x,y) in ImageJ", pos=(50, 15))
        wx.StaticText(self.MainPanel, -1, "x:", pos=(250, 15))
        self.xpos = wx.TextCtrl(self.MainPanel, -1, "100", pos=(270, 13))
        wx.StaticText(self.MainPanel, -1, "y:", pos=(400, 15))
        self.ypos = wx.TextCtrl(self.MainPanel, -1, "100", pos=(420, 13))
        '''
        
 
        #btn1 = wx.Button(MainPanel, pos=(50,400), label="Open Data File")
        #btn1.Bind(wx.EVT_BUTTON, self.onOpenFile1)
        
        #btn1.SetSize(80,80)
        
        wx.StaticText(self.MainPanel, -1, 'Grayscale Tiff Images', pos=(50, 20))
        
        self.btn_OpenImage = wx.Button(self.MainPanel, pos=(50,50), label="Open Image File")
        self.ImageFilePath = self.btn_OpenImage.Bind(wx.EVT_BUTTON, self.onOpenImageFile)
        self.btn_OpenImage_Text = wx.StaticText(self.MainPanel, -1, str(self.ImageFilePath), pos=(170, 53))
        
        
        
        
        wx.StaticText(self.MainPanel, -1, "Choose a Frame #:", pos=(170, 103))
        #self.Nframe = wx.TextCtrl(self.MainPanel, -1, "0", pos=(280, 100), size=(40,-1))
        self.Nframe = wx.SpinCtrl(self.MainPanel, -1, pos=(280, 100), size=(50,-1))
        self.Nframe.SetValue(0)
              
        
        self.btn_FindingMolecules = wx.Button(self.MainPanel, pos=(50,100), label="Finding Molecules")
        self.btn_FindingMolecules.Hide()
        self.btn_FindingMolecules.Bind(wx.EVT_BUTTON, self.FindingFeatureFunc)
        
        wx.StaticText(self.MainPanel, -1, "Average # frames:", pos=(370, 103))
        self.NofFramesForAveraging = wx.SpinCtrl(self.MainPanel, -1, pos=(480, 100), size=(50,-1))
        self.NofFramesForAveraging.SetValue(5)
        
        wx.StaticText(self.MainPanel, -1, "Feature Size:", pos=(580, 103))
        self.FeatureSize = wx.SpinCtrl(self.MainPanel, -1,  pos=(660, 100), size=(50,-1))
        self.FeatureSize.SetValue(7)
        
        wx.StaticText(self.MainPanel, -1, "Min Intensity:", pos=(400, 133))
        self.MinIntensity = wx.TextCtrl(self.MainPanel, -1, str(ThrethholdIntensity), pos=(480, 130), size=(70,-1))
        
        wx.StaticText(self.MainPanel, -1, "Number of Features:", pos=(170, 133))
        self.NofFeatures = wx.StaticText(self.MainPanel, -1, "None", pos=(300, 133), style=5)
        


        self.btn_ShowRawImages = wx.Button(self.MainPanel, pos=(750,100), label="Show Raw Images")
        self.btn_ShowRawImages.Bind(wx.EVT_BUTTON, self.ShowRawImages)
        self.btn_ShowRawImages.Show()



        self.btn_SelectData = wx.Button(self.MainPanel, pos=(880,100), label="Select Images")
        self.btn_SelectData.Bind(wx.EVT_BUTTON, self.SelectData)
        self.btn_SelectData.Show()
        
        












        
        self.btn_IntensityTrajectoryPlots = wx.Button(self.MainPanel, pos=(50,180), label="Intensity Trajectory Plots")
        self.btn_IntensityTrajectoryPlots.Hide()
        self.btn_IntensityTrajectoryPlots.Bind(wx.EVT_BUTTON, self.IntensityTrajectoryPlots)
        

        self.btn_DriftPlots = wx.Button(self.MainPanel, pos=(250,180), label="x,y drift Plots")
        self.btn_DriftPlots.Show()
        self.btn_DriftPlots.Bind(wx.EVT_BUTTON, self.DriftPlots)
        
        self.btn_DriftPlotsSave = wx.Button(self.MainPanel, pos=(450,180), label="Save drift Fig, Data")        
        self.btn_DriftPlotsSave.Bind(wx.EVT_BUTTON, self.DriftPlotsSave)
        self.btn_DriftPlotsSave.Show()
        


        
                
        
        self.btn_2Dfits = wx.Button(self.MainPanel, pos=(550,200), label="2D Gaussian Fit")
        self.btn_2Dfits.Hide()
        self.btn_2Dfits.Bind(wx.EVT_BUTTON, self.twoDGaussFit)
        wx.StaticText(self.MainPanel, -1, "Fitting 2D Gaussian for all features", pos=(550, 154))
        
        
        
        self.btn_gshrimp = wx.Button(self.MainPanel, pos=(50,220), label="AF gSHRImP")       
        self.btn_gshrimp.Hide()
        self.btn_gshrimp.Bind(wx.EVT_BUTTON, self.onRungshrimp)
        
        self.cb_show_afSHRImP_resultFigures_yn = wx.CheckBox(self.MainPanel, -1, 'Show figures', (50, 250))
        self.cb_show_afSHRImP_resultFigures_yn.SetValue(True)
        
        
        wx.StaticText(self.MainPanel, -1, "Start frame #", pos=(200, 222))
        self.frameStart = wx.TextCtrl(self.MainPanel, -1, "0", pos=(280, 220))
        wx.StaticText(self.MainPanel, -1, "End frame #", pos=(450, 222))
        self.frameEnd = wx.TextCtrl(self.MainPanel, -1, "4", pos=(520, 220))
        wx.StaticText(self.MainPanel, -1, "Adjacent Frame #s", pos=(670, 222))
        self.additionalFramesN = wx.SpinCtrl(self.MainPanel, -1, str(additionalFramesN), pos=(780, 220), size=(50,-1))
        self.additionalFramesN.SetValue(additionalFramesN)
        
        
        wx.StaticText(self.MainPanel, -1, "Min S(ave)/B(ave) ", pos=(180, 252))
        self.SvsB_Tolerance = wx.TextCtrl(self.MainPanel, -1, "1.1", pos=(280, 250))
        wx.StaticText(self.MainPanel, -1, "Max eccentricity", pos=(420, 252))
        self.MaxEccentricity = wx.TextCtrl(self.MainPanel, -1, "0.65", pos=(520, 250))
        wx.StaticText(self.MainPanel, -1, "Max fit err %: A", pos=(180, 282))
        self.MaxPercentage_Tolerance = wx.TextCtrl(self.MainPanel, -1, "1000", pos=(280, 280))        
        wx.StaticText(self.MainPanel, -1, "Max fit err px: x, y", pos=(410, 282))
        self.MaxPixel_Tolerance = wx.TextCtrl(self.MainPanel, -1, "0.1", pos=(520, 280))        
        
                
        
        

        
        wx.StaticText(self.MainPanel, -1, "Data and Figure Files Prefix", pos=(50, 322))
        self.FilePrefix = wx.TextCtrl(self.MainPanel, -1, FilePrefixInput, pos=(200, 320))
        
         
        self.btn_Save_figures = wx.Button(self.MainPanel, pos=(50,350), label="Save Figures")  
        self.btn_Save_figures.Hide()
        self.btn_Save_figures.Bind(wx.EVT_BUTTON, self.SaveFigures)
        self.text_SaveFig = wx.StaticText(self.MainPanel, -1, ' ', pos=(50, 380)) 
        
        self.cb_Save_figures_AllInfoOnly = wx.CheckBox(self.MainPanel, -1, 'All Info Only', (150, 355))
        self.cb_Save_figures_AllInfoOnly.SetValue(True)
        
        print self.cb_Save_figures_AllInfoOnly.GetValue()
        
        
        self.btn_Save_gdata = wx.Button(self.MainPanel, pos=(250,350), label="Save Data")   
        self.btn_Save_gdata.Hide()
        self.btn_Save_gdata.Bind(wx.EVT_BUTTON, self.Save_gdata)
        
        self.text_SaveGdata = wx.StaticText(self.MainPanel, -1, " ", pos=(250, 380)) 
        
 
        self.btn_SaveIntensityTrajectoryPlots = wx.Button(self.MainPanel, pos=(450,350), label="Save Intensity Trajectory Plots Only")
        self.btn_SaveIntensityTrajectoryPlots.Hide()
        self.btn_SaveIntensityTrajectoryPlots.Bind(wx.EVT_BUTTON, self.SaveIntensityTrajectoryPlotsOnly)
        
        
        


        self.MainPanel.Bind(wx.EVT_PAINT, self.on_paint)
        
        
        
        
       
        
        self.btn_PlotSingle = wx.Button(self.MainPanel, pos=(50,450), label="Plot Raw Image Data")        
        self.btn_PlotSingle.Bind(wx.EVT_BUTTON, self.PlotSingleFeatureImages)
        self.btn_PlotSingle.Hide()
        
        
        wx.StaticText(self.MainPanel, 1, "Min < Color range < Max", pos=(220, 405))
        self.singleImageColorMin = wx.TextCtrl(self.MainPanel, -1, "1000", pos=(230, 420), size=(40,-1))
        self.singleImageColorMax = wx.TextCtrl(self.MainPanel, -1, "6000", pos=(290, 420), size=(40,-1))
        


        self.btn_PlotSingleImage_colorRangeUpdate = wx.Button(self.MainPanel, pos=(350,410), label="Update color range")        
        self.btn_PlotSingleImage_colorRangeUpdate.Bind(wx.EVT_BUTTON, self.PlotSingleFeatureImagesUpdate)
        self.btn_PlotSingleImage_colorRangeUpdate.Show()
        
        
        
        
        
        wx.StaticText(self.MainPanel, 1, "ImageJ (x,y)", pos=(220, 453))
        #wx.StaticText(self.MainPanel, -1, "x:", pos=(200, 552))
        self.xposSingle1 = wx.TextCtrl(self.MainPanel, -1, str(SingleFeature_x), pos=(300, 450), size=(40,-1))
        #wx.StaticText(self.MainPanel, -1, "y:", pos=(300, 552))
        self.yposSingle1 = wx.TextCtrl(self.MainPanel, -1, str(SingleFeature_y), pos=(350, 450), size=(40,-1))

        wx.StaticText(self.MainPanel, 1, "- from frame #         0", pos=(420, 453), style=wx.TE_READONLY)
        self.xposRawImageStartFrameN = wx.TextCtrl(self.MainPanel, -1, "0", pos=(510, 450), size=(40,-1))
        self.xposRawImageStartFrameN.Hide()
        wx.StaticText(self.MainPanel, 1, "to frame #", pos=(570, 453))
        self.xposRawImageEndFrameN = wx.TextCtrl(self.MainPanel, -1, "79", pos=(640, 450), size=(40,-1))
        self.xposSingle1FrameTotalText2 = wx.StaticText(self.MainPanel, -1, '', pos=(600, 430)) 



        self.radio_NAvergeFrames_rawImage = wx.RadioBox(self.MainPanel, -1, choices=['Single (1 frame)', 'Adjacent (3 frames)'], label='# of Average Frames', pos=(700, 415), size=wx.Size(140, 60), style=wx.RA_SPECIFY_ROWS)
        #self.radio_ColorMap_Custom = wx.TextCtrl(self.MainPanel, -1, " ", pos=(820, 535), size=(60,-1))
        self.radio_NAvergeFrames_rawImage.SetSelection(0)
        

        
        
        wx.StaticText(self.MainPanel, 1, "Noise Amplitude:", pos=(810, 493))
        self.xposSingle1NoiseAmp = wx.TextCtrl(self.MainPanel, -1, "0", pos=(890, 490), size=(40,-1))


        self.btn_PlotSingleSave = wx.Button(self.MainPanel, pos=(850,450), label="Save Image Figures")        
        self.btn_PlotSingleSave.Bind(wx.EVT_BUTTON, self.PlotSingleFeatureImagesSave)
        self.btn_PlotSingleSave.Hide()
        





        self.btn_SingleTwoFrameShrimp = wx.Button(self.MainPanel, pos=(50,490), label="Two-Image Shrimp")        
        self.btn_SingleTwoFrameShrimp.Bind(wx.EVT_BUTTON, self.SingleTwoFrameShrimp)
        self.btn_SingleTwoFrameShrimp.Show()
        wx.StaticText(self.MainPanel, 1, "Image 1 frames:", pos=(220, 493))
        self.xposSingle1FrameStart1 = wx.TextCtrl(self.MainPanel, -1, "4", pos=(320, 490), size=(40,-1))       
        self.xposSingle1FrameStart2 = wx.TextCtrl(self.MainPanel, -1, "9", pos=(370, 490), size=(40,-1))       
        wx.StaticText(self.MainPanel, 1, "Image 2 frames:", pos=(450, 493))
        self.xposSingle1FrameEnd1 = wx.TextCtrl(self.MainPanel, -1, "19", pos=(550, 490), size=(40,-1))       
        self.xposSingle1FrameEnd2 = wx.TextCtrl(self.MainPanel, -1, "35", pos=(600, 490), size=(40,-1))               
            
        
        #wx.StaticText(self.MainPanel, -1, "Sc", pos=(850, 460))
        self.radio_ColorMap = wx.RadioBox(self.MainPanel, -1, choices=['jet', 'gray', 'rainbow' ], label='Color Map', pos=(800, 520), size=wx.Size(120, 80), style=wx.RA_SPECIFY_ROWS)
        #self.radio_ColorMap_Custom = wx.TextCtrl(self.MainPanel, -1, " ", pos=(820, 535), size=(60,-1))
        self.radio_ColorMap.SetSelection(0)
        
        self.btn_SingleTwoFrameShrimpAdjacent = wx.Button(self.MainPanel, pos=(50,530), label="Two-Image Shrimp2")        
        self.btn_SingleTwoFrameShrimpAdjacent.Bind(wx.EVT_BUTTON, self.SingleTwoFrameShrimpAdjacent)
        self.btn_SingleTwoFrameShrimpAdjacent.Show()
        wx.StaticText(self.MainPanel, 1, "Image 1 frame:", pos=(220, 533))
        self.xposSingle1FrameStart1Adjacent = wx.TextCtrl(self.MainPanel, -1, "4", pos=(320, 530), size=(40,-1))       
        wx.StaticText(self.MainPanel, 1, "Image 2 frame:", pos=(400, 533))
        self.xposSingle1FrameEnd1Adjacent = wx.TextCtrl(self.MainPanel, -1, "9", pos=(500, 530), size=(40,-1))   
        wx.StaticText(self.MainPanel, 1, "+- # of adjacent frames:", pos=(580, 533))
        self.xposSingle1AdjacentFrameN = wx.TextCtrl(self.MainPanel, -1, "1", pos=(730, 530), size=(40,-1))               
            
        
        
                
        
        
        
        
        



        self.btn_IntensitySingle = wx.Button(self.MainPanel, pos=(50,630), label="Intensity Trajectory Singles for afSHRImP")       
        self.btn_IntensitySingle.Hide()
        self.btn_IntensitySingle.Bind(wx.EVT_BUTTON, self.IntensityTrajectoryPlotsSingle)   
        
        
        
        self.btn_DriftPlots_with_xy_features = wx.Button(self.MainPanel, pos=(450,630), label="x,y drift Plots")
        self.btn_DriftPlots_with_xy_features.Show()
        self.btn_DriftPlots_with_xy_features.Bind(wx.EVT_BUTTON, self.DriftPlots_with_xy_features)
        
        self.btn_DriftPlotsSave_with_xy_features = wx.Button(self.MainPanel, pos=(550,630), label="Save drift Fig, Data")        
        self.btn_DriftPlotsSave_with_xy_features.Bind(wx.EVT_BUTTON, self.DriftPlotsSave_with_xy_features)
        self.btn_DriftPlotsSave_with_xy_features.Show()
        




        wx.StaticText(self.MainPanel, 1, "ImageJ (x,y)", pos=(50, 653))
        self.xposSingleM = [[] for _ in xrange(10)]      
        self.yposSingleM = [[] for _ in xrange(10)]    
        
        
        self.xposSingleM[0] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(150, 670), size=(40,-1))
        self.yposSingleM[0] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(200, 670), size=(40,-1))
        
        self.xposSingleM[1] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(300, 670), size=(40,-1))
        self.yposSingleM[1] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(350, 670), size=(40,-1))        
        
        self.xposSingleM[2] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(450, 670), size=(40,-1))
        self.yposSingleM[2] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(500, 670), size=(40,-1))        
        
        self.xposSingleM[3] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(600, 670), size=(40,-1))
        self.yposSingleM[3] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(650, 670), size=(40,-1))        
        
        self.xposSingleM[4] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(750, 670), size=(40,-1))
        self.yposSingleM[4] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(800, 670), size=(40,-1))
        
    
                
        self.xposSingleM[5] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(150, 720), size=(40,-1))
        self.yposSingleM[5] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(200, 720), size=(40,-1))
        
        self.xposSingleM[6] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(300, 720), size=(40,-1))
        self.yposSingleM[6] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(350, 720), size=(40,-1))        
        
        self.xposSingleM[7] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(450, 720), size=(40,-1))
        self.yposSingleM[7] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(500, 720), size=(40,-1))        
        
        self.xposSingleM[8] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(600, 720), size=(40,-1))
        self.yposSingleM[8] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(650, 720), size=(40,-1))        
        
        self.xposSingleM[9] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(750, 720), size=(40,-1))
        self.yposSingleM[9] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(800, 720), size=(40,-1))
        
    
                          
                
                        
        

    #----------------------------------------------------------------------#
    #----------------------------------------------------------------------#


    def SwitchToAFShrimp(self, event):
        frameSM2Tracking = MainTracking()
        frameSM2Tracking.Show()
        self.Close()
         

    def SwitchToDataFileAnalysis(self, event):
        frameDataFileAnalysis = MainTrackingDataFileAnalysis()
        frameDataFileAnalysis.Show()
        self.Close()
        

    def SwitchToIntensityTrajectoryPlots(self, event):
        IntensityTrajectoryPlots = MainTrackingIntensityOnlyPlots()
        IntensityTrajectoryPlots.Show()
        self.Close()
        




        
        
    def on_paint(self, event):
        dc = wx.PaintDC(event.GetEventObject())
        dc.Clear()
        dc.SetPen(wx.Pen("BLACK", 4))
        dc.DrawLine(0, 400, 1000, 400)  
        dc.SetPen(wx.Pen("BLACK", 4))
        dc.DrawLine(0, 620, 1000, 620)  
        
            
        
     #----------------------------------------------------------------------
 
 
 
    def onOpenFile1(self, event):
        plt.close("all")

        """
        Create and show the Open FileDialog
        """
        dlg = wx.FileDialog(
            self, message="Choose a file",
            defaultFile="",
            wildcard=wildcard,defaultDir=os.getcwd(),
            style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
            )
        if dlg.ShowModal() == wx.ID_OK:
            paths = dlg.GetPaths()
            print "You chose the following file(s):"
            for path in paths:
                print path
        print dlg
        print 'path = ', dlg.GetPath()
        print 'filename = ', dlg.GetFilename()
        #dlg.Destroy()
        
        self.Show()
        
        
    def onOpenImageFile(self, event):
        
        self.btn_FindingMolecules.Hide()
        self.btn_IntensityTrajectoryPlots.Hide()
        self.btn_SaveIntensityTrajectoryPlots.Hide()
        self.btn_gshrimp.Hide()
        self.btn_Save_figures.Hide()
        self.btn_Save_gdata.Hide()
        self.btn_IntensitySingle.Hide()
        
        
        self.text_SaveFig.SetLabel(' ')
        self.text_SaveGdata.SetLabel(" ")
        
        
        

        tp_ver = tp.__version__[:3]
        if float(tp_ver[:3]) != 0.3:
            MessageBox(self, 'The trackpy version is not 0.3.#\nUpdate trackpy', 'Error')



        
        plt.close('Data Image Frames')
        
        '''
        try:
            plt.close(self.fig_ShowAfewFrames)
        except:
            pass

        '''



          
        """
        Create and show the Open FileDialog
        """
        plt.close()
        dlg2 = wx.FileDialog(
            self, message="Choose a file",
            defaultFile="",
            wildcard=wildcard,
            style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
            )
        if dlg2.ShowModal() == wx.ID_OK:
            paths = dlg2.GetPaths()
            print "You chose the following file(s):"
            for path in paths:
                print path
        print ''
        print 'dlg2 = ', dlg2
        print 'path = ', dlg2.GetPath()
        print 'path type', type(dlg2.GetPath())
        print 'filename = ', dlg2.GetFilename()
        #dlg.Destroy()
        
        self.Show()
        
        filenameTemp = dlg2.GetFilename()
        self.btn_OpenImage_Text.SetLabel(str(filenameTemp))

        
        self.filepath = dlg2.GetPath()
        self.filenameOnly = dlg2.GetFilename()
        self.filedirectoryOnly = dlg2.GetDirectory()
        
        self.save_fig_ShowAfewFrames = self.ShowAFewFrames()
        
        self.btn_FindingMolecules.Show()
        self.btn_IntensitySingle.Show()
        self.btn_PlotSingle.Show()
        
        
    
    def ShowAFewFrames(self):
        frames = pims.TiffStack(self.filepath)
        self.LoadedFrames = frames
        #Nframes = frames._count
        
        Nframes = len(frames) #trackpy 0.3
            
            
            
        self.LastFrameN = Nframes - 1
        self.xposSingle1FrameTotalText2.SetLabel('Max End # '+str(Nframes-1))
        
        if self.LastFrameN < int(self.xposRawImageEndFrameN.GetValue()):
            self.xposRawImageEndFrameN.SetLabel(str(self.LastFrameN))
            
            
        
        print 'Nframes = ', Nframes
        self.fig_ShowAfewFrames = plt.figure('Data Image Frames')
        #fig_ShowAfewFrames.suptitle('test title', fontsize=20)
        self.fig_ShowAfewFrames.subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.55, hspace = 0.2)
        
        
        plt.figtext(0.5, 0.95, str(self.filedirectoryOnly)+'\n' + str(os.path.basename(self.filepath))
            ,ha='center', color='black', weight='normal', size='small')
            
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(1010,30,800, 800)
        
        if Nframes > 15:
            k = 16
        else:
            k = Nframes
            
        for n in range(k):
            print 'frame # = ', n
            plt.subplot(4,4,n+1)
            plt.title('frame # ' + str(n))
            #plt.hold(True)
            plt.imshow(frames[n], origin='upper')  # for trackpy 3.0
            #plt.imshow(np.rot90(frames[n] ,3), origin='upper')
            #plt.imshow(frames[n], vmin=1000, vmax=16000)
            
        plt.ion()
        plt.show()

        self.figFileInfo_text = str(self.filedirectoryOnly)+'\n' + str(os.path.basename(self.filepath))
        
        return self.fig_ShowAfewFrames
        
        
        
        
        
        

    def FindingFeatureFunc(self, event):
        plt.close('Finding Features')
        featureSize = int(self.FeatureSize.GetValue())
        minIntensity = int(self.MinIntensity.GetValue())
        
        #frames = tp.TiffStack(self.filepath)
        frames = self.LoadedFrames
        
        self.N_xPixel = frames.frame_shape[0]
        self.N_yPixel = frames.frame_shape[1]

        NframeTemp = int(self.Nframe.GetValue())
        NFramesForAveraging = int(self.NofFramesForAveraging.GetValue())
        
        AnalysisConditions = ( 'Chosen F # ' + str(self.Nframe.GetValue()) + '   Ave # F: ' + str(self.NofFramesForAveraging.GetValue()) + 
        '   f size: '+ str(self.FeatureSize.GetValue()) +  '   Min I: ' + str(self.MinIntensity.GetValue()) )
        
        print '\nFingding features from frame # ', NframeTemp
        
        
        #framesSummed = np.rot90(frames[0], -1) * 0
        framesSummed = frames[0] * 0 # for trackpy 3.0
        
        
        for n in range(NframeTemp, NframeTemp + NFramesForAveraging):
            print 'Summing frames for feature finding: frame # ', n
            #framesSummed += frames[n]       
            #framesSummed += np.rot90(frames[n], -1)
            framesSummed += frames[n] # for trackpy 3.0
        framesSummed /= NFramesForAveraging
        
        f = tp.locate(framesSummed, featureSize, minmass=minIntensity, invert=False)
        
        #NFeatures = f.values.size/8
        NFeatures = len(f.index) #trackpy 0.3


        
        
        if NFeatures == 0:
            MessageBox(self, '0 features. Adjust the Min Intensity', 'Error')
            return
            
            
            
            
            
        
        print 'NFeatures = ', NFeatures
        self.NofFeatures.SetLabel(str(NFeatures))
        
        #plt.close()
        self.fig_FindingFeatures = plt.figure('Finding Features')
        plt.figtext(0.5, 0.92, str(self.filedirectoryOnly)+'\n' + str(os.path.basename(self.filepath)) + '\n' + AnalysisConditions + 
            '\n From frame # ' + str(NframeTemp) + '  with average of ' + str(NFramesForAveraging) + ' frames' + ' ,    Total # of features: '+ str(NFeatures)
            ,ha='center', color='black', weight='normal', size='small')

        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(1010,30,700, 700)
        #tp.annotate(f, frames[NframeTemp])
        tp.annotate(f, framesSummed)
        
        self.xpos = []
        self.ypos = []
        for n in range(NFeatures):
            x_temp = f.values[n][0]
            y_temp = f.values[n][1]
            self.xpos.append(x_temp)
            self.ypos.append(y_temp)

        #N_xPixel = frames.frame_shape[0]
        #N_yPixel = frames.frame_shape[1]
         
        self.xposImageJ = np.array(self.xpos) + 1
        #self.yposImageJ = np.array(self.ypos)
        self.yposImageJ = np.array(self.ypos) + 1    #trackpy 0.3
        
        
        print 'self.xpos ', self.xpos
        print 'self.ypos ', self.ypos
        #print 'xposImageJ', self.xposImageJ
        #print 'yposImageJ', self.yposImageJ
        
        plt.ion()
        plt.show()
        
        self.btn_IntensityTrajectoryPlots.Show()
        frameStartValueTemp =  self.Nframe.GetValue()
        self.frameStart.SetValue(str(frameStartValueTemp))
        frameEndValueTemp = str(int(frameStartValueTemp) + 4)
        self.frameEnd.SetValue(frameEndValueTemp)
        
        #self.frameStartSingle.SetValue(frameStartValueTemp)
        #frameEndValueTemp = str(int(frameStartValueTemp) + 99)
        #self.frameEndSingle.SetValue(frameEndValueTemp)        
        
     
     
    def IntensityTrajectoryPlots(self, event):
        print '\nIntensityTrajectoryPlots'
        
        self.IntensityTrajectoryAve, self.IntensityTrajectoryMax, self.IntensityTrajectoryMax5Ave = IntensityTrajectoryPlotsFunc(self.filepath, self.xposImageJ, self.yposImageJ)
        
        
        AnalysisConditions = ( 'Chosen F # ' + str(self.Nframe.GetValue()) + '   Ave # F: ' + str(self.NofFramesForAveraging.GetValue()) + 
        '   f size: '+ str(self.FeatureSize.GetValue()) +  '   Min I: ' + str(self.MinIntensity.GetValue()) )
        
        


        
        NFeatures = len(self.IntensityTrajectoryMax)
        Nfig = int(  math.ceil((NFeatures/30.0)))
        
        TotalFeatureN = 'Total # ' + str(NFeatures)
        
        #self.fig_IntensityTrajectory = [[]] * Nfig
        self.fig_IntensityTrajectory = [[] for _ in xrange(Nfig)]
        
        k = 0
        fn = 0
        for n in range(Nfig):
            print 'Intensity Trajectory Nfig n = ', n
            self.fig_IntensityTrajectory[n] = plt.figure('IntensityTrajectory_'+ str(n), figsize = (18, 9))
            #fig_ShowAfewFrames.suptitle('test title', fontsize=20)
            self.fig_IntensityTrajectory[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.25, hspace = 0.35)
            plt.figtext(0.5, 0.95, str(self.filedirectoryOnly)+'\n' + str(os.path.basename(self.filepath) + '\n' + AnalysisConditions   )
                ,ha='center', color='black', weight='normal', size='small')
            self.fig_IntensityTrajectory[n].text(0.05, 0.96, TotalFeatureN + '    Max 5 ave', ha="left", va="bottom", size="medium",color="red")
            
            
            for m in range(30):
                #print 'subplot # = ', m
                plt.subplot(6,5, m+1)
                plt.tick_params(labelsize=9)
                plt.title('# ' + str(m + fn) + ',  (' + str(int(self.xposImageJ[m+fn])) + ', ' + str(int(self.yposImageJ[m+fn]))+')'  , fontsize=9)
                #plt.plot(self.IntensityTrajectoryAve[m + fn])
                #plt.plot(self.IntensityTrajectoryMax[m + fn])
                plt.plot(self.IntensityTrajectoryMax5Ave[m + fn])
                
                k += 1
                if k%30 ==0:
                    fn += 30
                if k == NFeatures:
                    break
                
                

        
        
        #print 'self.fig_IntensityTrajectory[0]', self.fig_IntensityTrajectory[0]
        #print 'self.fig_IntensityTrajectory[1]', self.fig_IntensityTrajectory[1]
        plt.show()
        
        
        self.btn_gshrimp.Show()
        self.btn_SaveIntensityTrajectoryPlots.Show()
        self.btn_DriftPlots.Show()
        
        
        
        
        
    def twoDGaussFit(self, event):
        plt.close('Finding Features')
        
        print 'running twoDGaussianFit'
        
        xpos = self.xposImageJ
        ypos = self.yposImageJ
        
        ManyTwoDFits(self.filepath, xpos, ypos, int(self.Nframe.GetValue()))
        
        print 'Nframe = ', int(self.Nframe.GetValue())
        
        self.Close()
        





    def ShowRawImages(self, event):
        print '\nstarting ShowRawImages'
        
        try:
            for n in range(len(self.fig_AllDataCenterImagesOnly)):
                plt.close(self.fig_AllDataCenterImagesOnly[n])
        except: pass
    
        
        xpos = self.xposImageJ
        ypos = self.yposImageJ
         
        ################################################
        x_ImageJ, y_ImageJ = xpos, ypos
        
        #xc = y_ImageJ
        #yc = self.N_xPixel - x_ImageJ +1 
        
        xc = x_ImageJ - 1 # for trackpy 3.0
        yc = y_ImageJ # for trackpy 3.0
            
        ################################################
        
        
        frames = self.LoadedFrames
        
        TotalNframes = len(frames)
        print 'TotalNframes: ', TotalNframes
        

        NCenterSubPixel = 7  # number of pixels from the center

        #frameTemp0 = frame0[yc-NCenterSubPixel:yc+NCenterSubPixel + 1, xc-NCenterSubPixel:xc+NCenterSubPixel+1]
        x3 = np.linspace(0, 2 * NCenterSubPixel, 2* NCenterSubPixel + 1)
        y3 = np.linspace(0, 2 * NCenterSubPixel, 2*NCenterSubPixel+1)
        x3, y3 = np.meshgrid(x3, y3)
        


        centerImage = [[] for _ in xrange(len(xc))] # centerImage[Molecule#][Frame#]
    
        frameNTemp = int(self.Nframe.GetValue())
        print 'center image data for frame # ', frameNTemp
        for k in range(len(xpos)):     
            print 'frame #: ', k
            
            frameLargeTemp = frames[frameNTemp][yc[k]-NCenterSubPixel:yc[k]+NCenterSubPixel+1, xc[k]-NCenterSubPixel:xc[k]+NCenterSubPixel+1]
            #frameCenterTemp = frames[frameTemp][yc[k]-2:yc[k]+3, xc[k]-2:xc[k]+3]
            
            centerImage[k].append(frameLargeTemp) # it's not matching to the actual frame #, only first a few are saved








        NFeatures = len(xpos)
        
        Nrow = 5
        Ncol = 10
        NimagePerFig = Nrow*Ncol
  
       
        NCfig = int(  math.ceil((NFeatures/(float(NimagePerFig)))))
        self.fig_AllDataCenterImagesOnly = [[] for _ in xrange(NCfig)]

        x3 = np.linspace(0, 15, 16)
        y3 = np.linspace(0, 15, 16)
        x3, y3 = np.meshgrid(x3, y3)
        
        fsn = int(self.Nframe.GetValue()) # frame start number
        k = 0
        fn = 0
        for n in range(NCfig):
            #print 'n = ', n
            self.fig_AllDataCenterImagesOnly[n] = plt.figure('Raw_Data_Images'+ str(n), figsize = (18, 9))
            plt.clf()
            self.fig_AllDataCenterImagesOnly[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.02, right = 0.99, wspace = 0.36, hspace = 0.4)
            plt.figtext(0.5, 0.94, str(self.filedirectoryOnly) + '\n'+ 'F#: '+ str(frameNTemp) + ',      ' + str(self.filenameOnly)
                ,ha='center', color='black', weight='normal', size='small')

            for m in range(NimagePerFig):
                try:
                    plt.subplot(Nrow, Ncol, m+1)
                    plt.title('M# ' + str(m+fn) + ',(' + str(int(x_ImageJ[m+fn])) + ',' + str(int(y_ImageJ[m+fn]))+')'  , fontsize=8, color = 'black')
                    plt.tick_params(labelsize=7)
                    plt.imshow(centerImage[m+fn][0].reshape(15, 15), interpolation = 'None', cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                except:
                    print '\nNo center image for ', 'M# = ', str(m + fn), '  F# ', str(fsn)

                k += 1
                if k%NimagePerFig ==0:
                    fn += NimagePerFig
                if k == NFeatures:
                    break





     
    def SelectData(self, event):
        #self.SelectData_window = wx.Frame(None, title='SelectData', size=(800, 1000), pos=(800,0))
        #self.SelectData_window.Show()
        

        self.SelectData_frame = wx.Frame(self, title='SelectData', size=(850, 1000), pos=(800,0))
        self.SelectData_window = wx.ScrolledWindow(self.SelectData_frame, -1)
        self.SelectData_window.SetScrollbars(1, 1, 800, int(len(self.xposImageJ)/10)*25 + 200)
        self.SelectData_frame.Show()
        
        
        
        

        wx.StaticText(self.SelectData_window, -1, "Data File Prefix", pos=(160, 35))
        self.FilePrefixForSelected = wx.TextCtrl(self.SelectData_window, -1, str(FilePrefixInput), pos=(250, 30), size=(150, -1))
        
        
        btn_ShowFiguresSelectedData = wx.Button(self.SelectData_window, pos=(50,30), label="Show Figures")
        btn_ShowFiguresSelectedData.Bind(wx.EVT_BUTTON, self.SelectData_ShowFigures)
        
        
        self.btn_DriftPlots_Selected = wx.Button(self.SelectData_window, pos=(420,30), label="Drift Plots Selected")
        self.btn_DriftPlots_Selected.Bind(wx.EVT_BUTTON, self.DriftPlots_Selected)
        self.btn_DriftPlots_Selected.Show()
        
        
        
        
        wx.StaticText(self.SelectData_window, -1, 'Select Data', pos=(50, 110)) 
                
        #self.xpos_ImageJ = [[] for _ in xrange(150)]
        self.cb_Data_yn = [[] for _ in xrange(len(self.xposImageJ))]
        rn = 0
        cn = 0
        countertemp = 0
        for n in range(len(self.xposImageJ)):
            self.cb_Data_yn[n] = wx.CheckBox(self.SelectData_window, -1, '# '+str(n), (50 + 70*cn, 145 + 25*rn))
            self.cb_Data_yn[n].SetValue(True)
            countertemp += 1
            cn += 1
            if countertemp % 10 == 0:
                rn += 1
                cn = 0
            
        self.btn_SelectAll = wx.Button(self.SelectData_window, pos=(150,105), label="Select All")
        self.btn_SelectAll.Bind(wx.EVT_BUTTON, self.SelectDataAll)
        self.btn_SelectAll.Show()
        
        self.btn_SelectNone = wx.Button(self.SelectData_window, pos=(250,105), label="Deselect All")
        self.btn_SelectNone.Bind(wx.EVT_BUTTON, self.SelectDataNone)
        self.btn_SelectNone.Show()
        
        
        
        
        
    def SelectDataAll(self, event):
        print '\nSelect All '
        for n in range(len(self.cb_Data_yn)):
            self.cb_Data_yn[n].SetValue(True)
    def SelectDataNone(self, event):
        print '\nDeselect All '
        for n in range(len(self.cb_Data_yn)):
            self.cb_Data_yn[n].SetValue(False)
            
                        
    def SelectData_ShowFigures(self, event):
        print '\n SelectData_ShowFigures '

        try:
            for n in range(len(self.fig_AllDataCenterImagesOnly_Selected)):
                plt.close(self.fig_AllDataCenterImagesOnly_Selected[n])
        except: pass
    
        xpos = []
        ypos = []
        for n in range(len(self.xposImageJ)):
            if self.cb_Data_yn[n].GetValue():
                xpos.append(self.xposImageJ[n])
                ypos.append(self.yposImageJ[n])
        
         
         
         
        ################################################
        x_ImageJ, y_ImageJ = np.array(xpos), np.array(ypos)
        
        #xc = y_ImageJ
        #yc = self.N_xPixel - x_ImageJ + 1
        
        xc = x_ImageJ - 1 # for trackpy 3.0
        yc = y_ImageJ # for trackpy 3.0
            
        ################################################
        
        
        frames = self.LoadedFrames
        
        TotalNframes = len(frames)
        print 'TotalNframes: ', TotalNframes
        

        NCenterSubPixel = 7  # number of pixels from the center

        #frameTemp0 = frame0[yc-NCenterSubPixel:yc+NCenterSubPixel + 1, xc-NCenterSubPixel:xc+NCenterSubPixel+1]
        x3 = np.linspace(0, 2 * NCenterSubPixel, 2* NCenterSubPixel + 1)
        y3 = np.linspace(0, 2 * NCenterSubPixel, 2*NCenterSubPixel+1)
        x3, y3 = np.meshgrid(x3, y3)
        


        centerImage = [[] for _ in xrange(len(xc))] # centerImage[Molecule#][Frame#]
    
        frameNTemp = int(self.Nframe.GetValue())
        print 'center image data for frame # ', frameNTemp
        for k in range(len(xpos)):            
            
            frameLargeTemp = frames[frameNTemp][yc[k]-NCenterSubPixel:yc[k]+NCenterSubPixel+1, xc[k]-NCenterSubPixel:xc[k]+NCenterSubPixel+1]
            #frameCenterTemp = frames[frameTemp][yc[k]-2:yc[k]+3, xc[k]-2:xc[k]+3]
            
            centerImage[k].append(frameLargeTemp) # it's not matching to the actual frame #, only first a few are saved


        NFeatures = len(xpos)
        
        Nrow = 5
        Ncol = 10
        NimagePerFig = Nrow*Ncol
  
       
        NCfig = int(  math.ceil((NFeatures/(float(NimagePerFig)))))
        self.fig_AllDataCenterImagesOnly_Selected = [[] for _ in xrange(NCfig)]

        x3 = np.linspace(0, 15, 16)
        y3 = np.linspace(0, 15, 16)
        x3, y3 = np.meshgrid(x3, y3)
        
        fsn = int(self.Nframe.GetValue()) # frame start number
        k = 0
        fn = 0
        for n in range(NCfig):
            #print 'n = ', n
            self.fig_AllDataCenterImagesOnly_Selected[n] = plt.figure('Raw_Data_Images_Selected_'+ str(n), figsize = (18, 9))
            self.fig_AllDataCenterImagesOnly_Selected[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.02, right = 0.99, wspace = 0.36, hspace = 0.4)
            plt.figtext(0.5, 0.94, str(self.filedirectoryOnly) + '\n'+ 'F#: '+ str(frameNTemp) + ',      ' + str(self.filenameOnly)
                ,ha='center', color='black', weight='normal', size='small')

            for m in range(NimagePerFig):
                try:
                    plt.subplot(Nrow, Ncol, m+1)
                    plt.title('M# ' + str(m+fn) + ',(' + str(int(x_ImageJ[m+fn])) + ',' + str(int(y_ImageJ[m+fn]))+')'  , fontsize=8, color = 'black')
                    plt.tick_params(labelsize=7)
                    plt.imshow(centerImage[m+fn][0].reshape(15, 15), interpolation = 'None', cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                except:
                    print '\nNo center image for ', 'M# = ', str(m + fn), '  F# ', str(fsn)

                k += 1
                if k%NimagePerFig ==0:
                    fn += NimagePerFig
                if k == NFeatures:
                    break


        self.xpos_Selected = xpos
        self.ypos_Selected = ypos




    def DriftPlots_Selected(self, event):
        print '\n\nDriftPlots\n' 
        
        frames = self.LoadedFrames
 
        N_xPixel = frames.frame_shape[0]
        N_yPixel = frames.frame_shape[1]

        xpos = []
        ypos = []
        
        for n in range(len(self.xpos_Selected)):
            if self.xpos_Selected[n] >= 10 and self.xpos_Selected[n] <= N_xPixel - 10 and self.ypos_Selected[n] >= 10 and self.ypos_Selected[n] <= N_yPixel - 10:
                xpos.append(self.xpos_Selected[n])
                ypos.append(self.ypos_Selected[n])





        
        print '# of molecules: xpos len = ', len(xpos)
        

        print '\nMaxEccentricity', float(self.MaxEccentricity.GetValue())        
        FrameN_fxy2, fx2, fy2, Fdata, FeccentricityData, fx2_dx, fy2_dy, FrameN_fxy2_driftAve, fx2_dx_ave, fy2_dy_ave = FIONA_drift(
            self.filepath, xpos, ypos, int(self.frameStart.GetValue()), int(self.frameEnd.GetValue())
            , float(self.SvsB_Tolerance.GetValue()), float(self.MaxEccentricity.GetValue())
            , float(self.MaxPercentage_Tolerance.GetValue()), float(self.MaxPixel_Tolerance.GetValue()), int(self.additionalFramesN.GetValue()))
        

        print 'FrameN_fxy2 ', FrameN_fxy2
        print 'fx2', fx2
        print 'fx2_dx', fx2_dx
        
        
        try:
            plt.close('drift')
        except: pass
        
        plt.figure('drift', figsize=(8,12))
        plt.subplot(211)
        plt.title('drift_x')
        for n in range(len(fx2)):
            plt.plot(FrameN_fxy2[n], fx2_dx[n])
        plt.plot(FrameN_fxy2_driftAve, fx2_dx_ave, 'ro')
        
        plt.subplot(212)
        plt.title('drift_y')
        for n in range(len(fy2)):
            plt.plot(FrameN_fxy2[n], fy2_dy[n])
        plt.plot(FrameN_fxy2_driftAve, fy2_dy_ave, 'ro')
        
        
            
            

    

    def DriftPlots(self, event):
        print '\n\nDriftPlots\n' 
        
        frames = self.LoadedFrames
 
        N_xPixel = frames.frame_shape[0]
        N_yPixel = frames.frame_shape[1]

        xpos = []
        ypos = []
        
        for n in range(len(self.xposImageJ)):
            if self.xposImageJ[n] >= 10 and self.xposImageJ[n] <= N_xPixel - 10 and self.yposImageJ[n] >= 10 and self.yposImageJ[n] <= N_yPixel - 10:
                xpos.append(self.xposImageJ[n])
                ypos.append(self.yposImageJ[n])





        
        print '# of molecules: xpos len = ', len(xpos)
        

        print '\nMaxEccentricity', float(self.MaxEccentricity.GetValue())        
        FrameN_fxy2, fx2, fy2, Fdata, FeccentricityData, fx2_dx, fy2_dy, FrameN_fxy2_driftAve, fx2_dx_ave, fy2_dy_ave = FIONA_drift(
            self.filepath, xpos, ypos, int(self.frameStart.GetValue()), int(self.frameEnd.GetValue())
            , float(self.SvsB_Tolerance.GetValue()), float(self.MaxEccentricity.GetValue())
            , float(self.MaxPercentage_Tolerance.GetValue()), float(self.MaxPixel_Tolerance.GetValue()), int(self.additionalFramesN.GetValue()))
        

        print 'FrameN_fxy2 ', FrameN_fxy2
        print 'fx2', fx2
        print 'fx2_dx', fx2_dx
        
        
        try:
            plt.close('drift')
        except: pass
        
        
        
        self.fig_drift = plt.figure('drift', figsize=(8,12))
        plt.figtext(0.5, 0.94, str(self.filenameOnly) +'\n# of total molecules: ' + str(len(FrameN_fxy2))
                ,ha='center', color='black', weight='normal', size='small')
        plt.subplot(211)
        plt.title('drift_x')
        for n in range(len(fx2)):
            plt.plot(FrameN_fxy2[n], fx2_dx[n])
        plt.plot(FrameN_fxy2_driftAve, fx2_dx_ave, 'ro')
        plt.ylim(-0.2, 0.2)
        
        plt.subplot(212)
        plt.title('drift_y')
        for n in range(len(fy2)):
            plt.plot(FrameN_fxy2[n], fy2_dy[n])
        plt.plot(FrameN_fxy2_driftAve, fy2_dy_ave, 'ro')
        plt.ylim(-0.2, 0.2)
        
        self.FrameN_fxy2_driftAve = FrameN_fxy2_driftAve
        self.fx2_dx_ave = fx2_dx_ave
        self.fy2_dy_ave = fy2_dy_ave
        
   




    def DriftPlotsSave(self, event):
        print '\n\nSaving drift data'


        figFileName = self.filenameOnly[:30]
        todayDate = time.strftime("%Y%m%d_%Hh%Mm")

        self.fig_drift.savefig(figFileName + '_' + todayDate + '_Drift' + '.png')
        plt.close(self.fig_drift)

        ff = open(figFileName + '_' + todayDate + '_Drift_FN_dx_dy' + '.txt','w')
        for n in range(len(self.FrameN_fxy2_driftAve)):
            ff.write(str(self.FrameN_fxy2_driftAve[n]) + ' ' + str(self.fx2_dx_ave[n]) + ' '  + str(self.fy2_dy_ave[n]) + '\n'   )    
                
        ff.close()       

        
        










 

    def DriftPlots_with_xy_features(self, event):
        print '\n\nDriftPlots with some x,y data\n' 
        
        frames = self.LoadedFrames
 
        N_xPixel = frames.frame_shape[0]
        N_yPixel = frames.frame_shape[1]

        


        xposImageJtemp = []
        yposImageJtemp = []
        
        for n in range(len(self.xposSingleM)):
            xtemp = int(self.xposSingleM[n].GetValue())
            ytemp = int(self.yposSingleM[n].GetValue())
            if (xtemp != 0) & (ytemp != 0):
                xposImageJtemp.append(xtemp)
                yposImageJtemp.append(ytemp)
                
        xpos = np.array(xposImageJtemp)
        ypos = np.array(yposImageJtemp)
        
        
        print 'xpos: ', xpos
        print 'ypos: ', ypos
        
        





        
        print '# of molecules: xpos len = ', len(xpos)
        

        print '\nMaxEccentricity', float(self.MaxEccentricity.GetValue())        
        FrameN_fxy2, fx2, fy2, Fdata, FeccentricityData, fx2_dx, fy2_dy, FrameN_fxy2_driftAve, fx2_dx_ave, fy2_dy_ave = FIONA_drift(
            self.filepath, xpos, ypos, int(self.frameStart.GetValue()), int(self.frameEnd.GetValue())
            , float(self.SvsB_Tolerance.GetValue()), float(self.MaxEccentricity.GetValue())
            , float(self.MaxPercentage_Tolerance.GetValue()), float(self.MaxPixel_Tolerance.GetValue()), int(self.additionalFramesN.GetValue()))
        

        print 'FrameN_fxy2 ', FrameN_fxy2
        print 'fx2', fx2
        print 'fx2_dx', fx2_dx
        
        
        try:
            plt.close('drift')
        except: pass
        
        
        
        self.fig_drift = plt.figure('drift', figsize=(8,12))
        plt.figtext(0.5, 0.94, str(self.filenameOnly) +'\n# of total molecules: ' + str(len(FrameN_fxy2))
                ,ha='center', color='black', weight='normal', size='small')
        plt.subplot(211)
        plt.title('drift_x')
        for n in range(len(fx2)):
            plt.plot(FrameN_fxy2[n], fx2_dx[n])
        plt.plot(FrameN_fxy2_driftAve, fx2_dx_ave, 'ro')
        plt.ylim(-0.2, 0.2)
        
        plt.subplot(212)
        plt.title('drift_y')
        for n in range(len(fy2)):
            plt.plot(FrameN_fxy2[n], fy2_dy[n])
        plt.plot(FrameN_fxy2_driftAve, fy2_dy_ave, 'ro')
        plt.ylim(-0.2, 0.2)
        
        self.FrameN_fxy2_driftAve = FrameN_fxy2_driftAve
        self.fx2_dx_ave = fx2_dx_ave
        self.fy2_dy_ave = fy2_dy_ave
   






    def DriftPlotsSave_with_xy_features(self, event):
        print '\n\nSaving drift data'


        figFileName = self.filenameOnly[:30]
        todayDate = time.strftime("%Y%m%d_%Hh%Mm")

        self.fig_drift.savefig(figFileName + '_' + todayDate + '_Drift' + '.png')
        plt.close(self.fig_drift)

        ff = open(figFileName + '_' + todayDate + '_Drift_FN_dx_dy' + '.txt','w')
        for n in range(len(self.FrameN_fxy2_driftAve)):
            ff.write(str(self.FrameN_fxy2_driftAve[n]) + ' ' + str(self.fx2_dx_ave[n]) + ' '  + str(self.fy2_dy_ave[n]) + '\n'   )    
                
        ff.close()       

        
        













        
        
        
        
    
    def onRungshrimp(self, event): 
        print '\nonRungshrimp'
        
        self.text_SaveFig.SetLabel('processing ... ')
        self.text_SaveGdata.SetLabel(' ')
        
    
        
         


        try:
            for n in range(len(self.fig_IntensityTrajectory)):
                plt.close(self.fig_IntensityTrajectory[n])
        except: pass
            
            
        try:        
            for n in range(len(self.fig_FionaDataAll)):    
                plt.close(self.fig_FionaDataAll[n])
        except: pass
            

        try:        
            for n in range(len(self.fig_FionaData1)):    
                plt.close(self.fig_FionaData1[n])
        except: pass
            
        try:        
            for n in range(len(self.fig_FionaData2)):    
                plt.close(self.fig_FionaData2[n])
        except: pass
            
            
        try:        
            for n in range(len(self.fig_gshrimp1)):    
                plt.close(self.fig_gshrimp1[n])
                plt.close(self.fig_gshrimp2[n])
        except: pass
        
        
        try:        
            for n in range(len(self.fig_geccentricityData)):    
                plt.close(self.fig_geccentricityData[n])
        except: pass

                    
    
    
    

        startTime = time.strftime("%Y-%m%d, %Hh %Mm")
        startTimeHour = int(time.strftime("%H"))
        startTimeMinute = int(time.strftime("%M"))
        startTimeDay = int(time.strftime("%d"))
        
        
        AnalysisConditions = ( 'Chosen F # ' + str(self.Nframe.GetValue()) + ',   Ave # F: ' + str(self.NofFramesForAveraging.GetValue()) + 
        ',   f size: '+ str(self.FeatureSize.GetValue()) +  ',   Min I: ' + str(self.MinIntensity.GetValue())  +
        ',   start F # '+ str(self.frameStart.GetValue()) + ',   end F # ' + str(self.frameEnd.GetValue()) +
        ',   S/B: ' + str(self.SvsB_Tolerance.GetValue()) + ',   Max Ecc: ' + str(self.MaxEccentricity.GetValue()) +  
        ',   Max A err%: ' + str(self.MaxPercentage_Tolerance.GetValue()) +  ',   Max Px x,y err: ' + str(self.MaxPixel_Tolerance.GetValue()) +
        ', #AveFrames: ' + str(int(self.additionalFramesN.GetValue())*2+1)  )
        
        g1Conditions = ('[AF g1]   S/B : ' + str(self.SvsB_Tolerance.GetValue()) + ',   Max Ecc: ' + str(self.MaxEccentricity.GetValue()) )
        g2Conditions = ('[AF g2]   S/B : ' + str(self.SvsB_Tolerance.GetValue()) + ',   Max Ecc: ' + str(self.MaxEccentricity.GetValue()) +  
        '\n, Max A err% : ' + str(self.MaxPercentage_Tolerance.GetValue()) +  ',   Max Px x,y err: ' + str(self.MaxPixel_Tolerance.GetValue()) )

        f0Conditions = ('[FIONA 0]   S/B : ' + str(self.SvsB_Tolerance.GetValue())  )
        f1Conditions = ('[FIONA 1]   S/B : ' + str(self.SvsB_Tolerance.GetValue()) + ',   Max Ecc: ' + str(self.MaxEccentricity.GetValue())   )  
        f2Conditions = ('[FIONA 2]   S/B : ' + str(self.SvsB_Tolerance.GetValue()) + ',   Max Ecc: ' + str(self.MaxEccentricity.GetValue()) +  
        '\n, Max A err% : ' + str(self.MaxPercentage_Tolerance.GetValue()) +  ',   Max Px x,y err: ' + str(self.MaxPixel_Tolerance.GetValue())  )
        
        
        
        xpos = self.xposImageJ
        ypos = self.yposImageJ
        print '# of molecules: xpos len = ', len(xpos)
        

        print '\nMaxEccentricity', float(self.MaxEccentricity.GetValue())        
        gx1, gy1, CenterImages, self.gdata, geccentricityData, gx0, gy0, fx0, fy0, self.data_FIONA, gx2, gy2, fx1,fy1, fx2, fy2, FeccentricityData = gshrimpAllFramesByEachFrame(
            self.filepath, xpos, ypos, int(self.frameStart.GetValue()), int(self.frameEnd.GetValue())
            , float(self.SvsB_Tolerance.GetValue()), float(self.MaxEccentricity.GetValue())
            , float(self.MaxPercentage_Tolerance.GetValue()), float(self.MaxPixel_Tolerance.GetValue()), int(self.additionalFramesN.GetValue()))
        
 
        
        NFeatures = len(xpos)
        Nfig = int(  math.ceil((NFeatures/60.0)))
        

    
    
             
        self.fig_FionaDataAll = [[] for _ in xrange(Nfig)]
        k = 0
        fn = 0
        for n in range(Nfig):
            #print 'n = ', n
            self.fig_FionaDataAll[n] = plt.figure('ALL_FIONA_'+ str(n), figsize = (14.4, 9))
            self.fig_FionaDataAll[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.4, hspace = 0.4)
            plt.figtext(0.5, 0.94, str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly) +'\nFIONA-0  _ ' + AnalysisConditions
                ,ha='center', color='black', weight='normal', size='small')
            self.fig_FionaDataAll[n].text(0.03, 0.96, f0Conditions, ha="left", va="bottom", size=9, color="red")
            
                


            for m in range(60):
                #print 'FIONA subplot # = ', m
                plt.subplot(6,10, m+1)
                plt.title('M# ' + str(m + fn) + ',  (' + str(int(self.xposImageJ[m+fn])) + ', ' + str(int(self.yposImageJ[m+fn]))+')'  , fontsize=9)
                plt.tick_params(labelsize=7)
                
                plt.scatter(fx0[m + fn], fy0[m + fn], marker = 'x', c = range(len(fx0[m + fn])),  s=35, vmin=0, vmax= len(fx0[m + fn]))
                plt.xlim( (np.around(xpos[k]) - 2, np.around(xpos[k]) + 3 ) )
                plt.ylim( (np.around(ypos[k]) - 2, np.around(ypos[k]) + 3 ) )        
                
                k += 1
                if k%60 ==0:
                    fn += 60
                if k == NFeatures:
                    break
        
        
        
             
        self.fig_FionaData1 = [[] for _ in xrange(Nfig)]
        k = 0
        fn = 0
        for n in range(Nfig):
            #print 'n = ', n
            self.fig_FionaData1[n] = plt.figure('FIONA1_'+ str(n), figsize = (14.4, 9))
            self.fig_FionaData1[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.4, hspace = 0.4)
            plt.figtext(0.5, 0.94, str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly) +'\nFIONA-1  _ ' + AnalysisConditions
                ,ha='center', color='black', weight='normal', size='small')
            self.fig_FionaData1[n].text(0.03, 0.96, f1Conditions, ha="left", va="bottom", size=9, color="red")
            
                


            for m in range(60):
                #print 'FIONA subplot # = ', m
                plt.subplot(6,10, m+1)
                plt.title('M# ' + str(m + fn) + ',  (' + str(int(self.xposImageJ[m+fn])) + ', ' + str(int(self.yposImageJ[m+fn]))+')'  , fontsize=9)
                plt.tick_params(labelsize=7)
                
                plt.scatter(fx1[m + fn], fy1[m + fn], marker = 'x', c = range(len(fx1[m + fn])),  s=35, vmin=0, vmax= len(fx1[m + fn]))
                plt.xlim( (np.around(xpos[k]) - 2, np.around(xpos[k]) + 3 ) )
                plt.ylim( (np.around(ypos[k]) - 2, np.around(ypos[k]) + 3 ) )        
                
                k += 1
                if k%60 ==0:
                    fn += 60
                if k == NFeatures:
                    break
        
             
             
        
        
             
        self.fig_FionaData2 = [[] for _ in xrange(Nfig)]
        k = 0
        fn = 0
        for n in range(Nfig):
            #print 'n = ', n
            self.fig_FionaData2[n] = plt.figure('FIONA2_'+ str(n), figsize = (14.4, 9))
            self.fig_FionaData2[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.4, hspace = 0.4)
            plt.figtext(0.5, 0.94, str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly) +'\nFIONA-2  _ ' + AnalysisConditions
                ,ha='center', color='black', weight='normal', size='small')
            self.fig_FionaData2[n].text(0.03, 0.96, f2Conditions, ha="left", va="bottom", size=9, color="red")
            
                


            for m in range(60):
                #print 'FIONA subplot # = ', m
                plt.subplot(6,10, m+1)
                plt.title('M# ' + str(m + fn) + ',  (' + str(int(self.xposImageJ[m+fn])) + ', ' + str(int(self.yposImageJ[m+fn]))+')'  , fontsize=9)
                plt.tick_params(labelsize=7)
                
                plt.scatter(fx2[m + fn], fy2[m + fn], marker = 'x', c = range(len(fx2[m + fn])),  s=35, vmin=0, vmax= len(fx2[m + fn]))
                plt.xlim( (np.around(xpos[k]) - 2, np.around(xpos[k]) + 3 ) )
                plt.ylim( (np.around(ypos[k]) - 2, np.around(ypos[k]) + 3 ) )        
                
                k += 1
                if k%60 ==0:
                    fn += 60
                if k == NFeatures:
                    break
        
             
             
            



        self.fig_gshrimp1 = [[] for _ in xrange(Nfig)]
        k = 0
        fn = 0
        for n in range(Nfig):
            #print 'n = ', n
            self.fig_gshrimp1[n] = plt.figure('AF_shrimp1_'+ str(n), figsize = (14.4, 9))
            self.fig_gshrimp1[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.4, hspace = 0.4)
            plt.figtext(0.5, 0.94, str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly) +'\nAF shrimp1  _  ' + AnalysisConditions
                ,ha='center', color='black', weight='normal', size='small')
            self.fig_gshrimp1[n].text(0.03, 0.96, g1Conditions, ha="left", va="bottom", size=9, color="red")
                


            for m in range(60):
                #print 'AF gshrimp subplot # = ', m
                plt.subplot(6,10, m+1)
                plt.title('M# ' + str(m + fn) + ',  (' + str(int(self.xposImageJ[m+fn])) + ', ' + str(int(self.yposImageJ[m+fn]))+')'  , fontsize=9)
                plt.tick_params(labelsize=7)
     
                plt.scatter(gx1[m + fn], gy1[m + fn], marker = 'x', c = range(len(gx1[m + fn])),  s=35, vmin=0, vmax= len(gx1[m + fn]))
                plt.xlim( (np.around(xpos[k]) - 2, np.around(xpos[k]) + 3 ) )
                plt.ylim( (np.around(ypos[k]) - 2, np.around(ypos[k]) + 3 ) )        


                k += 1
                if k%60 ==0:
                    fn += 60
                if k == NFeatures:
                    break
         


        self.fig_gshrimp2 = [[] for _ in xrange(Nfig)]
        k = 0
        fn = 0
        for n in range(Nfig):
            #print 'n = ', n
            self.fig_gshrimp2[n] = plt.figure('AF_shrimp2_'+ str(n), figsize = (14.4, 9))
            self.fig_gshrimp2[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.4, hspace = 0.4)
            plt.figtext(0.5, 0.94, str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly) +'\nAF shrimp2  _  ' + AnalysisConditions
                ,ha='center', color='black', weight='normal', size='small')
            self.fig_gshrimp2[n].text(0.03, 0.96, g2Conditions, ha="left", va="bottom", size=9, color="red")
                


            for m in range(60):
                #print 'AF gshrimp2 subplot # = ', m
                plt.subplot(6,10, m+1)
                plt.title('M# ' + str(m + fn) + ',  (' + str(int(self.xposImageJ[m+fn])) + ', ' + str(int(self.yposImageJ[m+fn]))+')'  , fontsize=9)
                plt.tick_params(labelsize=7)
     
                plt.scatter(gx2[m + fn], gy2[m + fn], marker = 'x', c = range(len(gx2[m + fn])),  s=35, vmin=0, vmax= len(gx2[m + fn]))
                plt.xlim( (np.around(xpos[k]) - 2, np.around(xpos[k]) + 3 ) )
                plt.ylim( (np.around(ypos[k]) - 2, np.around(ypos[k]) + 3 ) )        


                k += 1
                if k%60 ==0:
                    fn += 60
                if k == NFeatures:
                    break
         



  
             
             
        self.fig_geccentricityData = [[] for _ in xrange(Nfig)]
        k = 0
        fn = 0
        for n in range(Nfig):
            #print 'n = ', n
            self.fig_geccentricityData[n] = plt.figure('eccentricity_shrimp_'+ str(n), figsize = (14.4, 9))
            self.fig_geccentricityData[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.4, hspace = 0.4)
            plt.figtext(0.5, 0.94, str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly) +'\nAll AF shrimp eccentricity  _  ' + AnalysisConditions
                ,ha='center', color='black', weight='normal', size='small')
                


            for m in range(60):
                #print 'e subplot # = ', m
                plt.subplot(6,10, m+1)
                plt.title('M# ' + str(m + fn) + ',  (' + str(int(self.xposImageJ[m+fn])) + ', ' + str(int(self.yposImageJ[m+fn]))+')'  , fontsize=9)
                plt.tick_params(labelsize=7)
                
                try:
                    plt.hist(geccentricityData[m + fn])
                    plt.xlim(0, 1.2)
                    
                except:
                    print 'hist error'
                    print ' geccentricityData[m + fn]: ', geccentricityData[m + fn]
                
                
                k += 1
                if k%60 ==0:
                    fn += 60
                if k == NFeatures:
                    break
        
             
             

             
                          
             

        
        NCfig = int(  math.ceil((NFeatures/6.0)))
        self.fig_AllDataCenterImages = [[] for _ in xrange(NCfig)]

        x3 = np.linspace(0, 15, 16)
        y3 = np.linspace(0, 15, 16)
        x3, y3 = np.meshgrid(x3, y3)
        
        fsn = int(self.frameStart.GetValue()) # frame start number
        k = 0
        fn = 0
        for n in range(NCfig):
            #print 'n = ', n
            self.fig_AllDataCenterImages[n] = plt.figure('AF_shrimp_&_center_images_All'+ str(n), figsize = (18, 9))
            self.fig_AllDataCenterImages[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.02, right = 0.99, wspace = 0.36, hspace = 0.4)
            plt.figtext(0.5, 0.94, str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly) +'\n' + AnalysisConditions
                ,ha='center', color='black', weight='normal', size='small')
                


            for m in range(6):
                #print 'subplot # = ', m      
                

                plt.subplot(6,14, 14*m+1)
                plt.title('M# ' + str(m + fn) + ',  (' + str(int(self.xposImageJ[m+fn])) + ', ' + str(int(self.yposImageJ[m+fn]))+')'  , fontsize=9, color = 'red')
                plt.tick_params(labelsize=6)
                plt.plot(self.IntensityTrajectoryMax5Ave[m + fn])


               
                #print len(CenterImages[0])
                for cfn in range(5 ): #  images for center frames, only a few frames from the start frame are saved, so not matching to actual frame #.
                    #print 'center frame # cfn = ', cfn + fsn
                    try:
                            
                        plt.subplot(6,14, 14*m+cfn+2)
                        plt.title('F# ' + str(fsn + cfn) + ',  (' + str(int(self.xposImageJ[m+fn])) + ', ' + str(int(self.yposImageJ[m+fn]))+')'  , fontsize=9, color = 'black')
                        plt.tick_params(labelsize=7)
                        plt.imshow(np.rot90(CenterImages[m+fn][cfn].reshape(15, 15) ,3), interpolation = 'None', cmap=plt.cm.jet, origin='upper',
                                   extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                    except:
                        print '\nNo center image for ', 'M# = ', str(m + fn), '  F# ', str(cfn + fsn)


                

                plt.subplot(6,14, 14*m+7)
                plt.title('F '+str(len(fx0[m + fn]))+':'+str(len(fx1[m + fn]))+':'+str(len(fx2[m + fn])) , fontsize=8)
                plt.tick_params(labelsize=6)
                
                try:
                    plt.hist(FeccentricityData[m + fn])
                    plt.xlim(0, 1.2)
                    
                except:
                    print 'FIONA hist error'
                    print 'FeccentricityData[m + fn]: ', FeccentricityData[m + fn]
                    
                    
                
                plt.subplot(6,14, 14*m+8)
                plt.title('F0 S/B' , fontsize=8, color='DeepPink')
                plt.tick_params(labelsize=7)
                plt.scatter(fx0[m + fn], fy0[m + fn], marker = 'x', c = range(len(fx0[m + fn])),  s=35, vmin=0, vmax= len(fx0[m + fn]))
                plt.xlim( (np.around(xpos[k]) - 2, np.around(xpos[k]) + 3 ) )
                plt.ylim( (np.around(ypos[k]) - 2, np.around(ypos[k]) + 3 ) )  


                plt.subplot(6,14, 14*m+9)
                plt.title('F1 S/B MaxE' , fontsize=8, color='DeepPink')
                plt.tick_params(labelsize=7)
                plt.scatter(fx1[m + fn], fy1[m + fn], marker = 'x', c = range(len(fx1[m + fn])),  s=35, vmin=0, vmax= len(fx1[m + fn]))
                plt.xlim( (np.around(xpos[k]) - 2, np.around(xpos[k]) + 3 ) )
                plt.ylim( (np.around(ypos[k]) - 2, np.around(ypos[k]) + 3 ) )  

                plt.subplot(6,14, 14*m+10)
                plt.title('F2 + Max%, MaxPx ' , fontsize=8, color='DeepPink', weight='bold')
                plt.tick_params(labelsize=7)
                plt.scatter(fx2[m + fn], fy2[m + fn], marker = 'x', c = range(len(fx2[m + fn])),  s=35, vmin=0, vmax= len(fx2[m + fn]))
                plt.xlim( (np.around(xpos[k]) - 2, np.around(xpos[k]) + 3 ) )
                plt.ylim( (np.around(ypos[k]) - 2, np.around(ypos[k]) + 3 ) )  


                
                plt.subplot(6,14, 14*m+11)
                plt.title(' AFg2 + Max%, MaxPx'  , fontsize=8, color="ForestGreen", weight='bold')
                plt.tick_params(labelsize=7)
                plt.scatter(gx2[m + fn], gy2[m + fn], marker = 'x', c = range(len(gx2[m + fn])),  s=35, vmin=0, vmax= len(gx2[m + fn]))
                plt.xlim( (np.around(xpos[k]) - 2, np.around(xpos[k]) + 3 ) )
                plt.ylim( (np.around(ypos[k]) - 2, np.around(ypos[k]) + 3 ) )  
                
                
                              
                plt.subplot(6,14, 14*m+12)
                plt.title('AFg1 S/B MaxE'  , fontsize=8, color="ForestGreen", weight='normal')
                plt.tick_params(labelsize=7)
                plt.scatter(gx1[m + fn], gy1[m + fn], marker = 'x', c = range(len(gx1[m + fn])),  s=35, vmin=0, vmax= len(gx1[m + fn]))
                plt.xlim( (np.around(xpos[k]) - 2, np.around(xpos[k]) + 3 ) )
                plt.ylim( (np.around(ypos[k]) - 2, np.around(ypos[k]) + 3 ) )  
                

                

                plt.subplot(6,14, 14*m+13)
                plt.title('AFg0  S/B'  , fontsize=8, color = 'ForestGreen')
                plt.tick_params(labelsize=7)
                plt.scatter(gx0[m + fn], gy0[m + fn], marker = 'x', c = range(len(gx0[m + fn])),  s=35, vmin=0, vmax= len(gx0[m + fn]))
                plt.xlim( (np.around(xpos[k]) - 2, np.around(xpos[k]) + 3 ) )
                plt.ylim( (np.around(ypos[k]) - 2, np.around(ypos[k]) + 3 ) )
                
                
                

                
                plt.subplot(6,14, 14*m+14)
                Ngx1temp = len(gx1[m + fn])
                Ngx2temp = len(gx2[m + fn])
                Ngx0temp = len(gx0[m + fn])
                Ngecctemp = len(geccentricityData[m + fn])                
                try:
                    R1 = np.around(100 * float(Ngx1temp)/float(Ngecctemp), 0)
                except:
                    R1 = 999                    
                try:
                    R2 = np.around(100 * float(Ngx2temp)/float(Ngx0temp), 0)
                except:
                    R2 = 999                
                plt.title('%d:%d:%d:%d,  %d:%d' %(Ngx2temp, Ngx1temp, Ngx0temp, Ngecctemp, R2, R1)  , fontsize=6)
                plt.tick_params(labelsize=6)
                
                try:
                    plt.hist(geccentricityData[m + fn])
                    plt.xlim(0, 1.2)
                    
                except:
                    print 'hist error'
                    print 'geccentricityData[m + fn]: ', geccentricityData[m + fn]
                    
                    
                    
                self.fig_AllDataCenterImages[n].text(0.88, 0.95, 'Histogram info\n #g2 : #g1 : #g0 : #TotalEcc,  #g2/#TEcc[%] : #g1/#TEcc[%]', ha="center", va="bottom", size=8,color="red")
                






                
                k += 1
                if k%6 ==0:
                    fn += 6
                if k == NFeatures:
                    break
         
         
         
         
        self.text_SaveFig.SetLabel('done')          
         
        print '\n AF gshrimp started'
        print startTime

        print '\n AF gshrimp done' 
        print time.strftime("%Y-%m%d, %Hh %Mm")
        
        
        endTimeDay = int(time.strftime("%d"))
        endTimeHour = int(time.strftime("%H"))
        endTimeMinute = int(time.strftime("%M"))
        
        totalDay = endTimeDay - startTimeDay
        totalHour = endTimeHour - startTimeHour
        totalMinute = endTimeMinute - startTimeMinute
        
        totalTimeInHour = totalDay * 24.0 + totalHour + totalMinute / 60.0
        
        print '\nTotal Time in Hour = ', totalTimeInHour
        
        
        
        if self.cb_show_afSHRImP_resultFigures_yn.GetValue():
            plt.ion()
            plt.show()
        else: pass
  


        self.btn_Save_figures.Show()
        self.btn_Save_gdata.Show()



    def SaveIntensityTrajectoryPlotsOnly(self, event):
        prefix = self.FilePrefix.GetValue()
        figFileName = prefix + self.filenameOnly[:30]
        todayDate = time.strftime("%Y%m%d")


        try:
            self.fig_FindingFeatures.savefig(figFileName +'_' + todayDate +  '_IntensityPlots_00Features.png')
            for n in range(len(self.fig_IntensityTrajectory)):
                self.fig_IntensityTrajectory[n].savefig(figFileName + '_' + todayDate + '_IntensityPlots_' + str(n)+ '.png')
                #plt.close(self.fig_IntensityTrajectory[n])
            print '\nIntensity Trajectory figures are saved'
            wx.StaticText(self.MainPanel, -1, 'Intensity Trajectory figures saved', pos=(450, 380)) 
            
        except:
            print '\nIntensity Trajectory figures are NOT saved'     
            wx.StaticText(self.MainPanel, -1, 'Error', pos=(450, 380)) 
            




      
    def SaveFigures(self, event):
        prefix = self.FilePrefix.GetValue()
        if prefix != '':
            prefix += '_'
        figFileName = prefix + self.filenameOnly[:35]
        todayDate = time.strftime("%Y%m%d_%Hh")
        print "\nToday's Date: ", time.strftime("%Y-%m-%d_%Hh")
        
        
             
        save_fig_error = 0
            
            
        if self.cb_Save_figures_AllInfoOnly.GetValue() == False:
            
            self.save_fig_ShowAfewFrames.savefig(figFileName + '_' + todayDate + '_1_RawImageFrames.png')
            print 'Raw Frames figures are saved'
            
            try:
                self.fig_FindingFeatures.savefig(figFileName +'_' + todayDate +  '_2_Molecules.png')
                print 'Feature finding figure is saved'
            except:
                print 'No finding feature figure'
                save_fig_error += 1
             
             

    
            try:
                for n in range(len(self.fig_IntensityTrajectory)):
                    self.fig_IntensityTrajectory[n].savefig(figFileName + '_' + todayDate + '_3_Intensity_' + str(n)+ '.png')
                    plt.close(self.fig_IntensityTrajectory[n])
                print 'Intensity Trajectory figures are saved'
            except:
                print 'Intensity Trajectory figures are NOT saved'      
                save_fig_error += 1
                
                
            try:        
                for n in range(len(self.fig_FionaDataAll)):    
                    self.fig_FionaDataAll[n].savefig(figFileName + '_' + todayDate + '_4_FIONA0_' + str(n) + '.png')
                    plt.close(self.fig_FionaDataAll[n])
                print 'Fiona 0 Data figures are saved'
                
            except:
                print 'Fiona 0 Data figures are NOT saved'  
                save_fig_error += 1
                

            try:        
                for n in range(len(self.fig_FionaData1)):    
                    self.fig_FionaData1[n].savefig(figFileName + '_' + todayDate + '_4_FIONA1_' + str(n) + '.png')
                    plt.close(self.fig_FionaData1[n])
                print 'Fiona 1 Data figures are saved'
                
            except:
                print 'Fiona 1 Data figures are NOT saved'  
                save_fig_error += 1
                
                
            try:        
                for n in range(len(self.fig_FionaData2)):    
                    self.fig_FionaData2[n].savefig(figFileName + '_' + todayDate + '_4_FIONA2_' + str(n) + '.png')
                    plt.close(self.fig_FionaData2[n])
                print 'Fiona 2 Data figures are saved'
                
            except:
                print 'Fiona 2 Data figures are NOT saved'  
                save_fig_error += 1
                
                
            try:        
                for n in range(len(self.fig_gshrimp1)):    
                    self.fig_gshrimp1[n].savefig(figFileName + '_' + todayDate + '_5_gshrimp1_' + str(n)+ '.png')
                    self.fig_gshrimp2[n].savefig(figFileName + '_' + todayDate + '_5_gshrimp2_' + str(n)+ '.png')
                    plt.close(self.fig_gshrimp1[n])
                    plt.close(self.fig_gshrimp2[n])
                print 'gshrimp figures are saved'
                
            except:
                print 'gshrimp figures are NOT saved'   
                save_fig_error += 1
            
            
            
            try:        
                for n in range(len(self.fig_geccentricityData)):    
                    self.fig_geccentricityData[n].savefig(figFileName + '_' + todayDate + '_6_Eccen_' + str(n)+ '.png')
                    plt.close(self.fig_geccentricityData[n])
                print 'geccentricityData figures are saved'
                
            except:
                print 'geccentricityData figures are NOT saved'   
                save_fig_error += 1
                

            
        if self.cb_Save_figures_AllInfoOnly.GetValue() == True:
            
            try:
                for n in range(len(self.fig_IntensityTrajectory)):
                    plt.close(self.fig_IntensityTrajectory[n])
            except:
                save_fig_error += 1
                
                
            try:        
                for n in range(len(self.fig_FionaDataAll)):    
                    plt.close(self.fig_FionaDataAll[n])
                                
            except:
                save_fig_error += 1
                

            try:        
                for n in range(len(self.fig_FionaData1)):    
                    plt.close(self.fig_FionaData1[n])
                                
            except:
                    save_fig_error += 1
                
                
            try:        
                for n in range(len(self.fig_FionaData2)):    
                    plt.close(self.fig_FionaData2[n])
                                
            except:
                print 'Fiona 2 Data figures are NOT saved'  
                save_fig_error += 1
                
                
            try:        
                for n in range(len(self.fig_gshrimp1)):    
                    plt.close(self.fig_gshrimp1[n])
                    plt.close(self.fig_gshrimp2[n])
                                
            except:
                save_fig_error += 1
            
            
            
            try:        
                for n in range(len(self.fig_geccentricityData)):    
                    plt.close(self.fig_geccentricityData[n])
                                
            except:
                    save_fig_error += 1
                



        try:        
            for n in range(len(self.fig_AllDataCenterImages)):    
                self.fig_AllDataCenterImages[n].savefig(figFileName + '_' + todayDate + '_7_ALLInfo_' + str(n)+ '.png')
                plt.close(self.fig_AllDataCenterImages[n])
            print 'All-Info figures are saved'
            
        except:
            print 'All-Info figures are NOT saved'       
            save_fig_error += 1
            
        
        
        
        
        
        if save_fig_error == 0 :
            save_fig_text = 'All figures are saved'
            print '\n', save_fig_text, '\n'
        else:
            save_fig_text = 'NOT All figures are saved'
            print '\n NOT All figures are saved\n # error: ', save_fig_error
            
        
        self.text_SaveFig.SetLabel(save_fig_text)
        
     
     
     
    
    
    
    def Save_gdata(self, event):
        prefix = self.FilePrefix.GetValue()
        if prefix != '':
            prefix += '_'
        
        todayDate = time.strftime("_%Y%m%d_%Hh")
        print "\nToday's Date: ", time.strftime("%Y-%m-%d_%Hh")
        
        gdataFileName = prefix + self.filenameOnly[:35] + todayDate + '_AFgshrimp' 

        for n in range(len(self.gdata)):
            #print 'len(self.gdata)[n] ', n, '  ', len(self.gdata[n])
            tempStr = self.gdata[n][0].split()
            #xytemp = '(' + str(int(np.around(float(tempStr[3]) ))) +',' + str(int(np.around(float(tempStr[4]) ))) +')'
            xytemp = '(' + str(int(float(tempStr[3]) )) +',' + str(int(float(tempStr[4]) )) +')'
            
            ff = open(gdataFileName + str(n) + '_' + xytemp + '.txt','w')
            #ff = open(gdataFileName + '_M#' +str(n) + '_' + xytemp + '.txt','w')
            for i in range(len(self.gdata[n])):
                
                ff.write(str(self.gdata[n][i]) + '\n'   )    
                
            ff.close()



        
        MdataFionaFileName = prefix + self.filenameOnly[:35] + todayDate + '_FIONA' 

        for n in range(len(self.data_FIONA)):
            tempStr = self.gdata[n][0].split()
            xytemp = '(' + str(int(np.around(float(tempStr[3]) ))) +',' + str(int(np.around(float(tempStr[4]) ))) +')'
               
            ff = open(MdataFionaFileName  + str(n) + '_' + xytemp + '.txt','w')
            #ff = open(MdataFionaFileName  + '_M#' +str(n) + '_' + xytemp + '.txt','w')
            
            for i in range(len(self.data_FIONA[n])):
                
                ff.write(str(self.data_FIONA[n][i]) + '\n'   )    
                
            ff.close()       

        
       


        MdataIntensityFileName = prefix + self.filenameOnly[:35] + todayDate + '_Intensity' 

        for n in range(len(self.IntensityTrajectoryMax)):
            tempStr = self.gdata[n][0].split()
            xytemp = '(' + str(int(np.around(float(tempStr[3]) ))) +',' + str(int(np.around(float(tempStr[4]) ))) +')'
            
            
            ff = open(MdataIntensityFileName  + str(n) + '_' + xytemp + '.txt','w')
            #ff = open(MdataIntensityFileName  + '_M#' +str(n) + '_' + xytemp + '.txt','w')
            
            for i in range(len(self.IntensityTrajectoryMax[n])):
                
                ff.write(str(self.IntensityTrajectoryMax[n][i]) + '\n'   )    
                
            ff.close()       

        

        
        
        print '\nAll data files are saved\n'     
        
                
        self.text_SaveGdata.SetLabel("Data Saved")
        
        
         
         
    def PlotSingleFeatureImages(self, event):     
        self.btn_PlotSingleImage_colorRangeUpdate.Show()
                
        frames = self.LoadedFrames
        
        N_xPixel = frames.frame_shape[0]
        #N_yPixel = frames.frame_shape[1]
                 
        
        #Nframes = frames._count
        Nframes = len(frames) #trackpy 0.3
        
        self.Nframes = Nframes
        
        ################################################
        x_ImageJ, y_ImageJ = int(self.xposSingle1.GetValue()), int(self.yposSingle1.GetValue()  )
        
        xc = y_ImageJ
        yc = N_xPixel - x_ImageJ + 1
        ################################################        
        
        
        
        NofSubpixel = 7  # number of pixels from the center
        
        x3 = np.linspace(0, 2 * NofSubpixel, 2* NofSubpixel + 1)
        y3 = np.linspace(0, 2 * NofSubpixel, 2*NofSubpixel+1)
        x3, y3 = np.meshgrid(x3, y3)
        
        
        if int(self.xposRawImageEndFrameN.GetValue()) > self.LastFrameN:
            self.xposRawImageEndFrameN.SetValue(str(self.LastFrameN))
        
        
        Nimages = int(self.xposRawImageEndFrameN.GetValue()) - int(self.xposRawImageStartFrameN.GetValue()) + 1
        
        
        NfeaturesPerFig = 160
        
        NFeatures = Nimages
        Nfig = int(  math.ceil((NFeatures/float(NfeaturesPerFig))))

        self.fig_RawImageData = [[] for _ in xrange(Nfig)]
        
        
        
        self.singleImageFrameN = []
        self.singleImagePlotData = []
        self.singleImagePlotData_imshowObject = []
        
        k = 0
        fn = 0
        for n in range(Nfig):
            print 'Raw_Image_Data_ Fig # = ', n
            self.fig_RawImageData[n] = plt.figure('Raw_Image_Data_'+ str(n), figsize = (18, 9))
            plt.clf()
            #fig_ShowAfewFrames.suptitle('test title', fontsize=20)
            #self.fig_IntensityTrajectory[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.28, hspace = 0.4)
            plt.figtext(0.5, 0.97, str(self.filedirectoryOnly)+'\n' + str(os.path.basename(self.filepath) + ',     Total frame # ' + str(Nframes))
                ,ha='center', color='black', weight='normal', size='small')
                
            self.fig_RawImageData[n].text(0.05, 0.96, 'x, y : (' + str(x_ImageJ)+', '+str(y_ImageJ)+')', ha="left", va="bottom", size="medium",color="red")
            
            colorMin = int(self.singleImageColorMin.GetValue())
            colorMax = int(self.singleImageColorMax.GetValue())

            for m in range(NfeaturesPerFig):
                
                print 'frame #: ', fn+m



                plt.subplot(8, 20,m + 1, aspect='equal')
                plt.axis('off')
                plt.subplots_adjust(top = 0.92, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.2, hspace = 0.4)
                plt.title('F # ' + str(fn+m),  fontsize=10)
                self.singleImageFrameN.append(fn+m)
                #plt.hold(True)
                
                
                if self.radio_NAvergeFrames_rawImage.GetSelection() == 0:
                    ftemp  = frames[fn+m][yc-NofSubpixel:yc+NofSubpixel+1, xc-NofSubpixel:xc+NofSubpixel+1]
                
                elif self.radio_NAvergeFrames_rawImage.GetSelection() == 1:
                    ftemp  = (frames[fn+m][yc-NofSubpixel:yc+NofSubpixel+1, xc-NofSubpixel:xc+NofSubpixel+1] +\
                            frames[fn+m-1][yc-NofSubpixel:yc+NofSubpixel+1, xc-NofSubpixel:xc+NofSubpixel+1] +\
                            frames[fn+m+1][yc-NofSubpixel:yc+NofSubpixel+1, xc-NofSubpixel:xc+NofSubpixel+1])/3
                    
                self.singleImagePlotData.append(ftemp)
                
                
                imTemp = plt.imshow(np.rot90(ftemp.reshape(2*NofSubpixel+1, 2*NofSubpixel+1) ,3) , interpolation = 'None', cmap=plt.cm.jet, origin='upper', vmin=colorMin, vmax=colorMax,
                    extent=(x3.min(), x3.max() +1, y3.min(), y3.max()+1 ))
                    
                self.singleImagePlotData_imshowObject.append(imTemp)

                                  
                k += 1
                if k%NfeaturesPerFig ==0:
                    fn += NfeaturesPerFig
                if k == NFeatures:
                    break
        
        self.btn_PlotSingleSave.Show()
        plt.show()
    





         
    def PlotSingleFeatureImagesUpdate(self, event):     
        print '\nPlotSingleFeatureImagesUpdate'
 
        colorMin = int(self.singleImageColorMin.GetValue())
        colorMax = int(self.singleImageColorMax.GetValue())
        
        print '\ncolor Min ', colorMin
        print 'color Max ', colorMax
        
      
        for n in range(len(self.singleImagePlotData_imshowObject)):
            self.singleImagePlotData_imshowObject[n].set_clim(vmin=colorMin, vmax=colorMax)
        
        for n in range(len(self.fig_RawImageData)):    
            self.fig_RawImageData[n].canvas.draw()
        






    def PlotSingleFeatureImagesSave(self, event):
        prefix = self.FilePrefix.GetValue()
        if prefix != '':
            prefix += '_'
        figFileName = prefix + self.filenameOnly[:40]
        #todayDate = time.strftime("%Y%m%d")

        xyStr = '(' + self.xposSingle1.GetValue() +','+ self.yposSingle1.GetValue() +')'

        try:
            for n in range(len(self.fig_RawImageData)):
                self.fig_RawImageData[n].savefig(figFileName + '_' + 'Image_'+ xyStr + '_' + str(n)+ '.png')
                #self.fig_RawImageData[n].savefig(figFileName + '_' + todayDate +xyStr+ '_RawImage' + str(n)+ '.png')
                #plt.close(self.fig_IntensityTrajectory[n])
            print '\nfig_RawImageData figures are saved'
            wx.StaticText(self.MainPanel, -1, 'Raw Image Data figures saved', pos=(750, 430)) 
            
        except:
            print '\nRaw Image Data figures are NOT saved'     
            wx.StaticText(self.MainPanel, -1, 'Error', pos=(750, 430)) 
            







    def SingleTwoFrameShrimp(self, event):
        print '\n Starting SingleTwoFrameShrimp'
        
        if self.radio_ColorMap.GetSelection() == 0:
            colormap = cm.jet
        elif self.radio_ColorMap.GetSelection() == 1:
            colormap = cm.gray
        elif self.radio_ColorMap.GetSelection() == 2:
            colormap = cm.rainbow        
        
        noiseAmp = int(self.xposSingle1NoiseAmp.GetValue())
        
        frame1 = np.array(self.LoadedFrames[int(self.xposSingle1FrameStart1.GetValue())])*0.0
        k = 0
        for n in range(int(self.xposSingle1FrameStart1.GetValue()), int(self.xposSingle1FrameStart2.GetValue())+1):
            frame1 += np.array(self.LoadedFrames[n])
            k += 1
        frame1 /= k
        
        image1Frames_str = 'F#: ' + self.xposSingle1FrameStart1.GetValue() + ' - ' + self.xposSingle1FrameStart2.GetValue()
            
        
        
        frame2 = np.array(self.LoadedFrames[int(self.xposSingle1FrameEnd1.GetValue())])*0.0
        k = 0
        for n in range(int(self.xposSingle1FrameEnd1.GetValue()), int(self.xposSingle1FrameEnd2.GetValue())+1):
            frame2 += np.array(self.LoadedFrames[n])
            k += 1
        frame2 /= k
        image2Frames_str = 'F#: ' + self.xposSingle1FrameEnd1.GetValue() + ' - ' + self.xposSingle1FrameEnd2.GetValue()
            

        
        SingleTwoFrameShrimp_Function(int(self.xposSingle1.GetValue()), int(self.yposSingle1.GetValue()),  frame1 , frame2, self.figFileInfo_text, image1Frames_str, image2Frames_str , colormap, noiseAmp )







    def SingleTwoFrameShrimpAdjacent(self, event):
        print '\n Starting SingleTwoFrameShrimp Adjacent'
        
        if self.radio_ColorMap.GetSelection() == 0:
            colormap = cm.jet
        elif self.radio_ColorMap.GetSelection() == 1:
            colormap = cm.gray
        elif self.radio_ColorMap.GetSelection() == 2:
            colormap = cm.rainbow
        
        
        noiseAmp = int(self.xposSingle1NoiseAmp.GetValue())
        
        
        print 'self.radio_ColorMap.GetSelection() = ', self.radio_ColorMap.GetSelection()
        print 'colormap = ', colormap
        
        Nadjacent = int(self.xposSingle1AdjacentFrameN.GetValue())
        
        frame1 = np.array(self.LoadedFrames[int(self.xposSingle1FrameStart1Adjacent.GetValue())])*0.0
        k = 0
        for n in range(int(self.xposSingle1FrameStart1Adjacent.GetValue())-Nadjacent, int(self.xposSingle1FrameStart1Adjacent.GetValue())+Nadjacent+1):
            print 'frame1 # = ', n
            frame1 += np.array(self.LoadedFrames[n])
            k += 1
        frame1 /= k
        
        image1Frames_str = 'F#: ' + str(int(self.xposSingle1FrameStart1Adjacent.GetValue())-Nadjacent) + ' - ' + str(int(self.xposSingle1FrameStart1Adjacent.GetValue())+Nadjacent)
            
        
        
        
        
        frame2 = np.array(self.LoadedFrames[int(self.xposSingle1FrameEnd1Adjacent.GetValue())])*0.0
        k = 0
        for n in range(int(self.xposSingle1FrameEnd1Adjacent.GetValue())-Nadjacent, int(self.xposSingle1FrameEnd1Adjacent.GetValue())+Nadjacent+1):
            print 'frame2 # = ', n
            frame2 += np.array(self.LoadedFrames[n])
            k += 1
        frame2 /= k
        
        image2Frames_str = 'F#: ' + str(int(self.xposSingle1FrameEnd1Adjacent.GetValue())-Nadjacent) + ' - ' + str(int(self.xposSingle1FrameEnd1Adjacent.GetValue())+Nadjacent)
            
        
         



        
        SingleTwoFrameShrimp_Function(int(self.xposSingle1.GetValue()), int(self.yposSingle1.GetValue()),  frame1 , frame2, self.figFileInfo_text, image1Frames_str, image2Frames_str , colormap , noiseAmp)









    def IntensityTrajectoryPlotsSingle(self, event):
        
        try:
            for n in range(len(self.fig_IntensityTrajectory)):
                plt.close(self.fig_IntensityTrajectory[n])
    
            for n in range(len(self.btn_FindingMolecules.Hide())):
                plt.close(self.btn_FindingMolecules[n])
            
            for n in range(len(self.fig_FionaDataAll)):
                plt.close(self.fig_FionaDataAll[n])
                
            for n in range(len(self.fig_FionaData1)):
                plt.close(self.fig_FionaData1[n])      
                
            for n in range(len(self.fig_FionaData2)):
                plt.close(self.fig_FionaData2[n])
                
            for n in range(len(self.fig_gshrimp1)):
                plt.close(self.fig_gshrimp1[n])
                
            for n in range(len(self.fig_gshrimp2)):
                plt.close(self.fig_gshrimp2[n])
                
            for n in range(len(self.fig_geccentricityData)):
                plt.close(self.fig_geccentricityData[n])
                
            for n in range(len(self.fig_AllDataCenterImages)):
                plt.close(self.fig_AllDataCenterImages[n])
            
        except:
            pass
        
        
        print '\nIntensityTrajectoryPlots Singles'
        
        xposImageJtemp = []
        yposImageJtemp = []
        
        for n in range(len(self.xposSingleM)):
            xtemp = int(self.xposSingleM[n].GetValue())
            ytemp = int(self.yposSingleM[n].GetValue())
            if (xtemp != 0) & (ytemp != 0):
                xposImageJtemp.append(xtemp)
                yposImageJtemp.append(ytemp)
                
        self.xposImageJ = np.array(xposImageJtemp)
        self.yposImageJ = np.array(yposImageJtemp)
        
        
        print 'self.xposImageJSingle: ', self.xposImageJ
        print 'self.yposImageJSingle: ', self.yposImageJ
        
        
        
        
        self.IntensityTrajectoryAve, self.IntensityTrajectoryMax, self.IntensityTrajectoryMax5Ave = IntensityTrajectoryPlotsFunc(self.filepath, self.xposImageJ, self.yposImageJ)
        
        
        AnalysisConditions = ( 'Chosen F # ' + str(self.Nframe.GetValue()) + '   Ave # F: ' + str(self.NofFramesForAveraging.GetValue()) + 
        '   f size: '+ str(self.FeatureSize.GetValue()) +  '   Min I: ' + str(self.MinIntensity.GetValue()) )
        
        


        
        NFeatures = len(self.IntensityTrajectoryMax)
        Nfig = int(  math.ceil((NFeatures/30.0)))
        
        TotalFeatureN = 'Total # ' + str(NFeatures)
        
        #self.fig_IntensityTrajectory = [[]] * Nfig
        self.fig_IntensityTrajectory = [[] for _ in xrange(Nfig)]
        
        k = 0
        fn = 0
        for n in range(Nfig):
            print 'Intensity Trajectory Nfig n = ', n
            self.fig_IntensityTrajectory[n] = plt.figure('IntensityTrajectory_'+ str(n), figsize = (18, 9))
            #fig_ShowAfewFrames.suptitle('test title', fontsize=20)
            self.fig_IntensityTrajectory[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.25, hspace = 0.35)
            plt.figtext(0.5, 0.95, str(self.filedirectoryOnly)+'\n' + str(os.path.basename(self.filepath) + '\n' + AnalysisConditions   )
                ,ha='center', color='black', weight='normal', size='small')
            self.fig_IntensityTrajectory[n].text(0.05, 0.96, TotalFeatureN + '    Max 5 ave', ha="left", va="bottom", size="medium",color="red")
            
            
            for m in range(30):
                #print 'subplot # = ', m
                plt.subplot(6,5, m+1)
                plt.tick_params(labelsize=9)
                plt.title('# ' + str(m + fn) + ',  (' + str(int(self.xposImageJ[m+fn])) + ', ' + str(int(self.yposImageJ[m+fn]))+')'  , fontsize=9)
                #plt.plot(self.IntensityTrajectoryAve[m + fn])
                #plt.plot(self.IntensityTrajectoryMax[m + fn])
                plt.plot(self.IntensityTrajectoryMax5Ave[m + fn])
                
                k += 1
                if k%30 ==0:
                    fn += 30
                if k == NFeatures:
                    break
                
                

        
        
        #print 'self.fig_IntensityTrajectory[0]', self.fig_IntensityTrajectory[0]
        #print 'self.fig_IntensityTrajectory[1]', self.fig_IntensityTrajectory[1]
        plt.show()
        
        
        self.btn_gshrimp.Show()
        self.btn_SaveIntensityTrajectoryPlots.Show()
        
        
        




    
    
# End.














'''
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
'''







def DataProcessing_FIONA_AFgshrimp(file_path, NDatalinesToSkip, driftData):
    print '\n\nLoading FIONA & AFGSHRIMP data files'
    print 'file_path: ', file_path
    





    AFshrimp_filepath = file_path
    
    fdata = open(AFshrimp_filepath, 'r')
    linestemp = fdata.readlines()
    fdata.close()
    
    
    #print 'linestemp \n', linestemp
    
    print ' linestemp[0] ', linestemp[0]
    
    sptemp = linestemp[0].split()
    
    print 'sptemp[3] = ', sptemp[3]
    MoleculeXpos = int(float(sptemp[3]))
    MoleculeYpos = int(float(sptemp[4]))
    
    try:
        AdjacentFramesN = int(float(sptemp[6]))
    except:
        AdjacentFramesN = 0
        
    
    
    print 'MoleculeXpos ,  MoleculeYpos', MoleculeXpos, MoleculeYpos
    
    
    DataLines = linestemp[NDatalinesToSkip:]
    
    #print 'DataLines \n', DataLines
    
    #f1n	f2n	centerIntensityMax[f1n]	centerIntensityMax[f2n]	Afit	x0fitImageJ	y0fitImageJ	z0fit	sigma_x	sigma_y	theta	eccentricity	:::	AfitErr	x0fitImageJErr	y0fitImageJErr	z0fitErr	sigma_xErr	sigma_yErr	thetaErr

    
    AFgshrimp_f1n = []
    AFgshrimp_f2n = []
    AFgshrimp_centerIntensityMaxf1n = []
    AFgshrimp_centerIntensityMaxf2n = []
    AFgshrimp_Afit = []
    AFgshrimp_x0fitImageJ = []
    AFgshrimp_y0fitImageJ = []
    AFgshrimp_z0fit = []
    AFgshrimp_sigma_x = []
    AFgshrimp_sigma_y = []
    AFgshrimp_theta = []
    AFgshrimp_eccentricity = []
    AFgshrimp_AfitErr	 = []
    AFgshrimp_x0fitImageJErr	 = []
    AFgshrimp_y0fitImageJErr	 = []
    AFgshrimp_z0fitErr	 = []
    AFgshrimp_sigma_xErr	 = []
    AFgshrimp_sigma_yErr	 = []
    AFgshrimp_thetaErr = []
    

    

    for line in DataLines:
        p = line.split()
        AFgshrimp_f1n.append(float(p[0]))
        AFgshrimp_f2n.append(float(p[1]))
        AFgshrimp_centerIntensityMaxf1n.append(float(p[2]))
        AFgshrimp_centerIntensityMaxf2n.append(float(p[3]))
        AFgshrimp_Afit.append(float(p[4]))
        AFgshrimp_x0fitImageJ.append(float(p[5]))
        AFgshrimp_y0fitImageJ.append(float(p[6]))
        AFgshrimp_z0fit.append(float(p[7]))
        AFgshrimp_sigma_x.append(float(p[8]))
        AFgshrimp_sigma_y.append(float(p[9]))
        AFgshrimp_theta.append(float(p[10]))
        AFgshrimp_eccentricity.append(float(p[11]))
        
        
        AFgshrimp_AfitErr.append(float(p[13]))
        AFgshrimp_x0fitImageJErr.append(float(p[14]))        
        AFgshrimp_y0fitImageJErr.append(float(p[15]))
        AFgshrimp_z0fitErr.append(float(p[16]))
        AFgshrimp_sigma_xErr.append(float(p[17]))
        AFgshrimp_sigma_yErr.append(float(p[18]))
        AFgshrimp_thetaErr.append(float(p[19]))
        
        
    '''
    print 'AFgshrimp_f1n: ',    AFgshrimp_f1n 
    print 'AFgshrimp_f2n: ',     AFgshrimp_f2n
    print 'AFgshrimp_centerIntensityMaxf1n: ',     AFgshrimp_centerIntensityMaxf1n
    print 'AFgshrimp_centerIntensityMaxf2n: ',     AFgshrimp_centerIntensityMaxf2n
    print 'AFgshrimp_x0fitImageJErr: ',     AFgshrimp_x0fitImageJErr
    print 'AFgshrimp_y0fitImageJErr: ',  AFgshrimp_y0fitImageJErr  
    print ': ',     
    print ': ',     
    '''





    FIONA_filepath =  file_path.replace('AFgshrimp', 'FIONA')
    
    fdata = open(FIONA_filepath, 'r')
    linestemp = fdata.readlines()
    fdata.close()
    
    
    #print 'linestemp \n', linestemp
    
    
    DataLines = linestemp[NDatalinesToSkip:]
    
    #print 'DataLines \n', DataLines
    
    #frame#,	x0fitFionaImageJ,	y0fitFionaImageJ,	Afit,	z0fit,	sigma_x,	sigma_y,	theta,	eccentricity,	BG_IntensityAve,	centerIntensityAve,	centerIntensityMax,	SBratioAve,	SBratioMax	:::	AfitErr	x0fitImageJErr	y0fitImageJErr	z0fitErr	sigma_xErr	sigma_yErr	thetaErr

    
    FIONA_FrameN = []
    FIONA_x0fit = []
    FIONA_y0fit = []
    FIONA_Afit = []
    FIONA_z0fit = []
    FIONA_sigma_x = []
    FIONA_sigma_y = []
    FIONA_theta = []
    FIONA_eccentricity = []
    FIONA_BG_IntensityAve = []
    FIONA_centerIntensityAve = []
    FIONA_centerIntensityMax = []
    FIONA_SBratioAve = []
    FIONA_SBratioMax = []
    FIONA_AfitErr = []
    FIONA_x0fitImageJErr = []
    FIONA_y0fitImageJErr = []
    FIONA_z0fitErr = []
    FIONA_sigma_xErr = []
    FIONA_sigma_yErr = []
    FIONA_thetaErr = []
    
    
    
    

    for line in DataLines:
        p = line.split()
        FIONA_FrameN.append(float(p[0]))
        FIONA_x0fit.append(float(p[1]))
        FIONA_y0fit.append(float(p[2]))
        FIONA_Afit.append(float(p[3]))
        FIONA_z0fit.append(float(p[4]))
        FIONA_sigma_x.append(float(p[5]))
        FIONA_sigma_y.append(float(p[6]))
        FIONA_theta.append(float(p[7]))
        FIONA_eccentricity.append(float(p[8]))
        FIONA_BG_IntensityAve.append(float(p[9]))
        FIONA_centerIntensityAve.append(float(p[10]))
        FIONA_centerIntensityMax.append(float(p[11]))
        FIONA_SBratioAve.append(float(p[12]))
        FIONA_SBratioMax.append(float(p[13]))
        FIONA_AfitErr.append(float(p[15]))
        FIONA_x0fitImageJErr.append(float(p[16]))
        FIONA_y0fitImageJErr.append(float(p[17]))
        FIONA_z0fitErr.append(float(p[18]))
        FIONA_sigma_xErr.append(float(p[19]))
        FIONA_sigma_yErr.append(float(p[20]))
        FIONA_thetaErr.append(float(p[21]))
        
        

    '''
    print 'FIONA_FrameN: ',    FIONA_FrameN 
    print 'FIONA_x0fit: ',     FIONA_x0fit
    print 'FIONA_y0fit: ',     FIONA_y0fit
    print 'FIONA_eccentricity: ',     FIONA_eccentricity
    print 'FIONA_AfitErr: ',     FIONA_AfitErr
    print 'FIONA_x0fitImageJErr: ',  FIONA_x0fitImageJErr  
    print ': ',     
    print ': ',     
    '''
    
    




    TotalIntensity_filepath = file_path.replace('AFgshrimp', 'Intensity')
    
    TotalIntensity_Trajectory = []
    
    try:            
        fdata = open(TotalIntensity_filepath, 'r')
        linestemp = fdata.readlines()
        fdata.close()
        DataLines = linestemp
    
        
        
        
        for line in DataLines:
            p = line.split()
            TotalIntensity_Trajectory.append(float(p[0]))
    except:
        print '\nNo intensity file'
        pass
    
    


    ####################################################################################################################
    ### applying dedrift
    
    if driftData == 'na':
        pass
    else:
        drift_Frame_N, drift_dx, drift_dy = driftData[0], driftData[1], driftData[2]

        print '\n\napplyting dedrift'



        fx_dedrifted = []
        fy_dedrifted = []
        for n in range(len(FIONA_x0fit)):
            FN_temp = FIONA_FrameN[n]
            f_index_temp = drift_Frame_N.index(FN_temp)
            dx_temp = drift_dx[f_index_temp]
            dy_temp = drift_dy[f_index_temp]
            
            
            fx_dedrifted.append(FIONA_x0fit[n] - dx_temp)
            fy_dedrifted.append(FIONA_y0fit[n] - dy_temp)
    
        '''
        plt.figure()
        plt.plot(drift_Frame_N, drift_dx)
        plt.plot(drift_Frame_N, drift_dy)
        
    
        plt.figure()
        plt.plot(FIONA_FrameN, FIONA_x0fit)
        plt.plot(FIONA_FrameN, fx_dedrifted)
        plt.show()
        exit()
        '''
        
    
    
        FIONA_x0fit = fx_dedrifted
        FIONA_y0fit = fy_dedrifted
        
        
        
        
        
        
        gx_dedrifted = []
        gy_dedrifted = []
        for n in range(len(AFgshrimp_f1n)):
            FN_temp1 = AFgshrimp_f1n[n]
            FN_temp2 = AFgshrimp_f2n[n]
            f1_index_temp = drift_Frame_N.index(FN_temp1)
            f2_index_temp = drift_Frame_N.index(FN_temp2)
            dx_temp1 = drift_dx[f1_index_temp]
            dy_temp1 = drift_dy[f1_index_temp]
            dx_temp2 = drift_dx[f2_index_temp]
            dy_temp2 = drift_dy[f2_index_temp]
            
            dx_ave_temp = (dx_temp1 + dx_temp2)/2.0
            dy_ave_temp = (dy_temp1 + dy_temp2)/2.0
            
            gx_dedrifted.append(AFgshrimp_x0fitImageJ[n] - dx_ave_temp)
            gy_dedrifted.append(AFgshrimp_y0fitImageJ[n] - dy_ave_temp)
    
    
        AFgshrimp_x0fitImageJ = gx_dedrifted
        AFgshrimp_y0fitImageJ = gy_dedrifted
        
    ### applying dedrift
    ####################################################################################################################
                    

    
    return (MoleculeXpos, MoleculeYpos, FIONA_FrameN, FIONA_x0fit, FIONA_y0fit, FIONA_Afit, FIONA_z0fit, FIONA_sigma_x, FIONA_sigma_y, FIONA_theta, FIONA_eccentricity
        , FIONA_BG_IntensityAve, FIONA_centerIntensityAve, FIONA_centerIntensityMax, FIONA_SBratioAve, FIONA_SBratioMax
        , FIONA_AfitErr, FIONA_x0fitImageJErr, FIONA_y0fitImageJErr, FIONA_z0fitErr, FIONA_sigma_xErr, FIONA_sigma_yErr, FIONA_thetaErr
        , AFgshrimp_f1n, AFgshrimp_f2n, AFgshrimp_centerIntensityMaxf1n, AFgshrimp_centerIntensityMaxf2n, AFgshrimp_Afit
        , AFgshrimp_x0fitImageJ, AFgshrimp_y0fitImageJ, AFgshrimp_z0fit, AFgshrimp_sigma_x, AFgshrimp_sigma_y, AFgshrimp_theta, AFgshrimp_eccentricity
        , AFgshrimp_AfitErr, AFgshrimp_x0fitImageJErr, AFgshrimp_y0fitImageJErr, AFgshrimp_z0fitErr, AFgshrimp_sigma_xErr, AFgshrimp_sigma_yErr, AFgshrimp_thetaErr
        , AFshrimp_filepath, TotalIntensity_Trajectory, AdjacentFramesN)
    
        






'''
#############################################################################################################################
#############################################################################################################################
'''


  
class MainTrackingDataFileAnalysis(wx.Frame):
    #ImageFilePath = ''
 
    #----------------------------------------------------------------------
    def __init__(self):
        wx.Frame.__init__(self, None, wx.ID_ANY,
                          "FIONA_AFSHRIMP_DATA_ANALYSIS_Plot", size=(1000, 1000), pos=(0,0))
                          
        self.MainPanel = wx.Panel(self, wx.ID_ANY)
        
        

        

        self.btn_SwitchToAFShrimp = wx.Button(self.MainPanel, pos=(300,10), label="Switch To afShrimp")
        self.btn_SwitchToAFShrimp.Bind(wx.EVT_BUTTON, self.SwitchToAFShrimp)
        #self.btn_SwitchToDataFileAnalysis = wx.Button(self.MainPanel, pos=(500,10), label="Switch To Data File Analysis")
        #self.btn_SwitchToDataFileAnalysis.Bind(wx.EVT_BUTTON, self.SwitchToDataFileAnalysis)
        self.btn_SwitchToIntensityTrajectoryPlots = wx.Button(self.MainPanel, pos=(700,10), label="Switch To Intensity Trajectory Plots with Images")
        self.btn_SwitchToIntensityTrajectoryPlots.Bind(wx.EVT_BUTTON, self.SwitchToIntensityTrajectoryPlots)



        #self.btn_restart_program = wx.Button(self.MainPanel, pos=(600,50), label="Restart Program")
        #self.btn_restart_program.Bind(wx.EVT_BUTTON, self.restart_program)



        #self.cb_JumpSizeOnly = wx.CheckBox(self.MainPanel, -1, 'Jump Size to File Only', (550, 55))
        #self.cb_JumpSizeOnly.SetValue(False)

        #self.cb_1st_last_distance = wx.CheckBox(self.MainPanel, -1, '1st_last pair distance, RMS', (550, 75))
        #self.cb_1st_last_distance.SetValue(False)


        self.radio_options_autoMultiPleAnalysis = wx.RadioBox(self.MainPanel, -1, choices=['Individual Plots', 'Jump Size to File Only', '1st_last pair distance, RMS', 'Nearest distance histogram'], label='Options for Auto-Analysis', pos=(560, 35), size=wx.Size(160, 100), style=wx.RA_SPECIFY_ROWS)
        self.radio_options_autoMultiPleAnalysis.SetSelection(0)
        self.cb_autoMultiPleAnalysisSaveData = wx.CheckBox(self.MainPanel, -1, 'Auto_Data_Save', (565, 140))
        self.cb_autoMultiPleAnalysisSaveData.SetValue(False)



        self.cb_apply_dedrift = wx.CheckBox(self.MainPanel, -1, 'Apply dedrift', (310, 165))
        self.cb_apply_dedrift.SetValue(False)

        self.btn_Load_drift_file = wx.Button(self.MainPanel, pos=(400,160), label="Load drift file")
        self.btn_Load_drift_file.Bind(wx.EVT_BUTTON, self.Load_drift_file)
        self.text_Load_drift_file = wx.StaticText(self.MainPanel, -1, 'na ', pos=(490, 165))

     
        
        
        

        



        self.btn_AutomaticMultipleDataAnalysis = wx.Button(self.MainPanel, pos=(730,50), label="Automatic Multiple Data Analysis")
        self.btn_AutomaticMultipleDataAnalysis.Bind(wx.EVT_BUTTON, self.AutomaticMultipleDataAnalysis)

        
        self.btn_Save_figures_AutomaticMultipleDataAnalysis = wx.Button(self.MainPanel, pos=(730,80), label="Save Auto-Adjacent Figures & Data")  
        self.btn_Save_figures_AutomaticMultipleDataAnalysis.Show()
        self.btn_Save_figures_AutomaticMultipleDataAnalysis.Bind(wx.EVT_BUTTON, self.SaveFigures_AutomaticMultipleDataAnalysis)
        self.text_Save_figures_AutomaticMultipleDataAnalysis = wx.StaticText(self.MainPanel, -1, ' ', pos=(730, 110))

     


        
        wx.StaticText(self.MainPanel, -1, 'Open a AF-Shrimp data file', pos=(50, 20))
        
        self.btn_OpenData = wx.Button(self.MainPanel, pos=(50,50), label="Open Data File")
        self.DataFilePath = self.btn_OpenData.Bind(wx.EVT_BUTTON, self.onOpenDataFile)
        wx.StaticText(self.MainPanel, -1, "# of lines to skip:", pos=(170, 53))
        self.NLinesToSkip = wx.TextCtrl(self.MainPanel, -1, "2", pos=(280, 50), size=(40,-1))
        self.btn_OpenData_Text = wx.StaticText(self.MainPanel, -1, str(self.DataFilePath), pos=(40, 85))
        
        self.cb_ShowPublicationFigures = wx.CheckBox(self.MainPanel, -1, 'Show Publication Figures', (50, 105))
        self.cb_ShowPublicationFigures.SetValue(False)
        
        





        self.btn_PlotData = wx.Button(self.MainPanel, pos=(50,150), label="Plot Data")
        self.btn_PlotData.Hide()
        self.btn_PlotData.Bind(wx.EVT_BUTTON, self.FIONA_AFSHRIMP_DATA_ANALYSIS_Plot)
        
        wx.StaticText(self.MainPanel, -1, "<<     Plot the data again to change the constraints for the emission sites analysis    >>", pos=(150, 190))
        
        
        wx.StaticText(self.MainPanel, -1, "Min S(ave)/B(ave) ", pos=(180, 212))
        self.SvsB_Tolerance = wx.TextCtrl(self.MainPanel, -1, "1.1", pos=(290, 210), size=(40,-1))
        wx.StaticText(self.MainPanel, -1, "Max eccentricity", pos=(420, 212))
        self.MaxEccentricity = wx.TextCtrl(self.MainPanel, -1, "0.55", pos=(520, 210), size=(40,-1))
        wx.StaticText(self.MainPanel, -1, "Max fit err %: A", pos=(180, 242))
        self.MaxPercentage_Tolerance = wx.TextCtrl(self.MainPanel, -1, "1000", pos=(290, 240), size=(40,-1))        
        wx.StaticText(self.MainPanel, -1, "Max fit err px: x, y", pos=(410, 242))
        self.MaxPixel_Tolerance = wx.TextCtrl(self.MainPanel, -1, "0.099", pos=(520, 240), size=(40,-1))   
        
        
        wx.StaticText(self.MainPanel, -1, "Frame sections", pos=(50, 272))
        
        
        self.F_SectionBegin = [[] for _ in xrange(10)]
        self.F_SectionEnd = [[] for _ in xrange(10)]
        
        self.F_SectionBegin[0] = wx.TextCtrl(self.MainPanel, -1, "7", pos=(50, 300), size=(40,-1))
        self.F_SectionEnd[0] = wx.TextCtrl(self.MainPanel, -1, "11", pos=(100, 300), size=(40,-1))
        
        self.F_SectionBegin[1] = wx.TextCtrl(self.MainPanel, -1, "14", pos=(200, 300), size=(40,-1))
        self.F_SectionEnd[1] = wx.TextCtrl(self.MainPanel, -1, "23", pos=(250, 300), size=(40,-1))
        
        self.F_SectionBegin[2] = wx.TextCtrl(self.MainPanel, -1, "24", pos=(350, 300), size=(40,-1))
        self.F_SectionEnd[2] = wx.TextCtrl(self.MainPanel, -1, "29", pos=(400, 300), size=(40,-1))
        
        self.F_SectionBegin[3] = wx.TextCtrl(self.MainPanel, -1, "31", pos=(500, 300), size=(40,-1))
        self.F_SectionEnd[3] = wx.TextCtrl(self.MainPanel, -1, "65", pos=(550, 300), size=(40,-1))
        
        self.F_SectionBegin[4] = wx.TextCtrl(self.MainPanel, -1, "81", pos=(650, 300), size=(40,-1))
        self.F_SectionEnd[4] = wx.TextCtrl(self.MainPanel, -1, "113", pos=(700, 300), size=(40,-1))
        
        
        self.F_SectionBegin[5] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(50, 330), size=(40,-1))
        self.F_SectionEnd[5] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(100, 330), size=(40,-1))
        
        self.F_SectionBegin[6] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(200, 330), size=(40,-1))
        self.F_SectionEnd[6] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(250, 330), size=(40,-1))
        
        self.F_SectionBegin[7] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(350, 330), size=(40,-1))
        self.F_SectionEnd[7] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(400, 330), size=(40,-1))
        
        self.F_SectionBegin[8] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(500, 330), size=(40,-1))
        self.F_SectionEnd[8] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(550, 330), size=(40,-1))
        
        self.F_SectionBegin[9] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(650, 330), size=(40,-1))
        self.F_SectionEnd[9] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(700, 330), size=(40,-1))
        
        
        
        
        
        wx.StaticText(self.MainPanel, -1, "Magnification for RMS (nm/pixel): ", pos=(200, 403))
        self.M_nm_pixel = wx.TextCtrl(self.MainPanel, -1, str(magnification), pos=(380, 400), size=(40,-1))

        
        
        self.btn_EmissionSitesAnalysis = wx.Button(self.MainPanel, pos=(50,400), label="Emission Sites Analysis")
        self.btn_EmissionSitesAnalysis.Hide()
        self.btn_EmissionSitesAnalysis.Bind(wx.EVT_BUTTON, self.EmissionSitesAnalysis)


        self.btn_Save_EmissionSitesAnalysisData = wx.Button(self.MainPanel, pos=(50,430), label="Save Emission Sites Analysis Data")  
        self.btn_Save_EmissionSitesAnalysisData.Show()
        self.btn_Save_EmissionSitesAnalysisData.Bind(wx.EVT_BUTTON, self.Save_EmissionSitesAnalysisData)
        self.btn_Save_EmissionSitesAnalysisData = wx.StaticText(self.MainPanel, -1, ' ', pos=(50, 470))
        
        

        
        wx.StaticText(self.MainPanel, -1, "Start frame #", pos=(570, 393))
        self.AdjacentShrimpStartFrameN = wx.TextCtrl(self.MainPanel, -1, str(AdjacentShrimp_Start_FrameN), pos=(660, 390), size=(40,-1))
        wx.StaticText(self.MainPanel, -1, "End frame #", pos=(570, 423))
        self.AdjacentShrimpEndFrameN = wx.TextCtrl(self.MainPanel, -1, str(AdjacentShrimp_End_FrameN), pos=(660, 420), size=(40,-1))
        
        wx.StaticText(self.MainPanel, -1, "Max # of adjacent \n    next frames", pos=(570, 463))
        self.MaxN_nextAdjacentFrames = wx.TextCtrl(self.MainPanel, -1, str(NadjacentFrameNs), pos=(680, 460), size=(40,-1))
        
        
        #self.cb_BG_ref_correction_yn = wx.CheckBox(self.MainPanel, -1, 'Intensity Correction by reference background', (570, 510))
        #self.cb_BG_ref_correction_yn.SetValue(False)
        
        wx.StaticText(self.MainPanel, -1, "Shrimp Ecc < FIONA Ecc", pos=(570, 500))
        self.radio_AdjacentEcc = wx.RadioBox(self.MainPanel, -1, choices=['None', 'Either', 'Both'], label='2nd ECC Constraint', pos=(570, 520), size=wx.Size(120, 80), style=wx.RA_SPECIFY_ROWS)
        self.radio_AdjacentEcc.SetSelection(0)
        
        #print 'self.radio_AdjacentEcc.GetSelection() ', self.radio_AdjacentEcc.GetSelection()
        
        
        self.radio_AdjacentBleachingType = wx.RadioBox(self.MainPanel, -1, choices=['Random', 'Random with the Sections', 'Monotonic', 'Monotonic - Sections', 'Monotonic - Sections (all point)'], label='Bleaching Type', pos=(570, 605), size=wx.Size(180, 80), style=wx.RA_SPECIFY_ROWS)
        self.radio_AdjacentBleachingType.SetSelection(2)
        
        
        
        
        wx.StaticText(self.MainPanel, -1, "# forward frames for Min Err", pos=(570, 723))
        self.Adjacent_N_foward_frames = wx.SpinCtrl(self.MainPanel, -1, pos=(720, 720), size=(50,-1))
        self.Adjacent_N_foward_frames.SetValue(3)
        
        
        wx.StaticText(self.MainPanel, -1, "Min # Frames: f1n diff", pos=(570, 753))
        self.Adjacent_f1n_min_diff = wx.SpinCtrl(self.MainPanel, -1, pos=(700, 750), size=(50,-1))
        self.Adjacent_f1n_min_diff.SetValue(3)
        
        wx.StaticText(self.MainPanel, -1, "Min # Frames: f1n f2n", pos=(570, 783))
        self.Adjacent_f1n_f2n_min_diff = wx.SpinCtrl(self.MainPanel, -1, pos=(700, 780), size=(50,-1))
        self.Adjacent_f1n_f2n_min_diff.SetValue(3)
        
        
        
  
        self.btn_AdjacentShrimpPlot = wx.Button(self.MainPanel, pos=(570,820), label="Plot Adjacent Shrimp Data")
        self.btn_AdjacentShrimpPlot.Show()
        self.btn_AdjacentShrimpPlot.Bind(wx.EVT_BUTTON, self.AdjacentShrimpPlot)
        self.text_AdjacentShrimpPlot = wx.StaticText(self.MainPanel, -1, ' ', pos=(570, 850))        
        
        
        self.btn_LastAdjacentShrimp_SaveFig = wx.Button(self.MainPanel, pos=(570,870), label="Save the Last\nAdjacent-Shrimp Figure")
        self.btn_LastAdjacentShrimp_SaveFig.Show()
        self.btn_LastAdjacentShrimp_SaveFig.Bind(wx.EVT_BUTTON, self.SaveFigures_LastAdjacentShrimp)
        self.text_LastAdjacentShrimp_Figsaved = wx.StaticText(self.MainPanel, -1, ' ', pos=(570, 910))        
        
        
        
        
        
        
        
        
        
        


        self.NTotalEmissionSites = 15
        self.EmissionSiteX = [[] for _ in xrange(self.NTotalEmissionSites)]
        self.EmissionSiteY = [[] for _ in xrange(self.NTotalEmissionSites)]
        self.EmissionSiteErrR = [[] for _ in xrange(self.NTotalEmissionSites)]
        EmissionSiteYaxisTemp = 520
        
        wx.StaticText(self.MainPanel, -1, "Enter Emission Sites \n x                  y", pos=(250, EmissionSiteYaxisTemp-30))
        wx.StaticText(self.MainPanel, -1, "\nErr radius", pos=(390, EmissionSiteYaxisTemp-30))
        for n in range(self.NTotalEmissionSites):
            wx.StaticText(self.MainPanel, -1, "# " + str(n+1), pos=(215, EmissionSiteYaxisTemp + 25*n))
            self.EmissionSiteX[n] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(250, EmissionSiteYaxisTemp + 25*n), size=(50,-1))
            self.EmissionSiteY[n] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(310, EmissionSiteYaxisTemp + 25*n), size=(50,-1))
            self.EmissionSiteErrR[n] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(390, EmissionSiteYaxisTemp + 25*n), size=(50,-1))





        self.btn_PlotFinalEmissionSites = wx.Button(self.MainPanel, pos=(50,650), label="Plot Final Emission Sites", size=(150,-1))
        self.btn_PlotFinalEmissionSites.Hide()
        self.btn_PlotFinalEmissionSites.Bind(wx.EVT_BUTTON, self.PlotFinalEmissionSites)
        
        
        





        
        wx.StaticText(self.MainPanel, -1, "Figure Files Prefix", pos=(50, 700))
        self.FilePrefix = wx.TextCtrl(self.MainPanel, -1, "", pos=(50, 720))
        
        
        
        self.btn_Save_figures = wx.Button(self.MainPanel, pos=(50,750), label="Save Figures")  
        self.btn_Save_figures.Hide()
        self.btn_Save_figures.Bind(wx.EVT_BUTTON, self.SaveFigures)
        self.text_Save_figures = wx.StaticText(self.MainPanel, -1, ' ', pos=(50, 780))



        self.driftData = 'na'

                     


    #----------------------------------------------------------------------#
    #----------------------------------------------------------------------#

    def restart_program(self, event):
        print 'Restarting Program'
        #restart_program()
        

    def SwitchToAFShrimp(self, event):
        frameSM2Tracking = MainTracking()
        frameSM2Tracking.Show()
        self.Close()
         

    def SwitchToDataFileAnalysis(self, event):
        frameDataFileAnalysis = MainTrackingDataFileAnalysis()
        frameDataFileAnalysis.Show()
        self.Close()
        

    def SwitchToIntensityTrajectoryPlots(self, event):
        IntensityTrajectoryPlots = MainTrackingIntensityOnlyPlots()
        IntensityTrajectoryPlots.Show()
        self.Close()
        



        

    def on_paint(self, event):
        dc = wx.PaintDC(event.GetEventObject())
        dc.Clear()
        dc.SetPen(wx.Pen("BLACK", 4))
        dc.DrawLine(0, 500, 800, 500)  
            
        
     #----------------------------------------------------------------------

    def SwitchToAFShrimp(self, event):
        frameSM2Tracking = MainTracking()
        frameSM2Tracking.Show()
        self.Close()
        
         
 
 

        
        
    def onOpenDataFile(self, event):
        
        self.text_LastAdjacentShrimp_Figsaved.SetLabel(' ')
        self.text_AdjacentShrimpPlot.SetLabel(' ')
        #self.MaxPixel_Tolerance.SetValue('0.12')
        
        plt.close('fxy2_distance')
        
        try:
            for n in range(len(self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_Plot)):
                plt.close(self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_Plot[n])
    
            for n in range(len(self.fig_EmissionSitesAnalysis)):
                plt.close(self.fig_EmissionSitesAnalysis[n])
        except:
            pass             
                
            
        try:
            for n in range(len(self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot)):
                plt.close(self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot[n])
        except: pass             


            
        try:
            for n in range(len(self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_Plot)):
                plt.close(self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_Plot[n])
        except: pass             
            
        try:
            for n in range(len(self.fig_EmissionSitesAnalysis)):
                plt.close(self.fig_EmissionSitesAnalysis[n])
        except: pass             
            
        try:
            for n in range(len(self.fig_EmissionSitesAnalysisFinal)):
                plt.close(self.fig_EmissionSitesAnalysisFinal[n])
        except: pass             


                
        
        try:
            plt.close(self.fig_afshrimpDensity)
        except:
            pass
        
        try:
            plt.close(self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_ACF_Plot)
        except:
            pass
                    
            
            

        dlg2 = wx.FileDialog(
            self, message="Choose a file",
            defaultFile="",
            wildcard = wildcardDataFile,
            style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
            )
        if dlg2.ShowModal() == wx.ID_OK:
            paths = dlg2.GetPaths()
            print "The following file(s) are chosen:"
            for path in paths:
                print path
        
        
        print ''
        print 'dlg2 = ', dlg2
        print 'dlg2.GetPaths() = ', dlg2.GetPaths()
        print ' path[0]', paths[0]
        #print ' path[1]', paths[1]
        
        print 'path type', type(dlg2.GetPath())
        print 'filename = ', dlg2.GetFilename()
        #dlg.Destroy()
        
        self.Show()
        
        #wx.StaticText(self.MainPanel, -1, str(dlg2.GetFilename()), pos=(50, 150))
        
        filenameTemp = dlg2.GetFilename()
        self.btn_OpenData_Text.SetLabel(str(filenameTemp))

        
        
        self.filepath = dlg2.GetPath()
        print 'self.filepath[0]: ', self.filepath[0]
        print 'self.filepath[1]: ', self.filepath[1]
        print 'dlg2.GetPath():  ', dlg2.GetPath()
        
        self.filenameOnly = dlg2.GetFilename()
        self.filedirectoryOnly = dlg2.GetDirectory()
        
        (self.MoleculeXpos, self.MoleculeYpos, self.FIONA_FrameN, self.FIONA_x0fit, self.FIONA_y0fit, self.FIONA_Afit, self.FIONA_z0fit, self.FIONA_sigma_x, self.FIONA_sigma_y, self.FIONA_theta, self.FIONA_eccentricity
        , self.FIONA_BG_IntensityAve, self.FIONA_centerIntensityAve, self.FIONA_centerIntensityMax, self.FIONA_SBratioAve, self.FIONA_SBratioMax
        , self.FIONA_AfitErr, self.FIONA_x0fitImageJErr, self.FIONA_y0fitImageJErr, self.FIONA_z0fitErr, self.FIONA_sigma_xErr, self.FIONA_sigma_yErr, self.FIONA_thetaErr
        , self.AFgshrimp_f1n, self.AFgshrimp_f2n, self.AFgshrimp_centerIntensityMaxf1n, self.AFgshrimp_centerIntensityMaxf2n, self.AFgshrimp_Afit
        , self.AFgshrimp_x0fitImageJ, self.AFgshrimp_y0fitImageJ, self.AFgshrimp_z0fit, self.AFgshrimp_sigma_x, self.AFgshrimp_sigma_y, self.AFgshrimp_theta, self.AFgshrimp_eccentricity
        , self.AFgshrimp_AfitErr, self.AFgshrimp_x0fitImageJErr, self.AFgshrimp_y0fitImageJErr, self.AFgshrimp_z0fitErr, self.AFgshrimp_sigma_xErr, self.AFgshrimp_sigma_yErr, self.AFgshrimp_thetaErr
        , self.AFshrimp_filepath, self.TotalIntensity_Trajectory, self.AdjacentFramesN) = DataProcessing_FIONA_AFgshrimp(paths[0], int(self.NLinesToSkip.GetValue()), self.driftData )
        
        #print 'FIONA_FrameN\n', self.FIONA_FrameN
        
        
        self.dataFile_path = paths[0]
        
        self.btn_PlotData.Show()
        
        
        try:
            self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_Plot[0].clf()
        except:
            print 'Failed: self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_Plot[0].clf()'
            pass
        
        

        
        
        self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_Plot = []
        self.fig_EmissionSitesAnalysis = []
        self.fig_EmissionSitesAnalysisFinal = []
        

        self.figCounter = 0
        self.AdjacentPlot_figCounter = 0
        self.fig_EmissionSitesAnalysisFinalCounter = 0
        

        self.fig_AllDataCenterImagesSingle = []
        self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot = []



        
        self.MoleculeXYtext = '(' + str(self.MoleculeXpos)+', '+str(self.MoleculeYpos) + ')' + '\n# of Averaged Frames: ' + str(self.AdjacentFramesN*2+1)
        
        
        self.FIONA_AFSHRIMP_DATA_ANALYSIS_Plot(event)
        
        






    def Load_drift_file(self, event):
        print 'Loading drift file'
        
        
        
        dlg2 = wx.FileDialog(
            self, message="Choose a file",
            defaultFile="",
            wildcard = wildcardDataFile,
            style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
            )
                    
        if dlg2.ShowModal() == wx.ID_OK:
            paths = dlg2.GetPaths()
            filesNames = dlg2.GetFilenames()
            print "The following file(s) are chosen:"
            print paths
            
            self.filename_short = filesNames[0][3:22] 
 
        self.text_Load_drift_file.SetLabel(filesNames[0])
        
        fdata = open(paths[0], 'r')
        linestemp = fdata.readlines()
        fdata.close()
        
        
        #print 'linestemp \n', linestemp
        for n in range(len(linestemp)):
            print linestemp[n]
            
        
        print ' linestemp[0] ', linestemp[0]
        
        sptemp = linestemp[0].split()
        
        print 'sptemp \n', sptemp
    

    
        self.drift_Frame_N = []
        self.drift_dx = []
        self.drift_dy = []
        for line in linestemp:
            p = line.split()
            self.drift_Frame_N.append(int(p[0]))
            self.drift_dx.append(float(p[1]))
            self.drift_dy.append(float(p[2]))
            
            
            
        print 'self.drift_Frame_N \n ', self.drift_Frame_N
        print 'self.drift_dx \n ', self.drift_dx
        print 'self.drift_dy \n ', self.drift_dy 



        self.driftData_Loaded = [self.drift_Frame_N, self.drift_dx, self.drift_dy]
        
        
        try:
            plt.close('drift-data')
        except: pass
    
        plt.figure('drift-data')
        
        plt.title('drift: dx, dy \n' + filesNames[0], size = 11)
        #plt.plot(self.drift_Frame_N, np.array(self.drift_dx)*float(self.M_nm_pixel.GetValue()), '-bo', label='dx')
        #plt.plot(self.drift_Frame_N, np.array(self.drift_dy)*float(self.M_nm_pixel.GetValue()), '-ro', label='dy')
        plt.plot(self.drift_Frame_N, np.array(self.drift_dx), '-bo', label='dx')
        plt.plot(self.drift_Frame_N, np.array(self.drift_dy), '-ro', label='dy')
        plt.xlabel('Frame #')
        plt.legend(loc='upper left', fontsize=11)
        #plt.ylim(-50,50)
        plt.ylabel('Distance (pixels)')
        
        plt.ion()
        plt.show()
    




        
    def AutomaticMultipleDataAnalysis(self, event):
        print '\n\n#######\nAutomatic Multiple Data Analysis \n#######'
        
        
        try:
            plt.close(self.Fig_rmsHistogram)
            plt.close('fxy2_distance')
            for n in range(len(self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot)):
                plt.close(self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot[n])
        except:
            pass
        
    
    
        
        self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_Plot = []
        self.fig_EmissionSitesAnalysis = []
        self.fig_EmissionSitesAnalysisFinal = []
        

        self.figCounter = 0
        self.AdjacentPlot_figCounter = 0
        self.fig_EmissionSitesAnalysisFinalCounter = 0
        

        self.fig_AllDataCenterImagesSingle = []
        self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot = []


        
        dlg2 = wx.FileDialog(
            self, message="Choose a file",
            defaultFile="",
            wildcard = wildcardDataFile,
            style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
            )
            


        self.rmsTotal = []
        self.rmsErrTotal = []
        self.fxy2_rmsTotal = []
        self.fxy2_rmsErrTotal = []
        self.Section_gxy2_rmsTotal = []
        self.Section_gxy2_rmsErrTotal = []
        self.AfewPairs_gxy2_rmsTotal = []
        self.AfewPairs_gxy2_rmsErrTotal = []
        
        excitonJumpDistance = []
        
        dx_1st_last = []
        dy_1st_last = []
        dr_1st_last = []
        
        dx_1st_lastN2 = []
        dy_1st_lastN2 = []
        
        
        dx_1st_last_dedrifted = []
        dy_1st_last_dedrifted = []
        dr_1st_last_dedrifted = []
        dx_1st_lastN2_dedrifted = []
        dy_1st_lastN2_dedrifted = []
        self.rmsTotal_dedrifted = []
        self.nearest_distances = []
                
        
        
        if self.cb_apply_dedrift.GetValue() == True:
            self.driftData = self.driftData_Loaded
        else:
            self.driftData = 'na'
            
            
            
        
        
        if dlg2.ShowModal() == wx.ID_OK:
            paths = dlg2.GetPaths()
            filesNames = dlg2.GetFilenames()
            filedirectoryOnly = dlg2.GetDirectory()
            print "The following file(s) are chosen:"
            
            self.filename_short = filesNames[0][3:22] 
 
        
            for n in range(len(paths)):
                print '#', n, ' : ', paths[n]
                
            for n in range(len(paths)):
                print '\n processing #', n, ' : ', paths[n]
                
                (self.MoleculeXpos, self.MoleculeYpos, self.FIONA_FrameN, self.FIONA_x0fit, self.FIONA_y0fit, self.FIONA_Afit, self.FIONA_z0fit, self.FIONA_sigma_x, self.FIONA_sigma_y, self.FIONA_theta, self.FIONA_eccentricity
                , self.FIONA_BG_IntensityAve, self.FIONA_centerIntensityAve, self.FIONA_centerIntensityMax, self.FIONA_SBratioAve, self.FIONA_SBratioMax
                , self.FIONA_AfitErr, self.FIONA_x0fitImageJErr, self.FIONA_y0fitImageJErr, self.FIONA_z0fitErr, self.FIONA_sigma_xErr, self.FIONA_sigma_yErr, self.FIONA_thetaErr
                , self.AFgshrimp_f1n, self.AFgshrimp_f2n, self.AFgshrimp_centerIntensityMaxf1n, self.AFgshrimp_centerIntensityMaxf2n, self.AFgshrimp_Afit
                , self.AFgshrimp_x0fitImageJ, self.AFgshrimp_y0fitImageJ, self.AFgshrimp_z0fit, self.AFgshrimp_sigma_x, self.AFgshrimp_sigma_y, self.AFgshrimp_theta, self.AFgshrimp_eccentricity
                , self.AFgshrimp_AfitErr, self.AFgshrimp_x0fitImageJErr, self.AFgshrimp_y0fitImageJErr, self.AFgshrimp_z0fitErr, self.AFgshrimp_sigma_xErr, self.AFgshrimp_sigma_yErr, self.AFgshrimp_thetaErr
                , self.AFshrimp_filepath, self.TotalIntensity_Trajectory, self.AdjacentFramesN) = DataProcessing_FIONA_AFgshrimp(paths[n], int(self.NLinesToSkip.GetValue()), self.driftData )
                
            
            
            
         


                fdata = open(paths[n], 'r')
                linestemp = fdata.readlines()
                fdata.close()
                #print 'linestemp \n', linestemp
                
                print ' linestemp[0] ', linestemp[0]
                
                sptemp = linestemp[0].split()
                
                print 'sptemp[3] = ', sptemp[3]
                MoleculeXpos = int(float(sptemp[3]))
                MoleculeYpos = int(float(sptemp[4]))
                
                try:
                    AdjacentFramesN = int(float(sptemp[6]))
                except:
                    AdjacentFramesN = 0
                    print 'except: AdjacentFramesN = ', AdjacentFramesN
                    
                
                
                print 'MoleculeXpos ,  MoleculeYpos', MoleculeXpos, MoleculeYpos
                
                
                #DataLines = linestemp[NDatalinesToSkip:]
                
                #print 'DataLines \n', DataLines
                
                AnalysisStartFrameN = int(float(sptemp[8]))                
                AnalysisEndFrameN = int(float(sptemp[9]))
                
                
                
                self.SectionStartFN_temp = []
                self.SectionEndFN_temp = []
                try:
                    for n2 in range(999):
                        self.SectionStartFN_temp.append(int(float(sptemp[10 + 3*n2+1])))
                        self.SectionEndFN_temp.append(int(float(sptemp[10 + 3*n2+2])))
                        
                except:
                    try:
                        if sptemp[10] == 'sections:':
                            pass
                    except:
                            self.SectionStartFN_temp = [10, 20, 30]
                            self.SectionEndFN_temp = [19, 29, 39]
                            print '\n\n ########################## \n No section information \n'
                        
                        
                
                print 'self.SectionStartFN_temp ', self.SectionStartFN_temp
                print 'self.SectionEndFN_temp ', self.SectionEndFN_temp
                
                
                
                self.AdjacentShrimpStartFrameN.SetValue(str(AnalysisStartFrameN))
                self.AdjacentShrimpEndFrameN.SetValue(str(AnalysisEndFrameN))
                
            
                self.MoleculeXYtext = '(' + str(self.MoleculeXpos)+', '+str(self.MoleculeYpos) + ')' + '\n# of Averaged Frames: ' + str(self.AdjacentFramesN*2+1)
                
                
                self.filenameOnly = filesNames[n]
                self.filedirectoryOnly = dlg2.GetDirectory()
     



           
                
                ###########################################        
                ## only for exciton jump calculation ###### Sections must be defined in each file correctly, otherwise wrong result file is made.
                if self.radio_options_autoMultiPleAnalysis.GetSelection() == 1:
                    excitonJumpDistance.append(self.AdjacentShrimpPlot(event) ) # for  exciton jump calculation
                ###########################################        
 
 
 
 
 
                ###########################################        
                ## only for distance change from the 1st and last pair
                elif self.radio_options_autoMultiPleAnalysis.GetSelection() == 2:
                    if self.cb_apply_dedrift.GetValue() == False:                            
                        dx_1st_last_temp, dy_1st_last_temp, rms_temp, dx_1st_lastN2_temp, dy_1st_lastN2_temp = self.AdjacentShrimpPlot(event)
                        dx_1st_last.append(dx_1st_last_temp)
                        dy_1st_last.append(dy_1st_last_temp)
                        dr_1st_last.append(dx_1st_last_temp**2 + dy_1st_last_temp**2)
                        dx_1st_lastN2.append(dx_1st_lastN2_temp)
                        dy_1st_lastN2.append(dy_1st_lastN2_temp)
                        self.rmsTotal.append(rms_temp)
                        
                    elif self.cb_apply_dedrift.GetValue() == True:



                        dx_1st_last_temp, dy_1st_last_temp, rms_temp, dx_1st_lastN2_temp, dy_1st_lastN2_temp = self.AdjacentShrimpPlot(event)
                        dx_1st_last_dedrifted.append(dx_1st_last_temp)
                        dy_1st_last_dedrifted.append(dy_1st_last_temp)
                        dr_1st_last_dedrifted.append(dx_1st_last_temp**2 + dy_1st_last_temp**2)
                        dx_1st_lastN2_dedrifted.append(dx_1st_lastN2_temp)
                        dy_1st_lastN2_dedrifted.append(dy_1st_lastN2_temp)
                        self.rmsTotal_dedrifted.append(rms_temp)
   



                        self.driftData = 'na'

                        (self.MoleculeXpos, self.MoleculeYpos, self.FIONA_FrameN, self.FIONA_x0fit, self.FIONA_y0fit, self.FIONA_Afit, self.FIONA_z0fit, self.FIONA_sigma_x, self.FIONA_sigma_y, self.FIONA_theta, self.FIONA_eccentricity
                        , self.FIONA_BG_IntensityAve, self.FIONA_centerIntensityAve, self.FIONA_centerIntensityMax, self.FIONA_SBratioAve, self.FIONA_SBratioMax
                        , self.FIONA_AfitErr, self.FIONA_x0fitImageJErr, self.FIONA_y0fitImageJErr, self.FIONA_z0fitErr, self.FIONA_sigma_xErr, self.FIONA_sigma_yErr, self.FIONA_thetaErr
                        , self.AFgshrimp_f1n, self.AFgshrimp_f2n, self.AFgshrimp_centerIntensityMaxf1n, self.AFgshrimp_centerIntensityMaxf2n, self.AFgshrimp_Afit
                        , self.AFgshrimp_x0fitImageJ, self.AFgshrimp_y0fitImageJ, self.AFgshrimp_z0fit, self.AFgshrimp_sigma_x, self.AFgshrimp_sigma_y, self.AFgshrimp_theta, self.AFgshrimp_eccentricity
                        , self.AFgshrimp_AfitErr, self.AFgshrimp_x0fitImageJErr, self.AFgshrimp_y0fitImageJErr, self.AFgshrimp_z0fitErr, self.AFgshrimp_sigma_xErr, self.AFgshrimp_sigma_yErr, self.AFgshrimp_thetaErr
                        , self.AFshrimp_filepath, self.TotalIntensity_Trajectory, self.AdjacentFramesN) = DataProcessing_FIONA_AFgshrimp(paths[n], int(self.NLinesToSkip.GetValue()), self.driftData )
                        
                                             
                        dx_1st_last_temp, dy_1st_last_temp, rms_temp, dx_1st_lastN2_temp, dy_1st_lastN2_temp = self.AdjacentShrimpPlot(event)
                        dx_1st_last.append(dx_1st_last_temp)
                        dy_1st_last.append(dy_1st_last_temp)
                        dr_1st_last.append(dx_1st_last_temp**2 + dy_1st_last_temp**2)
                        dx_1st_lastN2.append(dx_1st_lastN2_temp)
                        dy_1st_lastN2.append(dy_1st_lastN2_temp)
                        self.rmsTotal.append(rms_temp)
                        
                        
                        self.driftData = self.driftData_Loaded
                                                
                    
                ###########################################        
                
                
                
                ###########################################        
                ## only for nearest distance histogram
                elif self.radio_options_autoMultiPleAnalysis.GetSelection() == 3:
                    
                    nearest_distances_Temp = self.AdjacentShrimpPlot(event)
                    self.nearest_distances.extend(nearest_distances_Temp)
                    
                ## only for nearest distance histogram    
                ###########################################        
                    
                    
                    
                    
                    
                    
                    
                    
                    
               
                else:
                    rmsTemp, rmsErrTemp, fxy2_rmsTemp, fxy2_rms_ErrTemp, rms2_adjacent_section_gxy2, rms2_Err_adjacent_section_gxy2, rms2_AfewPairs, rms2_Err_AfewPairs = self.AdjacentShrimpPlot(event)
                    
                    self.rmsTotal.append(rmsTemp)
                    self.rmsErrTotal.append(rmsErrTemp)
                    
                    self.fxy2_rmsTotal.append(fxy2_rmsTemp)
                    self.fxy2_rmsErrTotal.append(fxy2_rms_ErrTemp)
                    
                    self.Section_gxy2_rmsTotal.append(rms2_adjacent_section_gxy2)
                    self.Section_gxy2_rmsErrTotal.append(rms2_Err_adjacent_section_gxy2)
                    
                    self.AfewPairs_gxy2_rmsTotal.append(rms2_AfewPairs)
                    self.AfewPairs_gxy2_rmsErrTotal.append(rms2_Err_AfewPairs)
                
        
        
        ###########################################        
        ## only for exciton jump calculation ######
        
        if self.radio_options_autoMultiPleAnalysis.GetSelection() == 1:
            todayDate = time.strftime("%Y%m%d_%Hh%Mm")
            ff = open(todayDate + '_Exciton_Jump_distance_pixel_Data_each.txt','w')
            
            for n in range(len(excitonJumpDistance)):
                for k in range(len(excitonJumpDistance[n])):
                   
                    ff.write(str(excitonJumpDistance[n][k]) + ' '   )    
                ff.write('\n')
            ff.close()       
            
            
            ff = open(todayDate + '_Exciton_Jump_distance_pixel_Data_all.txt','w')
            
            for n in range(len(excitonJumpDistance)):
                for k in range(len(excitonJumpDistance[n])):
                   
                    ff.write(str(excitonJumpDistance[n][k]) + '\n'   )    
                ff.write('\n')
            ff.close()       
            
            print 'exciton jump size calculation done. Files were created ', 
            
            return
        
        ## only for exciton jump calculation ######
        ###########################################'''





        ###########################################
        ## 1st-last pair distance and dedrift plots ##
        if self.radio_options_autoMultiPleAnalysis.GetSelection() == 2:
            dedriftText = ' '
            if self.cb_apply_dedrift.GetValue():
                dedriftText = '_dedrifted'
                        
            self.fig_dedrifted = plt.figure(self.filename_short + dedriftText, figsize = (12,8))
            plt.figtext(0.5, 0.98, self.filename_short + dedriftText
                ,ha='center', color='black', weight='normal', size='small')

            bins = np.arange(-100, 100, 5)
            xticks = np.arange(-100, 100, 20)
            plt.subplot(331)
            #plt.title()
            plt.hist([dx_1st_last, dx_1st_lastN2], label=['dx', 'dx_1st-lastN2'], bins = bins, color = ['blue', 'y'])
            #plt.hist(dx_1st_last, label='dx', bins = bins, color = 'blue')
            #plt.hist(dx_1st_lastN2, label='dx_1st-lastN2', bins = bins, color = 'y', alpha=0.6)
            plt.xlim(-100, 100)
            plt.xticks(xticks)
            plt.xlabel('1st-last distance (nm)')
            plt.legend(fontsize=11)
            
            plt.subplot(334)
            plt.hist([dy_1st_last, dy_1st_lastN2], bins = bins, alpha=1, color = ['r', 'y'],label=['dy', 'dy_1st-lastN2'])
            plt.xlim(-100, 100)
            plt.xlabel('1st-last distance (nm)')
            plt.legend(fontsize=11)
            
            plt.subplot(337)
            plt.hist(self.rmsTotal, bins = bins, color = 'g', label = 'rms')
            plt.xlim(-100, 100)
            plt.xlabel('rms distance (nm)')
            plt.legend(loc='upper left', fontsize=11)

            print '\n\n  self.rmsTotal', self.rmsTotal

            
            if self.cb_apply_dedrift.GetValue():
                plt.subplot(332)
                #plt.title()
                plt.hist([dx_1st_last_dedrifted, dx_1st_lastN2_dedrifted], label=['dx_dedrifted', 'dx_1st-lastN2_dedrifted'], bins = bins, color = ['blue','y'])
                
                #plt.hist(dx_1st_last_dedrifted, label='dx_dedrifted', bins = bins, color = 'blue')
                #plt.hist(dx_1st_lastN2_dedrifted, label='dx_1st-lastN2_dedrifted', bins = bins, color = 'y', alpha=0.6)
                plt.xlim(-100, 100)
                plt.xlabel('1st-last distance (nm) - dedrifted')
                plt.legend(fontsize=11)
                
                plt.subplot(335)
                plt.hist([dy_1st_last_dedrifted, dy_1st_lastN2_dedrifted], bins = bins, alpha=1, color = ['r', 'y'],label=['dy_dedrifted', 'dy_1st-lastN2_dedrifted'])
                plt.xlim(-100, 100)
                plt.xlabel('1st-last distance (nm) - dedrifted')
                plt.legend(fontsize=11)
                
                plt.subplot(338)
                plt.hist(self.rmsTotal_dedrifted, bins = bins, color = 'g', label = 'rms_dedrifted')
                plt.xlim(-100, 100)
                plt.xlabel('rms distance (nm) - dedrifted')
                plt.legend(loc='upper left', fontsize=11)
    



                plt.subplot(336)
                plt.title('drift: dx, dy (nm)')
                plt.plot(self.drift_Frame_N, np.array(self.drift_dx)*float(self.M_nm_pixel.GetValue()), '-bo', label='dx')
                plt.plot(self.drift_Frame_N, np.array(self.drift_dy)*float(self.M_nm_pixel.GetValue()), '-ro', label='dy')
                plt.xlabel('Frame #')
                plt.legend(loc='upper left', fontsize=11)
                plt.ylim(-50,50)
            
                
                print '\n\n  self.rmsTotal_dedrifted', self.rmsTotal_dedrifted
                
            
            plt.tight_layout()
                        

            self.fig_RMS_change_by_dedrift = plt.figure('RMS_change')
            RMS_change = np.array(self.rmsTotal_dedrifted) - np.array(self.rmsTotal)
            plt.hist(RMS_change)
            
            
            if self.cb_autoMultiPleAnalysisSaveData.GetValue():
                todayDate = time.strftime("%Y%m%d_%Hh%Mm")
                self.fig_dedrifted.savefig(todayDate + '_Before_After_dedrift.png', dpi = 80)
                self.fig_RMS_change_by_dedrift.savefig(todayDate + '_RMS_change_histogram_by_dedrift.png', dpi = 80)
                ff = open(todayDate + '_RMS_before_after_dedrift.txt','w')
                for n in range(len(self.rmsTotal_dedrifted)):
                    ff.write(str(self.rmsTotal[n]) + ' ' + str(self.rmsTotal_dedrifted[n]) + '\n'   )    
                ff.close()       



                        
            return
            
        ## 1st-last pair distance and dedrift plots ##               
        ###########################################
         
                
                
        ###########################################        
        ## only for nearest distance histogram
        elif self.radio_options_autoMultiPleAnalysis.GetSelection() == 3:
            
            totalN_molecules = len(paths)
            
            fig_nearestDistanceHistogram = plt.figure('nearest_distances')
            plt.title(filedirectoryOnly + '\nNearest distance between emission sites' + '\nTotal # of molecules: ' + str(totalN_molecules) + '     total # of distances: ' + str(len(self.nearest_distances)), size = 9)
            plt.hist(self.nearest_distances, bins = np.arange(0, 100, 5))
            plt.xlabel('distance (nm)')
            plt.show()
            
            
            
            todayDate = time.strftime("%Y%m%d_%Hh%Mm")
            
            fig_nearestDistanceHistogram.savefig(todayDate + '_nearest_distance_Histogram.png', dpi = 80)
            
            
            if self.cb_autoMultiPleAnalysisSaveData.GetValue():
                ff = open(todayDate + '_nearest_distance.txt','w')
                for n in range(len(self.nearest_distances)):
                    ff.write(str(self.nearest_distances[n]) + '\n'   )    
                ff.close()       
                
                print '\n\nNearest_distance files were created ', 
                            
            
            

            
            return
            

        ## only for nearest distance histogram    
        ###########################################        
            
            
                    






        
        self.fig_rmsHistogram = plt.figure('rms_hist', figsize=(12,7))
        plt.suptitle('RMS histogram,   Total # = ' + str(len(self.rmsTotal)) + '\n' + self.DataFileAnalysisConditions2, size = 11 )
        plt.subplots_adjust(top = 0.85, bottom = 0.05, left = 0.05, right = 0.98, wspace = 0.2, hspace = 0.4)
        

        if len(self.rmsTotal) == 1:
            self.rmsTotal = [1,2]
            self.Section_gxy2_rmsTotal = [1,2]
            self.AfewPairs_gxy2_rmsTotal = [1,2]
            self.fxy2_rmsTotal = [1,2]
            
        
        print 'self.rmsTotal ', self.rmsTotal
        print 'self.Section_gxy2_rmsTotal ', self.Section_gxy2_rmsTotal
        plt.subplot(121)
        plt.title('Adjacent shrimp RMS')
        plt.hist(self.rmsTotal, label='RMS')
        plt.hist(self.AfewPairs_gxy2_rmsTotal, label='A few pairs RMS', alpha=0.4)
        
        try:
            plt.hist(self.Section_gxy2_rmsTotal, alpha=0.4, label = 'Sectioned RMS')
        except:
            pass

        plt.legend(loc=1,prop={'size':11})
        
        
        
        plt.subplot(122)
        plt.title('FIONA RMS')
        plt.hist(self.fxy2_rmsTotal, label='FIONA RMS')
        
        plt.legend(loc=1,prop={'size':11})
        
        
                
 
                           




    def SaveFigures_AutomaticMultipleDataAnalysis(self, event):
        print 'save SaveFigures_AutomaticMultipleDataAnalysis'
        
        todayDate = time.strftime("%Y%m%d_%Hh%Mm")
        
        try:        
            for n in range(len(self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot)):    
                self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot[n].savefig(todayDate + '_Adjacent_' + str(n) + '.png', dpi = 80)
                plt.close(self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot[n])
            
            self.fig_rmsHistogram.savefig(todayDate + '_RMS_Histogram.png', dpi = 80)
            plt.close(self.fig_rmsHistogram)
            
            print 'fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot figures are saved'
            
        except:
            print 'fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot figures are NOT saved'  
        
        
        print 'self.rmsTotal', self.rmsTotal
        
        print 'self.Section_gxy2_rmsTotal', self.Section_gxy2_rmsTotal
        
    
        ff = open(todayDate + '_SHRImP_RMS_RMSErr_Section_RMS_RMSerr_Data.txt','w')
        
        for n in range(len(self.rmsTotal)):
            ff.write(str(self.rmsTotal[n]) + ' ' + str(self.rmsErrTotal[n]) + ' ' + str(self.Section_gxy2_rmsTotal[n]) + ' ' + str(self.Section_gxy2_rmsErrTotal[n]) + '\n'   )    
            
        ff.close()       
    

        '''   
        ff = open(todayDate + '_SHRImP_RMS_RMSErr_Data_' + str(self.NPairs_newRMS) + '_Pairs.txt','w')
        
        for n in range(len(self.AfewPairs_gxy2_rmsTotal)):
            ff.write(str(self.AfewPairs_gxy2_rmsTotal[n]) + ' ' + str(self.AfewPairs_gxy2_rmsErrTotal[n]) + '\n'   )    
            
        ff.close()       




        ff = open(todayDate + '_FIONA_RMS_RMSErr_Data.txt','w')
        
        for n in range(len(self.rmsTotal)):
            ff.write(str(self.fxy2_rmsTotal[n]) + ' ' + str(self.fxy2_rmsErrTotal[n]) + '\n'   )    
            
        ff.close()       
        '''
    

                
        

    def FIONA_AFSHRIMP_DATA_ANALYSIS_Plot(self, event):   
        
        
        if self.cb_apply_dedrift.GetValue() == True:
            print '\nFIONA_AFSHRIMP_DATA_ANALYSIS_Plot with dedrift'
            
            self.driftData = self.driftData_Loaded
            

            (self.MoleculeXpos, self.MoleculeYpos, self.FIONA_FrameN, self.FIONA_x0fit, self.FIONA_y0fit, self.FIONA_Afit, self.FIONA_z0fit, self.FIONA_sigma_x, self.FIONA_sigma_y, self.FIONA_theta, self.FIONA_eccentricity
            , self.FIONA_BG_IntensityAve, self.FIONA_centerIntensityAve, self.FIONA_centerIntensityMax, self.FIONA_SBratioAve, self.FIONA_SBratioMax
            , self.FIONA_AfitErr, self.FIONA_x0fitImageJErr, self.FIONA_y0fitImageJErr, self.FIONA_z0fitErr, self.FIONA_sigma_xErr, self.FIONA_sigma_yErr, self.FIONA_thetaErr
            , self.AFgshrimp_f1n, self.AFgshrimp_f2n, self.AFgshrimp_centerIntensityMaxf1n, self.AFgshrimp_centerIntensityMaxf2n, self.AFgshrimp_Afit
            , self.AFgshrimp_x0fitImageJ, self.AFgshrimp_y0fitImageJ, self.AFgshrimp_z0fit, self.AFgshrimp_sigma_x, self.AFgshrimp_sigma_y, self.AFgshrimp_theta, self.AFgshrimp_eccentricity
            , self.AFgshrimp_AfitErr, self.AFgshrimp_x0fitImageJErr, self.AFgshrimp_y0fitImageJErr, self.AFgshrimp_z0fitErr, self.AFgshrimp_sigma_xErr, self.AFgshrimp_sigma_yErr, self.AFgshrimp_thetaErr
            , self.AFshrimp_filepath, self.TotalIntensity_Trajectory, self.AdjacentFramesN) = DataProcessing_FIONA_AFgshrimp(self.dataFile_path, int(self.NLinesToSkip.GetValue()), self.driftData )
                
                

        else:
            print '\nFIONA_AFSHRIMP_DATA_ANALYSIS_Plot without dedrift'
            
            self.driftData = 'na'
            

            (self.MoleculeXpos, self.MoleculeYpos, self.FIONA_FrameN, self.FIONA_x0fit, self.FIONA_y0fit, self.FIONA_Afit, self.FIONA_z0fit, self.FIONA_sigma_x, self.FIONA_sigma_y, self.FIONA_theta, self.FIONA_eccentricity
            , self.FIONA_BG_IntensityAve, self.FIONA_centerIntensityAve, self.FIONA_centerIntensityMax, self.FIONA_SBratioAve, self.FIONA_SBratioMax
            , self.FIONA_AfitErr, self.FIONA_x0fitImageJErr, self.FIONA_y0fitImageJErr, self.FIONA_z0fitErr, self.FIONA_sigma_xErr, self.FIONA_sigma_yErr, self.FIONA_thetaErr
            , self.AFgshrimp_f1n, self.AFgshrimp_f2n, self.AFgshrimp_centerIntensityMaxf1n, self.AFgshrimp_centerIntensityMaxf2n, self.AFgshrimp_Afit
            , self.AFgshrimp_x0fitImageJ, self.AFgshrimp_y0fitImageJ, self.AFgshrimp_z0fit, self.AFgshrimp_sigma_x, self.AFgshrimp_sigma_y, self.AFgshrimp_theta, self.AFgshrimp_eccentricity
            , self.AFgshrimp_AfitErr, self.AFgshrimp_x0fitImageJErr, self.AFgshrimp_y0fitImageJErr, self.AFgshrimp_z0fitErr, self.AFgshrimp_sigma_xErr, self.AFgshrimp_sigma_yErr, self.AFgshrimp_thetaErr
            , self.AFshrimp_filepath, self.TotalIntensity_Trajectory, self.AdjacentFramesN) = DataProcessing_FIONA_AFgshrimp(self.dataFile_path, int(self.NLinesToSkip.GetValue()), self.driftData )
                
                                
                
                
                
                
                
                
                
                
        
        AFSLastFrameN = np.amax([np.amax(self.AFgshrimp_f1n),np.amax(self.AFgshrimp_f2n)])
        #print 'LastFramN: ' , AFSLastFrameN
        
        
        FIONA_LastFrameN = int(np.amax(self.FIONA_FrameN))
        #print 'FIONA_LastFrameN = ', FIONA_LastFrameN
        FrameIntensity = np.zeros(FIONA_LastFrameN + 1)
        #print 'FrameIntensity = ', FrameIntensity
        FrameIndex = np.arange(FIONA_LastFrameN + 1)
        #print 'FrameIndex ', FrameIndex
        
        for n in range(len(self.FIONA_FrameN)):
            print 'n = ', n
            FNtemp = int(self.FIONA_FrameN[n])
            #print 'FNtemp =  ', FNtemp
            FrameIntensity[FNtemp] = self.FIONA_centerIntensityMax[n]
            
        
        FrameCorMatrixEcc = np.zeros(shape=(AFSLastFrameN+1, AFSLastFrameN+1))
        
        for n in range(len(self.AFgshrimp_f1n)):
            FrameCorMatrixEcc[int(self.AFgshrimp_f1n[n]),int(self.AFgshrimp_f2n[n])] = self.AFgshrimp_eccentricity[n]
            FrameCorMatrixEcc[int(self.AFgshrimp_f2n[n]),int(self.AFgshrimp_f1n[n])] = self.AFgshrimp_eccentricity[n]
        
        
        
        print 'FrameCorMatrixEcc \n', FrameCorMatrixEcc
        
        
                
        
        
        
        AFgshrimp_AfitErrPer = ((np.array(self.AFgshrimp_AfitErr)) / (np.array(self.AFgshrimp_Afit)))*100
        
        
        FrameCorMatrixAfitErrPer = np.zeros(shape=(AFSLastFrameN+1, AFSLastFrameN+1))
        
        for n in range(len(self.AFgshrimp_f1n)):
            FrameCorMatrixAfitErrPer[int(self.AFgshrimp_f1n[n]),int(self.AFgshrimp_f2n[n])] = AFgshrimp_AfitErrPer[n]
            FrameCorMatrixAfitErrPer[int(self.AFgshrimp_f2n[n]),int(self.AFgshrimp_f1n[n])] = AFgshrimp_AfitErrPer[n]
        
        
        
        print 'FrameCorMatrixAfitErrPer \n', FrameCorMatrixAfitErrPer
        
        
        
        


        FrameCorMatrixXfitErr = np.zeros(shape=(AFSLastFrameN+1, AFSLastFrameN+1))
        
        for n in range(len(self.AFgshrimp_f1n)):
            FrameCorMatrixXfitErr[int(self.AFgshrimp_f1n[n]),int(self.AFgshrimp_f2n[n])] = self.AFgshrimp_x0fitImageJErr[n]
            FrameCorMatrixXfitErr[int(self.AFgshrimp_f2n[n]),int(self.AFgshrimp_f1n[n])] = self.AFgshrimp_x0fitImageJErr[n]

        print 'FrameCorMatrixXfitErr \n', FrameCorMatrixXfitErr
        
        
        


        FrameCorMatrixYfitErr = np.zeros(shape=(AFSLastFrameN+1, AFSLastFrameN+1))
        
        for n in range(len(self.AFgshrimp_f1n)):
            FrameCorMatrixYfitErr[int(self.AFgshrimp_f1n[n]),int(self.AFgshrimp_f2n[n])] = self.AFgshrimp_y0fitImageJErr[n]
            FrameCorMatrixYfitErr[int(self.AFgshrimp_f2n[n]),int(self.AFgshrimp_f1n[n])] = self.AFgshrimp_y0fitImageJErr[n]                
        
        
        
        
        
        
        


        FrameCorMatrixYfitSigma_x = np.zeros(shape=(AFSLastFrameN+1, AFSLastFrameN+1))
        
        for n in range(len(self.AFgshrimp_f1n)):
            FrameCorMatrixYfitSigma_x[int(self.AFgshrimp_f1n[n]),int(self.AFgshrimp_f2n[n])] = self.AFgshrimp_sigma_x[n]
            FrameCorMatrixYfitSigma_x[int(self.AFgshrimp_f2n[n]),int(self.AFgshrimp_f1n[n])] = self.AFgshrimp_sigma_x[n]                
        
        
 
        FrameCorMatrixYfitSigma_y = np.zeros(shape=(AFSLastFrameN+1, AFSLastFrameN+1))
        
        for n in range(len(self.AFgshrimp_f1n)):
            FrameCorMatrixYfitSigma_y[int(self.AFgshrimp_f1n[n]),int(self.AFgshrimp_f2n[n])] = self.AFgshrimp_sigma_y[n]
            FrameCorMatrixYfitSigma_y[int(self.AFgshrimp_f2n[n]),int(self.AFgshrimp_f1n[n])] = self.AFgshrimp_sigma_y[n]                
        
        
        
               
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        #print '\n\nFIONA_AFSHRIMP_DATA_ANALYSIS_Plot'
        ### startTime = time.strftime("%Y-%m%d, %Hh %Mm")
        
        EccentricityCondition = float(self.MaxEccentricity.GetValue())
        #print 'EccentricityCondition: ',EccentricityCondition
        
        
            
            
        self.DataFileAnalysisConditions = ('S/B: ' + str(self.SvsB_Tolerance.GetValue()) + ',   Max Ecc: ' + str(self.MaxEccentricity.GetValue()) +  
        ',   Max A err% : ' + str(self.MaxPercentage_Tolerance.GetValue()) +  ',   Max Px x,y err: ' + str(self.MaxPixel_Tolerance.GetValue()) )
        

    
        maxPer = float(self.MaxPercentage_Tolerance.GetValue())
        maxPx = float(self.MaxPixel_Tolerance.GetValue())


    
    
        fx0, fy0 = np.array(self.FIONA_x0fit), np.array(self.FIONA_y0fit)
        fxy0FrameN = np.array(self.FIONA_FrameN)
        fx0Err, fy0Err = np.array(self.FIONA_x0fitImageJErr), np.array(self.FIONA_y0fitImageJErr)
            

        
        fx1 = []
        fy1 = []
        fx2 = []
        fy2 = []
        fx2Err = []
        fy2Err = []
        fxy2FrameN = []           
        
        self.DataFileAnalysis_FIONA_fxy2_MaxIntensity = []
        
        self.DataFileAnalysis_FIONA_fxy2_FN = []
       
        
        for n in range(len(fx0)):
            FionaAfitErrPtemp = np.around(100.0 * self.FIONA_AfitErr[n] / self.FIONA_Afit[n], 0)     

            if self.FIONA_eccentricity[n] <= EccentricityCondition:
                fx1.append(fx0[n])
                fy1.append(fy0[n])
                
                if (FionaAfitErrPtemp <= maxPer) & (self.FIONA_x0fitImageJErr[n] <= maxPx) & (self.FIONA_y0fitImageJErr[n] <= maxPx):
                    fx2.append(fx0[n])
                    fy2.append(fy0[n])
                    fx2Err.append(fx0Err[n])
                    fy2Err.append(fy0Err[n])
                    fxy2FrameN.append(fxy0FrameN[n])
                                        
                    self.DataFileAnalysis_FIONA_fxy2_MaxIntensity.append(self.FIONA_centerIntensityMax[n])
                    self.DataFileAnalysis_FIONA_fxy2_FN.append(int(self.FIONA_FrameN[n]))
                    
         
        self.DataFileAnalysis_fx2 = fx2
        self.DataFileAnalysis_fy2 = fy2
        self.DataFileAnalysis_fx2Err = fx2Err
        self.DataFileAnalysis_fy2Err = fy2Err
        
        
        print 'self.DataFileAnalysis_FIONA_fxy2_FN, fx2, fy2 ', self.DataFileAnalysis_FIONA_fxy2_FN, fx2, fy2
        
        print 'self.DataFileAnalysis_fx2Err = fx2Err, self.DataFileAnalysis_fy2Err = fy2Err ', self.DataFileAnalysis_fx2Err, '\n',  self.DataFileAnalysis_fy2Err
        
        
        
        
        
        gx0, gy0 = self.AFgshrimp_x0fitImageJ, self.AFgshrimp_y0fitImageJ
        gx0Err, gy0Err = self.AFgshrimp_x0fitImageJErr, self.AFgshrimp_y0fitImageJErr
                
        
        gx1 = []
        gy1 = []
        gx2 = []
        gy2 = []
        gx2Err = []
        gy2Err = []
                
        
        self.DataFileAnalysis_AFshrimp_gxy2_F1N = []
        self.DataFileAnalysis_AFshrimp_gxy2_F2N = []
        
        
        for n in range(len(gx0)):
            AFshrimpAfitErrPtemp = np.around(100.0 * self.AFgshrimp_AfitErr[n] / self.AFgshrimp_Afit[n], 0)     
            
            if self.AFgshrimp_eccentricity[n] <= EccentricityCondition:
                gx1.append(gx0[n])
                gy1.append(gy0[n])
                
                if (AFshrimpAfitErrPtemp <= maxPer) & (self.AFgshrimp_x0fitImageJErr[n] <= maxPx) & (self.AFgshrimp_y0fitImageJErr[n] <= maxPx):
                    gx2.append(gx0[n])
                    gy2.append(gy0[n])
                    gx2Err.append(gx0Err[n])
                    gy2Err.append(gy0Err[n])
                    
                    self.DataFileAnalysis_AFshrimp_gxy2_F1N.append(int(self.AFgshrimp_f1n[n]))
                    self.DataFileAnalysis_AFshrimp_gxy2_F2N.append(int(self.AFgshrimp_f2n[n]))
        
        
        self.DataFileAnalysis_gx2 = gx2
        self.DataFileAnalysis_gy2 = gy2
        self.DataFileAnalysis_gx2Err = gx2Err
        self.DataFileAnalysis_gy2Err = gy2Err
        





        
        self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_Plot.append(plt.figure('FIONA_AFSHRIMP_DATA_ANALYSIS_Plot_' + str(self.figCounter), figsize = (18, 9)))
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(150,30,1500, 750) 
        
        print '\nself.figCounter = ', self.figCounter
        self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_Plot[self.figCounter].text(0.05, 0.93, '         x,y: '+self.MoleculeXYtext, ha="left", va="bottom", size=13, weight = 'normal', color="red")
        
        x3 = np.linspace(0, 9, 10)
        y3 = np.linspace(0, 9, 10)
        x3, y3 = np.meshgrid(x3, y3)
        
        
        
        
        
        gs = gridspec.GridSpec(12, 25)
        gs.update(bottom = 0.05, top=0.95, left=0.05, right=0.95, wspace=5.0, hspace=8.0)
        plt.suptitle(str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly) + '\n' + self.DataFileAnalysisConditions
            , size=11)
 
       


        
        axTotal_Intensity = plt.subplot(gs[1:3,0:9])
        plt.tick_params(labelsize=12)
        axTotal_Intensity.plot(self.TotalIntensity_Trajectory)
        #ax1.plot(self.FIONA_centerIntensityAve)
        axTotal_Intensity.set_title('Total Intensity Trajectory', size=14)
        plt.xlabel('Frame #')
        

        axAFSIntensity_0 = plt.subplot(gs[3:5,0:9])
        plt.tick_params(labelsize=12)
        axAFSIntensity_0.plot(FrameIndex, FrameIntensity)
        #ax1.plot(self.FIONA_centerIntensityAve)
        axAFSIntensity_0.set_title('Intensity Trajectory for AF Shrimp', size=12)
        plt.xlabel('Frame #')
      
      
    
        
 
        
        axFIONA_ECC = plt.subplot(gs[5:7,0:9])
        plt.tick_params(labelsize=12)
        axFIONA_ECC.plot(self.FIONA_FrameN, self.FIONA_eccentricity, 'o-')
        #ax1.plot(self.FIONA_centerIntensityAve)
        axFIONA_ECC.set_title('FIONA ECC Trajectory', size=12)
        plt.xlabel('Frame #')
 




        '''
        AfitErrPer = ((np.array(self.FIONA_AfitErr)) / (np.array(self.FIONA_Afit)))*100
 
        axAfitErrPer = plt.subplot(gs[7:9,0:9])
        plt.tick_params(labelsize=12)
        axAfitErrPer.set_title('FIONA Afit Err %', size=12)
        axAfitErrPer.plot(self.FIONA_FrameN, AfitErrPer, label='Afit Err (%)')
        plt.xlabel('Frame #')
        plt.legend(loc='upper left', prop={'size':10})
        plt.yticks(np.arange(0.0, 20, 5 ))
        plt.ylim(0.0, 20)
        '''



        axFitWidth = plt.subplot(gs[7:9,0:9])
        plt.tick_params(labelsize=12)
        axFitWidth.set_title('FIONA fit width', size=12)
        axFitWidth.plot(self.FIONA_FrameN, self.FIONA_sigma_x, label='x_width')
        axFitWidth.plot(self.FIONA_FrameN, self.FIONA_sigma_y, label='y_width')
        plt.xlabel('Frame #')
        plt.legend(loc='upper left', prop={'size':10})
        plt.yticks(np.arange(0.0, 4, .5 ))
        plt.ylim(0.0, 4)





        axFXY0err = plt.subplot(gs[9:11,0:9])
        plt.tick_params(labelsize=12)
        axFXY0err.set_title('FIONA x, y fit err', size=12)
        axFXY0err.plot(self.FIONA_FrameN, self.FIONA_x0fitImageJErr, label='x err')
        axFXY0err.plot(self.FIONA_FrameN, self.FIONA_y0fitImageJErr, label='y err')
        #plt.xlim( 0, np.max(fxy0FrameN) )
        plt.xlabel('Frame #')
        plt.legend(loc='upper left', prop={'size':10})
        plt.yticks(np.arange(0.0, 0.2, 0.05 ))
        plt.ylim(0.0, 0.2)
        
 
 
        
        
        #print '\nfx2 ', fx2 
        #print 'len(fx2) ', len(fx2)
        #exit()
 
        xposAve = int(np.mean(fx2) )
        yposAve = int(np.mean(fy2) )
                
 
 
                
        
        
        axFIOINA2 = plt.subplot(gs[2:6,9:13], aspect='equal')
        plt.tick_params(labelsize=10)
        axFIOINA2.set_title('FIONA2', size =14)
        axFIOINA2.scatter(fx2, fy2, marker = '.',edgecolors='None', c = range(len(fx2)),  s=35, vmin=0, vmax= len(fx2))
        plt.xlim( xposAve - 2, xposAve + 2 )
        plt.ylim( yposAve - 2, yposAve + 2 )
        plt.xticks(np.arange(xposAve - 2, xposAve + 3, 1 ))
        plt.xlabel('Pixel Position')
        
        
        
        
        

        axAFShrimp2 = plt.subplot(gs[7:11,9:13], aspect='equal')
        plt.tick_params(labelsize=10)
        axAFShrimp2.set_title('AF Shrimp2', size =14)
        axAFShrimp2.scatter(gx2, gy2, marker = '.', edgecolors='None', c = range(len(gx2)),  s=35, vmin=0, vmax= len(gx2))
        plt.xlim( (np.around(xposAve) - 2, np.around(xposAve) + 2 ) )
        plt.ylim( (np.around(yposAve) - 2, np.around(yposAve) + 2 ) )    
        plt.xticks(np.arange(xposAve - 2, xposAve + 3, 1 ))
        plt.xlabel('Pixel Position')
        
        





        gs2 = gridspec.GridSpec(12, 25)
        gs2.update(bottom = 0.05, top=0.95, left=0.085 , right=0.922 , wspace=5.0, hspace=2.0)
        
        axAFSIntensity = plt.subplot(gs2[2:4,13:20])
        plt.tick_params(labelsize=12)
        axAFSIntensity.plot(FrameIndex, FrameIntensity)
        #ax1.plot(self.FIONA_centerIntensityAve)
        axAFSIntensity.set_title('Intensity Trajectory for AF Shrimp', size=15)
        plt.xlabel('Frame #')
        plt.xlim(0, AFSLastFrameN+1)
      
      
      
        
        
        axFrameCorEcc = plt.subplot(gs[4:12,13:21], aspect='equal')
        plt.tick_params(labelsize=10)
        axFrameCorEcc.set_title("Frames' Eccentricity Correlation", size =16)
        plt.imshow(FrameCorMatrixEcc, interpolation = 'None', origin='bottom', vmin=0.2, vmax=1.2)
        plt.colorbar()     
        plt.xlabel('Frame #')
        plt.ylabel('Frame #')
        








        axFrameCorAfitErrPer = plt.subplot(gs[0:4,21:25], aspect='equal')
        plt.tick_params(labelsize=10)
        axFrameCorAfitErrPer.set_title("Frames' AfitErrPer Corr", size =10)
        plt.imshow(FrameCorMatrixAfitErrPer, interpolation = 'None', origin='bottom', vmin=0, vmax=40)
        plt.colorbar()     
        plt.xlabel('Frame #')
        plt.ylabel('Frame #')




        axFrameCorXfitErr = plt.subplot(gs[4:8,21:25], aspect='equal')
        plt.tick_params(labelsize=10)
        axFrameCorXfitErr.set_title("Frames' XfitErr Corr", size = 10)
        plt.imshow(FrameCorMatrixXfitErr, interpolation = 'None', origin='bottom', vmin=0, vmax=0.3)
        plt.colorbar()     
        plt.xlabel('Frame #')
        plt.ylabel('Frame #')







        axFrameCorYfitErr = plt.subplot(gs[8:12,21:25], aspect='equal')
        plt.tick_params(labelsize=10)
        axFrameCorYfitErr.set_title("Frames' YfitErr Corr", size = 10)
        plt.imshow(FrameCorMatrixYfitErr, interpolation = 'None', origin='bottom', vmin=0, vmax=0.3)
        plt.colorbar()     
        plt.xlabel('Frame #')
        plt.ylabel('Frame #')






        
        self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_ACF_Plot = plt.figure('FIONA_AFSHRIMP_DATA_ANALYSIS_ACF_Plot_supplement', figsize = (18, 9))
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(300,40,1280, 750) 
        
        
        self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_ACF_Plot.text(0.02, 0.91, '         x,y: '+self.MoleculeXYtext, ha="left", va="bottom", size=11, weight = 'normal', color="red")
        
        x3 = np.linspace(0, 9, 10)
        y3 = np.linspace(0, 9, 10)
        x3, y3 = np.meshgrid(x3, y3)
        
        
        
        
        
        gs = gridspec.GridSpec(11, 25)
        gs.update(bottom = 0.05, top=0.95, left=0.05, right=0.95, wspace=5.0, hspace=5.0)
        plt.suptitle(str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly) + '\n' + self.DataFileAnalysisConditions
            , size=11)
 
       

        ######################################
        # Calculating AutoCorrelationFunction

        if len(self.TotalIntensity_Trajectory) >= 2:
            TotalIntensityACF = AutoCorrFunc(self.TotalIntensity_Trajectory)
        else:
            TotalIntensityACF = [0,0]
            print 'TotalIntensityACF ', TotalIntensityACF
            
        
        #fx0ACF = AutoCorrFunc(fx0)
        #fy0ACF = AutoCorrFunc(fy0)
        
        FIONAdx0 = fx0-np.mean(fx0)
        FIONAdy0 = fy0-np.mean(fy0)
        FIONAdx2 = fx2-np.mean(fx2)
        FIONAdy2 = fy2-np.mean(fy2)
        
        Fdx0ACF = AutoCorrFunc(FIONAdx0)
        Fdy0ACF = AutoCorrFunc(FIONAdy0)
        Fdx2ACF = AutoCorrFunc(FIONAdx2)
        Fdy2ACF = AutoCorrFunc(FIONAdy2)
        
 

        # Calculating AutoCorrelationFunction#
        ######################################





    
        ### self.TotalIntensity_Trajectory        
        axTotal_Intensity = plt.subplot(gs[1:3,0:12])
        plt.tick_params(labelsize=12)
        axTotal_Intensity.plot(self.TotalIntensity_Trajectory)
        #ax1.plot(self.FIONA_centerIntensityAve)
        axTotal_Intensity.set_title('Total Intensity Trajectory', size=16)
        plt.xlabel('Frame #')
        plt.xlim(0, len(self.TotalIntensity_Trajectory))
 
        
        axFIONA_ECC = plt.subplot(gs[1:3,12:18])
        plt.tick_params(labelsize=10)
        axFIONA_ECC.plot(self.FIONA_FrameN, self.FIONA_eccentricity, 'o-')
        #ax1.plot(self.FIONA_centerIntensityAve)
        axFIONA_ECC.set_title('FIONA ECC Trajectory', size=12)
        plt.xlabel('Frame #')
 


        axFX0 = plt.subplot(gs[3:7,0:12])
        plt.tick_params(labelsize=10)
        axFX0.plot(fxy0FrameN, FIONAdx0, 'o-', label='Fdx0')
        axFX0.plot(fxy0FrameN, FIONAdy0, 'o-', label='Fdy0')
        axFX0.plot(fxy2FrameN, FIONAdx2, 'x-', label='Fdx2')
        axFX0.plot(fxy2FrameN, FIONAdy2, 'x-', label='Fdy2')
        #axFX0.plot(FrameIndex, FrameIntensity, 'o-')
        #ax1.plot(self.FIONA_centerIntensityAve)
        axFX0.set_title('dfx0 dfy0', size=12)
        plt.xlim(0, len(self.TotalIntensity_Trajectory))
        plt.xlabel('Frame #')
        plt.ylabel('pixel')
        plt.legend(loc='upper right', prop={'size':10})
        
 
        axACFxy = plt.subplot(gs[7:11,0:12])
        plt.tick_params(labelsize=10)
        axACFxy.plot(TotalIntensityACF, 'o-', label='Total I ACF')
        #axACFxy.plot(fx0ACF, 'o-', label='fx0 ACF')
        axACFxy.plot(Fdx0ACF, 'o-', label='Fdx0ACF')
        axACFxy.plot(Fdy0ACF, 'o-', label='Fdy0ACF')
        axACFxy.plot(Fdx2ACF, 'x-', label='Fdx2ACF')
        axACFxy.plot(Fdy2ACF, 'x-', label='Fdy2ACF')
        #ax1.plot(self.FIONA_centerIntensityAve)
        axACFxy.set_title('ACF', size=12)
        #plt.xlim( 0, np.max(fxy0FrameN) )
        plt.xlabel('Frame #')
        plt.xlim(0, len(self.TotalIntensity_Trajectory))
        plt.legend(loc='upper right', prop={'size':10})
        


        gs2 = gridspec.GridSpec(11, 25)
        gs2.update(bottom = 0.05, top=0.95, left=0.102, right=0.924, wspace=5.0, hspace=2.0)
        
        axAFSIntensity = plt.subplot(gs2[1:3,19:26])
        plt.tick_params(labelsize=12)
        axAFSIntensity.plot(FrameIndex, FrameIntensity)
        #ax1.plot(self.FIONA_centerIntensityAve)
        axAFSIntensity.set_title('Intensity Trajectory for AF Shrimp', size=15)
        plt.xlabel('Frame #')
      
      











        ### publication contour maps figures ###
        if self.cb_ShowPublicationFigures.GetValue():
            self.fig_contourmaps = plt.figure('countour maps', figsize = (22, 10))
            plt.subplots_adjust(top = 0.90, bottom = 0.10, left = 0.04, right = 0.98, wspace = 0.15, hspace = 0.25)
            
            plt.subplot(241, aspect='equal')
            plt.tick_params(labelsize=10)
            plt.title(r'$\delta$ A', size = 22)
            plt.imshow(FrameCorMatrixAfitErrPer, interpolation = 'None', origin='bottom', vmin=0, vmax=100)
            plt.colorbar()     
            plt.xlabel('Frame #', fontsize=16)
            plt.ylabel('Frame #', fontsize=16)
    
            
            
            

            plt.subplot(242, aspect='equal')
            plt.tick_params(labelsize=10)
            plt.title("Frames' X fit Err Corr", size = 22)
            plt.imshow(FrameCorMatrixXfitErr, interpolation = 'None', origin='bottom', vmin=0, vmax=0.3)
            plt.colorbar()     
            plt.xlabel('Frame #', fontsize=16)
            plt.ylabel('Frame #', fontsize=16)
    
    
    
            plt.subplot(243, aspect='equal')
            plt.tick_params(labelsize=10)
            plt.title("Frames' Y fit Err Corr", size = 22)
            plt.imshow(FrameCorMatrixXfitErr, interpolation = 'None', origin='bottom', vmin=0, vmax=0.3)
            plt.colorbar()     
            plt.xlabel('Frame #', fontsize=16)
            plt.ylabel('Frame #', fontsize=16)
            
            
            
            FrameCorMatrixRfitErr = (FrameCorMatrixXfitErr + FrameCorMatrixYfitErr)/2.0
    
            plt.subplot(245, aspect='equal')
            plt.tick_params(labelsize=10)
            plt.title(r'$\delta$ r', size = 22)
            #plt.imshow(FrameCorMatrixRfitErr, interpolation = 'None', origin='bottom', vmin=0, vmax=0.3)
            plt.imshow(FrameCorMatrixRfitErr*92.6, interpolation = 'None', origin='bottom', vmin=0, vmax=50)
            #plt.imshow(FrameCorMatrixRfitErr, interpolation = 'None', origin='bottom')
            plt.colorbar()     
            plt.xlabel('Frame #', fontsize=16)
            plt.ylabel('Frame #', fontsize=16)
    
    
    

            plt.subplot(246, aspect='equal')
            plt.tick_params(labelsize=10)
            plt.title(r'$\sigma$ x', size = 22)
            plt.imshow(FrameCorMatrixYfitSigma_x*92.6, interpolation = 'None', origin='bottom', vmin=0, vmax=500)
            plt.colorbar()     
            plt.xlabel('Frame #', fontsize=16)
            plt.ylabel('Frame #', fontsize=16)
    
            
    

            plt.subplot(247, aspect='equal')
            plt.tick_params(labelsize=10)
            plt.title(r'$\sigma$ y', size = 22)
            plt.imshow(FrameCorMatrixYfitSigma_y*92.6, interpolation = 'None', origin='bottom', vmin=0, vmax=500)
            plt.colorbar()     
            plt.xlabel('Frame #', fontsize=16)
            plt.ylabel('Frame #', fontsize=16)
    
            

            
            self.fig_contourmaps.savefig('tempContourMap.png', transparent=True, dpi=400)
            
            
            









        #plt.tight_layout()        
        plt.ion()
        plt.show()
        
        self.btn_EmissionSitesAnalysis.Show()
        
        
        
        self.figCounter += 1

        allX = np.append(fx2, gx2)
        allY = np.append(fy2, gy2)
        
        afshrimpDensity = TowD_PointDensityCalculation(allX, allY, 0.1, 2)
        #levels = np.arange(10.0, 20.0, 50.0)
        
        x3countour = np.linspace(0, len(afshrimpDensity)-1, len(afshrimpDensity)*5)
        y3countour = np.linspace(0, len(afshrimpDensity)-1, len(afshrimpDensity)*5)
        x3countour, y3countour = np.meshgrid(x3countour, y3countour)
    

        
        self.fig_afshrimpDensity = plt.figure('afdensity')
        plt.imshow(afshrimpDensity, interpolation = 'None', origin='bottom')
        plt.colorbar()  
        plt.contour(x3countour, y3countour, scipy.ndimage.zoom(afshrimpDensity,5), [np.max(afshrimpDensity)*0.1], colors='w')
        
        


        
        

    def EmissionSitesAnalysis(self, event):
        print '\n Starting Emission Sites Analysis'
        
        self.text_Save_figures.SetLabel(' ')
        
        
        SectionBeginFN = [] # starting frame number for each section, eg) ...[0] value is the starting frome number of the first section
        SectionEndFN = []   
        
        TotalNsection = 0
        
        print 'len(self.F_SectionBegin) = ', len(self.F_SectionBegin)
        
        for n in range(len(self.F_SectionBegin)):
            if (self.F_SectionBegin[n].GetValue() == '0') & (self.F_SectionEnd[n].GetValue() == '0' ):
                break
            else:
                TotalNsection += 1
                
        print 'TotalNsection = ',   TotalNsection 
        
        
        
        
        for n in range(TotalNsection):
            SectionBeginFN.append(int(self.F_SectionBegin[n].GetValue() ))
            SectionEndFN.append(int(self.F_SectionEnd[n].GetValue() ) )
        
        print 'SectionBeginFN ', SectionBeginFN
        print 'SectionEndFN ', SectionEndFN
        
        
        

        
        Section_Intensity = [[] for _ in xrange(TotalNsection)]
        Section_FIONA_xpos = [[] for _ in xrange(TotalNsection)]
        Section_FIONA_ypos = [[] for _ in xrange(TotalNsection)]
        Section_FIONA_xposErr = [[] for _ in xrange(TotalNsection)]
        Section_FIONA_yposErr = [[] for _ in xrange(TotalNsection)]
                
        
        
        print 'self.DataFileAnalysis_FIONA_fxy2_FN = ', self.DataFileAnalysis_FIONA_fxy2_FN
        print 'self.DataFileAnalysis_fx2 = ', self.DataFileAnalysis_fx2
        print '\n\nself.DataFileAnalysis_FIONA_fxy2_MaxIntensity', self.DataFileAnalysis_FIONA_fxy2_MaxIntensity
        
        for n in range(len(self.DataFileAnalysis_FIONA_fxy2_FN)):
            FNtemp = self.DataFileAnalysis_FIONA_fxy2_FN[n]
            print 'FNtemp = ', FNtemp
            
            for k in range(TotalNsection):
                if SectionBeginFN[k] <= FNtemp <= SectionEndFN[k]:
                    Section_FIONA_xpos[k].append(self.DataFileAnalysis_fx2[n])
                    Section_FIONA_ypos[k].append(self.DataFileAnalysis_fy2[n])
                    
                    Section_FIONA_xposErr[k].append(self.DataFileAnalysis_fx2Err[n])
                    Section_FIONA_yposErr[k].append(self.DataFileAnalysis_fy2Err[n])
                    
                    Section_Intensity[k].append(self.DataFileAnalysis_FIONA_fxy2_MaxIntensity[n])
                    
                    
        

        print 'Section_FIONA_xpos = ', Section_FIONA_xpos
        print 'Section_FIONA_ypos = ', Section_FIONA_ypos
        
        print 'Section_FIONA_xposErr = ', Section_FIONA_xposErr
        print 'Section_FIONA_yposErr = ', Section_FIONA_yposErr
        
        
        print '\nSection_Intensity = \n', Section_Intensity
            
            
            
            
        Section_Intensity_ave = []
        Section_FIONA_xpos_WeightedAve = []
        Section_FIONA_ypos_WeightedAve = []
        Section_FIONA_xpos_WeightedAveErr = []
        Section_FIONA_ypos_WeightedAveErr = []
        
        for n in range(TotalNsection):
            Section_Intensity_ave.append(np.mean(Section_Intensity[n]))
            xposWAveTemp, xposWAveErrTemp = Weighted_Ave_Err(Section_FIONA_xpos[n], Section_FIONA_xposErr[n])
            yposWAveTemp, yposWAveErrTemp = Weighted_Ave_Err(Section_FIONA_ypos[n], Section_FIONA_yposErr[n])
            
            
            Section_FIONA_xpos_WeightedAve.append(xposWAveTemp)
            Section_FIONA_xpos_WeightedAveErr.append(xposWAveErrTemp)
            
            Section_FIONA_ypos_WeightedAve.append(yposWAveTemp)            
            Section_FIONA_ypos_WeightedAveErr.append(yposWAveErrTemp)
            
                

        #Weighted_Ave_Err()


        
        Section_FIONA_xypos_WeightedAveErr = (np.array(Section_FIONA_xpos_WeightedAveErr) + np.array(Section_FIONA_ypos_WeightedAveErr))/2.0
        
        
        
        
        
        
        Section_AFshrimp_xpos = [[] for _ in xrange(TotalNsection)]
        Section_AFshrimp_ypos = [[] for _ in xrange(TotalNsection)]
        Section_AFshrimp_xposErr = [[] for _ in xrange(TotalNsection)]
        Section_AFshrimp_yposErr = [[] for _ in xrange(TotalNsection)]
        
        Section_AFshrimp_CorrespondingSectionN = [[] for _ in xrange(TotalNsection)]
        
        print '\n\n##############################################\nstarting AF shrimp'
        for n in range(TotalNsection):
            
            print '\n\n ################## section n = ', n
            gx2postemp = []
            gy2postemp = []
            gx2posErrtemp = []
            gy2posErrtemp = []
            CorrSectiontemp = []
            
            for k in range(TotalNsection):
                print '\n $$$$$$$$$$$$$ section k = ', k

                kxtemp = []
                kytemp = []
                kxErrtemp = []
                kyErrtemp = []
                
                if Section_Intensity_ave[n] > Section_Intensity_ave[k]:
                    for p in range(len(self.DataFileAnalysis_AFshrimp_gxy2_F1N)):
                        
                       
                        if (SectionBeginFN[n] <= self.DataFileAnalysis_AFshrimp_gxy2_F1N[p] <= SectionEndFN[n]) & (SectionBeginFN[k] <= self.DataFileAnalysis_AFshrimp_gxy2_F2N[p] <= SectionEndFN[k]):   
                            kxtemp.append(self.DataFileAnalysis_gx2[p])
                            kytemp.append(self.DataFileAnalysis_gy2[p])
                            kxErrtemp.append(self.DataFileAnalysis_gx2Err[p])
                            kyErrtemp.append(self.DataFileAnalysis_gy2Err[p])
                            

                        elif (SectionBeginFN[n] <= self.DataFileAnalysis_AFshrimp_gxy2_F2N[p] <= SectionEndFN[n]) & (SectionBeginFN[k] <= self.DataFileAnalysis_AFshrimp_gxy2_F1N[p] <= SectionEndFN[k]):   
                            kxtemp.append(self.DataFileAnalysis_gx2[p])
                            kytemp.append(self.DataFileAnalysis_gy2[p])
                            kxErrtemp.append(self.DataFileAnalysis_gx2Err[p])
                            kyErrtemp.append(self.DataFileAnalysis_gy2Err[p])
                            

                        else:
                            pass
                    gx2postemp.append(kxtemp)
                    gy2postemp.append(kytemp)
                    gx2posErrtemp.append(kxErrtemp)
                    gy2posErrtemp.append(kyErrtemp)
                    
                    
                    #if len(gx2postemp) != 0:
                    CorrSectiontemp.append(k)
                                            
                    
                else:
                    pass
                
            Section_AFshrimp_xpos[n] += (gx2postemp)
            Section_AFshrimp_ypos[n] += (gy2postemp)
            Section_AFshrimp_xposErr[n] += (gx2posErrtemp)
            Section_AFshrimp_yposErr[n] += (gy2posErrtemp)
            Section_AFshrimp_CorrespondingSectionN[n] += (CorrSectiontemp)
            

        
        Section_AFshrimp_index = []
        Section_AFshrimp_xposWeightedAve = []
        Section_AFshrimp_yposWeightedAve = []
        Section_AFshrimp_xposWeightedAveErr = []
        Section_AFshrimp_yposWeightedAveErr = []
        Section_AFshrimp_xyposWeightedAveErr_Mean = []
        
        for n in range(TotalNsection):
            
            for k in range(len(Section_AFshrimp_xpos[n])):
                
                xposWAveTemp, xposWAveErrTemp = Weighted_Ave_Err(Section_AFshrimp_xpos[n][k], Section_AFshrimp_xposErr[n][k])
                yposWAveTemp, yposWAveErrTemp = Weighted_Ave_Err(Section_AFshrimp_ypos[n][k], Section_AFshrimp_yposErr[n][k])
                
                
                Section_AFshrimp_index.append([n, Section_AFshrimp_CorrespondingSectionN[n][k] ])
                Section_AFshrimp_xposWeightedAve.append(xposWAveTemp)
                Section_AFshrimp_yposWeightedAve.append(yposWAveTemp)
                
                Section_AFshrimp_xposWeightedAveErr.append(xposWAveErrTemp)
                Section_AFshrimp_yposWeightedAveErr.append(yposWAveErrTemp)
                
                Section_AFshrimp_xyposWeightedAveErr_Mean.append( (xposWAveErrTemp + yposWAveErrTemp )/2.0 )
           
  
  
  
  
  
  
  
        
        self.fig_EmissionSitesAnalysis.append(plt.figure('Emission_Sites_ANALYSIS_Plot_' + str(self.fig_EmissionSitesAnalysisFinalCounter), figsize = (18, 9)))
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(50,30,1600, 800)        
        plt.suptitle(str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly) + '\n' + self.DataFileAnalysisConditions
            , size=11)
        self.fig_EmissionSitesAnalysis[self.fig_EmissionSitesAnalysisFinalCounter].text(0.05, 0.94, '         x,y: '+self.MoleculeXYtext, ha="left", va="bottom", size=13, weight = 'normal', color="red")   

        plt.subplots_adjust(top = 0.90, bottom = 0.03, left = 0.03, right = 0.98, wspace = 0.2, hspace = 0.4)
        markerColor = ['b' , 'g' , 'c', 'm', 'y', 'k', 'b', 'g', 'c', 'm', 'y', 'k']
        
        for n in range(TotalNsection):
                
            plt.subplot(2,5,n+1, aspect='equal')
            plt.tick_params(labelsize=10)
            plt.title('Section '+str(n+1) + ',   F# ' +str(SectionBeginFN[n])+':'+str(SectionEndFN[n])       , size =12)
            
            for k in range(len(Section_AFshrimp_xpos[n])):
                if len(Section_AFshrimp_xpos[n][k]) != 0:
                    plt.scatter(Section_AFshrimp_xpos[n][k], Section_AFshrimp_ypos[n][k], marker = 'x', edgecolors=markerColor[k], color = markerColor[k], label = 'S'+str(n+1) +'-' + str(Section_AFshrimp_CorrespondingSectionN[n][k] +1))
                
                                
            #plt.scatter(gx2, gy2, marker = '.', edgecolors='None', c = range(len(gx2)),  s=35, vmin=0, vmax= len(gx2)) 
            print 'len(Section_FIONA_xpos[n]) = ', len(Section_FIONA_xpos[n])
            
            if len(Section_FIONA_xpos[n]) != 0 :
                plt.scatter(Section_FIONA_xpos[n], Section_FIONA_ypos[n], marker = '.',s= 100, edgecolors='r', color = 'r', label = 'FIONA')
                plt.scatter(Section_FIONA_xpos[n], Section_FIONA_ypos[n], marker = '.', edgecolors='None', c = range(len(Section_FIONA_xpos[n])),  s=35, vmin=0, vmax= len(Section_FIONA_xpos[n]), label = 'FIONA')
            
                plt.legend(loc='upper right', prop={'size':10})
                plt.xlim( (np.around(Section_FIONA_xpos_WeightedAve[n]) - 1, np.around(Section_FIONA_xpos_WeightedAve[n]) + 1 ) )
                plt.ylim( (np.around(Section_FIONA_ypos_WeightedAve[n]) - 1, np.around(Section_FIONA_ypos_WeightedAve[n]) + 1 ) )   
                plt.xlabel('Pixel Position')
                







            
        self.fig_EmissionSitesAnalysisFinal.append(plt.figure('Emission_Sites_ANALYSIS_Final_'+str(self.fig_EmissionSitesAnalysisFinalCounter), figsize = (18, 6)) )
        self.fig_EmissionSitesAnalysisFinal[self.fig_EmissionSitesAnalysisFinalCounter].text(0.05, 0.92, '         x,y: '+self.MoleculeXYtext, ha="left", va="bottom", size=13, weight = 'normal', color="red")   
        
        
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(50,500,1600, 533)
        #mngr.window.move(50, 500)
        

        
        
        plt.suptitle(str(self.filedirectoryOnly) + ' \n ' + str(self.filenameOnly) + '\n' + self.DataFileAnalysisConditions, size=11)
        plt.subplots_adjust(top = 0.84, bottom = 0.05, left = 0.04, right = 0.98, wspace = 0.2, hspace = 0.4)  
        
        markerColorAFS = ['r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k','r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k','r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k','r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k']

        plt.subplot(1,3,1, aspect='equal')
        plt.title("FIONA")
        print 'Section_FIONA_xpos_WeightedAve = ', Section_FIONA_xpos_WeightedAve

        FIONA_Cir = [[] for _ in xrange(TotalNsection)]
        for n in range(len(Section_FIONA_xpos_WeightedAve)):
            print 'n = ', n
            if str(Section_FIONA_xpos_WeightedAve[n]) != 'nan':
                print 'Section_FIONA_xpos_WeightedAve[n] ', Section_FIONA_xpos_WeightedAve[n], type(Section_FIONA_xpos_WeightedAve[n])
                print 'run n = ', n
                plt.scatter(Section_FIONA_xpos_WeightedAve[n], Section_FIONA_ypos_WeightedAve[n], marker = '.',s= 100, edgecolors = markerColorAFS[n], c = markerColorAFS[n]
                , label = 'S' + str(n+1) + ' ('+str(np.around(Section_FIONA_xpos_WeightedAve[n],3))+', ' + str(np.around(Section_FIONA_ypos_WeightedAve[n],3)) + ')  E ' + str(np.around(Section_FIONA_xypos_WeightedAveErr[n],3))  )
            else:
                print 'pass'
                pass
            FIONA_Cir[n] = plt.Circle((Section_FIONA_xpos_WeightedAve[n], Section_FIONA_ypos_WeightedAve[n]), radius = Section_FIONA_xypos_WeightedAveErr[n], color = markerColorAFS[n], fill=False)
            self.fig_EmissionSitesAnalysisFinal[self.fig_EmissionSitesAnalysisFinalCounter].gca().add_artist(FIONA_Cir[n])
            

        plt.legend(loc='upper right', prop={'size':12})
        plt.xlim( (np.around(Section_FIONA_xpos_WeightedAve[0]) - 1, np.around(Section_FIONA_xpos_WeightedAve[0]) + 1 ) )
        plt.ylim( (np.around(Section_FIONA_ypos_WeightedAve[0]) - 1, np.around(Section_FIONA_ypos_WeightedAve[0]) + 1 ) )     
            
        
        
        self.Section_FIONA_xpos_WeightedAve = Section_FIONA_xpos_WeightedAve
        self.Section_FIONA_ypos_WeightedAve = Section_FIONA_ypos_WeightedAve
        self.Section_FIONA_xypos_WeightedAveErr = Section_FIONA_xypos_WeightedAveErr
        self.TotalNsection = TotalNsection

        

        plt.subplot(1,3,2, aspect='equal')
        plt.title("AF Shrimp")

        N_SectionAfshrimpIndex = len(Section_AFshrimp_index)
        AFshrimp_Cir = [[] for _ in xrange(N_SectionAfshrimpIndex)]
        
        
        print 'N_SectionAfshrimpIndex ', N_SectionAfshrimpIndex
        print 'Section_AFshrimp_xposWeightedAve = ' , Section_AFshrimp_xposWeightedAve
        
        for n in range(N_SectionAfshrimpIndex):
            print ' n = ', n

            if str(Section_AFshrimp_xposWeightedAve[n]) != 'nan':
                plt.scatter(Section_AFshrimp_xposWeightedAve[n], Section_AFshrimp_yposWeightedAve[n], marker = 'x', c = markerColorAFS[n]
                , label = 'S' + str(Section_AFshrimp_index[n][0]+1) + '-' + str(Section_AFshrimp_index[n][1]+1) + ' ('+str(np.around(Section_AFshrimp_xposWeightedAve[n],3))+', ' + str(np.around(Section_AFshrimp_yposWeightedAve[n],3)) + ')  E ' + str(np.around(Section_AFshrimp_xyposWeightedAveErr_Mean[n],3)) )
            else:
                print 'pass'
                pass
            
            AFshrimp_Cir[n] = plt.Circle((Section_AFshrimp_xposWeightedAve[n], Section_AFshrimp_yposWeightedAve[n]), radius = Section_AFshrimp_xyposWeightedAveErr_Mean[n], color = markerColorAFS[n], fill=False)
            self.fig_EmissionSitesAnalysisFinal[self.fig_EmissionSitesAnalysisFinalCounter].gca().add_artist(AFshrimp_Cir[n])
            

        plt.legend(loc='upper right', prop={'size':12})
        plt.xlim( (np.around(Section_FIONA_xpos_WeightedAve[0]) - 1, np.around(Section_FIONA_xpos_WeightedAve[0]) + 1 ) )
        plt.ylim( (np.around(Section_FIONA_ypos_WeightedAve[0]) - 1, np.around(Section_FIONA_ypos_WeightedAve[0]) + 1 ) )     
            
        
        
        

        plt.subplot(1,3,3, aspect='equal')
        plt.title("FIONA + AF Shrimp")

        N_SectionAfshrimpIndex = len(Section_AFshrimp_index)
        AFshrimp_Cir = [[] for _ in xrange(N_SectionAfshrimpIndex)]

        FIONA_Cir = [[] for _ in xrange(TotalNsection)]
        for n in range(len(Section_FIONA_xpos_WeightedAve)):
            print 'n = ', n
            if str(Section_FIONA_xpos_WeightedAve[n]) != 'nan':
                print 'Section_FIONA_xpos_WeightedAve[n] ', Section_FIONA_xpos_WeightedAve[n], type(Section_FIONA_xpos_WeightedAve[n])
                print 'run n = ', n
                plt.scatter(Section_FIONA_xpos_WeightedAve[n], Section_FIONA_ypos_WeightedAve[n], marker = '.',s= 100, edgecolors = markerColorAFS[n], c = markerColorAFS[n], label = 'FIONA ' + str(n+1))
            else:
                print 'pass'
                pass

            
            
        for n in range(N_SectionAfshrimpIndex):

            if str(Section_AFshrimp_xposWeightedAve[n]) != 'nan':
                plt.scatter(Section_AFshrimp_xposWeightedAve[n], Section_AFshrimp_yposWeightedAve[n], marker = 'x', c = markerColorAFS[n])
            else:
                print 'pass'
                pass
            

            
        plt.legend(loc='upper right', prop={'size':12})
        plt.xlim( (np.around(Section_FIONA_xpos_WeightedAve[0]) - 1, np.around(Section_FIONA_xpos_WeightedAve[0]) + 1 ) )
        plt.ylim( (np.around(Section_FIONA_ypos_WeightedAve[0]) - 1, np.around(Section_FIONA_ypos_WeightedAve[0]) + 1 ) )     
        
        
        self.fig_EmissionSitesAnalysisFinalCounter += 1
        
        print 'self.fig_EmissionSitesAnalysisFinalCounter = ', self.fig_EmissionSitesAnalysisFinalCounter
        
        self.btn_Save_figures.Show()   
        #self.btn_EmissionSitesAnalysis.Hide()
        self.btn_PlotFinalEmissionSites.Show()
        
        
        
        
        self.Section_FIONA_xpos = Section_FIONA_xpos
        self.Section_FIONA_ypos = Section_FIONA_ypos

        #Section_AFshrimp_index = []
        self.Section_AFshrimp_xpos = Section_AFshrimp_xpos
        self.Section_AFshrimp_ypos = Section_AFshrimp_ypos
        self.Section_AFshrimp_CorrespondingSectionN = Section_AFshrimp_CorrespondingSectionN
        






        
    def Save_EmissionSitesAnalysisData(self, enven):
        print '\nSave_EmissionSitesAnalysisData'
        
        prefix = self.FilePrefix.GetValue()
        todayDate = time.strftime("%Y%m%d_%Hh%Mm")
        figFileName = prefix + self.filenameOnly[:35]
        
        
        
        
        
    
        ff = open('SectionData_' + figFileName + todayDate + '_FIONA.txt','w')
        
        for n in range(len(self.Section_FIONA_xpos_WeightedAve)):
            for k in range(len(self.Section_FIONA_xpos[n])):
                ff.write(str(n+1) + ' ' + str(self.Section_FIONA_xpos[n][k]) + ' ' + str(self.Section_FIONA_ypos[n][k]) + '\n'   )    
            
        ff.close()       
    


        
    
        ff = open('SectionData_' + figFileName + todayDate + '_shrimp.txt','w')
        
        for n in range(len(self.Section_AFshrimp_xpos)):
            for k in range(len(self.Section_AFshrimp_xpos[n])):
                for p in range(len(self.Section_AFshrimp_xpos[n][k])):
                    ff.write(str(n+1) + ' ' + str(self.Section_AFshrimp_CorrespondingSectionN[n][k]+1) + ' ' + str(self.Section_AFshrimp_xpos[n][k][p]) + ' ' + str(self.Section_AFshrimp_ypos[n][k][p]) + '\n'   )    
            
        ff.close()       
    



        
        
        
        
        
        
        
        
        
        
        



    def PlotFinalEmissionSites(self, event):
        print '\nPlotFinalEmissionSites '
        

        finalEmissionSitesX = []
        finalEmissionSitesY = []
        finalEmissionSitesErr = []
        
        for n in range(self.NTotalEmissionSites):
            if (str(self.EmissionSiteX[n].GetValue()) != '0' ) & (str(self.EmissionSiteY[n].GetValue()) != '0' ):
                finalEmissionSitesX.append(float(self.EmissionSiteX[n].GetValue()))
                finalEmissionSitesY.append(float(self.EmissionSiteY[n].GetValue()))
                finalEmissionSitesErr.append(float(self.EmissionSiteErrR[n].GetValue()))
                        

        #rms = RMS_distance_XY(finalEmissionSitesX, finalEmissionSitesY)
        
        rms2, rms2_Err = RMS_distance_XY_Error(finalEmissionSitesX, finalEmissionSitesY, finalEmissionSitesErr)
        
        rms2_nm = rms2 * float(self.M_nm_pixel.GetValue())
        rms2_Err_nm = rms2_Err * float(self.M_nm_pixel.GetValue())



        self.fig_PlotFinalEmissionSites = plt.figure('PlotFinalEmissionSites', figsize = (16, 8))
        self.fig_PlotFinalEmissionSites.text(0.02, 0.92, '         x,y: '+self.MoleculeXYtext, ha="left", va="bottom", size=12, weight = 'normal', color="red")   
        
        plt.subplot(1,2,1, aspect='equal')
        #mngr = plt.get_current_fig_manager()
        #mngr.window.setGeometry(50,30,1600, 800)        
        plt.suptitle(str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly) + '\n' + self.DataFileAnalysisConditions
            , size=11)
            
        plt.title('Emission Sites , rms: ' + str(np.around(rms2,3)) + ' +- ' + str(np.around(rms2_Err,3))  + ' ,  ' + str(np.around(rms2_nm,1)) + ' +- ' + str(np.around(rms2_Err_nm,1)) + ' nm' )
        

        
        markerColorAFS = ['r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k','r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k','r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k','r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k']
        

  


        
        EmissionSite_Cir = [[] for _ in xrange(self.NTotalEmissionSites)]
        
        for n in range(self.NTotalEmissionSites):
            if (str(self.EmissionSiteX[n].GetValue()) != '0' ) & (str(self.EmissionSiteY[n].GetValue()) != '0' ):
                plt.scatter(float(self.EmissionSiteX[n].GetValue()), float(self.EmissionSiteY[n].GetValue()), marker = 'x',s= 100, edgecolors = markerColorAFS[n], c = markerColorAFS[n]
                , label = '# ' + str(n+1) + ' (' + str(float(self.EmissionSiteX[n].GetValue())) + ', ' + str(float(self.EmissionSiteY[n].GetValue())) + ')   E ' + str(float(self.EmissionSiteErrR[n].GetValue()) ))
                
                EmissionSite_Cir[n] = plt.Circle((float(self.EmissionSiteX[n].GetValue()), float(self.EmissionSiteY[n].GetValue())), radius = float(self.EmissionSiteErrR[n].GetValue())
                , color = markerColorAFS[n], fill=False)
                
                self.fig_PlotFinalEmissionSites.gca().add_artist(EmissionSite_Cir[n])
                
                
                
            else:
                print 'break'
                break
            
            
            

        plt.legend(loc='upper right', prop={'size':12})
        plt.xlim( np.mean(finalEmissionSitesX) - 2, np.mean(finalEmissionSitesX) + 3 ) 
        plt.ylim( np.mean(finalEmissionSitesY) - 2, np.mean(finalEmissionSitesY) + 3 )
            




        plt.subplot(1,2,2, aspect='equal')
        #mngr = plt.get_current_fig_manager()
        #mngr.window.setGeometry(50,30,1600, 800)        
        plt.title('Comparing FIONA-only & AF-shrimp with FIONA')
        plt.suptitle(str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly) + '\n' + self.DataFileAnalysisConditions
            , size=11)
        #self.fig_PlotFinalEmissionSites.text(0.02, 0.91, '         x,y : '+self.MoleculeXYtext, ha="left", va="bottom", size=12, weight = 'normal', color="red")   

        
        markerColorAFS = ['r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k','r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k','r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k','r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k']
        
        
        EmissionSite_Cir = [[] for _ in xrange(self.NTotalEmissionSites)]
        
        for n in range(self.NTotalEmissionSites):
            if (str(self.EmissionSiteX[n].GetValue()) != '0' ) & (str(self.EmissionSiteY[n].GetValue()) != '0' ):
                plt.scatter(float(self.EmissionSiteX[n].GetValue()), float(self.EmissionSiteY[n].GetValue()), marker = 'x',s= 100, edgecolors = markerColorAFS[n], c = markerColorAFS[n]
                , label = '# ' + str(n+1) + ' (' + str(float(self.EmissionSiteX[n].GetValue())) + ', ' + str(float(self.EmissionSiteY[n].GetValue())) + ')   E ' + str(float(self.EmissionSiteErrR[n].GetValue()) ))
                
                EmissionSite_Cir[n] = plt.Circle((float(self.EmissionSiteX[n].GetValue()), float(self.EmissionSiteY[n].GetValue())), radius = float(self.EmissionSiteErrR[n].GetValue())
                , color = markerColorAFS[n], fill=False)
                
                self.fig_PlotFinalEmissionSites.gca().add_artist(EmissionSite_Cir[n])
                
                
                
            else:
                print 'break'
                break
            
 


           

        FIONA_Cir = [[] for _ in xrange(self.TotalNsection)]
        for n in range(len(self.Section_FIONA_xpos_WeightedAve)):
            print 'n = ', n
            if str(self.Section_FIONA_xpos_WeightedAve[n]) != 'nan':
                print 'Section_FIONA_xpos_WeightedAve[n] ', self.Section_FIONA_xpos_WeightedAve[n], type(self.Section_FIONA_xpos_WeightedAve[n])
                print 'run n = ', n
                plt.scatter(self.Section_FIONA_xpos_WeightedAve[n], self.Section_FIONA_ypos_WeightedAve[n], marker = '.',s= 100, edgecolors = markerColorAFS[n], c = markerColorAFS[n]
                , label = 'S' + str(n+1) + ' ('+str(np.around(self.Section_FIONA_xpos_WeightedAve[n],3))+', ' + str(np.around(self.Section_FIONA_ypos_WeightedAve[n],3)) + ')  E ' + str(np.around(self.Section_FIONA_xypos_WeightedAveErr[n],3))  )
            else:
                print 'pass'
                pass
            FIONA_Cir[n] = plt.Circle((self.Section_FIONA_xpos_WeightedAve[n], self.Section_FIONA_ypos_WeightedAve[n]), radius = self.Section_FIONA_xypos_WeightedAveErr[n], color = markerColorAFS[n], fill=False)
            self.fig_PlotFinalEmissionSites.gca().add_artist(FIONA_Cir[n])
 

            

        plt.legend(loc='upper right', prop={'size':12})
        plt.xlim( np.mean(finalEmissionSitesX) - 2, np.mean(finalEmissionSitesX) + 3 ) 
        plt.ylim( np.mean(finalEmissionSitesY) - 2, np.mean(finalEmissionSitesY) + 3 )
            













    def AdjacentShrimpPlot(self, event):
        print '\n Starting AdjacentShrimpPlot'
        
        self.text_AdjacentShrimpPlot.SetLabel(' ')
        
        

        
        
        
        AFSLastFrameN = np.amax([np.amax(self.AFgshrimp_f1n),np.amax(self.AFgshrimp_f2n)])
        #print 'LastFramN: ' , AFSLastFrameN
        
        
        FIONA_LastFrameN = int(np.amax(self.FIONA_FrameN))
        #print 'FIONA_LastFrameN = ', FIONA_LastFrameN
        FrameIntensity = np.zeros(FIONA_LastFrameN + 1)
        #print 'FrameIntensity = ', FrameIntensity
        FrameIndex = np.arange(FIONA_LastFrameN + 1)
        #print 'FrameIndex ', FrameIndex
        
        for n in range(len(self.FIONA_FrameN)):
            #print 'n = ', n
            FNtemp = int(self.FIONA_FrameN[n])
            #print 'FNtemp =  ', FNtemp
            FrameIntensity[FNtemp] = self.FIONA_centerIntensityMax[n]
            
        
        FrameCorMatrixEcc = np.zeros(shape=(AFSLastFrameN+1, AFSLastFrameN+1))
        
        for n in range(len(self.AFgshrimp_f1n)):
            FrameCorMatrixEcc[int(self.AFgshrimp_f1n[n]),int(self.AFgshrimp_f2n[n])] = self.AFgshrimp_eccentricity[n]
            FrameCorMatrixEcc[int(self.AFgshrimp_f2n[n]),int(self.AFgshrimp_f1n[n])] = self.AFgshrimp_eccentricity[n]
        
        
        



        


    
    
        fx0, fy0 = np.array(self.FIONA_x0fit), np.array(self.FIONA_y0fit)
        fx0Err, fy0Err = np.array(self.FIONA_x0fitImageJErr), np.array(self.FIONA_y0fitImageJErr)
        #fxy0Ecc = np.array(self.FIONA_eccentricity)
            
        
        #xposAve = self.MoleculeXpos
        #yposAve = self.MoleculeYpos


        fxy0FrameN = np.array(self.FIONA_FrameN, dtype=np.uint32)
        
        
        AdjacentShrimpEndFrameN = int(self.AdjacentShrimpEndFrameN.GetValue())
        
        #AdjacentShrimpEndFrameNindex = np.where(fxy0FrameN == AdjacentShrimpEndFrameN)[0][0]
        #print 'fxy0FrameN ', fxy0FrameN
        #print 'AdjacentShrimpEndFrameNindex ', AdjacentShrimpEndFrameNindex
        
        
        
        
        #tempI = AdjacentShrimpEndFrameNindex
        
        #AdjacentShrimpEndFrameFIONAerrAve = np.around( (fx0Err[tempI-1] + fx0Err[tempI] + fx0Err[tempI+1] + fy0Err[tempI-1] + fy0Err[tempI] + fy0Err[tempI+1])*1.414/6.0 , 3)
        #print 'AdjacentShrimpEndFrameFIONAerrAve ', AdjacentShrimpEndFrameFIONAerrAve
        #self.MaxPixel_Tolerance.SetValue(str(AdjacentShrimpEndFrameFIONAerrAve))
        
        #AdjacentShrimpEndFrameECCAve = np.around( (fxy0Ecc[tempI-1] + fxy0Ecc[tempI] + fxy0Ecc[tempI+1])/3.0 , 3)
        #print 'AdjacentShrimpEndFrameECCAve ', AdjacentShrimpEndFrameECCAve
        #self.MaxEccentricity.SetValue(str(AdjacentShrimpEndFrameECCAve))
        
        
        
        
        
        
        
        EccentricityCondition = float(self.MaxEccentricity.GetValue())
        maxPer = float(self.MaxPercentage_Tolerance.GetValue())
        maxPx = float(self.MaxPixel_Tolerance.GetValue())


        
        
        
        self.DataFileAnalysisConditions = ('S/B: ' + str(self.SvsB_Tolerance.GetValue()) + ',   Max Ecc: ' + str(self.MaxEccentricity.GetValue()) +  
        ',   Max A err% : ' + str(self.MaxPercentage_Tolerance.GetValue()) +  ',   Max Px x,y err: ' + str(self.MaxPixel_Tolerance.GetValue()) +
        ',   Max Ajd Next Frames: ' + str(self.MaxN_nextAdjacentFrames.GetValue()) + ',   # foward frames: ' + str(self.Adjacent_N_foward_frames.GetValue()) + 
        ',   Min f1n diff: ' + str(self.Adjacent_f1n_min_diff.GetValue()) + ',   Min f1n-f2n diff: ' + str(self.Adjacent_f1n_f2n_min_diff.GetValue()) +
        ',   F# ' + str(self.AdjacentShrimpStartFrameN.GetValue()) + '-' + str(self.AdjacentShrimpEndFrameN.GetValue()))
        

        # two lines conditions below
        self.DataFileAnalysisConditions2 = ('S/B: ' + str(self.SvsB_Tolerance.GetValue()) + ',   Max Ecc: ' + str(self.MaxEccentricity.GetValue()) +  
        ',   Max A err% : ' + str(self.MaxPercentage_Tolerance.GetValue()) +  ',   Max Px x,y err: ' + str(self.MaxPixel_Tolerance.GetValue()) +
        '\nMax Ajd Next Frames: ' + str(self.MaxN_nextAdjacentFrames.GetValue()) +  ',   # foward frames: ' + str(self.Adjacent_N_foward_frames.GetValue()) +
        ',   Min f1n diff: ' + str(self.Adjacent_f1n_min_diff.GetValue()) + ',   Min f1n-f2n diff: ' + str(self.Adjacent_f1n_f2n_min_diff.GetValue()) +
        ',   F# ' + str(self.AdjacentShrimpStartFrameN.GetValue()) + '-' + str(self.AdjacentShrimpEndFrameN.GetValue()))
        




        fx1 = []
        fy1 = []
        fx2 = []
        fy2 = []
        fx2Err = []
        fy2Err = []
        fxy2FrameN = []
                   
        
        self.DataFileAnalysis_FIONA_fxy2_MaxIntensity = []
        
        self.DataFileAnalysis_FIONA_fxy2_FN = []
       
        
        for n in range(len(fx0)):
            FionaAfitErrPtemp = np.around(100.0 * self.FIONA_AfitErr[n] / self.FIONA_Afit[n], 0)     

            if self.FIONA_eccentricity[n] <= EccentricityCondition:
                fx1.append(fx0[n])
                fy1.append(fy0[n])
                
                if (FionaAfitErrPtemp <= maxPer) & (self.FIONA_x0fitImageJErr[n] <= maxPx) & (self.FIONA_y0fitImageJErr[n] <= maxPx):
                    fx2.append(fx0[n])
                    fy2.append(fy0[n])
                    fx2Err.append(fx0Err[n])
                    fy2Err.append(fy0Err[n])
                    fxy2FrameN.append(fxy0FrameN[n])
                                        
                    self.DataFileAnalysis_FIONA_fxy2_MaxIntensity.append(self.FIONA_centerIntensityMax[n])
                    self.DataFileAnalysis_FIONA_fxy2_FN.append(int(self.FIONA_FrameN[n]))
                    
         
        self.DataFileAnalysis_fx2 = fx2
        self.DataFileAnalysis_fy2 = fy2
        self.DataFileAnalysis_fx2Err = fx2Err
        self.DataFileAnalysis_fy2Err = fy2Err
        fxy2Err = (np.array(fx2Err) + np.array(fy2Err))/2.0
        
        
        print 'self.DataFileAnalysis_FIONA_fxy2_FN, fx2, fy2 ', self.DataFileAnalysis_FIONA_fxy2_FN, fx2, fy2
        
        print 'self.DataFileAnalysis_fx2Err = fx2Err, self.DataFileAnalysis_fy2Err = fy2Err ', self.DataFileAnalysis_fx2Err, '\n',  self.DataFileAnalysis_fy2Err
        
        

        fx2Ave = np.mean(fx2)
        fy2Ave = np.mean(fy2)
        fx2DistanceFromMean = np.array(fx2) - fx2Ave
        fy2DistanceFromMean = np.array(fy2) - fy2Ave
        fxy2_distance = np.sqrt(fx2DistanceFromMean**2 + fy2DistanceFromMean**2)
        fxy2_distance_STD = np.std(fxy2_distance)
        
        print 'fxy2_distance_STD = ', fxy2_distance_STD
        
        #plt.figure('fxy2_distance')
        #plt.title('fxy2_distance STD = ' + str(fxy2_distance_STD))
        #plt.plot(fxy2_distance)


        fxy2_rms, fxy2_rms_Err = RMS_distance_XY_Error(fx2, fy2, fxy2Err)
        print 'fxy2_rms with err = ', fxy2_rms, fxy2_rms_Err







        ##############################################################
        ## obtaining the section frames from the user input boxes.

        SectionBeginFN = [] # starting frame number for each section, eg) ...[0] value is the starting frome number of the first section
        SectionEndFN = []   
        SectionFramesText = ''
        
        TotalNsection = 0
        
        print '\nlen(self.F_SectionBegin) = ', len(self.F_SectionBegin)
        
        for n in range(len(self.F_SectionBegin)):
            if (self.F_SectionBegin[n].GetValue() == '0') & (self.F_SectionEnd[n].GetValue() == '0' ):
                break
            else:
                TotalNsection += 1
                
        print 'TotalNsection = ',   TotalNsection 
        
        
        
        
        for n in range(TotalNsection):
            SectionBeginFN.append(int(self.F_SectionBegin[n].GetValue() ))
            SectionEndFN.append(int(self.F_SectionEnd[n].GetValue() ) )
            SectionFramesText += str(int(self.F_SectionBegin[n].GetValue() )) + '-' + str(int(self.F_SectionEnd[n].GetValue() )) + ', '
        
        SectionFramesText = SectionFramesText[:-2]
        
        print 'SectionBeginFN ', SectionBeginFN
        print 'SectionEndFN ', SectionEndFN
        
        
        AllSectionFramesNs = []
        for n in range(TotalNsection):
            AllSectionFramesNs.extend(np.arange(SectionBeginFN[n], SectionEndFN[n]+1))
            
        print 'AllSectionFramesNs = ', AllSectionFramesNs
        print 'SectionFramesText ', SectionFramesText
        ## End of obtaining the section frames
        #####################################################'''
        
        
        


        
        
        gx0, gy0 = self.AFgshrimp_x0fitImageJ, self.AFgshrimp_y0fitImageJ
        gx0Err, gy0Err = self.AFgshrimp_x0fitImageJErr, self.AFgshrimp_y0fitImageJErr
        
        gxy0_f1n = self.AFgshrimp_f1n
        gxy0_f2n = self.AFgshrimp_f2n
                
        
        gx1 = []
        gy1 = []
        gx2 = []
        gy2 = []
        gx2Err = []
        gy2Err = []
        gx3 = []
        gy3 = []
                
        
        self.DataFileAnalysis_AFshrimp_gxy2_F1N = []
        self.DataFileAnalysis_AFshrimp_gxy2_F2N = []
        
        
        for n in range(len(gx0)):
            ##############################################################
            ### bleaching type dependent selection 1 for sections ########
            if (self.radio_AdjacentBleachingType.GetSelection() == 1) or (self.radio_AdjacentBleachingType.GetSelection() == 3) or (self.radio_AdjacentBleachingType.GetSelection() == 4): # selecting only within the Sections
                if (gxy0_f1n[n] in AllSectionFramesNs) and (gxy0_f2n[n] in AllSectionFramesNs):
                    pass
                else:
                    continue
            
            ### bleaching type dependent selection 1 for sections ########
            ##############################################################

                
                
            AFshrimpAfitErrPtemp = np.around(100.0 * self.AFgshrimp_AfitErr[n] / self.AFgshrimp_Afit[n], 0)     
            
            f1n_FIONA_indexN = np.where(fxy0FrameN == gxy0_f1n[n])[0][0]
            f2n_FIONA_indexN = np.where(fxy0FrameN == gxy0_f2n[n])[0][0]
            
            #print 'self.AFgshrimp_f1n[n], f2n = ', self.AFgshrimp_f1n[n], self.AFgshrimp_f2n[n]
            
            #print 'f1n_FIONA_indexN, f2n = ', f1n_FIONA_indexN, f2n_FIONA_indexN
            
            
            
            ###############################
            ### 2nd Ecc constraint#########
            if self.radio_AdjacentEcc.GetSelection() == 0 :
                SecondEccConstraint = True
                SecondEccConstraintText = ' : None'
                
            if self.radio_AdjacentEcc.GetSelection() == 1 :
                SecondEccConstraint = (self.AFgshrimp_eccentricity[n] < self.FIONA_eccentricity[f2n_FIONA_indexN]) or (self.AFgshrimp_eccentricity[n] < self.FIONA_eccentricity[f1n_FIONA_indexN])
                SecondEccConstraintText = ' : Either'
                            
            if self.radio_AdjacentEcc.GetSelection() == 2 :
                SecondEccConstraint = (self.AFgshrimp_eccentricity[n] < self.FIONA_eccentricity[f2n_FIONA_indexN]) and (self.AFgshrimp_eccentricity[n] < self.FIONA_eccentricity[f1n_FIONA_indexN])
                SecondEccConstraintText = ' : Both'
                
            
            #print 'SecondEccConstraint ', SecondEccConstraint
            
            
            if self.AFgshrimp_eccentricity[n] <= EccentricityCondition:
                gx1.append(gx0[n])
                gy1.append(gy0[n])
                
                #if (AFshrimpAfitErrPtemp <= maxPer) & (self.AFgshrimp_x0fitImageJErr[n] <= maxPx) & (self.AFgshrimp_y0fitImageJErr[n] <= maxPx):
                if (AFshrimpAfitErrPtemp <= maxPer) & (self.AFgshrimp_x0fitImageJErr[n] <= maxPx) & (self.AFgshrimp_y0fitImageJErr[n] <= maxPx) & SecondEccConstraint :
                    #print 'gx2 done '
                    gx2.append(gx0[n])
                    gy2.append(gy0[n])
                    gx2Err.append(gx0Err[n])
                    gy2Err.append(gy0Err[n])
                    
                    self.DataFileAnalysis_AFshrimp_gxy2_F1N.append(int(gxy0_f1n[n]))
                    self.DataFileAnalysis_AFshrimp_gxy2_F2N.append(int(gxy0_f2n[n]))
                    
                    
                    if (self.AFgshrimp_eccentricity[n] < self.FIONA_eccentricity[f1n_FIONA_indexN]) & (self.AFgshrimp_eccentricity[n] < self.FIONA_eccentricity[f2n_FIONA_indexN]):
                        #print 'gx3 done '
                        gx3.append(gx0[n])
                        gy3.append(gy0[n])
            
            ### 2nd Ecc constraint#########
            ###############################
        
        
        
        #print 'gx2 \n ', gx2
        #print 'gx3 \n', gx3
        
        

        


      
        self.DataFileAnalysis_gx2 = gx2
        self.DataFileAnalysis_gy2 = gy2
        self.DataFileAnalysis_gx2Err = gx2Err
        self.DataFileAnalysis_gy2Err = gy2Err
        
        gxy2ErrAve = ( np.array(gx2Err) + np.array(gy2Err) )/2.0





        #print 'DataFileAnalysis_AFshrimp_gxy2_F1N ', self.DataFileAnalysis_AFshrimp_gxy2_F1N
        #print 'DataFileAnalysis_AFshrimp_gxy2_F2N ', self.DataFileAnalysis_AFshrimp_gxy2_F2N





        #print 'gx1', gx1
        #print '\ngx2', gx2
        print 'self.DataFileAnalysis_AFshrimp_gxy2_F1N ', self.DataFileAnalysis_AFshrimp_gxy2_F1N
        print 'self.DataFileAnalysis_AFshrimp_gxy2_F2N ', self.DataFileAnalysis_AFshrimp_gxy2_F2N



        ####################################################################################################################
        ####  collecting adjacent afShrimp temp data with constraints.

        adjacent_gx2temp = []
        adjacent_gy2temp = []
        adjacent_gxy2Errtemp = []
        adjacentFrameN1temp = []
        adjacentFrameN2temp = []
        
        #adjacentNumberTemp = int(self.MaxN_nextAdjacentFrames.GetValue())        
        
        #F1Ntemp = -99999

        for n in range(len(self.DataFileAnalysis_AFshrimp_gxy2_F1N)): 
            #print 'AFshrimp_gxy2_F1N = ', n
            #print self.DataFileAnalysis_AFshrimp_gxy2_F1N[n], self.DataFileAnalysis_AFshrimp_gxy2_F2N[n]
            
            #if F1Ntemp == self.DataFileAnalysis_AFshrimp_gxy2_F1N[n]: # excluding double selection for the same f1n
            #    continue
            
            #if F1Ntemp + self.Adjacent_f1n_min_diff.GetValue() > self.DataFileAnalysis_AFshrimp_gxy2_F1N[n]  : # excluding selecting mext immediate a few frames for the same f1n
            #    continue
            
            
            if self.DataFileAnalysis_AFshrimp_gxy2_F1N[n] < int(self.AdjacentShrimpStartFrameN.GetValue()): # limiting f1n from the AdjacentStartFrame
                continue                
                
            if self.DataFileAnalysis_AFshrimp_gxy2_F2N[n] > int(self.AdjacentShrimpEndFrameN.GetValue()): # limiting f2n to the AdjacentEndFrame
                continue

            #print '##'
            #print 'self.DataFileAnalysis_AFshrimp_gxy2_F1N[n] = ', self.DataFileAnalysis_AFshrimp_gxy2_F1N[n]
            #print 'self.Adjacent_f1n_f2n_min_diff.GetValue() = ', self.Adjacent_f1n_f2n_min_diff.GetValue()            
            #print 'self.DataFileAnalysis_AFshrimp_gxy2_F2N[n] = ', self.DataFileAnalysis_AFshrimp_gxy2_F2N[n]
            
            #print 'self.DataFileAnalysis_AFshrimp_gxy2_F1N[n] + int(self.Adjacent_f1n_f2n_min_diff.GetValue()) ', self.DataFileAnalysis_AFshrimp_gxy2_F1N[n] + int(self.Adjacent_f1n_f2n_min_diff.GetValue())
            #print '- = ', self.DataFileAnalysis_AFshrimp_gxy2_F1N[n] + int(self.Adjacent_f1n_f2n_min_diff.GetValue()) - self.DataFileAnalysis_AFshrimp_gxy2_F2N[n]
            
            if self.DataFileAnalysis_AFshrimp_gxy2_F1N[n] + int(self.Adjacent_f1n_f2n_min_diff.GetValue()) > self.DataFileAnalysis_AFshrimp_gxy2_F2N[n]:
                continue  # excluding close f1n , f2n pairs. min separation is the user defined frame difference.
                
            if self.DataFileAnalysis_AFshrimp_gxy2_F1N[n] + int(self.MaxN_nextAdjacentFrames.GetValue()) < self.DataFileAnalysis_AFshrimp_gxy2_F2N[n]:
                continue  # limiting number of next adjacent frames.





                
            
            adjacent_gx2temp.append(gx2[n])
            adjacent_gy2temp.append(gy2[n])
            adjacent_gxy2Errtemp.append(gxy2ErrAve[n])
            adjacentFrameN1temp.append(self.DataFileAnalysis_AFshrimp_gxy2_F1N[n])
            adjacentFrameN2temp.append(self.DataFileAnalysis_AFshrimp_gxy2_F2N[n])


            '''
            if self.DataFileAnalysis_AFshrimp_gxy2_F1N[n] >= int(self.AdjacentShrimpStartFrameN.GetValue()):
                    
                for m in range(adjacentNumberTemp):
                    # below m+self.Adjacent_f1n_f2n_min_diff.GetValue() is for excluding close frame pairs
                    if (self.DataFileAnalysis_AFshrimp_gxy2_F1N[n] + m + self.Adjacent_f1n_f2n_min_diff.GetValue() >= int(self.AdjacentShrimpEndFrameN.GetValue())):
                        break
                    
                    if (gx2[n] < xposAve-2) or (gx2[n] > xposAve + 3) or (gy2[n] < yposAve -2) or (gy2[n] > yposAve +3): # excluding unrealistic values
                        continue
                  
                    if (self.DataFileAnalysis_AFshrimp_gxy2_F1N[n] + m+self.Adjacent_f1n_f2n_min_diff.GetValue() == self.DataFileAnalysis_AFshrimp_gxy2_F2N[n]):
                        print '\n ok => f1n, f2n : ', n, n+m+self.Adjacent_f1n_f2n_min_diff.GetValue() 
                        print 'F1N F2N = ', self.DataFileAnalysis_AFshrimp_gxy2_F1N[n], self.DataFileAnalysis_AFshrimp_gxy2_F2N[n]
                        adjacent_gx2temp.append(gx2[n])
                        adjacent_gy2temp.append(gy2[n])
                        adjacent_gxy2Errtemp.append(gxy2ErrAve[n])
                        adjacentFrameN1temp.append(self.DataFileAnalysis_AFshrimp_gxy2_F1N[n])
                        adjacentFrameN2temp.append(self.DataFileAnalysis_AFshrimp_gxy2_F2N[n])
                        
                        F1Ntemp = self.DataFileAnalysis_AFshrimp_gxy2_F1N[n] # excluding double selection for the same f1n by the 1st if statement above.
                        
               '''         

        
        print '\n\nadjacentFrameN1temp ', adjacentFrameN1temp
        print 'adjacentFrameN2temp ', adjacentFrameN2temp
        #print 'adjacent_gx2temp ', adjacent_gx2temp        
        
        print 'End of collecting adjacent afShrimp temp data. \n '

        ####  collecting adjacent afShrimp temp data.
        ####################################################################################################################












        ####################################################################################################################
        ####  collecting adjacent afShrimp temp data.



        adjacent_gx2 = []
        adjacent_gy2 = []
        adjacent_gxy2Err = []
        adjacentFrameN1 = []
        adjacentFrameN2 = []


        #############################################
        ### bleaching type dependent selection 2 ########

        if self.radio_AdjacentBleachingType.GetSelection() == 0 or  self.radio_AdjacentBleachingType.GetSelection() == 1: # for random bleaching
            if self.radio_AdjacentBleachingType.GetSelection() == 0:
                AdjacentBleachingTypeText = ' : Random'
            if self.radio_AdjacentBleachingType.GetSelection() == 1:
                AdjacentBleachingTypeText = ' : Random,  Sections: ' + SectionFramesText
                
            print 'AdjacentBleachingTypeText = ', AdjacentBleachingTypeText

            for n in range(len(adjacentFrameN2temp)): # excluding multiple selections for f2n, selecting nearest f1n, f2n.  This is for randome bleaching
                print '\n n = ', n
                print 'adjacentFrameN2temp.count(adjacentFrameN2temp[n]) ', adjacentFrameN2temp.count(adjacentFrameN2temp[n])
                

                if adjacentFrameN2temp.count(adjacentFrameN2temp[n]) == 1: #if f2n is only one, include it.
                    adjacent_gx2.append(adjacent_gx2temp[n])
                    adjacent_gy2.append(adjacent_gy2temp[n])
                    adjacent_gxy2Err.append(adjacent_gxy2Errtemp[n])
                    adjacentFrameN1.append(adjacentFrameN1temp[n])
                    adjacentFrameN2.append(adjacentFrameN2temp[n])
                    continue


                if adjacentFrameN2temp[n] in adjacentFrameN2: # if there is a same f2n in the already selceted list, choosing only min err one.
                    previous_f2nIndexTemp = adjacentFrameN2.index(adjacentFrameN2temp[n])
                    previous_f2nErrTemp = adjacent_gxy2Err[previous_f2nIndexTemp]
                    if adjacent_gxy2Errtemp[n] < previous_f2nErrTemp:
                        del adjacent_gx2[previous_f2nIndexTemp]
                        del adjacent_gy2[previous_f2nIndexTemp]
                        del adjacent_gxy2Err[previous_f2nIndexTemp]
                        del adjacentFrameN1[previous_f2nIndexTemp]
                        del adjacentFrameN2[previous_f2nIndexTemp]
                        
                        adjacent_gx2.append(adjacent_gx2temp[n])
                        adjacent_gy2.append(adjacent_gy2temp[n])
                        adjacent_gxy2Err.append(adjacent_gxy2Errtemp[n])
                        adjacentFrameN1.append(adjacentFrameN1temp[n])
                        adjacentFrameN2.append(adjacentFrameN2temp[n])                        
                        
                        
                    else:
                        continue
                

                
                
                ''' # not yet done
                else: #selecting nearest f1n, f2n pair if thre are more than two pairs with same f2n.
                    indexTemp = len(adjacentFrameN2temp) -1 - adjacentFrameN2temp[::-1].index(adjacentFrameN2temp[n])
                    
                    adjacent_gx2.append(adjacent_gx2temp[indexTemp])
                    adjacent_gy2.append(adjacent_gy2temp[indexTemp])
                    adjacent_gxy2Err.append(adjacent_gxy2Errtemp[indexTemp])
                    adjacentFrameN1.append(adjacentFrameN1temp[indexTemp])
                    adjacentFrameN2.append(adjacentFrameN2temp[indexTemp])
                    
                    
 
                    absFrameDiffTemp = np.abs(adjacentFrameN1temp[n] - adjacentFrameN2temp[n])
                    indexTemp = n
                    print ''
                    for k in range(len(adjacentFrameN2temp)):
                        if adjacentFrameN2temp[k] == adjacentFrameN2temp[n]:
                            absFrameDiffTemp2 = np.abs(adjacentFrameN1temp[k] - adjacentFrameN2temp[k])
                            if absFrameDiffTemp > absFrameDiffTemp2:
                                indexTemp = k
                    
                    

                if adjacentFrameN2temp[n] > adjacentFrameN2[-1] + 1: #excluding immediate next frame
                    adjacent_gx2.append(adjacent_gx2temp[n])
                    adjacent_gy2.append(adjacent_gy2temp[n])
                    adjacent_gxy2Err.append(adjacent_gxy2Errtemp[n])
                    adjacentFrameN1.append(adjacentFrameN1temp[n])
                    adjacentFrameN2.append(adjacentFrameN2temp[n])
                    
                else:
                    pass            
                
                '''


            for n in range(len(adjacentFrameN2)-1): # choosing only one when there are adjacent two points for f2n.  ex), f2n is 33,34,35 => choosing min Err pair.
                if adjacentFrameN2[n]-1 in adjacentFrameN2 and adjacentFrameN2[n]+1 in adjacentFrameN2:
                    PreviousOneIndexTemp = adjacentFrameN2.index(adjacentFrameN2[n]-1)
                    NextOneIndexTemp = adjacentFrameN2.index(adjacentFrameN2[n]+1)
                    errsTemp = [adjacent_gxy2Err[PreviousOneIndexTemp], adjacent_gxy2Err[n], adjacent_gxy2Err[NextOneIndexTemp]]
                    errsTempMinIndex = np.argmin(errsTemp)
                    
                    if errsTempMinIndex == 0:
                        del adjacent_gx2[NextOneIndexTemp]
                        del adjacent_gy2[NextOneIndexTemp]
                        del adjacent_gxy2Err[NextOneIndexTemp]
                        del adjacentFrameN1[NextOneIndexTemp]
                        del adjacentFrameN2[NextOneIndexTemp]
                    
                        del adjacent_gx2[n]
                        del adjacent_gy2[n]
                        del adjacent_gxy2Err[n]
                        del adjacentFrameN1[n]
                        del adjacentFrameN2[n]
                        
                        
                
                    elif errsTempMinIndex == 1:
                        del adjacent_gx2[NextOneIndexTemp]
                        del adjacent_gy2[NextOneIndexTemp]
                        del adjacent_gxy2Err[NextOneIndexTemp]
                        del adjacentFrameN1[NextOneIndexTemp]
                        del adjacentFrameN2[NextOneIndexTemp]

                        del adjacent_gx2[PreviousOneIndexTemp]
                        del adjacent_gy2[PreviousOneIndexTemp]
                        del adjacent_gxy2Err[PreviousOneIndexTemp]
                        del adjacentFrameN1[PreviousOneIndexTemp]
                        del adjacentFrameN2[PreviousOneIndexTemp]
                        
                    
                
                    elif errsTempMinIndex == 2:
                        del adjacent_gx2[n]
                        del adjacent_gy2[n]
                        del adjacent_gxy2Err[n]
                        del adjacentFrameN1[n]
                        del adjacentFrameN2[n]
           
                        del adjacent_gx2[PreviousOneIndexTemp]
                        del adjacent_gy2[PreviousOneIndexTemp]
                        del adjacent_gxy2Err[PreviousOneIndexTemp]
                        del adjacentFrameN1[PreviousOneIndexTemp]
                        del adjacentFrameN2[PreviousOneIndexTemp]
                        
     



            for n in range(len(adjacentFrameN2)-1): # choosing only one when there are adjacent another point for f2n.  ex), f2n is 33,34 => choosing min Err pair.
                if adjacentFrameN2[n]+1 in adjacentFrameN2:
                    NextOneIndexTemp = adjacentFrameN2.index(adjacentFrameN2[n]+1)
                    errsTemp = [adjacent_gxy2Err[n], adjacent_gxy2Err[NextOneIndexTemp]]
                    errsTempMinIndex = np.argmin(errsTemp)
                    
                    if errsTempMinIndex == 0:
                        del adjacent_gx2[NextOneIndexTemp]
                        del adjacent_gy2[NextOneIndexTemp]
                        del adjacent_gxy2Err[NextOneIndexTemp]
                        del adjacentFrameN1[NextOneIndexTemp]
                        del adjacentFrameN2[NextOneIndexTemp]
                    
                       
                
                    elif errsTempMinIndex == 1:
                        del adjacent_gx2[n]
                        del adjacent_gy2[n]
                        del adjacent_gxy2Err[n]
                        del adjacentFrameN1[n]
                        del adjacentFrameN2[n]











        
        if self.radio_AdjacentBleachingType.GetSelection() == 2: # for monotonic continuous bleaching
            AdjacentBleachingTypeText = ' : Monotonic'
            print 'AdjacentBleachingTypeText = ', AdjacentBleachingTypeText        
            
            #adjacent_gx2.append(adjacent_gx2temp[0])
            #adjacent_gy2.append(adjacent_gy2temp[0])
            #adjacent_gxy2Err.append(adjacent_gxy2Errtemp[0])
            #adjacentFrameN1.append(adjacentFrameN1temp[0])
            #adjacentFrameN2.append(adjacentFrameN2temp[0])
            
            for n in range(len(adjacentFrameN2temp)): # excluding multiple selections for f2n, this is for monotonic bleaching
                print '\n monotonic bleaching: n = ', n
                print 'adjacentFrameN2temp.count(adjacentFrameN2temp[n]) ', adjacentFrameN2temp.count(adjacentFrameN2temp[n])
                
                if adjacentFrameN1temp[n] in adjacentFrameN1:
                    continue
                
                try:
                    if adjacentFrameN1temp[n] - adjacentFrameN1[-1] < self.Adjacent_f1n_min_diff.GetValue():
                        continue
                except: pass
            
                try:
                    if adjacentFrameN2temp[n] - adjacentFrameN2[-1] < self.Adjacent_f1n_min_diff.GetValue():
                        continue
                except: pass
            
                    
                
                gx2temp = adjacent_gx2temp[n]
                gy2temp = adjacent_gy2temp[n]
                gxy2Errtemp = adjacent_gxy2Errtemp[n]
                f1ntemp = adjacentFrameN1temp[n]
                f2ntemp = adjacentFrameN2temp[n]
                
                gx2tempBest = gx2temp
                gy2tempBest = gy2temp
                gxy2ErrtempBest = gxy2Errtemp
                f1ntempBest = f1ntemp
                f2ntempBest = f2ntemp
                
                NadjacentForwardF = self.Adjacent_N_foward_frames.GetValue()
                
                
                for k in range(len(adjacentFrameN1temp)): 
                # finding a min Err pair within next 3 frames for both f1n and f2n. ex) if f1n=3 and f2n=10, a best pair for f1n=3-6 and f2n=10-13
                    if adjacentFrameN1temp[k] < f1ntemp or adjacentFrameN2temp[k] < f2ntemp or adjacentFrameN1temp[k] > f1ntemp + NadjacentForwardF -1 or adjacentFrameN2temp[k] > f2ntemp + NadjacentForwardF - 1:
                        continue
                    #elif adjacentFrameN1temp[k] <= f1ntemp + 3 and adjacentFrameN2temp[k] <= f2ntemp +3:
                    if adjacent_gxy2Errtemp[k] < gxy2Errtemp:
                        gx2tempBest = adjacent_gx2temp[k]
                        gy2tempBest = adjacent_gy2temp[k]
                        gxy2ErrtempBest = adjacent_gxy2Errtemp[k]
                        f1ntempBest = adjacentFrameN1temp[k]
                        f2ntempBest = adjacentFrameN2temp[k]
                
                    
                    

                adjacent_gx2.append(gx2tempBest)
                adjacent_gy2.append(gy2tempBest)
                adjacent_gxy2Err.append(gxy2ErrtempBest)
                adjacentFrameN1.append(f1ntempBest)
                adjacentFrameN2.append(f2ntempBest)
                    








            
        if self.radio_AdjacentBleachingType.GetSelection() == 3: # for monotonic stepwise bleaching within the sections, only one pair in neiboring sections
            AdjacentBleachingTypeText = ' : Monotonic,  Sections: ' + SectionFramesText
            print 'AdjacentBleachingTypeText = ', AdjacentBleachingTypeText        
            
            
            f1n_temp = -99999
            f1n_f2n_distance_temp = 99999
            for n in range(len(adjacentFrameN1temp)): # selecting only one f1n point in each section
                print '\n monotonic bleaching: n = ', n
                print 'adjacentFrameN2temp.count(adjacentFrameN2temp[n]) ', adjacentFrameN2temp.count(adjacentFrameN2temp[n])
                
                print '\n adjacentFrameN1temp[n], adjacentFrameN2temp[n]: ', adjacentFrameN1temp[n], adjacentFrameN2temp[n]




                f1n_nextf1n_temp = ' '
                for k in range(len(SectionBeginFN)): # excluding multiple f1n in a same section
                    print 'range(SectionBeginFN[k], SectionEndFN[k]+1) ', range(SectionBeginFN[k], SectionEndFN[k]+1) 
                    if ( adjacentFrameN1temp[n] in range(SectionBeginFN[k], SectionEndFN[k]+1) ) and ( f1n_temp in range(SectionBeginFN[k], SectionEndFN[k]+1) ) :
                        #if adjacentFrameN2temp[n] - adjacentFrameN1temp[n] < f1n_f2n_distance_temp:
                        if adjacent_gxy2Errtemp[n] < f1n_f2n_distance_temp:
                            f1n_nextf1n_temp = 'sameFrameShorterDistance'
                            
                            del adjacent_gx2[-1]
                            del adjacent_gy2[-1]
                            del adjacent_gxy2Err[-1]
                            del adjacentFrameN1[-1]
                            del adjacentFrameN2[-1]
                            
                        else:
                            f1n_nextf1n_temp = 'sameFrameNotShorterDistance'
                            
                            

                        
                if f1n_nextf1n_temp == 'sameFrameNotShorterDistance':
                    print 'sameFrameNotShorterDistance: skipped'
                    continue




                f1n_f2n_temp = ' '
                for k in range(len(SectionBeginFN)): # excluding same section f1n f2n
                    if ( adjacentFrameN1temp[n] in range(SectionBeginFN[k], SectionEndFN[k]+1) ) and ( adjacentFrameN2temp[n] in range(SectionBeginFN[k], SectionEndFN[k]+1) ) :
                        f1n_f2n_temp = 'sameFrame'
                        
                if f1n_f2n_temp == 'sameFrame':
                    print 'sameframe2'
                    continue
                
                else:
                    adjacent_gx2.append(adjacent_gx2temp[n])
                    adjacent_gy2.append(adjacent_gy2temp[n])
                    adjacent_gxy2Err.append(adjacent_gxy2Errtemp[n])
                    adjacentFrameN1.append(adjacentFrameN1temp[n])
                    adjacentFrameN2.append(adjacentFrameN2temp[n])
                    
                    f1n_temp = adjacentFrameN1temp[n]
                    #f1n_f2n_distance_temp = adjacentFrameN2temp[n] - adjacentFrameN1temp[n]
                    f1n_f2n_distance_temp = adjacent_gxy2Errtemp[n]
                    
                        






        if self.radio_AdjacentBleachingType.GetSelection() == 4: # for monotonic stepwise bleaching within the sections, only one pair in neiboring sections
            AdjacentBleachingTypeText = ' : Monotonic,  Sections all points: ' + SectionFramesText
            print 'AdjacentBleachingTypeText = ', AdjacentBleachingTypeText        
            
            
            #f1n_temp = -99999
            for n in range(len(adjacentFrameN1temp)): # selecting only one f1n point in each section
                print '\n monotonic bleaching: n = ', n
                print 'adjacentFrameN2temp.count(adjacentFrameN2temp[n]) ', adjacentFrameN2temp.count(adjacentFrameN2temp[n])
                
                print '\n adjacentFrameN1temp[n], adjacentFrameN2temp[n]: ', adjacentFrameN1temp[n], adjacentFrameN2temp[n]


                f1n_f2n_temp = ' '
                for k in range(len(SectionBeginFN)-1): # excluding multiple f1n in a same section
                    print 'range(SectionBeginFN[k], SectionEndFN[k]+1) ', range(SectionBeginFN[k], SectionEndFN[k]+1) 
                    if ( adjacentFrameN1temp[n] in range(SectionBeginFN[k], SectionEndFN[k]+1) ) and ( adjacentFrameN2temp[n] in range(SectionBeginFN[k+1], SectionEndFN[k+1]+1) ) :
                        f1n_f2n_temp = 'adjacentSections'
                        
                if f1n_f2n_temp == 'adjacentSections':
                    print 'adjacentSections '
                    pass
                else:
                    print 'Not adjacentSections'
                    continue




                f1n_f2n_temp = ' '
                for k in range(len(SectionBeginFN)): # excluding same section f1n f2n
                    if ( adjacentFrameN1temp[n] in range(SectionBeginFN[k], SectionEndFN[k]+1) ) and ( adjacentFrameN2temp[n] in range(SectionBeginFN[k], SectionEndFN[k]+1) ) :
                        f1n_f2n_temp = 'sameFrame'
                        
                if f1n_f2n_temp == 'sameFrame':
                    print 'sameframe2'
                    continue
                
                else:
                    adjacent_gx2.append(adjacent_gx2temp[n])
                    adjacent_gy2.append(adjacent_gy2temp[n])
                    adjacent_gxy2Err.append(adjacent_gxy2Errtemp[n])
                    adjacentFrameN1.append(adjacentFrameN1temp[n])
                    adjacentFrameN2.append(adjacentFrameN2temp[n])
                    
                    #f1n_temp = adjacentFrameN1temp[n]
                    
                        
                            

                            
                        
        ### bleaching type dependent selection 2 ########
        #############################################
        
        
        
       
       
       
       
        
        
        print 'fxy2FrameN = ',  fxy2FrameN     
        
        print 'int(self.AdjacentShrimpStartFrameN.GetValue()) = ', int(self.AdjacentShrimpStartFrameN.GetValue())
        print 'int(self.AdjacentShrimpEndFrameN.GetValue()) = ', int(self.AdjacentShrimpEndFrameN.GetValue())
        

        ######################################
        #### finding the last FIONA location .

        if self.radio_AdjacentBleachingType.GetSelection() == 0 or self.radio_AdjacentBleachingType.GetSelection() == 2:   # without sections     
        
            
            ### finding a fiona point at the last shrimp frame.
            try:
                adjacent_shrimp_final_FIONA_fxy2_indexN = fxy2FrameN.index(adjacentFrameN2[-1])
                LastFIONAframeN = int(fxy2FrameN[adjacent_shrimp_final_FIONA_fxy2_indexN])
                #adjacent_shrimp_final_FIONA_fxy2_ErrTemp = fxy2Err[adjacent_shrimp_final_FIONA_fxy2_indexN]
            except:
            
                for n2 in range(AdjacentShrimpEndFrameN - adjacentFrameN2[-1]):
                    try:
                        adjacent_shrimp_final_FIONA_fxy2_indexN = fxy2FrameN.index(adjacentFrameN2[-1]+n2+1)
                        LastFIONAframeN = int(fxy2FrameN[adjacent_shrimp_final_FIONA_fxy2_indexN])
                        break
                    except:
                        pass
                    
        
            
        
            '''
            ### finding a fiona point near the last frame defined by a user input.
            adjacent_shrimp_final_FIONA_fxy2_indexN = 0
            adjacent_shrimp_final_FIONA_fxy2_ErrTemp = 99999
            LastFIONAframeN = 99999            

            for n in range(len(fxy2FrameN)):  
                if int(self.AdjacentShrimpEndFrameN.GetValue()) <= fxy2FrameN[n] <= int(self.AdjacentShrimpEndFrameN.GetValue()) + int(self.MaxN_nextAdjacentFrames.GetValue()):
                    if fxy2Err[n] < adjacent_shrimp_final_FIONA_fxy2_ErrTemp:
                        adjacent_shrimp_final_FIONA_fxy2_indexN = n
                        LastFIONAframeN = int(fxy2FrameN[n])
                        adjacent_shrimp_final_FIONA_fxy2_ErrTemp = fxy2Err[n]
            '''            
                
        
        if self.radio_AdjacentBleachingType.GetSelection() == 1 or self.radio_AdjacentBleachingType.GetSelection() == 3 or self.radio_AdjacentBleachingType.GetSelection() == 4: # with sections
            lastSectionFrameIndex = len(SectionBeginFN) - 1
            lastSectionFrameBeginN = SectionBeginFN[lastSectionFrameIndex]
            lastSectionFrameEndN = SectionEndFN[lastSectionFrameIndex]

            fxy2ErrTemp = 100
            for n in range(len(fx2)):
                if fxy2FrameN[n] in range(lastSectionFrameBeginN,lastSectionFrameEndN+1):
                    if fxy2Err[n] < fxy2ErrTemp:
                        adjacent_shrimp_final_FIONA_fxy2_indexN = n
                        LastFIONAframeN = int(fxy2FrameN[n])                       
                        fxy2ErrTemp = fxy2Err[n]
                

        #### finding the last FIONA location .
        ######################################


   
        
        
        try: 
            print 'LastFIONAframeN = ', LastFIONAframeN
            print 'adjacent_shrimp_final_FIONA_fxy2_indexN = ', adjacent_shrimp_final_FIONA_fxy2_indexN
        except:
            self.text_AdjacentShrimpPlot.SetLabel('Error: End Frame')
            
        print '\nadjacent_gx2:\n ', adjacent_gx2
        print '\nadjacent_gy2:\n ', adjacent_gy2     
        print ' adjacentFrame N1 N2 = ', adjacentFrameN1, adjacentFrameN2
        
        print '\nfxy2FrameN ', fxy2FrameN
        
        
        adjacent_gx2.append(fx2[adjacent_shrimp_final_FIONA_fxy2_indexN])
        adjacent_gy2.append(fy2[adjacent_shrimp_final_FIONA_fxy2_indexN])
        adjacent_gxy2Err.append(fxy2Err[adjacent_shrimp_final_FIONA_fxy2_indexN])
        adjacentFrameN1.append(LastFIONAframeN)
        adjacentFrameN2.append(LastFIONAframeN)
         

                
        print '\nadjacent_gx2:\n ', adjacent_gx2
        print '\nadjacent_gy2:\n ', adjacent_gy2
        
        for n in range(len(adjacentFrameN1)):
            print ' adjacentFrame N1 N2 = ', adjacentFrameN1[n], adjacentFrameN2[n]
            
            
            
           
            
            
        NAdjacentPairs = len(adjacentFrameN2)    
            
        
        print 'N Adjacent Pairs = ', NAdjacentPairs 
        ####  collecting adjacent afShrimp data.
        ####################################################################################################################





        print '\nadjacent_gx2:\n ', adjacent_gx2
        print '\nadjacent_gy2:\n ', adjacent_gy2




    

        adjacent_gx2_ave = np.mean(adjacent_gx2)
        adjacent_gy2_ave = np.mean(adjacent_gy2)
        print 'adjacent_gx2_ave = ', adjacent_gx2_ave
        print 'adjacent_gy2_ave = ', adjacent_gy2_ave
        
        adjacent_gxy2_distance = np.sqrt( (np.array(adjacent_gx2)-adjacent_gx2_ave)**2 + (np.array(adjacent_gy2)-adjacent_gy2_ave)**2 )
        print 'adjacent_gxy2_distance = ', adjacent_gxy2_distance
        
        adjacent_gxy2_distance_rms = np.sqrt(np.mean(adjacent_gxy2_distance**2))
        print 'adjacent_gxy2_distance_rms = ', adjacent_gxy2_distance_rms




        adjacent_gxy2_distance_rms = RMS_distance_XY(adjacent_gx2, adjacent_gy2)
        print 'adjacent_gxy2_distance_rms = ', adjacent_gxy2_distance_rms


        rms2, rms2_Err = RMS_distance_XY_Error(adjacent_gx2, adjacent_gy2, adjacent_gxy2Err)
        print 'adjacent_gxy2_distance_rms with err = ', rms2, rms2_Err


        rms2_nm = rms2 * float(self.M_nm_pixel.GetValue())
        rms2_Err_nm = rms2_Err * float(self.M_nm_pixel.GetValue())



        
        #############################################################################
        ## only for 1st last distance.
        if self.radio_options_autoMultiPleAnalysis.GetSelection() == 2:
            dx_1st_last_pair = adjacent_gx2[-1] - adjacent_gx2[0]
            dy_1st_last_pair = adjacent_gy2[-1] - adjacent_gy2[0]
            
            dx_1st_lastN2_pair = adjacent_gx2[-2] - adjacent_gx2[0]
            dy_1st_lastN2_pair = adjacent_gy2[-2] - adjacent_gy2[0]
            
            
            dx_1st_last_pair_nm = dx_1st_last_pair * float(self.M_nm_pixel.GetValue())
            dy_1st_last_pair_nm = dy_1st_last_pair * float(self.M_nm_pixel.GetValue())
            
            
            dx_1st_lastN2_pair_nm = dx_1st_lastN2_pair * float(self.M_nm_pixel.GetValue())
            dy_1st_lastN2_pair_nm = dy_1st_lastN2_pair * float(self.M_nm_pixel.GetValue())
            
            
            print 'dx_1st_last_pair_nm ', dx_1st_last_pair_nm
            print 'dy_1st_last_pair_nm ', dy_1st_last_pair_nm
            
            return dx_1st_last_pair_nm, dy_1st_last_pair_nm, rms2_nm, dx_1st_lastN2_pair_nm, dy_1st_lastN2_pair_nm
            
        ## only for 1st last distance
        #############################################################################
        
    
        
    
    


        
        #############################################################################
        ## only for nearest distance
        if self.radio_options_autoMultiPleAnalysis.GetSelection() == 3:
            
            siteN1 = []
            siteN2 = []
            distanceN1N2 = []
            
            for n in range(len(adjacent_gx2)):
                n1_shortest = 9999999999999999
                n2_shortest = 9999999999999999
                dr_shortest = 9999999999999999
                for m in range(len(adjacent_gx2)):
                    if n == m:
                        continue
                    else:
                        dxtemp = adjacent_gx2[m] - adjacent_gx2[n]
                        dytemp = adjacent_gy2[m] - adjacent_gy2[n]
                        drtemp = np.sqrt(dxtemp**2 + dytemp**2)
                        
                        if drtemp < dr_shortest:
                            dr_shortest = drtemp
                            n1_shortest = n
                            n2_shortest = m
                            
                newPair = 'y'
                for k in range(len(siteN1)):
                    if ((n1_shortest == siteN1[k]) and (n2_shortest == siteN2[k])) or  ((n1_shortest == siteN2[k]) and (n2_shortest == siteN1[k])):
                        newPair = 'n'
                        break
                    
                if newPair == 'y':
                    siteN1.append(n1_shortest)
                    siteN2.append(n2_shortest)
                    distanceN1N2.append(dr_shortest)
                    
                        
            
            
            
            
            print 'siteN1 ', siteN1
            print 'siteN2 ', siteN2
            print 'distanceN1N2 ', distanceN1N2
            
            distanceN1N2_nm = np.array(distanceN1N2) * float(self.M_nm_pixel.GetValue())

            
            return distanceN1N2_nm
            
        ## only for nearest distance
        #############################################################################
        
    
        
        
        
        

        self.NPairs_newRMS = 5
        NPairs_newRMS = self.NPairs_newRMS

        rms2_AfewPairs, rms2_Err_AfewPairs = RMS_distance_XY_Error(adjacent_gx2[:NPairs_newRMS], adjacent_gy2[:NPairs_newRMS], adjacent_gxy2Err[:NPairs_newRMS])
        print 'rms2_AfewPairs with err, rms2_Err_AfewPairs = ', rms2_AfewPairs, rms2_Err_AfewPairs


        rms2_AfewPairs_nm = rms2_AfewPairs * float(self.M_nm_pixel.GetValue())
        rms2_Err_AfewPairs_nm = rms2_Err_AfewPairs * float(self.M_nm_pixel.GetValue())

        AfewPairsRMS_text = 'First ' + str(NPairs_newRMS) + ' pairs RMS: ' + str(np.round(rms2_AfewPairs_nm,1)) + '+- ' + str(np.round(rms2_Err_AfewPairs_nm)) + ' nm'
















        
        
        
        #############################################################################################################################################
        ######### adjacent shrimp with sections ######################################################

                
        SectionBeginFN = self.SectionStartFN_temp # starting frame number for each section, eg) ...[0] value is the starting frome number of the first section
        SectionEndFN = self.SectionEndFN_temp
        TotalNsection = len(self.SectionEndFN_temp)

        
        print 'SectionBeginFN ', SectionBeginFN
        print 'SectionEndFN ', SectionEndFN
        
        
        

        
        Section_Intensity = [[] for _ in xrange(TotalNsection)]
        Section_FIONA_xpos = [[] for _ in xrange(TotalNsection)]
        Section_FIONA_ypos = [[] for _ in xrange(TotalNsection)]
        Section_FIONA_xposErr = [[] for _ in xrange(TotalNsection)]
        Section_FIONA_yposErr = [[] for _ in xrange(TotalNsection)]
                
        
        
        #print 'self.DataFileAnalysis_FIONA_fxy2_FN = ', self.DataFileAnalysis_FIONA_fxy2_FN
        #print 'self.DataFileAnalysis_fx2 = ', self.DataFileAnalysis_fx2
        #print '\n\nself.DataFileAnalysis_FIONA_fxy2_MaxIntensity', self.DataFileAnalysis_FIONA_fxy2_MaxIntensity
        
        for n in range(len(self.DataFileAnalysis_FIONA_fxy2_FN)):
            FNtemp = self.DataFileAnalysis_FIONA_fxy2_FN[n]
            #print 'FNtemp = ', FNtemp
            
            for k in range(TotalNsection):
                if SectionBeginFN[k] <= FNtemp <= SectionEndFN[k]:
                    Section_FIONA_xpos[k].append(self.DataFileAnalysis_fx2[n])
                    Section_FIONA_ypos[k].append(self.DataFileAnalysis_fy2[n])
                    
                    Section_FIONA_xposErr[k].append(self.DataFileAnalysis_fx2Err[n])
                    Section_FIONA_yposErr[k].append(self.DataFileAnalysis_fy2Err[n])
                    
                    Section_Intensity[k].append(self.DataFileAnalysis_FIONA_fxy2_MaxIntensity[n])
                    
                    
        

        #print 'Section_FIONA_xpos = ', Section_FIONA_xpos
        #print 'Section_FIONA_ypos = ', Section_FIONA_ypos
        
        #print 'Section_FIONA_xposErr = ', Section_FIONA_xposErr
        #print 'Section_FIONA_yposErr = ', Section_FIONA_yposErr
        
        
        #print '\nSection_Intensity = \n', Section_Intensity
            
            
            
            
        Section_Intensity_ave = []
        Section_FIONA_xpos_WeightedAve = []
        Section_FIONA_ypos_WeightedAve = []
        Section_FIONA_xpos_WeightedAveErr = []
        Section_FIONA_ypos_WeightedAveErr = []
        
        for n in range(TotalNsection):
            Section_Intensity_ave.append(np.mean(Section_Intensity[n]))
            xposWAveTemp, xposWAveErrTemp = Weighted_Ave_Err(Section_FIONA_xpos[n], Section_FIONA_xposErr[n])
            yposWAveTemp, yposWAveErrTemp = Weighted_Ave_Err(Section_FIONA_ypos[n], Section_FIONA_yposErr[n])
            
            
            Section_FIONA_xpos_WeightedAve.append(xposWAveTemp)
            Section_FIONA_xpos_WeightedAveErr.append(xposWAveErrTemp)
            
            Section_FIONA_ypos_WeightedAve.append(yposWAveTemp)            
            Section_FIONA_ypos_WeightedAveErr.append(yposWAveErrTemp)
            
                

        #Weighted_Ave_Err()


        
        Section_FIONA_xypos_WeightedAveErr = (np.array(Section_FIONA_xpos_WeightedAveErr) + np.array(Section_FIONA_ypos_WeightedAveErr))/2.0
        
        
        
        print 'Section_Intensity_ave : ', Section_Intensity_ave 
        print 'Section_FIONA_xpos_WeightedAve : ', Section_FIONA_xpos_WeightedAve 
        print 'Section_FIONA_ypos_WeightedAve : ', Section_FIONA_ypos_WeightedAve 
        print 'Section_FIONA_xpos_WeightedAveErr : ', Section_FIONA_xpos_WeightedAveErr 
        print 'Section_FIONA_ypos_WeightedAveErr : ', Section_FIONA_ypos_WeightedAveErr
        print 'Section_FIONA_xypos_WeightedAveErr : ', Section_FIONA_xypos_WeightedAveErr



        
        
        
        Section_AFshrimp_xpos = [[] for _ in xrange(TotalNsection)]
        Section_AFshrimp_ypos = [[] for _ in xrange(TotalNsection)]
        Section_AFshrimp_xposErr = [[] for _ in xrange(TotalNsection)]
        Section_AFshrimp_yposErr = [[] for _ in xrange(TotalNsection)]
        
        Section_AFshrimp_CorrespondingSectionN = [[] for _ in xrange(TotalNsection)]
        
        print '\n\n##############################################\n starting AF shrimp with sections'
        for n in range(TotalNsection):
            
            #print '\n ################## section n = ', n
            gx2postemp = []
            gy2postemp = []
            gx2posErrtemp = []
            gy2posErrtemp = []
            CorrSectiontemp = []
            
            for k in range(TotalNsection):
                #print '\n $$$$$$$$$$$$$ section k = ', k

                kxtemp = []
                kytemp = []
                kxErrtemp = []
                kyErrtemp = []
                
                if Section_Intensity_ave[n] > Section_Intensity_ave[k]:
                    for p in range(len(self.DataFileAnalysis_AFshrimp_gxy2_F1N)):
                        
                       
                        if (SectionBeginFN[n] <= self.DataFileAnalysis_AFshrimp_gxy2_F1N[p] <= SectionEndFN[n]) & (SectionBeginFN[k] <= self.DataFileAnalysis_AFshrimp_gxy2_F2N[p] <= SectionEndFN[k]):   
                            kxtemp.append(self.DataFileAnalysis_gx2[p])
                            kytemp.append(self.DataFileAnalysis_gy2[p])
                            kxErrtemp.append(self.DataFileAnalysis_gx2Err[p])
                            kyErrtemp.append(self.DataFileAnalysis_gy2Err[p])
                            

                        elif (SectionBeginFN[n] <= self.DataFileAnalysis_AFshrimp_gxy2_F2N[p] <= SectionEndFN[n]) & (SectionBeginFN[k] <= self.DataFileAnalysis_AFshrimp_gxy2_F1N[p] <= SectionEndFN[k]):   
                            kxtemp.append(self.DataFileAnalysis_gx2[p])
                            kytemp.append(self.DataFileAnalysis_gy2[p])
                            kxErrtemp.append(self.DataFileAnalysis_gx2Err[p])
                            kyErrtemp.append(self.DataFileAnalysis_gy2Err[p])
                            

                        else:
                            pass
                    gx2postemp.append(kxtemp)
                    gy2postemp.append(kytemp)
                    gx2posErrtemp.append(kxErrtemp)
                    gy2posErrtemp.append(kyErrtemp)
                    
                    
                    #if len(gx2postemp) != 0:
                    CorrSectiontemp.append(k)
                                            
                    
                else:
                    pass
                
            Section_AFshrimp_xpos[n] += (gx2postemp)
            Section_AFshrimp_ypos[n] += (gy2postemp)
            Section_AFshrimp_xposErr[n] += (gx2posErrtemp)
            Section_AFshrimp_yposErr[n] += (gy2posErrtemp)
            Section_AFshrimp_CorrespondingSectionN[n] += (CorrSectiontemp)
            

        
        Section_AFshrimp_index = []
        Section_AFshrimp_xposWeightedAve = []
        Section_AFshrimp_yposWeightedAve = []
        Section_AFshrimp_xposWeightedAveErr = []
        Section_AFshrimp_yposWeightedAveErr = []
        Section_AFshrimp_xyposWeightedAveErr_Mean = []
        
        for n in range(TotalNsection):
            
            for k in range(len(Section_AFshrimp_xpos[n])):
                
                xposWAveTemp, xposWAveErrTemp = Weighted_Ave_Err(Section_AFshrimp_xpos[n][k], Section_AFshrimp_xposErr[n][k])
                yposWAveTemp, yposWAveErrTemp = Weighted_Ave_Err(Section_AFshrimp_ypos[n][k], Section_AFshrimp_yposErr[n][k])
                
                
                Section_AFshrimp_index.append([n, Section_AFshrimp_CorrespondingSectionN[n][k] ])
                Section_AFshrimp_xposWeightedAve.append(xposWAveTemp)
                Section_AFshrimp_yposWeightedAve.append(yposWAveTemp)
                
                Section_AFshrimp_xposWeightedAveErr.append(xposWAveErrTemp)
                Section_AFshrimp_yposWeightedAveErr.append(yposWAveErrTemp)
                
                Section_AFshrimp_xyposWeightedAveErr_Mean.append( (xposWAveErrTemp + yposWAveErrTemp )/2.0 )
           
  
  
  
  
        #print 'Section_AFshrimp_index : ', Section_AFshrimp_index 
        #print 'Section_AFshrimp_xposWeightedAve : ', Section_AFshrimp_xposWeightedAve 
        #print 'Section_AFshrimp_yposWeightedAve : ', Section_AFshrimp_yposWeightedAve 
        #print 'Section_AFshrimp_xposWeightedAveErr : ', Section_AFshrimp_xposWeightedAveErr 
        #print 'Section_AFshrimp_yposWeightedAveErr : ', Section_AFshrimp_yposWeightedAveErr 
        #print 'Section_AFshrimp_xyposWeightedAveErr_Mean : ', Section_AFshrimp_xyposWeightedAveErr_Mean 
  
  
  
  
        adjacent_section_gx2 = []
        adjacent_section_gy2 = []
        adjacent_section_gxy2_err = []
        
        for n in range(TotalNsection-1):
            for k in range(len(Section_AFshrimp_index)):
                if Section_AFshrimp_index[k][0] == n and Section_AFshrimp_index[k][1] == n+1:
                    adjacent_section_gx2.append(Section_AFshrimp_xposWeightedAve[k])
                    adjacent_section_gy2.append(Section_AFshrimp_yposWeightedAve[k])
                    adjacent_section_gxy2_err.append(Section_AFshrimp_xyposWeightedAveErr_Mean[k])
                        
        adjacent_section_gx2.append(Section_FIONA_xpos_WeightedAve[-1])
        adjacent_section_gy2.append(Section_FIONA_ypos_WeightedAve[-1])
        adjacent_section_gxy2_err.append(Section_FIONA_xypos_WeightedAveErr[-1]) 
        
        
        
        #print 'adjacent_section_gx2 : ', adjacent_section_gx2
        #print 'adjacent_section_gy2 : ', adjacent_section_gy2
        #print 'adjacent_section_gxy2_err : ', adjacent_section_gxy2_err
        
        
        rms2_adjacent_section_gxy2, rms2_Err_adjacent_section_gxy2 = RMS_distance_XY_Error(adjacent_section_gx2, adjacent_section_gy2, adjacent_section_gxy2_err)
        #print 'Sectioned: adjacent_gxy2_distance_rms with err = ', rms2_adjacent_section_gxy2, rms2_Err_adjacent_section_gxy2


        rms2_nm_adjacent_section_gxy2 = rms2_adjacent_section_gxy2 * float(self.M_nm_pixel.GetValue())
        rms2_Err_nm_adjacent_section_gxy2 = rms2_Err_adjacent_section_gxy2 * float(self.M_nm_pixel.GetValue())
        
        SectionAdjacentRMS_text = 'Sectioned SHRImP rms = ' + str(np.around(rms2_adjacent_section_gxy2,3)) + ' +- ' + str(np.around(rms2_Err_adjacent_section_gxy2,3))+ ' px,  ' + str(np.around(rms2_nm_adjacent_section_gxy2,1)) + ' +- ' + str(np.around(rms2_Err_nm_adjacent_section_gxy2,1)) + ' nm'
        
        #print 'rms2_nm_adjacent_section_gxy2 : ', rms2_nm_adjacent_section_gxy2
        #print 'rms2_Err_nm_adjacent_section_gxy2 : ', rms2_Err_nm_adjacent_section_gxy2
        
        
        ######### adjacent shrimp with sections ###################################################### 
        #############################################################################################################################################

     

        





        
        #############################################################################################################################################
        ######### Exciton jump calculation by weighted Fiona ave in sections ###################################################### 

        if self.radio_options_autoMultiPleAnalysis.GetSelection() == 1 :
                
            print '\n\nCalculating Exciton jump distance'
            
            print 'Section_FIONA_xpos_WeightedAve ', Section_FIONA_xpos_WeightedAve
    
    
    
            #Section_Intensity_ave = []
            #Section_FIONA_xpos_WeightedAve = []
            #Section_FIONA_ypos_WeightedAve = []
            #Section_FIONA_xpos_WeightedAveErr = []
            #Section_FIONA_ypos_WeightedAveErr = []
    
    
            excitonJumpDistance = []
            for n in range(len(Section_FIONA_xpos_WeightedAve)-1):
                x1temp = Section_FIONA_xpos_WeightedAve[n]
                y1temp = Section_FIONA_ypos_WeightedAve[n]
                x2temp = Section_FIONA_xpos_WeightedAve[n+1]
                y2temp = Section_FIONA_ypos_WeightedAve[n+1]
                
                rtemp = np.sqrt((x1temp - x2temp)**2 + (y1temp - y2temp)**2)
                print 'rtemp ' , rtemp
                excitonJumpDistance.append(rtemp)
    
    
    
            return excitonJumpDistance
            

        ######## Exciton jump calculation by weighted Fiona ave ###################################################### '''
        #############################################################################################################################################

















        
        self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot.append(plt.figure('Adjacent_AFSHRIMP_DATA_ANALYSIS_Plot_' + str(self.AdjacentPlot_figCounter), figsize = (18, 9)))
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(150,30,1500, 750) 
        
        if self.cb_apply_dedrift.GetValue() == True:
            dedrift_text = 'de-drifted: applied'
        else:
            dedrift_text = 'de-drifted: none'
            
        
        adjacentFigText = '         x,y: '+self.MoleculeXYtext + '\n# of adjacent next frames: ' + str(self.MaxN_nextAdjacentFrames.GetValue()) + '\nShrimp Ecc < Each FIONA Ecc' + SecondEccConstraintText + '\nBleaching' + AdjacentBleachingTypeText + '\n' + dedrift_text
        self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot[self.AdjacentPlot_figCounter].text(0.02, 0.87, adjacentFigText, ha="left", va="bottom", size=10, weight = 'normal', color="red")
        
        
        x3 = np.linspace(0, 9, 10)
        y3 = np.linspace(0, 9, 10)
        x3, y3 = np.meshgrid(x3, y3)
        
        
        
        
        
        gs = gridspec.GridSpec(11, 25)
        gs.update(bottom = 0.04, top=0.94, left=0.05, right=0.95, wspace=5.0, hspace=3.0)
        plt.suptitle(str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly) + '\n' + self.DataFileAnalysisConditions
            , size=10)
 
       

        
        
        ### self.TotalIntensity_Trajectory
        NofOmission = int( self.AdjacentShrimpStartFrameN.GetValue() )
        initial_guess = (5000, 30, 2000)
        y = self.TotalIntensity_Trajectory[NofOmission:]
        x = np.arange(0, len(y)) 
        
        try:
            popt, pcov = opt.curve_fit(SingleExp, x, y, p0 = initial_guess)
            TotalIntensity_Trajectory_fit_temp = SingleExp(x, *popt)
            TotalIntensity_Trajectory_fit = np.append(np.arange(NofOmission),TotalIntensity_Trajectory_fit_temp)
        except:
            TotalIntensity_Trajectory_fit_temp = (np.arange(len(x)) * 0 + 1)*np.mean(self.TotalIntensity_Trajectory)
            TotalIntensity_Trajectory_fit = np.append(np.arange(NofOmission),TotalIntensity_Trajectory_fit_temp)
        

        

              
        
        axTotal_Intensity = plt.subplot(gs[1:3,0:9])
        plt.tick_params(labelsize=12)
        axTotal_Intensity.plot(TotalIntensity_Trajectory_fit, 'rx', label='Exp fit')
        axTotal_Intensity.plot(self.TotalIntensity_Trajectory)
        
        #ax1.plot(self.FIONA_centerIntensityAve)
        axTotal_Intensity.set_title('Total Intensity Trajectory', size=14)
        #plt.xlabel('Frame #')
        plt.legend(loc=1,prop={'size':10})
        
        
        
        axFIONA_ECC = plt.subplot(gs[1:3,9:15])
        plt.tick_params(labelsize=10)
        axFIONA_ECC.plot(self.FIONA_FrameN, self.FIONA_eccentricity, 'o-')
        #ax1.plot(self.FIONA_centerIntensityAve)
        axFIONA_ECC.set_title('FIONA ECC Trajectory', size=12)
        plt.xlabel('Frame #')
 

        
        
        
        axFIOINA2 = plt.subplot(gs[3:11,0:8], aspect='equal')
        plt.tick_params(labelsize=10)
        axFIOINA2.set_title('FIONA2 & Adjacent Shrimp,  # ' + str(NAdjacentPairs) + ',  rms = ' + str(np.around(rms2,3)) + ' +- ' + str(np.around(rms2_Err,3))+ 
        ',  ' + str(np.around(rms2_nm,1)) + ' +- ' + str(np.around(rms2_Err_nm,1)) + ' nm' + '\n' + SectionAdjacentRMS_text + '\n' + AfewPairsRMS_text, size =11)
        axFIOINA2.scatter(fx2, fy2, marker = '.',edgecolors='None', c = range(len(fx2)),  s=35, vmin=0, vmax= len(fx2))
        axFIOINA2.plot(adjacent_gx2, adjacent_gy2, 'b-')
        axFIOINA2.scatter(adjacent_gx2, adjacent_gy2, marker = 'H', edgecolors='None', c = range(len(adjacent_gx2)),  s=100, vmin=0, vmax= len(adjacent_gx2), label=str(adjacentFrameN1)+'\n'+str(adjacentFrameN2))
        
        axFIOINA2.plot(adjacent_section_gx2, adjacent_section_gy2, 'r-')
        axFIOINA2.scatter(adjacent_section_gx2, adjacent_section_gy2, marker = 's', edgecolors='None', c = range(len(adjacent_section_gx2)),  s=100, vmin=0, vmax= len(adjacent_section_gx2), label=str(SectionBeginFN)+'\n'+str(SectionEndFN))
        
        plt.xlim( np.mean(fx2) - 1, np.mean(fx2) + 1 )
        plt.ylim( np.mean(fy2) - 1, np.mean(fy2) + 1 )   
        plt.xlabel('Pixel Position')
        plt.legend(loc='upper right', prop={'size':10})
        
        
        
        

        axAFShrimp2 = plt.subplot(gs[3:11,8:16], aspect='equal')
        plt.tick_params(labelsize=10)
        axAFShrimp2.set_title('AF Shrimp0 & Adjacent Shrimp', size =14)
        axAFShrimp2.scatter(gx0, gy0, marker = '.', edgecolors='None', c = range(len(gx0)),  s=35, vmin=0, vmax= len(gx0))

        axAFShrimp2.plot(adjacent_gx2, adjacent_gy2, '-')
        axAFShrimp2.scatter(adjacent_gx2, adjacent_gy2, marker = 'H', edgecolors='None', c = range(len(adjacent_gx2)),  s=100, vmin=0, vmax= len(adjacent_gx2))

        
        plt.xlim( np.mean(fx2) - 1, np.mean(fx2) + 1 )
        plt.ylim( np.mean(fy2) - 1, np.mean(fy2) + 1 )     
        plt.xlabel('Pixel Position')
        
        



        gs2 = gridspec.GridSpec(11, 25)
        gs2.update(bottom = 0.05, top=0.95, left=0.102, right=0.924, wspace=5.0, hspace=2.0)
        
        axAFSIntensity = plt.subplot(gs2[1:3,16:24])
        plt.tick_params(labelsize=12)
        axAFSIntensity.plot(FrameIndex, FrameIntensity)
        #axAFSIntensity.xlim(0, AFSLastFrameN+1)
        
        #ax1.plot(self.FIONA_centerIntensityAve)
        axAFSIntensity.set_title('Intensity Trajectory for AF Shrimp', size=15)
        plt.xlabel('Frame #')
        plt.xlim(0, AFSLastFrameN+1)
      
      
      
        
        
        axFrameCorEcc = plt.subplot(gs[3:11,16:25], aspect='equal')
        plt.tick_params(labelsize=10)
        axFrameCorEcc.set_title("Frames' Eccentricity Correlation", size =16)
        plt.imshow(FrameCorMatrixEcc, interpolation = 'None', origin='bottom', vmin= 0.2, vmax=1.2)
        plt.colorbar()     
        plt.xlabel('Frame #')
        plt.ylabel('Frame #')
        

        
        plt.ion()
        plt.show()
        
        self.btn_EmissionSitesAnalysis.Show()
        self.AdjacentPlot_figCounter += 1


        # End of def AdjacentShrimpPlot(self, event)#
        #############################################
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        ####################################################
        ##  for publication figure  
        
        if self.cb_ShowPublicationFigures.GetValue():
        
            weighted_x = [386.275, 386.029, 386.207, 386.325]        
            weighted_y = [81.358, 81.422, 81.75, 81.703]
            
            
            
            
            
            font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 14}
    
            matplotlib.rc('font', **font)
            
            majortick_size = 8
            majortick_width = 3
            
            plt.rc('axes', linewidth = 3)
            
            #matplotlib.rcParams['lines.linewidth'] = 4
            
            matplotlib.rcParams['xtick.major.size'] = majortick_size
            matplotlib.rcParams['xtick.major.width'] = majortick_width
            matplotlib.rcParams['xtick.minor.size'] = 5
            matplotlib.rcParams['xtick.minor.width'] = 4
            
            matplotlib.rcParams['ytick.major.size'] = majortick_size
            matplotlib.rcParams['ytick.major.width'] = majortick_width
            matplotlib.rcParams['ytick.minor.size'] = 5
            matplotlib.rcParams['ytick.minor.width'] = 4
            
            
            matplotlib.rcParams.update({'font.size': 22, 'family' : 'normal', 'weight' : 'bold',  'size': 8})
       
    
    
    
    
            fig_publication = plt.figure('fxy2_adjacentXY_for_publication', figsize = (9.3,9))
            #fig_publication.patch.set_visible(False) # making transparent for outside of the main plot.
                       
            
            plt.tick_params(labelsize = 10 )
            #plt.title('FIONA2 & Adjacent Shrimp,  # ' + str(NAdjacentPairs) + ',  rms = ' + str(np.around(rms2,3)) + ' +- ' + str(np.around(rms2_Err,3)), size =12)
            plt.scatter(fx2, fy2, marker = '.',edgecolors='None', c = range(len(fx2)),  s=60, alpha=1, vmin=0, vmax= len(fx2), label = 'FIONA')
            #plt.plot(adjacent_gx2, adjacent_gy2, 'b-', linewidth  = 4)
    
            #plt.plot(weighted_x, weighted_y, 'r-', linewidth  = 4)
     
    
    
    
            AFshrimpErr_Cir = [[] for _ in xrange(len(adjacent_gxy2Err))]     
            
            markerColorAFS = ['DarkBlue', 'DarkTurquoise' , 'lime' , 'r']
            markerColorAFS2 = ['navy', 'DodgerBlue' , 'darkgreen' , 'darkred']
            
            '''
            for n in range(len(adjacent_gxy2Err)):
                
                AFshrimpErr_Cir[n] = plt.Circle((adjacent_gx2[n], adjacent_gy2[n]), radius = adjacent_gxy2Err[n], color = markerColorAFS[n], linewidth=3, fill=False)
                fig_publication.gca().add_artist(AFshrimpErr_Cir[n])
                
            '''
    
    
    
            plt.scatter(adjacent_gx2, adjacent_gy2, marker = 'H', edgecolors='k', linewidth='3', c = range(len(adjacent_gx2)),  s=500, vmin=0, vmax= len(adjacent_gx2), label='Adjacent SHRImP')
            #plt.scatter(adjacent_gx2, adjacent_gy2, marker = 'H', edgecolors='k', linewidth='3', c = markerColorAFS,  s=500, vmin=0, vmax= len(adjacent_gx2), label='Adjacent SHRImP')
       
           
            #plt.scatter(weighted_x, weighted_y, marker = 'H', edgecolors='k', c = range(len(adjacent_gx2)),  s=500, vmin=0, vmax= len(weighted_x), label='Sectioned SHRImP')
            #plt.scatter(weighted_x, weighted_y, marker = 'x', edgecolor='k',linewidth='5', s=500, c = range(len(adjacent_gx2)), vmin=0, vmax= len(weighted_x), label='Sectioned SHRImP')
            #plt.scatter(weighted_x, weighted_y, marker = 'x', linewidth='6', s=500, c = markerColorAFS2, vmin=0, vmax= len(weighted_x), label='Sectioned SHRImP')
            
            
            
            #plt.xlim( np.mean(fx2) - 0.48, np.mean(fx2) + 0.48 )
            #plt.ylim( np.mean(fy2) - 0.48, np.mean(fy2) + 0.48 )   

            plt.xlim( np.mean(fx2) - 0.8, np.mean(fx2) + 0.8 )
            plt.ylim( np.mean(fy2) - 0.8, np.mean(fy2) + 0.8 )   

            
            
            #plt.xticks(np.arange(232.185, 233.141, 0.2)) # for CF
            #plt.yticks(np.arange(281.426, 282.382, 0.2)) # for CF
            
            #plt.xticks(np.arange(386.0, 386.0 + 0.956, 0.2)) # for toluene
            #plt.yticks(np.arange(81.3, 81.3 + 0.956, 0.2)) # for toluene
            
            plt.xlabel('Pixel', fontweight="bold", fontsize = 29)
            plt.ylabel('Pixel', fontweight="bold", fontsize = 29)
            #plt.legend(loc='upper right', prop={'size':20, 'weight' : 'bold'})
            
            
            fig_publication.savefig('temp.png', transparent=True, dpi=400)
            
            
            
      
      
      
      
          
    
            
            plt.figure('ecc', figsize=(11,9.2))
            axFrameCorEcc_publication = plt.subplot(111, aspect='equal')
            plt.tick_params(labelsize=26)
            axFrameCorEcc_publication.set_title("Frames' Eccentricity Correlation", size =26)
            plt.imshow(FrameCorMatrixEcc, interpolation = 'None', origin='bottom', vmin= 0.19, vmax=1.21)
            plt.colorbar(ticks=[0.2, 0.4, 0.6, 0.8, 1.0, 1.2])     
            plt.xlabel('Frame #', fontweight="bold", fontsize = 36)
            plt.ylabel('Frame #', fontweight="bold", fontsize = 36)
            plt.xlim(0,97)
            plt.ylim(0,97)
            
            
            
            #plt.tight_layout()
            
    
            
        
        ##  for publication figure  temp use. '''
        #########################################
        
        
        
        
        
           
        
        
        return rms2, rms2_Err, fxy2_rms, fxy2_rms_Err, rms2_adjacent_section_gxy2, rms2_Err_adjacent_section_gxy2, rms2_AfewPairs, rms2_Err_AfewPairs
        
        












        

    def AdjacentEmissionSitesAnalysisWithSections(self, event):
        print '\n Starting Adjacent Emission Sites Analysis With Sections'
        
        self.text_Save_figures.SetLabel(' ')
        
        
        
        '''###
        SectionBeginFN = [] # starting frame number for each section, eg) ...[0] value is the starting frome number of the first section
        SectionEndFN = []   
        
        TotalNsection = 0
        
        print 'len(self.F_SectionBegin) = ', len(self.F_SectionBegin)
        
        for n in range(len(self.F_SectionBegin)):
            if (self.F_SectionBegin[n].GetValue() == '0') & (self.F_SectionEnd[n].GetValue() == '0' ):
                break
            else:
                TotalNsection += 1
                
        print 'TotalNsection = ',   TotalNsection 
        
        
        
        
        for n in range(TotalNsection):
            SectionBeginFN.append(int(self.F_SectionBegin[n].GetValue() ))
            SectionEndFN.append(int(self.F_SectionEnd[n].GetValue() ) )
            
        ### '''
        
        
        SectionBeginFN = [0, 7, 51, 55] # starting frame number for each section, eg) ...[0] value is the starting frome number of the first section
        SectionEndFN = [4, 49, 54, 97]   
        
        TotalNsection = 4

        
        print 'SectionBeginFN ', SectionBeginFN
        print 'SectionEndFN ', SectionEndFN
        
        
        

        
        Section_Intensity = [[] for _ in xrange(TotalNsection)]
        Section_FIONA_xpos = [[] for _ in xrange(TotalNsection)]
        Section_FIONA_ypos = [[] for _ in xrange(TotalNsection)]
        Section_FIONA_xposErr = [[] for _ in xrange(TotalNsection)]
        Section_FIONA_yposErr = [[] for _ in xrange(TotalNsection)]
                
        
        
        print 'self.DataFileAnalysis_FIONA_fxy2_FN = ', self.DataFileAnalysis_FIONA_fxy2_FN
        print 'self.DataFileAnalysis_fx2 = ', self.DataFileAnalysis_fx2
        print '\n\nself.DataFileAnalysis_FIONA_fxy2_MaxIntensity', self.DataFileAnalysis_FIONA_fxy2_MaxIntensity
        
        for n in range(len(self.DataFileAnalysis_FIONA_fxy2_FN)):
            FNtemp = self.DataFileAnalysis_FIONA_fxy2_FN[n]
            print 'FNtemp = ', FNtemp
            
            for k in range(TotalNsection):
                if SectionBeginFN[k] <= FNtemp <= SectionEndFN[k]:
                    Section_FIONA_xpos[k].append(self.DataFileAnalysis_fx2[n])
                    Section_FIONA_ypos[k].append(self.DataFileAnalysis_fy2[n])
                    
                    Section_FIONA_xposErr[k].append(self.DataFileAnalysis_fx2Err[n])
                    Section_FIONA_yposErr[k].append(self.DataFileAnalysis_fy2Err[n])
                    
                    Section_Intensity[k].append(self.DataFileAnalysis_FIONA_fxy2_MaxIntensity[n])
                    
                    
        

        print 'Section_FIONA_xpos = ', Section_FIONA_xpos
        print 'Section_FIONA_ypos = ', Section_FIONA_ypos
        
        print 'Section_FIONA_xposErr = ', Section_FIONA_xposErr
        print 'Section_FIONA_yposErr = ', Section_FIONA_yposErr
        
        
        print '\nSection_Intensity = \n', Section_Intensity
            
            
            
            
        Section_Intensity_ave = []
        Section_FIONA_xpos_WeightedAve = []
        Section_FIONA_ypos_WeightedAve = []
        Section_FIONA_xpos_WeightedAveErr = []
        Section_FIONA_ypos_WeightedAveErr = []
        
        for n in range(TotalNsection):
            Section_Intensity_ave.append(np.mean(Section_Intensity[n]))
            xposWAveTemp, xposWAveErrTemp = Weighted_Ave_Err(Section_FIONA_xpos[n], Section_FIONA_xposErr[n])
            yposWAveTemp, yposWAveErrTemp = Weighted_Ave_Err(Section_FIONA_ypos[n], Section_FIONA_yposErr[n])
            
            
            Section_FIONA_xpos_WeightedAve.append(xposWAveTemp)
            Section_FIONA_xpos_WeightedAveErr.append(xposWAveErrTemp)
            
            Section_FIONA_ypos_WeightedAve.append(yposWAveTemp)            
            Section_FIONA_ypos_WeightedAveErr.append(yposWAveErrTemp)
            
                

        #Weighted_Ave_Err()


        
        Section_FIONA_xypos_WeightedAveErr = (np.array(Section_FIONA_xpos_WeightedAveErr) + np.array(Section_FIONA_ypos_WeightedAveErr))/2.0
        
        
        
        
        
        
        Section_AFshrimp_xpos = [[] for _ in xrange(TotalNsection)]
        Section_AFshrimp_ypos = [[] for _ in xrange(TotalNsection)]
        Section_AFshrimp_xposErr = [[] for _ in xrange(TotalNsection)]
        Section_AFshrimp_yposErr = [[] for _ in xrange(TotalNsection)]
        
        Section_AFshrimp_CorrespondingSectionN = [[] for _ in xrange(TotalNsection)]
        
        print '\n\n##############################################\nstarting AF shrimp'
        for n in range(TotalNsection):
            
            print '\n\n ################## section n = ', n
            gx2postemp = []
            gy2postemp = []
            gx2posErrtemp = []
            gy2posErrtemp = []
            CorrSectiontemp = []
            
            for k in range(TotalNsection):
                print '\n $$$$$$$$$$$$$ section k = ', k

                kxtemp = []
                kytemp = []
                kxErrtemp = []
                kyErrtemp = []
                
                if Section_Intensity_ave[n] > Section_Intensity_ave[k]:
                    for p in range(len(self.DataFileAnalysis_AFshrimp_gxy2_F1N)):
                        
                       
                        if (SectionBeginFN[n] <= self.DataFileAnalysis_AFshrimp_gxy2_F1N[p] <= SectionEndFN[n]) & (SectionBeginFN[k] <= self.DataFileAnalysis_AFshrimp_gxy2_F2N[p] <= SectionEndFN[k]):   
                            kxtemp.append(self.DataFileAnalysis_gx2[p])
                            kytemp.append(self.DataFileAnalysis_gy2[p])
                            kxErrtemp.append(self.DataFileAnalysis_gx2Err[p])
                            kyErrtemp.append(self.DataFileAnalysis_gy2Err[p])
                            

                        elif (SectionBeginFN[n] <= self.DataFileAnalysis_AFshrimp_gxy2_F2N[p] <= SectionEndFN[n]) & (SectionBeginFN[k] <= self.DataFileAnalysis_AFshrimp_gxy2_F1N[p] <= SectionEndFN[k]):   
                            kxtemp.append(self.DataFileAnalysis_gx2[p])
                            kytemp.append(self.DataFileAnalysis_gy2[p])
                            kxErrtemp.append(self.DataFileAnalysis_gx2Err[p])
                            kyErrtemp.append(self.DataFileAnalysis_gy2Err[p])
                            

                        else:
                            pass
                    gx2postemp.append(kxtemp)
                    gy2postemp.append(kytemp)
                    gx2posErrtemp.append(kxErrtemp)
                    gy2posErrtemp.append(kyErrtemp)
                    
                    
                    #if len(gx2postemp) != 0:
                    CorrSectiontemp.append(k)
                                            
                    
                else:
                    pass
                
            Section_AFshrimp_xpos[n] += (gx2postemp)
            Section_AFshrimp_ypos[n] += (gy2postemp)
            Section_AFshrimp_xposErr[n] += (gx2posErrtemp)
            Section_AFshrimp_yposErr[n] += (gy2posErrtemp)
            Section_AFshrimp_CorrespondingSectionN[n] += (CorrSectiontemp)
            

        
        Section_AFshrimp_index = []
        Section_AFshrimp_xposWeightedAve = []
        Section_AFshrimp_yposWeightedAve = []
        Section_AFshrimp_xposWeightedAveErr = []
        Section_AFshrimp_yposWeightedAveErr = []
        Section_AFshrimp_xyposWeightedAveErr_Mean = []
        
        for n in range(TotalNsection):
            
            for k in range(len(Section_AFshrimp_xpos[n])):
                
                xposWAveTemp, xposWAveErrTemp = Weighted_Ave_Err(Section_AFshrimp_xpos[n][k], Section_AFshrimp_xposErr[n][k])
                yposWAveTemp, yposWAveErrTemp = Weighted_Ave_Err(Section_AFshrimp_ypos[n][k], Section_AFshrimp_yposErr[n][k])
                
                
                Section_AFshrimp_index.append([n, Section_AFshrimp_CorrespondingSectionN[n][k] ])
                Section_AFshrimp_xposWeightedAve.append(xposWAveTemp)
                Section_AFshrimp_yposWeightedAve.append(yposWAveTemp)
                
                Section_AFshrimp_xposWeightedAveErr.append(xposWAveErrTemp)
                Section_AFshrimp_yposWeightedAveErr.append(yposWAveErrTemp)
                
                Section_AFshrimp_xyposWeightedAveErr_Mean.append( (xposWAveErrTemp + yposWAveErrTemp )/2.0 )
           
  
  
  
  
  
  
  
  
  
  
  
  
  
  
        
        self.fig_EmissionSitesAnalysis.append(plt.figure('Emission_Sites_ANALYSIS_Plot_' + str(self.fig_EmissionSitesAnalysisFinalCounter), figsize = (18, 9)))
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(50,30,1600, 800)        
        plt.suptitle(str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly) + '\n' + self.DataFileAnalysisConditions
            , size=11)
        self.fig_EmissionSitesAnalysis[self.fig_EmissionSitesAnalysisFinalCounter].text(0.05, 0.94, '         x,y: '+self.MoleculeXYtext, ha="left", va="bottom", size=13, weight = 'normal', color="red")   

        plt.subplots_adjust(top = 0.90, bottom = 0.03, left = 0.03, right = 0.98, wspace = 0.2, hspace = 0.4)
        markerColor = ['b' , 'g' , 'c', 'm', 'y', 'k', 'b', 'g', 'c', 'm', 'y', 'k']
        
        for n in range(TotalNsection):
                
            plt.subplot(2,5,n+1, aspect='equal')
            plt.tick_params(labelsize=10)
            plt.title('Section '+str(n+1) + ',   F# ' +str(SectionBeginFN[n])+':'+str(SectionEndFN[n])       , size =12)
            
            for k in range(len(Section_AFshrimp_xpos[n])):
                if len(Section_AFshrimp_xpos[n][k]) != 0:
                    plt.scatter(Section_AFshrimp_xpos[n][k], Section_AFshrimp_ypos[n][k], marker = 'x', edgecolors=markerColor[k], color = markerColor[k], label = str(Section_AFshrimp_CorrespondingSectionN[n][k] +1))
                
                                
            #plt.scatter(gx2, gy2, marker = '.', edgecolors='None', c = range(len(gx2)),  s=35, vmin=0, vmax= len(gx2)) 
            print 'len(Section_FIONA_xpos[n]) = ', len(Section_FIONA_xpos[n])
            
            if len(Section_FIONA_xpos[n]) != 0 :
                plt.scatter(Section_FIONA_xpos[n], Section_FIONA_ypos[n], marker = '.',s= 100, edgecolors='r', color = 'r', label = 'FIONA')
                plt.scatter(Section_FIONA_xpos[n], Section_FIONA_ypos[n], marker = '.', edgecolors='None', c = range(len(Section_FIONA_xpos[n])),  s=35, vmin=0, vmax= len(Section_FIONA_xpos[n]), label = 'FIONA')
            
                plt.legend(loc='upper right', prop={'size':10})
                plt.xlim( (np.around(Section_FIONA_xpos_WeightedAve[n]) - 1, np.around(Section_FIONA_xpos_WeightedAve[n]) + 1 ) )
                plt.ylim( (np.around(Section_FIONA_ypos_WeightedAve[n]) - 1, np.around(Section_FIONA_ypos_WeightedAve[n]) + 1 ) )   
                plt.xlabel('Pixel Position')
                







            
        self.fig_EmissionSitesAnalysisFinal.append(plt.figure('Emission_Sites_ANALYSIS_Final_'+str(self.fig_EmissionSitesAnalysisFinalCounter), figsize = (18, 6)) )
        self.fig_EmissionSitesAnalysisFinal[self.fig_EmissionSitesAnalysisFinalCounter].text(0.05, 0.92, '         x,y: '+self.MoleculeXYtext, ha="left", va="bottom", size=13, weight = 'normal', color="red")   
        
        
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(50,500,1600, 533)
        #mngr.window.move(50, 500)
        

        
        
        plt.suptitle(str(self.filedirectoryOnly) + ' \n ' + str(self.filenameOnly) + '\n' + self.DataFileAnalysisConditions, size=11)
        plt.subplots_adjust(top = 0.84, bottom = 0.05, left = 0.04, right = 0.98, wspace = 0.2, hspace = 0.4)  
        
        markerColorAFS = ['r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k','r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k','r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k','r', 'b' , 'g' , 'c', 'm', 'y', 'k','r', 'b', 'g', 'c', 'm', 'y', 'k']

        plt.subplot(1,3,1, aspect='equal')
        plt.title("FIONA")
        print 'Section_FIONA_xpos_WeightedAve = ', Section_FIONA_xpos_WeightedAve

        FIONA_Cir = [[] for _ in xrange(TotalNsection)]
        for n in range(len(Section_FIONA_xpos_WeightedAve)):
            print 'n = ', n
            if str(Section_FIONA_xpos_WeightedAve[n]) != 'nan':
                print 'Section_FIONA_xpos_WeightedAve[n] ', Section_FIONA_xpos_WeightedAve[n], type(Section_FIONA_xpos_WeightedAve[n])
                print 'run n = ', n
                plt.scatter(Section_FIONA_xpos_WeightedAve[n], Section_FIONA_ypos_WeightedAve[n], marker = '.',s= 100, edgecolors = markerColorAFS[n], c = markerColorAFS[n]
                , label = 'S' + str(n+1) + ' ('+str(np.around(Section_FIONA_xpos_WeightedAve[n],3))+', ' + str(np.around(Section_FIONA_ypos_WeightedAve[n],3)) + ')  E ' + str(np.around(Section_FIONA_xypos_WeightedAveErr[n],3))  )
            else:
                print 'pass'
                pass
            FIONA_Cir[n] = plt.Circle((Section_FIONA_xpos_WeightedAve[n], Section_FIONA_ypos_WeightedAve[n]), radius = Section_FIONA_xypos_WeightedAveErr[n], color = markerColorAFS[n], fill=False)
            self.fig_EmissionSitesAnalysisFinal[self.fig_EmissionSitesAnalysisFinalCounter].gca().add_artist(FIONA_Cir[n])
            

        plt.legend(loc='upper right', prop={'size':12})
        plt.xlim( (np.around(Section_FIONA_xpos_WeightedAve[0]) - 1, np.around(Section_FIONA_xpos_WeightedAve[0]) + 1 ) )
        plt.ylim( (np.around(Section_FIONA_ypos_WeightedAve[0]) - 1, np.around(Section_FIONA_ypos_WeightedAve[0]) + 1 ) )     
            
        
        
        self.Section_FIONA_xpos_WeightedAve = Section_FIONA_xpos_WeightedAve
        self.Section_FIONA_ypos_WeightedAve = Section_FIONA_ypos_WeightedAve
        self.Section_FIONA_xypos_WeightedAveErr = Section_FIONA_xypos_WeightedAveErr
        self.TotalNsection = TotalNsection

        

        plt.subplot(1,3,2, aspect='equal')
        plt.title("AF Shrimp")

        N_SectionAfshrimpIndex = len(Section_AFshrimp_index)
        AFshrimp_Cir = [[] for _ in xrange(N_SectionAfshrimpIndex)]
        
        
        print 'N_SectionAfshrimpIndex ', N_SectionAfshrimpIndex
        print 'Section_AFshrimp_xposWeightedAve = ' , Section_AFshrimp_xposWeightedAve
        
        for n in range(N_SectionAfshrimpIndex):
            print ' n = ', n

            if str(Section_AFshrimp_xposWeightedAve[n]) != 'nan':
                plt.scatter(Section_AFshrimp_xposWeightedAve[n], Section_AFshrimp_yposWeightedAve[n], marker = 'x', c = markerColorAFS[n]
                , label = 'S' + str(Section_AFshrimp_index[n][0]+1) + '-' + str(Section_AFshrimp_index[n][1]+1) + ' ('+str(np.around(Section_AFshrimp_xposWeightedAve[n],3))+', ' + str(np.around(Section_AFshrimp_yposWeightedAve[n],3)) + ')  E ' + str(np.around(Section_AFshrimp_xyposWeightedAveErr_Mean[n],3)) )
            else:
                print 'pass'
                pass
            
            AFshrimp_Cir[n] = plt.Circle((Section_AFshrimp_xposWeightedAve[n], Section_AFshrimp_yposWeightedAve[n]), radius = Section_AFshrimp_xyposWeightedAveErr_Mean[n], color = markerColorAFS[n], fill=False)
            self.fig_EmissionSitesAnalysisFinal[self.fig_EmissionSitesAnalysisFinalCounter].gca().add_artist(AFshrimp_Cir[n])
            

        plt.legend(loc='upper right', prop={'size':12})
        plt.xlim( (np.around(Section_FIONA_xpos_WeightedAve[0]) - 1, np.around(Section_FIONA_xpos_WeightedAve[0]) + 1 ) )
        plt.ylim( (np.around(Section_FIONA_ypos_WeightedAve[0]) - 1, np.around(Section_FIONA_ypos_WeightedAve[0]) + 1 ) )     
            
        
        
        

        plt.subplot(1,3,3, aspect='equal')
        plt.title("FIONA + AF Shrimp")

        N_SectionAfshrimpIndex = len(Section_AFshrimp_index)
        AFshrimp_Cir = [[] for _ in xrange(N_SectionAfshrimpIndex)]

        FIONA_Cir = [[] for _ in xrange(TotalNsection)]
        for n in range(len(Section_FIONA_xpos_WeightedAve)):
            print 'n = ', n
            if str(Section_FIONA_xpos_WeightedAve[n]) != 'nan':
                print 'Section_FIONA_xpos_WeightedAve[n] ', Section_FIONA_xpos_WeightedAve[n], type(Section_FIONA_xpos_WeightedAve[n])
                print 'run n = ', n
                plt.scatter(Section_FIONA_xpos_WeightedAve[n], Section_FIONA_ypos_WeightedAve[n], marker = '.',s= 100, edgecolors = markerColorAFS[n], c = markerColorAFS[n], label = 'FIONA ' + str(n+1))
            else:
                print 'pass'
                pass

            
            
        for n in range(N_SectionAfshrimpIndex):

            if str(Section_AFshrimp_xposWeightedAve[n]) != 'nan':
                plt.scatter(Section_AFshrimp_xposWeightedAve[n], Section_AFshrimp_yposWeightedAve[n], marker = 'x', c = markerColorAFS[n])
            else:
                print 'pass'
                pass
            

            
        plt.legend(loc='upper right', prop={'size':12})
        plt.xlim( (np.around(Section_FIONA_xpos_WeightedAve[0]) - 1, np.around(Section_FIONA_xpos_WeightedAve[0]) + 1 ) )
        plt.ylim( (np.around(Section_FIONA_ypos_WeightedAve[0]) - 1, np.around(Section_FIONA_ypos_WeightedAve[0]) + 1 ) )     
        
        
        self.fig_EmissionSitesAnalysisFinalCounter += 1
        
        print 'self.fig_EmissionSitesAnalysisFinalCounter = ', self.fig_EmissionSitesAnalysisFinalCounter
        
        self.btn_Save_figures.Show()   
        #self.btn_EmissionSitesAnalysis.Hide()
        self.btn_PlotFinalEmissionSites.Show()
        
        
        
        
        self.Section_FIONA_xpos = Section_FIONA_xpos
        self.Section_FIONA_ypos = Section_FIONA_ypos

        #Section_AFshrimp_index = []
        self.Section_AFshrimp_xpos = Section_AFshrimp_xpos
        self.Section_AFshrimp_ypos = Section_AFshrimp_ypos
        self.Section_AFshrimp_CorrespondingSectionN = Section_AFshrimp_CorrespondingSectionN
        






        
    def SaveFigures(self, event):
        prefix = self.FilePrefix.GetValue()
        figFileName = prefix + self.filenameOnly[:35]
        todayDate = time.strftime("%Y%m%d_%Hh%Mm")
        print "\nToday's Date: ", time.strftime("%Y-%m-%d_%Hh%Mm")
        
        
        MoleculeXY = '_(' + str(self.MoleculeXpos)+','+str(self.MoleculeYpos) + ')'
                
        
         
        save_fig_error = 0
        
        
        
        
        
        print 'len(fig_FIONA_AFSHRIMP_DATA_ANALYSIS_Plot) = ', len(self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_Plot)
        
        
        try:
            for n in range(len(self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_Plot)):
                self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_Plot[n].savefig(figFileName + '_' + todayDate +MoleculeXY+ '_1_FIONA_AFshrimp' + str(n)+ '.png')
                plt.close(self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_Plot[n])
            print 'fig_FIONA_AFSHRIMP_DATA_ANALYSIS_Plot figures are saved'
        except:
            print 'fig_FIONA_AFSHRIMP_DATA_ANALYSIS_Plot figures are NOT saved'      
            save_fig_error += 1
            



        

        print 'len()self.fig_EmissionSitesAnalysis ', len(self.fig_EmissionSitesAnalysis)
                

        try:
            for n in range(len(self.fig_EmissionSitesAnalysis)):
                self.fig_EmissionSitesAnalysis[n].savefig(figFileName + '_' + todayDate +MoleculeXY+ '_2_EmissionSites1_' + str(n)+ '.png', dpi = 95)
                plt.close(self.fig_EmissionSitesAnalysis[n])
            print 'fig_Emission Sites Analysis figures are saved'
        except:
            print 'fig_Emission Sites Analysis figures are NOT saved'      
            save_fig_error += 1
            


         
        
        print 'len()self.fig_EmissionSitesAnalysisFinal ', len(self.fig_EmissionSitesAnalysisFinal)
        
        
        try:        
            for n in range(len(self.fig_EmissionSitesAnalysisFinal)):    
                self.fig_EmissionSitesAnalysisFinal[n].savefig(figFileName + '_' + todayDate +MoleculeXY+ '_3_EmissionSites2_' + str(n) + '.png', dpi = 80)
                plt.close(self.fig_EmissionSitesAnalysisFinal[n])
            print 'fig_EmissionSitesAnalysisFinal figures are saved'
            
        except:
            print 'fig_EmissionSitesAnalysisFinal figures are NOT saved'  
            save_fig_error += 1
            
            


        try:        
            
            self.fig_PlotFinalEmissionSites.savefig(figFileName + '_' + todayDate +MoleculeXY+ '_4_EmissionSites_Final' + '.png', dpi = 80)
            plt.close(self.fig_PlotFinalEmissionSites)
                
            print 'fig_PlotFinalEmissionSites figures are saved'
            
        except:
            print 'fig_PlotFinalEmissionSites figures are NOT saved'  
            save_fig_error += 1
            
            

        
        
        if save_fig_error == 0 :
            save_fig_text = 'All figures are saved'
            print '\n', save_fig_text, '\n'
        else:
            save_fig_text = 'NOT All figures are saved'
            print '\n NOT All figures are saved\n # error: ', save_fig_error
            
        
         
     
     
        self.text_Save_figures.SetLabel(save_fig_text)
        
        #self.btn_Save_figures.Hide()
        
        
            












        
    def SaveFigures_LastAdjacentShrimp(self, event): # saving only the last adjacent figure.
        prefix = self.FilePrefix.GetValue()
        figFileName = prefix + 'Adjacent_' + self.filenameOnly[:40]
        todayDate = time.strftime("%Y%m%d_%Hh%Mm")
        print "\nToday's Date: ", time.strftime("%Y-%m-%d_%Hh%Mm")
        
        
        MoleculeXY = '(' + str(self.MoleculeXpos)+','+str(self.MoleculeYpos) + ')_'
                
        
         
        save_fig_error = 0
        
        
        

        '''
        try:        
            for n in range(len(self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot)):    
                self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot[n].savefig(figFileName + '_' + todayDate +MoleculeXY+ '_3_EmissionSites2_' + str(n) + '.png', dpi = 80)
                #plt.close(self.fig_EmissionSitesAnalysisFinal[n])
            print 'fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot figures are saved'
            
        except:
            print 'fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot figures are NOT saved'  
            save_fig_error += 1
            
        '''    


        try:        
            adjacentFigN = len(self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot)
            self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot[adjacentFigN-1].savefig(figFileName + '_' + MoleculeXY + todayDate  + '.png', dpi = 80)
            plt.close(self.fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot[adjacentFigN-1])
                
            print 'fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot figures are saved'
            
        except:
            print 'fig_FIONA_AFSHRIMP_DATA_ANALYSIS_AdjacentPlot figures are NOT saved'  
            save_fig_error += 1
            
            

        
        
        if save_fig_error == 0 :
            save_fig_text = 'Adjacent Plots are saved'
            print '\n', save_fig_text, '\n'
        else:
            save_fig_text = 'Error'
            print '\n NOT All figures are saved\n # error: ', save_fig_error
            
        
        self.text_LastAdjacentShrimp_Figsaved.SetLabel(save_fig_text)
        
        plt.close('fxy2_distance')
     
     
        
        
        #self.btn_Save_figures.Hide()
        
        
            











'''
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
'''








def IntensityTrajectoryImagePlotsFunc(loadedFrames, filename, xpos, ypos, imageFrameN1, imageFrameN2):

    print '\n Running  IntensityTrajectoryImagePlotsFunc' 
    print 'filename = ', filename
    frames = loadedFrames
    
    frameStart = 0
    frameEnd = len(frames) - 1
    
    print 'len: frames ', len(frames)    
    
    
    
    ########################################################
    # mismatching correction between python numpy and ImageJ

    #N_xPixel = frames.frame_shape[0]
    #N_yPixel = frames.frame_shape[1]
    
    x_ImageJ, y_ImageJ = np.array(xpos), np.array(ypos)
    
    #xc = y_ImageJ - 1
    #yc = N_xPixel - x_ImageJ + 1
    
    xc = x_ImageJ - 1 # for trackpy 3.0
    yc = y_ImageJ # for trackpy 3.0
            
    ########################################################




    
    # calculating each frame for all features => 
    
    featureImages = [[] for _ in xrange(len(xc))] 
    featureImagesFrameN = [[] for _ in xrange(len(xc))] 
    
    
    centerIntensityAveTrajectory = [[] for _ in xrange(len(xc))] 
    centerIntensityMaxTrajectory = [[] for _ in xrange(len(xc))]
    centerIntensityMax5AveTrajectory = [[] for _ in xrange(len(xc))]
    centerIntensityIntegratedTrajectory = [[] for _ in xrange(len(xc))] 
    
    
    for k in range(frameStart,frameEnd+1):
        print '\n frame # ', k
        frametemp = frames[k]
        for n in range(len(xc)):
            
            frameLargeTemp = frametemp[yc[n]-7:yc[n]+8, xc[n]-7:xc[n]+8]
            frameCenterTemp = frametemp[yc[n]-3:yc[n]+4, xc[n]-3:xc[n]+4]
            #frameCenterTemp2 = frametemp[yc[n]-4:yc[n]+5, xc[n]-4:xc[n]+5]
            frameCenterTemp3 = frametemp[yc[n]-6:yc[n]+7, xc[n]-6:xc[n]+7]
            
            
            '''
            frameCenterTemp4 = frametemp[yc[n]-1:yc[n]+2, xc[n]-1:xc[n]+2]
            frameCenterTemp5 = frametemp[yc[n]-2:yc[n]+3, xc[n]-2:xc[n]+3]
            BGave3 = ( np.sum(frameCenterTemp5.ravel()) - np.sum(frameCenterTemp4.ravel()) ) / 16
            IntegratedIntesityTemp3 = np.sum(frameCenterTemp4.ravel()) - BGave3*9
            '''
            
            
            #print 'frameCenterTemp4 \n', frameCenterTemp4
            #print 'frameCenterTemp5 \n', frameCenterTemp5
            
            
            #print 'np.sum(frameCenterTemp4.ravel() ', np.sum(frameCenterTemp4.ravel() )
            #print 'np.sum(frameCenterTemp5.ravel() ', np.sum(frameCenterTemp5.ravel() )
            
            
            
            FiveMaxAveTemp = int(np.mean( np.sort(frameCenterTemp.ravel())[-5:]) )
            
            #BGave = ( np.sum(frameCenterTemp2.ravel()) - np.sum(frameCenterTemp.ravel()) ) / 32
            
            BGave2 = ( np.sum(frameLargeTemp.ravel()) - np.sum(frameCenterTemp3.ravel()) ) / 56
            
            #print 'BGave = ', BGave, BGave2, BGave3
            
            
            #IntegratedIntesityTemp = np.sum(frameCenterTemp2.ravel()) - BGave*81
            #print 'IntegratedIntesityTemp ', IntegratedIntesityTemp
            
            IntegratedIntesityTemp2 = np.sum(frameCenterTemp3.ravel()) - BGave2*169
            
            #print 'IntegratedIntesityTemp = ', IntegratedIntesityTemp, IntegratedIntesityTemp2, IntegratedIntesityTemp3
            
            
            centerIntensityAveTrajectory[n].append(np.mean(frameCenterTemp))
            centerIntensityMaxTrajectory[n].append(np.amax(frameCenterTemp)) 
            centerIntensityMax5AveTrajectory[n].append(FiveMaxAveTemp)
            centerIntensityIntegratedTrajectory[n].append(IntegratedIntesityTemp2)
            
            
            
            if k in [imageFrameN1, imageFrameN2]:
                featureImages[n].append(frameLargeTemp)
                
                
            





    print '\n Endding  IntensityTrajectoryImagePlotsFunc\n\n'       
    
    return centerIntensityIntegratedTrajectory, centerIntensityAveTrajectory, centerIntensityMaxTrajectory, centerIntensityMax5AveTrajectory, featureImages, featureImagesFrameN
    



  
class MainTrackingIntensityOnlyPlots(wx.Frame):
    #ImageFilePath = ''
 
    #----------------------------------------------------------------------
    def __init__(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, "SM2 Tracking Intensity Trajectory", size=(1000, 1000), pos=(0,0))
                                  
        self.MainPanel = wx.Panel(self, wx.ID_ANY)
                          
        
        self.MainPanel.Bind(wx.EVT_PAINT, self.on_paint)
        
                
        
        

        self.btn_SwitchToAFShrimp = wx.Button(self.MainPanel, pos=(300,10), label="Switch To afShrimp")
        self.btn_SwitchToAFShrimp.Bind(wx.EVT_BUTTON, self.SwitchToAFShrimp)
        self.btn_SwitchToDataFileAnalysis = wx.Button(self.MainPanel, pos=(500,10), label="Switch To Data File Analysis")
        self.btn_SwitchToDataFileAnalysis.Bind(wx.EVT_BUTTON, self.SwitchToDataFileAnalysis)
        #self.btn_SwitchToIntensityTrajectoryPlots = wx.Button(self.MainPanel, pos=(700,10), label="Switch To Intensity Trajectory Plots with Images")
        #self.btn_SwitchToIntensityTrajectoryPlots.Bind(wx.EVT_BUTTON, self.SwitchToIntensityTrajectoryPlots)



        
        wx.StaticText(self.MainPanel, -1, 'Grayscale Tiff Images', pos=(50, 20))
        
        self.btn_OpenImage = wx.Button(self.MainPanel, pos=(50,50), label="Open Image File")
        self.ImageFilePath = self.btn_OpenImage.Bind(wx.EVT_BUTTON, self.onOpenImageFile)
        self.btn_OpenImage_Text = wx.StaticText(self.MainPanel, -1, str(self.ImageFilePath), pos=(170, 53))
        
        
        
        
        wx.StaticText(self.MainPanel, -1, "Choose a Frame #:", pos=(170, 103))
        #self.Nframe = wx.TextCtrl(self.MainPanel, -1, "0", pos=(280, 100), size=(40,-1))
        self.Nframe = wx.SpinCtrl(self.MainPanel, -1, pos=(280, 100), size=(60,-1))
        self.Nframe.SetValue(0)
                
        
        self.btn_FindingMolecules = wx.Button(self.MainPanel, pos=(50,100), label="Finding Molecules")
        self.btn_FindingMolecules.Hide()
        self.btn_FindingMolecules.Bind(wx.EVT_BUTTON, self.FindingFeatureFunc)
        
        wx.StaticText(self.MainPanel, -1, "Average # frames:", pos=(370, 103))
        self.NofFramesForAveraging = wx.SpinCtrl(self.MainPanel, -1, pos=(480, 100), size=(50,-1))
        self.NofFramesForAveraging.SetValue(5)
        
        wx.StaticText(self.MainPanel, -1, "Feature Size:", pos=(580, 103))
        self.FeatureSize = wx.SpinCtrl(self.MainPanel, -1,  pos=(660, 100), size=(50,-1))
        self.FeatureSize.SetValue(7)
        
        wx.StaticText(self.MainPanel, -1, "Min Intensity:", pos=(400, 133))
        self.MinIntensity = wx.TextCtrl(self.MainPanel, -1, "40000", pos=(480, 130), size=(70,-1))
        
        wx.StaticText(self.MainPanel, -1, "Number of Features:", pos=(170, 133))
        self.NofFeatures = wx.StaticText(self.MainPanel, -1, "None", pos=(300, 133), style=5)
        
        
        self.btn_IntensityTrajectoryPlots = wx.Button(self.MainPanel, pos=(50,180), label="Image and Intensity Trajectory Plots")
        self.btn_IntensityTrajectoryPlots.Hide()
        self.btn_IntensityTrajectoryPlots.Bind(wx.EVT_BUTTON, self.IntensityTrajectoryPlots)
        
        wx.StaticText(self.MainPanel, -1, "Two image frame numbers: ", pos=(50, 213))
        self.imageFrameN1 = wx.TextCtrl(self.MainPanel, -1, "0", pos=(200, 210), size=(40,-1))
        self.imageFrameN2 = wx.TextCtrl(self.MainPanel, -1, "49", pos=(250, 210), size=(40,-1))
                




        self.btn_FIONA_Plots = wx.Button(self.MainPanel, pos=(450,180), label="FIONA Plots for A Few Bright Features")
        self.btn_FIONA_Plots.Hide()
        self.btn_FIONA_Plots.Bind(wx.EVT_BUTTON, self.FIONA_Plots)
        
        wx.StaticText(self.MainPanel, -1, "# of Frames: ", pos=(450, 213))
        self.NfastFIONAframe = wx.TextCtrl(self.MainPanel, -1, "100", pos=(525, 210), size=(40,-1))
                
        wx.StaticText(self.MainPanel, -1, "# of Features: \n (not yet) ", pos=(450, 233))
        self.N_FIONA = wx.TextCtrl(self.MainPanel, -1, "14", pos=(525, 230), size=(40,-1))
                
                
        
  
   
        wx.StaticText(self.MainPanel, -1, "Data and Figure Files Prefix", pos=(50, 322))
        self.FilePrefix = wx.TextCtrl(self.MainPanel, -1, "", pos=(200, 320))
         

        self.btn_SaveIntensityTrajectoryPlots = wx.Button(self.MainPanel, pos=(50,350), label="Save Intensity Trajectory & Image Plots , data files")
        self.btn_SaveIntensityTrajectoryPlots.Hide()
        self.btn_SaveIntensityTrajectoryPlots.Bind(wx.EVT_BUTTON, self.SaveIntensityTrajectoryPlotsData)
        self.text_SaveIntensityTrajectoryPlots = wx.StaticText(self.MainPanel, -1, ' ', pos=(50, 380)) 
        
        
        self.btn_SaveFIONA_Plots = wx.Button(self.MainPanel, pos=(450,350), label="Save FIONA Plots for Bright Features")
        self.btn_SaveFIONA_Plots.Hide()
        self.btn_SaveFIONA_Plots.Bind(wx.EVT_BUTTON, self.SaveFIONA_Plots)
        self.text_FIONA_Figsaved = wx.StaticText(self.MainPanel, -1, ' ', pos=(450, 380))
  

        self.cb_save_data_files_yn = wx.CheckBox(self.MainPanel, -1, 'save data files', (350, 323))
        self.cb_save_data_files_yn.SetValue(False)






        title2 = wx.StaticText(self.MainPanel, -1, "Plotting intensity trajectories from Image (x,y)", pos=(50, 505))
        font = wx.Font(12, wx.DECORATIVE, wx.NORMAL, wx.NORMAL)
        title2.SetFont(font)     

        
        
        
        
        self.N_ImageJ_features = wx.SpinCtrl(self.MainPanel, -1, pos=(50, 540), size=(50,-1))
        self.N_ImageJ_features.SetRange(1, 50)  
        self.N_ImageJ_features.SetValue(1)
        wx.StaticText(self.MainPanel, -1, ": Number of ImageJ (x,y) features", pos=(110, 543))
        
        

        self.btn_ImageJ_xy_inputs = wx.Button(self.MainPanel, pos=(330,540), label="Create ImageJ x,y inputs")
        self.btn_ImageJ_xy_inputs.Bind(wx.EVT_BUTTON, self.Create_ImageJ_xy_inputs)
        
        self.Create_ImageJ_xy_inputs('')
        
        
        

        
        
        wx.StaticText(self.MainPanel, -1, 'Enter (x, y) from ImageJ', pos=(50, 575))
        



        

        self.btn_PlotIntensity_from_ImageJ_xy_inputs = wx.Button(self.MainPanel, pos=(50,880), label="Plot intesnity trajectories from ImageJ x,y inputs")
        self.btn_PlotIntensity_from_ImageJ_xy_inputs.Bind(wx.EVT_BUTTON, self.PlotIntensity_from_ImageJ_xy_inputs)
        self.btn_PlotIntensity_from_ImageJ_xy_inputs.Hide()
        
        







        
    #----------------------------------------------------------------------#


    def on_paint(self, event):
        dc = wx.PaintDC(event.GetEventObject())
        dc.Clear()
        dc.SetPen(wx.Pen("BLACK", 4))
        dc.DrawLine(0, 500, 1000, 500)  
            
        

    def SwitchToAFShrimp(self, event):
        frameSM2Tracking = MainTracking()
        frameSM2Tracking.Show()
        self.Close()
         

    def SwitchToDataFileAnalysis(self, event):
        frameDataFileAnalysis = MainTrackingDataFileAnalysis()
        frameDataFileAnalysis.Show()
        self.Close()
        

    def SwitchToIntensityTrajectoryPlots(self, event):
        IntensityTrajectoryPlots = MainTrackingIntensityOnlyPlots()
        IntensityTrajectoryPlots.Show()
        self.Close()
        



         
        
    def onOpenImageFile(self, event):
        
        self.btn_FindingMolecules.Hide()
        self.btn_IntensityTrajectoryPlots.Hide()
        self.btn_SaveIntensityTrajectoryPlots.Hide()
        
        self.btn_FIONA_Plots.Hide()
        self.btn_SaveFIONA_Plots.Hide()
        
        
        self.text_SaveIntensityTrajectoryPlots.SetLabel(' ')
        self.text_FIONA_Figsaved.SetLabel(" ")



        tp_ver = tp.__version__[:3]
        if float(tp_ver[:3]) != 0.3:
            MessageBox(self, 'The trackpy version is not 0.3.#\nUpdate trackpy', 'Error')

        
        
        plt.close('Data Image Frames')
        


          
        """
        Create and show the Open FileDialog
        """
        
        dlg2 = wx.FileDialog(
            self, message="Choose a file",
            defaultFile="",
            wildcard=wildcard,
            style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
            )
        if dlg2.ShowModal() == wx.ID_OK:
            paths = dlg2.GetPaths()
            print "You chose the following file(s):"
            for path in paths:
                print path
        print ''
        print 'dlg2 = ', dlg2
        print 'path = ', dlg2.GetPath()
        print 'path type', type(dlg2.GetPath())
        print 'filename = ', dlg2.GetFilename()
        #dlg.Destroy()
        
        self.Show()
        
        #wx.StaticText(self.MainPanel, -1, str(dlg2.GetFilename()), pos=(50, 150))
        
        filenameTemp = dlg2.GetFilename()
        self.btn_OpenImage_Text.SetLabel(str(filenameTemp))
        
        
        self.filepath = dlg2.GetPath()
        self.filenameOnly = dlg2.GetFilename()
        self.filedirectoryOnly = dlg2.GetDirectory()
        
        self.save_fig_ShowAfewFrames = self.ShowAFewFrames()
        
        self.btn_FindingMolecules.Show()
        self.btn_PlotIntensity_from_ImageJ_xy_inputs.Show()

        
        
        
        
    
    def ShowAFewFrames(self):
        frames = pims.TiffStack(self.filepath) # trackpy 0.3
        self.LoadedFrames = frames


        Nframes = len(frames) # trackpy 0.3





        print 'Nframes = ', Nframes
        self.frameLength = Nframes        
        
        self.fig_ShowAfewFrames = plt.figure('Data Image Frames')
        #fig_ShowAfewFrames.suptitle('test title', fontsize=20)
        self.fig_ShowAfewFrames.subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.55, hspace = 0.2)
        plt.figtext(0.5, 0.95, str(self.filedirectoryOnly)+'\n' + str(os.path.basename(self.filepath))
            ,ha='center', color='black', weight='normal', size='small')
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(810,30,800, 800)
        
        if Nframes > 15:
            k = 16
        else:
            k = Nframes
            
        for n in range(k):
            print 'frame # = ', n
            plt.subplot(4,4,n+1)
            plt.title('frame # ' + str(n))
            #plt.hold(True)
            plt.imshow(frames[n], origin='upper')  # for trackpy 3.0
            #plt.imshow(np.rot90(frames[n] ,3), origin='upper')
            #plt.imshow(frames[n], vmin=1000, vmax=16000)
            
        plt.ion()
        plt.show()
        
        return self.fig_ShowAfewFrames
        
        
        
        
        
        

    def FindingFeatureFunc(self, event):
        plt.close('Finding Features')
        featureSize = int(self.FeatureSize.GetValue())
        minIntensity = int(self.MinIntensity.GetValue())
        #frames = tp.TiffStack(self.filepath)
        frames = self.LoadedFrames
        
        NframeTemp = int(self.Nframe.GetValue())
        NFramesForAveraging = int(self.NofFramesForAveraging.GetValue())
        
        AnalysisConditions = ( 'Chosen F # ' + str(self.Nframe.GetValue()) + '   Ave # F: ' + str(self.NofFramesForAveraging.GetValue()) + 
        '   f size: '+ str(self.FeatureSize.GetValue()) +  '   Min I: ' + str(self.MinIntensity.GetValue()) )
        
        print '\nFingding features from frame # ', NframeTemp
        
        #framesSummed = frames[0] * 0
        #self.framesSummed = np.rot90(frames[0].T, 2) * 0
        self.framesSummed = frames[0] * 0 # for trackpy 3.0
        
        
        
        
        for n in range(NframeTemp, NframeTemp + NFramesForAveraging):
            print 'Summing frames for feature finding: frame # ', n
            #framesSummed += frames[n]       
            #self.framesSummed += np.rot90(frames[n].T, 2)
            self.framesSummed += frames[NframeTemp + n] # for trackpy 3.0
            
        self.framesSummed /= NFramesForAveraging
        
        self.f = tp.locate(self.framesSummed, featureSize, minmass=minIntensity, invert=False)
        
        #NFeatures = self.f.values.size/8
        NFeatures = len(self.f.index) #trackpy 0.3
        
        print 'NFeatures = ', NFeatures
        
        
        if NFeatures == 0:
            MessageBox(self, '0 features. Adjust the Min Intensity', 'Error')
            return
            
            
                    
        
        
        self.NofFeatures.SetLabel(str(NFeatures))
        

        
        self.N_xPixel = frames.frame_shape[0]
        self.N_yPixel = frames.frame_shape[1]
        print 'self.N_xPixel: ', self.N_xPixel
        print 'self.N_yPixel: ', self.N_yPixel      


        
        self.xpos = []
        self.ypos = []
        pixel_dn = 15
        for n in range(NFeatures):
            x_temp = self.f.values[n][0]
            y_temp = self.f.values[n][1]
            
            print 'x_temp ', x_temp
            print 'y_temp ', y_temp
            if x_temp >=pixel_dn and x_temp <= self.N_xPixel - pixel_dn and y_temp >=pixel_dn and y_temp <= self.N_yPixel - pixel_dn:
                self.xpos.append(x_temp)
                self.ypos.append(y_temp)



        self.xposImageJ = np.array(self.xpos) + 1                      # the last - 1 is for correction between ImageJ and the Crock's Grier (x,y)
        #self.yposImageJ = self.N_xPixel - np.array(self.ypos) - 0    # the last subtraction is for correction between ImageJ and Crock's Grier (x,y)
        self.yposImageJ = np.array(self.ypos) + 1    #trackpy 0.3
        
        
        
        
        
        print 'self.xposImageJ ', self.xposImageJ
        print 'self.yposImageJ ', self.yposImageJ
        
        
        plt.ion()
        plt.show()

        
        self.btn_IntensityTrajectoryPlots.Show()
        self.btn_FIONA_Plots.Show()
        


        

        intensityMax5Ave = FastIntensityHistogramFunc(np.rot90(self.framesSummed, 1), self.xpos, self.ypos)
        


        self.fig_FindingFeatures = plt.figure('Finding Features')
        plt.subplots_adjust(top = 0.88, bottom = 0.07, left = 0.04, right = 0.98, wspace = 0.2, hspace = 0.4)  
        plt.figtext(0.5, 0.92, str(self.filedirectoryOnly)+'\n' + str(os.path.basename(self.filepath)) + '\n' + AnalysisConditions + 
            '\n From frame # ' + str(NframeTemp) + '  with average of ' + str(NFramesForAveraging) + ' frames' + ' ,    Total # of features: '+ str(NFeatures)
            ,ha='center', color='black', weight='normal', size='small')

        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(510,30,1400, 700)
        #tp.annotate(f, frames[NframeTemp])
        plt.subplot(121)
        tp.annotate(self.f, self.framesSummed)
        plt.subplot(122)
        plt.hist(intensityMax5Ave)
        plt.xlabel('Max 5 Ave Intensity')




        self.Condition_text = ( 'Finding features:  Chosen F # ' + str(self.Nframe.GetValue()) + '   Ave # F: ' + str(self.NofFramesForAveraging.GetValue()) + 
        '   f size: '+ str(self.FeatureSize.GetValue()) +  '   Min I: ' + str(self.MinIntensity.GetValue()) )
        
        

        '''
        plt.figure()
        plt.subplots_adjust(top = 0.84, bottom = 0.05, left = 0.04, right = 0.98, wspace = 0.2, hspace = 0.4)  
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(10,430,1600, 600)
        plt.subplot(131)
        plt.figtext(0.5, 0.91, str(self.filedirectoryOnly)+'\n' + str(os.path.basename(self.filepath)) + '\n' + AnalysisConditions + 
            '\n From frame # ' + str(NframeTemp) + '  with average of ' + str(NFramesForAveraging) + ' frames' + ' ,    Total # of features: '+ str(NFeatures)
            ,ha='center', color='black', weight='normal', size='small')
        tp.annotate(self.f, self.framesSummed)
        plt.subplot(132)
        plt.plot(intensityMax5Ave)
        
        plt.subplot(133)
        plt.hist(intensityMax5Ave)
        '''









      
     
     
    def IntensityTrajectoryPlots(self, event):
        print '\nIntensityTrajectoryPlots'
        
        
        try:
            for n in range(len(self.fig_IntensityTrajectory)):
                plt.close(self.fig_IntensityTrajectory[n])
        except: pass
                
        
        #self.xposImageJ = [184, 250]
        #self.yposImageJ = [373, 438]
        
        ImageFrameN1 = int(self.imageFrameN1.GetValue())
        ImageFrameN2 = int(self.imageFrameN2.GetValue())
        
        
        
        (self.IntensityTrajectoryInt, self.IntensityTrajectoryAve, self.IntensityTrajectoryMax, self.IntensityTrajectoryMax5Ave,
        self.featureImages, self.featureImagesFrameN) = IntensityTrajectoryImagePlotsFunc(self.LoadedFrames, self.filepath, self.xposImageJ, self.yposImageJ, ImageFrameN1, ImageFrameN2)
        
        
        
        
        self.IntensityTrajectoryInt = np.array(self.IntensityTrajectoryInt, dtype=float)
        self.IntensityTrajectoryAve = np.array(self.IntensityTrajectoryAve, dtype=float)
        self.IntensityTrajectoryMax  = np.array(self.IntensityTrajectoryMax, dtype=float)
        self.IntensityTrajectoryMax5Ave = np.array(self.IntensityTrajectoryMax5Ave, dtype=float)
        
        

        
        
        
        
        AnalysisConditions = self.Condition_text
        
        
        
        NfeaturesPerFig = 24
        
        NFeatures = len(self.IntensityTrajectoryMax)
        Nfig = int(  math.ceil((NFeatures/float(NfeaturesPerFig))))
        TotalFeatureN = 'Total # ' + str(NFeatures)

        

        self.fig_IntensityTrajectory = [[] for _ in xrange(Nfig)]
 

       
        k = 0
        fn = 0
        for n in range(Nfig):
            print 'Intensity Trajectory Nfig n = ', n
            self.fig_IntensityTrajectory[n] = plt.figure('Images_&_IntensityTrajectory_'+ str(n), figsize = (18, 9))
            plt.clf()
            #fig_ShowAfewFrames.suptitle('test title', fontsize=20)
            #self.fig_IntensityTrajectory[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.28, hspace = 0.4)
            plt.figtext(0.5, 0.95, str(self.filedirectoryOnly)+'\n' + str(os.path.basename(self.filepath) + '\n' + AnalysisConditions   )
                ,ha='center', color='black', weight='normal', size='small')
            self.fig_IntensityTrajectory[n].text(0.05, 0.94, TotalFeatureN + '    Max 5 Ave' + '\nImage frame #: ' + str(ImageFrameN1) + ', ' + str(ImageFrameN2), ha="left", va="bottom", size="medium",color="red")
            
            
            
            
            gs = gridspec.GridSpec(8, 15)
            gs.update(bottom = 0.02, top=0.91, left=0.01, right=0.99, wspace=0.45, hspace=0.5)
            #plt.suptitle(str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly) + '\n' + AnalysisConditions, size=11)
            
            RowN = 0
            ColN = 0
            
            for m in range(NfeaturesPerFig):
                print 'm = ', m
                print 'fn + m ', fn + m
                print 'k = ', k
                print ' RowN = ', RowN
                
                if ColN%3 == 0:
                    ColN = 0
                ColN += 1
                


                axImages = plt.subplot(gs[RowN, (ColN-1)*5+0], aspect='equal')
                plt.tick_params(labelsize=8)
                axImages.set_title("M# "+str(k), size =9)
                plt.axis('off')
                plt.imshow(np.rot90(self.featureImages[fn+m][0].reshape(15, 15),3), interpolation = 'None', cmap=plt.cm.jet, origin='upper')
                                   #extent=(x3.min(), x3.max(), y3.min(), y3.max()) )

                axImages = plt.subplot(gs[RowN, (ColN-1)*5+1], aspect='equal')
                plt.tick_params(labelsize=8)
                axImages.set_title("(" + str(int(self.xposImageJ[k])) + ', ' + str(int(self.yposImageJ[k])) + ')', size =9)
                plt.axis('off')
                plt.imshow(np.rot90(self.featureImages[fn+m][1].reshape(15, 15),3), interpolation = 'None', cmap=plt.cm.jet, origin='upper')
                                   #extent=(x3.min(), x3.max(), y3.min(), y3.max()) )
            

                axTotal_Intensity = plt.subplot(gs[RowN, (ColN-1)*5+2:(ColN-1)*5+5])
                plt.tick_params(labelsize=8)
                
                axTotal_Intensity.plot(self.IntensityTrajectoryMax5Ave[fn + m])
                


                
                
                if (m+1)%3 == 0:
                    RowN += 1
                                    
                k += 1
                if k%NfeaturesPerFig ==0:
                    fn += NfeaturesPerFig
                if k == NFeatures:
                    break
                
                
        plt.ion()
        plt.show()
        
        
        
        self.btn_SaveIntensityTrajectoryPlots.Show()
        
        
        
        
        
        

    def SaveIntensityTrajectoryPlotsData(self, event):
        prefix = self.FilePrefix.GetValue()
        if prefix != '':
            prefix += '_'
           
            
            
        figFileName = prefix +  self.filenameOnly[:30]
        todayDate = time.strftime("%Y%m%d")

        try:
            self.fig_FindingFeatures.savefig(figFileName +'_' + todayDate +  '_Image&IntensityPlots1_Features.png')
        except: pass


        try:
            for n in range(len(self.fig_IntensityTrajectory)):
                self.fig_IntensityTrajectory[n].savefig(figFileName + '_' + todayDate + '_Image&IntensityPlots2_' + str(n)+ '.png', dpi = 100)
                plt.close(self.fig_IntensityTrajectory[n])
                
                #plt.close(self.fig_IntensityTrajectory[n])
            print '\nImage & Intensity Trajectory figures (data) are saved'
            #wx.StaticText(self.MainPanel, -1, 'Image & Intensity Trajectory figures saved', pos=(50, 380)) 
            self.text_SaveIntensityTrajectoryPlots.SetLabel('Image & Intensity Trajectory figures (data) saved')
            
        except:
            print '\nImage & Intensity Trajectory figures (data) are NOT saved'     
            #wx.StaticText(self.MainPanel, -1, 'Error', pos=(50, 380)) 
            self.text_SaveIntensityTrajectoryPlots.SetLabel('Errors')
            
            
            
        plt.close('Finding Features')
        plt.close('Data Image Frames')




        if self.cb_save_data_files_yn.GetValue():
            IntensityFileName = prefix + prefix2 + self.filenameOnly[:30] + '_' + todayDate + '_Intensity' 
            for n in range(len(self.IntensityTrajectoryMax5Ave)):
                xytemp = '(' + str(int(self.xposImageJ[n])) + ',' + str(int(self.yposImageJ[n])) +')'
                ff = open(IntensityFileName  + '_' + str(n) + '_' + xytemp + '.txt','w')
                for i in range(len(self.IntensityTrajectoryMax[n])):
                    ff.write(str(self.IntensityTrajectoryMax[n][i]) + '\n'   )    
                ff.close()       
    
    
    
    




    
     
    def FIONA_Plots(self, event):
        print '\nFIONA Plots'


        featureSize = int(self.FeatureSize.GetValue())
        minIntensity = int(self.MinIntensity.GetValue())
        #frames = tp.TiffStack(self.filepath)
        frames = self.LoadedFrames
        
        NframeTemp = int(self.Nframe.GetValue())
        NFramesForAveraging = int(self.NofFramesForAveraging.GetValue())
        

        #framesSummed = np.rot90(frames[0], -1) * 0
        framesSummed = frames[0] * 0 # for trackpy 3.0

        
        for n in range(NframeTemp, NframeTemp + NFramesForAveraging):
            print 'Summing frames for feature finding: frame # ', n
            #framesSummed += frames[n]       
            #framesSummed += np.rot90(frames[n], -1)
            framesSummed += frames[n] # for trackpy 3.0
            
        framesSummed /= NFramesForAveraging
        
        NFeatures = 0
        
        print 'minIntensity ', minIntensity

        
        self.N_xPixel = frames.frame_shape[0]
        self.N_yPixel = frames.frame_shape[1]
        print 'self.N_xPixel: ', self.N_xPixel
        print 'self.N_yPixel: ', self.N_yPixel        
        

        
        while 14 > NFeatures or NFeatures > 30:
                
            f = tp.locate(framesSummed, featureSize, minmass=minIntensity, invert=False)
            #NFeaturesTemp = f.values.size/8
            NFeaturesTemp = len(self.f.index) #trackpy 0.3
            
            print 'NFeaturesTemp', NFeaturesTemp
            

            xposTemp = []
            yposTemp = []
            
            pixel_dn = 30
            
            
            for n in range(NFeaturesTemp):
                print n
                x_temp = f.values[n][0]
                y_temp = f.values[n][1]
                
                print 'x_temp ', x_temp
                print 'y_temp ', y_temp
                if x_temp >=pixel_dn and x_temp <= self.N_xPixel - pixel_dn and y_temp >=pixel_dn and y_temp <= self.N_yPixel - pixel_dn:
                    xposTemp.append(x_temp)
                    yposTemp.append(y_temp)


        
            NFeatures = len(xposTemp)
            print 'NFeatures = ', NFeatures
            
            
            if NFeatures<14:
                minIntensity *= 0.9
                print 'minIntensity ', minIntensity
                
            if 30 < NFeatures < 50:
                minIntensity *= 1.1
                print 'minIntensity ', minIntensity
                
            if 50 <= NFeatures:
                minIntensity *= 2
                print 'minIntensity ', minIntensity  
            

        




        
        AnalysisConditions = ( 'Finding features:  Chosen F # ' + str(self.Nframe.GetValue()) + '   Ave # F: ' + str(self.NofFramesForAveraging.GetValue()) + 
        '   f size: '+ str(self.FeatureSize.GetValue()) +  '   Min I: ' + str(self.MinIntensity.GetValue()) )
        
        



        
        
        xpos = np.array(xposTemp[:14]) + 1
        ypos = np.array(yposTemp[:14]) + 1
 
        
        NframesFiona = int(self.NfastFIONAframe.GetValue())
        
        CenterImages, fx0, fy0, FIONA_frameN, centerIntensityMax5Ave = FastFIONA(self.LoadedFrames, self.filepath, xpos, ypos, int(self.Nframe.GetValue()), int(self.Nframe.GetValue()) + NframesFiona -1, 1.1, 0.6, 10, 0.2, 0)
        
        
        
        #print 'CenterImages[0] ', CenterImages[0]
        print 'fx0 ', fx0

        
        NfeaturesPerFig = 14
        NfeaturesPerRow = 2

        
        NFeatures = len(xpos)
        Nfig = int(  math.ceil((NFeatures/float(NfeaturesPerFig))))
        TotalFeatureN = 'Total # ' + str(NFeatures)

        

        self.fig_FIONA_Plots = [[] for _ in xrange(Nfig)]
 

       
        k = 0
        fn = 0
        for n in range(Nfig):
            print 'Intensity Trajectory Nfig n = ', n
            self.fig_FIONA_Plots[n] = plt.figure('FIONA_Images_&_IntensityTrajectory_'+ str(n), figsize = (18, 9))
            plt.clf()
            #fig_ShowAfewFrames.suptitle('test title', fontsize=20)
            #self.fig_FIONA_Plots[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.28, hspace = 0.4)
            plt.figtext(0.5, 0.95, str(self.filedirectoryOnly)+'\n' + str(os.path.basename(self.filepath) + '\n' + AnalysisConditions   )
                ,ha='center', color='black', weight='normal', size='small')
            self.fig_FIONA_Plots[n].text(0.05, 0.94, TotalFeatureN + '    Max 5 ave' + '\nImage frame #: ' + str(int(self.Nframe.GetValue())) + ', ' + str(int(self.Nframe.GetValue())+1), ha="left", va="bottom", size="medium",color="red")
            
            
            
            
            gs = gridspec.GridSpec(7, 14)
            gs.update(bottom = 0.05, top=0.91, left=0.01, right=0.99, wspace=0.5, hspace=0.7)
            #plt.suptitle(str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly) + '\n' + AnalysisConditions, size=11)
            
            RowN = 0
            ColN = 0
            
            for m in range(NfeaturesPerFig):
                print 'm = ', m
                print 'fn + m ', fn + m
                print 'k = ', k
                print ' RowN = ', RowN
                
                if ColN%NfeaturesPerRow == 0:
                    ColN = 0
                ColN += 1
                


 
                axImages = plt.subplot(gs[RowN, (ColN-1)*7+0], aspect='equal')
                plt.tick_params(labelsize=8)
                axImages.set_title("M# "+str(k), size =9)
                plt.axis('off')
                plt.imshow(np.rot90(CenterImages[fn+m][0].reshape(15, 15),3), interpolation = 'None', cmap=plt.cm.jet, origin='bottom')
                                   #extent=(x3.min(), x3.max(), y3.min(), y3.max()) )

                axImages = plt.subplot(gs[RowN, (ColN-1)*7+1], aspect='equal')
                plt.tick_params(labelsize=8)
                axImages.set_title("(" + str(int(xpos[k])) + ', ' + str(int(ypos[k])) + ')', size =9)
                plt.axis('off')
                plt.imshow(np.rot90(CenterImages[fn+m][1].reshape(15, 15),3), interpolation = 'None', cmap=plt.cm.jet, origin='bottom')
                                   #extent=(x3.min(), x3.max(), y3.min(), y3.max()) )
            
                axTotal_Intensity = plt.subplot(gs[RowN, (ColN-1)*7+2:(ColN-1)*7+5])
                plt.tick_params(labelsize=8)
                plt.locator_params(nbins=4)
                axTotal_Intensity.plot(FIONA_frameN, centerIntensityMax5Ave[fn + m])
                
                
                
                axFX0 = plt.subplot(gs[RowN, (ColN-1)*7+5], aspect='equal')
                plt.tick_params(labelsize=8)
                axFX0.set_title("(" + str(int(xpos[k])) + ', ' + str(int(ypos[k])) + ')', size =9)
                #plt.axis('off')
                plt.scatter(fx0[m + fn], fy0[m + fn], marker = 'x', c = range(len(fx0[m + fn])),  s=35, vmin=0, vmax= len(fx0[m + fn]))
                plt.xlim( (np.around(fx0[k][0]) - 3, np.around(fx0[k][0]) + 3 ) )
                plt.ylim( (np.around(fy0[k][0]) - 3, np.around(fy0[k][0]) + 3 ) )
                plt.xticks(np.arange(np.around(fx0[k][0]) - 2, np.around(fx0[k][0]) + 4, 1 ), rotation='vertical')
                plt.yticks(np.arange(np.around(fy0[k][0]) - 2, np.around(fy0[k][0]) + 4, 1 ))



                
                
                if (m+1)%NfeaturesPerRow == 0:
                    RowN += 1
                                    
                k += 1
                if k%NfeaturesPerFig ==0:
                    fn += NfeaturesPerFig
                if k == NFeatures:
                    break
                
                
        
        plt.show()
        
        
        
        self.btn_SaveFIONA_Plots.Show()
        
        
  





    def SaveFIONA_Plots(self, event):
        prefix = self.FilePrefix.GetValue()
        if prefix != '':
            prefix += '_'

        figFileName = prefix + self.filenameOnly[:50]
        todayDate = time.strftime("%Y%m%d")


        try:
            
            for n in range(len(self.fig_FIONA_Plots)):
                self.fig_FIONA_Plots[n].savefig('FIONA_' + figFileName + '_' + todayDate + '_' + str(n)+ '.png', dpi = 100)
                plt.close(self.fig_FIONA_Plots[n])
                #plt.close(self.fig_IntensityTrajectory[n])
            print '\nImages_Intensity_FIONA figures are saved'
            #wx.StaticText(self.MainPanel, -1, 'Images_Intensity_FIONA_ figures saved', pos=(450, 380)) 
            self.text_FIONA_Figsaved.SetLabel("Images_Intensity_FIONA_ figures saved")
            
            
            
        except:
            print '\nImages_Intensity_FIONA figures are NOT saved'     
            #wx.StaticText(self.MainPanel, -1, 'Error', pos=(450, 380)) 
            self.text_FIONA_Figsaved.SetLabel("Error")
            




        plt.close('Finding Features')
        plt.close('Data Image Frames')






    def Create_ImageJ_xy_inputs(self, event):
        print ''
        
        
        try:
            for n in range(len(self.ImageJ_x_for_intensity)):
                self.ImageJ_xy_for_intensity_N[n].Destroy()
                self.ImageJ_x_for_intensity[n].Destroy()
                self.ImageJ_y_for_intensity[n].Destroy()
        except: pass        
            
        
        
        
        
        self.N_ImageJ_xy_inputs = int(self.N_ImageJ_features.GetValue())

        self.ImageJ_xy_for_intensity_N = [[] for _ in xrange(self.N_ImageJ_xy_inputs)] 
        self.ImageJ_x_for_intensity = [[] for _ in xrange(self.N_ImageJ_xy_inputs)] 
        self.ImageJ_y_for_intensity = [[] for _ in xrange(self.N_ImageJ_xy_inputs)]
        
        
        
        rn = 0
        cn = 0
        countertemp = 0
        for n in range(self.N_ImageJ_xy_inputs):
            self.ImageJ_xy_for_intensity_N[n] = wx.StaticText(self.MainPanel, -1, '# ' + str(n) + '.', pos=(50 + 120*cn, 603 + 25*rn), size=(-1,-1))
            self.ImageJ_x_for_intensity[n] = wx.TextCtrl(self.MainPanel, -1, '0', pos=(80 + 120*cn, 600 + 25*rn), size=(30,-1))
            self.ImageJ_y_for_intensity[n] = wx.TextCtrl(self.MainPanel, -1, '0', pos=(113 + 120*cn, 600 + 25*rn), size=(30,-1))
            countertemp += 1
            cn += 1
            if countertemp % 5 == 0:
                rn += 1
                cn = 0
            
        
        

    def PlotIntensity_from_ImageJ_xy_inputs(self, event):
        print ''
        
        self.xposImageJ = []
        self.yposImageJ = []
        for n in range(self.N_ImageJ_xy_inputs):
            xtemp = int(self.ImageJ_x_for_intensity[n].GetValue())
            ytemp = int(self.ImageJ_y_for_intensity[n].GetValue())
            
            if xtemp != 0 or ytemp != 0:
                self.xposImageJ.append(xtemp)
                self.yposImageJ.append(ytemp)
        

        if len(self.xposImageJ) == 0:
            MessageBox(self, 'No ImageJ x,y data', 'Error')
            return
        

        
        self.Condition_text = 'Directly from ImageJ (x,y)'
        
        self.IntensityTrajectoryPlots('')
        
        








      

#----------------------------------------------------------------------
# Run the program
if __name__ == "__main__":
    app = wx.App(False)
    if StartingPanel == 1:
        frame = MainTracking()
        frame.Show()
    elif StartingPanel == 2:
        frame = MainTrackingDataFileAnalysis()
        frame.Show()
        
    elif StartingPanel == 3:
        frame = MainTrackingIntensityOnlyPlots()
        frame.Show()
        
    else:
        print 'Nothing'
        pass
    
    #for test
    
    '''
    x = np.arange(0.0, 31, 0.01) 
    y = np.cos(x)
    y2 = np.exp(-x)
    
    acf = AutoCorrFunc(y)    
    acf2 = AutoCorrelationFunction(y)
    acf3 = AutoCorrFunc(x * 0 +100)
    acf22 =   AutoCorrFunc(y2)
    
    plt.figure()
    plt.plot(x,y)
    plt.plot(x, acf, '--')
    plt.plot(x, acf2)
    plt.plot(x, acf3)
    plt.plot(x, y2, 'o')
    plt.plot(x, acf22)
    
    plt.show()
    '''
    
    
    
    app.MainLoop()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
