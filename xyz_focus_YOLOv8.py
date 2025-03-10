# -*- coding: utf-8 -*-
"""
@author: FGLab (An Mei) also see 
https://github.com/TheImagingSource/IC-Imaging-Control-Samples/tree/master/Python/tisgrabber/samples

Last edited: Thu Jun 27 2024

This code images droplets in a timelapse.

Focal plane is found using a deep learning model, and slices are then taken
above and below the focal plane (starting x um below and moving up). The focal
plane is found dynamically, and therefore can correct for z drift.

Hough transform is used to find the centre of the drop, and allows correcting
for xy drift.

****V4 notes: using the recording history, can detect if the focus has been going
              back and forth between slightly above and below for a while (26 
              iterations) and break it. Can also determine if it's been going in
              one direction for a while (45 iterations) and still hasn't found
              the focus. If so, corrects back the other way (manually sets focus 
              for the other direction and returns back the other way to try again. 
              Also gets current time imaged (normalised to T=0s) in addition to timepoint.
"""



### IMPORT NECESSARY MODULES

import os
import serial
from datetime import datetime as dt
import time as t
import ctypes, re
import numpy as np
import pandas as pd

import cv2 as cv2
import tisgrabber as IC
from proscan import PriorStage

import ultralytics
from ultralytics import YOLO
ultralytics.checks()



### PARAMETERS TO ADJUST
Width=640 #change image width here
Height=480

total_nb=15000 #total number of time points
time_interval=45 #time interval between time points
restart_time=18 #change this number if you have restarted it to the last timepoint you had before restarting

recording = True #change for whether or not the focus history from autofocus is being recorded
focus_step=7 #step size for focus adjustments in 0.1 microns
step_size=7 #step size in 0.1 microns
nb_slices=12 #number of slices above and below focal plane (total slices would be this * 2)

exposure_time=1/9800#1/6803
gain=14 #5
tolerance=0.15 #chi-square:40000 // KL=4000 //Bhatta =0.05
############################



### DIRECTORIES & FILES

main_dir = os.getcwd()

# Read number of droplets
with open('positions.txt') as f:
    text = f.readlines()
    Nb_drops = len(text) # number of droplets in the positions.txt file
    print(f'Number of drops: {str(Nb_drops)}')

# Create folders and store droplet positions
x_drop = [0 for i in range(Nb_drops)]  #stores x coordinates for centres of drops
y_drop = [0 for i in range(Nb_drops)]  #stores y coordinates for centres of drops
for n in range(Nb_drops): #loops through each of the drops
    # Creates individual folders for drops
    newpath = f'{main_dir}\drop{str(n + 1)}' 
    if not os.path.exists(newpath): os.makedirs(newpath)
    # Finds the xy positions of the drops
    pos = re.findall(r'-?\d+', str(text[n]))
    x_drop[n] = float(pos[0]) # Position of drops multiply by 100
    y_drop[n] = float(pos[1])
     
# Load the pre-trained YOLO model
model = YOLO('bacteria-autofocus-model.pt')
remapping = [5,4,1,2,3] # remap classes to have in order of focus level to determine weight
if restart_time == 0: record_focus = pd.DataFrame(columns=['Timepoint', 'Time (s)', 'Drop', 'Focus Level', 'Confidence', 'Weight', 'Movement']) #dataframe to record focus history
else: record_focus = pd.read_csv('focus_history.csv')



### CONNECT TO CAMERA AND INITIALISE LIVESTREAM

def send_command(command, sleep=False):
    """A function to automatically send commands to the serial port for imaging."""
    ser.write(command.encode())
    ser.flush() #make sure all is sent
    if sleep==True: t.sleep(0.1) #wait for the command to be processed

ser = serial.Serial(
     port='COM3',
     baudrate=9600,
     timeout=1,
     bytesize=8,
     parity=serial.PARITY_NONE,
     stopbits=serial.STOPBITS_ONE,
     rtscts=False,
     dsrdtr=False,
     xonxoff=False
     )

print(ser.get_settings())
print(ser.name) 
print("OPENING port.")

if not ser.isOpen():
    ser.open()

send_command("\r") #return carriage to initialise comm
send_command("\r") #second return carriage
serialString = ser.readlines(16) #read output
print(serialString)
send_command("\r")

serialString = ser.read(16)
position = bytes.decode(serialString) #convert bytes to string
position_ini = int(position[4:10]) #crop useful number as it comes in form 0,0,number
position_ini2 = bytes(str(position_ini), 'utf-8') #convert bytes to string
print(position_ini)
         
step = bytes(str(step_size), 'utf-8') #convert to string
send_command(f"C {step}\r") #sets step size for focus motor

ic = ctypes.cdll.LoadLibrary("./tisgrabber_x64.dll")

print('ok')

IC.declareFunctions(ic)
ic.IC_InitLibrary(0)

hGrabber = IC.openDevice(ic)
ic.IC_SetPropertySwitch(hGrabber, IC.T("Exposure"), IC.T("Auto"), 0)
ic.IC_SetPropertyAbsoluteValue(hGrabber, IC.T("Exposure"), IC.T("Value"),
                                ctypes.c_float(exposure_time))
ic.IC_SetPropertySwitch(hGrabber, IC.T("Gain"), IC.T("Auto"), 0)
ic.IC_SetPropertyAbsoluteValue(hGrabber, IC.T("Gain"), IC.T("Value"),
                                ctypes.c_float(gain))

ic.IC_SetVideoFormat(hGrabber, IC.T(f"Y800 ({str(Width)}x{str(Height)})"))
    
print('Camera properties initialized \n\n\n')

# Start live video feed..
if ic.IC_IsDevValid(hGrabber):
    ic.IC_StartLive(hGrabber, 1)  

im_local = np.zeros((640,480), float)

# Communication with Prior Stage
p = PriorStage("COM4")



### DEFINE FUNCTIONS FOR THE MAIN LOOP

def capture_image(image_path):
    """Function to capture and save an image."""
    if ic.IC_SnapImage(hGrabber, 2000) == IC.IC_SUCCESS:
        ic.IC_SaveImage(hGrabber, IC.T(image_path), IC.ImageFileTypes['JPEG'], 100)

def Hough_T(image_path,x0,y0,n):
    """Function that uses the Hough transform to determine xy drift and correct it."""
    capture_image(image_path) #acquire one image
    # Use cv2 to perform Hough Transform      
    img = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)   
    circles = cv2.HoughCircles(img,cv2.HOUGH_GRADIENT,1,200,
    param1=80,param2=40,minRadius=150,maxRadius=280) #param2=40 
    new_x_drop=x0 #in case no new centre is found
    new_y_drop=y0 #in case no new centre is found    
    if circles is not None:
        circles = np.uint16(np.around(circles))
        #find new position of centre based on list of circles
        for i in circles[0,:]:
            if abs(Width/2-i[0]) <= 75 and abs(Height/2 - i[1]) <= 75:
                print(f'Correction found with x={str(i[0])}, y={str(i[1])}')
                cimg = cv2.cvtColor(img,cv2.COLOR_GRAY2BGR)
                # draw the outer circle
                cv2.circle(cimg,(i[0],i[1]),i[2],(0,255,0),2)
                # draw the center of the circle
                cv2.circle(cimg,(i[0],i[1]),2,(0,0,255),3)
                cv2.imwrite(image_dir+'\Hough.jpeg',cimg)
                new_x_drop=x_drop[n]+(Width/2-i[0])/344*60*100 #x_drop[n] is the top left absolute position, 60um =344 pixels, scale =0.1 microns
                new_y_drop=y_drop[n]-(Height/2-i[1])/344*60*100
                break
    return (new_x_drop,new_y_drop)

def determine_focus_level(model, image_path, remapping):
    """Function to determine the focus of the image using the deep learning model."""
    outputs = model(image_path, verbose=False, device=0)
    results = pd.DataFrame({'Class': [outputs[0].cpu().names[x] for x in outputs[0].cpu().probs.top5],
                            'Class ID': [remapping[x] for x in outputs[0].cpu().probs.top5],
                            'Confidence': outputs[0].cpu().probs.top5conf})
    weight = np.sum(results['Confidence'] * results['Class ID']) / 10
    focus_level = results['Class'].iloc[0]
    confidence = results['Confidence'].iloc[0]
    return focus_level, confidence, weight

def adjust_focus_based_on_yolo(focus_level, weight, step_size):
    """Function to adjust Z-axis based on YOLO focus level."""
    if round(weight, 3)>=0.28 and round(weight, 3)<=0.35: 
        focus_level='InFocus' #use weights to determine focus as well
    elif focus_level == 'AboveGreater3' or focus_level == 'AboveLess3':
        send_command(f'D {step_size}\r', sleep=True)  #move down by step size um
    elif focus_level == 'BelowGreater3' or focus_level == 'BelowLess3':
        send_command(f'U {step_size}\r', sleep=True)  #move up by step size um
    return focus_level 

def record_focus_history(record_focus, timepoint, time, drop, focus_level, confidence, weight, step_size):
    """Function to record the adjustments made to detect focus level."""
    movement=0
    if round(weight, 3)>=0.28 and round(weight, 3)<=0.35: movement=0
    elif focus_level == 'AboveGreater3' or focus_level == 'AboveLess3': movement = 0-step_size
    elif focus_level == 'BelowGreater3' or focus_level == 'BelowLess3': movement = step_size
    record_focus.loc[len(record_focus)] = [timepoint, time, drop, focus_level, confidence, weight, movement]    

def save_focus_history(record_focus, main_dir, recording):
    """Function to save the recorded focus data from the autofocus in a csv."""
    if recording==True:
        record_focus.to_csv(f'{main_dir}\\focus_history.csv', index=False) #save data
        # Create a summary of all the focus data
        record_focus = record_focus.groupby(['Timepoint', 'Drop']).agg({'Time (s)':'first', 'Movement': ['count', 'sum']})
        record_focus = record_focus.droplevel(0, axis=1).rename({'first':'Time (s)', 'count':'Number of Movements', 'sum':'Total Movement (in 0.1um)'}, axis=1)
        record_focus = record_focus.reset_index()
        record_focus.to_csv(f'{main_dir}\\focus_summary.csv', index=False)
    

        

### MAIN LOOP

for timepoint in range(total_nb):
    timepoint = timepoint + restart_time #will update if restarted
    print(f'Timepoint: {str(timepoint)} ({timepoint * time_interval}s)')
    
    for drop_idx in range(Nb_drops):
    
        p.move_to({'x': x_drop[drop_idx], 'y': y_drop[drop_idx]}) #move to next drop
        t.sleep(1)
        print(f'\n\n Droplet number: {str(drop_idx+1)} \n')
        image_dir = f"{main_dir}\drop{str(drop_idx+1)}" #new image folder for each drop
        
        # Focus correction
        focus_level='OutofFocus'
        wrong_direction = False
        rounds = 0
        if timepoint==0 or timepoint==restart_time: focus_start_time = dt.now()
        time = (dt.now() - focus_start_time).total_seconds()
        print('Finding focus to correct drift.')
        while focus_level != 'InFocus':
            # Detect focus using machine learning
            image_path = f'{image_dir}\\focus_test.jpeg'
            capture_image(image_path)
            focus_level, confidence, weight = determine_focus_level(model, image_path, remapping)
            # Correct detection - generally based on autofocus (may be very out-of-focus, but generally only drifts a little)
            if rounds == 26 and (np.sum(record_focus.tail(25)['Movement'])<focus_step*3 and np.sum(record_focus.tail(25)['Movement'])>-focus_step*3): #if gone back and forth a lot
                break #end it early
            if rounds == 45 and (np.sum(record_focus.tail(44)['Movement'])>focus_step*5 or np.sum(record_focus.tail(44)['Movement'])<-focus_step*5): #therefore if its taken >30 rounds going in one direction, its probably going in the wrong direction
                # When very out-of-focus, it can't always tell the difference between above or below the plane, so manually tell it to go the opposite way
                wrong_direction = True
                focus_this_step = record_focus.tail(44)
                adjust_step = (len(focus_this_step) + 2)*focus_step #first fix to go back the other way a bit
                if focus_level=='BelowGreater3': focus_level='AboveGreater3'; weight=np.nan; confidence=np.nan #manually adjust values for recording
                if focus_level=='AboveGreater3': focus_level='BelowGreater3'; weight=np.nan; confidence=np.nan  
                print(f'The image is {focus_level} ({weight:.3}), with {confidence*100:.3}% confidence.)')
                focus_level = adjust_focus_based_on_yolo(focus_level, weight, adjust_step)
                record_focus_history(record_focus, timepoint, time, drop_idx+1, focus_level, confidence, weight, step_size) #save history of focus correction
            elif rounds > 46 and wrong_direction==True: #keep going the other direction if still detecting below the plane after manual correcting
                if focus_level == 'BelowGreater3': focus_level='AboveGreater3'; weight=np.nan; confidence=np.nan #manually adjust values for recording
                if focus_level=='AboveGreater3': focus_level='BelowGreater3'; weight=np.nan; confidence=np.nan 
                print(f'The image is {focus_level} ({weight:.3}), with {confidence*100:.3}% confidence.)')
                focus_level = adjust_focus_based_on_yolo(focus_level, weight, focus_step)
                record_focus_history(record_focus, timepoint, time, drop_idx+1, focus_level, confidence, weight, focus_step) #save history of focus correction
            else: #else just do a normal adjustment based on the detected focus level
                print(f'The image is {focus_level} ({weight:.3}), with {confidence*100:.3}% confidence.)')
                focus_level = adjust_focus_based_on_yolo(focus_level, weight, focus_step)
                record_focus_history(record_focus, timepoint, time, drop_idx+1, focus_level, confidence, weight, focus_step) #save history of focus correction
            save_focus_history(record_focus, main_dir, recording)
            if rounds > 65: break #if gets trapped in a loop just break
            rounds+=1
            
        # XY drift correction
        print('\n XY drift correction')
        image_path = f'{image_dir}\Hough_test.jpeg'
        (new_x_drop,new_y_drop)=Hough_T(image_path,x_drop[drop_idx],y_drop[drop_idx],drop_idx)
        x_drop[drop_idx]=new_x_drop
        y_drop[drop_idx]=new_y_drop       
        p.move_to({'x': x_drop[drop_idx], 'y': y_drop[drop_idx]}) #moves to next drop
        t.sleep(1)
        print('\n')
        
        # Take image slices
        send_command(f'D {step_size*nb_slices}\r', sleep=True) #go below focal plane to image up
        for x in range(nb_slices*2):               
            image_path = f'{image_dir}\\foo{str(timepoint)}_{str(x)}.jpeg' #change j+ 'X' here if needed on restart
            capture_image(image_path)
            send_command(f'U {step_size}\r', sleep=True)
            print(f"Imaged slice {str(x)}")
        
        # Return to roughly focal plane for next drop
        send_command(f'D {step_size*nb_slices}\r', sleep=True)
        
    # Wait for the next time point
    print(f"\n Waiting {time_interval}s for next time step. \n\n")
    t.sleep(time_interval)
        

# CLOSE SERIAL PORT AND STOP LIVE FEED TO END
ser.close()
if ic.IC_IsDevValid(hGrabber):
    ic.IC_StopLive(hGrabber)