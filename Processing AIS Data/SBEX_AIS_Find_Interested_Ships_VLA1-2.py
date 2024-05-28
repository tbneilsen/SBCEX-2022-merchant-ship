# -*- coding: utf-8 -*-
"""
Created on Jun 14 2022
Last Modified: Sept 6 2023

@author: Alexandra Hopps McDaniel

Purpose: 
Sorts through a download of AIS data and finds all ships of interest for further analysis

Output: 
Formatted spreadsheet with all interested ships and their pertinent information for each VLA

Note: This code assume that the AIS file downloaded has been sorted in alphebetic order according to the vessel names
"""
import csv
import numpy as np
import os
import matplotlib
import pandas as pd
import datetime

# Reads in the AIS csv file 
# You can download the AIS file from: https://marinecadastre.gov/accessais/
# The file will be emailed a few days after a request has been submitted
# Once downloaded, make sure the file is saved wiith the list sorted alphabetically according to vessel name  

# Update the path to the AIS download
file=r'C:\Users\alexh\Documents\UW_Acoustics\Codes\sbc-vla-data\processing_AIS_data\AIS_download_2022_VLA1-2.csv' 
df = pd.read_csv(file)

# Creating arrays with the information from the headers of the CSV file 
# Check to make sure all arrays are the same size 
# If the square is blank in the AIS file, that place in the array will be filled with 'nan'
Vessel=df['VesselName']
date_time=df['BaseDateTime']
imo=df['IMO']
SOGknots=df['SOG'] # Ship's Speed over Ground in knots 
latitude=df['LAT']
longitude=df['LON']
length=df['Length']
width=df['Width']
draft=df['Draft']



def formatTimeStamp(date_time):
    
    '''
    This function takes the format of date and time in the downloadable AIS CSV 
    and converts it to the special time stamp format. 
    The format is as follows: YYDDDHHMMSS, wheree DDD is the three digit Julian date.

    The function returns and array with the entire column of formatted date/time stamps. 
    '''

    time_stamp=np.zeros_like(date_time) # creates an array same length as date_time 
   
    for i in range(len(date_time)):
        day=int(date_time[i][8:10])
        month=int(date_time[i][5:7])
        year=int(date_time[i][0:4])
        # python function for converting year, month, and day to julian day
        date_obj=datetime.date(year, month, day)
        julian_day=date_obj.timetuple().tm_yday
        
        # The following line formats the time step information
        time_stamp[i]=date_time[i][2:4]+str(julian_day)+date_time[i][11:13]+date_time[i][14:16]+date_time[i][17:19] 
        
    return time_stamp

def findInterestedShips(latVLA,longVLA,startTime,endTime,radius,formattedTime):

    '''
    This function goes through the complete array of AIS data and finds the ships that fall within the
    time period and are found less than 15 km from the arrays. 
    It tracks the distance of each ship to the arrays and saves only the time stamp and pertinent ship information
    at CPA (closest point of approach). 
    '''

    distance=[]
    vessel=[]
    imoShip=[]
    time=[]
    latShip=[]
    longShip=[]
    sogShip=[]
    lengthShip=[]
    widthShip=[]
    draftShip=[]
    
    r=6.3781e3 # radius of spherical Earth in km
   
    # Delete irrelevant data and find the distance between each ship and VLA
    for i in range(len(Vessel)):
        if formattedTime[i]>startTime and formattedTime[i]<endTime: # Filters out ships outside of time frame 
            
            # Use the haversine formula to compute the distance of each ship to the VLAs
            # Can find information on the haversine formula: https://en.wikipedia.org/wiki/Haversine_formula
            # Note: Latitudes and Longitudes must be in radians for the equation to work
            haversine=2*r*np.arcsin(np.sqrt(np.sin(np.radians(latVLA-latitude[i])/2)**2+np.cos(np.radians(latVLA))*np.cos(np.radians(latitude[i]))*np.sin(np.radians(longVLA-longitude[i])/2)**2))
                    
            if int(haversine)<radius: # filters out all ships that are farther than the set radius from the VLA
                distance.append(haversine)
                vessel.append(Vessel[i])
                time.append(formattedTime[i])
                latShip.append(latitude[i])
                longShip.append(longitude[i])
                sogShip.append(SOGknots[i])
                lengthShip.append(length[i])
                widthShip.append(width[i])
                draftShip.append(draft[i])
                imoShip.append(imo[i])
            
    # The list now has repeats of each ship but filtered 
    # This list saves each ship name that meets the criteria without repeats
    interestedVessels=[]
    interestedVessels.append(vessel[0]) 
    for i in range(len(vessel)-1):
        if vessel[i+1]!= vessel[i]:
            interestedVessels.append(vessel[i+1])
            
    cpa=[]
    timeStamp=[]
    sog=[]
    Length=[]
    Width=[]
    Draft=[]
    IMOship=[]
    # This loop sorts through the filtered list of ships, finds the CPA and time stamp, and saves the information in array in the same order as the interested ships array
    for i in range(len(interestedVessels)):
        d=[]
        js=[]
        for j in range(len(distance)):
            if vessel[j]==interestedVessels[i]:
                d.append(distance[j])
                js.append(j)
                
        index=d.index(min(d))
        cpa.append(d[index])
        timeStamp.append(time[js[index]])
        sog.append(sogShip[js[index]])
        Length.append(lengthShip[js[index]])
        Width.append(widthShip[js[index]])
        Draft.append(draftShip[js[index]])
        IMOship.append(imoShip[js[index]])
        
        
    return interestedVessels,IMOship,timeStamp,cpa,sog,Length,Width,Draft

# Information regarding the  deployment of VLA 1 and 2

# LOCATION 1
# In water from 23-May to 26-May
lat_VLA1_location1=40.47
long_VLA1_location1=-70.5971
lat_VLA2_location1=40.4418
long_VLA2_location1=-70.5272
# Start and end times 
startTime_VLA1_location1=22143132500 # 23-May-22 start time was 1319Z, but rounded up to 1325Z
endTime_VLA1_location1=22146134000 # 26-May-22 end time was 1354Z, but rounded down to 1340Z
startTime_VLA2_location1=22143161000 # 26-May-22 start time was 1603Z, but rounded up to 1610Z
endTime_VLA2_location1=22146121000 # 26-May-22 end time was 1220Z, but rounded down to 1210Z

# LOCATION 2
# In water from
lat_VLA1_location2=40.0485
long_VLA1_location2=-70.8842
lat_VLA2_location2=39.954
long_VLA2_location2=-70.7708
# Start and end times 
startTime_VLA1_location2=22147130000 # 27-May-22 start time was 1250Z, but rounded up to 1300Z
endTime_VLA1_location2=22150124000 # 30-May-22 end time was 1250Z, but rounded down to 1240Z
startTime_VLA2_location2=22147151000 # 27-May-22 start time was 1455Z, but rounded up to 1510Z
endTime_VLA2_location2=22151144000 # 31-May-22 end time was 1454Z, but rounded down to 1440Z

radius=15 # In km this is the max ship distance from each VLA that we are interested in 

# Convert the date/time location in the AIS file to a new format 
formattedTime=formatTimeStamp(date_time)

# Calls the findInterestedShips() function for both VLAs at both location
interestedVessels,IMOship,timeStamp,cpa,SOGShip,Length,Width,Draft=findInterestedShips(lat_VLA1_location1,long_VLA1_location1,startTime_VLA1_location1,endTime_VLA1_location1,radius,formattedTime)
VLA1_location1=pd.DataFrame({'Vessel':interestedVessels,'IMO':IMOship,'Time Stamp':timeStamp,'CPA':cpa,'SOG (knots)':SOGShip,'Length':Length,'Width':Width,'Draft':Draft})

del interestedVessels,IMOship,timeStamp,cpa,SOGShip,Length,Width,Draft
interestedVessels,IMOship,timeStamp,cpa,SOGShip,Length,Width,Draft=findInterestedShips(lat_VLA2_location1,long_VLA2_location1,startTime_VLA2_location1,endTime_VLA2_location1,radius,formattedTime)
VLA2_location1=pd.DataFrame({'Vessel':interestedVessels,'IMO':IMOship,'Time Stamp':timeStamp,'CPA':cpa,'SOG (knots)':SOGShip,'Length':Length,'Width':Width,'Draft':Draft})

del interestedVessels,IMOship,timeStamp,cpa,SOGShip,Length,Width,Draft
interestedVessels,IMOship,timeStamp,cpa,SOGShip,Length,Width,Draft=findInterestedShips(lat_VLA1_location2,long_VLA1_location2,startTime_VLA1_location2,endTime_VLA1_location2,radius,formattedTime)
VLA1_location2=pd.DataFrame({'Vessel':interestedVessels,'IMO':IMOship,'Time Stamp':timeStamp,'CPA':cpa,'SOG (knots)':SOGShip,'Length':Length,'Width':Width,'Draft':Draft})

del interestedVessels,IMOship,timeStamp,cpa,SOGShip,Length,Width,Draft
interestedVessels,IMOship,timeStamp,cpa,SOGShip,Length,Width,Draft=findInterestedShips(lat_VLA2_location2,long_VLA2_location2,startTime_VLA2_location2,endTime_VLA2_location2,radius,formattedTime)
VLA2_location2=pd.DataFrame({'Vessel':interestedVessels,'IMO':IMOship,'Time Stamp':timeStamp,'CPA':cpa,'SOG (knots)':SOGShip,'Length':Length,'Width':Width,'Draft':Draft})

# Creates a .xlsx file with the list of interested ships for each VLA and location separated onto different sheets and labeled accordingly 
with pd.ExcelWriter(r'C:\Users\underwater\Desktop\UW_Research\AIS_Processing_Alex\SBEX_2022_AIS_InterestedShips.xlsx') as writer:
    VLA1_location1.to_excel(writer,sheet_name='VLA1_Location1',index=False)
    VLA2_location1.to_excel(writer,sheet_name='VLA2_Location1',index=False)
    VLA1_location2.to_excel(writer,sheet_name='VLA1_Location2',index=False)
    VLA2_location2.to_excel(writer,sheet_name='VLA2_Location2',index=False)
 



       