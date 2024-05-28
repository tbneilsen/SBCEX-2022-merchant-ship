# -*- coding: utf-8 -*-
"""
Created on Jun 14 2022
Last Modified: Sept 6 2023

@author: Alexandra Hopps McDaniel

Purpose: 
Sorts through a download of AIS data and finds all ships of interest for further analysis

Output: 
Formatted spreadsheet with all interested ships and their pertinent information

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
file=r'C:\Users\alexh\Documents\UW_Acoustics\Codes\sbc-vla-data\processing_AIS_data\AIS_download_2022_Proteus.csv' 
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

size_date=len(Vessel) # Length of each array of data 


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

def findInterestedShips(lat,long,startTime,endTime,radius,TimeStamp):

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
    for i in range(size_date):
        # Filters out ships outside of time frame 
        if int(TimeStamp[i])>startTime and int(TimeStamp[i])<endTime: 
            
            # Use the haversine formula to compute the distance of each ship to the VLAs
            # Can find information on the haversine formula: https://en.wikipedia.org/wiki/Haversine_formula
            # Note: Latitudes and Longitudes must be in radians for the equation to work
            haversine=2*r*np.arcsin(np.sqrt(np.sin(np.radians(lat-latitude[i])/2)**2+np.cos(np.radians(lat))*np.cos(np.radians(latitude[i]))*np.sin(np.radians(long-longitude[i])/2)**2))
                    
            if int(haversine)<radius: # filters out all ships that are farther than the set radius from the VLA
                distance.append(haversine)
                vessel.append(Vessel[i])
                time.append(TimeStamp[i])
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
    
    # This initializes what will be the final list of information regarding each ship, no repeats
    cpa=[]
    timeStamp=[]
    sog=[]
    Length=[]
    Width=[]
    Draft=[]
    IMOship=[]
    # This loop sorts through the filtered list of ships, finds the CPA and time stamp, and saves the 
    # information in an array in the same order as the interested ships array
    for i in range(len(interestedVessels)):
        d=[]
        js=[]
        for j in range(len(distance)):
            # For each vessel in the time ranges, find the distance according to the haversine formula
            if vessel[j]==interestedVessels[i]:
                # Appends these distances to an array and keeps track of the indexes for the entire array of AIS ships
                d.append(distance[j])
                js.append(j)
                
        index=d.index(min(d)) # Finds the index from the complete AIS data array of the CPA
        cpa.append(d[index]) # Appends the cpa associated with that index
        
        # The following use the index from the entire AIS data set to add the ship information from the row at CPA
        timeStamp.append(time[js[index]])
        sog.append(sogShip[js[index]])
        Length.append(lengthShip[js[index]])
        Width.append(widthShip[js[index]])
        Draft.append(draftShip[js[index]])
        IMOship.append(imoShip[js[index]])
        
        
    return interestedVessels,IMOship,timeStamp,cpa,sog,Length,Width,Draft

# Information regarding the time and location of Proteus deployment 
lat=40.459
long=-70.5635
# Use the special format
startTime=22127025421 # 07-May-22, first data file is at time 025421
endTime=22154125500 # 03-June-22, first data file is at time 125500

radius=15 # In km this this the max radius of interested ships from the arrays that we care about
TimeStamp=formatTimeStamp(date_time) # Makes an array of correctly formatted time stamps for AIS data 

# Save to a excel file
interestedVessels,IMOship,timeStamp,cpa,SOGShip,Length,Width,Draft=findInterestedShips(lat,long,startTime,endTime,radius,TimeStamp)
data=pd.DataFrame({'Vessel':interestedVessels,'IMO':IMOship,'Time Stamp':timeStamp,'CPA':cpa,'SOG (knots)':SOGShip,'Length':Length,'Width':Width,'Draft':Draft})


with pd.ExcelWriter('SBEX_2022_InterestedShips_Proteus_jul.xlsx') as writer:
     data.to_excel(writer,sheet_name='Proteus',index=False)
     



    