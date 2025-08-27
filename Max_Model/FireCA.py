from Cell import Cell
from SimpleCell import SimpleCell
from Stats import *
import numpy as np
#import pandas as pd
import rasterio as rio
#import os, glob
import datetime
#import subprocess
import re
#import shutil
#import collections
from itertools import chain
from matplotlib import colors
from matplotlib import pyplot as plt

# assumptions:
# if there is a barrier raster the ignition is not on a barrier

class fireCA:

    def __init__(self, path_to_data=False):
        self.plotting = False
        self.show_plot = False
        self.hasBarriers = False
        self.isSimpleCell = False
    
        self.consts = {'fireID':None,
                       'growthWeight':1,
                       'maxTimeStep':150,
                       'n_burnt':False, #stop model if hits a limit (1000) of cells burnt
                       'spatial_resolution':250,
                       's_base':2**0, #timestep adjustment, right now = 1
                       'weight':0.32,
                       'weight_d':1.6,
                       'weight_phi':0.7,
                       'ign_loc':[(93,50)],
                       'fname':'',
                       'start_time':datetime.datetime.strptime('2017-07-24','%Y-%m-%d')}
        
        if path_to_data:
            self.path_to_data = path_to_data
    
    def get_data(self):
        print("Getting data...")
        path_to_data = self.path_to_data

        with rio.open(self.path_to_data+self.consts['fireID']+'/'+self.consts['fireID']+'_firegrowth_ha_'+str(self.consts['spatial_resolution'])+'m.tif') as src:
            self.RoS_cube = src.read()
            self.crs = src.crs
            self.transform = src.transform

        # topography
        self.slope = rio.open(path_to_data+self.consts['fireID']+'/'+self.consts['fireID']+'_slope_radians_'+str(self.consts['spatial_resolution'])+'m.tif').read()[0]
        self.aspect = rio.open(path_to_data+self.consts['fireID']+'/'+self.consts['fireID']+'_aspect_radians_'+str(self.consts['spatial_resolution'])+'m.tif').read()[0]
        self.hillshade = rio.open(path_to_data+self.consts['fireID']+'/'+self.consts['fireID']+'_hillshade_'+str(self.consts['spatial_resolution'])+'m.tif').read()[0]
        
        # wind
        self.wfrom = rio.open(path_to_data+self.consts['fireID']+'/'+self.consts['fireID']+'_winddir_deg_'+str(self.consts['spatial_resolution'])+'m.tif').read()
        self.wspd = rio.open(path_to_data+self.consts['fireID']+'/'+self.consts['fireID']+'_windspeed_'+str(self.consts['spatial_resolution']) + 'm.tif').read()
                        
        # burnday
        self.burnday = rio.open(path_to_data+self.consts['fireID']+'/'+self.consts['fireID']+'_VIIRS_burnday_'+str(self.consts['spatial_resolution'])+'m.tif').read()[0]
        
        # start date
        self.consts['jday'] = np.min(self.burnday[~np.isnan(self.burnday)]) #gives day of year first point is burned

        if self.hasBarriers:
            self.barriers = rio.open(path_to_data+self.consts['fireID']+'/'+self.consts['fireID']+"_burn_noburn_" + str(self.consts['spatial_resolution'])+'m.tif').read()[0]
        else:
            self.barriers = np.ones_like(self.RoS_cube[0])
        
        try:
            self.consts['start_time'] = datetime.datetime.strptime(f"{self.consts['fireID'][-8:-4]}-{int(self.consts['jday'])}", "%Y-%j")
            #old datetime.datetime(int(re.findall(r'\d{3,}',path_to_data)[0][-8:-4]), 1, 1) + datetime.timedelta(int(self.consts['jday'] - 1)) 
            #should use datetime.strptime(f"{fireID[-8:-4]}-{int(jday)}", "%Y-%j")
        except:
            self.consts['start_time'] = datetime.datetime.strptime('2017-07-23', '%Y-%m-%d')
            print("No start time found in directory name. Start time set to default")
        
    def preprocess(self):
        print("Starting Preprocessing...")
        t,self.m,self.n = self.RoS_cube.shape
        #t*24=# hours, pick smallest amount of data we have
        #now t is in hours
        self.t = min(self.wfrom.shape[0], t*24)

        #convert from hectacres/day to m/hour
        self.RoS_cube = ((self.RoS_cube*10000)/24)*self.consts['growthWeight'] 
        self.RoS_cube[self.RoS_cube<0] = np.nan

        # arrival map is an array that holds which cells were ignited at which times
        self.arrivalMap = np.full(self.RoS_cube[0].shape, np.nan)
        
        self.VIIRS_burnt_number = np.zeros(self.t)
        for i in range(self.t):
            j = i + self.consts['jday']*24
            self.VIIRS_burnt_number[i] = np.sum(self.burnday*24 <= j)
        
        self.upslope = ((((3*np.pi)/2) - self.aspect)+ np.pi) % (2*np.pi) - np.pi #???
        
        self.wfrom[self.wfrom<0] = np.nan
        self.wdir = (self.wfrom-90)*(np.pi/180) #horizontal flip for wind direction, rotation to change to math coords
        self.wspd[self.wspd<0] = 0
        
        # maximum rate and direction of spread
        # some bullshit i made up, just graph it
        # based off some article i read about the 10% rule
        wm = (self.wspd*3.6) # m s-1 -> km hr-1
        wm = np.exp(0.05*wm)/(1+np.exp(0.6*(10-wm))) # in km hr-1
        # slope conversion from radians to degrees
        self.slope = self.slope * (180/np.pi)
        # slope from Butler, Anderson, Catchpole 2007
        sm = 0.2*np.exp(0.13*self.slope) + 0.6

        # TODO: figure out what is going on with k
        k = np.sqrt(self.RoS_cube)
        self.r_max_cube = np.sqrt(np.sqrt(np.square(wm) + np.square(sm) + 2*wm*sm*np.cos(self.wdir-self.upslope)))*np.repeat(k,24,axis=0)[:wm.shape[0]]
        self.phi_max_cube = (self.upslope + np.arctan2(wm*np.sin(self.wdir-self.upslope), sm + wm*np.cos(self.wdir-self.upslope)))
        
        self.consts['ign_loc'] = set([(x[1], self.m-1-x[0]) for x in np.argwhere(self.burnday == self.consts['jday'])])
        print("Preprocessing Complete")
        self.reset()
        
    def reset(self):
        print("Resetting...")
        self.graph_burning = self.consts['ign_loc']
        self.extinguished = set()
        self.predicted = set()
                
        self.CA_burnt_number = np.zeros(self.t)
        self.FP = np.zeros(self.t)
        self.FN = np.zeros(self.t)
        self.TP = np.zeros(self.t)
        self.TN = np.zeros(self.t)
        self.accuracy = np.zeros(self.t)
        self.specificity = np.zeros(self.t)
        self.precision = np.zeros(self.t)
        self.recall = np.zeros(self.t)
        self.kappa = np.zeros(self.t)
        self.quantityDisagreement = np.zeros(self.t)
        self.allocationDisagreement = np.zeros(self.t)
        self.totalDisagreement = np.zeros(self.t)
        self.n_prev_burnt = 0
        self.i_prev = 0
        
        self.phi_max = self.phi_max_cube[0]
        self.r_max = self.r_max_cube[0]
        self.RoS = self.RoS_cube[0]
        self.burning = np.array([])
        for cell in self.consts['ign_loc']:
            sx = cell[0]
            sy = cell[1]
            col = sx
            row = self.m - sy - 1
            if self.isSimpleCell:
                newCell = SimpleCell('c', sx, sy, self.RoS[row][col], self.consts['spatial_resolution'])
            else:
                newCell = Cell('c', sx, sy, self.RoS[row][col], self.consts['spatial_resolution'])
            v = [self.r_max[row][col], self.phi_max[row][col]]
            newCell.calcRoS(v, 0.0001)
            self.burning = np.append(self.burning, newCell)
            self.predicted.add((row, col))
            self.arrivalMap[row][col] = 0
        print("Resetting Complete")
    
    def plot(self, path=False):
        state_mtx = np.ones_like(self.RoS).astype(int)*5

        for pixel in self.graph_burning:
            truePixel = (self.m-1 - pixel[1], pixel[0])
            state_mtx[truePixel] = 4
        for pixel in self.extinguished:
            truePixel = (self.m-1 - pixel[1], pixel[0])
            state_mtx[truePixel] = 2
        
        relevant_pixels = list(chain(self.graph_burning, self.extinguished))
        x = [pixel[0] for pixel in relevant_pixels]
        y = [self.m-1-pixel[1] for pixel in relevant_pixels]
        top = np.min(y)
        bottom = np.max(y)
        left = np.min(x)
        right = np.max(x)
        
        xpad = max((right-left)/3, 20)
        ypad = max((bottom-top)/3, 20)
        top = int(max(0, top-ypad))
        bottom = int(min(self.m, bottom+ypad))
        left = int(max(0, left-xpad))
        right = int(min(self.n, right+xpad))
        
        j = int(self.i//self.consts['s_base'])
        
        cmap = colors.ListedColormap([(0,0,0,0), 'black', 'dimgray', 'red', 'gold', (0,0,0,0)])
        norm = colors.BoundaryNorm([0,1,2,3,4,5,6], cmap.N)
        '''
        fig,ax = plt.subplots(2,2,figsize=(15,7))
        ax[0,0].imshow(self.hillshade[top:bottom,left:right], alpha=1, cmap='gray', interpolation='none', zorder=0)
        ax[0,0].imshow(state_mtx[top:bottom,left:right], cmap=cmap, alpha=0.6, norm=norm, interpolation='none', zorder=10)
        viirs = np.argwhere(self.burnday[top:bottom,left:right]<=self.consts['jday']+(j/24)).T
        ax[0,0].scatter(viirs[1], viirs[0], c='r', s=1, alpha=0.4, zorder=5)
        ax[0,0].tick_params(left=False, right=False, top=False, bottom=False, labelleft=False, labelright=False, labeltop=False, labelbottom=False)
        ax[0,1].plot(self.CA_burnt_number[:j], label='CA')
        ax[0,1].plot(self.VIIRS_burnt_number[:j], label='VIIRS')
        ax[0,1].plot(self.FP[:j], ':r', label='false positive')
        ax[0,1].plot(self.FN[:j], ':b', label='false negative')
        ax[0,1].plot(self.TP[:j], ':g', label='true positive')
        ax[0,1].set_ylabel('# cells')
        ax[0,1].legend()
        ax[1,0].plot(self.accuracy[:j], 'r', label='accuracy')
        ax[1,0].plot(self.specificity[:j], 'b', label='specificity')
        ax[1,0].plot(self.precision[:j], 'g', label='precision')
        ax[1,0].plot(self.recall[:j], 'o', label='recall')
        ax[1,0].plot(self.kappa[:j], '-.g', label='kappa')
        ax[1,0].legend()
        ax[1,1].plot(self.quantityDisagreement[:j], 'r', label='quatity disagreement')
        ax[1,1].plot(self.allocationDisagreement[:j], 'b', label='allocation disagreement')
        ax[1,1].plot(self.totalDisagreement[:j], 'g', label='total disagreement')
        ax[1,1].legend()
        '''
        
        # TODO: remove
        fig, ax = plt.subplots(2, figsize=(15,7))
        ax[0].imshow(self.hillshade[top:bottom,left:right], alpha=1, cmap='gray', interpolation='none', zorder=0)
        ax[0].imshow(state_mtx[top:bottom,left:right], cmap=cmap, alpha=0.6, norm=norm, interpolation='none', zorder=10)
        viirs = np.argwhere(self.burnday[top:bottom,left:right]<=self.consts['jday']+(j/24)).T
        ax[0].scatter(viirs[1], viirs[0], c='r', s=1, alpha=0.4, zorder=5)
        ax[0].tick_params(left=False, right=False, top=False, bottom=False, labelleft=False, labelright=False, labeltop=False, labelbottom=False)
        im1 = ax[1].imshow(self.arrivalMap[top:bottom, left:right], interpolation='none')
        fig.colorbar(im1, ax=ax[1], orientation='vertical')
        
        fig.suptitle(f"{str((np.datetime64(self.consts['start_time'],'s')+((self.hour+24*self.day)*60+self.minute)*60+self.second)).replace('T',' ')}\n"\
                  f"area burned: {((self.consts['spatial_resolution']/1000)**2*(len(self.graph_burning))):.2f} km2")
        if path:
            plt.savefig(path)
        if self.show_plot:
            plt.show()
        plt.close()
        self.saveArrivalMap()

    # returns the ignition point for the new cell
    def findNeighbIgnPt(self, ignPt):
        if np.array_equal(ignPt, [0,0]):
            return 'ne'
        elif np.array_equal(ignPt, [0, 0.5]):
            return 'e'
        elif np.array_equal(ignPt, [0,1]):
            return 'se'
        elif np.array_equal(ignPt, [0.5, 1]):
            return 's'
        elif np.array_equal(ignPt, [1, 1]):
            return 'sw'
        elif np.array_equal(ignPt, [1, 0.5]):
            return 'w'
        elif np.array_equal(ignPt, [1, 0]):
            return 'nw'
        elif np.array_equal(ignPt, [0.5, 0]):
            return 'n'
        else:
            return "Error: not a direction"
    
    def burn(self):
        # progress fire in all burning cells
        for cell in self.burning:
            newBurning = cell.updatePhi(self.s)
            for nb in newBurning:
                trueDir = (nb*2)-1
                # the x and y values of the cell
                tx = int(trueDir[0] + cell.x)
                ty = int(trueDir[1] + cell.y)
                # the row and column where the data for the cell is located 
                col = tx
                row = self.m-1-ty
                if 0<tx<self.n and 0<ty<self.m and self.barriers[row][col] == 1:
                    if self.isSimpleCell:
                        newCell = SimpleCell(self.findNeighbIgnPt(nb), tx, ty, self.RoS[row][col], self.consts['spatial_resolution'])
                    else:
                        newCell = Cell(self.findNeighbIgnPt(nb), tx, ty, self.RoS[row][col], self.consts['spatial_resolution'])
                    if (newCell not in self.burning) and ((newCell.x, newCell.y) not in self.extinguished): #should change to (tx,ty) not in graph_burning and ?? not in extinguished ??
                        v = [self.r_max[row][col], self.phi_max[row][col]]
                        newCell.calcRoS(v, 0.0001)
                        self.burning = np.append(self.burning, newCell)
                        self.graph_burning.add((tx, ty))
                        self.predicted.add((row, col))
                        self.arrivalMap[row][col] = self.i # + self.s

    def extinguish(self):
        removeFromBurning = np.array([])
        removeFromGraphBurning = set()
        for cell in self.burning:
            # if the cell is surrounded by burning cells remove that cell from burning
            neighbors = [(x,y) for x in range(max(cell.x-1, 0),min(cell.x+2, self.n)) 
                         for y in range(max(cell.y-1, 0),min(cell.y+2, self.m))]
            if self.allIn(neighbors, self.graph_burning.union(self.extinguished)):
                removeFromBurning = np.append(removeFromBurning, cell)
                removeFromGraphBurning.add((cell.x, cell.y))
                self.extinguished.add((cell.x, cell.y))
        self.burning = np.setdiff1d(self.burning, removeFromBurning)
        self.graph_burning.difference_update(removeFromGraphBurning)

    def saveArrivalMap(self):
        
        metadata = {
            'driver':'GTiff',
            'height':self.arrivalMap.shape[0],
            'width':self.arrivalMap.shape[1],
            'count':1,
            'dtype':self.arrivalMap.dtype,
            'crs':self.crs,
            'transform':self.transform
        }

        with rio.open(self.consts['gif_dir']+'arrivalMap_' + str(self.i).zfill(2) + '.tiff', 'w', **metadata) as dest:
            dest.write(self.arrivalMap, 1)
    
    def allIn(self, list1, list2):
        a = True
        for l in list1:
            if l not in list2:
                a = False
        return a
    
    def run(self):
        # take t and multiply it by 32 (s_base)
        # keep track of how many 32ths have passed
        # begin a while loop
        # update weather conditions
        # determine step size
        # transition smoldering to burnt
        # progress fire in all burning cells
        # if any burning cells reach the threshold, transition them to smoldering
        # ignite susceptible cells
        # increment 32ths
        
        self.i = 0
        self.j = 0
        self.hour = 0
        self.day = 0
        self.s = self.consts['s_base']
        #while self.i < (self.t-1)*self.consts['s_base']:
        while self.i < self.consts['maxTimeStep']:
            
            self.n_burnt = len(self.burning)+len(self.extinguished)

            
            if self.n_burnt > 1000 and self.consts['n_burnt']:
                print("Number of burning cells exceeds 1000")
                self.plot(self.consts['gif_dir']+f"{self.day:02}_{self.hour:02}_{self.minute:02}_{self.second:02}.png")
                break
            
            
            # update weather conditions
            if self.i > self.consts['s_base']*(24*self.day+self.hour):
                self.j = int(self.i//self.consts['s_base'])
                self.CA_burnt_number[self.j] = self.n_burnt
                #VIIRS_pixels = set([tuple(x) for x in np.argwhere(24*self.burnday<=self.consts['jday']*24+self.j)])
                VIIRS_pixels = set([tuple(x) for x in np.argwhere(self.burnday<=self.consts['jday']+(self.j/24))])
                CA_pixels = self.graph_burning.union(self.extinguished)
                self.FP[self.j] = calcFalsePositive(self.predicted, VIIRS_pixels)
                self.FN[self.j] = calcFalseNegative(self.predicted, VIIRS_pixels)
                self.TP[self.j] = calcTruePositive(self.predicted, VIIRS_pixels)
                self.accuracy[self.j] = calcAccuracy(self.FP[self.j], self.FN[self.j], self.TP[self.j])
                self.precision[self.j] = calcPrecision(self.FP[self.j], self.TP[self.j])
                self.recall[self.j] = calcRecall(self.FN[self.j], self.TP[self.j])
                self.quantityDisagreement[self.j] = calcQuantityDisagreement(self.predicted, VIIRS_pixels)
                self.allocationDisagreement[self.j] = calcAllocationDisagreement(self.FP[self.j], self.FN[self.j], self.quantityDisagreement[self.j])
                self.totalDisagreement[self.j] = self.quantityDisagreement[self.j] + self.allocationDisagreement[self.j]
                self.kappa[self.j] = calcKappa(self.FP[self.j], self.TP[self.j], self.FN[self.j])
                if self.hour != self.j%24:
                    print(self.i,end=' ')
                    self.hour = self.j%24
                    self.r_max = self.r_max_cube[self.j]
                    self.phi_max = self.phi_max_cube[self.j]
                    if self.hour == 0:
                        self.day = int(self.j//24)
                        self.RoS = self.RoS_cube[self.day]
                    for cell in self.burning:
                        row = self.m-1-cell.y
                        col = cell.x
                        cell.aos = self.RoS[row][col]
                        v = [self.r_max[row][col], self.phi_max[row][col]]
                        cell.calcRoS(v, 0.0001)

            # plotting
            if ((self.plotting) and ((self.n_burnt - self.n_prev_burnt > 0.05*self.n_prev_burnt) or (self.i - self.i_prev > 6*self.consts['s_base']))):
                self.n_prev_burnt = self.n_burnt
                self.i_prev = self.i
                self.minute = int(((self.i/self.consts['s_base'])%1)*60)
                self.second = int(((((self.i/self.consts['s_base'])%1)*60)%1)*60)
                self.plot(self.consts['gif_dir']+f"{self.day:02}_{self.hour:02}_{self.minute:02}_{self.second:02}.png")

            self.s = self.consts['s_base']
            # update step size
            for cell in self.burning:
                traveltimes = cell.getTravelTimes()
                if not len(traveltimes) == 0:
                    tmin = np.min(cell.getTravelTimes())
                    if tmin < self.s and tmin>0.0001: self.s=tmin
            #print('s =',self.s)
            
            # burn cells based on level set
            # this is where we update phi for the neighboring cells
            self.burn()
            self.extinguish()
            
            # increment 32ths
            self.i += self.s
            
    
