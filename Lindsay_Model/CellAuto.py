import os
import shutil
import time
import numpy as np
import rasterio as rio
from datetime import datetime, timedelta
from Cell import Cell
from utils import save_init_info, add_info, initial_data_plots, arrival_plots, plot_time_series #plot_daily_time_series, plot_hourly_time_series, plot_pred_vs_model_growth

class CellAuto:
    def __init__(self, 
                 path_to_data, 
                 name,
                 fireID,
                 out_dir = '',
                 num_plots = 0,
                 has_barriers = False,
                 max_hours = 100, 
                 res = 90,
                 burn_out_time = 50
                ):

        self.NAME = name
        self.FIRE_ID = fireID
        self.IN_DIR = f'{path_to_data}/raster_inputs_daily/{self.FIRE_ID}'   #directory from which to pull data
        self.OUT_DIR = f'{path_to_data}/lw_model_outputs/{self.NAME}/{out_dir}'   #directory to save data in. 
        # SUGGESTION: create a new output folder in place of lw_model_outputs and change in path above
        if os.path.exists(self.OUT_DIR):
            shutil.rmtree(self.OUT_DIR)
        os.makedirs(self.OUT_DIR)
        self.NUM_PLOTS = num_plots   #number of arrival plots to make
        self.HAS_BARRIERS = has_barriers   #whether or not there is barrier data available
        self.MAX_HOURS = max_hours   #max num of hours simulation is allowed to run
        self.RESOLUTION = res   #cell length
        self.BURN_OUT_TIME = burn_out_time
        self.burning_cells = set()   #set of burning cells
        self.burning_coords = set()   #set of coords of burning cells
        self.extinguished_cells = set()   #set of extinguished cells
        self.extinguished_coords = set()   #set of coords of extinguished cells
        self.curr_hour = 0   #time since simulation began (# hours)
        self.int_day = 0   #need integer days for indexing fire growth
        self.int_hour = 0   #need integer hours for indexing wind
        self.TIME_STEP_BASE = 1 #largest time step allowed to catch changes in underlying data
        self.time_step = 1 #measured in hours
        self._process_data()
        # plot initial data and save sim info
        save_init_info(self)
        
        initial_data_plots(self)
        


    def _get_data_file_name(self, data_name):
        return f'{self.IN_DIR}/{self.FIRE_ID}_{data_name}_{self.RESOLUTION}m.tif'

    def _process_data(self):
        # ---- LOAD DATA ----
        # key = class_name, value = file_name
        DATA_DICT = {
            'fire_growth': 'firegrowth_ha',
            'slope_rads' : 'slope_radians',
            'aspect' : 'aspect_radians',
            'hillshade' : 'hillshade',
            'wdir_deg' : 'winddir_deg',
            'wspd' : 'windspeed',
            'burnday' : 'VIIRS_burnday'
        }

        if self.HAS_BARRIERS:
            DATA_DICT['barriers'] = 'burn_noburn'
        
        CONSTANT_DATA = ['slope_rads','aspect','hillshade','burnday','barriers']
        
        # load data
        for key, value in DATA_DICT.items():
            file_name = self._get_data_file_name(value)
            with rio.open(file_name) as src:
                if key in CONSTANT_DATA:
                    setattr(self, key, src.read()[0])
                else:
                    setattr(self, key, src.read())

        # ---- EXTRACT INFO ----
        # determine shape of rasters
        self.NUM_ROWS = self.fire_growth.shape[1]
        self.NUM_COLS = self.fire_growth.shape[2]

        # if there was not any barrier data, allow the fire to burn everywhere
        if not self.HAS_BARRIERS:
            self.barriers = np.ones((self.NUM_ROWS,self.NUM_COLS))
        
        # determine the number of hours/days for which we have data
        self.TOTAL_HOURS = min(self.wdir_deg.shape[0],self.fire_growth.shape[0]*24)
        self.TOTAL_DAYS = self.TOTAL_HOURS//24

        tot_hours = min(self.MAX_HOURS, self.TOTAL_HOURS)
        self.num_burning = np.zeros(tot_hours+1)

        # create an arrival map to keep track of predictions
        self.arrival_map = np.full((tot_hours,self.NUM_ROWS,self.NUM_COLS), np.nan)

        # determine plotting frequency
        if self.NUM_PLOTS > 0: 
            freq = tot_hours//self.NUM_PLOTS
            self.PLOT_HOURS = np.arange(freq,tot_hours+1,freq)
        else:
            self.PLOT_HOURS = np.array([])
                    
        # from VIIRS data, get ignition date/time info
        ign_day_time_of_year = np.min(self.burnday[~np.isnan(self.burnday)])
        ign_year = self.FIRE_ID[-8:-4]
        ign_date = datetime.strptime(f"{ign_year}-{int(ign_day_time_of_year)}", "%Y-%j")
        ign_time = ign_day_time_of_year - int(ign_day_time_of_year)
        ign_secs = float(ign_time*24*60*60)
        time_delta = timedelta(seconds=ign_secs)
        self.START_DATE_TIME = ign_date + time_delta
        self.hour_of_day = self.START_DATE_TIME.hour

        # Edit VIIRS times to hours from ignition
        self.burnday -= ign_day_time_of_year
        self.burnhour = self.burnday*24
        
        # determine ignition point(s)
        self.IGN_LOC = set([(int(x[0]), int(x[1])) for x in np.argwhere(self.burnday == 0)])

        # ---- WIND ----
        # shape: num_hours X num_rows X num_cols
        self.w_mag = self.wspd # I believe units are meters/second 
        self.w_dir_rads = self.wdir_deg * np.pi / 180

        # ---- SLOPE ----
        # shape: num_rows X num_cols
        self.s_mag = self.slope_rads
        # aspect gives angle of downslope, want angle of upslope
        self.s_dir_rads = np.remainder(self.aspect+np.pi,2*np.pi)

        # ---- GROWTH ----
        # fire growth conversion 
        # hectares -> m^2: *10000
        # cancel area unit: /res^2
        # shape: num_days X num_rows X num_cols
        self.fire_growth[self.fire_growth<0] = np.nan
        self.fire_growth = self.fire_growth*10000/self.RESOLUTION**2 

        # now get hourly distribution
        # visualize: https://www.desmos.com/calculator/zg5nxw5eza
        hours = (np.arange(24) + self.hour_of_day) % 24
        hourly_weights = 1/24 *(np.sin(np.pi / 12 * hours - 1) + 1)
        hourly_growth = self.fire_growth[:, np.newaxis, :, :]
        hourly_weights = hourly_weights[np.newaxis, :, np.newaxis, np.newaxis]
        self.hourly_growth = (hourly_growth*hourly_weights).reshape(self.fire_growth.shape[0]*24, self.NUM_ROWS, self.NUM_COLS)
        
        # shape: num_hours X num_rows X num_cols
        self.growth_mag = self.hourly_growth
        
        # ---- FIRST CELL ---
        # create the first cell, add to burning, and indicate on arrival map 
        n = len(self.IGN_LOC)
        for row, col in self.IGN_LOC:
            self.burning_coords.add((row,col))
            self._add_burning_cell(row, col, [0,0], n)

            

    def _add_burning_cell(self, row, col, ign_loc, n):
        if self.barriers[row][col]: # barriers=0 means can't burn there. barriers=1 means we can burn
            # new cell: Cell(row, col, raster_shape, ign_loc, 
            #                wspd, w_ang, 
            #                slope, s_ang, 
            #                growth, 
            #                res, start_time, burn_out_time)
            newCell = Cell(row, col, [self.NUM_ROWS, self.NUM_COLS], ign_loc,
                           self.w_mag[self.int_hour][row][col], self.w_dir_rads[self.int_hour][row][col],
                           self.s_mag[row][col], self.s_dir_rads[row][col],
                           # divide the growth prediction by the number of cells contributing to that growth
                           self.growth_mag[self.int_hour][row][col]/n, 
                           self.RESOLUTION, self.curr_hour, self.BURN_OUT_TIME)
            self.burning_cells.add(newCell)
    
    # ---- DAILY PERIMETER UPDATE OPTION ---
    def _update_perimeter(self):
        # empty all cell sets
        self.burning_cells = set()  
        self.burning_coords = set()  
        self.extinguished_cells = set()   
        self.extinguished_coords = set()
        # limit the VIIRS data to only cells that have burned prior to current time
        current_VIIRS = self.burnhour.copy()
        current_VIIRS[self.burnhour >= self.int_hour] = np.nan
        rows, cols = current_VIIRS.shape
        # pad with nans for window extraction
        current_VIIRS = np.pad(current_VIIRS, pad_width=1, mode='constant', constant_values=np.nan)

        # loop through all the original cells
        for i in range(1, rows+1):
            for j in range(1, cols+1):
                center_val = current_VIIRS[i, j]
                # consider all the cells that are burning
                if not np.isnan(center_val):
                    # Extract 3x3 window around center
                    neighborhood = current_VIIRS[i-1:i+2, j-1:j+2].copy()
                    if np.isnan(neighborhood).any():
                        # if any of the neighbors are not burning, add the point to the perimeter (burning cells)
                        self.burning_coords.add((i-1, j-1))  # offset for padding
                    else:
                        # otherwise, all cells around it are burnt, so consider the cell to be extinguished
                        self.extinguished_coords.add((i-1, j-1))
        n = len(self.burning_coords)
        for row, col in self.burning_coords:
            self._add_burning_cell(row, col, [0,0], n)
        
    
    def _update_time(self):
        # update the current hour and day
        self.int_hour = int(self.curr_hour)
        self.int_day = int(self.curr_hour//24)
        if self.int_hour % 24 == 23:
            arrival_plots(self)
        '''
        # use to update the perimeter to VIIRS data at the end of each day
        if self.int_hour % 24 == 0:
            self._update_perimeter()
        '''
        curr_secs = self.curr_hour*60
        time_delta = timedelta(seconds=curr_secs)
        curr_date_time = self.START_DATE_TIME + time_delta
        self.hour_of_day = curr_date_time.hour

        self.num_burning[self.int_hour] = len(self.burning_cells)
    
        # plot if int_hour in PLOT_HOURS
        if self.int_hour in self.PLOT_HOURS:
            arrival_plots(self)

    def _update_growth(self):
        # update growth in all cells
        n = max(1, len(self.burning_cells))
        for cell in self.burning_cells:
            hour = self.int_hour
            row = cell.ROW
            col = cell.COL
            # cell.update_cell_growth(wspd, w_ang, 
            #                         slope, s_ang, growth)
            cell.update_cell_growth(self.w_mag[hour][row][col], self.w_dir_rads[hour][row][col], 
                                     self.s_mag[row][col], self.s_dir_rads[row][col], self.growth_mag[hour][row][col]/n)

    def _burnt_neighbors(self, cell):
        # check if all neighbors of a cell are burnt
        for n in cell.NEIGHBORS.values():
            if (n not in self.burning_coords) and (n not in self.extinguished_coords):
                # if one cell is not burning and not extinguished we can burn it
                return False
        # all neighbors were either burning or extinguished
        return True

    def extinguish(self, cell):
        cell.is_extinguished = True
        self.burning_cells.remove(cell)
        self.burning_coords.remove((cell.ROW,cell.COL))
        self.extinguished_cells.add(cell)
        self.extinguished_coords.add((cell.ROW,cell.COL))
    
    def burn_cells(self):
        new_cells, new_ign_locs = [], []
        cells_to_extinguish = set()
        
        for cell in self.burning_cells:
            if cell.is_extinguished:
                # move to extinguished if burning too long
                cells_to_extinguish.add(cell)
            elif self._burnt_neighbors(cell):
                # extinguish cell if has no where to burn
                cell.burnt_neighbors = True
                cells_to_extinguish.add(cell)
            else:
                # otherwise progress the fire within the cell
                cell.burn(self.time_step)
                if cell.min_time < 0.001:
                    # if fire is close enough to another ignition point, then ignite a new cell
                    next_cells, next_ign_locs = cell.get_next_burn()
                    new_cells += next_cells
                    new_ign_locs += next_ign_locs
        
        for ext_cell in cells_to_extinguish:
            self.extinguish(ext_cell)

        n = max(1, len(self.burning_cells))
        for i in range(len(new_cells)):
            # if cell has already been ignited, don't make a new one
            if (new_cells[i] in self.burning_coords) or (new_cells[i] in self.extinguished_coords):
                continue 

            # otherwise, make a new burning cell
            row = new_cells[i][0]
            col = new_cells[i][1]
            ign_loc = new_ign_locs[i]
            self.burning_coords.add((row,col))
            self._add_burning_cell(row, col, ign_loc, n)
            self.arrival_map[self.int_day][row][col] = self.curr_hour

    def run(self):
        # time the run
        start_time = time.time()

        # run the simulation until we hit the max time or run out of data
        while self.curr_hour <= min(self.MAX_HOURS, self.TOTAL_HOURS):
            # every hour we must update the wind and growth data
            if self.curr_hour >= self.int_hour+1:
                print(np.round(self.curr_hour,2),end=' ')
                self._update_time()
                if self.int_hour >= self.TOTAL_HOURS:
                    break
                self._update_growth()

            # get minimum time required to ignite new cell
            min_times = np.array([cell.min_time for cell in self.burning_cells])
            # want to catch changes in underlying data every hour
            time_to_next_hour = self.int_hour+1 - self.curr_hour 
            # progress time forward the min amount of these two times
            self.time_step = min(time_to_next_hour, np.nanmin(min_times)) 
            if np.isinf(self.time_step):
                print('time step is infinity')
                break
            self.curr_hour += self.time_step

            self.burn_cells()
            
            if len(self.burning_cells) < 1:
                print('empty burning cells set')
                break
        # run time logging
        elapsed_time = time.time() - start_time
        add_info(self, f'Time to run simulation: {elapsed_time:.4f} seconds')
        plot_time_series(self)

        




        
        