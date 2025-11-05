import os
import shutil
import time
import numpy as np
import rasterio as rio
from datetime import datetime, timedelta
from ToyCell import Cell
from toy_utils import save_init_info, add_info, initial_data_plots, arrival_plots, plot_time_series, plot_pred_vs_model_growth

class CellAuto:
    def __init__(self, 
                 path_to_data='data', 
                 fireID = 'TOY',
                 out_dir = '',
                 num_plots = 0,
                 has_barriers = False,
                 max_hours = 36, 
                 res = 5,
                 burn_out_time = 100
                ):
        
        self.FIRE_ID = fireID
        self.IN_DIR = f'{path_to_data}/raster_inputs_daily/{self.FIRE_ID}/'   #directory from which to pull data
        self.OUT_DIR = f'{path_to_data}/model_outputs/{self.FIRE_ID}/{out_dir}'   #directory to save data in
        if os.path.exists(self.OUT_DIR):
            shutil.rmtree(self.OUT_DIR)
        os.makedirs(self.OUT_DIR)
        self.NUM_PLOTS = num_plots   #number of arrival plots to make
        self.HAS_BARRIERS = has_barriers   #whether or not there is burn_noburn data available
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
        '''
        initial_data_plots(self)
        '''


    def _get_data_file_name(self, data_name):
        return f'{self.IN_DIR}/{self.FIRE_ID}_{data_name}_{self.RESOLUTION}m.tif'

    def _process_data(self):
        # ---- LOAD DATA ----
        # key = class_name, value = file_name
        DATA_DICT = {
            'fire_growth': 'firegrowth_ha'
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

        # create an arrival map to keep track of predictions
        self.arrival_map = np.full((self.NUM_ROWS,self.NUM_COLS), np.nan)
        
        # determine the number of hours/days for which we have data
        #self.TOTAL_HOURS = min(self.fire_growth.shape[0]*24 ,self.wdir_deg.shape[0])
        self.TOTAL_HOURS =  min(self.fire_growth.shape[0]*24 ,self.MAX_HOURS)
        self.TOTAL_DAYS = self.TOTAL_HOURS//24

        tot_hours = min(self.MAX_HOURS, self.TOTAL_HOURS)

        # determine plotting frequency
        if self.NUM_PLOTS > 0: 
            freq = tot_hours//self.NUM_PLOTS
            self.PLOT_HOURS = np.arange(freq,tot_hours+1,freq)
        else:
            self.PLOT_HOURS = np.array([])
                    
        # get ignition date/time info
        self.START_DATE_TIME = datetime.now()
        self.hour_of_day = self.START_DATE_TIME.hour
        
        # determine ignition point(s)
        self.IGN_LOC = set([(self.NUM_ROWS//2, self.NUM_COLS//2)])

        # ---- WIND ----
        # NO EDITS
        self.w_mag = np.full((self.TOTAL_HOURS, self.NUM_ROWS,self.NUM_COLS),30)  
        self.w_dir_rads = np.full((self.TOTAL_HOURS, self.NUM_ROWS,self.NUM_COLS),np.pi)

        # ---- SLOPE ----
        # NO EDITS
        self.s_mag = np.full((self.NUM_ROWS,self.NUM_COLS),0) 
        # aspect gives angle of downslope, want angle of upslope
        self.s_dir_rads = np.full_like(self.fire_growth[0],np.nan)

        # ---- GROWTH ----
        # fire growth conversion 
        # hectares -> m^2: *10000 
        # cancel area unit: /res^2
        # shape: num_days X num_rows X num_cols
        self.fire_growth[self.fire_growth<0] = np.nan
        self.og_growth = self.fire_growth
        growth_adjust = self.fire_growth*10000/self.RESOLUTION**2
        self.fire_growth = np.full((self.TOTAL_DAYS,self.NUM_ROWS,self.NUM_COLS),1)
        #self.fire_growth = 0.5*np.sqrt(self.fire_growth)
        '''
        days = np.arange(1,len(self.fire_growth)+1)
        self.fire_growth = self.fire_growth*np.sqrt(days)[:, np.newaxis, np.newaxis]
        self.fire_growth[1:,0,0] = self.fire_growth[1:,0,0] - self.fire_growth[:-1,0,0]
        
        #self.fire_growth = np.full_like(self.fire_growth, 18)
        '''

        # now get hourly distribution
        # visualize: https://www.desmos.com/calculator/zg5nxw5eza
        hours = (np.arange(24) + self.hour_of_day) % 24
        hourly_weights = 1/24 *(np.sin(np.pi / 12 * hours - 1) + 1)
        hourly_weights = np.full_like(hours,1/24, dtype='float')
        hourly_growth = self.fire_growth[:, np.newaxis, :, :]
        hourly_weights = hourly_weights[np.newaxis, :, np.newaxis, np.newaxis]
        self.hourly_growth = (hourly_growth*hourly_weights).reshape(self.fire_growth.shape[0]*24, self.NUM_ROWS, self.NUM_COLS)
        
        # NO EDITS
        self.growth_mag = self.hourly_growth
        
        # ---- FIRST CELL ---
        # create the first cell, add to burning, and indicate on arrival map 
        for row, col in self.IGN_LOC:
            # new cell: Cell(row, col, raster_shape, ign_loc, 
            #                wspd, w_ang, 
            #                slope, s_ang, 
            #                growth, res, start_time, burn_out_time)
            newCell = Cell(row, col, [self.NUM_ROWS, self.NUM_COLS], [0,0],
                           self.w_mag[self.int_hour][row][col], self.w_dir_rads[self.int_hour][row][col],
                           self.s_mag[row][col], self.s_dir_rads[row][col],
                           self.growth_mag[self.int_hour][row][col], 
                           self.RESOLUTION, 0, self.BURN_OUT_TIME)
            self.burning_cells.add(newCell)
            self.burning_coords.add((row,col))
            self.arrival_map[row][col] = 0

    
    def _update_time(self):
        # update the current hour and day
        self.int_hour = int(self.curr_hour)
        self.int_day = int(self.curr_hour//24)
        curr_secs = self.curr_hour*60
        time_delta = timedelta(seconds=curr_secs)
        curr_date_time = self.START_DATE_TIME + time_delta
        self.hour_of_day = curr_date_time.hour

        # plot if int_hour in PLOT_HOURS
        if self.int_hour in self.PLOT_HOURS:
            arrival_plots(self)

    def _update_growth(self):
        # update growth in all cells
        n_burn = len(self.burning_cells)
        for cell in self.burning_cells:
            hour = self.int_hour
            row = cell.ROW
            col = cell.COL
            #cell._update_cell_growth(wspd, w_ang, slope, s_ang, growth)
            cell.update_cell_growth(self.w_mag[hour][row][col], self.w_dir_rads[hour][row][col], 
                                     self.s_mag[row][col], self.s_dir_rads[row][col], self.growth_mag[hour][row][col])
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
                # otherwise burn the cell
                cell.burn(self.time_step)
                if cell.min_time < 0.001:
                    #if fire is close enough to another ignition point, then ignite a new cell
                    next_cells, next_ign_locs = cell.get_next_burn()
                    new_cells += next_cells
                    new_ign_locs += next_ign_locs
        
        for ext_cell in cells_to_extinguish:
            self.extinguish(ext_cell)

        n_burn = len(self.burning_cells)
        for i in range(len(new_cells)):
            # if cell has already been ignited, don't make a new one
            if (new_cells[i] in self.burning_coords) or (new_cells[i] in self.extinguished_coords):
                continue 

            # otherwise, make a new burning cell
            row = new_cells[i][0]
            col = new_cells[i][1]
            ign_loc = new_ign_locs[i]
            if not self.barriers[row][col]:
                continue # barriers=0 means can't burn there.
            # new cell: Cell(row, col, raster_shape, ign_loc, 
            #                wspd, w_ang, 
            #                slope, s_ang, 
            #                growth, res, start_time, burn_out_time)
            newCell = Cell(row, col, [self.NUM_ROWS, self.NUM_COLS], ign_loc,
                           self.w_mag[self.int_hour][row][col], self.w_dir_rads[self.int_hour][row][col],
                           self.s_mag[row][col], self.s_dir_rads[row][col],
                           self.growth_mag[self.int_hour][row][col], 
                           self.RESOLUTION, self.curr_hour, self.BURN_OUT_TIME)
            self.burning_cells.add(newCell)
            self.burning_coords.add((row,col))
            self.arrival_map[row][col] = self.curr_hour
            

    def run(self):
        start_time = time.time()
        
        while self.curr_hour <= min(self.MAX_HOURS, self.TOTAL_HOURS):
            
            if self.curr_hour >= self.int_hour+1:
                print(np.round(self.curr_hour,2),end=' ')
                self._update_time()
                if self.int_hour >= self.TOTAL_HOURS:
                    break
                self._update_growth()

            # get minimum time required to ignite new cell (or use time step base)
            min_times = np.array([cell.min_time for cell in self.burning_cells])
            time_to_next_hour = self.int_hour+1 - self.curr_hour 
            self.time_step = min(time_to_next_hour, np.nanmin(min_times)) # want to catch changes in underlying data every hour
            if np.isinf(self.time_step):
                print('time step is infinity')
                break
            self.curr_hour += self.time_step

            self.burn_cells()
            
            if len(self.burning_cells) < 1:
                print('empty burning cells set')
                break
        elapsed_time = time.time() - start_time
        add_info(self, f'Time to run simulation: {elapsed_time:.4f} seconds')
        arrival_plots(self)
        '''
        plot_time_series(self)
        '''
        plot_pred_vs_model_growth(self)
        
        




        
        