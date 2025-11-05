import numpy as np

class Cell:
    def __init__(self, row, col, raster_shape, ign_loc, wspd, w_ang, slope, s_ang, growth, res, start_time, burn_out_time):
        self.ROW = row
        self.COL = col
        self.MAX_ROWS = raster_shape[0]
        self.MAX_COLS = raster_shape[1]
        self.IGN_LOC = np.array(ign_loc)   #place where fire started within cell
        self.START_TIME = start_time  #time cell was ignited measured in hours since simulation began
        self.BURN_OUT_TIME = burn_out_time  #length of time cell is allowed to burn before being extinguished (hours)
        self.RESOLUTION = res   #length of cell
        self._set_arrays()
        # project growth onto angles
        self.update_cell_growth(wspd, w_ang, slope, s_ang, growth)
        self.burning_time = 0   #amount of time cell has been burning (hours)
        self.is_extinguished = False
        self.burnt_neighbors = False

    def __str__(self):
        return f"({self.ROW},{self.COL})"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return (self.ROW==other.ROW and self.COL==other.COL)

    def __hash__(self):
        return hash((self.ROW, self.COL))

    def print_cell_info(self):
        print(f"{self}    is_extinguished: {self.is_extinguished}  ")
        print(f"Ign location: {self.IGN_LOC}")
        print(f"min_time:{self.min_time:.2f}   burning_time: {self.burning_time:.2f}   burnt_neighbors: {self.burnt_neighbors}")
        print('angles:', np.round(self.ANGLES,2))
        print('distances:', np.round(self.distances,2))
        print('ros:', np.round(self.ros,2))
        print('times:', np.round(self.times,2))
        
    def _set_arrays(self):
        # directions of neighboring cells
        self.DIRS = np.array(['s','sw','w','nw','n','ne','e','se'])
        
        # ignition points of next cell (if you burn north, you will be on the south side of that cell)
        self.NEXT_IGN = {
            'n': [0,-1],
            'ne': [-1,-1],
            'e': [-1,0],
            'se': [-1,1],
            's': [0,1],
            'sw': [1,1],
            'w': [1,0],
            'nw': [1,-1] 
        }
        # find angles to all other cell ignition points
        ign_pts = np.array([[0,-1],[-1,-1],[-1,0],[-1,1],
                            [0,1],[1,1],[1,0],[1,-1]])
        vectors = ign_pts - self.IGN_LOC
        angles = np.array([np.arctan2(vect[1],vect[0]) for vect in vectors])
        self.ANGLES = (3*np.pi/2-angles)%(2*np.pi)
        
        # calc distances to other cell ign points
        dist = np.array([np.linalg.norm(vect)/2 for vect in vectors])
        dist[dist<=0] = np.nan
        self.distances = dist*self.RESOLUTION
        
        # find coords of neighbors
        self._get_neighbors()
        

    def _get_neighbors(self):
        neighbors = {
            'n': (self.ROW-1, self.COL),
            'ne': (self.ROW-1, self.COL+1),
            'e': (self.ROW, self.COL+1),
            'se': (self.ROW+1, self.COL+1),
            's': (self.ROW+1, self.COL),
            'sw': (self.ROW+1, self.COL-1),
            'w': (self.ROW, self.COL-1),
            'nw': (self.ROW-1, self.COL-1)
        }
        # remove any neighbors off the grid
        keys_to_remove = []
        
        for key, value in neighbors.items():
            if value[0]<0 or value[1]<0 or value[0]>=self.MAX_ROWS or value[1]>=self.MAX_COLS:
                keys_to_remove.append(key)

        self.NEIGHBORS = {key: value for key,value in neighbors.items() if key not in keys_to_remove}
        
        # remove them from burning calculation considerations
        indices = []
        for key in keys_to_remove:
            indices.append(np.argwhere(self.DIRS == key)[0,0])
        self.distances[indices] = np.nan

    
    def _get_wind_factor(self, wspd, w_ang):
        '''
        # CONSTANT 
        return np.full_like(self.ANGLES, 1)
        
        wspd = 10
        w_ang = 5/4*np.pi
        '''
        # visualize: https://www.desmos.com/calculator/bmcsrmhfn3
        w_projs = wspd*np.cos(self.ANGLES - w_ang)
        w_neg_ids = np.where(w_projs<=0)
        w_pos_ids = np.where(w_projs>0)
        '''
        # NO EDITS
        w_projs[w_neg_ids] = 0
        return w_projs
        '''
        # PIECEWISE
        w_pos = 50 / (1 + 49*np.exp(-0.25*w_projs))
        w_neg = 2 / (1 + np.exp(-0.1*w_projs))
        
        w_factor = np.ones_like(w_projs)
        w_factor[w_neg_ids] = w_neg[w_neg_ids]
        w_factor[w_pos_ids] = w_pos[w_pos_ids]
        return w_factor
        
        
        
    def _get_slope_factor(self, slope, s_ang):
        '''
        # CONSTANT 
        return np.full_like(self.ANGLES, 1)
        
        # use slope as ratio instead of angle
        slope = np.tan(slope)
        
        slope = np.pi/12
        s_ang = 5/4*np.pi
        '''
        # visualize: https://www.desmos.com/calculator/iuza8zwfnz
        s_projs = slope*np.cos(self.ANGLES - s_ang)
        s_neg_ids = np.where(s_projs<=0)
        s_pos_ids = np.where(s_projs>0)
        '''
        # NO EDITS
        s_projs[s_neg_ids<=0] = 0
        return s_projs
        '''
        # PIECEWISE
        s_pos = 5 / (1 + 4*np.exp(-10*s_projs))
        s_neg = 2 / (1 + np.exp(-1*s_projs))
    
        s_factor = np.ones_like(s_projs)
        s_factor[s_neg_ids] = s_neg[s_neg_ids]
        s_factor[s_pos_ids] = s_pos[s_pos_ids]
        return s_factor
        
             
    def update_cell_growth(self, wspd, w_ang, slope, s_ang, growth):
        
        w_factor = self._get_wind_factor(wspd, w_ang)
        s_factor = self._get_slope_factor(slope, s_ang)
        g_factor = growth
        
        self.ros = g_factor*w_factor*s_factor
        
        # require positive ros
        self.ros[self.ros<=0] = np.nan
        
        # return ros to meaningful units
        self.ros = self.ros*self.RESOLUTION

        # calculate time required to reach other ign pts
        self._get_times()
        
    def _get_times(self):
        # calc time it takes to get to other cell ign points
        self.times = self.distances/self.ros
        self.times[self.times<0] = np.nan

        # if times contains any values, set the min time required to ignite another cell to the min time
        if not np.all(np.isnan(self.times)):
            self.min_time = np.nanmin(self.times)
        else:
        # otherwise all times are nan means that for all other ign pts, 
        # either the ros does not point towards that point (ros = nan)
        # or ign pt has been ignited (distance = nan)
            self.min_time = np.inf
        
    def get_next_burn(self):
        # find where to burn next
        min_indices = np.argwhere(self.times == self.min_time).flatten()
        next_dirs = self.DIRS[min_indices]

        # take them out of consideration for future burning
        self.distances[min_indices] = np.nan
        self._get_times()

        # return new cell details
        next_cells, next_ign_locs = [], []
        for next_dir in next_dirs:
            next_cells.append(self.NEIGHBORS[next_dir])
            next_ign_locs.append(self.NEXT_IGN[next_dir]) 
        return next_cells, next_ign_locs

    def burn(self, time_step): 
        if self.is_extinguished: 
            print(f'Trying to burn cell {self}, but it is extinguished.')
        else:
            # progress forward in time within cell
            self.burning_time += time_step

        
        # if the cell has been burning too long, extinguish it
        if self.burning_time > self.BURN_OUT_TIME:
            self.is_extinguished = True
        else:
            # otherwise, determine how far to move fire forward 
            dist_steps = self.ros*time_step
            
            # if the ros is nan (meaning the fire isn't moving in that direction), then treat those distances as 0
            dist_steps = np.nan_to_num(dist_steps,nan=0) 
            
            # change distance to ign pts
            self.distances -= dist_steps
            
            # prevent negative distances
            self.distances[self.distances<=0] = 0

            # update times required to get to next ign pts
            self._get_times()
            
            
        
        



            
    
        
        