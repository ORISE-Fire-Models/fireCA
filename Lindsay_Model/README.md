# CellAuto: Cellular Automaton Fire Spread Model

This folder contains the code as it was when Lindsay left the project. The code should be well commented.

---

## File Structure

- **runCA.py**  
  Entry point for running simulations. Responsible for:
  - Assigning values for the `CellAuto` class arguments.
  - Creating a `CellAuto` instance.  
  - Running the cellular automaton (`CellAuto.run()`) simulation.  
  - Currently set to run through many fires
 
- **CellAuto.py**  
  Core cellular automaton logic. Handles data ingestion, preprocessing, model initialization, simulation, and plotting.  

- **Cell.py**  
  Defines cell-level behavior, including ignition, burning, extinguishment, and rate-of-spread vector calculations.  

- **utils.py**  
  Contains functions for logging simulation information along with calculating and plotting statistics to assess model performance.  

---

## Key Model Components
### Preprocessing (`CellAuto._process_data()`)
- Initializes the fire simulation, then reads in data and performs calculations to create rasters of wind, slope, and growth prediction for each location and time throughout the fire. 

 
### ROS Calcuation (`Cell.update_cell_growth()`)
- From the iginition point within a cell, there is a ROS vector to each of the eight surrounding cells.
- This ROS in each direction is calculated by $ ROS = g*w_f*s_f $
- where:
    - $g$ is the growth prediction from the growth prediction raster at the current time and location. To deal with runaway exponential growth, the growth prediction is divided by the number of currently burning cells before it is passed into this function.
    - $w_f$ is the wind factor and $s_f$ is the slope factor. 
    - Each of these factors is calculated by first projecting the wind vector (from the windspeed and wind angle rasters at the current time and location)  or the slope vector (from the slope and slope angle rasters at the current location) onto the current angle. 
    - Note that these projected values ($w_p$ or $s_p$) are postive if the angle between the wind (or slope) and the direction to the next cell is less than $90\degree$. If the angle is more than $90\degree$, then the projected value is negative and a backwind (or downslope) fire is being considered. 
    - These projetions are then passed through the functions $w_f$ and $s_f$ given below.  
    - The graph of $w_f$ can be visualized [here](https://www.desmos.com/calculator/bmcsrmhfn3) while the graph of $s_f$ can be visualized [here](https://www.desmos.com/calculator/iuza8zwfnz).
$$
w_f(w_p) = 
\begin{cases} 
\frac{2}{1 + e^{-0.1w_p}} &  w_p < 0 \\[2mm]
\frac{15}{1 + 14e^{-0.25w_p}} &  w_p > 0
\end{cases}
\quad
s_f(s_p) = 
\begin{cases} 
\frac{2}{1 + e^{-s_p}} &  s_p < 0 \\[2mm]
\frac{5}{1 + 4e^{-10s_p}} &  s_p > 0
\end{cases}
$$


- For both the wind and slope factors, the goal is to slow the growth prediction in the opposite direction (backwind or downslope) and increase the growth in the direction that aligns with the factor. Therefore, the wind and slope factor functions output values less than one for directions opposite the factors and values bigger than one in directions that align with the factors.
- The model was repeatedly tested while varying values within the factor functions to identify the parameters that yielded the best performance.  


### Running the Automaton (`CellAuto.run()`)
- Advances time and burning within the cells. 
- The time step is determined by finding the smallest amount of time required to ignite another cell or to get to the next hour (because then the ROS must be updated). 
- Each hour, all cells are updated to reflect updated weather conditions and growth predictions. 
- Cells are extinguished after all neighbors have burned or after they have burned longer than the burn-out time.

---
## Model Testing and Development

Throughout development, several model variations were tested to explore sensitivity and improve fire spread realism. Examples include:

- **Local normalization:** Instead of dividing by the total number of burning cells, dividing by the number of burning cells within a given radius was attempted. No matter what this radius was, the growth eventually grew exponentially bigger than the fire was predicted to grow.  
- **Altering ignition start times:** Shifting ignition timing by a set number of hours to test sensitivity to environmental conditions and growth predictions. This seemed to not have much effect.
- **ROS factor tuning:** As mentioned above, the coefficients in the wind and slope factor functions systematically adjusted to tune the function parameters.

These experiments guided parameter selection and informed refinements to improve model stability and performance across different fire events.

---

## Getting Started

### Requirements
- Python 3.9+
- numpy
- scipy
- rasterio
- matplotlib
- pandas
- scikit-image
- scikit-learn

### Example Run
First, edit runCA.py to assign all necessary values. Then, run the line below from the terminal. 
```bash
python runCA.py 
```

