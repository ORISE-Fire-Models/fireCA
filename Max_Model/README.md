# FireCA: Cellular Automaton Fire Spread Model

This folder contains the code as it was when Max left the project. There have been attempts at comments, but some of the calcuations lack clear reasoning.

---

## File Structure

- **[Main.py](Main.py)**  
  Entry point for running simulations. Responsible for:
  - Creating a `fireCA` object.  
  - Running data collection from the other folders (`fireCA.get_data()`).  
  - Preprocessing inputs (`fireCA.preprocess()`), including slope/aspect/wind calculations.  
  - Initializing the fire ignition point and RoS vectors (`fireCA.reset()`).  
  - Running the cellular automaton (`fireCA.run()`).  
  - Creates GeoTiff of arrival times (`fireCA.saveArrivalMap()`).

- **fireCA.py**  
  Core cellular automaton logic. Handles data ingestion, preprocessing, model initialization, and simulation.  

- **Cell.py**  
  Defines cell-level behavior, including ignition, burning, extinguishment, and rate-of-spread vector calculations.  
  - Uses **gradient descent** to calibrate RoS vectors.  

- **SimpleCell.py**  
  Inherits from `Cell`, but uses **bisection** instead of gradient descent to calculate RoS.  

- **Stats.py**  
  Contains functions for calculating statistics to assess model performance.  

- **max_model_report.pdf**  
  Max’s explanation of the process and reasoning behind model design choices.  

- **influence_of_slope_on_fire_spread_rate.pdf**  
  Paper used as justification for the slope magnitude calculation. 

---

## Key Model Components
### Preprocessing (`fireCA.preprocess()`)
Performs calculations on data to create a spread vector for each location and time throughout the fire.
- **ROS:** Converts from hectares/day to square meters/hour and multiplies by a growth weight which was used to scale the growth
- **Wind Direction:** Projects onto a "math-style" axis system (0° = +x-axis, 90° = +y-axis) and converts to radians
- **Wind Magnitude:** Converts windspeed ($w$) from meters/second to km/hour and then determimes the magnitude by calculating $\frac{e^{0.05w}}{1+e^{0.6(10-w)}}$
- **Slope Direction:** Calculates upslope direction from aspect and projects onto "math-style" axis system
- **Slope Magnitude:** Converts slope ($s$) to degrees and determines the magnitude by calculating $0.2e^{0.13s}+0.6$
- The wind and slope vectors are combined to create a spread vector where `r_max_cube` gives the magnitude and `phi_max_cube` gives the direction. 
  - This [graph](https://www.desmos.com/calculator/d4ca117200) gives a good representation of these calculations.  
- **NOTES** 
  - The spread magnitude (`r_max_cube`) is also multiplied by the square root of the predicted fire growth (ROS) for an unclear reason
  - Converting to mathematical axes seems like an unnecessary complication and may be done incorrectly
  - The paper mentioned in the slope comment is included here: [Influence of Slope on Fire Spread Rate (PDF)](influence_of_slope_on_fire_spread_rate.pdf).  

---
 
### Initialization (`fireCA.reset()`)
- Creates cells at the ignition locations based off the VIIRS burnday raster
- Uses `Cell.calcRoS()` to determine the spread within the cells. This is done by:
  - Projecting the wind/slope spread vector onto vectors to neighboring cells.
    - If this projection is negative (because the angle between the two vectors is greater than 90°), then it is replaced with 0.
  - Calculating the area of the polygon created by these vectors.
  - Running gradient descent (or bisection) to minimize the squared difference between the area of this polygon and the predicted growth amount.
  - Throughout the minimization process, the polygon is forced to stay convex and have no negative vectors.
  - **NOTE:** There are calculations within the functions that make up this process (like f, df, and convexify) that are only made on certain angles if the ignition point is the center of a cell (which is really only the first cells). I think this is an issue. 

### Running the Automaton (`fireCA.run()`)
- Advances time and burning within the cells. 
- The time step is determined by finding the smallest amount of time required to ignite another cell. 
- Each hour, all cells are updated to reflect updated weather conditions. 
- Each day, the growth predictions are also updated within all the cells.
- Cells are extinguished after all neighbors have burned.

---

## Getting Started

### Requirements
- Python 3.9+
- numpy
- scipy
- matplotlib

### Example Run
```bash
python Main.py --indir=/mnt/DataDrive1/data/zholden/VIIRS/CA_model_test/raster_inputs_daily/ --outdir=/mnt/DataDrive1/data/zholden/VIIRS/CA_model_test/<USER_OUTPUT_FOLDER>/ --fireID=MT4703711487920170717 --res=90 --max_iteration=400 --n_burnt=False --plotting=True --hasBarrier=False --isSimpleCell=False --growthWeight=1
```

### Explanation of Example Run Arguments

The example run executes the **`Main.py`** script and passes various parameters that control the model's input, output, and behavior. Each argument is a **flag** prefixed with `--`.

| Argument | Description | Example Value |
| :--- | :--- | :--- |
| **`--indir`** | **Input Directory:** The file path to the folder containing the input data. | `/mnt/DataDrive1/.../raster_inputs_daily/` |
| **`--outdir`** | **Output Directory:** The file path to the folder where all model outputs will be saved. **Best to create a new folder** | `/mnt/DataDrive1/.../<USER_OUTPUT_FOLDER>/` |
| **`--fireID`** | **Fire Identification:** The ID of the fire to be simulated. | `MT4703711487920170717` |
| **`--res`** | **Spatial Resolution:** Measured in **meters**. | `90` |
| **`--max_iteration`** | **Maximum Iterations:** The maximum number of simulation time steps. | `400` |
| **`--n_burnt`** | **Number Burnt:** A boolean flag that determines if the model will stop after 1000 cells have burned. | `False` |
| **`--plotting`** | **Enable Plotting:** A boolean flag that determines if simulation plots are saved. | `True` |
| **`--hasBarrier`** | **Barrier Inclusion:** A boolean flag to indicate if the simulation uses the final perimeter of the VIIRS data to restrit the simulation. | `False` |
| **`--isSimpleCell`** | **Simple Cell Model:** A boolean flag that determines whether gradient descent (`False`) or bisection (`True`) is used. | `False` |
| **`--growthWeight`**** | **Growth Weight Factor:** Factor to adjust the growth rate. Max seemed to have success with approx `0.1`. | `1` |

