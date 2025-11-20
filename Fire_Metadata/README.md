# Fire Metadata

This folder contains the compiled fire-level metadata used to align ignition timing, environmental data availability, and suppression activity records across datasets.

---

## File Overview

**fire_info.csv**  
A summary table containing metadata for each fire with available data.
 
**Columns**  
- **fireID**: Unique fire identifier.
- **fire_name**: The fire name. 
- **fire_start**:Estimated ignition date and time based on VIIRS detections.
- **wind_start**:First available timestamp in the wind dataset (derived from NetCDF metadata).
- **hours_of_data** and **days_of_data**: Number of hours and days for which wind and growth prediction data are available.
- **peak_aerial** and **peak_personnel**: Date and time of peak aerial and personnel suppression activity derived from ICS-209-PLUS.
- **first_suppression**: The earliest of the two peak suppression times (aerial or personnel).
- **hours_til_first, hours_til_aerial, hours_til_personnel, days_til_first, days_til_aerial,** and **days_til_personnel**: Time since ignition to each of the suppression metrics.
- **fs209_indices**: The ICS-209-PLUS record IDs from which the suppression information was derived.

---

## Supporting Scripts

The following scripts and notebooks were used to compile and cross-reference the information in `fire_info.csv`. These files are not fully organized or well-documented, but they outline the general data collection process.

- **getFireInfo.py**
    Script used to merge ignition, environmental, and suppression metadata across sources.
- **getSuppressionInfo.ipynb**
    Notebook used to extract and summarize suppression timing (aerial and personnel peaks) from ICS-209-PLUS data.
- **windMeta.ipynb**
    Notebook for reading and summarizing metadata from wind NetCDF files.
