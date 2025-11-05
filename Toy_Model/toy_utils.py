import numpy as np
import ToyCell as Cell
import ToyCellAuto as ca
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

# maybe create function that opens the sim info and adds timing??

def save_init_info(CA):
    file_path = CA.OUT_DIR + '/simulation_info.txt'
    info_lines = [
        f'Fire ID: {CA.FIRE_ID}',
        f'Number of plots: {CA.NUM_PLOTS}',
        f'Has barrier data? {CA.HAS_BARRIERS}',
        f'Max hours allowed to run: {CA.MAX_HOURS}',
        f'Cell resolution: {CA.RESOLUTION}',
        f'Raster dimensions: {CA.NUM_ROWS} rows by {CA.NUM_COLS} cols',
        f'Burn out time: {CA.BURN_OUT_TIME}',
        f'Fire start date and time: {CA.START_DATE_TIME}',
        f'Data for {CA.TOTAL_HOURS} hours = {CA.TOTAL_DAYS} days',
        f'Ignition ROW, COL = {CA.IGN_LOC}',
    ]
    with open(file_path, "w") as f:
       f.writelines([line + "\n" for line in info_lines])
        
        
def add_info(CA, info):
    file_path = CA.OUT_DIR + '/simulation_info.txt'
    with open(file_path, "a") as f:
       f.write(info + '\n')


def initial_data_plots(CA):
    # plot a bunch of the initial data
    file_path = CA.OUT_DIR + f'/initial_plots.png'
    
    burnhour_masked = CA.burnhour.copy()
    burnhour_masked[CA.burnhour >= CA.MAX_HOURS] = np.nan
    
    min_col, max_col, min_row, max_row = get_plot_bounds([burnhour_masked])

    i = 0  # frame index
    
    num_rows = 4
    num_cols = 4
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_rows*6, num_cols*6))
    
    # === Row 1 ===
    # Barriers
    im = axes[0, 0].imshow(CA.barriers)
    axes[0, 0].set_title('Barriers')
    cbar = fig.colorbar(im, ax=axes[0, 0])
    cbar.ax.set_visible(False) 

    # Wind direction
    im = axes[0, 1].imshow(CA.w_dir_rads[i])
    axes[0, 1].set_title(f'Wind Direction i={i} (rads)')
    fig.colorbar(im, ax=axes[0, 1])

    # Slope direction
    im = axes[0, 2].imshow(CA.s_dir_rads)
    axes[0, 2].set_title('Slope Direction (rads)')
    fig.colorbar(im, ax=axes[0, 2])
    
    # Blank slot
    axes[0, 3].axis('off')

    
    # === Row 2 ===
    # Burnhour all
    im = axes[1, 0].imshow(CA.burnhour)
    axes[1, 0].set_title('Hour of Arrival')
    fig.colorbar(im, ax=axes[1, 0])
    
    # Add extent rectangle
    rect = mpatches.Rectangle((min_col, min_row), width=max_col - min_col, height=max_row - min_row,
                              edgecolor='gray', facecolor='gray', linewidth=2, alpha=0.65)
    axes[1, 0].add_patch(rect)

    # Wind speed
    im = axes[1, 1].imshow(CA.wspd[i])
    axes[1, 1].set_title(f'Windspeed i={i}')
    fig.colorbar(im, ax=axes[1, 1])

    # Slope (rad)
    im = axes[1, 2].imshow(CA.slope_rads)
    axes[1, 2].set_title('Slope (rads)')
    fig.colorbar(im, ax=axes[1, 2])

    # Hourly growth
    im = axes[1, 3].imshow(CA.hourly_growth[i])
    axes[1, 3].set_title(f'Hourly Growth i={i}')
    fig.colorbar(im, ax=axes[1, 3])

    # === Row 3 ===
    # Burnhour masked
    im = axes[2, 0].imshow(burnhour_masked)
    axes[2, 0].set_title('Hour of Arrival')
    axes[2, 0].set_xlim(min_col, max_col)
    axes[2, 0].set_ylim(max_row, min_row)
    fig.colorbar(im, ax=axes[2, 0])

    # Wind magnitude
    im = axes[2, 1].imshow(CA.w_mag[i])
    axes[2, 1].set_title(f'Wind Magnitudei={i}')
    fig.colorbar(im, ax=axes[2, 1])
    
    # Slope magnitude
    im = axes[2, 2,].imshow(CA.s_mag)
    axes[2, 2].set_title('Slope Magnitude')
    fig.colorbar(im, ax=axes[2, 2])
    
    # Growth magnitde
    im = axes[2, 3].imshow(CA.growth_mag[i])
    axes[2, 3].set_title(f'Growth Magnitude i={i}')
    fig.colorbar(im, ax=axes[2, 3])

    # === Row 4 ===
    # Blank slot
    axes[3, 0].axis('off')

    # Histogram: Wind Speed
    ax = axes[3, 1]
    data = CA.wspd.flatten()
    ax.hist(data[~np.isnan(data)], bins=75, color='steelblue')
    ax.set_title('Wind Speed Histogram')
    ax.axvline(np.nanmean(data), color='firebrick', linestyle='--', label=f'Mean = {np.nanmean(data):.2f}', linewidth=2)
    ax.axvline(np.nanmedian(data), color='orange', linestyle='--', label=f'Median = {np.nanmedian(data):.2f}', linewidth=2)
    ax.legend()
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_visible(False) 

    # Histogram: Slope
    ax = axes[3, 2]
    data = CA.slope_rads.flatten()
    ax.hist(data[~np.isnan(data)], bins=50, color='steelblue')
    ax.set_title('Slope Histogram')
    ax.axvline(np.nanmean(data), color='firebrick', linestyle='--', label=f'Mean = {np.nanmean(data):.2f}', linewidth=2)
    ax.axvline(np.nanmedian(data), color='orange', linestyle='--', label=f'Median = {np.nanmedian(data):.2f}', linewidth=2)
    ax.legend()
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_visible(False) 

    # Histogram: Hourly Growth
    ax = axes[3, 3]
    data = CA.hourly_growth.flatten()
    ax.hist(data[~np.isnan(data)], bins=50, color='steelblue')
    ax.set_title('Hourly Growth Histogram')
    ax.axvline(np.nanmean(data), color='firebrick', linestyle='--', label=f'Mean = {np.nanmean(data):.2f}', linewidth=2)
    ax.axvline(np.nanmedian(data), color='orange', linestyle='--', label=f'Median = {np.nanmedian(data):.2f}', linewidth=2)
    ax.legend()
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_visible(False) 
    
    # === Figure title and layout ===
    fig.suptitle(f'{CA.FIRE_ID} \n burning for {CA.MAX_HOURS:.0f} hours', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.savefig(file_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    

def get_plot_bounds(data_arrays): 
    min_col, min_row = data_arrays[0].shape
    max_col, max_row = 0,0
    
    for data in data_arrays:
        non_nan_indices = np.argwhere(~np.isnan(data))
        row_indices = non_nan_indices[:, 0]
        col_indices = non_nan_indices[:, 1]
        temp_min_row = np.min(row_indices)
        temp_max_row = np.max(row_indices)
        temp_min_col = np.min(col_indices)
        temp_max_col = np.max(col_indices)
        min_col = min(min_col, temp_min_col)
        max_col = max(max_col, temp_max_col)
        min_row = min(min_row, temp_min_row)
        max_row = max(max_row, temp_max_row)
        
    min_col -= 5
    max_col += 5
    min_row -= 5
    max_row += 5
    return min_col, max_col, min_row, max_row


def create_burn_ext(CA):
    burn_ext = np.full_like(CA.arrival_map, np.nan)
    for cell in CA.burning_cells:
        burn_ext[cell.ROW][cell.COL] = 2
    for cell in CA.extinguished_cells:
        if cell.burnt_neighbors:
            burn_ext[cell.ROW][cell.COL] = 1
        else:
            burn_ext[cell.ROW][cell.COL] = 0
    return burn_ext

def arrival_plots(CA):
    file_path = CA.OUT_DIR + f'/plot_{CA.curr_hour}_hrs.png'

    data_arrays = [CA.arrival_map]
    min_col, max_col, min_row, max_row = get_plot_bounds(data_arrays)

    burn_ext = create_burn_ext(CA)

    values = [0, 1, 2]
    labels = ['Burned Out', 'Burnt Neighbors', 'Burning']
    colors = ['midnightblue', 'dimgrey', 'firebrick']
    
    cmap = plt.matplotlib.colors.ListedColormap(colors)
    bounds = [v - 0.5 for v in range(4)]  # for 3 values: [-0.5, 0.5, 1.5, 2.5]
    norm = plt.matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 4)) 
    
    
    im = axes[0].imshow(CA.arrival_map, alpha=0.7, cmap='viridis', origin='upper')
    axes[0].set_title('Hour of Arrival')
    g_pixels = np.sum(~np.isnan(CA.arrival_map))
    axes[0].set_xlabel(f'guess: {g_pixels} pixels')
    fig.colorbar(im, ax=axes[0], orientation='vertical')
    
    axes[1].imshow(burn_ext, cmap=cmap, norm=norm, alpha=0.7)
    axes[1].set_xlabel(f'Num burning = {len(CA.burning_cells)}')
    patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(3)]
    axes[1].legend(handles=patches, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    
    for i in range(len(axes)):
        axes[i].set_xlim(min_col,max_col)
        axes[i].set_ylim(max_row,min_row)
    
    fig.suptitle(f'{CA.FIRE_ID} burning for {CA.curr_hour:.0f} hours = {CA.curr_hour/24:.1f} days')
    
    plt.tight_layout(rect=[0, 0, 1, 0.98])

    plt.savefig(file_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

def plot_time_series(CA):
    file_path = CA.OUT_DIR + f'/area_burned_time_series.png'
    days = np.arange(len(CA.fire_growth))
    avg, first, med, third, actual_growth = (np.zeros_like(days) for _ in range(5))
    model_growth = np.full_like(days, np.nan, dtype=float)
    for day in days:
        first[day], med[day], third[day] = np.nanpercentile(CA.fire_growth[day],[25,50,75])
        avg[day] = np.nanmean(CA.fire_growth[day]) 
        hr = day*24
        actual_mask = np.logical_and(CA.burnhour<hr, CA.burnhour>=(hr-24))
        value = np.nansum(actual_mask)
        if value > 0:
            actual_growth[day] = value
        model_mask = np.logical_and(CA.arrival_map<hr, CA.arrival_map>=(hr-24))
        value = np.nansum(model_mask) 
        if value > 0:
            model_growth[day] = value
    
    hours = np.arange(1, CA.MAX_HOURS+1)
    true_area_burned, pred_area_burned, IOU, RMSE = (np.zeros(CA.MAX_HOURS) for _ in range(4))
    for hr in hours[:-1]:
        true_mask = CA.burnhour<hr
        pred_mask = CA.arrival_map<hr
        true = CA.burnhour[true_mask]
        pred = CA.arrival_map[pred_mask] 
        true_masked_pred = CA.burnhour[pred_mask]
        true_mask = true_mask.astype(int)
        pred_mask = pred_mask.astype(int)
        true_area_burned[hr] = true_mask.sum()*CA.RESOLUTION**2/1000000
        pred_area_burned[hr] = pred_mask.sum()*CA.RESOLUTION**2/1000000
        TP = ((pred_mask + true_mask) > 1).sum()
        FP = ((pred_mask - true_mask) > 0).sum()
        FN = ((pred_mask - true_mask) < 0).sum()
        if (TP + FP + FN) == 0:
            IOU[hr] = 0
            RMSE[hr] = 0
        else:
            IOU[hr] = TP /(TP + FP + FN)
            RMSE[hr] = np.sqrt((np.power((true_masked_pred - pred),2)).sum()/pred_mask.sum())
            
    
    fig, axes = plt.subplots(4, 2, figsize=(12, 16), sharex='col')
    
    # --- Plot 1: Truth vs. Guess ---
    ax = axes[0,0]
    ax.plot(hours, true_area_burned, label='Truth')
    ax.plot(hours, pred_area_burned, label='Guess')
    ax.set_title('Total Area Burned (sq km)')
    ax.legend()
    ax.grid(True,alpha = 0.35)
    
    # --- Plot 2: Difference ---
    ax = axes[1,0]
    ax.plot(hours, pred_area_burned - true_area_burned, label='Guess - Truth')
    #ax.set_xlabel('Hour')
    ax.set_title('Area Difference (sq km)')
    ax.legend()
    ax.grid(True,alpha = 0.35)
    
    
    # --- Plot 3: IOU ---
    ax = axes[2,0]
    ax.plot(hours, IOU, label='Intersection over Union')
    ax.set_title('Intersection over Union')
    #ax.legend()
    ax.grid(True,alpha = 0.35)
    
    # --- Plot 4: RMSE ---
    ax = axes[3,0]
    ax.plot(hours, RMSE, label='RMSE')
    ax.set_xlabel('Hour')
    ax.set_title('RMSE')
    #ax.legend()
    ax.grid(True,alpha = 0.35)
    
    # --- Plot 5: Actual Growth ---
    ax = axes[0,1]
    ax.plot(days, actual_growth, label='model')
    ax.set_title('Actual Fire Growth')
    #ax.set_ylim(0,750000)
    #ax.set_yscale('log')
    ax.grid(True,alpha = 0.35)
    
    # --- Plot 6: Predicted Growth ---
    ax = axes[1,1]
    ax.plot(days, avg, label='avg')
    ax.plot(days, first, label='Q1')
    ax.plot(days, med, label='med')
    ax.plot(days, third, label='Q3')
    ax.set_title('Predicted Fire Growth')
    #ax.set_yscale('log')
    ax.legend()
    ax.grid(True,alpha = 0.35)
    
    # --- Plot 7: Model Growth ---
    ax = axes[2,1]
    ax.plot(days, model_growth, label='model')
    ax.set_title('Model Fire Growth')
    ax.grid(True,alpha = 0.35)
    ax.set_xlabel('Day')
    axes[2,1].tick_params(labelbottom=True)
    
    # --- Blank  ---
    axes[3,1].axis('off')
    
    fig.suptitle(f'{CA.FIRE_ID}')
    plt.tight_layout(rect=[0, 0, 1, 0.98])

    plt.savefig(file_path, dpi=300, bbox_inches="tight")
    plt.close()

def plot_pred_vs_model_growth(CA):
    file_path = CA.OUT_DIR + f'/pred_vs_model_growth.png'
    hours = np.arange(CA.MAX_HOURS)
    days = np.arange(0, CA.TOTAL_DAYS)
    
    pred_growth = np.zeros(len(hours))
    model_growth = np.zeros(len(hours))
    for hr in hours:
        pred_growth[hr] = np.nanmean(CA.hourly_growth[hr])
        model_mask = np.logical_and(CA.arrival_map>=hr, CA.arrival_map<(hr+1))
        model_growth[hr] = np.nansum(model_mask)
    
    daily_pred_growth = np.zeros(len(days))
    daily_model_growth = np.zeros(len(days))
    for day in days:
        daily_model_growth[day] = np.nansum(model_growth[day*24:(day+1)*24])
        daily_pred_growth[day] = np.nansum(pred_growth[day*24:(day+1)*24])
        
    cum_pred = np.nancumsum(daily_pred_growth)
    cum_model = np.nancumsum(daily_model_growth)
    
    
    fig, axes = plt.subplots(2, 1, figsize=(6,8), sharex='col')
    
    # --- Plot 1: Daily Growth ---
    ax = axes[0]
    ax.plot(days, daily_pred_growth, label='predicted')
    ax.plot(days, daily_model_growth, label='model')
    ax.set_ylabel('Daily Growth (hectares)')
    ax.legend()
    ax.grid(True,alpha = 0.35)
    
    # --- Plot 2: Cumulative Growth ---
    ax = axes[1]
    ax.plot(days+1,cum_pred,label='predicted')
    ax.plot(days+1, cum_model,label='model')
    ax.set_ylabel('Total Area of Fire (hectares)')
    ax.legend()
    ax.grid(True,alpha = 0.35)
    ax.set_xlabel('Day')
    
    fig.suptitle(f'AOS: {CA.fire_growth[0,0,0]*CA.RESOLUTION**2/10000:0.2f} hectares/day')
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    
    plt.savefig(file_path, dpi=300, bbox_inches="tight")
    plt.close()