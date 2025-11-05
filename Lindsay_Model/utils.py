import numpy as np
import Cell as Cell
import CellAuto as ca
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from skimage.metrics import structural_similarity as ssim
from scipy.stats import pearsonr
from sklearn.metrics import r2_score

def save_init_info(CA):
    # saves simulation information to a text file
    file_path = f'{CA.OUT_DIR}/{CA.NAME}_simulation_info.txt'
    info_lines = [
        f'Fire Name: {CA.NAME}',
        f'Fire ID: {CA.FIRE_ID}',
        f'Number of plots: {CA.NUM_PLOTS}',
        f'Has barrier data? {CA.HAS_BARRIERS}',
        f'Max time allowed to run: {CA.MAX_HOURS} hours = {CA.MAX_HOURS/24:0.2f} days',
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
    # adds simulation information to text file
    file_path = f'{CA.OUT_DIR}/{CA.NAME}_simulation_info.txt'
    with open(file_path, "a") as f:
       f.write(info + '\n')


def write_summary_info(CA, AREA, IOU, TVER, SSIM, RMSE, R, R2):
    # writes simulation summary to text file
    file_path = ('/').join(CA.OUT_DIR.split('/')[:-1])+'/'+CA.NAME+'_summary.txt'
    lines = [
        f'Run: {CA.OUT_DIR.split('/')[-1]}',
        f'==== Model vs VIIRS Area Size Stats =====',
        f'RMSE = {RMSE:0.4f}',
        f'Pearson Correlation = {R:0.4f}',
        f'R2 = {R2:0.4f}',
        f'==== Spatial Stats =====',
        f'Average IOU: {np.nanmean(IOU):0.4f}',
        f'Area Weighted Average IOU: {np.nansum(AREA*IOU)/np.nansum(AREA):0.4f}',
        f'Average Tversky: {np.nanmean(TVER):0.4f}',
        f'Area Weighted Average Tversky: {np.nansum(AREA*TVER)/np.nansum(AREA):0.4f}',
        f'Average SSIM: {np.nanmean(SSIM):0.4f}',
    ]
    with open(file_path, "a") as f:
        f.writelines([line + "\n" for line in lines])
        f.write('\n' + '\n')

    
def initial_data_plots(CA):
    # plot a bunch of the initial data
    file_path = f'{CA.OUT_DIR}/{CA.NAME}_initial_plots.png'
    
    burnhour_masked = CA.burnhour.copy()
    burnhour_masked[CA.burnhour >= CA.MAX_HOURS] = np.nan
    
    min_col, max_col, min_row, max_row = get_plot_bounds([burnhour_masked])
    
    hours = np.arange(CA.MAX_HOURS)
    days = np.arange(0, CA.MAX_HOURS//24)
    
    actual_growth, pred_growth, pred_adj_growth, cum_pred, cum_pred_adj = (np.zeros(len(hours)) for _ in range(5))
    
    for hr in hours:
        pred_growth[hr] = np.nanmean(CA.hourly_growth[hr])*CA.RESOLUTION**2/10000
        pred_adj_growth[hr] = np.nanmean(CA.growth_mag[hr])*CA.RESOLUTION**2/10000
        actual_mask = np.logical_and(CA.burnhour>=hr, CA.burnhour<(hr+1))
        actual_growth[hr] = np.nansum(actual_mask)*CA.RESOLUTION**2/10000
        
        if hr%24 == 0:
            cum_pred[hr] =  np.nansum(actual_growth[:hr+1])
            cum_pred_adj[hr] =  np.nansum(actual_growth[:hr+1])
        else:
            cum_pred[hr] = cum_pred[hr-1] + pred_growth[hr]
            cum_pred_adj[hr] = cum_pred_adj[hr-1] + pred_adj_growth[hr]
    
    cum_actual = np.nancumsum(actual_growth)
    
    daily_pred_growth, daily_pred_adj_growth, daily_actual_growth = (np.zeros(len(days)) for _ in range(3))
    
    for day in days:
        daily_pred_growth[day] = np.nansum(pred_growth[day*24:(day+1)*24])
        daily_pred_adj_growth[day] = np.nansum(pred_adj_growth[day*24:(day+1)*24])
        daily_actual_growth[day] = np.nansum(actual_growth[day*24:(day+1)*24])     
    
    i = 0  # frame index
    
    num_rows = 4
    num_cols = 4
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_rows*6, num_cols*6))
    
    # === Row 1 ===
    # Barriers
    ax = axes[0, 0]
    im = ax.imshow(CA.barriers)
    ax.set_title('Barriers')
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_visible(False) 
    
    # Wind direction
    ax = axes[0, 1]
    im = ax.imshow(CA.w_dir_rads[i])
    ax.set_title(f'Wind Direction hr={i} (rads)')
    fig.colorbar(im, ax=ax)
    
    # Slope direction
    ax = axes[0, 2]
    im = ax.imshow(CA.s_dir_rads)
    ax.set_title('Slope Direction (rads)')
    fig.colorbar(im, ax=ax)
    
    # Daily growth
    ax = axes[0, 3]
    im = ax.imshow(CA.fire_growth[i]*CA.RESOLUTION**2/10000)
    ax.set_title(f'Daily Growth day={i} (hectares/day)')
    fig.colorbar(im, ax=ax)
    
    
    # === Row 2 ===
    # Burnhour all
    ax = axes[1, 0]
    im = ax.imshow(CA.burnhour)
    ax.set_title('Hour of Arrival')
    fig.colorbar(im, ax=ax)
    
    # Add extent rectangle
    rect = mpatches.Rectangle((min_col, min_row), width=max_col - min_col, height=max_row - min_row,
                              edgecolor='gray', facecolor='gray', linewidth=2, alpha=0.65)
    ax.add_patch(rect)
    
    # Histogram: Wind Direction
    ax = axes[1, 1]
    data = CA.w_dir_rads[:CA.MAX_HOURS+1].flatten()
    ax.hist(data[~np.isnan(data)], bins=75, color='steelblue')
    ax.set_title('Wind Direction Histogram (rads)')
    ax.axvline(np.nanmean(data), color='firebrick', linestyle='--', label=f'Mean = {np.nanmean(data):.2f}', linewidth=2)
    ax.axvline(np.nanmedian(data), color='orange', linestyle='--', label=f'Median = {np.nanmedian(data):.2f}', linewidth=2)
    ax.legend()
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_visible(False)
    
    # Histogram: Slope Direction
    ax = axes[1, 2]
    data = CA.s_dir_rads.flatten()
    ax.hist(data[~np.isnan(data)], bins=50, color='steelblue')
    ax.set_title('Slope Direction Histogram (rads)')
    ax.axvline(np.nanmean(data), color='firebrick', linestyle='--', label=f'Mean = {np.nanmean(data):.2f}', linewidth=2)
    ax.axvline(np.nanmedian(data), color='orange', linestyle='--', label=f'Median = {np.nanmedian(data):.2f}', linewidth=2)
    ax.legend()
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_visible(False)
    
    # Histogram: Daily Growth
    ax = axes[1, 3]
    data = (CA.fire_growth[:CA.MAX_HOURS+1]*CA.RESOLUTION**2/10000).flatten()
    ax.hist(data[~np.isnan(data)], bins=50, color='steelblue')
    ax.set_title('Daily Growth Histogram (hectares/day)')
    ax.axvline(np.nanmean(data), color='firebrick', linestyle='--', label=f'Mean = {np.nanmean(data):.2f}', linewidth=2)
    ax.axvline(np.nanmedian(data), color='orange', linestyle='--', label=f'Median = {np.nanmedian(data):.2f}', linewidth=2)
    ax.legend()
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_visible(False) 
    
    
    # === Row 3 ===
    # Burnhour masked
    ax = axes[2, 0]
    im = ax.imshow(burnhour_masked)
    ax.set_title('Hour of Arrival')
    ax.set_xlim(min_col, max_col)
    ax.set_ylim(max_row, min_row)
    fig.colorbar(im, ax=ax)
    
    # Wind speed
    ax = axes[2, 1]
    im = ax.imshow(CA.wspd[i])
    ax.set_title(f'Windspeed hr={i} (m/sec?)')
    fig.colorbar(im, ax=ax)
    
    # Slope (rad)
    ax = axes[2, 2]
    im = ax.imshow(CA.slope_rads)
    ax.set_title('Slope (rads)')
    fig.colorbar(im, ax=ax)
    
    # Time Series: Predicted/Actual Daily Growth 
    ax = axes[2, 3]
    ax.plot(days, daily_pred_growth, label='avg predicted')
    ax.plot(days, daily_pred_adj_growth, label='avg predicted adjusted', color='tab:blue')
    ax.plot(days, daily_actual_growth, label='VIIRS', color='tab:orange')
    ax.set_title('Daily Growth (hectares)')
    ax.legend()
    ax.grid(True,alpha = 0.35)
    ax.set_xlabel('Day')
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_visible(False)

    
    # === Row 4 ===
    # Blank slot
    axes[3, 0].axis('off')
    
    # Histogram: Wind Speed
    ax = axes[3, 1]
    data = CA.wspd[:CA.MAX_HOURS+1].flatten()
    ax.hist(data[~np.isnan(data)], bins=75, color='steelblue')
    ax.set_title('Windspeed Histogram (m/sec?)')
    ax.axvline(np.nanmean(data), color='firebrick', linestyle='--', label=f'Mean = {np.nanmean(data):.2f}', linewidth=2)
    ax.axvline(np.nanmedian(data), color='orange', linestyle='--', label=f'Median = {np.nanmedian(data):.2f}', linewidth=2)
    ax.legend()
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_visible(False) 
    
    # Histogram: Slope
    ax = axes[3, 2]
    data = CA.slope_rads.flatten()
    ax.hist(data[~np.isnan(data)], bins=50, color='steelblue')
    ax.set_title('Slope Histogram (rads)')
    ax.axvline(np.nanmean(data), color='firebrick', linestyle='--', label=f'Mean = {np.nanmean(data):.2f}', linewidth=2)
    ax.axvline(np.nanmedian(data), color='orange', linestyle='--', label=f'Median = {np.nanmedian(data):.2f}', linewidth=2)
    ax.legend()
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_visible(False) 
    
    # Time Series: Predicted/Actual Total Growth
    ax = axes[3, 3]
    ax.plot(hours/24,cum_pred,label='avg predicted')
    ax.plot(hours/24,cum_pred_adj,label='avg predicted adjusted', color='tab:blue')
    ax.plot(hours/24, cum_actual,label='VIIRS', color='tab:orange')
    ax.set_title('Total Area of Fire (hectares)')
    ax.legend()
    ax.grid(True,alpha = 0.35)
    ax.set_xlabel('Day')
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_visible(False) 
    
    # === Figure title and layout ===
    fig.suptitle(f'{CA.NAME}', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.savefig(file_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    

def get_plot_bounds(data_arrays): 
    min_col, min_row = data_arrays[0].shape
    max_col, max_row = 0,0
    
    for data in data_arrays:
        if not np.isnan(data).all():
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
    # creates a raster that demonstrates a cell's state
    # 2 = burning, 1 = extinguished beause all neighbors are burning
    # 0 = extinguished because burned longer than the burn out time
    burn_ext = np.full_like(CA.burnhour, np.nan)
    for cell in CA.burning_cells:
        burn_ext[cell.ROW][cell.COL] = 2
    for cell in CA.extinguished_cells:
        if cell.burnt_neighbors:
            burn_ext[cell.ROW][cell.COL] = 1
        else:
            burn_ext[cell.ROW][cell.COL] = 0
    return burn_ext

def create_burn_ext2(CA):
    # creates a raster that demonstrates a cell's state
    # 2 = burning, 1 = extinguished 
    burn_ext = np.full_like(CA.burnhour, np.nan)
    for row, col in CA.burning_coords:
        burn_ext[row][col] = 2
    for row, col in CA.extinguished_coords:
        burn_ext[row][col] = 1
    return burn_ext

def arrival_plots(CA):
    # plot the model and actual arrival maps at any time during the simulation
    file_path = f'{CA.OUT_DIR}/{CA.NAME}_plot_{CA.curr_hour}_hrs.png'

    max_burnhour = CA.curr_hour
    burnhour_masked = CA.burnhour.copy()
    burnhour_masked[CA.burnhour >= max_burnhour] = np.nan
    burnhour_masked[CA.burnhour < CA.int_day*24] = np.nan

    burn_ext = create_burn_ext2(CA)
    
    data_arrays = [burnhour_masked, CA.arrival_map[CA.int_day],burn_ext]
    min_col, max_col, min_row, max_row = get_plot_bounds(data_arrays)
    
    values = [0, 1, 2]
    labels = ['Burned Out', 'Burnt Neighbors', 'Burning']
    colors = ['midnightblue', 'dimgrey', 'firebrick']
    
    cmap = plt.matplotlib.colors.ListedColormap(colors)
    bounds = [v - 0.5 for v in range(4)]  # for 3 values: [-0.5, 0.5, 1.5, 2.5]
    norm = plt.matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    
    true_false = np.full_like(CA.burnhour, np.nan)
    true_pos = np.logical_and(~np.isnan(burnhour_masked), ~np.isnan(CA.arrival_map[CA.int_day]))
    false_pos = np.logical_and(np.isnan(burnhour_masked), ~np.isnan(CA.arrival_map[CA.int_day]))
    false_neg = np.logical_and(~np.isnan(burnhour_masked), np.isnan(CA.arrival_map[CA.int_day]))
    true_false[true_pos] = 1
    true_false[false_pos] = 2
    true_false[false_neg] = 0
    
    values1 = [0, 1, 2]
    labels1 = ['False Neg', 'True Pos', 'False Pos']
    colors1 = ['gold', 'forestgreen', 'firebrick']
    
    cmap1 = plt.matplotlib.colors.ListedColormap(colors1)
    norm1 = plt.matplotlib.colors.BoundaryNorm(bounds, cmap1.N)
    
    num_rows = 2
    num_cols = 2
    
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(6*num_cols, 4*num_rows)) 
    
    # === ROW 1 ===
    # --- Model ---
    ax = axes[0, 0]
    im = ax.imshow(CA.arrival_map[CA.int_day], alpha=0.7, cmap='viridis', origin='upper')
    g_pixels = np.sum(~np.isnan(CA.arrival_map[CA.int_day]))
    ax.set_title(f'Model Hour of Arrival {g_pixels} pixels')
    fig.colorbar(im, ax=ax, orientation='vertical')
    
    # --- Burn/Ext ---
    ax = axes[0, 1]
    ax.imshow(CA.hillshade, alpha=0.6, cmap='gray', origin='upper')
    ax.imshow(burn_ext, cmap=cmap, norm=norm, alpha=0.7)
    patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(3)]
    ax.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_visible(False)
    
    
    # === ROW 2 ===
    # --- Actual ---
    ax = axes[1,0]
    im = ax.imshow(burnhour_masked, alpha= 0.7, cmap='viridis', origin='upper')
    t_pixels = np.sum(~np.isnan(burnhour_masked))
    ax.set_title(f'VIIRS Hour of Arrival {t_pixels} pixels')
    #ax.set_xlabel(f'truth: {t_pixels} pixels')
    fig.colorbar(im, ax=ax, orientation='vertical')
    
    # --- True/False ---
    ax = axes[1,1]
    ax.imshow(CA.hillshade, alpha=0.6, cmap='gray', origin='upper')
    ax.imshow(true_false, cmap=cmap1, norm=norm1, alpha=0.5)
    patches1 = [mpatches.Patch(color=colors1[i], label=labels1[i]) for i in range(3)]
    ax.legend(handles=patches1, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_visible(False)
    
    for ax in axes.flatten():
        ax.set_xlim(min_col,max_col)
        ax.set_ylim(max_row,min_row)
    
    fig.suptitle(f'{CA.NAME} Fire burning for {CA.curr_hour:.0f} hours = {CA.curr_hour/24:0.1f} days')
    
    plt.tight_layout(rect=[0, 0, 1, 0.98])

    plt.savefig(file_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def normalize(array, mask):
    array[~mask] = np.nan
    rows = np.any(mask, axis=1)
    cols = np.any(mask, axis=0)
    rmin, rmax = np.where(rows)[0][[0, -1]]
    if rmax <= rmin+1:
        if rmax < array.shape[0] - 3:
            rmax += 3
        else:
            rmin -= 3
    cmin, cmax = np.where(cols)[0][[0, -1]]
    if cmax <= cmin+1:
        if cmax < array.shape[1] - 3:
            cmax += 3
        else:
            cmin -= 3
    array /= np.nanmax(array) #if all nan, just gives back all nans
    array = np.nan_to_num(array[rmin:rmax+1, cmin:cmax+1], nan=0) #if all nan, just changes to all 0
    return array

    
def plot_time_series(CA):
    # plot stats comparing the model, the actual, and predicted throughout the simulation
    hours = np.arange(CA.MAX_HOURS+1)
    
    actual_growth, pred_growth, pred_adj_growth, model_growth, = (np.full(len(hours), np.nan) for _ in range(4))
    cum_model, cum_pred, cum_pred_adj = (np.full(len(hours), np.nan) for _ in range(3))
    
    for hr in hours:
        pred_growth[hr] = np.nanmean(CA.hourly_growth[hr])*CA.RESOLUTION**2/10000
        pred_adj_growth[hr] = np.nanmean(CA.growth_mag[hr])*CA.RESOLUTION**2/10000
        model_mask = np.logical_and(CA.arrival_map[hr//24]>=hr, CA.arrival_map[hr//24]<(hr+1))
        model_growth[hr] = np.nansum(model_mask)*CA.RESOLUTION**2/10000
        actual_mask = np.logical_and(CA.burnhour>=hr, CA.burnhour<(hr+1))
        actual_growth[hr] = np.nansum(actual_mask)*CA.RESOLUTION**2/10000
        if hr%24 == 0:
            cum_model[hr] =  np.nansum(actual_growth[:hr+1])
            cum_pred[hr] =  np.nansum(actual_growth[:hr+1])
            cum_pred_adj[hr] =  np.nansum(actual_growth[:hr+1])
        else:
            cum_model[hr] = cum_model[hr-1] + model_growth[hr]
            cum_pred[hr] = cum_pred[hr-1] + pred_growth[hr]
            cum_pred_adj[hr] = cum_pred_adj[hr-1] + pred_adj_growth[hr]
    
    #cum_pred = np.nancumsum(pred_growth)
    #cum_pred_adj = np.nancumsum(pred_adj_growth)
    cum_actual = np.nancumsum(actual_growth)

    days = np.arange(0, CA.MAX_HOURS//24)
    daily_model_growth, daily_actual_growth = (np.zeros(len(days)) for _ in range(2))
    
    for day in days:
        daily_model_growth[day] = np.nansum(model_growth[day*24:(day+1)*24])
        daily_actual_growth[day] = np.nansum(actual_growth[day*24:(day+1)*24])
    
    RMSE = np.sqrt(np.mean((daily_model_growth - daily_actual_growth)**2))
    if len(days) > 1:
        R, _ = pearsonr(daily_model_growth, daily_actual_growth)
        R2 = r2_score(daily_model_growth, daily_actual_growth)
    else:
        R = 0
        R2 = -1
    
    TP, FP, FN, IOU, TVER, SSIM = (np.full(len(hours), np.nan) for _ in range(6))
    
    for hr in hours:
        if hr%24 == 0:
            continue
        burnhour_masked = CA.burnhour.copy()
        burnhour_masked[CA.burnhour >= hr] = np.nan
        burnhour_masked[CA.burnhour >= hr] = np.nan
        burnhour_masked[CA.burnhour < hr//24*24] = np.nan
        arrival_masked = CA.arrival_map[hr//24].copy()
        arrival_masked[CA.arrival_map[hr//24] >= hr] = np.nan
        TP[hr] = np.logical_and(~np.isnan(burnhour_masked), ~np.isnan(arrival_masked)).sum()
        FP[hr] = np.logical_and(np.isnan(burnhour_masked), ~np.isnan(arrival_masked)).sum()
        FN[hr] = np.logical_and(~np.isnan(burnhour_masked), np.isnan(arrival_masked)).sum()
        if (TP[hr] + FP[hr] + FN[hr]) == 0:
            IOU[hr] = 0
            TVER[hr] = 0
        else:
            IOU[hr] = TP[hr] /(TP[hr] + FP[hr] + FN[hr])
            TVER[hr] = TP[hr] /(TP[hr] + 0.3*FP[hr] + 0.7*FN[hr])
        if hr == 0:
            continue
        ssim_mask = np.logical_or(CA.burnhour < hr, CA.arrival_map[hr//24] < hr)
        burnhour_ssim = normalize(CA.burnhour.copy(), ssim_mask)
        arrival_ssim = normalize(CA.arrival_map[hr//24].copy(), ssim_mask)
        win = min(burnhour_ssim.shape[0], burnhour_ssim.shape[1], 7)
        if win%2 == 0:
            win -= 1
        score, ssim_map = ssim(
            burnhour_ssim,
            arrival_ssim,
            data_range=np.nanmax([burnhour_ssim, arrival_ssim]),
            win_size=win,
            full=True,
            gaussian_weights=True
        )
        SSIM[hr] = score

    write_summary_info(CA, cum_actual, IOU, TVER, SSIM, RMSE, R, R2)
    
    num_rows = 3
    num_cols = 2

    # === Plot Hourly ===
    file_path = f'{CA.OUT_DIR}/{CA.NAME}_area_burned_hourly_time_series.png'
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(6*num_cols,4*num_rows), sharex='col')
    
    # --- Plot 1: Hourly Growth ---
    ax = axes[0,0]
    ax.plot(hours, pred_growth, label='avg predicted')
    ax.plot(hours, pred_adj_growth, label='avg predicted adjusted')
    ax.plot(hours, actual_growth, label='VIIRS')
    ax.plot(hours, model_growth, label='model')
    ax.set_title('Hourly Growth (hectares)')
    ax.legend()
    ax.grid(True,alpha = 0.35)
    all_vals = np.concatenate([pred_growth, pred_growth*10, model_growth])
    max_cutoff = min(1.5*np.nanmax(all_vals),np.max(actual_growth))
    ax.set_ylim(0,max_cutoff)
    
    # --- Plot 2: Cumulative Growth ---
    ax = axes[1,0]
    ax.plot(hours, cum_pred,label='avg predicted')
    ax.plot(hours, cum_pred_adj,label='avg predicted adjusted', color='tab:blue')
    ax.plot(hours, cum_actual,label='VIIRS', color='tab:orange')
    ax.plot(hours, cum_model,label='model', color='tab:red')
    ax.set_title('Total Area of Fire (hectares)')
    ax.legend()
    ax.grid(True,alpha = 0.35)
    
    # --- Plot 3: Area Difference ---
    ax = axes[2,0]
    ax.plot(hours,cum_model-cum_actual,label='model - VIIRS', color='tab:red')
    #ax.plot(hours,cum_model-cum_pred_adj,label='model - avg predicted adjusted', color='tab:gray')
    #ax.plot(hours, cum_pred_adj-cum_actual,label='avg predicted adjusted - VIIRS', color='tab:orange')
    ax.set_title('Difference (hectares)')
    ax.legend()
    ax.grid(True,alpha = 0.35)
    ax.set_xlabel('Hour')

    # --- Plot 4: Confusion Matrix ---
    ax = axes[0,1]
    tot_pix = CA.NUM_ROWS*CA.NUM_COLS
    ax.plot(hours, TP/tot_pix, label='True Positive', color='tab:green')
    ax.plot(hours, FP/tot_pix, label='False Positive', color='tab:red')
    ax.plot(hours, FN/tot_pix, label='False Negative', color='tab:orange')
    ax.set_title('Confusion Matrix')
    ax.legend()
    ax.grid(True,alpha = 0.35)

    # --- Plot 5: Accuracy ---
    ax = axes[1,1]
    ax.plot(hours, IOU, label='IOU')
    ax.plot(hours, TVER, label='Tversky')
    ax.plot(hours, SSIM, label='SSIM')
    ax.set_title('Accuracy')
    ax.legend()
    ax.grid(True,alpha = 0.35)
    
    # --- Plot 6: Num Burning ---
    ax = axes[2,1]
    ax.plot(hours, CA.num_burning)
    ax.set_title('Number Cells Burning')
    ax.grid(True,alpha = 0.35)
    ax.set_xlabel('Hour')
    
    fig.suptitle(f'{CA.NAME} Fire burning for {CA.MAX_HOURS:0.0f} hours = {CA.MAX_HOURS/24:0.1f} days')
    plt.tight_layout(rect=[0, 0, 1, 0.98])

    plt.savefig(file_path, dpi=300, bbox_inches="tight")
    plt.close()

    '''
    # === Plot Daily ===
    file_path = f'{CA.OUT_DIR}/{CA.NAME}_area_burned_daily_time_series.png'

    days = np.arange(0, CA.MAX_HOURS//24)
    daily_pred_growth, daily_pred_adj_growth, daily_model_growth, daily_actual_growth = (np.zeros(len(days)) for _ in range(4))
    
    for day in days:
        daily_model_growth[day] = np.nansum(model_growth[day*24:(day+1)*24])
        daily_pred_growth[day] = np.nansum(pred_growth[day*24:(day+1)*24]) 
        daily_pred_adj_growth[day] = np.nansum(pred_adj_growth[day*24:(day+1)*24]) 
        daily_actual_growth[day] = np.nansum(actual_growth[day*24:(day+1)*24]) 
    
    cum_pred = np.nancumsum(daily_pred_growth)
    cum_pred_adj = np.nancumsum(daily_pred_adj_growth)
    cum_model = np.nancumsum(daily_model_growth)
    cum_actual = np.nancumsum(daily_actual_growth)
    
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(6*num_cols,4*num_rows), sharex='col')
    
    # --- Plot 1: Daily Growth ---
    ax = axes[0,0]
    ax.plot(days, daily_pred_growth, label='avg predicted')
    ax.plot(days, daily_pred_adj_growth, label='avg predicted adjusted')
    ax.plot(days, daily_actual_growth, label='VIIRS')
    ax.plot(days, daily_model_growth, label='model')
    ax.set_title('Daily Growth (hectares)')
    ax.legend()
    ax.grid(True,alpha = 0.35)
    
    # --- Plot 2: Cumulative Growth ---
    ax = axes[1,0]
    ax.plot(days+1,cum_pred,label='avg predicted')
    ax.plot(days+1,cum_pred_adj,label='avg predicted adjusted')
    ax.plot(days+1, cum_actual,label='VIIRS')
    ax.plot(days+1, cum_model,label='model')
    ax.set_title('Total Area of Fire (hectares)')
    ax.legend()
    ax.grid(True,alpha = 0.35)
    
    # --- Plot 3: Area Difference ---
    ax = axes[2,0]
    ax.plot(days+1,cum_model-cum_actual,label='model - VIIRS', color='tab:red')
    ax.plot(days+1,cum_model-cum_pred_adj,label='model - avg predicted_adj', color='tab:gray')
    ax.plot(days+1, cum_pred_adj-cum_actual,label='avg predicted_adj - VIIRS', color='tab:orange')
    ax.set_title(' Difference (hectares)')
    ax.legend()
    ax.grid(True,alpha = 0.35)
    ax.set_xlabel('Day')

    # --- Plot 4: Confusion Matrix ---
    ax = axes[0,1]
    tot_pix = CA.NUM_ROWS*CA.NUM_COLS
    ax.plot(hours/24, TP/tot_pix, label='True Positive')
    ax.plot(hours/24, FP/tot_pix, label='False Positive')
    ax.plot(hours/24, FN/tot_pix, label='False Negative')
    ax.set_title('Confusion Matrix')
    ax.legend()
    ax.grid(True,alpha = 0.35)

    # --- Plot 5: Accuracy ---
    ax = axes[1,1]
    ax.plot(hours/24, IOU, label='IOU')
    ax.plot(hours/24, TVER, label='Tversky')
    ax.plot(hours/24, SSIM, label='SSIM')
    ax.set_title('Accuracy')
    ax.legend()
    ax.grid(True,alpha = 0.35)
    
    # --- Plot 6: Num Burning ---
    ax = axes[2,1]
    ax.plot(hours/24, CA.num_burning)
    ax.set_title('Number Cells Burning')
    ax.grid(True,alpha = 0.35)
    ax.set_xlabel('Day')
    
    fig.suptitle(f'{CA.NAME} Fire burning for {CA.MAX_HOURS:0.0f} hours = {CA.MAX_HOURS/24:0.1f} days')
    plt.tight_layout(rect=[0, 0, 1, 0.98])

    plt.savefig(file_path, dpi=300, bbox_inches="tight")
    plt.close()
    '''
    
    
    
