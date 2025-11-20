import os
import pandas as pd
from datetime import datetime

root_dir = '/mnt/DataDrive1/data/zholden/VIIRS/CA_model_test'
out_file = root_dir + '/fire_info.csv'
res = 90

fireDF = pd.read_csv(f'{root_dir}/fire_list.csv')
new_cols = ['fire_start','first_supression','days_til_sup','hours_of_data','days_of_data']

for i,row in fireDF.iterrows():
	fire = row.fire_name
	print(fire)
	summary_file = f'{root_dir}/lw_model_outputs/{fire}_{res}m/combo_50_1_10_10_1_25_0.25_0.1/{fire}_simulation_info.txt'
	with open(summary_file, "r") as file:
		text = file.readlines()
		
	start = text[8].split(': ',1)[1].strip().split('.')[0]
	start = datetime.strptime(start,'%Y-%m-%d %H:%M:%S')
	
	sup = text[9].split(': ',1)[1].strip().split('.')[0]
	sup = datetime.strptime(sup,'%Y-%m-%d %H:%M:%S')
	
	days_til_sup = int(float(text[4].split('=')[1].split('days')[0].strip()))
	
	hours_of_data = int(text[10].split('hours')[0].split('for')[1].strip())
	days_of_data = int(text[10].split('=')[1].split('days')[0].strip())
	
	fireDF.loc[i,new_cols] = [start, sup, days_til_sup, hours_of_data, days_of_data]
    
fireDF.to_csv(out_file, index=False)
