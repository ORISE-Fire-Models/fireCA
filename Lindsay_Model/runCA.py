import CellAuto as ca
import itertools
import pandas as pd

def main():
	root_dir = '/mnt/DataDrive1/data/zholden/VIIRS/CA_model_test'
	# open DF with info for all fires
	fireDF = pd.read_csv(f'{root_dir}/fire_info.csv')
	# limit DF to fires that burn longer than 24 hours before suppression activities
	fireDF = fireDF.loc[fireDF.hours_til_first > 24]

	# for each fire, create and run a simulation          
	for i,row in fireDF.iterrows():

		print(f'{row.fire_name} {row.fireID}')

		CA = ca.CellAuto(root_dir, # sets the path to the current directory
						 row.fire_name, # sets the fire name
						 row.fireID, # sets the fire ID
						 out_dir = 'no_perim_update', # adds a model name within fire folder
						 num_plots = row.hours_til_first//12, # indicates the number of arrival plots (this value does one a day)
						 has_barriers = False, # indicates whether or not to use barriers to fire
						 max_hours = row.hours_til_first, # max hours the simulation will run (set to hours until first suppression)
						 res = 90, # spatial resolution in meters
						 burn_out_time = 72, # number of hours a cell is allowed to burn before it burns out
			)
		
		CA.run()

if __name__ == "__main__":
    main() 
