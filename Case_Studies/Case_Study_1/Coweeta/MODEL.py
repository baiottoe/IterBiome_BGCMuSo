# IMPORT PACKAGES
import xarray as xr
import scipy
import pandas as pd
import re
from datetime import datetime, timedelta, date
import numpy as np
import matplotlib.pyplot as plt
from csv import reader
import math
import pyproj
from pyproj import CRS, Transformer, Proj, transform
import multiprocessing
from time import time
import os, shutil
import sys

# ---------------------------------------------------------------------------------------------------------------------------------
# PARAMETERS FROM COMMAND LINE
# !python MODEL.py {PROJECT_NAME} {CURRENT_DIRECTORY} {OUTPUT_TYPE} {CELL_SIZE}
# args: ['MODEL.py', 'Test', '/scratch/tmb4ee/MODEL_62', 'daily', '30']

args = list(sys.argv)
PROJECT_NAME = str(args[1])
CURRENT_DIRECTORY = str(args[2])
output_type = str(args[3])
CELLSIZE = int(args[4])

os.chdir(CURRENT_DIRECTORY)
# ---------------------------------------------------------------------------------------------------------------------------------
# IMPORTANT DIRECTORIES
CURRENT_PATH = os.getcwd()
PROJECT_DIR = os.path.join(CURRENT_PATH, PROJECT_NAME)
RAWGIS_DIR = os.path.join(PROJECT_DIR, "gis_data")
RAWOBS_DIR = os.path.join(PROJECT_DIR, "obs")
RAWSOIL_DIR = os.path.join(RAWGIS_DIR, "soil")
MODEL_DIR = os.path.join(PROJECT_DIR, 'model')
DEF_DIR = os.path.join(MODEL_DIR, 'defs')
INI_DIR = os.path.join(MODEL_DIR, 'ini_files')
EPC_DIR = os.path.join(MODEL_DIR, 'epc_files')
OUTPUT_DIR = os.path.join(MODEL_DIR, 'output')
CO2_DIR = os.path.join(MODEL_DIR, 'co2')
NDEP_DIR = os.path.join(MODEL_DIR, 'ndep')
ENDPOINT_DIR = os.path.join(MODEL_DIR, 'endpoint_files')
SPINUP_DIR = os.path.join(MODEL_DIR, 'spinup')
NORMAL_DIR = os.path.join(MODEL_DIR, 'normal')
MODEL_RAST_DIR = os.path.join(MODEL_DIR, 'raster_inputs')
IMAGE = os.path.join(PROJECT_DIR, 'image_map')

# ---------------------------------------------------------------------------------------------------------------------------------
# Give system permission to use executable
os.system(f"chmod u+x {os.path.join(MODEL_DIR, 'muso')}")
     
# ---------------------------------------------------------------------------------------------------------------------------------
# INTERNAL MODEL VARIABLES AND INFORMATION
# Coordinate System
COORDINATE_SYSTEM = CRS("+proj=lcc +ellps=WGS84 +a=6378137 +b=6356752.314245 +lat_1=25 +lat_2=60 +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +units=m +no_defs")
crs_wgs = CRS("WGS84")
transformer = Transformer.from_crs(COORDINATE_SYSTEM, crs_wgs)

# Dictionary indicating the file extension given based on input above
output_suffix = {'daily': '.dayout', 'monthly_agg': '.monavgout', 'annual_agg': '.annavgout'}


# ---------------------------------------------------------------------------------------------------------------------------------
# READ DEF AND INPUT FILES
# DEF FILES
# Read main.def
main_def = os.path.join(DEF_DIR, "main.def")
with open(main_def, 'r') as file:
    # read a list of lines into data
    main_def_lines = file.readlines()
main_d = {}
for i in main_def_lines:
    line = list(filter(None, re.split(r'\s', i)))
    if len(line) == 2:
        main_d[line[0]] = line[1]  

# Read lu_def
lu_def = str(DEF_DIR) + "/lu.def"
with open(lu_def, 'r') as file:
    # read a list of lines into data
    lu_def_lines = file.readlines()
lu_dict = {}
for i in lu_def_lines:
    line = list(filter(None, re.split(r'\s', i)))
    lu_dict[line[0]] = line[1:]

# EPC file info
epc_file_list = []
skip_first = 0
for i in lu_dict:
    if skip_first == 0:
        skip_first = 1
    else:
        if lu_dict[i][0] not in epc_file_list and lu_dict[i][0] != 'none':
            epc_file_list.append(lu_dict[i][0])  
        if lu_dict[i][1] not in epc_file_list and lu_dict[i][1] != 'none':
            epc_file_list.append(lu_dict[i][1])  
epc_file_dict = {}
for i in epc_file_list:
    with open(os.path.join(EPC_DIR, i), 'r') as file:
        epc = file.readlines()
        file.close()
    epc_file_dict[i] = epc
    
# Soil texture information
texture_mukey_dict = {}
with open(os.path.join(RAWSOIL_DIR, "soil_mukey_texture.csv"), 'r') as file:
    texture_mukey_raw = reader(file)
    texture_mukey_LOL = list(texture_mukey_raw)
    for i in texture_mukey_LOL:
        texture_mukey_dict[i[0]] = i[1:]
    file.close()

# Soil mukey information
soil_mukey_dict = {}
texture_dict = {}
with open(os.path.join(MODEL_DIR, "soil_mukey_rhessys.csv"), 'r') as file:
    soil_mukey_raw = reader(file)
    soil_mukey_LOL = list(soil_mukey_raw)
    for i in soil_mukey_LOL:
        soil_mukey_dict[i[1]] = i[2:]
        if i[2] not in texture_dict.keys():
            texture_dict[i[2]] = i[3:]
    file.close()

# Read output file
with open(os.path.join(DEF_DIR, "output.def"), 'r') as file:
    output_defs = file.readlines()
    file.close()
output_def = list(filter(None, output_defs))


# INPUT RASTER FILES
# Read All Necessary Files
with open(os.path.join(MODEL_RAST_DIR, "nlcd.txt"), 'r') as file:
    lu_gis_raw = file.readlines()
    file.close()

with open(os.path.join(MODEL_RAST_DIR, "dem.txt"), 'r') as file:
    dem_gis_raw = file.readlines()
    file.close()
    
with open(os.path.join(MODEL_RAST_DIR, "ymap.txt"), 'r') as file:
    lat_gis_raw = file.readlines()
    file.close()

with open(os.path.join(MODEL_RAST_DIR, "xmap.txt"), 'r') as file:
    long_gis_raw = file.readlines()  
    file.close()

with open(os.path.join(MODEL_RAST_DIR, "soil_texture.txt"), 'r') as file:
    soil_texture_gis_raw = file.readlines()  
    file.close()

# Read SSURGO GIS data
with open(os.path.join(MODEL_RAST_DIR, "ssurgo.txt"), 'r') as file:
    ssurgo_raw = file.readlines()
    file.close()
ssurgo_LOL = []
for i in ssurgo_raw:
    line = re.split(r'\s', i)
    line = list(filter(None, line))
    ssurgo_LOL.append(line)
    file.close()
    
# Read Daymet Cell GIS data
with open(os.path.join(MODEL_RAST_DIR, "daymet.txt"), 'r') as file:
    daymet_cells_raw = file.readlines()
    file.close()
    
# CLEAN FILES UP - MAKE INTO LIST OF LISTS WITH NO WHITESPACE
def clean_file(LOL):
    clean_LOL = []
    line_num = 0
    for i in LOL:
        clean_LOL.append(list(filter(None, re.split(r'\s', i))))
        line_num += 1
    return clean_LOL


lu_gis = clean_file(lu_gis_raw)
dem_gis = clean_file(dem_gis_raw)
lat_gis = clean_file(lat_gis_raw)
long_gis = clean_file(long_gis_raw)
soil_texture_gis = clean_file(soil_texture_gis_raw)
daymet_cells_gis = clean_file(daymet_cells_raw)
    
# ---------------------------------------------------------------------------------------------------------------------------------
# GENERATE INI FILES

#Generate s_ini
with open(os.path.join(INI_DIR, "ini_template.ini"), 'r') as file:
    s_ini = file.readlines()
    file.close()
    
# Restart Info
s_ini[8] = "0" + s_ini[8]
s_ini[9] = "1" + s_ini[9]

# Met Years
s_ini[14] = str(int(main_d["s_end_date"].split('-')[0]) - int(main_d["s_start_date"].split('-')[0]) + 1) + s_ini[14]
s_ini[15] = str(main_d["s_start_date"].split('-')[0]) + s_ini[15]

# Spinup?
s_ini[16] = "1" + s_ini[16]

# CO2 control
s_ini[20] = "0" + s_ini[20]
s_ini[22] = str(os.path.join(CO2_DIR, main_d["co2"])) + s_ini[22]

s_ini[25] = "0" + s_ini[25]
s_ini[27] = str(os.path.join(NDEP_DIR, main_d["ndep"])) + s_ini[27]

s_ini[48] = '1' + s_ini[48]

# Output prefix
# s_ini[83] = str(os.path.join(OUTPUT_DIR, PROJECT_NAME + "_S")) + s_ini[83]

# Output type
s_ini[106] = '2' + s_ini[106]
s_ini[107] = '0' + s_ini[107]
s_ini[108] = '0' + s_ini[108]

# Daily Outputs
ind_daily = 0
for output_line in output_def:
    s_ini.insert(113 + ind_daily, output_line)
    ind_daily += 1

# Generate n_ini
with open(os.path.join(INI_DIR, "ini_template.ini"), 'r') as file:
    n_ini = file.readlines()
    file.close()        
    
# Restart Info
n_ini[8] = "1" + n_ini[8]
n_ini[9] = "0" + n_ini[9]
# n_ini[11] = PROJECT_NAME + n_ini[11]
# n_ini[12] = PROJECT_NAME + n_ini[12]

# Met Years
n_ini[14] = str(int(main_d["n_end_date"].split('-')[0]) - int(main_d["n_start_date"].split('-')[0]) + 1) + n_ini[14]
n_ini[15] = str(main_d["n_start_date"].split('-')[0]) + n_ini[15]

# Spinup?
n_ini[16] = "0" + n_ini[16]

# CO2 control
n_ini[20] = "1" + n_ini[20]
n_ini[22] = str(os.path.join(CO2_DIR, main_d["co2"])) + n_ini[22]

n_ini[25] = "1" + n_ini[25]
n_ini[27] = str(os.path.join(NDEP_DIR, main_d["ndep"])) + n_ini[27]

# Output type
if output_type == 'daily':
    n_ini[106] = '2' + n_ini[106]
    n_ini[107] = '0' + n_ini[107]
    n_ini[108] = '0' + n_ini[108]
elif output_type == 'monthly_agg':
    n_ini[106] = '0' + n_ini[106]
    n_ini[107] = '2' + n_ini[107]
    n_ini[108] = '0' + n_ini[108]
elif output_type == 'annual_agg':
    n_ini[106] = '0' + n_ini[106]
    n_ini[107] = '0' + n_ini[107]
    n_ini[108] = '2' + n_ini[108]

# Daily Outputs
ind_daily = 0
for output_line in output_def:
    n_ini.insert(113 + ind_daily, output_line)
    ind_daily += 1
    
# ---------------------------------------------------------------------------------------------------------------------------------
# GENERATE SOIL MAP TEMPLATE FILE
with open(os.path.join(INI_DIR, "soil.SOI"), 'r') as file:
    soil_soi = file.readlines()
    file.close()
    
# ---------------------------------------------------------------------------------------------------------------------------------
# OPEN METEOROLOGICAL DATASET
met_data = xr.open_dataset(os.path.join(MODEL_DIR, 'met_data', 'aoi_met.nc'))
daymet_proj = Proj("+proj=lcc +ellps=WGS84 +a=6378137 +b=6356752.314245 +lat_1=25 +lat_2=60 +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +units=m +no_defs")

x_coords = list(met_data.coords['x'].values)
y_coords = list(met_data.coords['y'].values)

met_cells = str(DEF_DIR) + "/met_cells.csv"
with open(met_cells, 'r') as file:
    # read a list of lines into data
    met_cells_lines = file.readlines()

met_cells_dict = {}
for i in met_cells_lines:
    line = [int(float(x.strip())) for x in i.split(',')]
    met_cells_dict[line[0]] = tuple(line[1:]) 
    
# ---------------------------------------------------------------------------------------------------------------------------------
# CREATE PARALLELIZATION DICTIONARY
# Dictionary defines where in the input files each process is responsible for running and combining. This is based on boundaries of the Daymet meteorological input data cells. If the cell center for a simulation cell is within a Daymet cell, that is where it is assigned to.

north_bound = daymet_cells_gis[0][1]
south_bound = daymet_cells_gis[1][1]
east_bound = daymet_cells_gis[2][1]
west_bound = daymet_cells_gis[3][1]
rows = daymet_cells_gis[4][1]
columns = daymet_cells_gis[5][1]

west_cell_center = int(west_bound) + int((CELLSIZE/2))
east_cell_center = int(east_bound) + int((CELLSIZE/2))
north_cell_center = int(north_bound) - int((CELLSIZE/2))
south_cell_center = int(south_bound) - int((CELLSIZE/2))

west_daymet, east_daymet, north_daymet, south_daymet = (round((west_cell_center+250)/1000)*1000 - 250), (round((east_cell_center+250)/1000)*1000 - 250), (round((north_cell_center)/1000)*1000), (round((south_cell_center)/1000)*1000)

# print(west_cell_center, east_cell_center, north_cell_center, south_cell_center)
# print(west_daymet, east_daymet, north_daymet, south_daymet)

num_1km_columns = int((east_daymet - west_daymet)/1000) + 1
num_1km_rows = int((north_daymet - south_daymet)/1000) + 1


row_breaks = []
i_num = 0
for i in range(num_1km_rows):
    if i_num == 0:
        row_breaks.append([6, math.ceil(((abs(north_daymet) + 500)-abs(north_cell_center))/30) + 6])
    elif i_num != 0 and i_num != (num_1km_rows-1):
        row_breaks.append([row_breaks[-1][1] + 1, math.ceil(((abs(north_daymet) + 500 + 1000*i_num)-abs(north_cell_center))/30) + 6])
    elif i_num == (num_1km_rows - 1):
        row_breaks.append([row_breaks[-1][1] + 1, int(rows) + 5])
    i_num += 1
    
column_breaks = []
i_num = 0
for i in range(num_1km_columns):
    if i_num == 0:
        column_breaks.append([0, math.floor(((west_daymet + 500)-west_cell_center)/30)])
    elif i_num != 0 and i_num != (num_1km_columns-1):
        column_breaks.append([column_breaks[-1][1] + 1, math.floor(((west_daymet + 500 + 1000*i_num)-west_cell_center)/30)])
    elif i_num == (num_1km_columns-1):
        column_breaks.append([column_breaks[-1][1] + 1, int(columns) - 1])
    
    i_num += 1

for key, value in met_cells_dict.items():
    if value == (west_daymet, north_daymet):
        northwest_cell = key
    if value == (west_daymet, south_daymet):
        southwest_cell = key

row_step = (southwest_cell - northwest_cell)/(num_1km_rows - 1)

parallization_dict = {}
cell_value_list = []
i_row = 0
for i in range(num_1km_rows):
    i_column = 0
    row_list = []
    for j in range(num_1km_columns):
        parallization_dict[northwest_cell + i_row*row_step + i_column] = {'rows': row_breaks[i_row], 'columns': column_breaks[i_column], 'value': (northwest_cell + i_row*row_step + i_column), 'row_bounds': [int(north_bound)-(row_breaks[i_row][0]-6)*30 -15, int(north_bound)-(row_breaks[i_row][1]-6)*30 -15], 'column_bounds': [int(west_bound)+(column_breaks[i_column][0])*30 + 15, int(west_bound)+(column_breaks[i_column][1])*30 + 15]}
        row_list.append(int(northwest_cell + i_row*row_step + i_column))
        i_column += 1
    cell_value_list.append(row_list)
    i_row += 1

# print(west_daymet, east_daymet, north_daymet, south_daymet)

cell_values_df = pd.DataFrame(cell_value_list)
cell_values_df.to_csv(os.path.join(DEF_DIR, 'cell_values.csv'), header=False, index=False)

# ---------------------------------------------------------------------------------------------------------------------------------
# DECIDE HOW MANY PROCESSES TO ASSIGN TO MULTIPROCESSING
# Based on the number of daymet grid cells in which we are iterating over. If >40, then we use the maximum number allowable of 40. If under 40, we use the number of grid cells we are iterating over in the AOI for maximum model speed/efficiency.

if len(parallization_dict) >= multiprocessing.cpu_count():
    num_processes = multiprocessing.cpu_count()
else:
    num_processes = len(parallization_dict)


# ---------------------------------------------------------------------------------------------------------------------------------
# MODEL FUNCTIONS
def iter_cells(cell_data):
    patch_ID = int(cell_data['value'])
    no_data = '*'
    unique_num = 0
    
    patch_ID_list = {}
    ID_first = {}
    
    ncdf_long_list = ncdf_list()
    
    # Create output directory for this patch ID
    if os.path.exists(os.path.join(NORMAL_DIR, str(patch_ID))) is False:
        os.makedirs(os.path.join(NORMAL_DIR, str(patch_ID)))
    
    line_num = cell_data['rows'][0]
    for row in range(cell_data['rows'][0], cell_data['rows'][1]+1):
        column_num = cell_data['columns'][0]
        for j in range(cell_data['columns'][0], cell_data['columns'][1]+1):
            if dem_gis[line_num][column_num] != no_data and lu_gis[line_num][column_num] != no_data and lat_gis[line_num][column_num] != no_data and long_gis[line_num][column_num] != no_data and soil_texture_gis[line_num][column_num] != no_data and daymet_cells_gis[line_num][column_num] != no_data:
                x_daymet, y_daymet = met_cells_dict[patch_ID][0], met_cells_dict[patch_ID][1]
                dem_midpoint = math.floor(float(dem_gis[line_num][column_num])/200)*200 + 100
                characteristic_tuple = (x_daymet, y_daymet, dem_midpoint, (lu_dict[lu_gis[line_num][column_num]][0]).split('.')[0], str(soil_texture_gis[line_num][column_num]))
                if characteristic_tuple not in patch_ID_list.keys():
                    # Add site info to ID list to prevent future duplicate spinups
                    patch_ID_list[characteristic_tuple] = f"{patch_ID}_{unique_num}"
                    ID_first[f"{patch_ID}_{unique_num}"] = characteristic_tuple

                    # Spinup run
                    try:
                        mean_annual_t, mean_annual_t_range = spinup(f"{patch_ID}_{unique_num}", ID_first)
                    except:
                        print(f"Spinup of {characteristic_tuple} failed")

                    unique_num += 1

                # Normal Run
                if lu_dict[lu_gis[line_num][column_num]][1] != "none":
                    x_daymet, y_daymet = met_cells_dict[patch_ID][0], met_cells_dict[patch_ID][1]
                    dem_midpoint = math.floor(float(dem_gis[line_num][column_num])/200)*200 + 100
                    characteristic_tuple = (x_daymet, y_daymet, dem_midpoint, (lu_dict[lu_gis[line_num][column_num]][0]).split('.')[0], str(soil_texture_gis[line_num][column_num]))

                    try:
                        change_setup(n_ini, line_num, column_num, 'n', patch_ID_list[characteristic_tuple], patch_ID, mean_annual_t, mean_annual_t_range)

                        run_model('n', patch_ID)
                    except:
                        pass
                        
            column_num += 1
        line_num += 1
    print(f"DONE WITH {patch_ID}!!!!!!!!!!!!!!!")

    ds_output = setup_output(cell_data, ncdf_long_list)
    ds_output = update_output(ds_output, patch_ID)
    
def get_met(long, lat):
    daymet_pt_data = met_data.sel(x=long, y=lat)
    daymet_pt_data = daymet_pt_data.to_dataframe()
    daymet_pt_data['year'] = pd.DatetimeIndex(daymet_pt_data.index).year
    daymet_pt_data['yday'] = daymet_pt_data.index.dayofyear
    daymet_pt_data['prcp'] = daymet_pt_data['prcp']/10
    daymet_pt_data['tday'] = 0.606*daymet_pt_data['tmax'] + 0.394*daymet_pt_data['tmin']
    daymet_pt_data['vpsat'] = 0.61094 * 2.71828**((17.725*daymet_pt_data['tday'])/(daymet_pt_data['tday']+243.04))*1000
    daymet_pt_data['vpd'] = abs(daymet_pt_data['vpsat'] - daymet_pt_data['vp'])
    daymet_final = daymet_pt_data[['year', 'yday', 'tmax', 'tmin', 'tday', 'prcp', 'vpd', 'srad', 'dayl']].copy()
    daymet_final.to_csv(os.path.join(MODEL_DIR, 'met_data', f'{os.getpid()}_metdata.mtc43'), index = False, sep="\t")
    
    t_mean_ann_avg = round(daymet_pt_data['tday'].mean(), 2)
    t_range_ann_avg = round(max(daymet_pt_data.groupby(daymet_pt_data.index.year)['tday'].transform('mean').unique()) - min(daymet_pt_data.groupby(daymet_pt_data.index.year)['tday'].transform('mean').unique()), 2)
    
    return t_mean_ann_avg, t_range_ann_avg

def spinup(ID, ID_first):

    if ID_first[ID][3] != "none":
        ini_file = s_ini.copy()
        soil_file = soil_soi.copy()
        
        # Change met data file
        mean_annual_t, mean_annual_t_range = get_met(int(ID_first[ID][0]), int(ID_first[ID][1]))
        ini_file[3] = str(os.path.join(MODEL_DIR, 'met_data', f'{os.getpid()}_metdata.mtc43')) + ini_file[3]
        
        # Endpoint File
        ini_file[10] = os.path.join(ENDPOINT_DIR, str(ID)) + ini_file[10]
        ini_file[11] = os.path.join(ENDPOINT_DIR, str(ID)) + ini_file[11]

        # Change Site Elevation
        ini_file[30] = str(ID_first[ID][2]) + ini_file[30]

        # Change Site Latitude
        ini_lat, ini_lon = transform(daymet_proj, crs_wgs, int(ID_first[ID][0]), int(ID_first[ID][1]))
        ini_file[31] = str(ini_lat) + ini_file[31]
        
        # Change Site Mean Annual Temp
        ini_file[33] = str(mean_annual_t) + ini_file[33]
        ini_file[34] = str(mean_annual_t_range) + ini_file[34]
        
        # Change Management File, Model does not simulate management in spinup anyways
        ini_file[44] = 'none' + ini_file[44]
        
        ini_file[105] = str(os.path.join(OUTPUT_DIR, str(ID) + "_S")) + ini_file[105]

        # Change Soil Info
        silt = round(float(texture_dict[ID_first[ID][4]][6])*100, 1)
        sand = round(float(texture_dict[ID_first[ID][4]][4])*100, 1)
        rooting_depth = texture_dict[ID_first[ID][4]][3]
        # sand
        line_pre = str(sand) + " " + str(sand) + " " + str(sand) + " " + str(sand) + " " + str(sand) + " " + str(sand) + " " + str(sand) + " " + str(sand) + " " + str(sand) + " " + str(sand)
        soil_file[81] = line_pre + soil_file[81]
        # silt
        line_pre = str(silt) + " " + str(silt) + " " + str(silt) + " " + str(silt) + " " + str(silt) + " " + str(silt) + " " + str(silt) + " " + str(silt) + " " + str(silt) + " " + str(silt)
        soil_file[82] = line_pre + soil_file[82]
        
        # Change soil File
        soil_path = str(os.path.join(SPINUP_DIR, f"{os.getpid()}soil"))
        ini_file[38] = soil_path + ini_file[38]
        
        
        # Change EPC File
        epc = ID_first[ID][3] + '.epc'
        epc_path = str(os.path.join(SPINUP_DIR, f"{os.getpid()}{epc}"))
        ini_file[41] = epc_path + ini_file[41]

        epc_file = epc_file_dict[epc].copy()

        # EPC - Change rooting depth
        epc_file[63] = str(rooting_depth) + epc_file[63]

        with open(str(os.path.join(SPINUP_DIR, f"{os.getpid()}{epc}")), 'w') as file:
            for item in epc_file:
                file.write(item)
        
        with open(str(os.path.join(SPINUP_DIR, f"{os.getpid()}soil")), 'w') as file:
            for item in soil_file:
                file.write(item)

        with open(str(os.path.join(SPINUP_DIR, f"s{os.getpid()}.ini")), 'w') as file:
            for item in ini_file:
                file.write(item)
        
        ini_file_name = f"s{os.getpid()}.ini"
        os.system(f"{os.path.join(MODEL_DIR, 'muso')} {os.path.join(SPINUP_DIR, ini_file_name)} > /dev/null")
        
        return mean_annual_t, mean_annual_t_range
        
#         os.remove(os.path.join(SPINUP_DIR, f"{os.getpid()}{epc}"))
#         os.remove(os.path.join(SPINUP_DIR, f"{os.getpid()}{soil}"))
#         os.remove(os.path.join(SPINUP_DIR, f"s{os.getpid()}.ini"))
    
def ncdf_list():
    ncdf_long_list = []
    step = (int(long_gis[2][1]) - int(long_gis[3][1]))/int(long_gis[5][1])
    for i in range(int(long_gis[5][1])):
        ncdf_long_list.append(str(int(int(long_gis[3][1]) + i*step + step/2)))
    return(ncdf_long_list)   

def change_setup(ini_template, line_num, column_num, Type, ID, patch_ID, mean_annual_t, mean_annual_t_range):
    ini_file = ini_template.copy()
    soil_file = soil_soi.copy()
    
    # Change met data
    ini_file[3] = str(os.path.join(MODEL_DIR, 'met_data', f'{os.getpid()}_metdata.mtc43')) + ini_file[3]
    
    # Change endpoint file
    ini_file[10] = os.path.join(ENDPOINT_DIR, str(ID)) + ini_file[10]
    ini_file[11] = os.path.join(ENDPOINT_DIR, str(ID)) + ini_file[11]
    
    # Change Site Elevation
    ini_file[30] = dem_gis[line_num][column_num] + ini_file[30]
    
    # Change Site Latitude
    ini_file[31] = str(round(transformer.transform(long_gis[line_num][column_num], lat_gis[line_num][column_num])[0], 4)) + ini_file[31]
    
    # Change Site Mean Annual Temp
    ini_file[33] = str(mean_annual_t) + ini_file[33]
    ini_file[34] = str(mean_annual_t_range) + ini_file[34]
    
    # Change Soil Info
    clay, silt, sand, rooting_depth = update_soil(line_num, column_num)
    # sand
    line_pre = str(sand[0]) + " " + str(sand[1]) + " " + str(sand[2]) + " " + str(sand[3]) + " " + str(sand[4]) + " " + str(sand[5]) + " " + str(sand[6]) + " " + str(sand[7]) + " " + str(sand[8]) + " " + str(sand[9])
    soil_file[81] = line_pre + soil_file[81]
    # silt
    line_pre = str(silt[0]) + " " + str(silt[1]) + " " + str(silt[2]) + " " + str(silt[3]) + " " + str(silt[4]) + " " + str(silt[5]) + " " + str(silt[6]) + " " + str(silt[7]) + " " + str(silt[8]) + " " + str(silt[9])
    soil_file[82] = line_pre + soil_file[82]
    
    # Change soil File
    soil_path = str(os.path.join(NORMAL_DIR, str(patch_ID), f"{os.getpid()}soil"))
    ini_file[38] = soil_path + ini_file[38]
    
    # Change EPC File
    epc = lu_dict[lu_gis[line_num][column_num]][1]
    epc_path = str(os.path.join(NORMAL_DIR, str(patch_ID), f"{os.getpid()}{epc}"))
    ini_file[41] = epc_path + ini_file[41]
    
    epc_file = epc_file_dict[epc].copy()
    
    # Change Management File
    if lu_dict[lu_gis[line_num][column_num]][2] != "none":
        ini_file[44] = str(os.path.join(MODEL_DIR, 'mgm_files', lu_dict[lu_gis[line_num][column_num]][2])) + ini_file[44]
        
        with open(str(os.path.join(MODEL_DIR, 'mgm_files', lu_dict[lu_gis[line_num][column_num]][2])), 'r') as file:
            mgm_iter = 0
            mgm_line_num = 9999
            for line in file:
                if line.strip('\s\n') == 'PLANTING':
                    mgm_line_num = mgm_iter
                if mgm_iter == mgm_line_num+1:
                    if list(filter(None, re.split(r'\s', line)))[0] == '1':
                        ini_file[48] = '0' + ini_file[48]
                    else:
                        ini_file[48] = '1' + ini_file[48]
                mgm_iter += 1
    else:
        ini_file[44] = 'none' + ini_file[44]
        ini_file[48] = '1' + ini_file[48]

    # EPC - Change rooting depth
    epc_file[63] = str(rooting_depth) + epc_file[63]
    
    # Change Output Prefix
    ini_file[105] = str(os.path.join(NORMAL_DIR, str(patch_ID), f"{lat_gis[line_num][column_num]}_{long_gis[line_num][column_num]}")) + ini_file[105]
    
    # Create files    
    with open(str(os.path.join(NORMAL_DIR, str(patch_ID), f"{os.getpid()}{epc}")), 'w') as file:
        for item in epc_file:
            file.write(item)
            
    with open(str(os.path.join(NORMAL_DIR, str(patch_ID), f"{os.getpid()}soil")), 'w') as file:
        for item in soil_file:
            file.write(item)
                
    with open(str(os.path.join(NORMAL_DIR, str(patch_ID), f"{os.getpid()}{Type}.ini")), 'w') as file:
        for item in ini_file:
            file.write(item)
    
def update_soil(line_num, column_num):
    depth_start = [0, 2, 5, 10, 20, 50, 100, 150, 200, 400] # cm
    depth_end = [2, 5, 10, 20, 50, 100, 150, 200, 400, 1000] # cm
    mean_depth = [1, 3.5, 7.5, 15, 35, 75, 125, 175, 300, 700] # cm

    mukey = ssurgo_LOL[line_num][column_num]

    clay = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    sand = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    silt = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    if texture_mukey_dict[mukey][1] != "NA":
        depth = 0
        horizons = 0
        horizon_thickness_list = [1, 16, 31]
        for horizon in horizon_thickness_list:
            if texture_mukey_dict[mukey][horizon] != "NA":
                horizons += 1

        horizon_num = 0
        while horizon_num <= (horizons-1):
            top_depth = float(texture_mukey_dict[mukey][int(horizon_num*15 + 3)]) * 100 # cm
            bottom_depth = float(texture_mukey_dict[mukey][int(horizon_num*15 + 4)]) * 100 # cm
            for depth in mean_depth:
                if top_depth <= depth < bottom_depth:
                    sand[mean_depth.index(depth)] = round(float(texture_mukey_dict[mukey][int(horizon_num*15 + 5)])*100, 1) # %
                    silt[mean_depth.index(depth)] = round(float(texture_mukey_dict[mukey][int(horizon_num*15 + 6)])*100, 1) # %
                    clay[mean_depth.index(depth)] = round(100 - sand[mean_depth.index(depth)] - silt[mean_depth.index(depth)], 1)
            horizon_num += 1

    for depth in mean_depth:
        if sand[mean_depth.index(depth)] + silt[mean_depth.index(depth)] + clay[mean_depth.index(depth)] != 1:
            sand[mean_depth.index(depth)] = round(float(soil_mukey_dict[mukey][5])*100, 1)
            silt[mean_depth.index(depth)] = round(float(soil_mukey_dict[mukey][7])*100, 1)
            clay[mean_depth.index(depth)] = round((100 - sand[mean_depth.index(depth)] - silt[mean_depth.index(depth)]), 1)

    #Is this the best way????
    soil_depth = soil_mukey_dict[mukey][1]   
    rooting_depth = soil_mukey_dict[mukey][4]

    return clay, silt, sand, rooting_depth

def run_model(Type, patch_ID):
    ini_file_name = f"{os.getpid()}{Type}.ini"
    os.system(f"{os.path.join(MODEL_DIR, 'muso')} {os.path.join(NORMAL_DIR, str(patch_ID), ini_file_name)} > /dev/null")
      
def setup_output(cell_data, ncdf_long_list):
    # Generate list of dates, prep for xarray coords
    if output_type == 'daily':
        dr = pd.date_range(start=main_d["n_start_date"], end=main_d["n_end_date"], freq='D')
        dates_list = dr[(dr.day != 29) | (dr.month != 2)]
    elif output_type == 'monthly_agg':
        dates_list = pd.date_range(start=main_d["n_start_date"], end=main_d["n_end_date"], freq='MS')
    elif output_type == 'annual_agg':
        dates_list = pd.date_range(start=main_d["n_start_date"], end=main_d["n_end_date"], freq='Y')

    # Read ini file, utilize output variable names in xarray
    output_variables = ['year', 'month', 'day', 'yday']
    number_outputs = int(output_def[0].split(' ', 1)[0])
    current_output = 1
    while current_output <= number_outputs:
        var = list(filter(None, re.split(r'\s', output_def[current_output])))[1]
        output_variables.append(var)
        current_output += 1
    
    # Latitude and Longitude list for uniform dimensions
    latitude_list = list(reversed(range(cell_data["row_bounds"][1], cell_data["row_bounds"][0]+30, 30)))
    longitude_list = list(range(cell_data["column_bounds"][0], cell_data["column_bounds"][1]+30, 30))

    # Create dict for data in xarray
    data_dict = {}
    na_dict = {}
    for i in output_variables:
        empty = np.empty((len(latitude_list), len(longitude_list), len(dates_list)))
        empty.fill(np.nan)
        data_dict[str(i)] = (["lat", "long", "time"], empty)

    # Generate xarray dataset with coords and data, but empty values
    ds_output = xr.Dataset(data_dict, coords={"time": dates_list, "long": longitude_list, "lat": latitude_list})
    
    return ds_output
        
def update_output(output_empty, patch_ID):
    ds_run = output_empty
    
    patch_output_files = [f for f in os.listdir(os.path.join(NORMAL_DIR, str(patch_ID))) if f.endswith(f"{output_suffix[output_type]}")]
    
    for model_file in patch_output_files:
        try:
            # Change values for first cell
            with open(os.path.join(NORMAL_DIR, str(patch_ID), model_file), 'r') as file:
                n_out_lines = file.readlines()
                file.close()

            output_LOL = []
            for i in n_out_lines:
                line = re.split(r'\s', i)
                line = list(filter(None, line))
                output_LOL.append(line)

            output_array = np.array(output_LOL)
            transposed_output = output_array.T

            for i in transposed_output:
                ds_run[str(i[0])].loc[dict(lat=int(model_file.split('.')[0].split('_')[0]), long=int(model_file.split('.')[0].split('_')[1]))] = i[1:]
        except:
            print(f"ERROR with {model_file}")

    ds_run.to_netcdf(os.path.join(OUTPUT_DIR, str(patch_ID) + ".nc"), unlimited_dims=["lat"])
    
    print(f"OUTPUT {str(patch_ID)}.nc DONE")
    
    # Open a file with access mode 'a'
    gridcell_complete_file = open(os.path.join(INI_DIR, 'gridcells_complete.txt'), 'a')
    # Append 'hello' at the end of file
    gridcell_complete_file.write(f'{patch_ID}\n')
    # Close the file
    gridcell_complete_file.close()
    
    os.system(f"rm -r {os.path.join(NORMAL_DIR, str(patch_ID))}")
    
    return ds_run

def main(parallization_dict):
    start = time()
    line_num = 0
    nodata_value = '*'
    c1 = time()
    
    if os.path.exists(os.path.join(INI_DIR, 'gridcells_complete.txt')) is False:
        # Creating a file at specified location
        with open(os.path.join(INI_DIR, 'gridcells_complete.txt'), 'w') as fp:
            pass
        
        pool = multiprocessing.Pool(processes=num_processes)
        normal_map = pool.map(iter_cells, parallization_dict.values())
        # for item in parallization_dict.values():
        #     iter_cells(item)
        
    else:
        with open(os.path.join(INI_DIR, 'gridcells_complete.txt'), 'r') as file:
            # read a list of lines into data
            completed_ids = file.readlines()
        try:
            completed_ids.remove('\n')
        except:
            pass

        num = 0
        for i in completed_ids:
            completed_ids[num] = round(float(i.strip()), 1)
            num += 1

        incomplete_gridcells = []
        for i in parallization_dict:
            if i not in completed_ids:
                incomplete_gridcells.append(parallization_dict[i])
        
        pool = multiprocessing.Pool(processes=num_processes)
        normal_map = pool.map(iter_cells, incomplete_gridcells)
        
#     normal_map = pool.map(iter_cells, error_dict.values())

    
    c2 = time()
    print(f"Iterative Model Run: {c2 - c1}")
    

# ---------------------------------------------------------------------------------------------------------------------------------
# RUN MODEL  
main(parallization_dict)