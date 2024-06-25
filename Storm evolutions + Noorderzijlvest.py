# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 08:30:16 2024

@author: jilde
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def generate_dataset_NZV(timestep_hours=0.1, duration_hours=35, water_level_min=-0.1, 
                         water_level_max=2.72, wave_height_min=0.1, wave_height_max=2.65, Tp_peak=5.75):
    # Calculate the total number of steps
    total_steps = int(duration_hours / timestep_hours) + 1
    
    # Create time arrays
    time_hours = np.linspace(0, duration_hours, total_steps)
    time_seconds = time_hours * 3600
    
    # Initialize the water level array
    water_level = np.zeros(total_steps)
    wave_height = np.zeros(total_steps)
    
    # Define middle duration for peak
    peak_start = 15.5
    peak_end = 19.5
    peak_start_idx = int(peak_start / timestep_hours)
    peak_end_idx = int(peak_end / timestep_hours)

    # Linearly interpolate water level
    water_level[:peak_start_idx] = np.linspace(water_level_min, water_level_max, peak_start_idx)
    water_level[peak_start_idx:peak_end_idx] = water_level_max
    water_level[peak_end_idx:] = np.linspace(water_level_max, water_level_min, total_steps - peak_end_idx)

    # Linearly interpolate significant wave height
    wave_height[:peak_start_idx] = np.linspace(wave_height_min, wave_height_max, peak_start_idx)
    wave_height[peak_start_idx:peak_end_idx] = wave_height_max
    wave_height[peak_end_idx:] = np.linspace(wave_height_max, wave_height_min, total_steps - peak_end_idx)

    # Calculate peak period Tp based on the significant wave height Hs
    Tp = np.zeros(total_steps)
    for i in range(1, total_steps):
        Tp[i] = (wave_height[i] * (Tp_peak**2 / wave_height_max))**0.5
    Tp[0] = Tp[1]  # Initial peak period
    
    # Define constants for erosion volume calculation
    m2 = 1.0
    s_op = 0.05
    height_berm = 2.52
    B_berm = 3.0
    tan_alpha = 0.333
    a1 = np.radians(18.4)  # Original slope angle in radians
    a2 = np.radians(5.7)   # Terrace angle in radians
    a3 = np.radians(63.4)  # Cliff angle in radians
    
    # Initialize the erosion volume array
    Ve_NZV = np.zeros(total_steps)
    dt = np.zeros(total_steps)
    do = np.zeros(total_steps)
    de = np.zeros(total_steps)
    
    # Calculate Ve_NZV for each time step starting from 15.12 hours
    # start_time = 15.12
    start_time = 0
    
    start_time_idx = int(start_time / timestep_hours)
    for i in range(start_time_idx, total_steps):
        h_overgang = water_level[i] - height_berm
        Hs = wave_height[i]
        
        # Check conditions
        condition1 = (height_berm / Hs) > -0.05
        condition2 = (h_overgang / Hs) > -0.7
        condition3 = Hs <= 2.8
        
        if condition1 and condition2 and condition3:
            # Calculate erosion volume
            erosion_factor = (0.068 * m2 * ((Hs - 0.25)**2 / np.sqrt(s_op)) *
                              min(2.4, max(0, 1.2 + 1.6 * (h_overgang / Hs))) *
                              (13 / ((B_berm / Hs) * min(1, max(0, 2 - height_berm / Hs)) + 8) - 0.5) *
                              (1.69 - 0.17 / tan_alpha))
            Ve_NZV[i] = Ve_NZV[i - 1] + timestep_hours * erosion_factor
        else:
            Ve_NZV[i] = 0
        
        # Calculate dt
        dt[i] = min(0.4 * Ve_NZV[i]**0.25 / Hs**1.5 + 0.7, 2)
        
        # Calculate do
        do[i] = (Hs * max(0, min(h_overgang / Hs + 1.5, 1)) * 
                 max(0, min((dt[i] - h_overgang) / (Hs * np.sin(a1)) * np.tan(a1 - a2), 0.18)))
        
        # Calculate de
        de[i] = np.sqrt((2 * Ve_NZV[i] * np.tan(a1 - a2) + do[i]**2) / 
                        (1 + (np.tan(a1 - a2) / np.tan(a3 - a1))))
    
    # Create DataFrame
    df = pd.DataFrame({
        'time (seconds)': time_seconds,
        'time (hours)': time_hours,
        'water level h (m+NAP)': water_level,
        'sign. wave height (in m)': wave_height,
        'peak period (sec.)': Tp,
        'Ve_NZV': Ve_NZV,
        'dt_NZV': dt,
        'd0_NZV': do,
        'de_NZV': de
    })
    
    return df

# Generate the dataset
dataset_NZV = generate_dataset_NZV()

# Filter out rows where Ve_NZV is zero or constant over three consecutive time steps
Ve_NZV = dataset_NZV['Ve_NZV']
mask_zero = Ve_NZV != 0
mask_constant = ~((Ve_NZV.shift(1) == Ve_NZV) & (Ve_NZV.shift(-1) == Ve_NZV))

# Combine masks
mask = mask_zero & mask_constant

# Apply mask to DataFrame
filtered_data = dataset_NZV[mask]

# Plot Ve_NZV against time
plt.figure(figsize=(10, 5))
plt.plot(filtered_data['time (hours)'], filtered_data['Ve_NZV'], label='Ve_NZV', color='purple')
plt.xlabel('Time (hours)')
plt.ylabel('Erosion volume (in m3/m/h)')
plt.title('Cumulative erosion volume (Klein Breteler, 2022) at Ketelmeerdijk HM 13.4')
#plt.legend()
plt.xlim(10,35)
plt.ylim(bottom=0)
plt.grid(True)
plt.show()

# Plot de against time
plt.figure(figsize=(10, 5))
plt.plot(filtered_data['time (hours)'], filtered_data['de_NZV'], label='de', color='purple')
plt.xlabel('Time (hours)')
plt.ylabel('Cumulative erosion volume (in m3/m)')
plt.title('Erosion depth development following Klein Breteler (2022) at Ketelmeerdijk HM 13.4')
plt.grid(True)
plt.show()

# Find the index of timestep T=15.5 hours
index_14_5 = int(15.5 / dataset_NZV['time (hours)'][1])

# Find the index of timestep T=19.5 hours
index_20_5 = int(19.5 / dataset_NZV['time (hours)'][1])

# Calculate the erosion volume difference between timestep T=19.5 hours and T=15.5 hours
erosion_volume_difference = dataset_NZV['Ve_NZV'][index_20_5] - dataset_NZV['Ve_NZV'][index_14_5]

# Calculate the erosion depth (de) at timestep T=15.5 hours
de_14_5 = dataset_NZV['de_NZV'][index_14_5]

# Calculate the erosion depth (de) at timestep T=19.5 hours
de_20_5 = dataset_NZV['de_NZV'][index_20_5]

# Calculate the value of de related to the erosion volume difference
de_related_to_erosion_difference = de_20_5 - de_14_5

# Save to CSV
dataset_NZV.to_csv('dataset_NZV.csv', index=False)


print("Erosion volume at 4-hour peak is", "%.2f" % erosion_volume_difference, "m3")
print("Erosion depth at 4-hour peak is", "%.2f" % de_related_to_erosion_difference, "m")
