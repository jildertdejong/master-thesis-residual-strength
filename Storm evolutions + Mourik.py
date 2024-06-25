# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 07:50:51 2024

@author: jilde
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def generate_hydraulic_table(timestep_hours=0.001, duration_hours=35, water_level_min=-0.1, 
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
    
    # Create DataFrame
    df = pd.DataFrame({
        'time (seconds)': time_seconds,
        'time (hours)': time_hours,
        'water level h (m+NAP)': water_level,
        'sign. wave height (in m)': wave_height,
        'peak period (sec.)': Tp
    })
    
    return df

# Generate the dataset
hydraulic_table = generate_hydraulic_table()

# Save to CSV
hydraulic_table.to_csv('hydraulic_conditions_Mourik.csv', index=False)

# Function to add the Ve_Mourik, dt, and de columns with the specified conditions
def add_ve_mourik_and_dt_columns(df, Ce=0.44, alpha=0.333, s_op=0.05, start_time=0):
    timestep_hours = (df['time (hours)'][1] - df['time (hours)'][0])
    start_index = int(start_time / timestep_hours)
    Ve_Mourik = np.zeros(len(df))
    dt_Mourik = np.zeros(len(df))
    de_Mourik = np.zeros(len(df))
    
    for i in range(1, len(df)):
        if i < start_index:
            Ve_Mourik[i] = 0
            dt_Mourik[i] = 0
            de_Mourik[i] = 0
            continue
        
        Hs = df['sign. wave height (in m)'][i]
        if Hs <= 0.4:
            Ve_Mourik[i] = 0
        else:
            delta_Ve = (timestep_hours * Ce * 
                        (1.32 - 0.079 * (Ve_Mourik[i-1] / Hs**2)) *
                        (16.4 * np.tan(alpha) / Hs**2) *
                        min(3.6, 0.0061 / s_op**1.5) *
                        (1.7 * (Hs - 0.4)**2))
            if delta_Ve < 0:
                delta_Ve = 0
            Ve_Mourik[i] = Ve_Mourik[i-1] + delta_Ve
            if Ve_Mourik[i] < Ve_Mourik[i-1]:
                Ve_Mourik[i] = Ve_Mourik[i-1]

        # Calculate dt using the given formula
        dt_Mourik[i] = min(0.4 * Ve_Mourik[i]**0.25 / np.sqrt(Hs) + 0.7, 2 * Hs)

        # Calculate de using the given formula (if Ve_Mourik[i] >= 0.75)
        if Ve_Mourik[i] >= 0.75:
            de_Mourik[i] = np.sqrt(Ve_Mourik[i] * np.tan(alpha)) - 0.14
    
    df['Ve_Mourik'] = Ve_Mourik
    df['dt_Mourik'] = dt_Mourik
    df['de_Mourik'] = de_Mourik
    return df

# Function to plot Ve_Mourik against time
def plot_ve_mourik(df):
    # Filter out rows where Ve_Mourik is zero or constant
    mask = (df['Ve_Mourik'] != 0) & (df['Ve_Mourik'] != df['Ve_Mourik'].shift(1)) & (df['Ve_Mourik'] != df['Ve_Mourik'].shift(-1))
    filtered_df = df[mask]
    
    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(filtered_df['time (hours)'], filtered_df['Ve_Mourik'], linestyle='-', color='green')
    plt.xlabel('Time (in hours)')
    plt.ylabel('Erosion volume (in m3/m)')
    plt.title('Cumulative erosion volume (Mourik, 2020) for Ketelmeerdijk HM13.4')
    plt.grid(True)
    plt.xlim(10,35)
    plt.ylim(bottom=0)
    plt.show()

    # Plot de against time
    plt.figure(figsize=(10, 5))
    plt.plot(filtered_df['time (hours)'], filtered_df['de_Mourik'], label='de', color='purple')
    plt.xlabel('Time (hours)')
    plt.ylabel('Erosion depth (in meters)')
    plt.title('Erosion depth development (Mourik, 2020) at Ketelmeerdijk HM 13.4')
    plt.grid(True)
    plt.show()

# Generate the original dataset
hydraulic_table = generate_hydraulic_table()

# Copy the original dataset
dataset_Mourik = hydraulic_table.copy()

# Add the Ve_Mourik, dt, and de columns to the dataset_Mourik
dataset_Mourik = add_ve_mourik_and_dt_columns(dataset_Mourik)

# Find the index of timestep T=15.5 hours
index_15_5 = int(15.5 / dataset_Mourik['time (hours)'][1])

# Find the index of timestep T=19.5 hours
index_19_5 = int(19.5 / dataset_Mourik['time (hours)'][1])

# Calculate the erosion volume difference between timestep T=19.5 hours and T=15.5 hours
erosion_volume_difference = dataset_Mourik['Ve_Mourik'][index_19_5] - dataset_Mourik['Ve_Mourik'][index_15_5]

# Calculate the erosion depth (de) at timestep T=15.5 hours
de_15_5 = dataset_Mourik['de_Mourik'][index_15_5]

# Calculate the erosion depth (de) at timestep T=19.5 hours
de_19_5 = dataset_Mourik['de_Mourik'][index_19_5]

# Calculate the value of de related to the erosion volume difference
de_related_to_erosion_difference = de_19_5 - de_15_5

print("Erosion volume at 4-hour peak is", "%.2f" % erosion_volume_difference, "m3")
print("Erosion depth at 4-hour peak is", "%.2f" % de_related_to_erosion_difference, "m")

dataset_Mourik.to_csv('dataset_Mourik.csv', index=False, sep=',', decimal='.')

# Display the plot
plot_ve_mourik(dataset_Mourik)
