###### Gilad Tsfaty #####
######## GEO 2 NED ##########


##### USER INPUT #####

# ### written by GILAD TSFATY for SEE MAPPING - TM ###
# ### all rights reserved ###
#
# import numpy as np
# import csv
#
# # Constants for WGS-84
# a = 6378137.0  # Semi-major axis (meters)
# e2 = 6.69437999014e-3  # Square of eccentricity
#
#
# def geodetic_to_ecef(lat, lon, alt):
#     """
#     Converts geodetic coordinates to ECEF coordinates.
#     """
#     lat_rad = np.radians(lat)
#     lon_rad = np.radians(lon)
#     N = a / np.sqrt(1 - e2 * np.sin(lat_rad) ** 2)
#     X = (N + alt) * np.cos(lat_rad) * np.cos(lon_rad)
#     Y = (N + alt) * np.cos(lat_rad) * np.sin(lon_rad)
#     Z = (N * (1 - e2) + alt) * np.sin(lat_rad)
#     return np.array([X, Y, Z])
#
#
# def ecef_to_ned_matrix(lat_ref, lon_ref):
#     """
#     Computes the rotation matrix to convert ECEF to NED.
#     """
#     lat_rad = np.radians(lat_ref)
#     lon_rad = np.radians(lon_ref)
#     R = np.array([
#         [-np.sin(lat_rad) * np.cos(lon_rad), -np.sin(lat_rad) * np.sin(lon_rad), np.cos(lat_rad)],
#         [-np.sin(lon_rad), np.cos(lon_rad), 0],
#         [-np.cos(lat_rad) * np.cos(lon_rad), -np.cos(lat_rad) * np.sin(lon_rad), -np.sin(lat_rad)]
#     ])
#     return R
#
#
# def geodetic_to_ned(lat, lon, alt, lat_ref, lon_ref, alt_ref):
#     """
#     Converts geodetic coordinates to NED coordinates relative to a reference point.
#     """
#     ecef_ref = geodetic_to_ecef(lat_ref, lon_ref, alt_ref)
#     ecef_target = geodetic_to_ecef(lat, lon, alt)
#     delta_ecef = ecef_target - ecef_ref
#     R_ecef_to_ned = ecef_to_ned_matrix(lat_ref, lon_ref)
#     ned = R_ecef_to_ned @ delta_ecef
#     return ned
#
#
# def calculate_distance_and_azimuth(north, east, down):
#     """
#     Calculates the distance and azimuth from the origin (0,0,0) in NED coordinates.
#     """
#     # Distance
#     distance = np.sqrt(north**2 + east**2 + down**2)
#
#     # Azimuth (angle in degrees from North direction)
#     azimuth = np.degrees(np.arctan2(east, north))  # atan2 handles the quadrant
#     azimuth = azimuth % 360  # Ensure azimuth is in the range [0, 360)
#
#     return distance, azimuth
#
#
# def calculate_horizontal_distance_and_elevation(north, east, down):
#     """
#     Calculates horizontal distance and elevation angle.
#     """
#     horizontal_distance = np.sqrt(north**2 + east**2)
#     elevation_angle = np.degrees(np.arctan2(-down, horizontal_distance))  # -down because down is positive in NED
#     return horizontal_distance, elevation_angle
#
#
# if __name__ == "__main__":
#     print("Enter Reference Point (latitude, longitude, altitude):")
#     lat_ref = float(input("Latitude (degrees): "))
#     lon_ref = float(input("Longitude (degrees): "))
#     alt_ref = float(input("Altitude (meters): "))
#
#     target_points = []
#     target_names = []
#     ned_results = []
#
#     while True:
#         print("\nEnter Target Point:")
#         target_name = input("Target Name: ").strip()  # Prompt for target name
#         lat = float(input("Latitude (degrees): "))
#         lon = float(input("Longitude (degrees): "))
#         alt = float(input("Altitude (meters): "))
#
#         # Compute NED coordinates for this target point
#         ned_coords = geodetic_to_ned(lat, lon, alt, lat_ref, lon_ref, alt_ref)
#         ned_results.append(ned_coords)
#         target_points.append((lat, lon, alt))
#         target_names.append(target_name)
#
#         # Calculate distance, azimuth, horizontal distance, and elevation angle
#         distance, azimuth = calculate_distance_and_azimuth(ned_coords[0], ned_coords[1], ned_coords[2])
#         horizontal_distance, elevation_angle = calculate_horizontal_distance_and_elevation(
#             ned_coords[0], ned_coords[1], ned_coords[2]
#         )
#
#         # Print results
#         print(f"Target Name: {target_name}")
#         print(
#             f"NED Coordinates (meters): North: {ned_coords[0]:.3f}, East: {ned_coords[1]:.3f}, Down: {ned_coords[2]:.3f}")
#         print(f"Distance (m): {distance:.3f}, Azimuth (°): {azimuth:.3f}")
#         print(f"Horizontal Distance (m): {horizontal_distance:.3f}, Elevation Angle (°): {elevation_angle:.3f}")
#
#         # Ask if user wants to add another target point
#         add_more = input("\nDo you want to add another Target Point? (y/n): ").strip().lower()
#         if add_more != "y":
#             break
#
#     # Save data to a CSV file
#     file_name = "ned_coordinates.csv"
#     try:
#         with open(file_name, mode='w', newline='') as file:
#             writer = csv.writer(file)
#
#             # Write the header
#             writer.writerow(["Reference Point", f"Lat: {lat_ref:.3f}, Lon: {lon_ref:.3f}, Alt: {alt_ref:.3f}"])
#             writer.writerow([
#                 "Target #", "Target Name", "Latitude", "Longitude", "Altitude",
#                 "North (m)", "East (m)", "Down (m)",
#                 "Distance (m)", "Azimuth (°)",
#                 "Horizontal Distance (m)", "Elevation Angle (°)"
#             ])
#
#             # Write each Target Point and its data
#             for i, (name, target, ned) in enumerate(zip(target_names, target_points, ned_results), 1):
#                 distance, azimuth = calculate_distance_and_azimuth(ned[0], ned[1], ned[2])
#                 horizontal_distance, elevation_angle = calculate_horizontal_distance_and_elevation(
#                     ned[0], ned[1], ned[2]
#                 )
#                 writer.writerow([
#                     f"Target {i}", name,
#                     f"{target[0]:.8f}",  # Latitude with 8 decimal places
#                     f"{target[1]:.8f}",  # Longitude with 8 decimal places
#                     round(target[2], 3),  # Altitude
#                     round(ned[0], 3),  # North
#                     round(ned[1], 3),  # East
#                     round(ned[2], 3),  # Down
#                     round(distance, 3),  # Distance
#                     round(azimuth, 3),  # Azimuth
#                     round(horizontal_distance, 3),  # Horizontal Distance
#                     round(elevation_angle, 3)  # Elevation Angle
#                 ])
#
#         print(f"\nData has been successfully saved to '{file_name}'")
#     except Exception as e:
#         print(f"Error: Could not save data to '{file_name}'.\n{e}")













##### FILE INPUT #####


### written by GILAD TSFATY for SEE MAPPING - TM ###
### all rights reserved ###

import numpy as np
import csv

# Constants for WGS-84
a = 6378137.0  # Semi-major axis (meters)
e2 = 6.69437999014e-3  # Square of eccentricity

# Functions remain unchanged
def geodetic_to_ecef(lat, lon, alt):
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)
    N = a / np.sqrt(1 - e2 * np.sin(lat_rad) ** 2)
    X = (N + alt) * np.cos(lat_rad) * np.cos(lon_rad)
    Y = (N + alt) * np.cos(lat_rad) * np.sin(lon_rad)
    Z = (N * (1 - e2) + alt) * np.sin(lat_rad)
    return np.array([X, Y, Z])

def ecef_to_ned_matrix(lat_ref, lon_ref):
    lat_rad = np.radians(lat_ref)
    lon_rad = np.radians(lon_ref)
    R = np.array([
        [-np.sin(lat_rad) * np.cos(lon_rad), -np.sin(lat_rad) * np.sin(lon_rad), np.cos(lat_rad)],
        [-np.sin(lon_rad), np.cos(lon_rad), 0],
        [-np.cos(lat_rad) * np.cos(lon_rad), -np.cos(lat_rad) * np.sin(lon_rad), -np.sin(lat_rad)]
    ])
    return R

def geodetic_to_ned(lat, lon, alt, lat_ref, lon_ref, alt_ref):
    ecef_ref = geodetic_to_ecef(lat_ref, lon_ref, alt_ref)
    ecef_target = geodetic_to_ecef(lat, lon, alt)
    delta_ecef = ecef_target - ecef_ref
    R_ecef_to_ned = ecef_to_ned_matrix(lat_ref, lon_ref)
    ned = R_ecef_to_ned @ delta_ecef
    return ned

def calculate_distance_and_azimuth(north, east, down):
    distance = np.sqrt(north**2 + east**2 + down**2)
    azimuth = np.degrees(np.arctan2(east, north)) % 360
    return distance, azimuth

def calculate_horizontal_distance_and_elevation(north, east, down):
    horizontal_distance = np.sqrt(north**2 + east**2)
    elevation_angle = np.degrees(np.arctan2(-down, horizontal_distance))
    return horizontal_distance, elevation_angle

# Main program
if __name__ == "__main__":
    input_file = "input_data.csv"
    output_file = "ned_coordinates.csv"

    try:
        with open(input_file, mode='r') as file:
            reader = csv.reader(file)
            rows = list(reader)

            # First line is the reference point
            ref_lat, ref_lon, ref_alt = map(float, rows[0][1:])
            targets = rows[1:]

        ned_results = []

        with open(output_file, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Reference Point", f"Lat: {ref_lat:.8f}, Lon: {ref_lon:.8f}, Alt: {ref_alt:.3f}"])
            writer.writerow([
                "Target #", "Target Name", "Latitude", "Longitude", "Altitude",
                "North (m)", "East (m)", "Down (m)",
                "Distance (m)", "Azimuth (°)",
                "Horizontal Distance (m)", "Elevation Angle (°)"
            ])

            for i, (name, lat, lon, alt) in enumerate(targets, 1):
                lat, lon, alt = map(float, (lat, lon, alt))
                ned = geodetic_to_ned(lat, lon, alt, ref_lat, ref_lon, ref_alt)
                distance, azimuth = calculate_distance_and_azimuth(ned[0], ned[1], ned[2])
                horizontal_distance, elevation_angle = calculate_horizontal_distance_and_elevation(
                    ned[0], ned[1], ned[2]
                )
                writer.writerow([
                    f"Target {i}", name.strip(),
                    f"{lat:.8f}", f"{lon:.8f}", f"{alt:.3f}",
                    round(ned[0], 3), round(ned[1], 3), round(ned[2], 3),
                    round(distance, 3), round(azimuth, 3),
                    round(horizontal_distance, 3), round(elevation_angle, 3)
                ])
        print(f"Data successfully saved to '{output_file}'.")

    except Exception as e:
        print(f"Error: {e}")

