import os
from PIL import Image
import pillow_heif
import piexif
from natsort import natsorted
from datetime import datetime
from pyproj import Proj, transform
import matplotlib.pyplot as plt


# Register HEIF format for Pillow
pillow_heif.register_heif_opener()

# Define the Photo data structure
class Photo:
    def __init__(self, shot_number, time, lat, lon, x=None, y=None):
        self.shot_number = shot_number
        self.time = time
        self.lat = lat
        self.lon = lon
        self.x = x
        self.y = y

    def minutes_elapsed(self, start_time):
        return (self.time - start_time).total_seconds() / 60.0

    def __repr__(self):
        return f"{self.shot_number},{self.time.isoformat()},{self.lat},{self.lon},{self.x},{self.y}"

# Extract EXIF data
def extract_exif_data(image):
    exif_data = image.info.get('exif')
    if exif_data:
        exif_dict = piexif.load(exif_data)
        return exif_dict
    return None

# Extract time and coordinates
def extract_time_and_coordinates(image):
    exif_data = extract_exif_data(image)
    if not exif_data:
        return None, None
    
    time_taken = None
    try:
        time_taken = exif_data['Exif'][36867].decode('utf-8')  # DateTimeOriginal
        time_taken = datetime.strptime(time_taken, "%Y:%m:%d %H:%M:%S")
    except KeyError:
        pass  # Handle cases where DateTimeOriginal is not present

    gps_info = exif_data.get('GPS', None)
    
    if gps_info:
        def convert_to_degrees(value):
            d = value[0][0] / value[0][1]
            m = value[1][0] / value[1][1]
            s = value[2][0] / value[2][1]
            return d + (m / 60.0) + (s / 3600.0)
        
        lat = convert_to_degrees(gps_info[2])  # GPSLatitude
        lon = convert_to_degrees(gps_info[4])  # GPSLongitude
        
        lat_ref = gps_info[1].decode('utf-8')  # GPSLatitudeRef
        lon_ref = gps_info[3].decode('utf-8')  # GPSLongitudeRef

        if lat_ref != "N":
            lat = -lat
        if lon_ref != "E":
            lon = -lon
    else:
        lat, lon = None, None

    return time_taken, (lat, lon)

# Convert GPS coordinates to local Cartesian coordinates
def convert_to_cartesian(photos):
    if not photos:
        return

    # Use the first photo's coordinates as the origin
    origin_lat = photos[0].lat
    origin_lon = photos[0].lon

    # Define a projection (use a local UTM zone or a suitable projection for your area)
    proj_latlon = Proj(proj='latlong', datum='WGS84')
    proj_xy = Proj(proj='utm', zone=33, datum='WGS84')  # Adjust zone as needed

    for photo in photos:
        x, y = transform(proj_latlon, proj_xy, photo.lon, photo.lat)
        origin_x, origin_y = transform(proj_latlon, proj_xy, origin_lon, origin_lat)

        # Calculate the local Cartesian coordinates relative to the origin
        photo.x = x - origin_x
        photo.y = y - origin_y

# Process images and create a list of Photo objects
def process_images(directory):
    photo_list = []
    heic_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.lower().endswith('.heic')]
    heic_files = natsorted(heic_files)
    
    for idx, image_path in enumerate(heic_files):
        try:
            image = Image.open(image_path)
            time_taken, coordinates = extract_time_and_coordinates(image)
            lat, lon = coordinates if coordinates else (None, None)
            photo = Photo(shot_number=idx+1, time=time_taken, lat=lat, lon=lon)
            photo_list.append(photo)
        except Exception as e:
            print(f"Error processing {image_path}: {e}")

    # Convert to local Cartesian coordinates
    convert_to_cartesian(photo_list)

    return photo_list

# Write the list of photos to a plain text file
def write_photos_to_file(photo_list, filename):
    with open(filename, 'w') as file:
        for photo in photo_list:
            file.write(f"{photo}\n")

# Read the list of photos from a plain text file
def read_photos_from_file(filename):
    photo_list = []
    with open(filename, 'r') as file:
        for line in file:
            shot_number, time_str, lat, lon, x, y = line.strip().split(',')
            time = datetime.fromisoformat(time_str)
            lat, lon = float(lat), float(lon)
            x, y = float(x), float(y)
            photo = Photo(shot_number=int(shot_number), time=time, lat=lat, lon=lon, x=x, y=y)
            photo_list.append(photo)
    return photo_list

# Plot the photos on a Cartesian plane
def plot_photos(photo_list):
    if not photo_list:
        print("No photos to plot.")
        return

    start_time = photo_list[0].time

    plt.figure(figsize=(10, 6))
    for photo in photo_list:
        plt.scatter(photo.x, photo.y, label=f"{photo.minutes_elapsed(start_time):.1f} min")
        plt.text(photo.x, photo.y, f"{photo.minutes_elapsed(start_time):.1f} min", fontsize=9, ha='right')

    plt.title("Photo Locations in Local Cartesian Coordinates")
    plt.xlabel("X (meters)")
    plt.ylabel("Y (meters)")
    plt.grid(True)
    plt.legend()
    plt.show()