import os
import obspy
from obspy import UTCDateTime
from obspy.clients.nrl import NRL
from obspy.core.inventory import Inventory, Network, Station, Channel, Site, Response
import matplotlib.pyplot as plt
import glob
import numpy as np
import pytz
from datetime import datetime
import shutil






class GeophoneArray:
    def __init__(self, recorders=None, channels=None, source_offset=0, geophone_spacing=1):
        """
        Initialize the GeophoneArray object.

        Parameters:
        - recorders (list): A list of recorder names (default is ['GECK1', 'GECK2', 'GECK3', 'GECK4', 'GECK5', 'GECK6']).
        - channels (list): A list of channel identifiers (default is ['E', 'N', 'Z', 'O']).
        - source_offset (float): The initial offset of the source in meters (default is 0).
        - geophone_spacing (float): The spacing between geophones in meters (default is 1).
        """
        self.recorders = recorders if recorders else ['GECK1', 'GECK2', 'GECK3', 'GECK4', 'GECK5', 'GECK6']
        self.channels = channels if channels else ['E', 'N', 'Z', 'O']
        self.source_offset = source_offset
        self.geophone_spacing = geophone_spacing

    def __iter__(self):
        """
        Returns an iterator that yields recorder-channel pairs with their corresponding distance.

        Yields:
        - tuple: A tuple (recorder, channel, source_distance) representing the current recorder, channel, and distance from the source.
        """
        for i, recorder in enumerate(self.recorders):
            for j, channel in enumerate(self.channels):
                source_distance = self.calculate_source_distance(i, j)
                yield (recorder, channel, source_distance)

    def iter_reverse(self):
        """
        Returns an iterator that yields recorder-channel pairs with their corresponding distance in reverse order.

        Yields:
        - tuple: A tuple (recorder, channel, source_distance) representing the current recorder, channel, and distance from the source.
        """
        for i, recorder in reversed(list(enumerate(self.recorders))):
            for j, channel in reversed(list(enumerate(self.channels))):
                source_distance = self.calculate_source_distance(i, j)
                yield (recorder, channel, source_distance)

    def calculate_source_distance(self, recorder_index, channel_index):
        """
        Calculate the distance for a given recorder and channel index relative to the source.

        Parameters:
        - recorder_index (int): The index of the recorder in the list.
        - channel_index (int): The index of the channel in the list.

        Returns:
        - float: The calculated distance from the source.
        """
        # Total number of geophones in the array
        total_geophones = len(self.recorders) * len(self.channels)
        # Calculate the current geophone index, from the first to the last
        geophone_index = recorder_index * len(self.channels) + channel_index
        # Calculate the corresponding source distance by reversing the index order
        reversed_index = total_geophones - geophone_index - 1
        return self.source_offset + reversed_index * self.geophone_spacing

    def __len__(self):
        """
        Returns the total number of recorder-channel pairs.

        Returns:
        - int: The number of recorder-channel pairs.
        """
        return len(self.recorders) * len(self.channels)

    def __repr__(self):
        """
        Returns a string representation of the GeophoneArray object.

        Returns:
        - str: A string representation of the object.
        """
        return (f"GeophoneArray(recorders={self.recorders}, channels={self.channels}, "
                f"source_offset={self.source_offset}, geophone_spacing={self.geophone_spacing})")


def convert_aus_time_to_utc(time_input, timezone_str):
    """
    Convert a datetime from any mainland Australian time zone to UTC.
    
    Parameters:
    time_input (str or datetime): A datetime string in format "YYYY-MM-DDTHH:MM:SS" or a datetime.datetime object.
    timezone_str (str): The Australian time zone as a string, one of "AEST", "ACST", "AWST".
    
    Returns:
    str: The corresponding datetime string in UTC in ISO format.
    """
    # Define the time zones
    timezones = {
        'AEST': 'Australia/Sydney',  # Eastern Standard Time (UTC+10:00)
        'ACST': 'Australia/Adelaide',  # Central Standard Time (UTC+9:30)
        'AWST': 'Australia/Perth',  # Western Standard Time (UTC+8:00)
    }

    # Get the corresponding timezone
    if timezone_str not in timezones:
        raise ValueError(f"Unsupported timezone: {timezone_str}. Supported timezones are: {', '.join(timezones.keys())}")
    
    aus_timezone = pytz.timezone(timezones[timezone_str])

    # Convert time_input to a datetime object if it's a string
    if isinstance(time_input, str):
        local_time = datetime.strptime(time_input, "%Y-%m-%dT%H:%M:%S")
    elif isinstance(time_input, datetime):
        local_time = time_input
    else:
        raise TypeError("time_input must be a string or a datetime.datetime object")

    # Localize the datetime object to the given Australian timezone
    local_time = aus_timezone.localize(local_time)

    # Convert the local time to UTC
    utc_time = local_time.astimezone(pytz.utc)

    # Return the UTC time as a string in ISO format
    return utc_time.strftime("%Y-%m-%dT%H:%M:%SZ")


def change_fho_channel_id(trace):
    """
    Change the FHO channel id from 01 to 00.
    """
    if 'FHO' in trace.id:
        parts = trace.id.split('.')
        if parts[2] == '01':
            parts[2] = '00'
        trace.id = '.'.join(parts)


def order_traces_by_pattern(stream, pattern):
    """
    Order channels within each recorder based on the user input pattern.
    """
    order_dict = {'Z': 'FHZ', 'O': 'FHO', 'N': 'FHN', 'E': 'FHE'}
    pattern_list = pattern.split(',')

    # Group traces by geophone set and retain original order
    geophone_groups = {}
    ordered_geophones = []
    for trace in stream:
        change_fho_channel_id(trace)
        geophone_id = trace.id.split('.')[1]  # Assuming geophone id is the second element
        if geophone_id not in geophone_groups:
            geophone_groups[geophone_id] = []
            ordered_geophones.append(geophone_id)
        geophone_groups[geophone_id].append(trace)
    
    ordered_traces = []
    for geophone_id in ordered_geophones:
        traces = geophone_groups[geophone_id]
        for char in pattern_list:
            if char in order_dict:
                ordered_traces.extend([trace for trace in traces if order_dict[char] in trace.id])
    
    return obspy.Stream(ordered_traces)




def plot_custom_section(stream, start_time, end_time, save_path=None, pad_dist=2):
    """
    Plot the stack of seismograms with the waves filled on one side of the zero axis.

    Parameters:
    - stream (Stream): The ObsPy Stream containing the traces.
    - start_time (UTCDateTime): The start time for plotting.
    - end_time (UTCDateTime): The end time for plotting.
    - save_path (str): The path to save the plot (optional).
    - pad_dist (float): The padding distance for the x-axis limits.
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    distances = []
    total_duration = end_time - start_time  # Calculate total duration in seconds

    for trace in stream:
        times = trace.times("relative")  # Relative times for plotting
        data = trace.data
        distance = trace.stats.distance
        distances.append(distance)

        # Normalize data for better visualization
        data = data / max(abs(data))

        # Shift times to be relative to the start_time and convert to seconds
        relative_times = times + (trace.stats.starttime - start_time)

        # Create an array for distance with the same length as times and data
        distance_array = np.full_like(relative_times, distance)

        # Plot the trace with filled positive areas
        ax.fill_betweenx(relative_times, distance_array, distance_array + data, where=(data > 0), facecolor='blue', alpha=0.5)
        ax.plot(distance_array + data, relative_times, color="black", linewidth=0.5)

        # Add text annotation for the recorder and channel
        recorder_channel = f"{trace.stats.station}/{trace.stats.channel}"
        ax.text(distance, 0, recorder_channel, fontsize=8, ha='center', va='bottom', rotation=90)

    # Set the y-axis limit from 0 to total duration
    ax.set_ylim(0, total_duration)

    # Set x-axis limits with padding
    ax.set_xlim(np.min(distances) - pad_dist, np.max(distances) + pad_dist)
    ax.set_ylabel("Time (s)")
    ax.set_xlabel("Distance (m)")
    ax.set_title(f"Time: {start_time.strftime('%Y-%m-%d %H:%M:%S')} - {end_time.strftime('%Y-%m-%d %H:%M:%S')}")

    if save_path:
        plt.savefig(save_path, format="png")
    plt.show()







def retrieve_and_combine_streams(base_directory, recorders, start_time, end_time, file_extension=".ms",
                                 reverse_order=False, source_offset=None, sensor_spacing=None, freqmin=None, freqmax=None,
                                 pattern=None, buffer=0, trim=True, scale=1, file_duration=60, verbose=False):
    """
    Retrieve files with the specified extension for all recorders from a nested directory structure,
    and combine them into a single stream with one trace per channel per recorder.

    Parameters:
    - base_directory (str): The base directory to search for files.
    - recorders (list or str): List of recorder names (e.g., ["GECK1", "GECK2"]) or a single recorder name.
    - start_time (UTCDateTime): The start time for the trimming interval.
    - end_time (UTCDateTime): The end time for the trimming interval.
    - file_extension (str): The file extension to look for (default is ".ms").
    - reverse_order (bool): If True, process recorders in reverse order.
    - source_offset (int or None): The initial offset of the source in meters. If None, distances are not applied.
    - sensor_spacing (int or None): The spacing between sensors in meters. If None, distances are not applied.
    - freqmin (float or None): The minimum frequency for the bandpass filter in Hz. If None, no filtering is applied.
    - freqmax (float or None): The maximum frequency for the bandpass filter in Hz. If None, no filtering is applied.
    - pattern (str): The pattern to use for ordering the channels (e.g., "O,Z,N,E"). If None, no reordering is done.
    - buffer (int): The buffer time in seconds to add on either side of the start and end times (default is 0).
    - trim (bool): If True, trim the stream to the exact start and end times.
    - scale (float): Scaling factor to apply to the data.
    - file_duration (int): Duration of each file in seconds (default is 60).
    - verbose (bool): If True, print informational messages (default is False).

    Returns:
    - Stream: A combined ObsPy Stream with the relevant data.
    """
    # Ensure recorders is a list, even if a single string is provided
    if isinstance(recorders, str):
        recorders = [recorders]

    combined_stream = obspy.Stream()
    geophone_index = 0  # Initialize geophone index globally

    if reverse_order:
        recorders = recorders[::-1]

    # Add a buffer on either side of the start and end times
    buffered_start_time = start_time - buffer
    buffered_end_time = end_time + buffer

    for recorder in recorders:
        # Create temporary streams for each channel
        temp_streams = {
            'FHN': obspy.Stream(),
            'FHE': obspy.Stream(),
            'FHZ': obspy.Stream(),
            'FHO': obspy.Stream()
        }

        # Construct the search pattern
        search_pattern = os.path.join(base_directory, "**", f"*{recorder}*{file_extension}")
        # Retrieve all files that match the pattern
        files = glob.glob(search_pattern, recursive=True)
        files = sorted(files)

        for file_path in files:
            # Skip files with zero bytes
            if os.path.getsize(file_path) == 0:
                if verbose:
                    print(f"Skipping zero-byte file: {file_path}")
                continue

            # Extract the timestamp from the filename
            file_name = os.path.basename(file_path)
            file_time_str = file_name.split()[0] + file_name.split()[1]
            try:
                file_time = UTCDateTime(file_time_str)
            except Exception as e:
                if verbose:
                    print(f"Error parsing time from filename {file_name}: {e}")
                continue

            # Calculate the file's end time
            file_end_time = file_time + file_duration

            # Check if there is overlap between the file's time range and the buffered interval
            if not (file_end_time <= buffered_start_time or file_time >= buffered_end_time):
                # Read the MiniSEED file
                stream = obspy.read(file_path)
                
                # Add the traces to the corresponding temporary stream
                for trace in stream:
                    channel = trace.stats.channel
                    if channel in temp_streams:
                        temp_streams[channel] += trace

            # Stop processing if the current file's timestamp is after the buffered end_time
            if file_time >= buffered_end_time:
                break

        # Combine the temporary streams into the final stream, ensuring only one trace per channel
        for channel, temp_stream in temp_streams.items():
            if len(temp_stream) > 0:
                # Merge the traces in each temporary stream
                temp_stream.merge(method=1, fill_value='interpolate')
                # Apply the bandpass filter if both freqmin and freqmax are provided
                if freqmin is not None and freqmax is not None:
                    temp_stream.filter("bandpass", freqmin=freqmin, freqmax=freqmax)
                # Scale and add the processed trace to the combined stream
                for trace in temp_stream:
                    trace.data = scale * trace.data
                    combined_stream += trace

    # Trim the combined stream to the exact start and end times
    if trim:
        combined_stream.trim(starttime=UTCDateTime(start_time), endtime=UTCDateTime(end_time))

    # Order the combined stream based on the pattern, if provided
    if pattern:
        combined_stream = order_traces_by_pattern(combined_stream, pattern)

    # Apply distances based on source_offset and sensor_spacing, if both are provided
    if source_offset is not None and sensor_spacing is not None:
        geophone_index = 0
        for trace in combined_stream:
            trace.stats.distance = source_offset + geophone_index * sensor_spacing
            geophone_index += 1  # Increment the index for each channel

    return combined_stream

def group_detections(detections, min_shots=4, max_shots=6, max_gap=10):
    """
    Groups detections into lists representing groups of sources based on time intervals.

    Parameters:
    - detections (list): List of detection dictionaries with 'time' and 'similarity'.
    - min_shots (int): Minimum number of shots required to form a group.
    - max_shots (int): Maximum number of shots allowed in a group.
    - max_gap (float): Maximum time gap (in seconds) allowed between consecutive detections within a group.

    Returns:
    - List of lists, where each sublist contains a group of detections.
    """
    groups = []
    current_group = []
    
    for i, detection in enumerate(detections):
        if not current_group:
            current_group.append(detection)
            continue

        # Calculate time gap in seconds between consecutive detections
        time_gap = detection['time'] - current_group[-1]['time']

        if time_gap <= max_gap:
            current_group.append(detection)
        else:
            if min_shots <= len(current_group) <= max_shots:
                groups.append(current_group)
            current_group = [detection]

        if i == len(detections) - 1:
            if min_shots <= len(current_group) <= max_shots:
                groups.append(current_group)

    return groups


def filter_detections_by_threshold(detections, threshold):
    """
    Filters detections based on a similarity threshold.

    Parameters:
    - detections (list): A list of detections, each detection is a dictionary with keys 'time' and 'similarity'.
    - threshold (float): The minimum similarity score a detection must have to be included in the output.

    Returns:
    - filtered_detections (list): A list of detections with similarity greater than or equal to the threshold.
    """
    # Filter detections based on the similarity threshold
    filtered_detections = [detection for detection in detections if detection['similarity'] >= threshold]
    
    return filtered_detections

def get_max_correlation_detection(group):
    """
    Returns the detection with the maximum correlation (similarity) in a group.

    Parameters:
    - group (list): A list of detections, where each detection is a dictionary with a 'similarity' key.

    Returns:
    - dict: The detection with the maximum similarity in the group.
    """
    if not group:
        return None  # Return None if the group is empty
    
    # Find the detection with the maximum similarity
    max_detection = max(group, key=lambda detection: detection['similarity'])
    
    return max_detection

def stack_channel(group, ref_time, recorder, channel, base_directory, order=0, endbuffer=1, startbuffer=0, scale=1):
    """
    Stack the traces for a specific channel across a group of detections.

    Parameters:
    - group (list): List of detection dictionaries.
    - ref_time (UTCDateTime): Reference time to align traces.
    - recorder (str): Recorder name.
    - channel (str): Channel name.
    - base_directory (str): Base directory where data is stored.
    - order (int): Order of the phase-weighted stack (default is 0).
    - endbuffer (float): Buffer time after the detection time (default is 1 second).
    - startbuffer (float): Buffer time before the detection time (default is 0 second).
    - scale (float): Scaling factor to apply to the data (default is 1).

    Returns:
    - Stream: Stacked ObsPy Stream object.
    """
    source_stream = obspy.Stream()

    for si in range(len(group)):
        # Get the waveform
        start_time = group[si]['time'] - startbuffer
        end_time = group[si]['time'] + endbuffer
        
        shot_stream = retrieve_and_combine_streams(
            base_directory=base_directory, 
            recorders=recorder, 
            start_time=start_time, 
            end_time=end_time,
            reverse_order=True, 
            source_offset=None, 
            sensor_spacing=None, 
            freqmin=None, 
            freqmax=None, 
            pattern=channel, 
            buffer=60, 
            trim=True,
            scale=scale
        )
          
        for trace in shot_stream:
            # Adjust the start time of the trace by subtracting the time shift
            time_shift = trace.stats.starttime - ref_time
            trace.stats.starttime -= time_shift
            source_stream.append(trace)
            
    # Detrend and stack the stream
    source_stream.detrend()
    pw_stack = source_stream.stack(
        group_by='all',
        stack_type=('pw', order),
        npts_tol=0,
        time_tol=0
    )
    
    return pw_stack

def save_stream_to_tab_delimited(stream, filename):
    """
    Save an ObsPy Stream object to a tab-delimited text file.

    Parameters:
    - stream (Stream): The ObsPy Stream object containing the traces.
    - filename (str): The path and name of the output file.
    """
    # Determine the number of traces (channels) and the length of the data
    n_channels = len(stream)
    n_samples = len(stream[0].data)

    # Create a 2D array to hold the data
    data_array = np.zeros((n_samples, n_channels))

    # Populate the array with data from each trace
    for i, trace in enumerate(stream):
        data_array[:, i] = trace.data

    # Save the data to a tab-delimited file
    header = "\t".join([f"Channel {i+1}" for i in range(n_channels)])
    np.savetxt(filename, data_array, delimiter="\t", header=header, comments='')


def process_and_save_stacks(data_dir, output_dir, groups, geophone_array,
                            overwrite=False, save_mseed=True, 
                            save_tab=False, save_su=False, save_segy=False,
                            startbuffer=1.01, endbuffer=1.5, save_in_single_dir=False, detrend=True):
    """
    Processes groups of shots, stacks the traces for each channel, and saves the results into a specified output directory.

    Parameters:
    - data_dir (str): The base directory containing the data.
    - output_dir (str): The directory where the stacked results will be saved.
    - groups (list or list of lists): A list of grouped shots to be processed, or a single group.
    - geophone_array (object): The geophone array object to iterate over.
    - overwrite (bool): If True, overwrite existing directories. If False, skip processing if the directory exists.
    - save_mseed (bool): If True, save the stacked stream as MiniSEED.
    - save_tab (bool): If True, save the stacked stream as a tab-delimited text file.
    - save_su (bool): If True, save the stacked stream as Seismic Unix (SU) format.
    - save_segy (bool): If True, save the stacked stream as SEGY format.
    - startbuffer (float): The time before the detection to start the extraction.
    - endbuffer (float): The time after the detection to end the extraction.
    - save_in_single_dir (bool): If True, save all stacked streams into a single directory.
    """
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Handle the case where a single group is passed instead of a list of groups
    if not isinstance(groups, list) or isinstance(groups[0], dict):
        groups = [groups]

    # Loop through each group and process
    for i, group in enumerate(groups):
        if save_in_single_dir:
            # If saving in a single directory, file name follows "s1.mseed", "s2.mseed", etc.
            sub_dir = output_dir
            base_filename = f"s{i+1}"
        else:
            # Define the subdirectory name based on the group index
            sub_dir = os.path.join(output_dir, f"s{i+1}")
            base_filename = "stacked_stream"
        
        # Check if the directory exists and handle based on the overwrite flag
        if not save_in_single_dir:
            if os.path.exists(sub_dir):
                if overwrite:
                    shutil.rmtree(sub_dir)  # Remove existing directory and its contents
                    os.makedirs(sub_dir)  # Create a new, empty directory
                else:
                    if os.listdir(sub_dir):
                        print(f"Skipping {sub_dir}: Directory already exists and contains files.")
                        continue
                    else:
                        print(f"Directory {sub_dir} exists but is empty. Proceeding with processing.")
            else:
                os.makedirs(sub_dir)  # Create the directory if it doesn't exist

        # Stacking and saving logic
        print(f"Processing group {i+1} and saving to {sub_dir}...")
        
        ref_time = get_max_correlation_detection(group)["time"]
        stacked_stream = obspy.Stream()
        geophone_index = 0
        for rec, chan, source_distance in geophone_array.iter_reverse():
            stacked = stack_channel(group, ref_time, recorder=rec, channel=chan, base_directory=data_dir, 
                                    order=1, endbuffer=endbuffer, startbuffer=startbuffer)
            trace = stacked[0]
            trace.stats.distance = source_distance
            if detrend is True:
                trace.detrend(type='linear')
            stacked_stream.append(trace)
        
        # Save the stacked stream based on the flags
        if save_mseed:
            mseed_file = os.path.join(sub_dir, f"{base_filename}.mseed")
            stacked_stream.write(mseed_file, format="MSEED")
            print(f"Saved stacked stream to {mseed_file}")

        if save_tab:
            tab_file = os.path.join(sub_dir, f"{base_filename}.txt")
            save_stream_to_tab_delimited(stacked_stream, tab_file)
            print(f"Saved stacked stream to {tab_file}")

        if save_su:
            # Convert data to float32 for SU and SEGY formats
            for trace in stacked_stream:
                trace.data = trace.data.astype(np.float32)
            
            su_file = os.path.join(sub_dir, f"{base_filename}.su")
            stacked_stream.write(su_file, format="SU")
            print(f"Saved stacked stream to {su_file}")

        if save_segy:
            segy_file = os.path.join(sub_dir, f"{base_filename}.sgy")
            stacked_stream.write(segy_file, format="SEGY", data_encoding=5)  # 5 is for IEEE float
            print(f"Saved stacked stream to {segy_file}")

    print("Processing complete.")


def save_individual_shots(data_dir, output_dir, groups, geophone_array,
                          target_number, startbuffer=1.01, endbuffer=1.5,
                          overwrite=False, save_mseed=True, detrend=True):
    """
    Processes groups of shots and saves each shot individually, ensuring each group has the target number of shots.

    Parameters:
    - data_dir (str): The base directory containing the data.
    - output_dir (str): The directory where the individual shots will be saved.
    - groups (list or list of lists): A list of grouped shots to be processed, or a single group.
    - geophone_array (object): The geophone array object to iterate over.
    - target_number (int): The target number of shots per group. Groups with fewer shots will have duplicates.
    - startbuffer (float): The time before the detection to start the extraction.
    - endbuffer (float): The time after the detection to end the extraction.
    - overwrite (bool): If True, overwrite existing files. If False, skip if files already exist.
    - save_mseed (bool): If True, save each shot as MiniSEED.
    """
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Handle the case where a single group is passed instead of a list of groups
    if not isinstance(groups, list) or isinstance(groups[0], dict):
        groups = [groups]

    file_index = 1  # For naming the output files like S1.mseed, S2.mseed, etc.

    # Loop through each group
    for group_index, group in enumerate(groups):
        # Sort the group by similarity to get the highest values first
        group_sorted = sorted(group, key=lambda x: x['similarity'], reverse=True)

        # If the group has more shots than the target_number, select only the top `target_number` shots
        if len(group_sorted) > target_number:
            selected_shots = group_sorted[:target_number]
        else:
            # If fewer shots, select the top shot and duplicate it as needed
            selected_shots = group_sorted[:]
            top_shot = selected_shots[0]  # Shot with the highest similarity
            while len(selected_shots) < target_number:
                selected_shots.append(top_shot)

        # Process each shot and save it
        for shot in selected_shots:
            ref_time = shot["time"]
            stacked_stream = obspy.Stream()
            for rec, chan, source_distance in geophone_array.iter_reverse():
                stacked = stack_channel([shot], ref_time, recorder=rec, channel=chan, base_directory=data_dir, 
                                        order=1, endbuffer=endbuffer, startbuffer=startbuffer)
                trace = stacked[0]
                trace.stats.distance = source_distance
                if detrend is True:
                    trace.detrend(type='linear')
                stacked_stream.append(trace)

            # Save the shot as MiniSEED
            if save_mseed:
                mseed_file = os.path.join(output_dir, f"S{file_index}.mseed")
                if overwrite or not os.path.exists(mseed_file):
                    stacked_stream.write(mseed_file, format="MSEED")
                    print(f"Saved shot {file_index} to {mseed_file}")
                file_index += 1

    print("Individual shot saving complete.")


def plot_average_correlation(groups, plot_title="Average Correlation for Each Group"):
    """
    Plots the average correlation value for each group of detections.

    Parameters:
    - groups (list): A list of groups, where each group is a list of detections.
    - plot_title (str): The title of the plot.
    """
    # Calculate the average correlation (similarity) for each group
    average_correlations = [np.mean([detection['similarity'] for detection in group]) for group in groups]
    
    # Plot the average correlations
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, len(average_correlations) + 1), average_correlations, marker='o', linestyle='-', color='b')
    plt.title(plot_title)
    plt.xlabel("Group Number")
    plt.ylabel("Average Correlation (Similarity)")
    plt.grid(True)
    plt.show()


def default_masw_channel_response(nrl_path):
    """
    Fetch the instrument response using a local NRL path.
    
    Parameters:
    nrl_path (str): Path to the local NRL directory.
    
    Returns:
    Response object: The response from the sensor and datalogger keys.
    """
    # Initialize the NRL with the given root path
    nrl = NRL(root=nrl_path)
    
    # Print the datalogger information
    datalogger_info = nrl.dataloggers['SeismologyResearchCentre']['Gecko']['128']['all']['1000 Hz']
    print(datalogger_info)
    
    # Print the sensor information
    sensor_info = nrl.sensors['HGSProducts']['HG-4']['8 Hz']['375']['None']['28.8']
    print(sensor_info)
    
    # Fetch the response object
    response = nrl.get_response(sensor_keys=['HGSProducts', 'HG-4', '8 Hz', '375', 'None', '28.8'],
                                datalogger_keys=['SeismologyResearchCentre', 'Gecko', '128', 'all', '1000 Hz'])
    
    return response

def build_inventory_from_trace_and_response(trace, channel_response):
    """
    Build an ObsPy Inventory from a trace and a given channel response.

    Parameters:
    - trace: ObsPy Trace object
    - channel_response: Channel Response object (as obtained from NRL or similar)
    
    Returns:
    - inventory: ObsPy Inventory object containing a single station and channel with the given response
    """
    # Create the Inventory object (empty at first)
    inv = Inventory(networks=[], source="Generated for Trace Response")

    # Create the Network object, taking network code from the trace
    net = Network(
        code=trace.stats.network or "",  # Use trace's network code or empty string
        stations=[]
    )

    # Create the Station object, using metadata from the trace
    sta = Station(
        code=trace.stats.station,  # Use trace's station code
        latitude=0.0,  # You can fill these values as needed or leave at default
        longitude=0.0,
        elevation=0.0,
        site=Site(name=trace.stats.station)
    )

    # Create the Channel object, using metadata from the trace
    cha = Channel(
        code=trace.stats.channel,  # Use trace's channel code
        location_code=trace.stats.location or "",  # Use location code from trace or empty string
        latitude=0.0,  # Again, you can fill this in with real values if you have them
        longitude=0.0,
        elevation=0.0,
        depth=0.0,
        azimuth=0.0,
        dip=0.0,
        sample_rate=trace.stats.sampling_rate  # Use sampling rate from trace
    )

    # Attach the provided channel response to the channel
    cha.response = channel_response

    # Add the channel to the station
    sta.channels.append(cha)

    # Add the station to the network
    net.stations.append(sta)

    # Add the network to the inventory
    inv.networks.append(net)

    return inv