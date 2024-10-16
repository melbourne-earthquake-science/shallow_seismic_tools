# shallow_seismic_tools

These Python functions, built on **ObsPy**, assist in preprocessing data from a geophone array at the University of Melbourne, School of Geography, Earth and Atmospheric Sciences. The array is mainly used in **MASW-type surveys**.

The geophone kit includes 24 single-component, 8 Hz geophones, connected to 6 **Gecko** digitizers. Instrument response information for both components is available through the **IRIS NRL (Nominal Response Library)**. Based on existing settings, you can construct the response as follows:

## Instrument Response

You can fetch the response using the following code:

```python
# Fetch the response object
response = nrl.get_response(sensor_keys=['HGSProducts', 'HG-4', '8 Hz', '375', 'None', '28.8'],
                            datalogger_keys=['SeismologyResearchCentre', 'Gecko', '128', 'all', '1000 Hz'])
```

## Picking Shots and Stacking

One of the main objectives of these functions is to (semi-)automate the **windowing** and **stacking** of seismic shots in a roll-along geophone array. The functionality allows you to select a **template shot** (ideally a clear one) and uses correlation to identify other shots within the array. The shots can then be grouped and stacked accordingly.

Look for an example in the notebooks directory.

## Data Formats and Inversion

The University of Melbourne's Earthquake Hazard Group uses **KGS Surfseis6**, a proprietary package that requires input data in **SEG-2** format. While ObsPy can read SEG-2 files, it cannot write them. For writing SEG-2 files, the open-source waveform viewer **Waves** can be used, but it requires manual interaction.

Alternatively, the **Geopsy** Python package can both read and write a variety of formats, including SEG-2. We currently use Geopsy via a Docker image (available at [Geopsy Docker](https://github.com/jpvantassel/geopsy-docker)). This works on Linux, but not on Mac/ARM.

### Converting MiniSEED Files to SEG-2

To convert MiniSEED files to SEG-2 using Geopsy, follow these steps:

1. Ensure that a local directory `/path/data` exists, containing a set of MiniSEED files.
2. Copy the example script `misc/batch_convert.sh` into `/path/`.

The `MINISEED_DIR` in the shell script must be relative to the mount point. By default, the local directory is mounted at `/mnt/path`, so you should set `MINISEED_DIR="/mnt/data/"` in the script.

3. Pull the Docker image and run it as follows:

```bash
# Pull the Geopsy Docker image
docker pull jpvantassel/geopsy-docker -t geopsyimage

# Run the Docker container
docker run -it --rm \
    --name geopsy-container \
    -v /path:/mnt/path \
    geopsyimage

# Inside the container, convert the files
cd /mnt/path
./batch_convert.sh
```

This will convert the MiniSEED files in `/path/data` into SEG-2 format with the `.dat` extension.

## Additional Notes

- **Linux Compatibility**: The Geopsy Docker image works on Linux but does not currently support Mac/ARM systems.
- **Waves** can be used as a backup option for manual conversion if Docker-based automation isn't suitable for your workflow.