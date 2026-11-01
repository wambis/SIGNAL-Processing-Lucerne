import geopandas as gpd
import matplotlib.pyplot as plt
import os
import rasterio
from rasterio.plot import show
from matplotlib.colors import LightSource
import numpy as np





#Path to the folder containing all your shapefiles
folder_path = "CAMEROONDATA"

lake_monoun = None

#Loop through all files in the folder
for file in os.listdir(folder_path):
    if file.endswith(".shp"):  # only load shapefiles
        shp_path = os.path.join(folder_path, file)
        gdf = gpd.read_file(shp_path)

        #Check available columns
        print(f"Inspecting {file}: {gdf.columns}")

        #Try filtering by name (adjust column name if needed)
        if "name" in gdf.columns:
            match = gdf[gdf["name"].str.contains("Monoun", case=False, na=False)]
            if not match.empty:
                lake_monoun = match
                print(f"Found Lake Monoun in {file}")
                break

#print(lake_monoun)
#get the tography file
# Example: SRTM DEM GeoTIFF file covering Cameroon
dem = rasterio.open("Lake_Monoun_Bathymetry.tif")
#exit()
# Read DEM data
#dem_data = dem.read(1)
dem_data = dem.read(1)
dem_data = dem_data.astype(float)
dem_data[dem_data == dem.nodata] = np.nan
# Plot if found
if lake_monoun is not None:
    fig, ax = plt.subplots(figsize=(8, 8),facecolor="white")
    #lake_monoun.plot(ax=ax, color="skyblue", edgecolor="black")
    lake_monoun.plot(ax=ax, color="skyblue", edgecolor="black", linewidth=3)
    #3.Reproject lake to DEM CRS
    #-----------------------------
    lake = lake_monoun.to_crs(dem.crs)

    #Overlay lake boundary
#    lake.plot(ax=ax, facecolor="none", edgecolor="blue", linewidth=2)
    #Plot DEM with terrain colormap
#    show(dem, ax=ax, cmap="terrain") # <-- change colormap here

#    #Create hillshade (shaded relief)
#    ax.imshow(dem_data, cmap="terrain", extent=(dem.bounds.left, dem.bounds.right, dem.bounds.bottom, dem.bounds.top))
#    ax.set_facecolor("white")
    #ax.imshow(dem_data, cmap="grey", extent=(dem.bounds.left, dem.bounds.right, dem.bounds.bottom, dem.bounds.top))
    #Extend the map
    #Get bounding box of the lake
    bounds = lake_monoun.total_bounds  # [minx, miny, maxx, maxy]
    minx, miny, maxx, maxy = bounds

    # Extend longitude/latitude by a margin (e.g., 0.05 degrees)
    margin = 0.005
    ax.set_xlim(minx - margin, maxx + margin)  # longitude (x-axis)
    ax.set_ylim(miny - margin, maxy + margin)  # latitude (y-axis)
    #####################################
    #Control tick spacing (e.g., every 0.02 degrees)
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.003))
    #ax.yaxis.set_major_locator(plt.MultipleLocator(0.02))


    ax.set_title("Lake Monoun", fontsize=14)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    figname   = "Map_LakeMonoun_TEST.pdf"
    #Set the font size of yticks
    plt.yticks(fontsize=13)
    # Set the font size of xticks
    plt.xticks(fontsize=14)
    ##Space between the subplots
    #Aligned all the ylabels of the figure
    #fig.align_ylabels()
    #Save the figure
    plt.savefig(figname, bbox_inches = 'tight', dpi = 300)
    #plt.show()
else:
    print("Lake Monoun not found in any shapefile.")

