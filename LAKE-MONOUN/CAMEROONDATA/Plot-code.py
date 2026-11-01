import geopandas as gpd
import matplotlib.pyplot as plt

# -----------------------------
# 1. Load the shapefile
# -----------------------------
# Replace with the path to your Lake Monoun shapefile (.shp file)
#shapefile_path = "Lake_Monoun.shp"
shapefile_path = "gis_osm_landuse_a_free_1.shp"
lake = gpd.read_file(shapefile_path)

# -----------------------------
# 2. Inspect attributes
# -----------------------------
#print(lake.head())   # See the attribute table (columns like 'name', 'type', etc.)

# -----------------------------
# 3. Filter for Lake Monoun
# -----------------------------
#If the shapefile contains multiple lakes, filter by name
#Adjust the column name if different (e.g., 'NAME', 'name', 'lake_name')
#lake_monoun = lake[lake['name'].str.contains("Monoun", case=False)]

lake_monoun = lake[lake['name'].str.contains("Monoun", case=False, na=False)]

print(lake_monoun)
exit()
##############################
print(lake_monoun)
print(lake.columns)   # see available attribute fields
print(lake[['name']].head())  # inspect the first few entries
exit()
# -----------------------------
# 4. Plot the lake
# -----------------------------
fig, ax = plt.subplots(figsize=(8, 8))
lake_monoun.plot(ax=ax, color="skyblue", edgecolor="black")

# Add title and labels
ax.set_title("Lake Monoun (OpenStreetMap Data)", fontsize=14)
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")

plt.show()

