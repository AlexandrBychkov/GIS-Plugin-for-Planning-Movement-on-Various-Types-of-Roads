from qgis.core import (
    QgsProject,
    QgsVectorLayer,
    QgsFeature,
    QgsGeometry,
    QgsField,
    QgsFields,
    QgsPointXY
)
from PyQt5.QtCore import QVariant
from collections import defaultdict

# Name of the source layer
input_layer_name = "Normalize_Data"
output_layer_name = "Nodal_points"

# Getting the input layer
input_layer = QgsProject.instance().mapLayersByName(input_layer_name)
if not input_layer:
    raise Exception(f"Layer {input_layer_name} not faund!")
input_layer = input_layer[0]

# Creating an output layer
fields = QgsFields()
fields.append(QgsField("fid1", QVariant.Int))
fields.append(QgsField("osm_id1", QVariant.Int))
fields.append(QgsField("vertex_index1", QVariant.Int))
fields.append(QgsField("fid2", QVariant.Int))
fields.append(QgsField("osm_id2", QVariant.Int))
fields.append(QgsField("vertex_index2", QVariant.Int))
fields.append(QgsField("X", QVariant.Double))
fields.append(QgsField("Y", QVariant.Double))

output_layer = QgsVectorLayer(
    "Point?crs=EPSG:3857", output_layer_name, "memory"
)
output_layer_data_provider = output_layer.dataProvider()
output_layer_data_provider.addAttributes(fields)
output_layer.updateFields()

# Dictionary for storing points
points_dict = {}

# Filling in the dictionary of points
for feature in input_layer.getFeatures():
    x = feature["X"]
    y = feature["Y"]
    point_key = (x, y)

    if point_key not in points_dict:
        points_dict[point_key] = []

    points_dict[point_key].append({
        "fid_point": feature["fid_point"],
        "osm_id": feature["osm_id"],
        "vertex_index": feature["vertex_index"]
    })

# Creating features for intersections
features_to_add = []
for point_key, point_features in points_dict.items():
    if len(point_features) > 1:
        for i in range(len(point_features)):
            for j in range(i + 1, len(point_features)):
                point1 = point_features[i]
                point2 = point_features[j]

                new_feature = QgsFeature(output_layer.fields())
                new_feature.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(*point_key)))
                new_feature.setAttributes([
                    point1["fid_point"], point1["osm_id"], point1["vertex_index"],
                    point2["fid_point"], point2["osm_id"], point2["vertex_index"],
                    point_key[0], point_key[1]
                ])
                features_to_add.append(new_feature)

output_layer_data_provider.addFeatures(features_to_add)
output_layer.updateExtents()

# Adding a layer to a project
QgsProject.instance().addMapLayer(output_layer)

print(f"Layer '{output_layer_name}' successfully created with {len(features_to_add)} intersections.")

# --- ADDING EXTREME POINTS OF ROADS ---
extreme_features_to_add = []

# Finding the minimum and maximum vertex_index for each osm_id
osm_groups = {}
for feature in input_layer.getFeatures():
    osm_id = feature["osm_id"]
    if osm_id not in osm_groups:
        osm_groups[osm_id] = []
    osm_groups[osm_id].append(feature)

for osm_id, features in osm_groups.items():
    # Sorting by vertex_index
    features.sort(key=lambda f: f["vertex_index"])
    min_feature = features[0]
    max_feature = features[-1]

    for extreme_feature in [min_feature, max_feature]:
        x = extreme_feature["X"]
        y = extreme_feature["Y"]
        point_key = (x, y)

        # Checking if a point exists in the source layer (to avoid duplication)
        exists = any(
            f["X"] == x and f["Y"] == y
            for f in output_layer.getFeatures()
        )
        
        if not exists:
            new_feature = QgsFeature(output_layer.fields())
            new_feature.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
            new_feature.setAttributes([
                extreme_feature["fid_point"],  # fid1
                extreme_feature["osm_id"],    # osm_id1
                extreme_feature["vertex_index"],  # vertex_index1
                None,  # fid2
                None,  # osm_id2
                None,  # vertex_index2
                x, y
            ])
            extreme_features_to_add.append(new_feature)

# Adding edge points to a layer
output_layer_data_provider.addFeatures(extreme_features_to_add)
output_layer.updateExtents()

print(f"Added {len(extreme_features_to_add)} extreme points of roads.")

#Point aggregation

# We get a layer with intersections
layer = QgsProject.instance().mapLayersByName('Nodal_points')[0]

# Dictionary for storing points by coordinates
points_dict = defaultdict(list)

# Overlay Counter
intersection_counter = 1

# We go through all the points of the layer and group them by coordinates
for feature in layer.getFeatures():
    x = feature['X']
    y = feature['Y']
    points_dict[(x, y)].append(feature)

# Switching the layer to edit mode
layer.startEditing()

# A list of existing fields
existing_fields = [field.name() for field in layer.fields()]

# We only process groups with overlays.
for (x, y), features in points_dict.items():
    if len(features) > 1:  # The overlay was found
        
        # A list for storing unique attribute combinations
        unique_combinations = []

        # Collecting all the unique parameters
        for feature in features:
            fid1 = feature['fid1']
            osm_id1 = feature['osm_id1']
            vertex_index1 = feature['vertex_index1']
            fid2 = feature['fid2']
            osm_id2 = feature['osm_id2']
            vertex_index2 = feature['vertex_index2']

            # Creating combinations of parameters
            combination1 = (fid1, osm_id1, vertex_index1)
            combination2 = (fid2, osm_id2, vertex_index2)

            # Adding unique combinations to the list
            if combination1 not in unique_combinations:
                unique_combinations.append(combination1)
            if combination2 not in unique_combinations:
                unique_combinations.append(combination2)

        # We check and add the missing columns
        for index in range(1, len(unique_combinations) + 1):
            if f'fid{index}' not in existing_fields:
                layer.addAttribute(QgsField(f'fid{index}', QVariant.Int))
                existing_fields.append(f'fid{index}')
            if f'osm_id{index}' not in existing_fields:
                layer.addAttribute(QgsField(f'osm_id{index}', QVariant.String))
                existing_fields.append(f'osm_id{index}')
            if f'vertex_index{index}' not in existing_fields:
                layer.addAttribute(QgsField(f'vertex_index{index}', QVariant.Int))
                existing_fields.append(f'vertex_index{index}')

        # Overwriting the parameters for all points in the group
        for feature in features:
            for index, (fid, osm_id, vertex_index) in enumerate(unique_combinations, start=1):
                layer.changeAttributeValue(feature.id(), layer.fields().indexOf(f'fid{index}'), fid)
                layer.changeAttributeValue(feature.id(), layer.fields().indexOf(f'osm_id{index}'), osm_id)
                layer.changeAttributeValue(feature.id(), layer.fields().indexOf(f'vertex_index{index}'), vertex_index)

        
        intersection_counter += 1

# Removing duplicate points
print("Removing duplicate points...")
unique_coords = set()  # A set for storing unique coordinates
features_to_delete = []  # A list for storing the IDs of points to delete

for feature in layer.getFeatures():
    coords = (feature['X'], feature['Y'])  # Coordinates of the point
    if coords in unique_coords:
        # If the coordinates are already there, add the point to the list for deletion.
        features_to_delete.append(feature.id())
    else:
        # Otherwise, we add the coordinates to the set
        unique_coords.add(coords)

# Removing duplicate points
layer.deleteFeatures(features_to_delete)

print(f"Deleted {len(features_to_delete)} duplicate points.")

# Saving the changes
layer.updateFields()
layer.commitChanges()

print("The layer has been successfully updated!")

