from qgis.core import (
    QgsProject, QgsRectangle, QgsCoordinateReferenceSystem,
    QgsCoordinateTransform, QgsFeature, QgsGeometry,
    QgsPointXY, QgsField, QgsVectorLayer
)
from PyQt5.QtCore import QVariant
from qgis.PyQt.QtWidgets import QInputDialog, QMessageBox
from qgis.utils import iface

# Default coordinates (EPSG:3857)
default_coords = {
    "north": 7306903.5433,
    "west": 4220660.2002,
    "south": 7294319.0911,
    "east": 4239298.3101
}

def get_coordinate(prompt, default_value):
    """A function for entering coordinates with a check for cancellation."""
    coord, ok = QInputDialog.getDouble(None, "Entering coordinates", prompt, value=default_value, decimals=6)
    if not ok:
        raise ValueError("The input was canceled.")
    return coord

def add_points_to_layer(points, layer_name="Normalize_Data"):
    """Adding points to a layer, even if the coordinates match the existing ones."""
    # Getting or creating a point layer
    existing_layer = QgsProject.instance().mapLayersByName(layer_name)
    if existing_layer:
        point_layer = existing_layer[0]
    else:
        # Creating a new layer with the necessary attributes
        point_layer = QgsVectorLayer("Point?crs=EPSG:3857", layer_name, "memory")
        provider = point_layer.dataProvider()
        provider.addAttributes([
            QgsField("fid", QVariant.Int),
            QgsField("osm_id", QVariant.Int),
            QgsField("X", QVariant.Double),
            QgsField("Y", QVariant.Double),
            QgsField("vertex_index", QVariant.Int)
        ])
        point_layer.updateFields()
        QgsProject.instance().addMapLayer(point_layer)

    # Converting coordinates from EPSG:4326 to EPSG:3857
    crs_project = QgsCoordinateReferenceSystem("EPSG:3857")
    transform = QgsCoordinateTransform(QgsCoordinateReferenceSystem("EPSG:4326"), crs_project, QgsProject.instance())

    provider = point_layer.dataProvider()

    features = []
    for point in points:
        feature = QgsFeature()
        transformed_point = transform.transform(point["point"])
        x_coord, y_coord = transformed_point.x(), transformed_point.y()

        # Adding a new point with attributes
        feature.setGeometry(QgsGeometry.fromPointXY(transformed_point))
        feature.setAttributes([
            point["fid"],          # fid
            point["osm_id"],       # osm_id
            x_coord,               # X
            y_coord,               # Y
            point["vertex_index"]  # vertex_index
        ])
        features.append(feature)

    provider.addFeatures(features)
    point_layer.updateExtents()

def extract_road_vertices_in_area(north, west, south, east):
    """Extracting the vertexes of roads in the specified area and adding them to the layer."""
    bounding_box = QgsRectangle(west, south, east, north)
    crs_project = QgsCoordinateReferenceSystem("EPSG:3857")
    layer_name = "road_layer"
    layer = QgsProject.instance().mapLayersByName(layer_name)

    if not layer:
        QMessageBox.critical(None, "Error", f"Layer '{layer_name}' not found.")
        return

    layer = layer[0]
    layer_crs = layer.crs()
    if layer_crs != crs_project:
        try:
            transform = QgsCoordinateTransform(crs_project, layer_crs, QgsProject.instance())
            bounding_box = transform.transform(bounding_box)
        except Exception as e:
            QMessageBox.critical(None, "Transformation error", f"Couldn't perform coordinate conversion: {e}")
            return

    points = []
    for feature in layer.getFeatures():
        geometry = feature.geometry()
        if geometry:
            for index, vertex in enumerate(geometry.vertices()):  # Добавляем индекс вершины
                point = QgsPointXY(vertex)
                if bounding_box.contains(point):
                    points.append({
                        "fid": feature.id(),           # fid from the original layer
                        "osm_id": feature["osm_id"],   # osm_id from the original layer
                        "point": point,                # The point itself
                        "vertex_index": index          # vertex index
                    })

    if points:
        add_points_to_layer(points)
    else:
        QMessageBox.information(None, "Result", "No road vertexes were found in the specified area.")

try:
    # Dialog for selecting the use of default coordinates
    use_defaults, ok = QInputDialog.getItem(
        None,
        "Choosing coordinates",
        "Should I use the default coordinates?",
        ["Yes", "No"],
        editable=False
    )

    if not ok:
        raise ValueError("The selection was canceled.")

    if use_defaults == "Yes":
        north = default_coords["north"]
        west = default_coords["west"]
        south = default_coords["south"]
        east = default_coords["east"]
    else:
        north = get_coordinate("Enter the north coordinate (EPSG:3857):", default_coords["north"])
        west = get_coordinate("Enter the western coordinate (EPSG:3857):", default_coords["west"])
        south = get_coordinate("Enter the southern coordinate (EPSG:3857):", default_coords["south"])
        east = get_coordinate("Enter the eastern coordinate (EPSG:3857):", default_coords["east"])

    # Extracting vertexes of roads and adding them to a layer
    extract_road_vertices_in_area(north, west, south, east)

except ValueError as e:
    QMessageBox.warning(None, "Cancel", str(e))


from qgis.core import (
    QgsProject,
    QgsVectorLayer,
    QgsFeature,
    QgsGeometry,
    QgsField,
    QgsRectangle,
    QgsCoordinateTransform,
    QgsCoordinateReferenceSystem,
    QgsWkbTypes
)
from qgis.PyQt.QtCore import QVariant
from qgis.PyQt.QtWidgets import QMessageBox, QInputDialog

# Default coordinates (EPSG:3857)
default_coords = {
    "north": 7306903.5433,
    "west": 4220660.2002,
    "south": 7294319.0911,
    "east": 4239298.3101
}

# A function for getting the coordinates of an area through dialog boxes
def get_user_coords():
    coords = default_coords.copy()
    try:
        for direction in ["north", "south", "west", "east"]:
            value, ok = QInputDialog.getDouble(
                None,
                f"Enter the coordinate {direction.capitalize()}",
                f"The coordinate {direction.capitalize()} (by default {coords[direction]}):",
                coords[direction],
                decimals=4
            )
            if ok:
                coords[direction] = value
            else:
                raise Exception("The input was canceled by the user.")
        return coords
    except Exception as e:
        QMessageBox.warning(None, "Error", f"Default coordinates are used.\n{e}")
        return default_coords

# Показать диалог выбора
use_default = QMessageBox.question(
    None,
    "Choosing the calculation area",
    "Should I use the default coordinates?",
    QMessageBox.Yes | QMessageBox.No,
    QMessageBox.Yes
) == QMessageBox.Yes

coords = default_coords if use_default else get_user_coords()

# Create a rectangular area for filtering
bounding_box = QgsRectangle(coords["west"], coords["south"], coords["east"], coords["north"])

# Name of the source layer
layer_name = "road_layer"

# Find a layer in the project
layer = QgsProject.instance().mapLayersByName(layer_name)
if not layer:
    QMessageBox.critical(None, "Error", f"Layer '{layer_name}' not faund!")
    raise Exception(f"Layer '{layer_name}' not faund!")
layer = layer[0]

# Make sure that the layer is linear
if layer.geometryType() != QgsWkbTypes.LineGeometry:
    QMessageBox.critical(None, "Error", "The selected layer must be linear!")
    raise Exception("The selected layer must be linear!")

# Create a vector layer for intersection points
intersection_layer = QgsVectorLayer(
    "Point?crs=EPSG:3857",
    "Intersection_Points", #Intersection_Points
    "memory"
)
provider = intersection_layer.dataProvider()

# Add margins to the points layer
provider.addAttributes([
    QgsField("line1_id", QVariant.String),  # ID of the first line
    QgsField("line2_id", QVariant.String),  # ID of the second line
    QgsField("x", QVariant.Double),         # X coordinate
    QgsField("y", QVariant.Double)          # Y coordinate
])
intersection_layer.updateFields()

# Find intersections
features = list(layer.getFeatures())
for i, feature1 in enumerate(features):
    for j, feature2 in enumerate(features):
        # Avoid duplication and checking the same line with yourself
        if i >= j:
            continue

        geom1 = feature1.geometry()
        geom2 = feature2.geometry()

        # Find the intersection points
        if geom1.intersects(geom2):
            intersection = geom1.intersection(geom2)

            # If the intersection is a point (or several points)
            if intersection.type() == QgsWkbTypes.PointGeometry:
                if intersection.isMultipart():
                    points = intersection.asMultiPoint()
                else:
                    points = [intersection.asPoint()]

                for point in points:
                    # Convert coordinates to EPSG:3857
                    crs_transform = QgsCoordinateTransform(
                        layer.crs(),
                        QgsCoordinateReferenceSystem("EPSG:3857"),
                        QgsProject.instance()
                    )
                    point_3857 = crs_transform.transform(point)

                    # Create a new point in the resulting layer
                    new_feature = QgsFeature()
                    new_feature.setGeometry(QgsGeometry.fromPointXY(point_3857))
                    new_feature.setAttributes([
                        feature1["osm_id"],  # ID of the first line
                        feature2["osm_id"],  # ID of the second line
                        point_3857.x(),      # X coordinate
                        point_3857.y()       # Y coordinate
                    ])
                    provider.addFeature(new_feature)

# Deleting points outside the specified rectangle
with edit(intersection_layer):
    for feature in intersection_layer.getFeatures():
        point = feature.geometry().asPoint()
        if not bounding_box.contains(point):
            intersection_layer.deleteFeature(feature.id())

# Добавить слой в проект
QgsProject.instance().addMapLayer(intersection_layer)
QMessageBox.information(None, "Success", "The intersection point layer was created successfully.")


from qgis.core import (
    QgsProject,
    QgsPointXY,
    QgsFeature,
    QgsGeometry,
    QgsField,
    QgsFeatureRequest
)
from PyQt5.QtCore import QVariant
from math import sqrt

# The first script: Adding new points from the "Intersection Points" layer
vertices_layer = QgsProject.instance().mapLayersByName("Normalize_Data")[0]
intersections_layer = QgsProject.instance().mapLayersByName("Intersection_Points")[0] #Intersection Points

if not vertices_layer or not intersections_layer:
    print("One of the layers was not found!")
    exit()

# Editing mode
vertices_layer.startEditing()

# Checking and adding a new 'fid_point' attribute if there is none
if "fid_point" not in [field.name() for field in vertices_layer.fields()]:
    vertices_layer.dataProvider().addAttributes([
        QgsField("fid_point", QVariant.Int)
    ])
    vertices_layer.updateFields()

# Adding new points from the "Intersection Points" layer
new_features = []
for intersection in intersections_layer.getFeatures():
    intersection_point = QgsPointXY(intersection["x"], intersection["y"])
    line1_id = intersection["line1_id"]
    line2_id = intersection["line2_id"]

    for line_id in [line1_id, line2_id]:
        # Checking the existence of a point
        exists = any(
            QgsPointXY(vertex.geometry().asPoint()) == intersection_point
            for vertex in vertices_layer.getFeatures()
        )
        if exists:
            continue

        # Creating a new point
        new_feature = QgsFeature(vertices_layer.fields())
        new_feature.setGeometry(QgsGeometry.fromPointXY(intersection_point))
        new_feature.setAttributes([
            None,  # A new fid will be automatically assigned.
            line_id,
            intersection["x"],
            intersection["y"],
            None,  # vertex_index will be set later
            None,  # fid_point will be set later
        ])
        new_features.append(new_feature)

vertices_layer.addFeatures(new_features)

# Saving changes
if not vertices_layer.commitChanges():
    print("Error saving changes!")
else:
    print("The points were successfully added and processed.")

# Second script: Updating the fid_point attributes
# A function for calculating the distance between two points
def calculate_distance(point1, point2):
    return sqrt((point1.x() - point2.x())**2 + (point1.y() - point2.y())**2)

# Getting a layer "Вершины_дорог_и_перекрёстки_vi"
layer = QgsProject.instance().mapLayersByName('Normalize_Data')[0]

if layer is None:
    print("The layer was not found!")
else:
    layer.startEditing()  # Starting layer editing

    try:
        expression = '"vertex_index" IS NULL'
        request = QgsFeatureRequest().setFilterExpression(expression)

        for feature in layer.getFeatures(request):
            # Coordinates of the current point and its osm_id
            current_point = feature.geometry().asPoint()
            current_osm_id = feature['osm_id']
            current_fid = feature.id()

            # Initialize variables to find the nearest point
            min_distance = float('inf')
            nearest_feature = None

            # Looking for the nearest point with the same osm_id
            for candidate in layer.getFeatures():
                if candidate['osm_id'] == current_osm_id and candidate.id() != feature.id():
                    candidate_point = candidate.geometry().asPoint()
                    distance = calculate_distance(current_point, candidate_point)

                    if distance < min_distance:
                        min_distance = distance
                        nearest_feature = candidate

            # If the nearest point is found, we are looking for a point forming a segment with the nearest
            if nearest_feature:
                nearest_point = nearest_feature.geometry().asPoint()
                nearest_vertex_index = nearest_feature['vertex_index']
                nearest_fid_point = nearest_feature['fid_point']
                nearest_fid = nearest_feature.id()

                if nearest_vertex_index is not None:
                    # Search for the second point forming the segment
                    segment_point = None
                    for candidate in layer.getFeatures():
                        if candidate['osm_id'] == current_osm_id and candidate.id() != nearest_feature.id():
                            candidate_vertex_index = candidate['vertex_index']
                            if candidate_vertex_index is not None:
                                candidate_point = candidate.geometry().asPoint()

                                # Check if the current point lies between the nearest one and the candidate
                                dist_to_current = calculate_distance(current_point, nearest_point)
                                dist_to_candidate = calculate_distance(current_point, candidate_point)
                                dist_between_candidate_and_nearest = calculate_distance(nearest_point, candidate_point)

                                if abs(dist_to_current + dist_to_candidate - dist_between_candidate_and_nearest) < 1e-9:
                                    segment_point = candidate
                                    break

                    if segment_point:
                        segment_vertex_index = segment_point['vertex_index']
                        segment_fid_point = segment_point['fid_point']
                        segment_fid = segment_point.id()

                        # Calculating a new vertex_index for the current point
                        new_vertex_index = max(nearest_vertex_index, segment_vertex_index)

                        # Setting the fid_point of the current point
                        current_fid_point = nearest_fid_point or segment_fid_point
                        layer.changeAttributeValue(current_fid, layer.fields().lookupField('fid_point'), current_fid_point)

                        # We increase the vertex_index for all points if it is greater than or equal to new_vertex_index.
                        for candidate in layer.getFeatures():
                            if candidate['osm_id'] == current_osm_id and candidate['vertex_index'] is not None:
                                if candidate['vertex_index'] >= new_vertex_index:
                                    layer.changeAttributeValue(candidate.id(), layer.fields().lookupField('vertex_index'), candidate['vertex_index'] + 1)

                        # Setting the vertex_index of the current point
                        layer.changeAttributeValue(current_fid, layer.fields().lookupField('vertex_index'), new_vertex_index)

        layer.commitChanges()  # Saving the changes in the layer
        print("The changes were saved successfully!")

        layer_name = "Normalize_Data"

        # Getting a layer from the project
        layer = None
        for lyr in QgsProject.instance().mapLayers().values():
            if lyr.name() == layer_name:
                layer = lyr
                break

        if not layer:
            raise Exception(f"The layer with the name '{layer_name}' not faund.")

        # We check that the layer is vector and has the necessary fields.
        required_fields = ["fid", "fid_point", "osm_id"]
        for field in required_fields:
            if not layer.fields().indexFromName(field) >= 0:
                raise Exception(f"The layer must contain the attribute '{field}'.")

        # Creating a dictionary to store the largest fid for each osm_id
        osm_id_to_max_fid = {}
        for feature in layer.getFeatures():
            osm_id = feature["osm_id"]
            fid = feature["fid"]
            if osm_id is not None and fid is not None:
                # If osm_id is already in the dictionary, update the value to the maximum value.
                if osm_id not in osm_id_to_max_fid:
                    osm_id_to_max_fid[osm_id] = fid
                else:
                    osm_id_to_max_fid[osm_id] = max(osm_id_to_max_fid[osm_id], fid)

        # Starting editing the layer
        layer.startEditing()

        # Updating the 'fid_point' values
        for feature in layer.getFeatures():
            osm_id = feature["osm_id"]
            if osm_id is not None and osm_id in osm_id_to_max_fid:
                feature["fid_point"] = osm_id_to_max_fid[osm_id]
                layer.updateFeature(feature)

        # Saving the changes
        if layer.commitChanges():
            print(f"The 'fid' values have been successfully updated in the 'fid_point', taking into account the largest 'fid' for each 'osm_id' on the layer '{layer_name}'.")
        else:
            print(f"Couldn't save changes to the layer '{layer_name}'.")


    except Exception as e:
        layer.rollBack()
        print(f"Error: {e}")


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

from qgis.core import (
    QgsProject,
    QgsVectorLayer,
    QgsFeature,
    QgsField,
    QgsCoordinateReferenceSystem,
    QgsGeometry, 
    QgsPointXY,
    QgsDistanceArea, 
    QgsCoordinateTransform
)
from PyQt5.QtCore import QVariant

import math

# Loading the intersection layer
layer_crossings = QgsProject.instance().mapLayersByName("Nodal_points")[0]

# Determining the maximum number of points dynamically
field_names = [field.name() for field in layer_crossings.fields()]
max_points = max([int(name[6:]) for name in field_names if name.startswith("osm_id")], default=1)

# Creating a new layer
crs = QgsCoordinateReferenceSystem("EPSG:3857")
layer_output = QgsVectorLayer("Point?crs=EPSG:3857", "Graph_Of_Roads", "memory")
provider = layer_output.dataProvider()

# Adding fields
def get_crossing_data():
    fields = []
    for i in range(1, max_points + 1):
        fields += [
            QgsField(f"fid{i}", QVariant.Int),
            QgsField(f"osm_id{i}", QVariant.Int),
            QgsField(f"vertex_index{i}", QVariant.Int),
            QgsField(f"fid{i}_pre", QVariant.Int),
            QgsField(f"osm_id{i}_pre", QVariant.Int),
            QgsField(f"vertex_index{i}_pre", QVariant.Int),
            QgsField(f"fid{i}_pos", QVariant.Int),
            QgsField(f"osm_id{i}_pos", QVariant.Int),
            QgsField(f"vertex_index{i}_pos", QVariant.Int),
        ]
    return fields

provider.addAttributes(get_crossing_data())
layer_output.updateFields()

# Uploading data to the dictionary for quick search
osm_id_dict = {}

for f in layer_crossings.getFeatures():
    for i in range(1, max_points + 1):
        osm_id = f[f"osm_id{i}"]
        vertex_index = f[f"vertex_index{i}"]
        
        if osm_id is None or vertex_index is None:
            continue
        
        if osm_id is not None and osm_id != NULL:
            osm_id = int(osm_id)
        else:
            continue  # Skip it if the value is missing.
                                    # Conversion to int to eliminate discrepancies
        if osm_id not in osm_id_dict:
            osm_id_dict[osm_id] = []
        
        osm_id_dict[osm_id].append((vertex_index, f.id(), f))  # Storing (vertex_index, fid, feature)

# Sorting indexed data by vertex_index
for osm_id in osm_id_dict:
    osm_id_dict[osm_id].sort()

# Neighbor search function
def find_neighbors(osm_id, vertex_index):
    if osm_id not in osm_id_dict:
        return None, None
    
    neighbors = osm_id_dict[osm_id]
    pre, pos = None, None
    
    for i, (v_index, fid, feature) in enumerate(neighbors):
        if v_index == vertex_index:
            if i > 0:
                pre = neighbors[i - 1]
            if i < len(neighbors) - 1:
                pos = neighbors[i + 1]
            break

    return pre, pos

# Processing of intersection points
output_features = []
for crossing in layer_crossings.getFeatures():
    new_feature = QgsFeature(layer_output.fields())
    new_feature.setGeometry(crossing.geometry())

    for i in range(1, max_points + 1):
        fid = crossing[f"fid{i}"]
        osm_id = crossing[f"osm_id{i}"]
        vertex_index = crossing[f"vertex_index{i}"]

        if fid is None or osm_id is None or vertex_index is None:
            continue

        if osm_id is not None and osm_id != NULL:
            osm_id = int(osm_id)
        else:
            continue  # Skip it if the value is missing.

        pre, pos = find_neighbors(osm_id, vertex_index)

        new_feature.setAttribute(f"fid{i}", fid)
        new_feature.setAttribute(f"osm_id{i}", osm_id)
        new_feature.setAttribute(f"vertex_index{i}", vertex_index)

        if pre:
            new_feature.setAttribute(f"fid{i}_pre", fid)
            new_feature.setAttribute(f"osm_id{i}_pre", osm_id)
            new_feature.setAttribute(f"vertex_index{i}_pre", pre[0])
        else:
            new_feature.setAttribute(f"fid{i}_pre", None)
            new_feature.setAttribute(f"osm_id{i}_pre", None)
            new_feature.setAttribute(f"vertex_index{i}_pre", None)

        if pos:
            new_feature.setAttribute(f"fid{i}_pos", fid)
            new_feature.setAttribute(f"osm_id{i}_pos", osm_id)
            new_feature.setAttribute(f"vertex_index{i}_pos", pos[0])
        else:
            new_feature.setAttribute(f"fid{i}_pos", None)
            new_feature.setAttribute(f"osm_id{i}_pos", None)
            new_feature.setAttribute(f"vertex_index{i}_pos", None)

    output_features.append(new_feature)

provider.addFeatures(output_features)
layer_output.updateExtents()

QgsProject.instance().addMapLayer(layer_output)

def get_layer_by_name(name):
    layers = QgsProject.instance().mapLayers().values()
    for layer in layers:
        if layer.name() == name:
            return layer
    return None

def transform_point_to_msk72(point):
    """ Converts a point to the MSK-72 coordinate system (EPSG:28411) """
    crs_src = QgsProject.instance().crs()
    crs_dest = QgsCoordinateReferenceSystem("EPSG:28411")  # MSK 72
    transform = QgsCoordinateTransform(crs_src, crs_dest, QgsProject.instance())
    return transform.transform(point)

def calculate_distance(point1, point2):
    """ Calculates the distance between two points in MSK-72"""
    d = QgsDistanceArea()
    d.setEllipsoid('WGS84')
    d.setSourceCrs(QgsCoordinateReferenceSystem("EPSG:28411"), QgsProject.instance().transformContext())
    p1 = transform_point_to_msk72(point1)
    p2 = transform_point_to_msk72(point2)
    return d.measureLine(p1, p2)

def get_sorted_points(fid):
    """ Gets a list of road vertex points sorted by vertex_index """
    points = []
    for vertex in layer_vertices.getFeatures():
        if vertex.attribute("fid_point") == fid:
            index = vertex.attribute("vertex_index")
            point = vertex.geometry().asPoint()
            points.append((index, point))
    return sorted(points, key=lambda x: x[0])  # Sorting by vertex_index

layer_graph = get_layer_by_name("Graph_Of_Roads")
layer_vertices = get_layer_by_name("Normalize_Data")

if not layer_graph or not layer_vertices:
    raise ValueError("One of the layers was not found! Check the names.")

# We'll check and add attributes if there aren't any.
layer_graph.startEditing()
provider = layer_graph.dataProvider()
existing_fields = {field.name() for field in provider.fields()}

for i in range(1, 4):
    for suffix in ["_pre", "_pos"]:
        field_name = f"distance{i}{suffix}"
        if field_name not in existing_fields:
            provider.addAttributes([QgsField(field_name, QVariant.Double)])
            existing_fields.add(field_name)

layer_graph.updateFields()

# Filling in the distance attributes
for feature in layer_graph.getFeatures():
    distances = {}
    for i in range(1, 4):  
        fid = feature.attribute(f"fid{i}")
        vertex_start = feature.attribute(f"vertex_index{i}")
        vertex_end_pre = feature.attribute(f"vertex_index{i}_pre")
        vertex_end_pos = feature.attribute(f"vertex_index{i}_pos")
        fid_pre = feature.attribute(f"fid{i}_pre")
        fid_pos = feature.attribute(f"fid{i}_pos")

        sorted_points = get_sorted_points(fid)  # Getting a list of road points

        # A function for summing distances between consecutive points
        def sum_distances(start, end):
            if start is None or end is None:
                return math.inf
            distance = 0
            for j in range(len(sorted_points) - 1):
                index1, p1 = sorted_points[j]
                index2, p2 = sorted_points[j + 1]
                if index1 >= start and index2 <= end:
                    distance += calculate_distance(p1, p2)
            return distance

        # Calculating distances with a NULL check
        distances[f"distance{i}_pre"] = sum_distances(vertex_end_pre, vertex_start) if fid_pre is not None else math.inf
        distances[f"distance{i}_pos"] = sum_distances(vertex_start, vertex_end_pos) if fid_pos is not None else math.inf
        if distances[f"distance{i}_pos"] is 0:
            distances[f"distance{i}_pos"] = float("inf")
        if distances[f"distance{i}_pre"] is 0:
            distances[f"distance{i}_pre"] = float("inf")
            
    for key, value in distances.items():
        feature.setAttribute(key, value)

    layer_graph.updateFeature(feature)

layer_graph.commitChanges()
print("The update is complete!")


import re
from collections import deque
from qgis.core import (
    QgsProject,
    QgsFeature,
    QgsVectorLayer,
    QgsField,
    QgsGeometry,
    QgsPointXY,
    QgsVectorDataProvider
)
from PyQt5.QtCore import QVariant

# A function for safely obtaining a value.
# If the value is missing (None or the string "NULL"), "inf" is returned
def safe_value(val):
    if val is None or str(val).strip().upper() == "NULL":
        return "inf"
    return val

###########################
# Input data for points
###########################

# Setting the input data for the start and end lines
# (# (the comparison is based on fid 1, osm_id 1, vertex_index1, as local parameters)
start_fid1 = 641 
start_osm_id1 = 320825085
start_vertex_index1 = 0 

end_fid1 = 2526 
end_osm_id1 = 727279923
end_vertex_index1 = 0

# We set a list of intermediate points (if there are none, leave an empty list [])
# Each point is set as a tuple: (fid, osm_id, vertex_index)
intermediate_points = [
    # Example: (500, 123456789, 3)
    # Add the necessary intermediary points here
]

#########################################
# We get a layer and collect all the objects
#########################################

layers = QgsProject.instance().mapLayersByName("Graph_Of_Roads")
if not layers:
    raise Exception('Layer "Graph_Of_Roads" not faund')
layer = layers[0]

# Collecting all the objects of the layer in a list
features = list(layer.getFeatures())

#########################################
# Step 1. Create global nodes.
#########################################
# Each row (object) is a global intersection.
# For each object, we collect all local points (groups with the index {i})
# and create the reverse mapping: for each local point (tuple) -> global node id.
global_nodes = {}    # key: global_id (for example, feat.id ()), value: dict with the 'feature' and 'locals' keys
local_to_global = {} # key: (fid, osm_id, vertex_index), value: global_id

for feat in features:
    gid = feat.id()  # the unique identifier of the global node
    global_nodes[gid] = {'feature': feat, 'locals': []}
    field_names = feat.fields().names()
    for field in field_names:
        m = re.match(r"^fid(\d+)$", field)
        if m:
            i = m.group(1)
            local_tuple = (
                safe_value(feat["fid" + i]),
                safe_value(feat["osm_id" + i]),
                safe_value(feat["vertex_index" + i])
            )
            # If at least one value is missing, we skip such a local point.
            if "inf" in local_tuple:
                continue
            global_nodes[gid]['locals'].append(local_tuple)
            local_to_global[local_tuple] = gid

#########################################
# Шаг 2. Defining global nodes for the start and end points.
#########################################
start_global = None
end_global = None
for gid, data in global_nodes.items():
    for local_tuple in data['locals']:
        if local_tuple == (start_fid1, start_osm_id1, start_vertex_index1):
            start_global = gid
        if local_tuple == (end_fid1, end_osm_id1, end_vertex_index1):
            end_global = gid
    if start_global is not None and end_global is not None:
        break

if start_global is None or end_global is None:
    raise Exception("No global nodes found for the start or end point")

print("A global node has been found to start:", start_global, "and for the finish:", end_global)

#########################################
#Step 2.5. Define global nodes for intermediate points.
#########################################
intermediate_globals = []
for local_tuple in intermediate_points:
    found = None
    for gid, data in global_nodes.items():
        for lt in data['locals']:
            if lt == local_tuple:
                found = gid
                break
        if found is not None:
            break
    if found is None:
        raise Exception("The global node for the midpoint was not found.: " + str(local_tuple))
    intermediate_globals.append(found)

print("Found global nodes for intermediate points:", intermediate_globals)

#########################################
# Step 3. Graph construction.
#########################################
# Graph nodes are global nodes.
# For each local point of the object, we look at the pointers to the next (_pos) and the previous (_pre) point.
# If the pointer gives a tuple, we search for it in local_to_global and create an edge with a weight from the distance field.
graph = {}  # key: global_id, value: list of tuples(target_global_id, weight)
for gid in global_nodes:
    graph[gid] = []

# Function for adding an edge to a graph with symmetric transition
def add_edge(graph, src, dst, weight):
    graph[src].append((dst, weight))
    # Adding a symmetrical edge, if it doesn't exist yet.
    if not any(neighbor == src for neighbor, _ in graph[dst]):
        graph[dst].append((src, weight))

for feat in features:
    gid = feat.id()
    field_names = feat.fields().names()
    for field in field_names:
        m = re.match(r"^fid(\d+)$", field)
        if not m:
            continue
        i = m.group(1)
        local_tuple = (
            safe_value(feat["fid" + i]),
            safe_value(feat["osm_id" + i]),
            safe_value(feat["vertex_index" + i])
        )
        if "inf" in local_tuple:
            continue

        # Wrapper function for processing pointers by suffix (for example, "_pos" or "_pre")
        def process_direction(suffix):
            key = "fid" + i + suffix
            if key in field_names and feat[key] is not None:
                neighbor_tuple = (
                    safe_value(feat["fid" + i + suffix]),
                    safe_value(feat["osm_id" + i + suffix]),
                    safe_value(feat["vertex_index" + i + suffix])
                )
                if "inf" in neighbor_tuple:
                    return
                if neighbor_tuple in local_to_global:
                    target_gid = local_to_global[neighbor_tuple]
                    weight_field = "distance" + i + suffix
                    weight = feat[weight_field] if weight_field in field_names and feat[weight_field] is not None else 0
                    add_edge(graph, gid, target_gid, weight)
        # We process both directions
        process_direction("_pos")
        process_direction("_pre")


########################################
# Step 4. Find the shortest path using Dijkstra's algorithm with mandatory intermediate points.
#########################################
# The shortest path search function (Dijkstra) between two nodes.
def dijkstra_path(graph, start, goal):
    distances = {node: float('inf') for node in graph}
    previous = {node: None for node in graph}
    distances[start] = 0
    unvisited = set(graph.keys())
    
    while unvisited:
        current = min(unvisited, key=lambda node: distances[node])
        if distances[current] == float('inf'):
            break  # The other nodes are unreachable
        if current == goal:
            break
        unvisited.remove(current)
        for neighbor, weight in graph.get(current, []):
            alt = distances[current] + weight
            if alt < distances[neighbor]:
                distances[neighbor] = alt
                previous[neighbor] = current
                
    if distances[goal] == float('inf'):
        return None, None
    path = []
    node = goal
    while node is not None:
        path.insert(0, node)
        node = previous[node]
    return path, distances[goal]

# If there are no intermediate points, just look for a way from start to finish.
# If there is, we search the path sequentially by segments.:
# start -> intermediate1, intermediate1 -> intermediate2, ..., intermediate_last -> finish.
all_points = [start_global] + intermediate_globals + [end_global]
full_path = []
total_distance = 0
for idx in range(len(all_points)-1):
    seg_start = all_points[idx]
    seg_end = all_points[idx+1]
    seg_path, seg_distance = dijkstra_path(graph, seg_start, seg_end)
    if seg_path is None:
        raise Exception("Couldn't find the path between the nodes {} and {}".format(seg_start, seg_end))
    # When merging, we remove the duplicate point at the junction of the segments.
    if full_path and full_path[-1] == seg_path[0]:
        full_path.extend(seg_path[1:])
    else:
        full_path.extend(seg_path)
    total_distance += seg_distance

print("Shortest path found (global nodes):", full_path)
print("Total route distance:", total_distance)

#########################################
# Step 5. Formation of the output layer "rout_obletchiki".
#########################################
# For each global point on the path, we use the geometry of the original feature.
crs = layer.crs().authid()
mem_layer = QgsVectorLayer("Point?crs=" + crs, "Dijkstra_path", "memory")
pr = mem_layer.dataProvider()
pr.addAttributes(layer.fields())
mem_layer.updateFields()

new_features = []
for gid in full_path:
    feat = global_nodes[gid]['feature']
    new_feat = QgsFeature()
    new_feat.setGeometry(feat.geometry())
    new_feat.setAttributes(feat.attributes())
    new_features.append(new_feat)

pr.addFeatures(new_features)
QgsProject.instance().addMapLayer(mem_layer)
