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
