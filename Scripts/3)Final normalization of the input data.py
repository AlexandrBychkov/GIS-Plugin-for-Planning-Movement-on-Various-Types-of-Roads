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
