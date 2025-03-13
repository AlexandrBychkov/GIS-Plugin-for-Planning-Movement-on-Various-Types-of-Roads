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
