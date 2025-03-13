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
