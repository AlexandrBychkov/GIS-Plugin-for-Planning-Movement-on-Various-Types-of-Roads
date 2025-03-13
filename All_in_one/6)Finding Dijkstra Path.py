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
