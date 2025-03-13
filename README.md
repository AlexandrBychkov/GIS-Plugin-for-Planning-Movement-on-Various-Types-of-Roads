# GIS-Plugin-for-Planning-Movement-on-Various-Types-of-Roads
GIS Plugin for Planning Movement on Various Types of Roads
# Connected Weighted Road Graph Implementation

This repository is open source and contains an example implementation of constructing a connected weighted road graph. This implementation is presented for the article **"GIS Plugin for Planning Movement on Various Types of Roads"**, which was showcased at the IEEE EDM 2025 conference.

## Repository Structure

- **All_in_one**: Contains all necessary scripts as separate files.
- **Example**: Contains everything required for a test run.

## Tested Environment

- **QGIS v3.22**

## Demonstration Layers

The following layers in the **Road_graph** group are provided as examples of the expected outcome:
- `Normalize_Data`
- `Intersection_Points`
- `Nodal_points`
- `Graph_Of_Roads`
- `Dijkstra_path`

**Note:** Before running individual scripts or the "Example" script, these demonstration layers must be removed.

## Maps for Spatial Orientation

For proper spatial orientation, it is recommended to add:
- OpenStreetMap
- Google Maps (satellite version)

For the test area, these maps are located in the **maps** folder.

## Example Script

The **Example** script is a combination of all the individual scripts and is intended for the simplest demonstration of the implementation.

## Usage

1. **Remove Demonstration Layers:** Delete the demonstration layers before executing the scripts.
2. **Add Map Layers:** Add the necessary OpenStreetMap and Google Maps (satellite) layers from the **maps** folder for correct spatial orientation.
3. **Run the Scripts:** Execute the individual scripts from the **All_in_one** folder or run the **Example** script for a full demonstration.

## License

This project is open source.
