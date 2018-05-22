# Operationalizing Streets GPS/AVL data

![](https://beta.phila.gov/media/20170531115900/oit-logo1.png)



[TOCM]


##Introduction

The city manages a fleet of garbage trucks that collect residential and commercial waste on a daily basis; these trucks are equipped with GPS tracking devices that collect information such as location and speed during the trucks' routes.  OIT aims to develop a real-time visual representation of a truck's location and route, providing assurance that a truck is efficiently collecting garbage and allowing departure/arrival time to be estimated.

The primary challenge to visualizing a truck's location is the limited quality of GPS data, which poses two problems.  First, GPS systems often produce an inaccurate coordinate location, and regardless of accuracy, some method must be implemented to determine which street a truck is on based on its GPS coordinates.  Such identification methods are known as map-matching algorithms.  Second, GPS data is often collected at such a low frequency that there is ambiguity regarding the streets taken in between two points.  Some method must be implemented to determine the most probable set of streets taken in these cases.

This guide documents the steps required to utilize OIT's map-matching algorithm.

##Raw Data

OIT collects GPS tracking information from active garbage trucks every 15 seconds using Verizon NetworkFleet.  The Esri GeoEvent Server receives the GPS data and this data is harvested every 5 minutes.

The data collected from GeoEvent Server is exported to a CSV file, which is processed using Python script(s).  After associations, the resulting data and street centerline data are input to ArcMap for visualization.

##Python Scripts

* **nf_tim1.py** is a preprocessing script, which entails sorting the data by truck ID, calculating distance and heading changes, and removing excess data points that represent idle time.
* **associations.py** is a reproduction of Esri's Near geoprocessing function, computing the shortest distance point-to-line association for each GPS point.  According to Esri:
>the shortest distance from a point to a line segment is the perpendicular to the line segment. If a perpendicular cannot be drawn within the end vertices of the line segment, then the distance to the closest end vertex is the shortest distance.

* **weight_based_model.py** selects a set of candidate points as the map-matched route.  Weight-based models can provide a simple means of customizing the set of properties included in a scoring, the way the subscore is calculated for each property, and the relative weighting of each subscore in calculation of the total score.  This model scores all possible paths of candidate points for a given route, because the likelihood of a certain link being correct is partially dependent upon its relation to the preceding and following candidate links chosen.

##Getting Started

Place the CSV file in the same folder as the Python script or change the script's path.  Ensure that the CSV column fields match those in the script's class instantiations.  Run the script and an output CSV file will be created.  Files can be viewed in ArcMap as follows:

1.	In ArcMap, Insert -> Data Frame.  Right click the data frame -> Add Data -> locate data file
2.	Right click data csv file -> Display XY Data -> X Field: lon, Y Field: lat -> Edit -> Geographic Coordinate Systems -> coordinate system of choice
3.	Include a street centerline Feature Class file
4.	Once data export is obtained as a new layer, delete the Event Data Layer
5.	Search for Points to Line geoprocess -> For Input Features, select the exported file.  Line Field is runid, Sort Field is runid_sequence
6.	Right click the Points to Line export -> Properties -> Symbology -> Categories -> Value Field: runid -> Add All Values
7.	Select specific routes by clicking the Select Features option in the toolbar