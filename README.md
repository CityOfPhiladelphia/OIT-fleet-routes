![](https://beta.phila.gov/media/20170531115900/oit-logo1.png)

# Operationalizing Streets GPS/AVL data


## Introduction

The city manages a fleet of garbage trucks that collect residential and commercial waste on a daily basis; these trucks are equipped with GPS tracking devices that collect information such as location and speed during the trucks' routes.  OIT aims to develop a real-time visual representation of a truck's location and route, providing assurance that a truck is efficiently collecting garbage and allowing departure/arrival time to be estimated.

The primary challenge to visualizing a truck's location is the limited quality of GPS data, which poses two problems.  First, GPS systems often produce an inaccurate coordinate location, and regardless of accuracy, some method must be implemented to determine which street a truck is on based on its GPS coordinates.  Such identification methods are known as map-matching algorithms.  Second, GPS data is often collected at such a low frequency that there is ambiguity regarding the streets taken in between two points.  Some method must be implemented to determine the most probable set of streets taken in these cases.

This guide documents the steps required to utilize OIT's map-matching algorithm.

## Data Inputs

OIT collects GPS tracking information from active garbage trucks every 15 seconds using Verizon NetworkFleet.  The Esri GeoEvent Server receives the GPS data and this data is harvested every 5 minutes.

The data collected from GeoEvent Server is exported to a CSV file, which is processed using Python script(s).  After associations, the resulting data and street centerline data are input to ArcMap for visualization.

Inputs include: vehicle id number (VIN), longitude/latitude, average speed, time stamp, street names, street segment topology

## Python Scripts

* **preprocess_fleetpts.py** is a preprocessing script, which entails sorting the data by truck ID, calculating distance and heading changes, and removing excess data points that represent idle time.
* **associations.py** is a reproduction of Esri's Near geoprocessing function, computing the shortest distance point-to-line association for each GPS point.  According to Esri:
>the shortest distance from a point to a line segment is the perpendicular to the line segment. If a perpendicular cannot be drawn within the end vertices of the line segment, then the distance to the closest end vertex is the shortest distance.

* **weight_based_model.py** selects a set of candidate points as the map-matched route.  Weight-based models can provide a simple means of customizing the set of properties included in a scoring, the way the subscore is calculated for each property, and the relative weighting of each subscore in calculation of the total score.  This model scores all possible paths of candidate points for a given route, because the likelihood of a certain link being correct is partially dependent upon its relation to the preceding and following candidate links chosen.
* 

## Guide to Processing Script (sliding_window.py)

1. **read_roadnetwork()** 
   - The centerline and topology files are read.  
   - A CENTERLINE object is created for each centerline, and these objects are stored in the centerline_map dictionary.  
   - A SEGPOINT object is created for each x,y coordinate that is an endpoint in a centerline segment, and these objects are stored in the coord_to_segpoints dictionary.
2. **extract_centerline_points()**
   - The LineStrings of all the centerlines in centerline_map are parsed to retrieve each centerline coordinate.  They are saved in centerline_endpoints_map so that we can have a mapping of coordinates to centerlines, and a complete list of these coordinates (the map keyset)
   - If any two points on a centerline LineString have a Euclidean distance of over 500 ft, we append the midpoint as an additional point on this centerline, so that the k-d tree will more accurately determine which gps points are close to this centerline.
3. A **k-d tree** is made, indexing all of the centerline points, including the additional ones.  This is all the preprocessing needed before gps points can be looked at.
4. **preprocess_fleetpoints.run()** is called, to do the initial preprocessing script written a few months ago
5. **read_fleetpoints()**
   - Read each row of the new fleetpoint csv and immediately process it, rather than reading and collecting all the rows and then processing them all together.
   - This was implemented as a way of simulating receiving the data in real-time.  We can imagine that read_fleetpoints will be called every 15 minutes or so when we receive a new csv.
   - a FLEETPOINT object is created for each data point, and stored in the ongoing list fleetpoints_list.  Then, we run identify_candidate_segments on the this point to determine the candidate associations to later choose from
   - **identify_candidate_segments()**
     * Queries the tree for the 6 centerline endpoints that are closest to the gps point (Euclidean distance).  For each of these endpoints, it considers each centerline that includes the endpoint.  It parses each centerline’s linestring and collects each sequential pair of linestring coordinates as a candidate segment.
     * For each of these candidate segments, it calculates the point on this segment line that is closest to the gps point (reproduction of Esri’s Near function).  It calculates the equation of a line perpendicular to the centerline segment that goes through the original gps point, and calculates where these two lines intersect as the association point.
     * These candidate points are stored as CANDIDATE objects in candidate_associations, where they are mapped to the id of the gps point.
     * If, for some gps point, any of its associations are within a Euclidean distance of 75 feet from it, then we will no longer consider any associations that are over 250 feet away from it, because it is highly unlikely that such points would be a correct match. (Saves time)
   - Once the list of candidate association points is generated for a gps point, we look at the current **sliding window** for its route.
     * If this point is the 1st, 2nd, 3rd, or 4th in its route, then we will wait for more points to be collected before selecting associations.  
     * If this point is the 5th in its route, we run init_path() to select associations for the first 3 gps points.  
     * If this point is 6th or greater in its route, we run continue_path() to select an association for the third to last (i.e. n-2) gps point in its route
     * sp() is a method used by both init_path() and continue_path() to generate the candidate scores.
5. **init_path()**
   - Begin by scoring all the candidate associations for the first gps point in the route, calling score_candidate.  The scoring is weight-based, meaning that each candidate is sub-scored on several factors, and the overall score is a sum of these sub-scores which are each given distinct weights based on their relative importance in making a correct selection.  The lower the overall score, the more likely a candidate is the correct selection.
   - For the first point in the route, the only factor considered is the candidate point’s distance from the gps point.  The scores are saved in path_dist.  
   - Run sp() for points 2, 3, 4, and 5:
     * For each candidate of the current gps point, we calculate a temporary score based on each of the candidates of the previous gps point: we consider the Euclidean distance between these two candidates, whether the two candidates are on the same street, and whether the current candidate’s centerline segment is directly reachable from the previous candidate’s centerline segment.  The candidate of the previous point that produces the lowest score for a current point candidate will be stored in the pred dictionary as this current point candidate’s predecessor, and this current point candidate’s path_dist value will be the cumulative sum of its score plus its predecessor’s score.
   - Finally, we identify the candidate of point 5 with the lowest path_dist value, and backtrack through its predecessors to select the candidates of points 1, 2, and 3 on its path as correct matchings.  The list of these points is added to associations_list.
6. **continue_path()**
   - This method works very similarly to init_path(); however, it begins by scoring the second-to-last point candidates using the third-to-last point candidates’ path_dist values that have been previously computed.  We score the last point candidates, select its candidate with the lowest path_dist value, and backtrack through its predecessors to select the candidate of the third-to-last point on its path as the correct matching.  This point is added to associations_list.
7. The process of reading a gps point from the csv and running init_path() or continue_path() accordingly continues until all rows have been read.  Then, in practice we would run finish_paths() on any unfinished routes that have not received new gps points since the previous csv file.  However, in this simulation we only have one file so we run finish_paths() to finish up all of the routes of the current csv file.  
8. **finish_paths()**
   - Call init_path() on a route if the route has fewer than 5 points and none of these points have yet been given matchings.  
   - Otherwise, look at the path_dist values that have been generated for the final two points that need to be matched, and give them matchings.
9. Finally, we call **write_associations()** to write all new matchings made from this batch of data to a new csv, and we clear associations_list so that these matchings will not be written again.


## Getting Started

Place the CSV file in the same folder as the Python script or change the script's path.  Ensure that the CSV column fields match those in the script's class instantiations.  Run the script and an output CSV file will be created.  Files can be viewed in ArcMap as follows:

1.	In ArcMap, Insert -> Data Frame.  Right click the data frame -> Add Data -> locate data file
2.	Right click data csv file -> Display XY Data -> X Field: lon, Y Field: lat -> Edit -> Geographic Coordinate Systems -> coordinate system of choice
3.	Include a street centerline Feature Class file
4.	Once data export is obtained as a new layer, delete the Event Data Layer
5.	Search for Points to Line geoprocess -> For Input Features, select the exported file.  Line Field is runid, Sort Field is runid_sequence
6.	Right click the Points to Line export -> Properties -> Symbology -> Categories -> Value Field: runid -> Add All Values
7.	Select specific routes by clicking the Select Features option in the toolbar
