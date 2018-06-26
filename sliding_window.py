from scipy import spatial
from collections import defaultdict
import os
import csv
import sys
from math import sqrt, pow
from shapely.wkt import loads

import preprocess_fleetpts

cwd = os.path.dirname(__file__)
centerlines_file = 'centerline_shape_2272'
topology_file = 'roadnetwork_topology'
original_fleetpoints_file = 'fleet_2018_06_01-00'
preprocessed_fleetpoints_file = '06_01-00_processed'
outfile = '6_01-00_assoc'

centerline_map = defaultdict(list)  # centerline id -> centerline object
centerline_endpoints_map = defaultdict(list)  # coordinates of point in centerline -> centerline object
id_to_segpoints = defaultdict(list)  # id of segpoint -> segpoint object
coord_to_segpoints = defaultdict(list)  # integer-casted coordinates (in case of rounding during processing) -> segpoint object

new_fleetpoint_list = []
saved_fleetpoints = []  # imagine this is a list of fleetpoints from the previous call to this script (retrieved from geoeventserver), which still have to be matched to centerlines
routes_list = defaultdict(list)  # route runid -> list of fleetpoint objects in route.  For real time this should be retrieved from geoeventserver and keys should be VIN
candidate_associations = defaultdict(list)  # orig fleetpoint id -> list of candidate point objects
path_dist = defaultdict(list)  # coordinate of original point -> lowest score among its candidate paths
pred = defaultdict(list)  # coordinate of original point -> candidate object of predecessor on its best candidate path
associations_list = defaultdict(list)  # route runid -> list of candidate objects that are chosen as matches for each point in route


# represents a centerline segment (one row of centerline csv)
class CENTERLINE:
    def __init__(self, row):
        self.linestring = loads(row[11])
        self.seg_id = row[9]
        self.street_name = row[1].strip()


# represents a point on one or more centerline segments, retrieved from csv's LINESTRING field
class SEGPOINT:
    def __init__(self, row):
        self.segpoint_id = int(row[0])
        self.lon = float(row[1])
        self.lat = float(row[2])
        self.in_segs = set(map(int, row[3].split(' ')))
        if row[4] != '':
            self.goes_to = list(map(int, row[4].split(' ')))
        else:
            self.goes_to = ''


# represents a data point from NetworkFleet
class FLEETPOINT:
    def __init__(self, row, id):
        self.id = id
        self.runid = int(row[1].strip())
        self.runid_sequence = int(row[2].strip())
        self.vin = row[3].strip()
        self.datetime = row[4].strip()
        self.lon = float(row[22])
        self.lat = float(row[23])


# represents a coordinate point that may be chosen as a matching for a fleetpoint
class CANDIDATE:
    def __init__(self, assocx, assocy, dist, seg, origx, origy):
        self.assocx = assocx
        self.assocy = assocy
        self.near_dist = dist
        self.seg_id = seg
        self.origx = origx
        self.origy = origy


def csv_path(file_name):
    return os.path.join(cwd, file_name + '.csv')


# creates lists of centerline and segpoint objects from the csv files
def read_roadnetwork():
    path = csv_path(centerlines_file)
    f = open(path, 'r')
    i = 0
    try:
        reader = csv.reader(f)
        for row in reader:
            if i == 0:  # first row contains the header
                i += 1
                continue
            r = CENTERLINE(row)
            centerline_map[r.seg_id] = r
            i += 1
    except IOError:
        print('Error opening ' + path, sys.exc_info()[0])
    f.close()

    path = csv_path(topology_file)
    f = open(path, 'r')
    i = 0
    try:
        reader = csv.reader(f)
        for row in reader:
            if i == 0:
                i += 1
                continue
            r = SEGPOINT(row)
            id_to_segpoints[r.segpoint_id] = r
            # if two points have same integer casting, merge their connectivity
            if coord_to_segpoints[int(r.lon), int(r.lat)]:
                coord_to_segpoints[int(r.lon), int(r.lat)].in_segs.union(r.in_segs)
                coord_to_segpoints[int(r.lon), int(r.lat)].goes_to.extend(r.goes_to)
            else:
                coord_to_segpoints[int(r.lon), int(r.lat)] = r
            i += 1
    except IOError:
        print('Error opening ' + path, sys.exc_info()[0])
    f.close()
    return


# maps each point contained in a centerline segment LineString to every centerline segment it belongs to, and returns a list of all such points
def extract_centerline_points():
    for seg_id in centerline_map:
        centerline = centerline_map[seg_id]
        centerline_coords = list(centerline.linestring.coords)
        l = 0
        centerline_endpoints_map[centerline_coords[l][0], centerline_coords[l][1]].append(centerline)
        # traverse each sequential pair of points in the LineString, as forming separate lines
        while l < len(centerline_coords) - 1:
            segx1 = centerline_coords[l][0]
            segy1 = centerline_coords[l][1]
            segx2 = centerline_coords[l + 1][0]
            segy2 = centerline_coords[l + 1][1]
            # if this piece of the centerline is longer than 500ft, create an intermediate point
            if sqrt(pow(segx2 - segx1, 2) + pow(segy2 - segy1, 2)) > 500:
                segx3 = min(segx1, segx2) + abs(segx1 - segx2) / 2
                segy3 = min(segy1, segy2) + abs(segy1 - segy2) / 2
                centerline_endpoints_map[segx3, segy3].append(centerline)
            l += 1
            centerline_endpoints_map[segx2, segy2].append(centerline)

    return list(centerline_endpoints_map.keys())


# simulates receiving a batch of data points
def read_fleetpoints():
    path = csv_path(preprocessed_fleetpoints_file)
    f = open(path, 'r')
    i = 0
    if len(saved_fleetpoints) == 0:
        id = 0
    # else:
        # retrieve id of most recent timestamp fleetpoint, and increment by 1
    try:
        reader = csv.reader(f)
        for row in reader:
            if i == 0:
                i += 1
                continue
            # create a gps point object from the csv row
            r = FLEETPOINT(row, id)
            new_fleetpoint_list.append(r)

            # identify potential centerline segment mappings based on lat/lon
            identify_candidate_segs(r.lon, r.lat, id)

            # perform sliding window function based on the number of points that have already been collected for this route
            if r.runid in routes_list:
                route = routes_list[r.runid]
                route.append(r)
                if len(route) == 5:
                    init_path(route)
                elif len(route) > 5:
                    continue_path(route)
            else:
                routes_list[r.runid].append(r)

            id += 1
            i += 1
    except IOError:
        print('Error opening ' + path, sys.exc_info()[0])
    f.close()
    return


# identify potential centerline segment matchings based on coordinate
def identify_candidate_segs(fleetptx, fleetpty, fleetpt_id):
    candidate_segs = []
    checked_centerlines = set()  # keep track of the segments we have already considered for this particular point

    # query the tree for the 6 closest points to the gps point, returns their index in centerline_endpoints_list
    closest_coords_indices = tree.query([fleetptx, fleetpty], k=6)[1]

    for coord_index in closest_coords_indices:
        if coord_index != len(centerline_endpoints_list):  # else indicates fewer than 6 points returned
            # get the list of centerlines that include this coordinate
            candidate_centerlines = centerline_endpoints_map[centerline_endpoints_list[coord_index]]

            for centerline in candidate_centerlines:
                centerline_id = centerline.seg_id
                if centerline_id not in checked_centerlines:
                    checked_centerlines.add(centerline_id)
                    centerline_coords = list(centerline.linestring.coords)
                    l = 0
                    while l < len(centerline_coords) - 1:  # each line within this centerline's LineString becomes a candidate
                        segx1 = centerline_coords[l][0]
                        segy1 = centerline_coords[l][1]
                        segx2 = centerline_coords[l + 1][0]
                        segy2 = centerline_coords[l + 1][1]
                        candidate_segs.append([segx1, segy1, segx2, segy2, centerline_id])
                        l += 1

    # Use Near to calculate lat/lon of candidate points.
    # If any of these are within 75 ft of original coordinate, get rid of candidates that are over 250 ft away from original
    within_75 = False
    for seg in candidate_segs:
        near_results = near(fleetptx, fleetpty, seg)
        if near_results.near_dist <= 75:
            within_75 = True
        candidate_associations[fleetpt_id].append(near_results)
    if within_75:
        j = 0
        while j < len(candidate_associations[fleetpt_id]):
            if candidate_associations[fleetpt_id][j].near_dist > 250:
                del candidate_associations[fleetpt_id][j]
            else:
                j += 1


# reproduction of Esri's Near function: “the shortest distance from a point to a line segment is the perpendicular to the line segment. If a
# perpendicular cannot be drawn within the end vertices of the line segment, then the distance to the closest end vertex is the shortest distance.”
def near(fleetptx, fleetpty, seg):
    segx1 = seg[0]
    segy1 = seg[1]
    segx2 = seg[2]
    segy2 = seg[3]
    seg_id = seg[4]

    # check if segment has the same lat or lon for both endpoints
    if segx2 == segx1:
        assocx = segx1
        assocy = fleetpty
    elif segy2 == segy1:
        assocx = fleetptx
        assocy = segy1
    else:
        segSlope = (segy2 - segy1) / (segx2 - segx1)
        segB = segy2 - segSlope * segx2
        perpendicularSlope = -1 / segSlope
        pB = fleetpty - perpendicularSlope * fleetptx
        assocx = (segB - pB) / (perpendicularSlope - segSlope)
        assocy = perpendicularSlope * assocx + pB
    # if a perpendicular cannot be drawn within the end vertices of the line segment:
    if not (((segx1 <= assocx <= segx2) & ((segy2 <= assocy <= segy1) | (segy1 <= assocy <= segy2))) |
            ((segx2 <= assocx <= segx1) & ((segy2 <= assocy <= segy1) | (segy1 <= assocy <= segy2)))):
        if sqrt(pow(fleetptx - segx1, 2) + pow(fleetpty - segy1, 2)) < sqrt(pow(fleetptx - segx2, 2) + pow(fleetpty - segy2, 2)):
            assocx = segx1
            assocy = segy1
        else:
            assocx = segx2
            assocy = segy2

    dist = sqrt(pow(fleetptx - assocx, 2) + pow(fleetpty - assocy, 2))
    return CANDIDATE(assocx, assocy, dist, seg_id, fleetptx, fleetpty)


# scores current candidate point on several factors, including its likelihood of coming after this previous candidate point (if not first in route)
def score_candidate(current, previous=None):
    # score Near distance of the current point
    if current.near_dist / 160 > 1:
        emiss = 1
    else:
        emiss = current.near_dist / 160

    # if first point in route, score based on Near distance only
    if previous is None:
        return emiss
    else:
        previous_centerline = centerline_map[previous.seg_id]
        current_centerline = centerline_map[current.seg_id]

        # score the selection of both the current and next point based on Euclidean distance
        trans_dist = sqrt(pow(previous.assocx - current.assocx, 2) + pow(previous.assocy - current.assocy, 2))
        orig_dist = sqrt(pow(previous.origx - current.origx, 2) + pow(previous.origy - current.origy, 2))
        if trans_dist / 750 > 1 or orig_dist > 500:
            trans_dist = 1
        else:
            trans_dist = trans_dist / 750

        # improve the score if the points are on the same street
        #if previous_centerline.street_name == current_centerline.street_name:
        #    trans_dist -= 0.25

        # improve the score if the next point is on a segment that can be reached within one or two hops from the current point
        previous_endpts = list(previous_centerline.linestring.coords)
        previous_endpts = [coord_to_segpoints[int(previous_endpts[0][0]), int(previous_endpts[0][1])],
                        coord_to_segpoints[int(previous_endpts[len(previous_endpts) - 1][0]), int(previous_endpts[len(previous_endpts) - 1][1])]]
        current_endpts = list(current_centerline.linestring.coords)
        current_endpts = [coord_to_segpoints[int(current_endpts[0][0]), int(current_endpts[0][1])].segpoint_id,
                       coord_to_segpoints[int(current_endpts[len(current_endpts) - 1][0]), int(current_endpts[len(current_endpts) - 1][1])].segpoint_id]
        connected = False
        for previous_pt in previous_endpts:
            for seg in previous_pt.goes_to:
                if seg in current_endpts:
                    connected = True
                else:
                    for next_seg in id_to_segpoints[seg].goes_to:
                        if next_seg in current_endpts:
                            connected = True
        if connected:
            trans_dist -= 0.25

        return emiss + trans_dist


# scores all candidates for a given point (time represented by fleetpoint id) trying all combinations of previous and current points.
# The best score for the last fleetpoint corresponds to the best path of candidates
def sp(start_time, end_time, scored=False):
    if not scored:
        # find shortest path for each candidate point from 2nd original to final original
        time = start_time
        while time < end_time:
            for previous in candidate_associations[time]:
                # print('CURRENT', previous.assocx, previous.assocy, path_dist[previous.assocx, previous.assocy, time])
                for current in candidate_associations[time + 1]:
                    score = score_candidate(current, previous)
                    # if the path to the next point using current point is the shortest thus far for this next point, use it
                    if path_dist[current.assocx, current.assocy, time + 1] > path_dist[previous.assocx, previous.assocy, time] + score:
                        path_dist[current.assocx, current.assocy, time + 1] = path_dist[previous.assocx, previous.assocy, time] + score
                        pred[current.assocx, current.assocy, time + 1] = previous
                    # print('FOR', current.assocx, current.assocy, ': ', pred[current.assocx, current.assocy, time + 1], score)
            time += 1

    # find shortest distance among all candidates of last original point; this corresponds to shortest overall path so far
    final_score = sys.float_info.max
    current_pred = []
    for candidate in candidate_associations[end_time]:
        score = path_dist[candidate.assocx, candidate.assocy, end_time]
        if score < final_score:
            current_pred = candidate
            final_score = score
    return current_pred


# scores points on a route for the first time
def init_path(route):
    start_time = route[0].id
    end_time = start_time + len(route) - 1
    associations = []

    # initialize all scores to high number
    i = start_time + 1
    while i <= end_time:
        for candidate in candidate_associations[i]:
            path_dist[candidate.assocx, candidate.assocy, i] = sys.float_info.max
        i += 1

    # candidate points for the first original are scored only by the Near distance
    for candidate in candidate_associations[start_time]:
        path_dist[candidate.assocx, candidate.assocy, start_time] = score_candidate(candidate)

    # score the rest of the candidate points and return best path of candidates
    current_pred = sp(start_time, end_time)
    time = end_time

    # retrieve all candidates on this shortest path
    if end_time - start_time == 4:  # if there are 5 points in the route, we will start by matching only the first 3
        while time > start_time:
            current_pred = pred[current_pred.assocx, current_pred.assocy, time]
            if time <= start_time + 3:
                associations.insert(0, current_pred)
            time -= 1
    else:  # if the route has ended with less than 5 points, we will match all of them now
        associations.insert(0, current_pred)
        while time > start_time:
            current_pred = pred[current_pred.assocx, current_pred.assocy, time]
            associations.insert(0, current_pred)
            time -= 1

    associations_list[route[0].runid] = associations


# if we have obtained 6 or more points for a route path, we use our current information to match the third to last point
def continue_path(route):
    start_time = route[len(route) - 1].id - 3
    end_time = route[len(route) - 1].id

    # reset scores for points 3 and 4, initialize score for point 5
    for i in [1, 2, 3]:
        for c in candidate_associations[start_time + i]:
            path_dist[c.assocx, c.assocy, start_time + i] = sys.float_info.max

    # score candidates for these points
    current_pred = sp(start_time, end_time)
    time = end_time

    # retrieve third-to-last candidate on this shortest path
    while time > start_time + 1:
        current_pred = pred[current_pred.assocx, current_pred.assocy, time]
        time -= 1

    associations_list[route[0].runid].append(current_pred)


# once we have determined that a route path is finished, complete all remaining map matches
# this would be a good place to check the time of the last point and decide whether to wait for more points on the route
def finish_paths():
    for run_id in routes_list:
        route = routes_list[run_id]

        if len(route) < 5:  # if we have not matched any of the points on this route path
            init_path(route)
        else:  # if we only have the final two points remaining to match on this route path
            end_time = route[len(route) - 1].id
            associations = []
            current_pred = sp(0, end_time, True)

            associations.insert(0, current_pred)  # found match for final point
            current_pred = pred[current_pred.assocx, current_pred.assocy, end_time]
            associations.insert(0, current_pred)  # found match for second to last point
            associations_list[route[0].runid].extend(associations)


def write_associations():
    path = csv_path(outfile)
    f = open(path, 'w+')
    header = 'id,runid,runid_sequence,vin,datetime,lon,lat,assoc_lon,assoc_lat,assoc_seg_id\n'
    f.write(header)
    run_id = 0
    j = 0
    while run_id < len(associations_list):
        route = associations_list[run_id]
        for assoc in route:
            s = '%i,%i,%i,%s,%s,%f,%f,%f,%f,%s\n' % (
                j,
                new_fleetpoint_list[j].runid,
                new_fleetpoint_list[j].runid_sequence,
                new_fleetpoint_list[j].vin,
                str(new_fleetpoint_list[j].datetime),
                new_fleetpoint_list[j].lon,
                new_fleetpoint_list[j].lat,
                assoc.assocx,
                assoc.assocy,
                assoc.seg_id,
            )
            f.write(s)
            j += 1
        run_id += 1
    f.close()


# EXECUTION
read_roadnetwork()
centerline_endpoints_list = extract_centerline_points()
tree = spatial.KDTree(centerline_endpoints_list)
print('Finished preprocessing road network.  Now preprocessing Verizon data')

preprocess_fleetpts.run(original_fleetpoints_file, preprocessed_fleetpoints_file)
print('Finished preprocessing Verizon data. Now simulating real-time association')


read_fleetpoints()
# at whatever point we determine it is time to "close off" a path, run finish_paths().  Here it is when there are no more data points to read.
# this should only finish paths that have not received new points since the previous csv.
finish_paths()
print('Finished simulation. Now writing')

write_associations()
