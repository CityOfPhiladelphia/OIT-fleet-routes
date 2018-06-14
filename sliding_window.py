from scipy import spatial
from collections import defaultdict
import os
import csv
import sys
from math import sqrt, pow
from shapely.wkt import loads
import preprocess_fleetpts


centerlines_file = 'centerline_shape_2272'
topology_file = '5_17_18_topology'
original_fleetpoints_file = '06_01'
preprocessed_fleetpoints_file = 'did' #'XY06_01_originals_feet'
outfile = '06_12_18_refactored_2'
cwd = os.path.dirname(__file__)

fleetpoint_list = []
routes_list = defaultdict(list)  # route runid -> list of fleetpoint objects in route
associations_list = defaultdict(list)
candidate_associations = defaultdict(list)  # orig fleetpoint id -> list of candidate objects
candidate_distances = defaultdict(list)
path_dist = defaultdict(list)  # coordinates of originals -> shortest dist
pred = defaultdict(list)  # coordinates of originals -> predecessor

# things that could be made into a csv
centerline_map = defaultdict(list)  # centerline id -> centerline object
centerline_endpoints_map = defaultdict(list)  # coordinates of point in centerline -> centerline object
coord_to_segpoints = defaultdict(list)  # integer-casted coordinates (in case of rounding during processing) -> segpoint object


class CENTERLINE:
    def __init__(self, row):
        self.linestring = loads(row[11])
        self.seg_id = row[9]
        self.street_name = row[1].strip()


class SEGPOINT:
    def __init__(self, row):
        self.segpoint_id = int(row[0])
        self.lon = float(row[1])
        self.lat = float(row[2])
        self.in_segs = list(map(int, row[3].split(' ')))
        if row[4] != '':
            self.goes_to = list(map(int, row[4].split(' ')))
        else:
            self.goes_to = ''


class FLEETPOINT:
    def __init__(self, row):
        self.id = int(row[0].strip())
        self.runid = int(row[1].strip())
        self.runid_sequence = int(row[2].strip())
        self.vin = row[3].strip()
        self.datetime = row[4].strip()
        self.lon = float(row[22])
        self.lat = float(row[23])


class CANDIDATE:
    def __init__(self, assocx, assocy, dist, seg, origx, origy):
        self.assocx = assocx
        self.assocy = assocy
        self.near_dist = dist
        self.seg = seg
        self.origx = origx
        self.origy = origy


def csv_path(file_name):
    return os.path.join(cwd, file_name + '.csv')


def read_centerlines():  # create a list of centerline objects from the csv file
    path = csv_path(centerlines_file)
    f = open(path, 'r')
    i = 0
    try:
        reader = csv.reader(f)
        for row in reader:
            if i == 0:  # because the first row contains the header
                i += 1
                continue
            r = CENTERLINE(row)
            centerline_map[r.seg_id] = r
            i += 1
    except IOError:
        print('Error opening ' + path, sys.exc_info()[0])
    f.close()
    return


def read_topology():
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
            coord_to_segpoints[int(r.lon), int(r.lat)] = r
            i += 1
    except IOError:
        print('Error opening ' + path, sys.exc_info()[0])
    f.close()
    return


# maps each point contained in a centerline segment LineString to every centerline segment it belongs to, and returns a list of all such points
def extract_centerline_points():
    points = []
    for seg_id in centerline_map:
        centerline = centerline_map[seg_id]
        centerline_coords = list(centerline.linestring.coords)
        l = 0
        centerline_endpoints_map[centerline_coords[l][0], centerline_coords[l][1]].append(centerline)
        while l < len(centerline_coords) - 1:  # consider each line that is a part of the current centerline LineString
            segx1 = centerline_coords[l][0]
            segy1 = centerline_coords[l][1]
            segx2 = centerline_coords[l + 1][0]
            segy2 = centerline_coords[l + 1][1]
            if sqrt(pow(segx2 - segx1, 2) + pow(segy2 - segy1,
                                                2)) > 500:  # if a piece of the centerline is longer than 500ft, include an intermediate point
                segx3 = min(segx1, segx2) + abs(segx1 - segx2) / 2
                segy3 = min(segy1, segy2) + abs(segy1 - segy2) / 2
                centerline_endpoints_map[segx3, segy3].append(centerline)
                points.extend([(segx3, segy3)])
            l += 1
            centerline_endpoints_map[segx2, segy2].append(centerline)
        points.extend(centerline_coords)
    return list(set(points))


# simulates receiving data points one at a time.
def read_fleetpoints():
    path = csv_path(preprocessed_fleetpoints_file)
    f = open(path, 'r')
    i = 0
    try:
        reader = csv.reader(f)
        for row in reader:
            if i == 0:
                i += 1
                continue
            # create a gps point object from the csv row
            if row[3].strip() == '1HTWGAAT96J346882':
                r = FLEETPOINT(row)
                fleetpoint_list.append(r)

                # identify potential centerline segment mappings based on coordinate
                identify_candidate_segs(r.lon, r.lat, r.id)

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
            i += 1
    except IOError:
        print('Error opening ' + path, sys.exc_info()[0])
    print('ROUTES LIST', routes_list)
    f.close()
    return


# identify potential centerline segment mappings based on coordinate
def identify_candidate_segs(fleetptx, fleetpty, fleetpt_id):
    candidate_segs = []
    closest_coords_indices = tree.query([fleetptx, fleetpty], k=6)[
        1]  # query the tree for the 6 closest points to the gps point, returns their index in centerline_points
    checked_centerlines = set()  # keep track of the segments we have already considered for this particular point
    for coord_index in closest_coords_indices:
        # print('new candidate coordinate')
        if coord_index != len(centerline_endpoints_list):  # else indicates fewer than 6 points returned
            candidate_centerlines = centerline_endpoints_map[
                centerline_endpoints_list[coord_index]]  # find the list of centerlines that each coordinate is mapped to
            for centerline in candidate_centerlines:
                seg_id = centerline.seg_id
                # print('current centerline:', linestring)
                if seg_id not in checked_centerlines:
                    # print('this centerline has not been checked yet')
                    checked_centerlines.add(seg_id)
                    centerline_coords = list(centerline.linestring.coords)
                    l = 0
                    while l < len(centerline_coords) - 1:  # consider each line that is a part of the current centerline LineString
                        segx1 = centerline_coords[l][0]
                        segy1 = centerline_coords[l][1]
                        segx2 = centerline_coords[l + 1][0]
                        segy2 = centerline_coords[l + 1][1]
                        # print('Segment:', [segx1, segy1, segx2, segy2])
                        candidate_segs.append([segx1, segy1, segx2, segy2, seg_id])
                        l += 1
                # print('already checked segment')

    # if we have potential mappings within 75 ft of original coordinate, do not consider potentials that are over 250 ft away from original
    within_75 = False
    for candidate in candidate_segs:
        near_results = near(fleetptx, fleetpty, candidate)
        if near_results[2] <= 75:
            within_75 = True
        print(fleetpt_id, candidate)
        candidate_associations[fleetpt_id].append(near_results)
    if within_75:
        j = 0
        while j < len(candidate_associations[fleetpt_id]):
            if candidate_associations[fleetpt_id][j][2] > 250:
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
    # if a perpendicular can be drawn within the end vertices of the line segment:
    if not (((segx1 <= assocx <= segx2) & ((segy2 <= assocy <= segy1) | (segy1 <= assocy <= segy2))) | \
            ((segx2 <= assocx <= segx1) & ((segy2 <= assocy <= segy1) | (segy1 <= assocy <= segy2)))):
        if sqrt(pow(fleetptx - segx1, 2) + pow(fleetpty - segy1, 2)) < sqrt(pow(fleetptx - segx2, 2) + pow(fleetpty - segy2, 2)):
            assocx = segx1
            assocy = segy1
        else:
            assocx = segx2
            assocy = segy2
    dist = sqrt(pow(fleetptx - assocx, 2) + pow(fleetpty - assocy, 2))
    return assocx, assocy, dist, seg_id, fleetptx, fleetpty


def score_candidate(current, next=None):
    current_near_dist = current[2]
    if next is None:
        if current_near_dist / 160 > 1:
            score = 1
        else:
            score = current_near_dist / 160
    else:
        current_lon = current[0]
        current_lat = current[1]
        current_centerline = centerline_map[current[3]]
        next_lon = next[0]
        next_lat = next[1]
        next_near_dist = next[2]
        next_centerline = centerline_map[next[3]]
        # score the next candidate point based on Near distance
        if next_near_dist / 160 > 1:
            emiss = 1
        else:
            emiss = next_near_dist / 160
        # score the selection of both the current and next point based on Euclidean distance
        trans_dist = sqrt(pow(current_lon - next_lon, 2) + pow(current_lat - next_lat, 2))
        orig_dist = sqrt(pow(current[4] - next[4], 2) + pow(current[5] - next[5], 2))
        if trans_dist / 750 > 1 or orig_dist > 500:
            trans_dist = 1
        else:
            trans_dist = trans_dist / 750
        if current_centerline.street_name == next_centerline.street_name:
            trans_dist -= 0.25
        segpts = list(current_centerline.linestring.coords)
        segpts = [coord_to_segpoints[int(segpts[0][0]), int(segpts[0][1])],
                  coord_to_segpoints[int(segpts[len(segpts) - 1][0]), int(segpts[len(segpts) - 1][1])]]
        nextpts = list(next_centerline.linestring.coords)
        nextpts = [coord_to_segpoints[int(nextpts[0][0]), int(nextpts[0][1])].segpoint_id,
                   coord_to_segpoints[int(nextpts[len(nextpts) - 1][0]), int(nextpts[len(nextpts) - 1][1])].segpoint_id]
        connected = False
        for segpt in segpts:
            for seg in segpt.goes_to:
                if seg in nextpts:
                    connected = True
        if connected:
            trans_dist -= 0.25
        # print('emiss',emiss,2*trans_dist)
        score = emiss + trans_dist
    return score


def sp(start_time, end_time, scored=False):
    if not scored:
        # find shortest path for each candidate point from 2nd original to final original
        time = start_time
        while time < end_time:
            # print('TIME', time, candidate_associations[time])
            for current in candidate_associations[time]:
                current_lon = current[0]
                current_lat = current[1]
                # print('CURRENT', current_lon, current_lat, path_dist[current_lon, current_lat, time])
                for next in candidate_associations[time + 1]:
                    next_lon = next[0]
                    next_lat = next[1]
                    score = score_candidate(current, next)
                    # if the path to the next point using current point is the shortest thus far for this next point, use it
                    if path_dist[next_lon, next_lat, time + 1] > path_dist[current_lon, current_lat, time] + score:
                        path_dist[next_lon, next_lat, time + 1] = path_dist[current_lon, current_lat, time] + score
                        pred[next_lon, next_lat, time + 1] = [current_lon, current_lat]
                    # print('FOR', next_lon, next_lat, ': ', pred[next_lon, next_lat, time + 1], path_dist[next_lon, next_lat, time + 1])
            time += 1

    # find shortest distance among all candidates of last original point; this corresponds to shortest overall path so far
    final_score = sys.float_info.max
    current_pred = []
    for c in candidate_associations[end_time]:
        sc = path_dist[c[0], c[1], end_time]
        if sc < final_score:
            current_pred = [c[0], c[1]]
            final_score = sc
    return current_pred


def init_path(route):
    # print('Init Path')
    # print(route)
    start_time = route[0].id
    end_time = start_time + len(route) - 1
    associations = []

    i = start_time + 1
    while i <= end_time:
        for c in candidate_associations[i]:
            path_dist[c[0], c[1], i] = sys.float_info.max
        i += 1

    # candidate points for the first original are scored only by the Near distance
    for c in candidate_associations[start_time]:
        lon = c[0]
        lat = c[1]
        path_dist[lon, lat, start_time] = score_candidate(c)

    current_pred = sp(start_time, end_time)
    time = end_time

    # retrieve all candidates on this shortest path
    if end_time - start_time == 4:  # if there are 5 points, we will start by matching the first 3
        while time > start_time:
            current_lon = current_pred[0]
            current_lat = current_pred[1]
            current_pred = pred[current_lon, current_lat, time]
            if time <= start_time + 3:
                associations.insert(0, current_pred)
            time -= 1
    else:  # if the route has ended with less than 5 points, we will match all of them now
        associations.insert(0, current_pred)
        while time > start_time:
            current_lon = current_pred[0]
            current_lat = current_pred[1]
            current_pred = pred[current_lon, current_lat, time]
            associations.insert(0, current_pred)
            time -= 1

    associations_list[route[0].runid] = associations
    # print('assoc', route[0].runid, associations_list[route[0].runid])


# if we have obtained 6 or more points for a route path, we use our current information to match the third to last point
def continue_path(route):
    # print('Continue Path')
    start_time = route[len(route) - 1].id - 3
    end_time = route[len(route) - 1].id

    for i in [1, 2, 3]:
        for c in candidate_associations[start_time + i]:
            path_dist[c[0], c[1], start_time + i] = sys.float_info.max

    current_pred = sp(start_time, end_time)
    time = end_time

    # retrieve third-to-last candidate on this shortest path
    while time > start_time + 1:
        current_lon = current_pred[0]
        current_lat = current_pred[1]
        current_pred = pred[current_lon, current_lat, time]
        time -= 1

    associations_list[route[0].runid].append(current_pred)
    # print('assoc', associations_list[route[0].runid])


# once we have determined that a route path is finished, complete all remaining map matches
def finish_paths():
    for run_id in routes_list:
        route = routes_list[run_id]

        if len(route) < 5:  # if we have not matched any of the points on this route path
            # print('route',route,'runid',run_id)
            init_path(route)
        else:  # if we only have the final two points remaining to match on this route path
            # print('Fin Path')
            end_time = route[len(route) - 1].id
            associations = []

            current_pred = sp(0, end_time, True)

            associations.insert(0, current_pred)  # found match for final point
            current_lon = current_pred[0]
            current_lat = current_pred[1]
            current_pred = pred[current_lon, current_lat, end_time]
            associations.insert(0, current_pred)  # found match for second to last point
            associations_list[route[0].runid].extend(associations)
            # print('assoc', route[0].runid, associations_list[route[0].runid])


def write_associations():
    path = csv_path(outfile)
    f = open(path, 'w+')
    header = 'id,runid,runid_sequence,vin,datetime,lon,lat,assoc_lon,assoc_lat\n'
    f.write(header)
    run_id = 0
    j = 0
    while run_id < len(associations_list):
        route = associations_list[run_id]
        for assoc in route:
            s = '%i,%i,%i,%s,%s,%f,%f,%f,%f\n' % (
                j,
                fleetpoint_list[j].runid,
                fleetpoint_list[j].runid_sequence,
                fleetpoint_list[j].vin,
                str(fleetpoint_list[j].datetime),
                fleetpoint_list[j].lon,
                fleetpoint_list[j].lat,
                assoc[0],
                assoc[1]
            )
            f.write(s)
            j += 1
        run_id += 1
    f.close()


read_centerlines()
read_topology()
# things that could be made into a csv
centerline_endpoints_list = extract_centerline_points()
tree = spatial.KDTree(centerline_endpoints_list)
print('Finished preprocessing road network.  Now preprocessing Verizon data')

#preprocess_fleetpts.run(original_fleetpoints_file, preprocessed_fleetpoints_file)
print('Finished preprocessing Verizon data. Now simulating real-time association')
read_fleetpoints()
# at whatever point we determine it is time to "close off" a path, run finish_paths().  Here it is when there are no more data points to read.
# this should only finish paths that have not received new points since the previous csv.
finish_paths()
print('Finished simulation. Now writing')

write_associations()
associations_list.clear()
