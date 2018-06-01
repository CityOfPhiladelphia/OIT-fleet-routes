from scipy import spatial
from collections import defaultdict
import os
import csv
import sys
from math import sqrt, pow
from shapely.geometry import Point, LineString
from shapely.wkt import loads


centerlines_file = 'centerline_shape_2272'
fleetpoints_file = 'originals_feet'
outfile = '6_1_18_2'
cwd = os.path.dirname(__file__)

centerline_map = defaultdict(list)
fleetpoint_list = []
routes_list = defaultdict(list)
associations_list = defaultdict(list)
candidate_associations = defaultdict(list)  # mapping is orig point -> candidate points, so orig index -> [lon, lat, perp dist]
candidate_distances = defaultdict(list)
path_dist = defaultdict(list)  # coordinates of originals -> shortest dist
pred = defaultdict(list)  # coordinates of originals -> predecessor


class CENTERLINE:
    def __init__(self, row):
        self.linestring = loads(row[11])
        self.seg_id = row[9]
        self.street_name = row[1].strip()


class FLEETPOINT:
    def __init__(self, row):
        self.id = int(row[0].strip())
        self.runid = int(row[1].strip())
        self.runid_sequence = int(row[2].strip())
        self.vin = row[3].strip()
        self.datetime = row[4].strip()
        self.lon = float(row[22])
        self.lat = float(row[23])
        self.point = Point(self.lon, self.lat)


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


# simulates receiving data points one at a time
def read_fleetpoints():  # create a list of gps point objects from the csv file
    path = csv_path(fleetpoints_file)
    f = open(path, 'r')
    i = 0
    try:
        reader = csv.reader(f)
        for row in reader:
            if i == 0:
                i += 1
                continue
            r = FLEETPOINT(row)
            fleetpoint_list.append(r)

            fleetptx = r.lon
            fleetpty = r.lat
            candidate_segs = identify_candidate_segs(fleetptx, fleetpty)
            for candidate in candidate_segs:
                near_results = near(fleetptx, fleetpty, candidate)
                candidate_associations[r.id].append(near_results)
            if r.runid in routes_list:
                route = routes_list[r.runid]
                route.append(r)
                if len(route) == 5:
                    init_path(route)
                elif len(route) > 5:
                    continue_path(route)
            else:
                routes_list[r.runid] = []
                routes_list[r.runid].append(r)

            i += 1
    except IOError:
        print('Error opening ' + path, sys.exc_info()[0])
    print('ROUTES LIST', routes_list[r.runid])
    f.close()
    finish_paths()
    return


# maps each point that is a endpoint for a centerline segment to every centerline segment it is an endpoint for, and returns a list of all the points
def extract_centerline_points():
    points = []
    for seg_id in centerline_map:
        centerline = centerline_map[seg_id]
        centerline_coords = list(centerline.linestring.coords)
        for coord in centerline_coords:
            centerline_endpoints_map[coord].append(centerline)
        # print(map[coord])
        points.extend(centerline_coords)
    return list(set(points))


def identify_candidate_segs(fleetptx, fleetpty):  # , distance_upper_bound=500
    candidate_segs = []
    closest_coords_indices = tree.query([fleetptx, fleetpty], k=4)[1]  # query the tree for the 4 closest points to the gps point, returns their index in centerline_points
    checked_centerlines = set()  # keep track of the segments we have already considered for this particular point
    k = 0  # useful in print statements for keeping track of the coordinates we are checking
    for coord_index in closest_coords_indices:
        # print('new candidate coordinate')
        if coord_index != len(centerline_endpoints_list):
            candidate_centerlines = centerline_endpoints_map[centerline_endpoints_list[coord_index]]  # find the list of LineStrings that each coordinate is mapped to
            for centerline in candidate_centerlines:
                seg_id = centerline.seg_id
                linestring = centerline.linestring
                # print(k)
                # print('current centerline:', linestring)
                if str(linestring) not in checked_centerlines:
                    # print('this centerline has not been checked yet')
                    checked_centerlines.add(str(linestring))
                    l = 0
                    centerline_coords = list(linestring.coords)
                    while l < len(centerline_coords) - 1:  # consider each line that is a part of the current centerline LineString
                        segx1 = centerline_coords[l][0]
                        segy1 = centerline_coords[l][1]
                        segx2 = centerline_coords[l + 1][0]
                        segy2 = centerline_coords[l + 1][1]
                        # print('Segment:', [segx1, segy1, segx2, segy2])
                        candidate_segs.append([segx1, segy1, segx2, segy2, seg_id])
                        l += 1
                # print('already checked segment')
                k += 1
    return candidate_segs


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
    if (((segx1 <= assocx <= segx2) & ((segy2 <= assocy <= segy1) | (segy1 <= assocy <= segy2))) | \
            ((segx2 <= assocx <= segx1) & ((segy2 <= assocy <= segy1) | (segy1 <= assocy <= segy2)))) == False:
        if sqrt(pow(fleetptx - segx1, 2) + pow(fleetpty - segy1, 2)) < sqrt(pow(fleetptx - segx2, 2) + pow(fleetpty - segy2, 2)):
            assocx = segx1
            assocy = segy1
        else:
            assocx = segx2
            assocy = segy2
    dist = sqrt(pow(fleetptx - assocx, 2) + pow(fleetpty - assocy, 2))
    return assocx, assocy, dist, seg_id


def score_candidate(current, next=None):
    current_perp = current[2]
    if next is None:
        if current_perp / 160 > 1:
            score = 1
        else:
            score = current_perp / 160
    else:
        current_lon = current[0]
        current_lat = current[1]
        current_centerline = centerline_map[current[3]]
        next_lon = next[0]
        next_lat = next[1]
        next_perp = next[2]
        next_centerline = centerline_map[next[3]]
        # score the next candidate point based on Near distance
        if next_perp / 160 > 1:
            emiss = 1
        else:
            emiss = next_perp / 160
        # score the selection of both the current and next point based on Euclidean distance
        trans_dist = sqrt(pow(current_lon - next_lon, 2) + pow(current_lat - next_lat, 2))
        if trans_dist / 750 > 1:
            trans_dist = 1
        else:
            trans_dist = trans_dist / 750
        if current_centerline.street_name == next_centerline.street_name:
            trans_dist -= 0.25
        score = emiss + 3*trans_dist
    return score


def init_path(route):
    # print('Init Path')
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

    # find shortest path for each candidate point from 2nd original to 4th original
    time = start_time
    while time < end_time:
        #print('TIME', time, candidate_associations[time])
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
                #print('FOR', next_lon, next_lat, ': ', pred[next_lon, next_lat, time + 1], path_dist[next_lon, next_lat, time + 1])
        time += 1

    # find shortest distance among all candidates of last original point; this corresponds to shortest overall path so far
    final_score = sys.float_info.max
    current_pred = []
    for c in candidate_associations[time]:
        sc = path_dist[c[0], c[1], time]
        if sc < final_score:
            current_pred = [c[0], c[1]]
            final_score = sc

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
    #print('Continue Path')
    start_time = route[len(route)-1].id - 3
    end_time = route[len(route)-1].id

    for i in [1, 2, 3]:
        for c in candidate_associations[start_time + i]:
            path_dist[c[0], c[1], start_time + i] = sys.float_info.max

    time = start_time
    while time < end_time:
        #print('TIME', time, candidate_associations[time])
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
                #print('FOR', next_lon, next_lat, ': ', pred[next_lon, next_lat, time + 1], path_dist[next_lon, next_lat, time + 1])
        time += 1

    final_assoc = sys.float_info.max
    time = end_time
    # find shortest distance among all candidates of last original point; this corresponds to shortest overall path
    current_pred = []
    for c in candidate_associations[time]:
        sc = path_dist[c[0], c[1], time]
        if sc < final_assoc:
            current_pred = [c[0], c[1]]
            final_assoc = sc
    # retrieve all candidates on this shortest path
    while time > start_time:
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
            init_path(route)
        else:  # if we only have the final two points remaining to match on this route path
            # print('Fin Path')
            end_time = route[len(route) - 1].id
            time = end_time
            associations = []

            # find shortest distance among all candidates of last original point; this corresponds to shortest overall path
            current_pred = []
            final_assoc = sys.float_info.max
            for c in candidate_associations[time]:
                sc = path_dist[c[0], c[1], time]
                if sc < final_assoc:
                    current_pred = [c[0], c[1]]
                    final_assoc = sc

            associations.insert(0, current_pred)  # found match for final point
            current_lon = current_pred[0]
            current_lat = current_pred[1]
            current_pred = pred[current_lon, current_lat, time]
            associations.insert(0, current_pred)  # found match for second to last point
            associations_list[route[0].runid].extend(associations)
            # print('assoc', route[0].runid, associations_list[route[0].runid])


def write_associations():
    path = csv_path(outfile)
    f = open(path, 'w+')
    header = 'id,runid,runid_sequence,vin,datetime,lon,lat,assoc_lon,assoc_lat\n'
    f.write(header)
    print(len(associations_list))
    run_id = 0
    j = 0
    while run_id < len(associations_list):
        print(run_id)
        route = associations_list[run_id]
        for assoc in route:
            # print('assoc', assoc)
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
print('Finished reading centerlines. Now indexing centerlines')

centerline_endpoints_map = defaultdict(list)
centerline_endpoints_list = extract_centerline_points()
tree = spatial.KDTree(centerline_endpoints_list)
print('Finished indexing.  Now simulating real-time association')

read_fleetpoints()
print('Finished simulation. Now writing')

write_associations()
