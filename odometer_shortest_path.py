from scipy import spatial
from collections import defaultdict
import os
import csv
import sys
from math import sqrt, pow
from shapely.geometry import Point, LineString
from shapely.wkt import loads
from sortedcontainers import SortedDict


centerlines_file = 'centerlines_joined'
topology_file = '5_17_18_topology'
fleetpoints_file = 'originals_feet'
outfile = '6_5_18_2'
cwd = os.path.dirname(__file__)

centerline_map = defaultdict(list)
fleetpoint_list = []
associations_list = []
candidate_associations = defaultdict(list)  # mapping is orig point -> candidate points, so orig index -> [lon, lat, perp dist]
candidate_distances = defaultdict(list)
id_to_segpoints = defaultdict(list)
coord_to_segpoints = defaultdict(list)


class CENTERLINE:
    def __init__(self, row):
        self.linestring = loads(row[24])
        self.seg_id = row[17]
        self.street_name = row[23].strip()
        self.dir = row[18].strip()


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
        self.odometer_dist = float(row[15])
        self.lon = float(row[22])
        self.lat = float(row[23])
        self.point = Point(self.lon, self.lat)


class CANDIDATE:
    def __init__(self, assocx, assocy, dist, seg):
        self.assocx = assocx
        self.assocy = assocy
        self.near_dist = dist
        self.seg = seg


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
            id_to_segpoints[r.segpoint_id] = r
            coord_to_segpoints[int(r.lon), int(r.lat)] = r
            i += 1
    except IOError:
        print('Error opening ' + path, sys.exc_info()[0])
    f.close()
    return


def read_fleetpoints(): # create a list of gps point objects from the csv file
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
            i += 1
    except IOError:
        print('Error opening ' + path, sys.exc_info()[0])
    f.close()
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
    return CANDIDATE(assocx, assocy, dist, seg)


def odometer_sp(current, time, next):
    odometer = fleetpoint_list[time + 1].odometer_dist * 5280
    visited_dist = SortedDict()
    coord_to_dist = defaultdict(list)
    visited_pred = defaultdict(list)

    # calculate the distance from current candidate to closest logical centerline endpoint
    start_seg = centerline_map[current.seg[4]] # obtain the complete centerline object for the current candidate's segment
    if start_seg.dir == 'TF':
        start_dist = sqrt(pow(current.seg[1] - current.assocy, 2) + pow(current.seg[0] - current.assocx, 2))
        visited_dist[start_dist] = coord_to_segpoints[int(current.seg[0]), int(current.seg[1])]
        coord_to_dist[current.seg[0], current.seg[1]] = start_dist
    elif start_seg.dir == 'FT':
        start_dist = sqrt(pow(current.seg[3] - current.assocy, 2) + pow(current.seg[2] - current.assocx, 2))
        visited_dist[start_dist] = coord_to_segpoints[int(current.seg[2]), int(current.seg[3])]
        coord_to_dist[current.seg[2], current.seg[3]] = start_dist
    else:
        start_dist = sqrt(pow(current.seg[1] - current.assocy, 2) + pow(current.seg[0] - current.assocx, 2))
        visited_dist[start_dist] = coord_to_segpoints[int(current.seg[0]), int(current.seg[1])]
        coord_to_dist[current.seg[0], current.seg[1]] = start_dist
        start_dist = sqrt(pow(current.seg[3] - current.assocy, 2) + pow(current.seg[2] - current.assocx, 2))
        visited_dist[start_dist] = coord_to_segpoints[int(current.seg[2]), int(current.seg[3])]
        coord_to_dist[current.seg[2], current.seg[3]] = start_dist

    end_seg = centerline_map[next.seg[4]]
    if end_seg.dir == 'TF':
        end_segpoint1 = next.seg[2], next.seg[3]
        end_dist1 = sqrt(pow(current.seg[3] - current.assocy, 2) + pow(current.seg[2] - current.assocx, 2))
        end_segpoint2 = None
    elif end_seg.dir == 'FT':
        end_segpoint1 = next.seg[0], next.seg[1]
        end_dist1 = sqrt(pow(current.seg[1] - current.assocy, 2) + pow(current.seg[0] - current.assocx, 2))
        end_segpoint2 = None
    else:
        end_segpoint1 = next.seg[0], next.seg[1]
        end_dist1 = sqrt(pow(current.seg[1] - current.assocy, 2) + pow(current.seg[0] - current.assocx, 2))
        end_segpoint2 = next.seg[2], next.seg[3]
        end_dist2 = sqrt(pow(current.seg[3] - current.assocy, 2) + pow(current.seg[2] - current.assocx, 2))

    current_segpoint = visited_dist.popitem(0)
    current_dist = current_segpoint[0]
    current_segpoint = current_segpoint[1]

    i = 0
    while (current_segpoint.lon, current_segpoint.lat) != end_segpoint1 and (current_segpoint.lon, current_segpoint.lat) != end_segpoint2 and current_dist < 2000:
        print(time, i)
        for nextpoint in current_segpoint.goes_to:
                nextpoint = id_to_segpoints[nextpoint]
                add_dist = sqrt(pow(nextpoint.lon - current_segpoint.lon, 2) + pow(nextpoint.lat - current_segpoint.lat, 2))
                # could identify seg id here by comparing current and next in_segs
                if (nextpoint.lon, nextpoint.lat) in coord_to_dist:
                    next_dist = coord_to_dist[nextpoint.lon, nextpoint.lat]
                    if next_dist > current_dist + add_dist:
                        coord_to_dist[nextpoint.lon, nextpoint.lat] = current_dist + add_dist
                        del visited_dist[next_dist]
                        visited_dist[current_dist + add_dist] = nextpoint
                        visited_pred[nextpoint.lon, nextpoint.lat] = current_segpoint.lon, current_segpoint.lat
                else:
                    coord_to_dist[nextpoint.lon, nextpoint.lat] = current_dist + add_dist
                    visited_dist[current_dist + add_dist] = nextpoint
                    visited_pred[nextpoint.lon, nextpoint.lat] = current_segpoint.lon, current_segpoint.lat
        current_segpoint = visited_dist.popitem(0)
        current_dist = current_segpoint[0]
        current_segpoint = current_segpoint[1]
        i += 1

    if current_dist >= 2000:
        return 1
    elif [current_segpoint.lon, current_segpoint.lat] != end_segpoint1:
        current_dist += end_dist1
        return current_dist / odometer
    else:
        current_dist += end_dist2
        return current_dist / odometer


def score_candidate(current, time=None, next=None):
    if next is None:
        if current.near_dist / 160 > 1:
            score = 1
        else:
            score = current.near_dist / 160
    else:
        current_centerline = centerline_map[current.seg[4]]
        next_centerline = centerline_map[next.seg[4]]
        # score the next candidate point based on Near distance
        if next.near_dist / 160 > 1:
            emiss = 1
        else:
            emiss = next.near_dist / 160
        # score the selection of both the current and next point based on Euclidean distance
        trans_dist = sqrt(pow(current.assocx - next.assocx, 2) + pow(current.assocy - next.assocy, 2))
        if trans_dist / 2000 > 1:
            trans_dist = 1
        else:
            trans_dist = odometer_sp(current, time, next)
        if current_centerline.street_name == next_centerline.street_name:
            trans_dist -= 0.25
        score = emiss + 3*trans_dist
    return score


# find the best score selection among all candidate points using a shortest path algorithm.  Only does 1 route per method call
def path(route):
    path_dist = defaultdict(list)  # coordinates of originals -> shortest dist
    pred = defaultdict(list)  # coordinates of originals -> predecessor

    for i in candidate_associations:
        for c in candidate_associations[i]:
            path_dist[c.assocx, c.assocy, i] = sys.float_info.max

    # candidate points for the first original are scored only by the Near distance
    for c in candidate_associations[route[0]]:
        path_dist[c.assocx, c.assocy, route[0]] = score_candidate(c)

    # find shortest path for each candidate point from 2nd original onward
    time = route[0]
    while time < route[0] + len(route) - 1:
        print('TIME', time, candidate_associations[time])
        for current in candidate_associations[time]:
            #print('CURRENT', current.assocx, current.assocy, path_dist[current.assocx, current.assocy, time])
            for next in candidate_associations[time+1]:
                score = score_candidate(current, time, next)
                # if the path to the next point using current point is the shortest thus far for this next point, use it
                if path_dist[next.assocx, next.assocy, time + 1] > path_dist[current.assocx, current.assocy, time] + score:
                    path_dist[next.assocx, next.assocy, time + 1] = path_dist[current.assocx, current.assocy, time] + score
                    pred[next.assocx, next.assocy, time + 1] = [current.assocx, current.assocy]
                print('FOR', next.assocx, next.assocy, ': ', pred[next.assocx, next.assocy, time+1], path_dist[next.assocx, next.assocy, time+1])
        time += 1

    final_assoc = sys.float_info.max
    time = route[0] + len(route) - 1
    # find shortest distance among all candidates of last original point; this corresponds to shortest overall path
    associations = []
    for c in candidate_associations[time]:
        sc = path_dist[c.assocx, c.assocy, time]
        print('final', time, sc)
        if sc < final_assoc:
            associations.clear()
            associations.append([c.assocx, c.assocy])
            final_assoc = sc
    print('final assoc', associations[0])
    current_pred = associations[0]
    # retrieve all candidates on this shortest path
    while time > route[0]:
        print(time)
        current_lon = current_pred[0]
        current_lat = current_pred[1]
        current_pred = pred[current_lon, current_lat, time]
        print(current_pred)
        associations.insert(0, current_pred)
        print('current list', associations)
        time -= 1

    print('ASSOCIATIONS')
    print(associations)
    associations_list.extend(associations)


def associate():
    j = 0
    run_id = 0
    route = []
    print('Making candidate Near associations for each original point')
    while j < len(fleetpoint_list):
        fleetptx = fleetpoint_list[j].lon
        fleetpty = fleetpoint_list[j].lat
        #print([fleetptx, fleetpty])
        candidate_segs = identify_candidate_segs(fleetptx, fleetpty)
        for candidate in candidate_segs:
            near_results = near(fleetptx, fleetpty, candidate)
            candidate_associations[j].append(near_results)
        if fleetpoint_list[j].runid != run_id:
            path(route)
            route.clear()
            run_id = fleetpoint_list[j].runid
        route.append(j)
        j += 1
    print('last route', route)
    path(route) # last route
    print('Finished making Near associations. Now finding best candidate path')


def write_associations():
    path = csv_path(outfile)
    f = open(path, 'w+')
    header = 'id,runid,runid_sequence,vin,datetime,lon,lat,assoc_lon,assoc_lat\n'
    f.write(header)
    print(len(associations_list))
    j=0
    while j < len(fleetpoint_list):
        assoc = associations_list[j]
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
    f.close()


read_fleetpoints()
read_centerlines()
read_topology()
print('Finished reading input files. Now indexing centerlines')

centerline_endpoints_map = defaultdict(list)
centerline_endpoints_list = extract_centerline_points()
tree = spatial.KDTree(centerline_endpoints_list)
print('Finished indexing.  Now associating')

associate()
print('Finished associating. Now writing')

write_associations()
