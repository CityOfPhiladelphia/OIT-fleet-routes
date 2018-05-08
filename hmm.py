from scipy import spatial
from collections import defaultdict
import shapely.geos
import os
import csv
import sys
from math import sqrt, pow
from shapely.geometry import Point, LineString
from shapely.wkt import loads


centerlines_file = 'centerline_shape_2272'
fleetpoints_file = 'hmm_test_route'
outfile = '5_8_18_2'
cwd = os.path.dirname(__file__)

centerline_list = []
fleetpoint_list = []
associations_list = []
candidate_associations = defaultdict(list)  # mapping is orig point -> candidate points, so orig index -> [lon, lat, perp dist]
candidate_distances = defaultdict(list)


class CENTERLINE:
    def __init__(self, row):
        self.linestring = loads(row[11])
        self.seg_id = row[9]


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
            centerline_list.append(r)
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


def extract_centerline_points():  # maps each point that is a endpoint for a centerline segment to every centerline segment it is an endpoint for, and returns a list of all the points
    points = []
    for centerline in centerline_list:
        centerline_coords = list(centerline.linestring.coords)
        for coord in centerline_coords:
            centerline_endpoints_map[coord].append(centerline.linestring)
        # print(map[coord])
        points.extend(centerline_coords)
    return list(set(points))


def near(fleetptx, fleetpty, segx1, segy1, segx2, segy2):
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
    if ((segx1 <= assocx <= segx2) & ((segy2 <= assocy <= segy1) | (segy1 <= assocy <= segy2))) | \
            ((segx2 <= assocx <= segx1) & ((segy2 <= assocy <= segy1) | (segy1 <= assocy <= segy2))):
        print('ok')
    else:
        if sqrt(pow(fleetptx - segx1, 2) + pow(fleetpty - segy1, 2)) < sqrt(pow(fleetptx - segx2, 2) + pow(fleetpty - segy2, 2)):
            assocx = segx1
            assocy = segy1
        else:
            assocx = segx2
            assocy = segy2
    dist = sqrt(pow(fleetptx - assocx, 2) + pow(fleetpty - assocy, 2))
    return assocx, assocy, dist


def euclid_dist(p1, p2):
    p1_candidates = candidate_associations[p1]
    p2_candidates = candidate_associations[p2]
    for c1 in p1_candidates:
        coord1 = [p1, c1[0], c1[1]]
        for c2 in p2_candidates:
            dist = sqrt(pow(c2[0] - coord1[1], 2) + pow(c2[1] - coord1[2], 2))
            candidate_distances[coord1].append([c2[0], c2[1], dist])


def viterbi():
    score = defaultdict(list)  # for each candidate, format is [lon, lat, time] and it stores score
    tscore = defaultdict(list)  # for each pair combination of candidates, format is [lon, lat, lon, lat] and it stores score
    pred = defaultdict(list) # for each matching, [current lon, current lat, current time] -> [previous lon, previous lat]

    print('INITIAL SCORING')
    # Scoring of candidates for first point
    for c in candidate_associations[0]:
        lon = c[0]
        lat = c[1]
        perp_dist = c[2]
        if 1 - perp_dist/160 > 0:
            score[lon, lat, 1] = 1 - perp_dist / 160
        else:
            score[lon, lat, 1] = 0
        print(perp_dist)
        print(score[lon, lat, 1])

    time = 2
    while time <= len(fleetpoint_list):
        print(time)
        print('NEW TIME')
        for current in candidate_associations[time-1]:  # time index starts at 1 while orig index starts at 0
            print('candidate', current)
            score[current[0], current[1], time] = -1
            if 1 - current[2]/160 > 0:
                emiss = 1 - current[2]/160
            else:
                emiss = 0
            for prev in candidate_associations[time-2]:
                trans_dist = sqrt(pow(current[0] - prev[0], 2) + pow(current[1] - prev[1], 2))
                print('dist', trans_dist)
                if 1 - trans_dist/1000 < 0:
                    trans_dist = 0
                else:
                    trans_dist = 1 - trans_dist/1000
                tscore[prev[0], prev[1], current[0], current[1]] = score[prev[0], prev[1], time-1] * pow(trans_dist, 2) * emiss
                print(tscore[prev[0], prev[1], current[0], current[1]])
                if tscore[prev[0], prev[1], current[0], current[1]] > score[current[0], current[1], time]:
                    score[current[0], current[1], time] = tscore[prev[0], prev[1], current[0], current[1]]
                    pred[current[0], current[1], time] = [prev[0], prev[1]]
                    print('PRED')
                    print(prev[0], prev[1])
        time += 1
    m = 0
    for s in score:
        print(s)
    final_assoc = -1
    for c in candidate_associations[len(fleetpoint_list) - 1]:
        print(c)
        sc = score[c[0], c[1], len(fleetpoint_list)]
        if sc > final_assoc:
            associations_list.clear()
            associations_list.append([c[0], c[1]])
            final_assoc = score[c[0], c[1], len(fleetpoint_list)]
    print('final assoc', associations_list[0])
    time = len(fleetpoint_list)
    current_pred = associations_list[0]
    while time >= 2:
        print(time)
        current_lon = current_pred[0]
        current_lat = current_pred[1]
        current_pred = pred[current_lon, current_lat, time]
        print(current_pred)
        associations_list.insert(0, current_pred)
        print('current list', associations_list)
        time -= 1
        print('ok?')

    print('ASSOCIATIONS')
    print(associations_list)



def associate():
    j = 0
    while j < len(fleetpoint_list):
        #if j == 103:
        fleetptx = fleetpoint_list[j].lon
        fleetpty = fleetpoint_list[j].lat
        print([fleetptx, fleetpty])
        closest_centerline_coords = tree.query([fleetptx, fleetpty], 4)[1]  # query the tree for the 4 closest points to the gps point, returns their index in centerline_points
        checked_centerlines = set()  # keep track of the segments we have already considered for this particular point
        shortestDist = sys.float_info.max
        shortestDist_seg = 0
        assocx = 0
        assocy = 0

        k = 0 #useful in print statements for keeping track of the coordinates we are checking
        for coord in closest_centerline_coords:
            print('new candidate coordinate')
            linestrings = centerline_endpoints_map[centerline_endpoints_list[coord]]  # find the list of LineStrings that each coordinate is mapped to
            for linestring in linestrings:
                print(k)
                print('current centerline:')
                print(linestring)
                if str(linestring) not in checked_centerlines:
                    print('this centerline has not been checked yet')
                    checked_centerlines.add(str(linestring))
                    l = 0
                    centerline_coords = list(linestring.coords)
                    while l < len(centerline_coords)-1: # consider each line that is a part of the current centerline LineString
                        segx1 = centerline_coords[l][0]
                        segy1 = centerline_coords[l][1]
                        segx2 = centerline_coords[l+1][0]
                        segy2 = centerline_coords[l+1][1]
                        print('Segment:')
                        print([segx1, segy1, segx2, segy2])
                        near_results = near(fleetptx, fleetpty, segx1, segy1, segx2, segy2)
                        candidate_associations[j].append(near_results)
                        l += 1
                else:
                    print('already checked segment')
                k += 1

        #print([assocx, assocy]) #print(shortestDist_seg)
        #associations_list.append([assocx, assocy])
        #if j > 0:
            #euclid_dist(j-1, j)
        j += 1
    print("VITERBI")
    viterbi()


def write_associations():
    path = csv_path(outfile)
    f = open(path, 'w+')
    header = 'id,runid,runid_sequence,vin,datetime,lon,lat,assoc_lon,assoc_lat\n'
    f.write(header)
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
print('Finished reading input files. Now indexing centerlines')

centerline_endpoints_map = defaultdict(list)
centerline_endpoints_list = extract_centerline_points()
tree = spatial.KDTree(centerline_endpoints_list)

print('Finished indexing.  Now associating')

associate()

print('Finished associating. Now writing')

write_associations()
