from scipy import spatial
from collections import defaultdict
import shapely.geos
import os
import csv
import sys
from math import sqrt, pow
from shapely.geometry import Point, LineString
from shapely.wkt import loads


centerlines_file = 'centerline_shape'
fleetpoints_file = 'test_fleet'
outfile = '3_20_18'
cwd = os.path.dirname(__file__)

centerline_list = []
fleetpoint_list = []


class CENTERLINE:
    def __init__(self, row):
        self.linestring = loads(row[11])
        self.seg_id = row[9]


class FLEETPOINT:
    def __init__(self, row):
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


read_fleetpoints()
read_centerlines()
print('Finished reading input files')


def extract_centerline_points(): # maps each point that is a endpoint for a centerline segment to every centerline segment it is an endpoint for, and returns a list of all the points
    points = []
    for centerline in centerline_list:
        centerline_coords = list(centerline.linestring.coords)
        for coord in centerline_coords:
            centerline_endpoints_map[coord].append(centerline.linestring)
        #print(map[coord])
        points.extend(centerline_coords)
    return list(set(points))


centerline_endpoints_map = defaultdict(list)
centerline_endpoints_list = extract_centerline_points()
tree = spatial.KDTree(centerline_endpoints_list)


def perpendicular(px, py, sx1, sy1, sx2, sy2): # returns the closest point on the input line to the input point
    segSlope = (sy2-sy1)/(sx2-sx1)
    segB = sy2 - segSlope*sx2
    pSlope = -1/segSlope
    pB = py - pSlope*px
    assocx = (segB-pB)/(pSlope-segSlope)
    assocy = pSlope*assocx + pB
    return assocx, assocy


associations_list = []


def associate():
    path = csv_path(outfile)
    f = open(path, 'w+')
    header = 'id,runid,runid_sequence,vin,datetime,lon,lat,assoc_lon,assoc_lat\n'
    f.write(header)

    associations_list.append([fleetpoint_list[0].lon, fleetpoint_list[0].lat]) # hard-coded first point, needs to be fixed
    j = 1
    while j < len(fleetpoint_list):
        previous_lon = associations_list[j-1][0]
        previous_lat = associations_list[j - 1][1]
        print([previous_lon, previous_lat])
        print([fleetpoint_list[j].lon, fleetpoint_list[j].lat])
        closest_centerline_coords = tree.query([fleetpoint_list[j].lon, fleetpoint_list[j].lat], 4)[1]  # query the tree for the 4 closest points to the gps point, returns their index in centerline_points
        centerlines = set()  # keep track of the segments we have already considered for this particular point
        shortestDist = 100
        shortestDist_seg = 0
        shortestDist_nonSlope = 100
        shortestDist_seg_nonSlope = 0
        assocx = 0
        assocy = 0
        if previous_lon == fleetpoint_list[j].lon: # either if the location has not changed at all or if the slope is vertical
            slope_from_prev = 100
        else:
            slope_from_prev = abs((previous_lat-fleetpoint_list[j].lat)/(previous_lon-fleetpoint_list[j].lon))
        dist_from_prev = sqrt(pow(previous_lat-fleetpoint_list[j].lat, 2) + pow(previous_lon-fleetpoint_list[j].lon, 2))
        k = 0 #useful in print statements for keeping track of the coordinates we are checking
        for coord in closest_centerline_coords:
            print('new candidate coordinate')
            segs = centerline_endpoints_map[centerline_endpoints_list[coord]]  # find the list of LineStrings that each coordinate is mapped to
            for seg in segs:
                print(k)
                print('current centerline:')
                print(seg)
                if str(seg) not in centerlines:
                    # print('considering new segment')
                    centerlines.add(str(seg))
                    l = 0
                    seg_coords = list(seg.coords)
                    px = fleetpoint_list[j].lon
                    py = fleetpoint_list[j].lat
                    while l < len(seg_coords)-1: # consider each line that is a part of the current centerline LineString
                        sx1 = seg_coords[l][0]
                        sy1 = seg_coords[l][1]
                        sx2 = seg_coords[l+1][0]
                        sy2 = seg_coords[l+1][1]
                        print('Segment:')
                        print([sx1, sy1, sx2, sy2])
                        if j > 0:
                            seg_slope = abs((sy2-sy1)/(sx2-sx1))
                        else:
                            seg_slope = 100
                        print(seg_slope)
                        assoc = ()
                        if abs(slope_from_prev - seg_slope) < 0.1: # and dist_from_prev < .001
                            # should really check the non-perpendicular distance every time as well
                            if sqrt(pow(px - sx1, 2) + pow(py - sy1, 2)) < sqrt(pow(px - sx2, 2) + pow(py - sy2, 2)):
                                dist = sqrt(pow(px - sx1, 2) + pow(py - sy1, 2))
                                assoc = (sx1, sy1)
                            else:
                                dist = sqrt(pow(px - sx2, 2) + pow(py - sy2, 2))
                                assoc = (sx2, sy2)
                            if (sx1 < px < sx2) | (sx2 < px < sx1) | (sy1 < py < sy2) | (sy2 < py < sy1):
                                assoc = perpendicular(px, py, sx1, sy1, sx2, sy2)
                                print(assoc)
                                dist2 = sqrt(pow(px-assoc[0],2) + pow(py-assoc[1],2))
                                if dist2 < dist:
                                    dist = dist2
                            if 0 < dist < shortestDist: # dist is 0 when the point is out of bounds of the line
                                shortestDist = dist
                                shortestDist_seg = [sx1, sy1, sx2, sy2]
                                assocx = assoc[0]
                                assocy = assoc[1]
                                print('new shortest dist:')
                                print(shortestDist)
                        else: # if the points don't have a similar slope, then consider its perpendicular separately
                            if sqrt(pow(px - sx1, 2) + pow(py - sy1, 2)) < sqrt(pow(px - sx2, 2) + pow(py - sy2, 2)):
                                dist = sqrt(pow(px - sx1, 2) + pow(py - sy1, 2))
                                assoc = (sx1, sy1)
                            else:
                                dist = sqrt(pow(px - sx2, 2) + pow(py - sy2, 2))
                                assoc = (sx2, sy2)
                            if (sx1 < px < sx2) | (sx2 < px < sx1) | (sy1 < py < sy2) | (sy2 < py < sy1):
                                assoc = perpendicular(px, py, sx1, sy1, sx2, sy2)
                                print(assoc)
                                dist2 = sqrt(pow(px-assoc[0],2) + pow(py-assoc[1],2))
                                if dist2 < dist:
                                    dist = dist2
                            if 0 < dist < shortestDist_nonSlope: # dist is 0 when the point is out of bounds of the line
                                shortestDist_nonSlope = dist
                                shortestDist_seg_nonSlope = [sx1, sy1, sx2, sy2]
                                assocx = assoc[0]
                                assocy = assoc[1]
                                print('new shortest dist:(non slope)')
                                print(shortestDist_nonSlope)
                        l += 1
                else:
                    print('already checked segment')
                k += 1
        chosen = shortestDist
        if chosen > 100: # there wasn't a centerline with a good slope
            chosen = shortestDist_nonSlope
        print('chosen shortest dist:')
        print([assocx, assocy]) #print(shortestDist_seg)
        associations_list.append([assocx, assocy])

        s = '%i,%i,%i,%s,%s,%f,%f,%f,%f\n' % (
            j,
            fleetpoint_list[j].runid,
            fleetpoint_list[j].runid_sequence,
            fleetpoint_list[j].vin,
            str(fleetpoint_list[j].datetime),
            fleetpoint_list[j].lon,
            fleetpoint_list[j].lat,
            assocx,
            assocy
        )
        f.write(s)
        j += 1
    f.close()


associate()


