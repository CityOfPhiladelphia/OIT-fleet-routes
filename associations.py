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
fleetpoints_file = 'XY06_01_originals_feet'
outfile = 'new_near_06_01'
cwd = os.path.dirname(__file__)

centerline_list = []
fleetpoint_list = []
associations_list = []


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
        self.lon = float(row[5])
        self.lat = float(row[6])
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


def extract_centerline_points():
    points = []
    for centerline in centerline_list:
        centerline_coords = list(centerline.linestring.coords)
        l = 0
        centerline_endpoints_map[centerline_coords[l][0],centerline_coords[l][1]].append(centerline.linestring)
        while l < len(centerline_coords) - 1:  # consider each line that is a part of the current centerline LineString
            segx1 = centerline_coords[l][0]
            segy1 = centerline_coords[l][1]
            segx2 = centerline_coords[l + 1][0]
            segy2 = centerline_coords[l + 1][1]
            if sqrt(pow(segx2 - segx1, 2) + pow(segy2 - segy1, 2)) > 500:
                segx3 = min(segx1, segx2) + abs(segx1-segx2)/2
                segy3 = min(segy1, segy2) + abs(segy1 - segy2) / 2
                centerline_endpoints_map[segx3, segy3].append(centerline.linestring)
                points.extend([(segx3, segy3)])
            l += 1
            centerline_endpoints_map[segx2, segy2].append(centerline.linestring)
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


def associate():
    j = 0
    while j < len(fleetpoint_list):
        #if j == 103:
        fleetptx = fleetpoint_list[j].lon
        fleetpty = fleetpoint_list[j].lat
        print([fleetptx, fleetpty])
        closest_centerline_coords = tree.query([fleetptx, fleetpty], 6)[1]  # query the tree for the 4 closest points to the gps point, returns their index in centerline_points
        checked_centerlines = set()  # keep track of the segments we have already considered for this particular point
        shortestDist = sys.float_info.max
        shortestDist_seg = 0
        assocx = 0
        assocy = 0

        k = 0 #useful in print statements for keeping track of the coordinates we are checking
        for coord in closest_centerline_coords:
            if j == 41:
                print('coord',centerline_endpoints_list[coord])
            #print('new candidate coordinate')
            linestrings = centerline_endpoints_map[centerline_endpoints_list[coord]]  # find the list of LineStrings that each coordinate is mapped to
            #print(linestrings)
            for linestring in linestrings:
                #print(k)
            #    print('current centerline:')
            #    print(linestring)
                if str(linestring) not in checked_centerlines:
            #        print('this centerline has not been checked yet')
                    checked_centerlines.add(str(linestring))
                    l = 0
                    centerline_coords = list(linestring.coords)
                    while l < len(centerline_coords)-1: # consider each line that is a part of the current centerline LineString
                        segx1 = centerline_coords[l][0]
                        segy1 = centerline_coords[l][1]
                        segx2 = centerline_coords[l+1][0]
                        segy2 = centerline_coords[l+1][1]
             #           print('Segment:')
                        #print([segx1, segy1, segx2, segy2])

                        near_results = near(fleetptx, fleetpty, segx1, segy1, segx2, segy2)

                        if near_results[2] < shortestDist: # decide whether to use the calculated Near value
             #               print('new shortest dist:')
                            #print(near_results)
                            shortestDist = near_results[2]
                            shortestDist_seg = [segx1, segy1, segx2, segy2]
                            assocx = near_results[0]
                            assocy = near_results[1]
                        else:
                            print('not shorter')
                        l += 1
                else:
                    print('already checked segment')
                k += 1

       #print('chosen shortest dist:')
        #print([assocx, assocy]) #print(shortestDist_seg)
        associations_list.append([assocx, assocy])
        j += 1


def write_associations():
    path = csv_path(outfile)
    f = open(path, 'w+')
    header = 'id,runid,runid_sequence,vin,datetime,lon,lat,assoc_lon,assoc_lat\n'
    f.write(header)
    j=0
    while j < len(fleetpoint_list):
        s = '%i,%i,%i,%s,%s,%f,%f,%f,%f\n' % (
            j,
            fleetpoint_list[j].runid,
            fleetpoint_list[j].runid_sequence,
            fleetpoint_list[j].vin,
            str(fleetpoint_list[j].datetime),
            fleetpoint_list[j].lon,
            fleetpoint_list[j].lat,
            associations_list[j][0],
            associations_list[j][1]
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
