from collections import defaultdict
import os
import csv
import sys
from shapely.wkt import loads


centerlines_file = 'centerlines_joined'
outfile = '5_17_18_topology'
cwd = os.path.dirname(__file__)

centerline_list = []
centerline_endpoints_map = defaultdict(list)  # a mapping of coordinates to segment point ids
seg_points = defaultdict(list)  # a mapping of segment point ids to the SEGPOINT objects


class CENTERLINE:
    def __init__(self, row):
        self.linestring = loads(row[24])
        self.seg_id = row[17]
        self.dir = row[18]


class SEGPOINT:
    def __init__(self, id, coords):
        self.segpoint_id = id
        self.lon = coords[0]
        self.lat = coords[1]
        self.in_segs = []
        self.goes_to = set([])


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


def extract_centerline_points():  # maps each point that is a endpoint for a centerline segment to every centerline segment it is an endpoint for, and returns a list of all the points
    points = []
    id = 1
    for centerline in centerline_list:
        centerline_coords = list(centerline.linestring.coords)
        for coord in centerline_coords:  # first, we make an object for each point in the centerline
            if coord not in centerline_endpoints_map:
                seg_points[id] = SEGPOINT(id, coord)
                centerline_endpoints_map[coord] = id
            segpoint = seg_points[centerline_endpoints_map[coord]]
            segpoint.in_segs.extend(centerline.seg_id)
        if centerline.dir == 'TF' or centerline.dir == 'B':
            i = 1
            while i < len(centerline_coords):  # if it is To-From, we assign goes_to accordingly
                segpoint = seg_points[centerline_endpoints_map[centerline_coords[i]]]
                prev_segpoint = seg_points[centerline_endpoints_map[centerline_coords[i-1]]]
                if prev_segpoint.segpoint_id not in segpoint.goes_to:
                    segpoint.goes_to.add(prev_segpoint.segpoint_id)
                i += 1
        if centerline.dir == 'FT' or centerline.dir == 'B':
            i = 0
            while i < len(centerline_coords) - 1:
                segpoint = seg_points[centerline_endpoints_map[centerline_coords[i]]]
                next_segpoint = seg_points[centerline_endpoints_map[centerline_coords[i+1]]]
                if next_segpoint not in segpoint.goes_to:
                    segpoint.goes_to.add(next_segpoint.segpoint_id)
                i += 1
        id += 1
        points.extend(centerline_coords)


def write_seg_topology():
    path = csv_path(outfile)
    f = open(path, 'w+')
    header = 'id,lon,lat,in_segs,goes_to\n'
    f.write(header)
    for pointid in seg_points:
        segpoint = seg_points[pointid]
        s = '%i,%f,%f,%s,%s\n' % (
            segpoint.segpoint_id,
            segpoint.lon,
            segpoint.lat,
            str(segpoint.goes_to),
            str(segpoint.in_segs)
        )
        f.write(s)
    f.close()


read_centerlines()
print('Finished reading input files. Now indexing centerlines')

extract_centerline_points()
print('Finished extracting. Now writing')

write_seg_topology()
