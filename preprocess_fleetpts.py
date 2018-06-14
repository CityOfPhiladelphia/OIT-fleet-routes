import os
import time
import csv
import json
import datetime
from pyproj import Proj, transform
from math import *
import math

__author__ = 'tom.swanson'
cwd = os.path.dirname(__file__) # the current directory
p1 = Proj(init='epsg:{}'.format(4326), preserve_units=True)
p2 = Proj(init='epsg:{}'.format(2272), preserve_units=True)

networfleet_list = []

def csv_path(file_name):
    return os.path.join(cwd, file_name + '.csv')


def lineDirAngle(lon1, lat1, lon2, lat2):
    # var y = Math.sin(dLon) * Math.cos(lat2);
    # var x = Math.cos(lat1)*Math.sin(lat2) -
    #        Math.sin(lat1)*Math.cos(lat2)*Math.cos(dLon);
    # var brng = Math.atan2(y, x).toDeg();
    return atan2(lon2 - lon1, lat2 - lat1) * 180 / math.pi
    deltax = lon2 - lon1
    deltay = lat2 - lat1
    angle_rad = atan2(deltay, deltax)
    angle_deg = angle_rad * 180.0 / pi
    # print "The angle is %.5f radians (%.5f degrees)." % (angle_rad,angle_deg)
    return angle_deg


# 7
def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    miles = 3961 * c
    # meters = (6367 * c) * 1000
    return miles


start = time.time()

# 3
class NETWORKFLEET:
    def __init__(self, row):
        # self.field1 = row[0].strip()
        # TODO error handle empty row
        self.vin = row[1].strip()
        self.fleetid = row[2].strip()
        self.messageid = row[3].strip()
        # self.datetime = row[4].strip()
        self.timeint = int(row[5].strip())
        self.status = row[6].strip()
        self.statusid = row[7].strip()
        ack = row[8].strip()
        self.coords = row[9].strip()
        temp = self.coords.split(',')
        x,y = transform(p1, p2, temp[0],temp[1])
        self.lon = x
        self.lat = y
        self.json_items = json.loads(ack)
        self.odometer = self.json_items['Odometer']
        self.ignition = str(self.json_items['Ignition'])
        # self.timeUTF = self.json_items['FixTimeUTF']
        dt_tmp = datetime.datetime.utcfromtimestamp(float(self.timeint) - 18000.0).strftime('%Y-%m-%dT%H:%M:%S')
        self.dt = datetime.datetime.strptime(dt_tmp, '%Y-%m-%dT%H:%M:%S')
        self.dayofweek = self.dt.weekday()
        self.year = self.dt.year
        self.month = self.dt.month
        self.day = self.dt.day
        self.hour = self.dt.hour
        self.minute = self.dt.minute
        self.ageMiles = self.json_items['AgeInMiles']
        speed = self.json_items['Speed']
        self.speedInst = speed[1]['Value']
        self.speedAvg = speed[0]['Value']
        self.speedMax = speed[2]['Value']
        if 'Heading' in self.json_items:
            self.heading = self.json_items['Heading']
        else:
            self.heading = None
        self.sort = "{}_{}".format(self.vin, dt_tmp)

# 2
def read_networfleet(infile):
    path = csv_path(infile)
    # We define f to be the file and make a reader for it
    f = open(path, 'r')
    i = 0
    try:
        reader = csv.reader(f)
        for row in reader:
            ''' Commented code causes reader to skip the first line of the input '''
            # if i == 0:
            #     i += 1
            #     continue
            # Create a NETWORKFLEET object from each row
            r = NETWORKFLEET(row)
            networfleet_list.append(r)
            i += 1
    except IOError:
        print('Error opening ' + path, sys.exc_info()[0])
    f.close()
    return


# 6
def write_networkfleet(newlist, outfile):
    i = 0
    path = csv_path(outfile)
    # now our file is the output file
    f = open(path, 'w+')
    header = 'id,runid,runid_sequence,vin,datetime,time_delta,year,month,day,hour,min,dayofweek,messageid,ignition,odometer,odometer_dist,xy_dist,speedInst,speedAvg,speedMax,heading,headingcalculated,lon,lat\n'
    f.write(header)
    runid = 0
    runid_sequence = 0
    # time_delta = 0
    # odometer_dist = 0.0
    # xy_dist = 0.0
    # heading = 0.0
    idle = False
    previous = {}
    for row in newlist:
        if i == 0:
            xy_dist = 0
            heading = 'None'
            odometer_dist = 0
            time_delta = 0
            runid_sequence = 0
            runid = 0
        # we want to see if it's idle only if it's the same truck so check the id
        elif row.vin == previous.vin:
            # So we start by seeing if the truck is currently off.  If so, it is the next data point in the current run's sequence, and idle is false
            if row.ignition == 'Off':
                # calculate the distance that has been traveled since the previous point
                xy_dist = haversine(float(row.lon), float(row.lat), float(previous.lon), float(previous.lat))
                if idle == False:
                    # calculate the angle degree of change that has occurred since the previous point
                    heading = lineDirAngle(float(row.lon), float(row.lat), float(previous.lon), float(previous.lat))
                else:
                    heading = 'None' # So the heading is only a thing if it's not idle
                odometer_dist = row.odometer - previous.odometer
                time_delta = (row.timeint - previous.timeint)
                runid_sequence += 1
                idle = False   
            # Truck is currently on but at last tracking it was off.  So we start a new run and this data point is the 0th in the sequence.
            elif previous.ignition == 'Off':
                xy_dist = 0
                heading = 'None'
                odometer_dist = 0
                time_delta = 0
                runid_sequence = 0
                runid += 1
            # Truck has not been on for the past two and has been marked idle
            elif idle == True:
                # I think this continue statement indicates that if it is idle and the odometer hasn't changed, the point is not written to the output file.
                if row.odometer == previous.odometer:
                    xy_dist = 0
                    heading = 'None'
                    odometer_dist = 0
                    time_delta = (row.timeint - previous.timeint)
                    continue
                else:
                    xy_dist = haversine(float(row.lon), float(row.lat), float(previous.lon), float(previous.lat))
                    heading = lineDirAngle(float(row.lon), float(row.lat), float(previous.lon), float(previous.lat))
                    odometer_dist = row.odometer - previous.odometer
                    time_delta = (row.timeint - previous.timeint)
                    runid_sequence += 1
                    idle = False
            # Truck has not been on for the past two and has not been marked idle
            elif idle == False:
                if row.odometer == previous.odometer:
                    idle = True
                    xy_dist = 0
                    heading = 'None'
                    odometer_dist = 0
                    time_delta = (row.timeint - previous.timeint)
                    runid_sequence += 1
                else:
                    xy_dist = haversine(float(row.lon), float(row.lat), float(previous.lon), float(previous.lat))
                    heading = lineDirAngle(float(row.lon), float(row.lat), float(previous.lon), float(previous.lat))
                    odometer_dist = row.odometer - previous.odometer
                    time_delta = (row.timeint - previous.timeint)
                    runid_sequence += 1
        # must be a new truck from the previous row so we need a new runid (and won't be false)
        else:
            xy_dist = 0
            heading = 'None'
            odometer_dist = 0
            time_delta = 0
            runid_sequence = 0
            runid += 1
            idle = False
        # when i=0, it will be false.  After doing whatever happens above (I would assume calculating if it's idle), you set what the csv row data should be.  
        # The difference is that if it's idle, you don't include ignition
        if idle == False:
            s = '%i,%i,%i,%s,%s,%s,%i,%i,%i,%i,%i,%i,%s,%s,%f,%f,%f,%i,%i,%i,%s,%s,%s,%s\n' % (
                i,
                runid,
                runid_sequence,
                row.vin,
                str(row.dt),
                time_delta,
                row.year,
                row.month,
                row.day,
                row.hour,
                row.minute,
                row.dayofweek,
                row.messageid,
                row.ignition,
                row.odometer,
                odometer_dist,
                xy_dist,
                row.speedInst,
                row.speedAvg,
                row.speedMax,
                row.heading,
                heading,
                row.lon,
                row.lat
            )
        else:
            s = '%i,%i,%i,%s,%s,%s,%i,%i,%i,%i,%i,%i,%s,%s,%f,%f,%f,%i,%i,%i,%s,%s,%s,%s\n' % (
                i,
                runid,
                runid_sequence,
                row.vin,
                str(row.dt),
                time_delta,
                row.year,
                row.month,
                row.day,
                row.hour,
                row.minute,
                row.dayofweek,
                row.messageid,
                'Idle',
                row.odometer,
                odometer_dist,
                xy_dist,
                row.speedInst,
                row.speedAvg,
                row.speedMax,
                row.heading,
                heading,
                row.lon,
                row.lat
            )
        f.write(s)
        i += 1
        # keep track of the previous data
        previous = row
    f.close()


def record_tracking(i):
    global intermediate, i_time, rec_per_sec
    i += 1
    if i % 100 == 0:
        intermediate = time.time()
        i_time = intermediate - start
        rec_per_sec = i / i_time
        # print('%i %f' % (i, rec_per_sec))
    return i


def convert_size(size_bytes):
    if size_bytes == 0:
        return '0'
    s = round((size_bytes / 1048576.0), 2)
    s = str(s)
    return s


def run(inf, outf):
    infile = inf
    outfile = outf

    # 1. STARTS HERE -- creates objects for each row of data in the csv and puts them in networfleet_list
    print('reading')
    read_networfleet(infile)


    # 4. Sort the networfleet_list data by the class's sort property, which takes in vin, dt_tmp
    print('sorting')
    networkfleet_list2 = sorted(networfleet_list, key=lambda nf: nf.sort)

    # 5. Make some calculations/cleaning for each data point and then write the finished product to the new file
    print('writing')
    write_networkfleet(networkfleet_list2, outfile)