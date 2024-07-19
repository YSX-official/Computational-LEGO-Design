import math
import cv2 as cv
import numpy as np
import graph as gh
import bricks as bk
import shapely.geometry as geo

def counter_clockwise(contour):
    max_right,index = -math.inf,0
    for point in contour:
        if point[0] > max_right:
            max_right = point[0]
            index = contour.index(point)

    index_pre = index - 1
    index_post = (index + 1) % len(contour)
    point_pre = contour[index_pre]
    point_post = contour[index_post]
    point = contour[index]
    product = (point[0] - point_pre[0]) * (point_post[1] - point[1]) - (point[1] - point_pre[1]) * (point_post[0] - point[0])
    
    if product > 0:
        return True
    elif product < 0:
        return False
    elif point[1] - point_pre[1] > 0:
        return False
    else:
        return True

def rotate(angle):
    if angle == 90:
        return [1,0]
    elif angle == 180:
        return [0,-1]
    elif angle == 270:
        return [-1,0]
    else:
        return [0,1]

def flip(s,direction):
    if direction == 1:
        return [s[0],-s[1]]
    elif direction == 2:
        return [-s[0],-s[1]]
    elif direction == 3:
        return [-s[0],s[1]]
    else:
        return [s[0],s[1]]

def opposite(start,end,contour):
    min1,idx1 = math.inf,0
    min2,idx2 = math.inf,0

    for i in range(len(contour)):
        point = contour[i]
        if abs(point[0] - start[0]) <= 1:
            dis = (point[0] - start[0]) ** 2 + (point[1] - start[1]) ** 2
            if dis < min1:
                min1 = dis
                idx1 = i
        if abs(point[0] - end[0]) <= 1:
            dis = (point[0] - end[0]) ** 2 + (point[1] - end[1]) ** 2
            if dis < min2:
                min2 = dis
                idx2 = i

    if (idx2 - idx1 >= 0 and idx2 - idx1 < 25) or idx1 - idx2 >= 100:
        return False
    else:
        return True

def infeasible(s,l):
    if s == (1,2) or s == (3,1):
        if l < 0:
            return True
    return False

object_name = 'Apple'
file_name = 'test/' + object_name + '.png'
img = cv.imread(file_name)
width,height = 640,540
img = cv.resize(img,(width,height))
sample = list()
contour = list()

gray = cv.cvtColor(img,cv.COLOR_BGR2GRAY)
gray = cv.GaussianBlur(gray,(5,5),0)
ret,thresh = cv.threshold(gray,0,255,cv.THRESH_BINARY_INV + cv.THRESH_OTSU)
contours,hierarchy = cv.findContours(thresh,cv.RETR_TREE,cv.CHAIN_APPROX_NONE)
cnt = contours[0]

for point in cnt:
    x,y = point.ravel()
    y = height - y
    x = round(x * 0.05,2)
    y = round(y * 0.05,2)
    sample.append((x,y))

sample_polygon = geo.LinearRing(sample)
length = int(sample_polygon.length)
sample_points = [sample_polygon.interpolate(i * 0.5) for i in range(2 * length + 2)]

for point in sample_points:
    x = point.x
    y = point.y
    contour.append((x,y))

if counter_clockwise(contour):
    contour.reverse()

grid_points = set()
angles = [0,90,180,270]
directions = [0,1,2,3]
thresh = 1

for i in range(len(contour)):
    x,y = contour[i]
    xi = int(x)
    yi = int(y)
    for l in range(xi - 1,xi + 2):
        for w in range(yi - 1,yi + 2):
            if ((x - l) ** 2 + (y - w) ** 2) <= thresh:
                if l >= 0 and l <= 32 and w >= 0 and w <= 27:
                    grid_points.add((l,w))

V = set()
S = set()
S = {(0,1),(1,1),(1,2),(1,3),(2,1),(3,1)}
region = geo.Polygon(contour)

for point in grid_points:
    for (x,y) in S:
        if x == 0:
            for angle in angles:
                [l,w] = rotate(angle)
                (l1,w1) = (point[0] + l,point[1] + w)
                if opposite(point,(l1,w1),contour):
                    continue
                if (l1,w1) in grid_points:
                    subcurve = gh.projection((point,(l1,w1)),contour)
                    if subcurve.length < 3 and 1 < 3 * subcurve.length:
                        V.add((point,(l1,w1)))
        else:
            for direction in directions:
                [l,w] = flip((x,y),direction)
                (l1,w1) = (point[0] + l,point[1] + w)
                if infeasible((x,y),l) or opposite(point,(l1,w1),contour):
                    continue
                if (l1,w1) in grid_points:
                    subcurve = gh.projection((point,(l1,w1)),contour)
                    length = np.sqrt(x ** 2 + y ** 2)
                    if subcurve.length < 3 * length and length < 3 * subcurve.length:
                        V.add((point,(l1,w1)))

C = set()
M = dict()

for s in V:
    D = bk.brick_occupy(s)
    for cell in D:
        if cell in M.keys():
            M[cell].add(s)
        else:
            edges_set = set()
            edges_set.add(s)
            M.update({cell:edges_set})
            
for key in M.keys():
    for s1 in M[key]:
        for s2 in M[key]:
            if s1 == s2:
                continue
            elif bk.is_slope(s1) and bk.is_slope(s2):
                C.add((s1,s2))
            elif bk.is_slope(s1):
                if bk.penetrate(s1,s2):
                    C.add((s1,s2))
            elif bk.is_slope(s2):
                if bk.penetrate(s2,s1):
                    C.add((s1,s2))

polygon = gh.graph_optimization(V,C,contour,grid_points)
bk.lego_construct(polygon,object_name)