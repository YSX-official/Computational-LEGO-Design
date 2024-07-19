import math
import numpy as np
import gurobipy as gp
import shapely.geometry as geo

def projection(s,contour):
    start,end = geo.Point(s[0]),geo.Point(s[1])
    curve = geo.LinearRing(contour)
    nearest_start = curve.interpolate(curve.line_locate_point(start))
    nearest_end = curve.interpolate(curve.line_locate_point(end))
    start_point = list(nearest_start.coords)[0]
    end_point = list(nearest_end.coords)[0]

    idx1,idx2 = 0,0
    min1,min2 = math.inf,math.inf
    for i in range(len(contour)):
        point = contour[i]
        dis1 = (start_point[0] - point[0]) ** 2 + (start_point[1] - point[1]) ** 2
        dis2 = (end_point[0] - point[0]) ** 2 + (end_point[1] - point[1]) ** 2
        if dis1 < min1:
            min1 = dis1
            idx1 = i
        if dis2 < min2:
            min2 = dis2
            idx2 = i

    if idx2 > idx1:
        sub_contour = contour[idx1:idx2+1]
    elif idx2 == idx1:
        sub_contour = [start_point,end_point]
    else:
        sub_contour = contour[idx2:idx1+1]
        
    if len(contour) < 2 * len(sub_contour):
        id1,id2 = max(idx1,idx2),min(idx1,idx2)
        before,after = contour[id1:],contour[0:id2+1]
        sub_contour = before + after
    subcurve = geo.LineString(sub_contour)

    return subcurve

def points_sampling(s,contour):
    subcurve = projection(s,contour)
    line = geo.LineString(s)
    d = 0.05
    
    length_s = line.length
    num_samples_s = int(length_s * 20)
    Ps = [line.interpolate(i * d) for i in range(num_samples_s + 2)]
    Ds = [subcurve.distance(p) for p in Ps]

    length_c = subcurve.length
    num_samples_c = int(length_c * 20)
    Pc = [subcurve.interpolate(i * d) for i in range(num_samples_c + 2)]
    Dc = [line.distance(c) for c in Pc]

    return Ds,Dc

def distance_deviation(Ds,Dc):
    Es = sum(Ds) / len(Ds)
    Ec = sum(Dc) / len(Dc)
    return Es + Ec
    
def distance_variation(Ds,Dc):
    Ds = np.array(Ds)
    Dc = np.array(Dc)
    Vs = np.std(Ds)
    Vc = np.std(Dc)
    return Vs + Vc

def decomposition(x1,y1,x2,y2):
    length = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    vector = ((x2 - x1) / length,(y2 - y1) / length)
    v = set()
    
    for i in range(int(length)):
        x = x1 + i * vector[0]
        y = y1 + i * vector[1]
        s = x + vector[0]
        t = y + vector[1]
        v.add(((x,y),(s,t)))

    return v

def graph_optimization(V,C,contour,grid_points):
    curve = geo.LinearRing(contour)
    Lc = curve.length
    coefficients = dict()
    wd,wv = 1,3

    for s in V:
        Ds,Dc = points_sampling(s,contour)
        Ldx = distance_deviation(Ds,Dc) * len(Ds) / Lc
        v = decomposition(s[0][0],s[0][1],s[1][0],s[1][1])
        Lvx = 0
        for l in v:
            Dl,Dc = points_sampling(l,contour)
            Lvx += distance_variation(Dl,Dc) * 20
        Lvx /= Lc
        coefficients.update({s:wd * Ldx + wv * Lvx})

    model = gp.Model()
    var = dict()

    for s in V:
        x = model.addVar(vtype = gp.GRB.BINARY)
        var.update({s:x})
    
    objective_function = 0
    for s in V:
        objective_function += var[s] * coefficients[s]
        
    model.setObjective(objective_function,gp.GRB.MINIMIZE)

    for s1 in V:
        for s2 in V:
            if s1 == s2:
                continue
            elif s1[0] == s2[1] and s1[1] == s2[0]:
                model.addConstr(var[s1] * var[s2] == 0)

    for (s1,s2) in C:
        model.addConstr(var[s1] * var[s2] == 0)
        
    model.addConstr(gp.quicksum(var.values()) >= 1)

    for point in grid_points:
        degree = gp.LinExpr()
        for s in V:
            if s[1] == point:
                degree += var[s]
            elif s[0] == point:
                degree -= var[s]
        model.addConstr(degree == 0)
        
    model.optimize()
    polygon = list()

    for s in V:
        if var[s].x == 1:
            polygon.append(s)
            
    return polygon