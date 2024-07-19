import math
import pymunk as pk
import pygame as pg
import gurobipy as gp
import shapely.geometry as geo

class Brick:
    def __init__(self,color,h,w,weight,dat_code,slope):
        self.color = color
        self.h = h
        self.w = w
        self.weight = weight
        self.dat_code = dat_code
        self.slope = slope 

    def occupy(self,s):
        (x,y),(x1,y1) = s
        dx,dy = x1 - x,y1 - y
        z,t = min(x,x1),min(y,y1)
        D = set()
        if self.slope[1] > 0:
            if dx > 0 and dy > 0:
                for i in range(z + 1,z + self.w + 1):
                    for j in range(t + 1,t + self.h + 1):
                        D.add((i,j))
            else:
                for i in range(z,z + self.w):
                    for j in range(t + 1,t + self.h + 1):
                        D.add((i,j))
        else:
            if dx < 0 and dy < 0:
                for i in range(z,z + self.w):
                    for j in range(t + 1,t + self.h + 1):
                        D.add((i,j))
            else:
                for i in range(z + 1,z + self.w + 1):
                    for j in range(t + 1,t + self.h + 1):
                        D.add((i,j))
        return D

# Brick
brick1 = Brick('9',1,1,0.00043,'3005.DAT',[])
brick2 = Brick('9',1,2,0.00081,'3004.DAT',[])
brick3 = Brick('9',1,3,0.00118,'3622.DAT',[])
brick4 = Brick('9',1,4,0.00157,'3010.DAT',[])
brick5 = Brick('9',1,6,0.00228,'3009.DAT',[])
brick6 = Brick('9',1,8,0.00303,'3008.DAT',[])
brick7 = Brick('9',1,10,0.00380,'6111.DAT',[])
brick8 = Brick('9',1,16,0.00620,'2465.DAT',[]) 

# Slope
brick9 = Brick('9',1,2,0.00069,'3040.DAT',[1,1])
brick10 = Brick('9',2,2,0.00116,'60481.DAT',[1,2])
brick11 = Brick('9',3,2,0.00168,'4460b.DAT',[1,3])
brick12 = Brick('9',1,3,0.00098,'4286.DAT',[2,1])
brick13 = Brick('9',1,4,0.00110,'60477.DAT',[3,1])
brick14 = Brick('9',1,2,0.00066,'3665.DAT',[1,-1])
brick15 = Brick('9',1,3,0.00094,'4287.DAT',[2,-1])
brick16 = Brick('9',3,2,0.00188,'2449.DAT',[1,-3])

def is_slope(s):
    if s[0][0] == s[1][0] or s[0][1] == s[1][1]:
        return False
    else:
        return True

def brick_associated(s):
    (x,y),(x1,y1) = s
    dx,dy = x1 - x,y1 - y
    if abs(dy) == 2:
        return brick10
    elif abs(dx) == 3:
        return brick13
    elif abs(dx) == 2:
        if dx > 0:
            return brick12
        else:
            return brick15
    elif abs(dy) == 1:
        if dx > 0:
            return brick9
        else:
            return brick14
    else:
        if dx > 0:
            return brick11
        else:
            return brick16

def brick_occupy(s):
    (x,y),(x1,y1) = s
    dx,dy = x1 - x,y1 - y
    D = set()
    if dx == 0:
        D.add((x + 1,max(y,y1)))
        D.add((x,max(y,y1)))
    elif dy == 0:
        D.add((max(x,x1),y))
        D.add((max(x,x1),y + 1))
    else:
        brick = brick_associated(s)
        D = brick.occupy(s)
    return D

def penetrate(s1,s2):
    D1 = brick_occupy(s1)
    D2 = brick_occupy(s2)
    if D2.issubset(D1):
        return True
    else:
        return False

def sort(s):
    return min(s[0][1],s[1][1]),s[0][0]

def intersect(s1,s2):
    x11,x12 = s1[0][0],s1[1][0]
    x21,x22 = s2[0][0],s2[1][0]
    if x11 <= x21 and x21 < x12:
        return True
    elif x21 <= x11 and x11 < x22:
        return True
    return False

def get_length(s):
    x1,x2 = s[0][0],s[1][0]
    return x2 - x1

def get_intersection(s1,s2):
    x11,x12 = s1[0][0],s1[1][0]
    x21,x22 = s2[0][0],s2[1][0]
    l1,l2 = get_length(s1),get_length(s2)
    start = min(x11,x21)
    end = max(x12,x22)
    return l1 + l2 - (end - start)

def pre_optimize(G,C,grid_cnt):
    objective_function = gp.LinExpr()
    model = gp.Model()
    var = dict()
    M = set()

    for s in G:
        x = model.addVar(vtype = gp.GRB.BINARY)
        var.update({s:x})
    
    for s in G:
        if get_length(s) == 1:
            objective_function += var[s]
    model.setObjective(objective_function,gp.GRB.MINIMIZE)

    for (s1,s2) in C:
        model.addConstr(var[s1] * var[s2] == 0)

    degree = gp.LinExpr()
    for s in G:
        dx = s[1][0] - s[0][0]
        degree += var[s] * dx
    model.addConstr(degree == grid_cnt)

    model.optimize()

    return model.ObjVal

def cross_detect(s1,s2):
    x11,x12 = s1[0][0],s1[1][0]
    x21,x22 = s2[0][0],s2[1][0]
    if x11 == x21 or x12 == x22:
        return True
    return False

def brick_select(s):
    x1,x2 = s[0][0],s[1][0]
    dx = x2 - x1
    if dx == 1:
        return brick1
    elif dx == 2:
        return brick2
    elif dx == 3:
        return brick3
    elif dx == 4:
        return brick4
    elif dx == 6:
        return brick5
    elif dx == 8:
        return brick6
    elif dx == 10:
        return brick7
    else:
        return brick8
        
def lego_construct(polygon,object_name):
    file_content = ''
    file_name = object_name + '.ldr'
    lego_connect = set()
    lego_slope = list()
    slope_set = set()
    N1,N2,min_height = 0,0,math.inf

    for s in polygon.copy():
        if is_slope(s):
            (x,y),(x1,y1) = s
            dx,dy = x1 - x,y1 - y
            brick = brick_associated(s)
            lego_slope.append([brick,s])
            N2 += 1
            if min_height > min(y,y1):
                min_height = min(y,y1) + 1
            if brick.slope[1] > 0:
                if dx > 0 and dy > 0:
                    z = (s[0][0] + brick.w) * 20
                    y = -(s[0][1] + brick.h) * 24
                    file_content += '0 Step\n' + '1 '
                    file_content += brick.color + ' '
                    file_content += str(0) + ' ' + str(y) + ' ' + str(z)
                    file_content += ' 1 0 0 0 1 0 0 0 1 '
                    file_content += brick.dat_code + '\n'
                    polygon.remove(s)
                    polygon.append(((x1 + 1,y1),(x1,y1)))
                    for i in range(x,x + brick.w):
                        polygon.append(((i,s[0][1]),(i + 1,s[0][1])))
                    for j in range(s[0][1],s[0][1] + brick.h):
                        polygon.append(((x1 + 1,j),(x1 + 1,j + 1)))

                    if brick.h == 1:
                        slope_set.add(((x,y1),(x + brick.w,y1)))
                    else:
                        slope_set.add(((x,s[0][1] + 1),(x + brick.w,s[0][1] + 1)))
                        slope_set.add(((x + 1,y1),(x + brick.w,y1)))

                elif dx > 0 and dy < 0:
                    z = s[0][0] * 20
                    y = -(s[1][1] + brick.h) * 24
                    file_content += '0 Step\n' + '1 '
                    file_content += brick.color + ' '
                    file_content += str(0) + ' ' + str(y) + ' ' + str(z)
                    file_content += ' -1 0 0 0 1 0 0 0 -1 '
                    file_content += brick.dat_code + '\n'
                    polygon.remove(s)
                    polygon.append(((x,s[0][1]),(x - 1,s[0][1])))
                    for i in range(x - 1,x1):
                        polygon.append(((i,y1),(i + 1,y1)))
                    for j in range(y1,s[0][1]):
                        polygon.append(((x - 1,j + 1),(x - 1,j)))

                    if brick.h == 1:
                        slope_set.add(((x - 1,s[0][1]),(x1,s[0][1])))
                    else:
                        slope_set.add(((x - 1,s[0][1]),(x,s[0][1])))
                        slope_set.add(((x - 1,y1 + 1),(x1,y1 + 1)))

            else:
                if dx < 0 and dy > 0:
                    z = (s[0][0] + 1) * 20
                    y = -(s[0][1] + brick.h) * 24
                    file_content += '0 Step\n' + '1 '
                    file_content += brick.color + ' '
                    file_content += str(0) + ' ' + str(y) + ' ' + str(z)
                    file_content += ' 1 0 0 0 1 0 0 0 1 '
                    file_content += brick.dat_code + '\n'
                    polygon.remove(s)
                    polygon.append(((x,s[0][1]),(x + 1,s[0][1])))
                    for i in range(x1,x + 1):
                        polygon.append(((i + 1,y1),(i,y1)))
                    for j in range(s[0][1],y1):
                        polygon.append(((x + 1,j),(x + 1,j + 1)))

                    if brick.h == 1:
                        slope_set.add(((x1,y1),(x1 + brick.w,y1)))
                    else:
                        slope_set.add(((x,s[0][1] + 1),(x + 1,s[0][1] + 1)))
                        slope_set.add(((x1,y1),(x1 + brick.w,y1)))

                elif dx < 0 and dy < 0:
                    z = (s[0][0] - brick.slope[0]) * 20
                    y = -(s[1][1] + brick.h) * 24
                    file_content += '0 Step\n' + '1 '
                    file_content += brick.color + ' '
                    file_content += str(0) + ' ' + str(y) + ' ' + str(z)
                    file_content += ' -1 0 0 0 1 0 0 0 -1 '
                    file_content += brick.dat_code + '\n'
                    polygon.remove(s)
                    polygon.append(((x1 - 1,y1),(x1,y1)))
                    for i in range(x1 - 1,x):
                        polygon.append(((i + 1,s[0][1]),(i,s[0][1])))
                    for j in range(y1,s[0][1]):
                        polygon.append(((x1 - 1,j + 1),(x1 - 1,j)))

                    if brick.h == 1:
                        slope_set.add(((x1 - 1,s[0][1]),(x,s[0][1])))
                    else:
                        slope_set.add(((x1 - 1,s[0][1]),(x,s[0][1])))
                        slope_set.add(((x1 - 1,y1 + 1),(x1,y1 + 1)))

    for s1 in polygon.copy():
        for s2 in polygon.copy():
            if s1 == s2:
                continue
            elif s1 in polygon and s2 in polygon and s1[0] == s2[1] and s1[1] == s2[0]:
                polygon.remove(s1)
                polygon.remove(s2)

    for s in polygon.copy():
        if s[0][1] == s[1][1]:
            polygon.remove(s)

    E = set()
    G = set()
    P = set()
    C = set()
    grid_cnt = 0
    E = {(1,0),(2,0),(3,0),(4,0),(6,0),(8,0),(10,0),(16,0)}
    polygon = sorted(polygon,key = sort)

    for i in range(0,len(polygon),2):
        s1,s2 = polygon[i],polygon[i+1]
        y = max(s1[0][1],s1[1][1])
        P.add(((s1[0][0],y),(s2[0][0],y)))
        grid_cnt += s2[0][0] - s1[0][0]
        
    for s in P:
        for start in range(s[0][0],s[1][0]):
            for v in E:
                if start + v[0] <= s[1][0]:
                    G.add(((start,s[0][1]),(start + v[0],s[0][1])))

    for s1 in G:
        for s2 in G:
            if s1 == s2:
                continue
            elif s1[0][1] == s2[0][1]:
                if intersect(s1,s2):
                    C.add((s1,s2))

    most = pre_optimize(G,C,grid_cnt)
    objective_function = 0
    model = gp.Model()
    var = dict()
    M1 = set()
    M2 = set()

    for s in G:
        x = model.addVar(vtype = gp.GRB.BINARY)
        var.update({s:x})
    
    for s1 in G:
        for s2 in G:
            if s1 == s2:
                continue
            elif (s1,s2) not in M1:
                if abs(s1[0][1] - s2[0][1]) == 1 and intersect(s1,s2):
                    potential = get_intersection(s1,s2) / (get_length(s1) + get_length(s2))
                    objective_function += var[s1] * var[s2] * potential
                    M1.add((s1,s2))
                    M1.add((s2,s1))
    
    for s1 in G:
        for s2 in slope_set:
            if abs(s1[0][1] - s2[0][1]) == 1 and intersect(s1,s2):
                potential = get_intersection(s1,s2) / get_length(s1)
                objective_function += var[s1] * potential

    for s1 in G:
        for s2 in G:
            if s1 == s2:
                continue
            elif (s1,s2) not in M2:
                if abs(s1[0][1] - s2[0][1]) == 1 and cross_detect(s1,s2):
                    if get_length(s1) == get_length(s2):
                        objective_function += var[s1] * var[s2] * 2
                    else:
                        objective_function += var[s1] * var[s2]
                    M2.add((s1,s2))
                    M2.add((s2,s1))

    for s1 in G:
        for s2 in slope_set:
            if abs(s1[0][1] - s2[0][1]) == 1 and cross_detect(s1,s2):
                if get_length(s1) == get_length(s2):
                    objective_function += var[s1] * 2
                else:
                    objective_function += var[s1]

    model.setObjective(objective_function,gp.GRB.MINIMIZE)

    one = gp.LinExpr()
    for s in G:
        if get_length(s) == 1:
            one += var[s]
    model.addConstr(one == most)

    for (s1,s2) in C:
        model.addConstr(var[s1] * var[s2] == 0)

    degree = gp.LinExpr()
    for s in G:
        dx = s[1][0] - s[0][0]
        degree += var[s] * dx
    model.addConstr(degree == grid_cnt)

    model.optimize()

    lego_bar = list()
    for s in G:
        if var[s].x == 1:
            brick = brick_select(s)
            lego_bar.append([brick,(s[0][0],max(s[0][1],s[1][1]))])
            if max(s[0][1],s[1][1]) < min_height:
                min_height = max(s[0][1],s[1][1])
            z = s[0][0] * 20 + (brick.w + 1) * 10
            y = -s[0][1] * 24
            file_content += '0 Step\n' + '1 '
            file_content += '14' + ' '
            file_content += str(0) + ' ' + str(y) + ' ' + str(z)
            file_content += ' 0 0 1 0 1 0 -1 0 0 '
            file_content += brick.dat_code + '\n'
            N1 += 1

    for s1 in G:
        for s2 in G:
            if s1 == s2:
                continue
            if var[s1].x == 1 and var[s2].x == 1:
                if abs(s1[0][1] - s2[0][1]) == 1 and intersect(s1,s2):
                    lego_connect.add(((s1[0][0],s1[0][1]),(s2[0][0],s2[0][1])))

    stability_simulation(N1,N2,min_height,lego_bar,lego_slope,lego_connect)

    with open(file_name,"w") as ldr_file:
        ldr_file.write(file_content)

def get_node(s1,s2,w1,w2):
    x11,x12 = s1[0],s1[0] + w1
    x21,x22 = s2[0],s2[0] + w2
    y = int(640 - min(s1[1],s2[1]) * 20)
    x = [x11,x12,x21,x22]
    x.sort() 
    return (x[1] * 20,y),(x[2] * 20,y)

def get_vertices(brick,point1,point2):
    (x,y),(x1,y1) = point1,point2
    dx,dy = x1 - x,y1 - y
    points = list()

    if dx > 0 and dy > 0:
        points.append((x1 * 20,640 - y1 * 20))
        points.append((x * 20,640 - y * 20))
        points.append(((x + brick.w) * 20,640 - y * 20))
        points.append(((x + brick.w) * 20,640 - y1 * 20))
    elif dx > 0 and dy < 0:
        points.append((x * 20,640 - y * 20))
        points.append(((x - 1) * 20,640 - y * 20))
        points.append(((x - 1) * 20,640 - y1 * 20))
        points.append((x1 * 20,640 - y1 * 20))
    elif dx < 0 and dy > 0:
        points.append((x1 * 20,640 - y1 * 20))
        points.append((x * 20,640 - y * 20))
        points.append(((x + 1) * 20,640 - y * 20))
        points.append(((x + 1) * 20,640 - y1 * 20))
    else:
        points.append((x * 20,640 - y * 20))
        points.append(((x - brick.w) * 20,640 - y * 20))
        points.append(((x - brick.w) * 20,640 - y1 * 20))
        points.append((x1 * 20,640 - y1 * 20))

    return points

def overlap(s1,s2,w):
    (x,y),(x1,y1) = s1
    dx,dy = x1 - x,y1 - y
    point1 = tuple()
    point2 = tuple()
    start,end,top = s2[0],s2[0] + w,s2[1]

    if dx > 0 and dy > 0:
        if top - y1 == 1:
            if (start <= x1 and x1 < end) or x1 == start:
                return get_node((x1,y1),s2,1,w)
        elif y == top:
            if (start <= x and x < end) or (x <= start and start < x1 + 1):
                return get_node((x,y + 1),s2,dx + 1,w)
    elif dx > 0 and dy < 0:
        if top - y == 1:
            if (start <= x - 1 and x - 1 < end) or x - 1 == start:
                return get_node((x - 1,y),s2,1,w)
        elif y1 == top:
            if (start <= x - 1 and x - 1 < end) or (x - 1 <= start and start < x1):
                return get_node((x - 1,y1 + 1),s2,dx + 1,w)
    elif dx < 0 and dy > 0:
        if top - y1 == 1:
            if (start <= x1 and x1 < end) or (x1 <= start and start < x + 1):
                return get_node((x1,y1),s2,1 - dx,w)
        elif y == top:
            if (start <= x and x < end) or x == start:
                return get_node((x,y + 1),s2,1,w)
    elif dx < 0 and dy < 0:
        if top - y == 1:
            if (start <= x1 - 1 and x1 - 1 < end) or (x1 - 1 <= start and start < x):
                return get_node((x1 - 1,y),s2,1 - dx,w)
        elif y1 == top:
            if (start <= x1 - 1 and x1 - 1 < end) or start == x1 - 1:
                return get_node((x1 - 1,y1 + 1),s2,1,w)

    return point1,point2

def stability_simulation(N1,N2,min_height,lego_bar,lego_slope,lego_connect):
    width,height = 640,640
    pg.init()
    screen = pg.display.set_mode((width,height))
    space = pk.Space()
    space.gravity = 0,9.8

    ground = pk.Body(body_type = pk.Body.STATIC)
    shape = pk.Segment(ground,(0,height - (min_height - 1) * 20),(width,height - (min_height - 1) * 20),0)
    space.add(ground,shape)
    slope_bodies = list()
    slope_shapes = list()
    shape_lines = list()
    bodies = list()
    shapes = list()

    for i in range(N1):
        brick = lego_bar[i][0]
        (x,y) = lego_bar[i][1]
        size = (brick.w * 20,20)
        moment = pk.moment_for_box(brick.weight,size)
        body = pk.Body(brick.weight,moment,body_type = pk.Body.DYNAMIC)
        body.position = ((x + brick.w / 2) * 20,height - (y - 0.5) * 20)
        shape = pk.Poly.create_box(body,size)
        shape.friction = 1
        space.add(body,shape)
        bodies.append(body)
        shapes.append(shape)
    
    for i in range(N2):
        brick = lego_slope[i][0]
        (x,y),(x1,y1) = lego_slope[i][1]
        vertices = get_vertices(brick,(x,y),(x1,y1))
        shape_lines.append(vertices)
        moment = pk.moment_for_poly(brick.weight,vertices)
        body = pk.Body(brick.weight,moment,body_type = pk.Body.DYNAMIC)
        shape = pk.Poly(body,vertices)
        (xc,yc) = shape.center_of_gravity
        vertices = [(point[0] - xc,point[1] - yc) for point in vertices]
        shape = pk.Poly(body,vertices)
        shape.friction = 1
        body.position = (xc,yc)
        slope_bodies.append(body)
        slope_shapes.append(shape)
        space.add(body,shape)

    for i in range(N1):
        w1 = lego_bar[i][0].w
        s1 = lego_bar[i][1]
        for j in range(i + 1,N1):
            w2 = lego_bar[j][0].w
            s2 = lego_bar[j][1]
            if (s1,s2) in lego_connect or (s2,s1) in lego_connect:
                point1,point2 = get_node(s1,s2,w1,w2)
                constraint1 = pk.PinJoint(bodies[i],bodies[j],point1)
                constraint2 = pk.PinJoint(bodies[i],bodies[j],point2)
                space.add(constraint1)
                space.add(constraint2)
    
    for i in range(N2):
        for j in range(N1):
            point1,point2 = overlap(lego_slope[i][1],lego_bar[j][1],lego_bar[j][0].w)
            if len(point1) > 0 and len(point2) > 0:
                constraint1 = pk.PinJoint(slope_bodies[i],bodies[j],point1)
                constraint2 = pk.PinJoint(slope_bodies[i],bodies[j],point2)
                space.add(constraint1)
                space.add(constraint2)

        for j in range(i + 1,N2):
            poly1 = geo.Polygon(shape_lines[i])
            poly2 = geo.Polygon(shape_lines[j])
            intersection = poly1.intersection(poly2)
            if intersection.geom_type == "LineString":
                point1 = intersection.coords[0]
                point2 = intersection.coords[-1]
                constraint1 = pk.PinJoint(slope_bodies[i],slope_bodies[j],point1)
                constraint2 = pk.PinJoint(slope_bodies[i],slope_bodies[j],point2)
                space.add(constraint1)
                space.add(constraint2)

    running = True
    while running:  
        for event in pg.event.get():
            if event.type == pg.QUIT:
                running = False
        dt = 1.0 / 240.0
        space.step(dt)
        screen.fill((255,255,255))
        pg.draw.line(screen,(255,0,0),(int(0),int(height - (min_height - 1) * 20)),(int(width),int(height - (min_height - 1) * 20)),2)

        for j in range(N1):
            points = shapes[j].get_vertices()
            points = [(int(bodies[j].position.x + point.x),int(bodies[j].position.y + point.y)) for point in points]
            pg.draw.polygon(screen,(0,255,0),points)

        for i in range(N2):
            (x,y) = (slope_bodies[i].position.x,slope_bodies[i].position.y)
            points = slope_shapes[i].get_vertices()
            points = [(int(point.x + x),int(point.y + y)) for point in points]
            pg.draw.polygon(screen,(0,255,0),points)

        pg.display.flip()
    pg.quit()