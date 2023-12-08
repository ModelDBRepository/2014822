from math import sqrt
from collections import deque
import neuron as nrn
import os
import random
import numpy as np

class Node:
    def __init__(self, item, soma=None, parent=None):
        if(soma == None):
            self.soma = self
        else:
            self.soma = soma
        self.syndist = 9999999
        self.level = 0
        self.parent = None
        self.children = []
        self.x = item.x
        self.item = item
        self.somadist = 0
        self.trunk = True
        if(parent != None):
            self.parent = parent
            self.level = parent.level+1
            self.parent.children.append(self)
            if(parent.trunk and parent.item.diam > 1.5):
                self.trunk = True
            else:
                self.trunk = False
            if(not parent.trunk):
                self.trunk = False
            
        self.nodeid = Node.x_to_id(self)
        self.dist_to_parent = self.calcdistancetoparent()
        
        if self.parent != None:
            self.somadist = self.parent.somadist + self.dist_to_parent
            
        self.dist_to_trunk = self.calc_dist_to_trunk()
        
        p1 = self.getpoint()
        self.point = p1
        p1y = 0
        p2y = 0
        if len(p1) < 1:
            p1 = 0
        else:
            p1y = p1[1]
        p2 = self.soma.getpoint()
        if len(p2) < 1:
            p2 = 0
        else:
            p2y = p2[1]
        self.soma_y_dist = p1y - p2y
        #self.weight_by_distance = self.somadist / 600
        self.weight_by_distance = self.soma_y_dist * 0.002
        
    def calc_neigh_diam(self):
        diam = 0
        diam += self.item.diam
        diam += self.parent.item.diam
        for ch in self.children:
            diam += ch.item.diam
        count = 2 + len(self.children)
        return diam/count
        
    def calc_dist_to_trunk(self):
        node = self
        dist = 0
        while(not node.trunk):
            dist += node.dist_to_parent
            node = node.parent
        return dist
        
    def eucl_dist_to_node(self, node):
        if len(node.point) < 1:
            return 99999
        return Node.eucliddistance(self.point, node.point)
            
    def in_SR_range(self):
        if self.soma_y_dist >= 100 and self.soma_y_dist <= 300:
            return True
        return False
        
    def x_to_id(obj):
        #points = obj.item.sec.psection()['morphology']['pts3d']
        n = obj.item.sec.n3d()
        if(n <= 0):
            return -1
        pid = min(round(n*obj.x), n-1)
        return pid
        
    def getstr(self, level=0):
        pre = " " * level
        pstring = f"{pre}{self.item}[{self.level}]"
        return pstring
    
    def printtree(self, level=0):
        if os.environ["MILEDIDEBUG"] == '1':
            print(self.getstr(level))
        for child in self.children:
            child.printtree(level+1)
    
    def getpoints(self):
        n = self.item.sec.n3d()
        points = []
        for i in range(n):
            point = [
                self.item.sec.x3d(i),
                self.item.sec.y3d(i),
                self.item.sec.z3d(i),
                self.item.sec.diam
            ]
            points.append(point)
        return points
        
    def getpoint(self):
        if(self.nodeid == -1):
            return []
        points = self.getpoints()
        p = points[self.nodeid]
        #return [p[0], p[1], p[2], self.level]
        return [p[0], p[1], p[2], p[3]]
    
    def gettreepoints(self):
        points = []
        cp = self.getpoint()
        if (len(cp) > 0):
            points = [cp]
        for child in self.children:
            cp = child.gettreepoints()
            if (len(cp) > 0):
                points += cp
        return points
    
    def getpointstox(self, x, y):
        points = self.getpoints()
        n = len(points)
        if(n <= 0):
            return []
        pid1 = min(round(n*y), n-1)
        pid2 = min(round(n*x), n-1)
        if pid1 > pid2:
            temp = pid1
            pid1 = pid2
            pid2 = temp
            
        p = []
        for i in range(pid1, pid2+1):
            p.append(points[i])
        return p
    
    def eucliddistance(p1, p2):
        x1 = p1[0]
        x2 = p2[0]
        y1 = p1[1]
        y2 = p2[1]
        z1 = p1[2]
        z2 = p2[2]
        
        ssum = (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2        
        return sqrt(ssum)
        
    def calcdistancefrompoints(points):
        dist = 0
        for pid in range(len(points)-1):
            dist += Node.eucliddistance(points[pid], points[pid+1])
        return dist
    
    def calcdistancetoparent(self):
        if(self.parent == None):
            return 0
        
        isself = self.parent.item.sec.same(self.item.sec)
        otherx = self.parent.x
        
        dist = 0
        if(not isself):
            otherx = 0
            points = self.parent.getpointstox(self.parent.x, 1)
            dist += Node.calcdistancefrompoints(points)
            
        points = self.getpointstox(otherx, self.x)
        dist += Node.calcdistancefrompoints(points)
        return dist
    
    def get_neighbour_dists(self, caller=None):
        nodes = []
        if caller != self.parent:
            if self.parent != None:
                nodes.append((self.parent, self.dist_to_parent))
        
        for child in self.children:
            if caller != child:
                if child != None:
                    nodes.append((child, child.dist_to_parent))
            
        return nodes
    
    def adjust_synapse_distances(self, caller=None):
        neighbours = self.get_neighbour_dists(caller)
        for neigh in neighbours:
            node = neigh[0]
            node_dist = neigh[1]
            node.syndist = min(node.syndist, node_dist + self.syndist)
            node.adjust_synapse_distances(self)
                
    def foreach(self, func, reverse=False):
        if(not reverse):
            func(self)
        for child in self.children:
            child.foreach(func)
        if(reverse):
            func(self)

    def function_for_neighbors_in_distance(self, function, dist, caller=None):
        neighbours = self.get_neighbour_dists(caller)
        for neigh in neighbours:
            node = neigh[0]
            node_dist = neigh[1]
            if node_dist < dist:
                function(node)
                node.function_for_neighbors_in_distance(function, dist-node_dist, self)

            
    def __str__(self):
        return str(self.item)

    def __repr__(self):
        return str(self.item)


def buildNodes(soma):
    que = deque()
    allitems = {}

    prev = None
    que.append(soma)
    root = Node(soma(0.5))
    allitems[str(soma(0.5))] = root

    while (len(que) > 0):
        item = que.pop()
        que.extend(item.children())
        for itemyi in list(item.allseg()):
            itemy = str(itemyi)
            tps = itemyi.sec.trueparentseg()
            if (tps == None):
                if(itemyi.x != 0.5):
                    allitems[itemy] = Node(itemyi, root, parent=allitems[str(itemyi.sec(0.5))])
            else:
                if itemyi.sec.same(prev.item.sec):
                    allitems[itemy] = Node(itemyi, root, parent = prev)
                else:
                    allitems[itemy] = Node(itemyi, root, parent = allitems[str(tps)])

            prev = allitems[itemy]
    return allitems

def get_nodes_in_radius(nodes, center, radius=20):
    nodes_in_radius = []
    for n1 in nodes:
        if center.eucl_dist_to_node(n1) < radius:
            nodes_in_radius.append(n1)
        
    return nodes_in_radius


def select_random_synapse(allnodes, somadist_L=None, trunkdist_L=None, syn_diam=None, seed=None, centerneuron=None):
    if seed == None:
        random.seed()
    else:
        random.seed(seed)

    ''' Select nodes that are:
    1. in Stratum Radiatum
    2. not trunk (diameter > x)
    3. not edge of segment (0 or 1), since it doesn't allow us to play synapse on there
    '''
    possible_synapse_nodes = [node for node in allnodes.values() if node.in_SR_range() and not node.trunk and not node.x == 0 and not node.x == 1]
    possible_synapse_nodes_distance = possible_synapse_nodes
    
    if somadist_L is not None:
        possible_synapse_nodes_distance = [x for x in possible_synapse_nodes_distance if abs(x.soma_y_dist - somadist_L) < 20]
        if len(possible_synapse_nodes_distance) < 1:
            possible_synapse_nodes_distance = sorted(possible_synapse_nodes, key=lambda x: abs(x.soma_y_dist - somadist_L))[:5]
        if (os.environ.get("NRN_DEBUG") == "1"):
            print("1. Selected distance from soma nodes: ")
            for x, dist in [(x, x.soma_y_dist) for x in possible_synapse_nodes_distance]:
                print(f" {f'{str(x)}:'.ljust(42)}  {f'{dist:.2f}'.ljust(5)}")
        #centernode = possible_synapse_nodes[0]
    #    else:


    if trunkdist_L is not None:
        possible_synapse_nodes_distance_trunk = [x for x in possible_synapse_nodes_distance if abs(x.dist_to_trunk - trunkdist_L) < 20]
        if len(possible_synapse_nodes_distance_trunk) < 1:
            possible_synapse_nodes_distance = sorted(possible_synapse_nodes_distance, key=lambda x: abs(x.dist_to_trunk - trunkdist_L))[:5]
        else:
            possible_synapse_nodes_distance = possible_synapse_nodes_distance_trunk
        if (os.environ.get("NRN_DEBUG") == "1"):
            print("2. Selected distance from trunk nodes: ")
            for x, dist in [(x, x.dist_to_trunk) for x in possible_synapse_nodes_distance]:
                print(f" {f'{str(x)}:'.ljust(42)}  {f'{dist:.2f}'.ljust(5)}")
        

    if syn_diam is not None:
        possible_synapse_nodes_distance_diam = [x for x in possible_synapse_nodes_distance if abs(x.calc_neigh_diam() - syn_diam) < 20]
        if len(possible_synapse_nodes_distance_diam) < 1:
            possible_synapse_nodes_distance = sorted(possible_synapse_nodes_distance, key=lambda x: abs(x.calc_neight_diam() - syn_diam))[:5]
        else:
            possible_synapse_nodes_distance = possible_synapse_nodes_distance_diam
        if (os.environ.get("NRN_DEBUG") == "1"):
            print("3. Selected dendrite diameters: ")
            for x, dist in [(x, x.calc_neigh_diam()) for x in possible_synapse_nodes_distance]:
                print(f" {f'{str(x)}:'.ljust(42)}  {f'{dist:.2f}'.ljust(5)}")


            
    centernodes = random.sample(possible_synapse_nodes_distance, 1)

    if centerneuron is not None:
        centernodes = [centerneuron]
    

    if os.environ["MILEDIDEBUG"] == '1':
        print("Final cluster: ")
        for x in [x for x in centernodes]:
            print(f" {x}")
        print("Center node: ")
        print(f" {centernodes[0]}")
    
    return (centernodes, centernodes[0])


def create_random_cluster(allnodes, syn_density, syn_count, somadist_L=None, trunkdist_L=None, syn_diam=None, seed=None, centerneuron=None):
    if syn_count == 1:
        return select_random_synapse(allnodes, somadist_L, trunkdist_L, syn_diam, seed, centerneuron)
    if seed == None:
        random.seed()
    else:
        random.seed(seed)


    possible_synapse_nodes = [node for node in allnodes.values() if node.in_SR_range() and not node.trunk and not node.x == 0 and not node.x == 1]
    possible_synapse_nodes_distance = possible_synapse_nodes
    
    if somadist_L is not None:
        possible_synapse_nodes_distance = [x for x in possible_synapse_nodes_distance if abs(x.soma_y_dist - somadist_L) < 20]
        if len(possible_synapse_nodes_distance) < 1:
            possible_synapse_nodes_distance = sorted(possible_synapse_nodes, key=lambda x: abs(x.soma_y_dist - somadist_L))[:5]


    if trunkdist_L is not None:
        possible_synapse_nodes_distance_trunk = [x for x in possible_synapse_nodes_distance if abs(x.dist_to_trunk - trunkdist_L) < 20]
        if len(possible_synapse_nodes_distance_trunk) < 1:
            possible_synapse_nodes_distance = sorted(possible_synapse_nodes_distance, key=lambda x: abs(x.dist_to_trunk - trunkdist_L))[:5]
        else:
            possible_synapse_nodes_distance = possible_synapse_nodes_distance_trunk
        

    if syn_diam is not None:
        possible_synapse_nodes_distance_diam = [x for x in possible_synapse_nodes_distance if abs(x.calc_neigh_diam() - syn_diam) < 20]
        if len(possible_synapse_nodes_distance_diam) < 1:
            possible_synapse_nodes_distance = sorted(possible_synapse_nodes_distance, key=lambda x: abs(x.calc_neight_diam() - syn_diam))[:5]
        else:
            possible_synapse_nodes_distance = possible_synapse_nodes_distance_diam


    centernodes = random.sample(possible_synapse_nodes_distance, 1)

    if centerneuron is not None:
        centernodes = [centerneuron]
        
    radius = 30
    selected_count = 4
    
    circle_nodes = get_nodes_in_radius(possible_synapse_nodes, centernodes[0], radius)
    #centernodes += random.sample(circle_nodes, selected_count)
    
    cluster = create_cluster(possible_synapse_nodes, centernodes, syn_density, syn_count, circle_nodes, somadist_L)
    return (cluster, centernodes[0])


def create_cluster(allnodes, centernodes, syn_density, syn_count, circle_nodes, somadist_L):
    for node in allnodes:
        node.syndist = 999999
    
    distance = 1/syn_density
    if type(centernodes) is not list:
        cluster_synapses = [centernodes]
    else:
        cluster_synapses = centernodes
        
    syn_count -= len(cluster_synapses)
    
    for node in cluster_synapses:
        node.syndist = 0

    for node in cluster_synapses:
        node.adjust_synapse_distances()
        
    while syn_count > 0:
        syndistances = [[node, node.syndist, centernodes[0].eucl_dist_to_node(node) + abs(node.soma_y_dist - somadist_L) * 5, abs(distance - node.syndist)] for node in allnodes if node.in_SR_range() and node.syndist > 0 and node.syndist > (0.9*distance)]
        #syndistances += [[node, 0, centernodes[0].eucl_dist_to_node(node), abs(distance - node.syndist)] for node in allnodes if node.syndist > 0]
        #print(syndistances)
        
        syndistances = sorted(syndistances, key=lambda x: x[3])
        if syndistances[0][1] > (distance * 2.0):
            syndistances = sorted(syndistances, key=lambda x: x[2])
        
        added_synapses = random.sample(syndistances[0:5], 1)
        for added_synapse in [x[0] for x in added_synapses]:
            if (syn_count <= 0):
                break
            added_synapse.syndist = 0
            cluster_synapses.append(added_synapse)
            syn_count -= 1
            cluster_synapses[-1].adjust_synapse_distances()
    
    return cluster_synapses

