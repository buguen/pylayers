from xml.etree.ElementTree import parse

class Worldsens(object):
    def __init__(self,Nnodesi=13,Tsim=50,Dx=2,Dy,2,Dz=2):
        self.Tsim = Tsim
        self.Dx = Dx
        self.Dy = Dy
        self.Dz = Dz
        self.Tsim = Tsim

    def write(self):

class Entities(object):
    def __init__(self):
        pass

class Environment(object):
    def __init__(self):
        pass

class Bundle(object):
    def __init__(self):
        pass

class Nodes(object):
    def __init__(self):
        pass


data = parse('ban_cea.xml')
for e in data.iter():
    print e.tag, e.attrib


