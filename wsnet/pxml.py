from xml.etree.ElementTree import parse
data = parse('ban_cea.xml')
for e in data.iter():
    print e.tag, e.attrib


