import xml.etree.ElementTree

#findall('search') gets all the children where the root's tag matches 'search' 
#get() gets the attribute 
root = xml.etree.ElementTree.parse('../sim_problem.xml').getroot()
for design in root.findall('design'):
    for layer in design.findall('layer'):
        for electrodes in layer.findall('electrode'):
            for property in electrodes:
                if property.tag == "dim":
                    # print property.attrib
                    print 'x1', property.get('x1')
                    print 'y1', property.get('y1')
                    print 'x2', property.get('x2')
                    print 'y2', property.get('y2')
                else:
                    print property.tag, property.text
                    
                    