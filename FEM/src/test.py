import xml.etree.ElementTree

#findall('search') gets all the children where the root's tag matches 'search' 
#get() gets the attribute 
root = xml.etree.ElementTree.parse('../sim_problem.xml').getroot()
for design in root.findall('design'):
    for layer in design.findall('layer'):
        for electrodes in layer.findall('electrode'):
            for property in electrodes:
                if property.tag == "dim":
                    print property.attrib
                    print property.get('x1')
                    print property.get('y1')
                    print property.get('x2')
                    print property.get('y2')
                else:
                    print property.tag, property.text
                    
                    