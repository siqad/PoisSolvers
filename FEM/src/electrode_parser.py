import xml.etree.ElementTree

#findall('search') gets all the children where the root's tag matches 'search' 
#get() gets the attribute 

def xml_parse(path):
    m_per_A = 1.0E-10
    root = xml.etree.ElementTree.parse(path).getroot()
    layer_props = {}
    for layer_prop in root.findall("layer_prop"):    
        desired = False
        for prop in layer_prop:
            if prop.tag == "name" and prop.text == "Metal":
                desired = True
            if desired == True:
                layer_props[prop.tag] = prop.text
    elec_list = []
    for design in root.findall('design'):
        for layer in design.findall('layer'):
            for electrodes in layer.findall('electrode'):
                curr_elec = {}
                for property in electrodes:
                    if property.tag == "dim":
                        curr_elec['x1'] = float(property.get('x1'))
                        curr_elec['y1'] = float(property.get('y1'))
                        curr_elec['x2'] = float(property.get('x2'))
                        curr_elec['y2'] = float(property.get('y2'))
                    else:
                        curr_elec[str(property.tag)] = property.text
                elec_list.append(curr_elec)
    for elec in elec_list:
        elec['x1'] *= m_per_A/float(elec['pixel_per_angstrom'])
        elec['y1'] *= m_per_A/float(elec['pixel_per_angstrom'])
        elec['x2'] *= m_per_A/float(elec['pixel_per_angstrom'])
        elec['y2'] *= m_per_A/float(elec['pixel_per_angstrom'])
    return elec_list, layer_props
    
def getBB(elec_list):
    x1s = [x['x1'] for x in elec_list]
    y1s = [x['y1'] for x in elec_list]
    x2s = [x['x2'] for x in elec_list]
    y2s = [x['y2'] for x in elec_list]

    min_x = min(x1s)
    min_y = min(y1s)
    max_x = max(x2s)
    max_y = max(y2s)
    
    min_x -= (max_x - min_x)*0.5
    max_x += (max_x - min_x)*0.5
    min_y -= (max_y - min_y)*0.5
    max_y += (max_y - min_y)*0.5
    return [min_x, max_x], [min_y, max_y]

def main():
    elec_list, layer_props = xml_parse("../sim_problem.xml")
    print elec_list
    print layer_props
    [min_x, max_x], [min_y, max_y] = getBB(elec_list)
    print min_x, max_x, min_y, max_y
    
# main()