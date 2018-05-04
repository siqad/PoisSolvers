import xml.etree.ElementTree

#findall('search') gets all the children where the root's tag matches 'search' 
#get() gets the attribute 

def xml_parse(path):
    m_per_A = 1.0E-10
    root = xml.etree.ElementTree.parse(path).getroot()
    sim_params = {}
    for sim_param in root.findall("sim_params"):
        for param in sim_param:
            sim_params[param.tag] = param.text
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
    return elec_list, layer_props, sim_params
    
def getBB(elec_list):
    x1s = [a['x1'] for a in elec_list]
    y1s = [a['y1'] for a in elec_list]
    x2s = [a['x2'] for a in elec_list]
    y2s = [a['y2'] for a in elec_list]

    min_x = min(x1s)
    min_y = min(y1s)
    max_x = max(x2s)
    max_y = max(y2s)
    
    xs = [min_x-2.0*(max_x-min_x), max_x+2.0*(max_x-min_x)]
    ys = [min_y-2.0*(max_y-min_y), max_y+2.0*(max_y-min_y)]
    return xs, ys
    
def getZparams(layer_props):
    return float(layer_props['zheight']), float(layer_props['zoffset'])
    
def getResolutionScale(sim_params):
    return float(sim_params['resolution'])
    
def getPixPerAngstrom(elec_list):
    return float(elec_list[0]['pixel_per_angstrom'])

def main():
    elec_list, layer_props, sim_params = xml_parse("../sim_problem.xml")
    print elec_list
    print layer_props
    [min_x, max_x], [min_y, max_y] = getBB(elec_list)
    print min_x, max_x, min_y, max_y
    print sim_params
    
# main()