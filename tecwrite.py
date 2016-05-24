'''
Collection of functions for writing Tecplot data files
Dr. John T. Hwang
November, 2015
'''

def _open_file(filename):
    return open(filename, 'w')

def _write_tec_header(tecfile, title):
    variables = ['x', 'y', 'z']
    tecfile.write('title = ' + title + '\n')
    tecfile.write('variables = ')
    for ivar in range(len(variables)):
        tecfile.write(variables[ivar] + ',')
    tecfile.write('\n')

def _write(tecfile, text):
    tecfile.write(text)

def _write_line(tecfile, data, label=''):
    tecfile.write(label)
    for ind in range(data.shape[0]):
        if data[ind] == data[ind]:
            tecfile.write(str(data[ind]) + ' ')
        else:
            tecfile.write(str(0.0) + ' ')
    tecfile.write('\n')

def _close_file(tecfile):
    tecfile.close()

def _write_surf(tecfile, surf):
    num_u, num_v = surf.shape[:2]
    _write(tecfile,
           'zone i='+str(num_u) + \
           ', j=' + str(num_v) + \
           ', DATAPACKING=POINT\n')
    for ind_v in xrange(num_v):
        for ind_u in xrange(num_u):
            _write_line(tecfile, surf[ind_u, ind_v, :])

def write_surf(filename, surf):
    tecfile = _open_file(filename)
    _write_tec_header(tecfile, 'Data_file')
    _write_surf(tecfile, surf)
    _close_file(tecfile)

def write_surf_multi(filename, surfs):
    tecfile = _open_file(filename)
    _write_tec_header(tecfile, 'Data_file')
    for surf in surfs:
        _write_surf(tecfile, surf)
    _close_file(tecfile)

def write_scatter(filename, data):
    tecfile = _open_file(filename)
    _write_tec_header(tecfile, 'Data_file')
    for ind in xrange(data.shape[0]):
        _write_line(tecfile, data[ind, :])
    _close_file(tecfile)
