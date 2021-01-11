import numpy as np
from matplotlib import pyplot as plt
from camfr import *


def get_field(section, field, mode=0,resolution = 0.001):
    dy = resolution
    dx = dy
    obj = section
    x = np.arange(0, obj.width(), dx)
    y = np.arange(0, obj.height(), dy)
    X,Y = np.meshgrid(x,y)

    F = np.zeros((np.size(y),np.size(x)), float)
    # print np.size(y), np.size(x)
    # print X.shape, Y.shape
    for i in range(0, np.size(y)):
        for j in range(0, np.size(x)):
                    
            F[i,j] =  \
                np.abs( 
                    getattr(    obj.mode(mode).field(Coord(X[i,j], Y[i,j], 0)) , field)()
                ) ** 2.0
    return F.astype(np.float)

def get_refr(section, resolution = 0.001):
    dy = resolution
    dx = dy
    obj = section
    x = np.arange(0, obj.width(), dx)
    y = np.arange(0, obj.height(), dy)
    X,Y = np.meshgrid(x,y)

    R = np.zeros((np.size(y),np.size(x)), float)
    print np.size(y), np.size(x)
    print X.shape, Y.shape
    for i in range(0, np.size(y)):
        for j in range(0, np.size(x)):
            R[i,j] = np.real(obj.n(Coord(X[i,j],Y[i,j],0)))
    return R.astype(np.float)
        
def max_inside(F, R, n_core):
    max_pos = unravel_index(F.argmax(), F.shape)
    return R[max_pos] == n_core

def plot_max(F, R, dim=None):
    (my, mx) = unravel_index(F.argmax(), F.shape)
    fig, ax = subplots(1)
    c = plt.Circle((mx, my),1,fill=False, color='green')
    
    if dim !=None:
        im = ax.imshow(F/F.max(),cmap='plasma', extent = dim) #/F.max()
        # ax.set_aspect(dim[3]/dim[1])
        ax.imshow(R,cmap='gray',alpha=0.2, extent = dim)
        ax.set_aspect(dim[3]/dim[1])
    else:
        im = ax.imshow(F/F.max(),cmap='plasma') #/F.max()
        ax.imshow(R, cmap='gray', alpha=0.2)
    # ax.add_patch(c)
    ax.set_xlabel('nm')
    ax.set_ylabel('nm')
    fig.colorbar(im)
    plt.show()

def plot_modes(section, modes, fields, normalize = False, resolution = 0.001):
    rows = len(modes)
    cols = len(fields)
    fig, axes = plt.subplots(rows, cols)
    images = []
    R = get_refr(section, resolution=resolution)
    for row in range(rows):
        for col in range(cols):
            mode = modes[row]
            field = fields[col]

            F = get_field(section, mode = mode, field=field, resolution=resolution)
            images.append(F)
    
    if normalize:
        vmin, vmax = 0, 0
        for im in images:
            if vmin > im.min():
                vmin = im.min()
            if vmax < im.max():
                vmax = im.max()

    for row in range(rows):
        for col in range(cols):
            mode = modes[row]
            field = fields[col]
            F = images[row*cols + col]
            ax = axes[row][col]
            if normalize:
                ax.imshow(F, cmap='hot', vmin=vmin, vmax=vmax)
            else:
                ax.imshow(F, cmap='hot')
            ax.imshow(R, cmap='gray',alpha=0.2)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title('{}-{}'.format(mode, field))

    plt.show()

def calc_section(n_core, n_clad, l, k, modes_count=5):
    wl = 0.01352
    n_clad += k
    set_lambda(wl)
    set_N(modes_count)
    core = Material(n_core)
    clad = Material(n_clad)
    set_section_solver(L)
    set_mode_correction(guided_only) #it can be also 'full'
    set_left_wall(E_wall)
    set_right_wall(E_wall)
    set_right_PML(-0.001)
    set_upper_PML(-0.001)
    set_lower_wall(slab_E_wall)
    set_upper_wall(slab_E_wall)
    wg = Slab(clad(l) + core(l)+clad(l))
    air = Slab(clad(wg.width()))
    all = air(l)+ wg(l) + air(l)
    s = Section(all,10,40)
    s.calc()
    return s


class FilterModes:
    def __init__(self,section, field, modes, resolution, threshold = 0.2):
        '''
            Attrs:
                fields_correct: list of correct fields
                n_effs: list of n_effs of correct fields
                fields_all: all computed fields
                bool_fields [bool]: list of correct and incorrect fields
        '''
        
        def maxdiff(f1, f2):
            f1 = abs(f1.copy())
            f2 = abs(f2.copy())
            f1 = f1 / float(f1.max())
            f2 = f2 / float(f2.max())

            res = abs(f2 - f1)
            return res.max()

        n_effs = []
        fields_all = []
        bool_fields = []
        for i in modes:
            F = get_field(section, mode = i, field=field, resolution=resolution)
            n_eff = section.mode(i).n_eff()
            fields_all.append(F)
            n_effs.append(n_eff)
        fields_correct = list(fields_all)
        inds_to_del = []
        similarities = {}
        for m in modes:
            try:
                loc = []
                for mc in modes:
                    if m != mc:
                        mxdf = maxdiff(fields_all[m], fields_all[mc])
                        if mxdf < threshold:
                            loc.append(mc,fields_all[m], fields_all[mc])
                similarities[m] = loc
            except:
                pass
        similars = []
        for mode in similarities:
            try:
                for sim in similarities[mode]:
                    loc = [mode, sim[0]]
                    loc.sort()
                    if loc not in similars:
                        similars.append(loc)
            except Exception as e:
                print(e)

        print('similaRS!!!!', similars)
        if len(similars) != 0:
            inds_to_del = list(array(similars)[:,1])
        for fi, field in enumerate(fields_all):
            bool_fields.append(fi not in inds_to_del)
        print('\n'*3)
        n_effs_correct = list(n_effs)
        inds_to_del.reverse()
        for i in inds_to_del:
            fields_correct.pop(i)
            n_effs_correct.pop(i)

        self.fields_correct = fields_correct
        self.n_effs_correct = n_effs_correct
        self.n_effs         = n_effs
        self.fields_all     = fields_all
        self.bool_fields    = bool_fields
        self.similarities   = similarities
        # return fields_correct, n_effs_correct, n_effs, fields_all, bool_fields