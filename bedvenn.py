import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')  # for saving image file
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

from pybedtools import BedTool


def get_num_bp(bt):
    return np.sum([i.length for i in bt])


def get_label(filepath):
    """Return filename without extension"""
    return os.path.splitext(os.path.basename(filepath))[0]


f1 = sys.argv[1]
f2 = sys.argv[2]
b1 = BedTool(f1)
b2 = BedTool(f2)

if len(sys.argv) == 4:
    out_fn = sys.argv[3]
    b12_intersect = b1.intersect(b2)
    n_12 = get_num_bp(b12_intersect)
    n_1 = get_num_bp(b1) - n_12
    n_2 = get_num_bp(b2) - n_12
    subsets = {'01': n_2,
               '10': n_1,
               '11': n_12}
    venn2(subsets=subsets, set_labels=(get_label(f1),
                                       get_label(f2)))
elif len(sys.argv) == 5:
    out_fn = sys.argv[4]
    f3 = sys.argv[3]
    b3 = BedTool(f3)

    b12_intersect = b1.intersect(b2)
    b13_intersect = b1.intersect(b3)
    b23_intersect = b2.intersect(b3)
    b123_intersect = b12_intersect.intersect(b3)

    n_123 = get_num_bp(b123_intersect)
    n_12 = get_num_bp(b12_intersect) - n_123
    n_13 = get_num_bp(b13_intersect) - n_123
    n_23 = get_num_bp(b23_intersect) - n_123
    n_1 = get_num_bp(b1) - n_12 - n_13 - n_123
    n_2 = get_num_bp(b2) - n_12 - n_23 - n_123
    n_3 = get_num_bp(b3) - n_13 - n_23 - n_123

    subsets = {'001': n_3,
               '010': n_2,
               '011': n_23,
               '100': n_1,
               '101': n_13,
               '110': n_12,
               '111': n_123}

    venn3(subsets=subsets, set_labels=(get_label(f1),
                                       get_label(f2),
                                       get_label(f3)))
else:
    raise Exception('Must provide 3 BED files.')

plt.savefig(out_fn)
