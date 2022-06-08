#!/usr/bin/env python
"""
Class and methods for getting the Green's function data
"""
import os
import glob
import optparse
import h5py
import utils


class ExtractData:
    def __init__(self, options=None):
        if options:
            self.verbose = options.verbose
        else:
            self.verbose = False

        #  Place holders to show what data will be filled
        self.root_name = None

    def _get_shotlabel(self, h5in):
        """
        Find the machine and shot to label the root directory
        """
        metagrp = h5in.get("dataset_description/data_entry")
        if "machine" in metagrp.keys():
 
    def read_green():
        nesum = 6
        nsilop = 44
        magpr2 = 76
        nfsum = 18
        nvsum = 24
        nw=129
        nh=129

        filedir ='/fusion/projects/codes/efit/efitai/efit_support_files/DIII-D//green/168191/'
        endian = '>'

        # Response of poloidal field (F) coils on grid, and grid on itself.
        filename = f'ec{nw}{nh}.ddd'
        file = filedir+filename
        f=scipy.io.FortranFile(file,'r',f'{endian}i4')
        [mw,mh] = f.read_ints(f'{endian}i4')
        grid = f.read_reals(f'{endian}f8')
        rgrid = grid[:mw]
        zgrid = grid[mh:]
        ggridfc = f.read_reals(f'{endian}f8').reshape([nfsum,mw*mh])
        gridpc = f.read_reals(f'{endian}f8').reshape(mh,mw*mh)

        # Response of Ohmic heating coil on flux loops, magnetic probes, and grid
        filename = f're{nw}{nh}.ddd'
        file = filedir+filename
        f=scipy.io.FortranFile(file,'r',f'{endian}i4')
        rsilec = f.read_reals(f'{endian}f8').reshape([nesum,nsilop])
        rmp2ec = f.read_reals(f'{endian}f8').reshape([nesum,magpr2])
        gridec = f.read_reals(f'{endian}f8').reshape([nesum,mw*mh])

        # Reponse of vessel on flux loops, probes, and grid
        filename = f'rv{nw}{nh}.ddd'
        file = filedir+filename
        f=scipy.io.FortranFile(file,'r',f'{endian}i4')
        gsilvs= f.read_reals(f'{endian}f8').reshape([nvsum,nsilop])
        gmp2vs = f.read_reals(f'{endian}f8').reshape([nvsum,magpr2])
        ggridvs = f.read_reals(f'{endian}f8').reshape([nvsum,mw*mh])

        # Reponse of magnetics flux loops and probes on grid
        filename = f'ep{nw}{nh}.ddd'
        file = filedir+filename
        f=scipy.io.FortranFile(file,'r',f'{endian}i4') 
        rsilpc = f.read_reals(f'{endian}f8').reshape([mw*mh, nsilop])
        rmp2pc = f.read_reals(f'{endian}f8').reshape([mw*mh, magpr2])


def parse_greenargs():
    """
    Routine for getting the options and arguments for extracting equilibrium.py
    It is it's own routine to allow other scripts (like plot_equilibrium.py)
    to use it.
    """
    parser = optparse.OptionParser(usage="%prog [options] omas_datafile(s)")
    parser.add_option(
        "-v", "--verbose", help="Verbose output", dest="verbose", action="store_true"
    )
    return parser


def main():
    """
    Parse arguments and options and act accordingly
    """
    parser = parse_eqargs()
    options, args = parser.parse_args()

    # Sanity checks
    if len(args) < 1:
        parser.print_usage()
        return

    input_directory = args[0]

    if not os.path.isdir(input_directory):
        print("Input directory ", input_directory, " must exist.")
        return
    
    # Show how it's used
    if options.verbose:
        print("Running code")

    read_green(input_directory)

if __name__ == "__main__":
    main()

