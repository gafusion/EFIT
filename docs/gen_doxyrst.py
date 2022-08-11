#!/usr/bin/env python
import os
import optparse


def gen_rst(directory):
    """Go through all fortran files in directory and generate the rst doxygen
    function calls.   Directory must be full path"""
    files = os.listdir(directory)
    files.sort()

    rstfile=os.path.join(directory,'../docs/gen_subroutines.rst')
    with open(rstfile, "w") as fh:
        for file in files:
            print(file)
            if "f90" in file or "F90" in file and "swp" not in file:
                fh.write(file)
                divstring = ""
                for i in range(len(file)):
                    divstring += "-"
                fh.write(divstring + "\n\n")
                with open(os.path.join(directory, file), "r") as f:
                    filelines = f.readlines()
                    for line in filelines:
                        linestrip = line.strip()
                        if linestrip[0:10].lower() == "subroutine":
                            fh.write(
                                ".. doxygenfunction::"+
                                line.strip().split()[1].split("(")[0].lower()+"\n",
                            )
                fh.write("\n")


def main():
    """
    Parse arguments and options and act accordingly
    """
    parser = optparse.OptionParser(usage="%prog [options] directory")
    options, args = parser.parse_args()

    # Sanity checks
    if not len(args) == 1:
        parser.print_usage()
        return

    directory = args[0]

    if not os.path.isdir(directory):
        print("Input directory ", directory, " must exist.")
        return

    gen_rst(directory)


if __name__ == "__main__":
    main()
