#!/usr/bin/env python3
import os
import optparse


def gen_rst(directory):
    """Go through all fortran files in directory and generate the rst doxygen
    function calls.   Directory must be full path"""
    files = os.listdir(directory)
    files.sort()

    rstfile = os.path.join(directory, "../docs/gen_subroutines.rst")
    with open(os.path.join(directory, "../docs/subroutines.header"), "r") as tfh:
        header = tfh.read()

    with open(rstfile, "w") as fh:
        fh.write(header)
        for file in files:
            if "modules-efit" in file or 'chkerr' in file:
                # contains subroutines inside of modules which currently breaks
                # doxygen (could possibly be fixed, but they aren't essential anyway)
                continue
            # print(file)
            if file.lower().endswith("f90"):
                first = 0
                with open(os.path.join(directory, file), "r") as f:
                    filelines = f.readlines()
                    for line in filelines:
                        linestrip = line.strip()
                        if linestrip[0:10].lower() == "subroutine":
                            if first == 0:
                                fh.write(file)
                                divstring = ""
                                for i in range(len(file)):
                                    divstring += "-"
                                fh.write(divstring + "\n\n")
                                first += 1
                            fh.write(
                                ".. doxygenfunction:: "
                                + line.strip().split()[1].split("(")[0].lower()
                                + "\n",
                            )
                fh.write("\n")

    rstfile = os.path.join(directory, "../docs/gen_functions.rst")
    with open(os.path.join(directory, "../docs/functions.header"), "r") as tfh:
        header = tfh.read()

    with open(rstfile, "w") as fh:
        fh.write(header)
        for file in files:
            # print(file)
            if file.lower().endswith("f90"):
                first = 0
                with open(os.path.join(directory, file), "r") as f:
                    filelines = f.readlines()
                    for line in filelines:
                        linestrip = line.strip()
                        if linestrip[0:8].lower() == "function":
                            if first == 0:
                                fh.write(file)
                                divstring = ""
                                for i in range(len(file)):
                                    divstring += "-"
                                fh.write(divstring + "\n\n")
                                first += 1
                            fh.write(
                                ".. doxygenfunction:: "
                                + line.strip().split()[1].split("(")[0].lower()
                                + "\n",
                            )
                fh.write("\n")

    rstfile = os.path.join(directory, "../docs/gen_modules.rst")
    with open(os.path.join(directory, "../docs/modules.header"), "r") as tfh:
        header = tfh.read()

    toublesome_modules = ['var_nio','error_control']
    with open(rstfile, "w") as fh:
        fh.write(header)
        for file in files:
            # print(file)
            if file.lower().endswith("f90"):
                first = 0
                with open(os.path.join(directory, file), "r") as f:
                    filelines = f.readlines()
                    for line in filelines:
                        linestrip = line.strip()
                        if linestrip[0:6].lower() == "module":
                            if first == 0:
                                fh.write(file)
                                divstring = ""
                                for i in range(len(file)):
                                    divstring += "-"
                                fh.write(divstring + "\n\n")
                                first += 1
                            fh.write(
                                ".. doxygennamespace:: "
                                + line.strip().split()[1].split("(")[0].lower()
                                + "\n",
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
