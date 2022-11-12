read_green = False
plot_mag = False
update_boundary = False
iprobe = 20

nesum = 6
nsilop = 44
magpr2 = 76
nfsum = 18
nvsum = 24
nw = 129
nh = 129

filedir = "/fusion/projects/codes/efit/efitai/efit_support_files/DIII-D//green/168191/"
endian = ">"
if read_green:
    # Response of poloidal field (F) coils on grid, and grid on itself.
    filename = f"ec{nw}{nh}.ddd"
    file = filedir + filename
    f = scipy.io.FortranFile(file, "r", f"{endian}i4")
    [mw, mh] = f.read_ints(f"{endian}i4")
    grid = f.read_reals(f"{endian}f8")
    rgrid = grid[:mw]
    zgrid = grid[mh:]
    ggridfc = f.read_reals(f"{endian}f8").reshape([nfsum, mw * mh])
    gridpc = f.read_reals(f"{endian}f8").reshape(mw * mh, mh)

    # Response of Ohmic heating coil on flux loops, magnetic probes, and grid
    filename = "rfcoil.ddd"
    file = filedir + filename
    f = scipy.io.FortranFile(file, "r", f"{endian}i4")
    rsilfc = f.read_reals(f"{endian}f8").reshape([nfsum, nsilop])
    rmp2fc = f.read_reals(f"{endian}f8").reshape([nfsum, magpr2])

    # Response of Ohmic heating coil on flux loops, magnetic probes, and grid
    filename = f"re{nw}{nh}.ddd"
    file = filedir + filename
    f = scipy.io.FortranFile(file, "r", f"{endian}i4")
    rsilec = f.read_reals(f"{endian}f8").reshape([nesum, nsilop])
    rmp2ec = f.read_reals(f"{endian}f8").reshape([nesum, magpr2])
    gridec = f.read_reals(f"{endian}f8").reshape([nesum, mw * mh])

    # Response of vessel on flux loops, probes, and grid
    filename = f"rv{nw}{nh}.ddd"
    file = filedir + filename
    f = scipy.io.FortranFile(file, "r", f"{endian}i4")
    gsilvs = f.read_reals(f"{endian}f8").reshape([nvsum, nsilop])
    gmp2vs = f.read_reals(f"{endian}f8").reshape([nvsum, magpr2])
    ggridvs = f.read_reals(f"{endian}f8").reshape([nvsum, mw * mh])

    # Response of magnetics flux loops and probes on grid
    filename = f"ep{nw}{nh}.ddd"
    file = filedir + filename
    f = scipy.io.FortranFile(file, "r", f"{endian}i4")
    rsilpc = f.read_reals(f"{endian}f8").reshape([mw * mh, nsilop])
    rmp2pc = f.read_reals(f"{endian}f8").reshape([mw * mh, magpr2])

if plot_mag:

    plt.figure()
    contourf(RR, ZZ, rmp2pc[:, 0].reshape([nw, nh]).T, levels=20)
    plt.title("probe 0")
    plt.xlabel("R grid")
    plt.ylabel("Z grid")
    colorbar()

    plt.figure()
    contourf(RR, ZZ, sum(rmp2pc[:, :], axis=1).reshape([nw, nh]).T, levels=20)
    plt.title("sum of all probes")
    plt.xlabel("R grid")
    plt.ylabel("Z grid")
    colorbar()


def isPointInPath(x, y, poly):
    """
    Function for deteriming if given x,y coord is inside a boundary defined
    by the polynomial poly

    Inputs:
    -------
      x, y : scalar variable, RZ coord. of the point of interest
      poly : boundary shape
    """
    num = len(poly)
    i = 0
    j = num - 1
    c = False

    for i in range(num):
        if ((poly[i][1] > y) != (poly[j][1] > y)) and (
            x
            < (poly[j][0] - poly[i][0]) * (y - poly[i][1]) / (poly[j][1] - poly[i][1])
            + poly[i][0]
        ):
            c = not c
        j = i

    return c


def get_vals(geqdsk):

    points = list(zip(geqdsk["RBBBS"], geqdsk["ZBBBS"]))

    vals = np.zeros(geqdsk["PSIRZ"].shape)
    for ir, R in enumerate(geqdsk["AuxQuantities"]["R"]):
        for iz, Z in enumerate(geqdsk["AuxQuantities"]["Z"]):

            vals[iz, ir] = isPointInPath(R, Z, points)

    return vals.T.flatten()


treeloc = OMFIT["EFIT"]["FILES"]
keqdsk = treeloc["kEQDSK"]
meqdsk = treeloc["mEQDSK"]
geqdsk = treeloc["gEQDSK"]

R = geqdsk["AuxQuantities"]["R"]
Z = geqdsk["AuxQuantities"]["Z"]
RR, ZZ = np.meshgrid(R, Z)
dr = R[1] - R[0]
dz = Z[1] - Z[0]
darea = dr * dz

expmp2 = keqdsk["IN1"]["EXPMP2"]
brsp = meqdsk["ccbrsp"]["data"][0, :]
ecurrt = meqdsk["cecurr"]["data"][0, :]

cm = sum(rmp2fc * brsp[:, None], axis=0)
cm += sum(rmp2ec * ecurrt[:, None], axis=0)
# no vessel model for vcurrt
# cm += sum(gmp2vs*vcurrt[:,None],axis=0)

if update_boundary:
    vals = get_vals(geqdsk)

jtorR = geqdsk["AuxQuantities"]["Jt"].T.flatten()

cm += sum(rmp2pc * (jtorR * vals)[:, None], axis=0) * darea
error = expmp2 - cm

print(error[iprobe], cm[iprobe], expmp2[iprobe])
