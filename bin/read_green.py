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
grid = f.read_reals(f'{endian}f8'
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