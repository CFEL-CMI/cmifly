import numpy as np
from scipy import interpolate

def readDfield(filename):
    rawdata = open(filename,'r')
    # Read the file contents and generate a list with each line
    lines = rawdata.readlines()

    # initial positions of x,y,z grid
    line = lines[0][:-2] # remove EOL and one space
    linesplit = line.split(" ")
    xfs, yfs, zfs = map(float, linesplit)

    # step size for each dimension
    line = lines[1][:-2] # remove EOL and one space
    linesplit = line.split(" ")
    dfx, dfy, dfz = map(float, linesplit)

    # number of steps for each dimension
    line = lines[2][:-2] # remove EOL and one space
    linesplit = line.split(" ")
    xfn, yfn, zfn = map(int, linesplit)

    # create a new array for import field grid values
    field = np.zeros((xfn,yfn),float)

    # create x,y mesh of field grid
    xf = np.linspace(xfs,xfs+dfx*(xfn-1),xfn)
    yf = np.linspace(yfs,yfs+dfy*(yfn-1),yfn)

    for xfi in range(xfn):
        for yfi in range(yfn):
            line = lines[3 + yfi + xfi * yfn][:-1] # remove EOL
            field[xfi][yfi] = float(line)
    rawdata.close()
    return xf, yf, field

def readDgradient(filename):
    rawdata = open(filename,'r')
    # Read the file contents and generate a list with each line
    lines = rawdata.readlines()

    # initial positions of x,y,z grid
    line = lines[0][:-2] # remove EOL and one space
    linesplit = line.split(" ")
    xgs, ygs, zgs = map(float, linesplit)

    # step size for each dimension
    line = lines[1][:-2] # remove EOL and one space
    linesplit = line.split(" ")
    dxg, dyg, dzg = map(float, linesplit)

    # number of steps for each dimension
    line = lines[2][:-2] # remove EOL and one space
    linesplit = line.split(" ")
    xgn, ygn, zgn = map(int, linesplit)

    # create a new array for import field values
    gradx = np.zeros((xgn,ygn),float)
    grady = np.zeros((xgn,ygn),float)
    gradz = np.zeros((xgn,ygn),float)

    # create x,y mesh of field griddd
    xg = np.linspace(xgs,xgs+dxg*(xgn-1),xgn)
    yg = np.linspace(ygs,ygs+dyg*(ygn-1),ygn)

    for xgi in range(xgn):
        for ygi in range(ygn):
            line = lines[3 + ygi + xgi * ygn][:-2] # remove EOL and one space
            linesplit = line.split(" ")
            gradx[xgi][ygi], grady[xgi][ygi], gradz[xgi][ygi] = map(float, linesplit)
    return xg, yg, gradx, grady

def Dfieldintp(x,y,vs,xf,yf,field):
    finterp = interpolate.interp2d(xf,yf,field)
    return finterp(x,y)[0] * vs

def Dgradintp(x,y,vs,xg,yg,grad):
    gradinterp = interpolate.interp2d(xg,yg,grad)
    return gradinterp(x,y)[0] * vs
