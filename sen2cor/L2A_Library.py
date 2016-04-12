#!/usr/bin/env python

from numpy import *
from scipy.signal import medfilt
from scipy.ndimage import map_coordinates
from scipy import interpolate as sp
from scipy.ndimage.filters import uniform_filter

import sys, os


def stdoutWrite(s):
    sys.stdout.write(s)
    sys.stdout.flush()
    

def stderrWrite(s):
    sys.stderr.write(s)
    sys.stderr.flush()  

def statistics(arr, comment = ''):
    if len(arr) == 0:
        return False
    s = 'object:' + str(comment) + '\n'
    s += '--------------------' + '\n'
    s += 'shape: ' + str(shape(arr)) + '\n'
    s += 'sum  : ' + str(arr.sum()) + '\n'
    s += 'mean : ' + str(arr.mean()) + '\n'
    s += 'std  : ' + str(arr.std())  + '\n'
    s += 'min  : ' + str(arr.min())  + '\n'
    s += 'max  : ' + str(arr.max())  + '\n'
    s += '-------------------' + '\n'
    return s


def showImage(arr):
    from PIL import Image
    if(arr.ndim) != 2:
        sys.stderr.write('Must be a two dimensional array.\n')
        return False

    arrmin = arr.mean() - 3*arr.std()
    arrmax = arr.mean() + 3*arr.std()
    arrlen = arrmax-arrmin
    arr = clip(arr, arrmin, arrmax)
    scale = 255.0
    scaledArr = (arr-arrmin).astype(float32) / float32(arrlen) * scale
    arr = (scaledArr.astype(uint8))
    img = Image.fromarray(arr)
    img.show()
    return True


def reverse(a): return a[::-1]


def interplin(vin, xin, uin):
    """
     NAME:
        interplin()

     PURPOSE:
        Perform 1-d linear interpolation.  Values outside the bounds are
        permitted unlike the scipy.interpolate.interp1d module. The are
        extrapolated from the line between the 0,1 or n-2,n-1 entries.
        This program is not as powerful as interp1d but it does provide
        this which makes it compatible with the IDL interpol() function.

     CALLING SEQUENCE:
        yint = interplin(y, x, u)

     INPUTS:
        y, x:  The y and x values of the data.
        u: The x-values to which will be interpolated.

     REVISION HISTORY:
        Created: 2006-10-24, Erin Sheldon, NYU, from:
        http://sdss.physics.nyu.edu/esheldon/python/code/astro_code-2009-05-15/es_util.py
    """
    # Make sure inputs are arrays.  Copy only made if they are not.
    v=array(vin, ndmin=1, copy=False)
    x=array(xin, ndmin=1, copy=False)
    u=array(uin, ndmin=1, copy=False)
    # Find closest indices
    xm = x.searchsorted(u) - 1

    # searchsorted returns size(array) when the input is larger than xmax
    # Also, we need the index to be less than the last since we interpolate
    # *between* points.
    w, = where(xm >= (x.size-1))
    if w.size > 0:
        xm[w] = x.size-2

    w, = where(xm < 0)
    if w.size > 0:
        xm[w] = 0

    xmp1 = xm+1
    return (u-x[xm])*(v[xmp1] - v[xm])/(x[xmp1] - x[xm]) + v[xm]


def interpol(vin, xin, uin):
    return interplin(vin, xin, uin)


def linear_interpolation(x, grid):
    ind = arange(x.size,dtype=int)
    lin1 = sp.interpolate.interp1d(ind, x)
    res = lin1(grid)
    return res


def bilinear_interpolation(points, x, y):
    if len(points) == 4:
        aux = asarray(points)
        lp = reshape(aux, (2,2))
    else:
        lp = points

    arr = interpol2d(lp, x, y)
    res = asarray(asmatrix(arr).H)
    return res


def interpol1d(x, num):
    ind = arange(x.size,dtype=int)
    lin1 = sp.interpolate.interp1d(ind, x)
    grid = linspace(0, x.size - 1, num)
    res = lin1(grid)
    return float32(res)


def rectBivariateSpline(xIn, yIn, zIn):
    x = arange(zIn.shape[0], dtype=float32)
    y = arange(zIn.shape[1], dtype=float32)

    f = sp.RectBivariateSpline(x,y,zIn)
    return f(xIn,yIn)


def interpol2d(lp, x, y):
    f = sp.interpolate.interp2d(arange(lp.shape[1], dtype=float32),arange(lp.shape[0], dtype=float32),lp)
    lph = f(x,y)
    return float32(lph)


def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return array(map(pointwise, array(xs)))

    return ufunclike


def smooth(data, filt, edge_truncate=True):
    if edge_truncate == True:
        res = uniform_filter(data,size=filt, mode='nearest')
    else:
        res = uniform_filter(data,size=filt)
    return float32(res)


def congrid(a, newdims, method='linear', centre=False, minusone=False):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).

    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                                scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    if not a.dtype in [float32, float]:
        a = cast[float32](a)
    m1 = cast[int](minusone)
    ofs = cast[int](centre) * 0.5
    old = array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        sys.stderr.write('[congrid] dimensions error. ' \
                  'This routine currently only supports ' \
                  'rebinning to the same number of dimensions.\n')
        return None
    newdims = asarray( newdims, dtype=float32 )
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                                     * (base + ofs) - ofs )
        cd = array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return float32(newa)

    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                                     * (base + ofs) - ofs )
        # specify old dims
        olddims = [arange(i, dtype = float32) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = intp.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + range( ndims - 1 )
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )
            mint = intp.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )

        if ndims > 1:
        # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return float32(newa)
    elif method in ['spline']:
        #oslices = [ slice(0,j) for j in old ]
        #oldcoords = ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = mgrid[nslices]

        newcoords_dims = range(rank(newcoords))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = map_coordinates(a, newcoords, order=1, mode='nearest')
        return float32(newa)
    else:
        sys.stderr.write('Congrid error: Unrecognized interpolation type.\n'\
                  'Currently only \'neighbour\', \'nearest\',\'linear\'\n'\
                  'and \'spline\' are supported.\n')
        return None

    return


def median_filter_2d(arr, width):
    if width == 1:
        medarray = arr
    else:
        medarray = medfilt(arr, kernel_size = width)
        #
        # medfilt doesn't handle edges the same way as IDL, so we
        # have to fix them up
        #
        istart = (width-1)/2
        iend = arr[0].size - (width+1)/2
        i = arange(arr[0].size)
        w = (i < istart) |  (i > iend)
        arr1d = ravel(arr[0])
        med1d = ravel(medarray[0])
        med1d[w] = arr1d[w]
        medarray[0] = med1d
        #
        iend = arr[1].size - (width+1)/2
        i = arange(arr[1].size)
        w = (i < istart) |  (i > iend)
        arr1d = ravel(arr[1])
        med1d = ravel(medarray[1])
        med1d[w] = arr1d[w]
        medarray[1] = med1d

    return float32(medarray)


def regress(x,y):
    A = array([x,ones(x.size)])
    w = linalg.lstsq(A.T,y)[0]
    return array([w[1],w[0]])


def fit_coeff(x, y):
# Input:
#     x = wv grid (usually 3 grid points, scaled wv values, e.g. [1000, 2000, 2900] )
#     y = radiative transfer function (transmittance, Lpath, Edir etc)
#
#    Function type: y  = exp( a + b*sqrt(x) )
#    Convert into linear regression:  log(y) = a + b*sqrt(x)
#
# Remarks:
#     The inverse function to "fit_coeff" is provided by function "polx"
#     an NaN usually happens in the 1.37-1.44 or 1.80-1.94 micron region of the LUT conversion
#     The numpy warning is disabled, the nan is replaced by a nan_to_num

    seterr(invalid='ignore', divide='ignore')
    xr = sqrt(x)  # make 2-D array for function regress
    yrr = log(y)
    yr = nan_to_num(yrr)
    return regress(xr, yr)


def rebin(a, *args):
    '''rebin ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)
    '''
    shape = a.shape
    lenShape = len(shape)
    factor = asarray(shape)/asarray(args)
    evList = ['a.reshape('] + \
                ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
                [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)] + \
                ['/factor[%d]'%i for i in range(lenShape)]

    return eval(''.join(evList))


#--------------------------------------------------------------------------------
def triang_interpol(visindex):
# triangular (x,y) interpolation with delny, not used for sentinel-2
    '''
    list_valid = nonzero(visindex > 0)
    list_gaps = nonzero(visindex == 0)

    n_valid = sum(list_valid)
    n_gaps = sum(list_gaps)

    if (n_valid < 3):
        sys.stderr.write('Not enough reference pixels for triangular interpolation.\n')
        return 1
    if (n_gaps == 0):
        return 0

    nx = array(visindex[0,:].size)
    ny = array(visindex[:,0].size)

    xcoor = array(list_valid % nx, copy=0).astype(int32)
    ycoor = array(list_valid / nx, copy=0).astype(int32)

    # calculate triangles

    triangles = delny([xcoor,ycoor])
    z = visindex[list_valid]

    visav = uint8(sum(visindex[list_valid] / n_valid))

    # IDL2PY_TBD find treegrid equivalent: interp_visindex = trigrid(xcoor, ycoor, z, triangles, array([1, 1]), array([0, 0, nx - 1, ny - 1]), missing=visav)
    # free memory
    interp_visindex = 0
    return uint8(interp_visindex + 0.5)
    '''
    return False

def unitmvector(matrixm, length=None):
    #;if matrixM is integer, then it is converted to
    #;float to operate with enough accuracy:
    if(matrixm.dtype != float32):
        matmf = array(matrixm, copy=0).astype(float32)
    else:
        matmf = matrixm

    length = sqrt(sum(matmf ** 2,2))

    for i in arange(0, 3):
        matmf[:,:,i] = matmf[:,:,i] / length

    return matmf

def jav_cellgradient(dem, dl):

    m = size(dem[0,:]) # number of columns
    n = size(dem[:,0]) # number of rows

    #******** COMPUTE NORMAL VECTOR    ************
    cellgr = zeros([n, m, 3], float32) * 0.0

    zr = zeros([n, m], float32)
    zd = zeros([n, m], float32)
    zrd = zeros([n, m], float32)

    zr[:,0:(m - 2)+1] = dem[:,1:(m - 1)+1]
    zd[0:(n - 2)+1,:] = dem[1:(n - 1)+1,:]
    zrd[0:(n - 2)+1,0:(m - 2)+1] = dem[1:(n - 1)+1,1:(m - 1)+1]

    tmp1 = dem + zd - zr - zrd
    tmp2 = dem - zd + zr - zrd

    # cellgr(0,*,*) = .5*dl*(DEM+Zd-Zr-Zrd)
    cellgr[:,:,0] = .5 * dl * tmp1
    # cellgr[1,*,*] = .5*dl*(DEM-Zd+Zr-Zrd)
    cellgr[:,:,1] = .5 * dl * tmp2
    cellgr[:,:,2] = (dl ** 2)

    cellgr[0:(n - 1)+1,m - 1,:] = cellgr[0:(n - 1)+1,m - 2,:]
    cellgr[n - 1,0:(m - 1)+1,:] = cellgr[n - 2,0:(m - 1)+1,:]
    cellgr[n - 1,m - 1,:] = cellgr[n - 2,m - 2,:]

    # norm_matvector() returns the matrix of unit vectors and the length of
    # every vector  (area, in this case)

    return unitmvector(cellgr)

def atc_skyview(dem, azimuthinterval, elevationinterval, dl, unders=None, group=None):

    rows0 = size(dem[:,0]) # number of rows
    cols0 = size(dem[0,:]) # number of columns
    radeg = 180.0 / pi

    if (unders is not None):
        unders = unders[0]
        dl = dl * unders
        cols = cols0 / unders
        rows = rows0 / unders
        if unders > 1:
            dem = congrid(dem, cols, rows, interp=True, center=True)
    else:
        cols = cols0
        rows = rows0
        unders = 1

    skyviewcell = zeros([rows, cols], float32) # integer appr
    sh = zeros([rows, cols])
    skyviewfazstep = zeros([rows, cols], float32)

    gr = jav_cellgradient(dem, dl)
    maxslope = arccos(array(gr[:,:,2], copy=0).min()) * radeg

    for azim in arange(0, (360 - azimuthinterval)+(azimuthinterval), azimuthinterval):

        skyviewfazstep[:,:] = maxslope
        for elev in arange(maxslope, (0)+(-elevationinterval), -elevationinterval):
            sun = array([azim, elev])
            sh[:,:] = doshadeskv(dem, sun, dl, gr)
            sunny = nonzero(ravel(sh == 1))
            if(sunny != 0): #IDL2PY_TBD: check this !!!
                skyviewfazstep[sunny] = elev

        skyviewcell = skyviewcell + ((cos(skyviewfazstep / radeg)) ** 2) * azimuthinterval / 360.0

    #; go back to original resolution
    if (unders is not None):
        if unders != 1:    #IDL2PY: congrid - attention !!!
            skyviewcell = congrid(skyviewcell, cols0, rows0, interp=True, center=True)

    return uint8(skyviewcell * 100)             # returns percentage


def jav_doshade(dem, sun, dl, shading=None, smfact=None):

    dl = array(dl, copy=0).astype(float32)

    #; define smoothing factor:
    if bitwise_not((smfact is not None)):
        smfact = 3  #; default smoothing!

    sunvector = sun / sqrt((sun ** 2).sum())

    if (size(sunvector, dimensions=True))[0] == 2:     # for the 2D case
        azimuth = sun[0] * (pi / 180.0)
        elevation = sun[1] * (pi / 180.0)
        sunvector = zeros([3], float32)

        sunvector[0] = sin(azimuth) * cos(elevation) #positive Eastwards
        sunvector[1] = -cos(azimuth) * cos(elevation) #positive Southwards
        sunvector[2] = sin(elevation)

    # If operated with a different reference system, eg: Y is positive northwards,
    #    simply make SunVector[1] = -SunVector[1]

    # solvector is opposite to sunvector and used for integer steps scanning
    #    along the solar direction
    # this makes the largest coordinate of solvector =  1

    solvector = -sunvector / max(abs(sunvector[0:2]))

    sunvectorpxy = sunvector         # horizontal pojection sunvector, step1
    sunvectorpxy[2] = 0.0e0            # horizontal pojection sunvector, step2

    # vector normal to sun in the horizontal plane
    sunvectorh = cross(sunvectorpxy, sunvector)
    sunvectorh = sunvectorh / sqrt(sum(sunvectorh ** 2)) #** Unit SunVectorH

    # unit vector normal to sun upwards
    sunvectoru = cross(sunvectorh, sunvector)

    # direction of view of sun
    # sunvectorn = -sunvector

    sz = size(dem, dimensions=True)
    cols = array(sz[0], copy=0).astype(int32)                  # number of columns
    rows = array(sz[1], copy=0).astype(int32)                  # number of rows

    #; Vector to Origin of coordinates,
    #; for projection onto plane perpendicular to sun
    #; It doesn't need to be (0,0,0) any arbitrary point is fine (0,0,0)
    #; is simpler
    if (dl % 1 != 0):  # suppixel resolution...
        xvalues = array(transpose(dot(transpose(arange(cols)), transpose((1)*ones([rows])))) * dl, copy=0).astype(float32)
        yvalues = array(transpose(dot(transpose((1)*ones([cols])), transpose(arange(rows)))) * dl, copy=0).astype(float32)
    else:
        xvalues = transpose(dot(transpose(arange(cols)), transpose((1)*ones([rows])))) * array(dl, copy=0).astype(int32)
        yvalues = transpose(dot(transpose((1)*ones([cols])), transpose(arange(rows)))) * array(dl, copy=0).astype(int32)
        zvalues = dem

    #VectortoOrigin[i] is [xvalues[i],yvalues[i],zvalues[i]
    #dot product [x,y,z] * SunVectorU

    zprojection = xvalues[:,:] * sunvectoru[0] + yvalues[:,:] * sunvectoru[1] + zvalues[:,:] * sunvectoru[2]

    # Determine origin of scanning lines in the direction of the sun
    _expr = 1
    if _expr == ((sunvector[0] <= 0.0)):
        xend = cols - 1 # sun is on the West
    elif _expr == ((sunvector[0] > 0.0)):
        xend = 0 # sun is on the East
    else:
        pass

    _expr = 1
    if _expr == ((sunvector[1] <= 0.0)):
        yend = rows - 1    # sun is on the North
    elif _expr == ((sunvector[1] > 0.0)):
        yend = 0 # sun is on the South
    else:
        pass


    # sombra is the array to store binary shadow values
    sombra = (1)*ones([rows, cols])

    #------------- SCANNING CASE SUN at  0, 90, 180, 270 degrees -----------

    # 0 or 180 degrees
    if (sunvector[0] == 0):
        lastloop = cols - 1
        lengthsc = rows - 1
        firsti = 0
        firstj = rows - yend - 1
        scanindex = concatenate([0, 1])
        scanvector = zeros(2, lengthsc)
        scanvector[:,1] = firstj + arange(lengthsc) * solvector[1]
        # IDL2PY_TBD goto, JumpNoZero

    # 90 or 270 degrees
    if (sunvector[1] == 0):
        lastloop = rows - 1
        lengthsc = cols - 1
        firsti = xend - cols - 1
        firstj = 0
        scanindex = concatenate([1, 0])
        scanvector = zeros(2, lengthsc)
        scanvector[:,0] = firsti + arange(lengthsc) * solvector[0]
        # IDL2PY_TBD goto, JumpNoZero

    #-------------------------  SCANNING BLOCK I ----------------------------------

    #***************** scanning along X-AXIS  ****************

    for i in arange(0, (cols)):
    # length of scanning line (vector)
        lengthx1 = abs(xend - i) / solvector[0] # solve for: l*Solvector=Xmax
        lengthx2 = (rows - 1) / solvector[1] # solve for: l*Sv=Ymax
        lengthsc = array(array(abs(concatenate([lengthx1, lengthx2])), copy=0).min(), copy=0).astype(int32)

        if (lengthsc == 0):
            continue # this means a corner = sunny

        vectorlength = arange(lengthsc)

        #X-origin of scan line
        firsti = i
        #* Y-origin of scan line
        firstj = rows - yend - 1

        # ScanVectorX: arrays of xvalues along scan line  ScanVectorY: Yvalues
        scanvectorx = i + vectorlength * solvector[0]

        # ScanVectorY: Yvalues
        scanvectory = firstj + vectorlength * solvector[1]

        # Array of Z_Projections along scan line
        #IDL2PY_TBD: check interp2d:
        zpdxdy = intp.interp2d(scanvectorx, scanvectory, zprojection)

        #; Some antialiasing device is required.  A temporary improvement is
        #;achieved by averaging the cell projections with the previous value
        #;along the scan line

        #IF (lengthsc GT 2) THEN $
        #     ZPdxdy[1:lengthsc-1] = 0.5*(ZPdxdy[0:lengthsc-2] + $
        #                                          ZPdxdy[1:lengthsc-1] )

        if (lengthsc > smfact):
            zpdxdy = smooth(zpdxdy, smfact, edge_truncate=True)
        #; initial zpcompare smaller than any possible zprojection
        zpcompare = -1e30

        for n in arange(0, (lengthsc)):
            if (zpdxdy[n] < zpcompare):
                sombra[scanvectory[n],scanvectorx[n]] = 0
            else:
                zpcompare = zpdxdy[n]

    xend = abs(xend - i) - 1

    #==============================================================================
    #-------------------------  SCANNING BLOCK J ----------------------------------

    #*** scanning    Along Y-AXIS
    lenarray = zeros([rows], float32)
    for j in arange(0, (rows)):
    #length of scanning line (vector)
        lengthy1 = abs(yend - j) / solvector[1] # solve for: l*Sv=Ymax
        lengthy2 = (cols - 1) / solvector[0] # solve for: l*Sv=Xmax
        lengthsc = array(array(abs(concatenate([lengthy1, lengthy2])), copy=0).min(), copy=0).astype(int32)

        if (lengthsc == 0):
            continue # this means a corner = sunny

        vectorlength = arange(lengthsc)
        lenarray[j] = lengthsc
        # X-origin of scan line
        firsti = xend
        # Y-origin of scan line
        firstj = j

        #; ScanVectorX: arrays of xvalues along scan line
        scanvectorx = xend + vectorlength * solvector[0]

        #; ScanVectorY: Yvalues
        scanvectory = firstj + vectorlength * solvector[1]

        #;Array of Z_Projections along scan line
        #ZPdxdy = zprojection[ScanVectorX,ScanVectorY]
        #IDL2PY_TBD: check interp2d:
        zpdxdy = intp.interp2d(scanvectorx, scanvectory, zprojection)

        #; Some antialiasing device is required.  A temporary improvement is
        #;achieved by averaging the cell projections with the previous value
        #;along the scan line

        if (lengthsc > smfact):
            zpdxdy = smooth(zpdxdy, smfact, edge_truncate=True)

        #IF (lengthsc GT 2) THEN $
        #     ZPdxdy[1:lengthsc-1] = 0.5*(ZPdxdy[0:lengthsc-2] + ZPdxdy[1:lengthsc-1] )

        #; initial zpcompare smaller than any possible zprojection
        zpcompare = -1e30

        for n in arange(0, (lengthsc)):
            if (zpdxdy[n] < zpcompare):
                sombra[scanvectory[n],scanvectorx[n]] = 0
            else:
                zpcompare = zpdxdy[n]


    #==============================================================================

    # IDL2PY_TBD goto, JumpZero

    #-------------------------  SCANNING BLOCK ZERO -----------------------------

    # jumpnozero:

    for i in arange(0, lastloop+1):

        scanvector[:,scanindex[0]] = i
        #; Array of Z_Projections along scan line
        #IDL2PY_TBD: check interp2d:
        zpdxdy = intp.interp2d(scanvector[:,scanindex[0]], scanvector[:,scanindex[1]], zprojection)

        #; Some antialiasing device is required.  A temporary improvement is
        #;achieved by averaging the cell projections with the previous value
        #;along the scan line

        if (lengthsc > smfact):
            zpdxdy = smooth(zpdxdy, smfact, edge_truncate=True)
        #ZPdxdy[1:lengthsc-1] = 0.5*(ZPdxdy[0:lengthsc-2] + ZPdxdy[1:lengthsc-1] )

        #; initial zpcompare smaller than any possible zprojection
        zpcompare = -1e30

        for n in arange(0, (lengthsc)):
            if (zpdxdy[n] < zpcompare):
                sombra[scanvector[n,scanindex[1]],scanvector[n,scanindex[0]]] = 0
            else:
                zpcompare = zpdxdy[n]

    #==============================================================================

    # jumpzero:


    # border columns and rows are set to 1st (or 2nd) inner column and row
    #;sombra[1,*] = sombra[2,*]
    sombra[:,0] = sombra[:,1]
    #;sombra[cols-2,*] = sombra[cols-3,*]

    sombra[:,cols - 1] = sombra[:,cols - 2]
    #;sombra[*,1] = sombra[*,2]

    sombra[0,:] = sombra[1,:]
    #;sombra[*,rows-2] = sombra[*,rows-3]

    sombra[rows - 1,:] = sombra[rows - 2,:]

    if (shading):
        gr = jav_cellgradient(dem, array(dl, copy=0).astype(float32))
        shading = zeros([rows, cols], float32)
        shading[:,:] = gr[:,:,0] * sunvector[0] + gr[:,:,1] * sunvector[1] + gr[:,:,2] * sunvector[2]
        shading = maximum(shading, 0.)

    return sombra


def doshadeskv(dem, sun, dl, gr, group=None):

    azimuth = (sun[0] * (pi / 180.0e0))
    elevation = (sun[1] * (pi / 180.0e0))
    sunvector = zeros([3], float32)    # RR (instead of dblarr)
    sunvector[0] = sin(azimuth) * cos(elevation) #positive Eastwards
    sunvector[1] = -cos(azimuth) * cos(elevation) #positive Southwards
    sunvector[2] = sin(elevation)

    solvector = -sunvector / max(abs(sunvector[0:2]))

    sunvectorpxy = sunvector         # horizontal pojection sunvector, step1
    sunvectorpxy[2] = 0.0 # RR, not 0.0D0 ; horizontal pojection sunvector, step2

    # vector normal to sun in the horizontal plane
    sunvectorh = cross(sunvectorpxy, sunvector)
    sunvectorh = sunvectorh / sqrt(sum(sunvectorh ** 2)) # Unit SunVectorH

    # unit vector normal to sun upwards
    sunvectoru = cross(sunvectorh, sunvector)

    # direction of view of sun
    #sunvectorn = -sunvector
    rows = size(dem[:,0]) # number of rows
    cols = size(dem[0,:]) # number of columns

    #Vector to Origin of coordinates, for projection onto
    #plane perpendicular to sun
    # It doesn't need to be (0,0,0)
    #any arbitrary point is fine (0,0,0) is simpler

    x1 = arange(cols)
    x2 = (ones(cols)).astype(int32)
    y1 = arange(rows)
    y2 = (ones(rows)).astype(int32)

    xvalues = outer(y2,x1)
    yvalues = outer(y1,x2)

    #VectortoOrigin[i] is [xvalues[i],yvalues[i],zvalues[i]
    #dot product [x,y,z] * SunVectorU

    z1 = xvalues * dl * sunvectoru[0]
    z2 = yvalues * dl * sunvectoru[1]
    z3 = dem * sunvectoru[2]

    zprojection = z1 + z2 + z3

    # Determine origin of scanning lines in the direction of the sun
    _expr = 1
    if _expr == ((sunvector[0] <= 0.0)):
        xend = cols - 1 # sun is on the West
    elif _expr == ((sunvector[0] > 0.0)):
        xend = 0 # sun is on the East
    else:
        pass

    _expr = 1
    if _expr == ((sunvector[1] <= 0.0)):
        yend = rows - 1    # sun is on the North
    elif _expr == ((sunvector[1] > 0.0)):
        yend = 0 # sun is on the South
    else:
        pass

    # sombra is the array to store binary shadow values
    sombra = (1)*ones([rows, cols])

    #------------- SCANNING CASE SUN at  0, 90, 180, 270 degrees -----------

    jump = False

    # 0 or 180 degrees
    if (sunvector[0] == 0):
        lastloop = cols - 1
        lengthsc = rows - 1
        firsti = 0
        firstj = rows - yend - 1
        scanindex = array([0, 1]).astype(int32)
        scanvector = zeros([lengthsc,2]).astype(int32)
        scanvector[:,1] = firstj + arange(lengthsc) * solvector[1]
        jump = True

    # 90 or 270 degrees
    if (sunvector[1] == 0):
        lastloop = rows - 1
        lengthsc = cols - 1
        firsti = xend - cols - 1
        firstj = 0
        scanindex = array([1, 0]).astype(int32)
        scanvector = zeros([lengthsc,2]).astype(int32)
        scanvector[:,0] = firsti + arange(lengthsc) * solvector[0]
        jump = True

    if(jump == False):
        #-------------------------  SCANNING BLOCK I ----------------------------------

        #***************** scanning along X-AXIS  ****************

        for i in arange(0, (cols)):
            # length of scanning line (vector)
            lengthx1 = abs(xend - i) / solvector[0] # solve for: l*Solvector=Xmax
            lengthx2 = (rows - 1) / solvector[1] # solve for: l*Sv=Ymax
            lengthsc = int(min(absolute(array([lengthx1, lengthx2]))))

            if (lengthsc == 0):
                continue # this means a corner = sunny

            vectorlength = zeros(lengthsc)

            # X-origin of scan line
            firsti = i
            # Y-origin of scan line
            firstj = rows - yend - 1

            # ScanVectorX: arrays of xvalues along scan line  ScanVectorY: Yvalues
            scanvectorx = array(i + vectorlength * solvector[0], copy=0).astype(int32)

            # ScanVectorY: Yvalues
            scanvectory = array(firstj + vectorlength * solvector[1], copy=0).astype(int32)

            # Array of Z_Projections along scan line
            zpdxdy = zprojection[scanvectory,scanvectorx]

            # Some antialiasing device is required.  A temporary improvement is
            # achieved by averaging the cell projections with the previous value
            # along the scan line

            if (lengthsc > 2):
                zpdxdy[1:(lengthsc - 1)+1] = 0.5 * (zpdxdy[0:(lengthsc - 2)+1] + zpdxdy[1:(lengthsc - 1)+1])

            #; initial zpcompare smaller than any possible zprojection
            zpcompare = -1e30

            for n in arange(0, (lengthsc)):
                if (zpdxdy[n] < zpcompare):
                    sombra[scanvectory[n],scanvectorx[n]] = 0
                else:
                    zpcompare = zpdxdy[n]

            xend = abs(xend - i) - 1

        #==============================================================================
        #-------------------------  SCANNING BLOCK J ----------------------------------

        # scanning    Along Y-AXIS
        lenarray = zeros([rows], float32)
        for j in arange(0, (rows)):
        #length of scanning line (vector)
            lengthy1 = abs(yend - j) / solvector[1] # solve for: l*Sv=Ymax
            lengthy2 = (cols - 1) / solvector[0] # solve for: l*Sv=Xmax
            lengthsc = int(min(absolute(array([lengthy1, lengthy2]))))

            if (lengthsc == 0):
                continue # this means a corner = sunny

            vectorlength = arange(lengthsc)
            lenarray[j] = lengthsc
            #; X-origin of scan line
            firsti = xend
            #; Y-origin of scan line
            firstj = j

            #; ScanVectorX: arrays of xvalues along scan line
            scanvectorx = array(xend + vectorlength * solvector[0], copy=0).astype(int32)

            #; ScanVectorY: Yvalues
            scanvectory = array(firstj + vectorlength * solvector[1], copy=0).astype(int32)

            #Array of Z_Projections along scan line
            zpdxdy = zprojection[scanvectory,scanvectorx]

            # Some antialiasing device is required.  A temporary improvement is
            #;achieved by averaging the cell projections with the previous value
            #;along the scan line

            if (lengthsc > 2):
                zpdxdy[1:(lengthsc - 1)+1] = 0.5 * (zpdxdy[0:(lengthsc - 2)+1] + zpdxdy[1:(lengthsc - 1)+1])

            # initial zpcompare smaller than any possible zprojection
            zpcompare = -1e30

            for n in arange(0, (lengthsc)):
                if (zpdxdy[n] < zpcompare):
                    sombra[scanvectory[n],scanvectorx[n]] = 0
                else:
                    zpcompare = zpdxdy[n]


    elif(jump == True):
    #-------------------------  SCANNING BLOCK ZERO -----------------------------
        for i in arange(0, lastloop+1):

            scanvector[:,scanindex[0]] = i
            # Array of Z_Projections along scan line
            zpdxdy = zprojection[scanvector[:,scanindex[1]],scanvector[:,scanindex[0]]]

            # Some antialiasing device is required.  A temporary improvement is
            # achieved by averaging the cell projections with the previous value
            # along the scan line

            zpdxdy[1:(lengthsc - 1)+1] = 0.5 * (zpdxdy[0:(lengthsc - 2)+1] + zpdxdy[1:(lengthsc - 1)+1])

            # initial zpcompare smaller than any possible zprojection
            zpcompare = -1e30

            for n in arange(0, (lengthsc)):
                if (zpdxdy[n] < zpcompare):
                    sombra[scanvector[n,scanindex[1]],scanvector[n,scanindex[0]]] = 0
                else:
                    zpcompare = zpdxdy[n]

    # border columns and rows are set to 1st (or 2nd) inner column and row
    sombra[:,0] = sombra[:,1]
    sombra[:,cols - 1] = sombra[:,cols - 2]
    sombra[0,:] = sombra[1,:]
    sombra[rows - 1,:] = sombra[rows - 2,:]

    #set to zero self-shading cells
    shading = zeros([rows, cols], float32)
    shading[:,:] = gr[:,:,0] * sunvector[0] + gr[:,:,1] * sunvector[1] + gr[:,:,2] * sunvector[2]
    sombra = sombra * (shading > 0)

    return sombra

#------------------------------------------------------------------------------------
def clear_line(datyp, bandm1x, bandm3x, li_clear):
# IDL2PY_TBD; IDL2PY_TBD @onerror.inc

# calculation of "clear_line" correlation between RED and BLUE band
# (for haze HOT transform)
#
# Input:    bandm1x = blue band
#             bandm3x = red band
#             li_clear = list of clear elements in image
#             datyp     = ENVI data type
# Output :
#             a_cl, b_cl linear regression y = a_cl + b_cl*x
#             corrcoeff = correlation coefficient

    siz = size(bandm1x)
    ncol = siz[1]
    nrow = siz[2]

    _expr = (datyp)
    if _expr == 1:
        clear_band1 = zeros(ncol, nrow)
        clear_band3 = zeros(ncol, nrow)
    elif _expr == 2:
        clear_band1 = zeros(ncol, nrow)
        clear_band3 = zeros(ncol, nrow)
    elif _expr == 12:
        clear_band1 = zeros(ncol, nrow)
        clear_band3 = zeros(ncol, nrow)
    elif _expr == 3:
        clear_band1 = zeros(ncol, nrow)
        clear_band3 = zeros(ncol, nrow)
    elif _expr == 4:
        clear_band1 = zeros([nrow, ncol],float32)
        clear_band3 = zeros([nrow, ncol],float32)
    else:
        raise RuntimeError('no match found for expression')

    clear_band1[li_clear] = bandm1x[li_clear]

    clear_band3[li_clear] = bandm3x[li_clear]

    vclear1 = clear_band1[li_clear]
    vclear3 = clear_band3[li_clear]

    #weights = (1.0)*ones(li_clear)
    #IDL2PY_TBD: http://www.scipy.org/Cookbook/LinearRegression
    return regress(vclear1, vclear3)


#--------------------------------------------------------------------------------
def check_if_finite(yfit, y0):
# IDL2PY_TBD; IDL2PY_TBD @onerror.inc

# Input:  yfit=fltarr(2) from regression
#            y0  = y for wv=0.4 cm (scalar, e.g. y=tdir, or y=e0t etc)
# Output: yfit is not altered if numeric data
#            yfit[0] = log(y0) if yfit contains NaN
#            yfit[1] =  1.0E-6                 "         "      "
#  Since  y(wv) = exp(yfit[0] + yfit[1]*sqrt(wv))  the NaN case provides
#            almost constant values of y independent of wv and y(wv=400) = y0
#            (wv is scaled water vapor [cm*1000])
# Called from "apdat_lut" routines

    lis = isfinite(yfit) # find equivalent for finite ...
    lix = where(ravel(lis == 0))[0]
    if (lix.size > 0):
        yfit[0] = log(y0)
        yfit[1] = 1.0e-6
    return yfit


#--------------------------------------------------------------------------------
def polx1(uu, coeff):
# IDL2PY_TBD; IDL2PY_TBD @onerror.inc

# Input: uu, coeff
#     uu is scaled water vapor (scalar, or vector)
#     coeff is 2 element vector
# Output: array (or vector, or scalar) of exp. function below
# This function is called from routines "apda_lut", "wv_retrieval", and coeff is already
# the appropriate 2-element vector adapted to the current of the 3 wv regions.
        return exp(coeff[0] + coeff[1] * sqrt(uu))      # exponential fit


#--------------------------------------------------------------------------------
def polx(uu1_interp, uu, coeff):
# IDL2PY_TBD; IDL2PY_TBD @onerror.inc

# Input:
#     uu1_interp : 4-, 5- or 6-element vector of scaled wv for the current elevation
#                      i.e., [400, 1000, 2000, 2900] or
#                              [400, 1000, 2000, 2900, 4000] or [400, 1000, 2000, 2900, 4000, 5000]
#     uu            : is scaled water vapor (scalar, or vector, or 2-D matrix)
#     coeff(n)    : n=2*2=4 element vector to fit 2 regions [400, 1000, 2000] and [ 1000, 2000, 2900]
#                    : n=2*3=6     "                        3    "      [400 - 4000]
#                    : n=2*4=8     "                        4    "      [400 - 5000]  (see apda1_lut_constvis)
# Output: fct = exp( a + b*sqrt(uu))    (where a=coeff[0], b=coeff[1] )
#            fct has same size as uu  (i.e. 2-matrix, vector, or scalar)
# Similar to function polx1, but coeff is the full n-element fit vector

    li1 = where(ravel(uu <= uu1_interp[1]))[0]
    if (uu1_interp.size == 4):
        li2 = where(ravel(uu > uu1_interp[1]))[0]
    else:
        li2 = where(ravel(bitwise_and(uu > uu1_interp[1], uu <= uu1_interp[3])))[0]
    if (uu1_interp.size >= 5):
        if (uu1_interp.size == 5):
            li3 = where(ravel(uu > uu1_interp[3]))[0]
        else:
            li3 = where(ravel(bitwise_and(uu > uu1_interp[3], uu <= uu1_interp[4])))[0]
        if (uu1_interp.size >= 6):
            li4 = where(ravel(uu > uu1_interp[4]))[0]

    siz = size(uu)
    if (siz[0] == 0):
        fct = zeros([1], float32)
    else:
        if (siz[0] == 1):
            fct = zeros([siz[1]], float32)
        else:
            fct = zeros([siz[2], siz[1]], float32)

    if (len(li1) > 0):
        fct[li1] = exp(coeff[0] + coeff[1] * sqrt(uu[li1]))  # coeff[0:1] for region 1
    if (len(li2) > 0):
        fct[li2] = exp(coeff[2] + coeff[3] * sqrt(uu[li2]))  # coeff[2:3]  "     "     2
    if (len(li3) > 0):
        fct[li3] = exp(coeff[4] + coeff[5] * sqrt(uu[li3]))
    if (len(li4) > 0):
        fct[li4] = exp(coeff[6] + coeff[7] * sqrt(uu[li4]))

    if (siz[0] == 0):
        fct = fct[0]

    return fct


#---------------------------------------------------------------------------------
def polx2(uu1_interp, uu, coeff):
# IDL2PY_TBD; IDL2PY_TBD @onerror.inc

# Input:
#     uu1_interp : 4-, 5- or 6-element vector of scaled wv for the current elevation
#                      (similar to polx routine)
#     uu            : is scaled water vapor (vector, or 2-D matrix)
#     coeff(n,*) or coeff(n,*,*)  where n=4, 6, or 8 to fit the 2, 3, 4 wv regions
#                      with two coefficients per region  (see apda1_lut_constvis)
# Output: fct = exp( a + b*sqrt(uu))    (where a=coeff[0], b=coeff[1] )
#            fct has same size as uu  (i.e. 2D-matrix or vector)
# Similar to polx, but with 2- or 3-D array coeff

    li1 = where(ravel(uu <= uu1_interp[1]))[0]
    if (uu1_interp.size == 4):
        li2 = where(ravel(uu > uu1_interp[1]))[0]
    else:
        li2 = where(ravel(bitwise_and(uu > uu1_interp[1], uu <= uu1_interp[3])))[0]
    if (uu1_interp.size >= 5):
        if (uu1_interp.size == 5):
            li3 = where(ravel(uu > uu1_interp[3]))[0]
        else:
            li3 = where(ravel(bitwise_and(uu > uu1_interp[3], uu <= uu1_interp[4])))[0]
        if (uu1_interp.size >= 6):
            li4 = where(ravel(uu > uu1_interp[4]))[0]

    siz = size(uu)
    if (siz[0] == 1):
        fct = zeros([siz[1]], float32)
    else:
        fct = zeros([siz[2], siz[1]], float32)

    if (len(li1) > 0):
        fct[li1] = exp(coeff[li1,0] + coeff[li1,1] * sqrt(uu[li1]))
    if (len(li2) > 0):
        fct[li2] = exp(coeff[li2,2] + coeff[li2,3] * sqrt(uu[li2]))
    if (len(li3) > 0):
        fct[li3] = exp(coeff[li3,4] + coeff[li3,5] * sqrt(uu[li3]))
    if (len(li4) > 0):
        fct[li4] = exp(coeff[6] + coeff[7] * sqrt(uu[li4]))

    return fct


#---------------------------------------------------------------------------------
def polx3(ju, uu, coeff):
# IDL2PY_TBD; IDL2PY_TBD @onerror.inc

# Input:
#     ju = start index for wv region based on scene average wv
#            coeff[ju:ju+1,*] is used (ju=0, 2, 4, or 6)
#            ju=0 : scaled wv region with grid points ( 400, 1000, 2000) is used
#            ju=2 :                                              (1000, 2000, 2900)
#            ju=4 :                                              (2000, 2900, 4000)
#            ju=6 :                                              (2900, 4000, 5000)
#     uu = scaled water vapor (vector or 2-D matrix)
#     coeff = fltarr(n,*) or fltarr(n,*,*)
# Output :
#     fct = exp( a + b*sqrt(uu))    (where a=coeff[0], b=coeff[1] )
# called from atcor3

    if(uu.ndim == 1):
        #fct = zeros([uu[0]], float32)
        fct = exp(coeff[:,ju] + coeff[:,ju + 1] * sqrt(uu))
    else:
        #fct = zeros([uu[0], uu[1]], float32)
        fct = exp(coeff[:,:,ju] + coeff[:,:,ju + 1] * sqrt(uu))

    return fct


#---------------------------------------------------------------------------------
def slice_histogram(slice_levels, hist_phi_unscl):

# NAME:
#    slice_histogram  (category: de-shadowing)
#
# PURPOSE:
#     slice histogram at different levels to get the number of intersection points
#     for each level (part of cloud shadow algorithm).
#
#     Case 1: high slice levels: 0.5 to 0.8, histogram hist_phi_unscl starts with
#                the main peak within the first 3 histogram points: skewed "strange" histogram
#                caused by (too) many shadow pixels representing the peak.
#                Processing continues, but results may have a poor quality.
#     Case 2: high slice levels: 0.5 to 0.8
#                If the number of intersection points "left" of the main histogram peak is greater than 3
#                then the histogram has large amplitude fluctuations, and it must to be smoothed
#                with a greater bin size.
#     Case 3: low slice levels: 0.05 to 0.40: the second smaller peak of shadow pixels is
#                searched for. In this case the is1, list_score_intersection are passed to the routine
#                threshold_shadow_histo to calculate the histogram threshold for the shadow pixels.
#
# CALLING SEQUENCE:
#     slice_histogram, slice_levels, list_score_intersection, hist_phi_unscl, is1, np_slice
#
# PARAMETERS:
#    Input:
#        slice_levels of normalized histogram, e.g. slice_levels=[0.8. 0.7, 0.6, 0.5]
#                If n_elements(slice_levels) = 4 then slice_histogram is called the first time
#                If n_elements(slice_levels) =10 then slice_histogram is called the second time
#                    with a finer grid of 10 slice levels (normal case).
#                IF n_elements(slice_levels) = 3 then a skewed strange histogram is encountered
#                    with mostly shadow pixels, similar to :
#
#                                                  !*
#                                                  ! *
#                                                  !  *
#                                                  !    *
#                                                  !     *
#                                                  !      *        *
#                                                  !         *  *    * * *
#                                                  ----------------------
#                In this case, slice_levels = [ 0.8, 0.6, 0.4 ], is1=1 (slice at 0.6)
#                and n_peak = 1 (no left shoulder of histogram peak) and de-shadowing results
#                may have a poor quality.
#
#        hist_phi_unscl : (normalized) histogram of unscaled shadow function
#
#    Output:
#      list_score_intersection = intarr(n, n_elements(slice_levels))
#                                         contains number of intersection scores per slice level
#                                         only up to n=3 scores are recorded
#      is1        : index where slice_levels[i]  was met
#      n_peak    : index at main peak of histogram
#                     only left part of histogram (index < n_peak) will be searched for
#                     intersection points. If n_peak=1 then the above "strange histogram"
#                     situation is encountered.
#      np_slice : number of slice points (LE 3, only up to 3 will be counted)
#
#
# KEYWORDS:
#    none
#
# COMMON BLOCKS:
#      none
#
# RR, May 2004
#      Jan. 2006 (new case 1 added)


    def _ret():  return (list_score_intersection, is1, n_peak, np_slice)

    list_score_intersection = zeros(3, slice_levels.size)

    li_max = where(ravel(hist_phi_unscl >= 0.95))[0]
    n_peak = li_max[0] + 1
    if (n_peak <= 3):
        # main peak is encountered with the first 3 histogram points
        # "strange" histogram, reset n_peak=1 for further processing
        n_peak = 1
        is1 = 1
        np_slice = 1
        list_score_intersection[is1] = 1
        return _ret()  #  strange histogram


    # make interpolated histogram for level slicing
    n20 = 20
    n_interp = n_peak * n20    # increase interpolated grid by factor n20
    xinterp = arange(n_interp, dtype = float32) / array(n20, copy=0).astype(float32)
    hist_interp = interpol(hist_phi_unscl[0:(n_peak - 1)+1], arange(n_peak), xinterp)
    slice_margin = 0.008

    islice = zeros(slice_levels.size)  # default is zero
    for j in arange(0, (slice_levels.size)):

        hmax1 = slice_levels[j]    # slicing of normalized histogram at level j

        adif = abs(hist_interp - (hmax1)*ones([n_interp]))  #  absolute differences
        # get list of elements within margin
        li_interp_margin = where(ravel(adif <= slice_margin))[0]

        if (len(li_interp_margin) >= 2):
            # remove neighboring elements from list, because 2-3 successive points were met within margin
            li_save = li_interp_margin
            for k in arange(0, (len(li_interp_margin) - 1)):
                if (bitwise_or(li_interp_margin[k] + 1 == li_interp_margin[k + 1], li_save[k] + 1 == li_interp_margin[k + 1])):
                    li_interp_margin[k + 1] = 0
            li_save = 0  # free memory

            # update li_interp_margin
            li_new = where(ravel(li_interp_margin > 0))[0]
            if (len(li_new) > 0):
                li_interp_margin = li_interp_margin[li_new]
            li_new = 0

        np_slice = len(li_new)
        if (np_slice >= 2):
            islice[j] = 1    # at least two slice points are met
            n3 = (minimum(np_slice, 3))
            list_score_intersection[j,0:(n3 - 1)+1] = li_interp_margin[0:(n3 - 1)+1] / n20

            if (n3 == 3):
                break

            # 3 intersection points were found :
            # (a) slicing at high levels (0.50 to 0.80) : high amplitude fluctuations  --> increase binsize
            # (b) slicing at low  levels (0.05 to 0.40) : secondary "shadow" peak was found
            #        --> see sketch below, left part for 3 intersection points
            #      If two points were found the next (lower) slice level will detect
            #      the local "shadow" peak and the slope leading to the main peak
            #      --> see sketch below, right part

        #
        #      3 intersection points                !    2 points
        #                  *                              !                 *
        #                 * *                             !              *  *
        #      *        *    *                            ! xx*xxxxx*x     *
        # __ * *___ *__    *                 or      !  * *    *          *
        #    *    * **         *                         !      ***             *
        # --------------------------              --------------------

    is1 = max(where(ravel(islice == 1))[0])

    return _ret()

#-------------------------------------------------------------------------------------
def threshold_shadow_histo(is1, slice_levels, bins, list_score_intersection, n_peak, hist_phi_unscl, percent_cloud_cover, phi_unscl_min):
# NAME:
#    threshold_shadow_histo  (category: de-shadowing)
#
# PURPOSE:
#        calculate the histogram threshold for the cloud shadow pixels
#
# CALLING SEQUENCE:
# result = threshold_shadow_histo (is1, bins, list_score_intersection, hist_phi_unscl, phi_unscl_min)
#
# PARAMETERS:
#  Input :
#     is1  = index for list_score_intersection(3,is1) for the slice level is1
#     slice_levels          : slice levels of normalized histogram
#     bins = current binsize of hist_phi_unscl
#     list_score_intersection(3,n)  : holds histogram indices of up to 3 intersection points
#     n_peak                  :  histogram index at main peak
#     hist_phi_unscl        :  normalized histogram of fshadow for the non_water pixels
#     phi_unscl_min         :  minimum of fshadow[li_nonwater_rd], possibly after increasing binsize
#     percent_cloud_cover :  float value
#  Output:
#     thr_shad : histogram threshold for core shadow areas
# KEYWORDS:
#    none
#
# RETURN VALUE:
#    histogram threshold for the core shadow areas: thr_shad
#
# COMMON BLOCKS:
#  none
#
# RR, May 2004

# IDL2PY_TBD; IDL2PY_TBD @onerror.inc

    if (is1 > 0):
        # level slicing was successful, at least one point was found.
        # If there is only one slice point then the point is repeated.
        # Determine the histogram "shadow" minimum between the last two slice points
        # Approximation: min is in the center of the last two slice points
        # If there is only one slice point the minimum is taken at this point.
        li = where(ravel(list_score_intersection[is1,:] > 0))[0]  # get number of nonzero elements
        cnt1 = len(list_score_intersection)
        if (cnt1 == 1):
            li1 = array([list_score_intersection[is1,li[0]], list_score_intersection[is1,li[0]]])
        else:
            li1 = array([list_score_intersection[is1,li[cnt1 - 2]], list_score_intersection[is1,li[cnt1 - 1]]])
        li = 0  # free memory
        cnt1 = 2
    else:
        # no intersection points were encountered
        # is1 = -1  --> use slice level 0.10 if intersection point is left of histogram max.
        #                    otherwise take threshold at (phi_unscl_max - phi_unscl_min) /2
        li1 = where(ravel(hist_phi_unscl[0:(n_peak - 1)+1] <= 0.10))[0]
        cnt1 = len(hist_phi_unscl)
        if (cnt1 == 0):
            # define threshold(core_shadow) depending on cloud cover
            # thr_shad = phi_unscl_min + cloud_cover* (phi_unscl_max - phi_unscl_min)
            # the following line is this equation on index basis
            ntt = array(0.01 * percent_cloud_cover * n_peak, copy=0).astype(int32)
            if (ntt < 3):
                ntt = 3
            li1 = array([ntt, ntt - 1])
            cnt1 = 2
        if (cnt1 == 1):
            # take successive two points
            li1 = array([li1[0], li1[0] + 1])
            cnt1 = 2

    i1 = li1[cnt1 - 1]
    i2 = li1[cnt1 - 2]

    # check whether i2 is left of local shadow peak
    # if yes then take intersection right of local shadow peak

    #  local                                          ! local
    #  shadow                                         ! shadow
    #    peak                                          !  peak
    #      .            main                          !  .              main
    #  i3 . i2    i1 peak                          !  . i2    i1    peak
    #    o . o     o  *                              !  . o     o     *
    #    o . o     o  *                              !  . o     o    * *
    #    o . o     o * *                             !  **o     o  *    *
    #    o **o     o*    *                            !     *____o*_      *
    # __o*  *___ *__    *                 or      !      *    *          *
    #    *     * **         *                         !        ***             *
    # --------------------------              --------------------


    fmax_shad = max(hist_phi_unscl[0:(i1)+1])
    li_max_shad = where(ravel(hist_phi_unscl[0:(i1)+1] == fmax_shad))[0]
    imax = li_max_shad[0]
    if (bitwise_and(imax > i2, is1 >= 0)):
        # search for i2 in the region right of imax
        li2 = where(ravel(hist_phi_unscl[imax:(li1[cnt1 - 1])+1] <= slice_levels[is1]))[0]
        i2 = imax + li2[0]
        li2 = 0
    li1 = 0  # free memory


    thr_shad = phi_unscl_min + (i1 + 1) * bins

    return thr_shad


#-------------------------------------------------------------------------------------------
def get_grid_pos(vector, val):

# calculate grid position of "val" in "vector" for later interpolation
# elements of "vector" can be spaced arbitrarily (not equidistant)

    index = where(ravel(vector > val))[0]

    _expr = index[0]
    if _expr == -1:
        g_pos = vector.size - 1.
    elif _expr == 0:
        g_pos = 0.
    else:
        g_pos = (val - vector[index[0] - 1]) / (vector[index[0]] - vector[index[0] - 1]) + index[0] - 1


    return g_pos    # (ranges between 0 and n_elements(vector)- 1)

# ---------------------------------------------------------------
def time_string(time_sec):
# Input : float number of the processing time [sec]
# Output: time as a string such as '1 hours 1 min.' or '10 min. 43 sec'
    m, s = divmod(time_sec, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)


#---------------------------------------------------------------------
def indexvis(visi1, nvisx, vis_reg):

# return the visibility index for given visibility visi1
# input : visi1
#            nvisx ; extended visib set (=55)
#            vis_reg : the nvisx vis. regions
# output: vis.index

    indvis = -1
    if (visi1 > vis_reg[0,1]):
        indvis = 0
    if (visi1 <= vis_reg[nvisx - 2,1]):
        indvis = nvisx - 1
    if (indvis >= 0):
        return indvis

    # now find visi1 in the vis_reg boundaries
    for j in arange(1, nvisx - 1):
        if (bitwise_and(visi1 <= vis_reg[j,0], visi1 >= vis_reg[j,1])):
            indvis = j
            break

    return indvis


#-------------------------------------------------------------------------
def extrapl(y):
# input : vector y[n]
# output: vector y[n+1], where last element is obtained with linear extrapolation
    def _ret():  return y

    n = y.size
    y = array([y, array([y[n - 1] + y[n - 1] - y[n - 2]])])

    return y


#-----------------------------------------------------------------------
def adjacency_weight(nadj_regions, adj_km, pixelsize):

# set zone-dpendent weight factors for adjacency correction

# input : nadj_regions : # adjacency regions (1:5) , intarr(5)
#            adj_km         : adjacency range [km]
#            pixelsize     : pixel size [m]
# output:
#            nadj_pix(*)  = adj. range for each region (pixels), intarr(5)
#            wadj(*)        = weight factor for each region, fltarr(5)

    nadj_pix = 50
    wadj = 1.0

    if (int(adj_km * 1000 / pixelsize + 0.5) > 255):
        nadj_regions = 1    # max adj. range is 255 pixels from center
        return(nadj_pix, wadj)      # i.e. max adj. box is 510 pixels

    _expr = (nadj_regions)
    if _expr == 1:
        adj_rangem = array([1.0]) * adj_km * 1000.0
        wadj = array([1.0])
    elif _expr == 2:
        adj_rangem = array([0.70, 1.0]) * adj_km * 1000.0
        wadj = array([0.5355, 0.4645])        # weight factor
    elif _expr == 3:
        adj_rangem = array([0.65, 0.90, 1.0]) * adj_km * 1000.0
        wadj = array([0.4699, 0.3757, 0.1545])
    elif _expr == 4:
        adj_rangem = array([0.52, 0.75, 0.90, 1.0]) * adj_km * 1000.0
        wadj = array([0.3150, 0.3101, 0.2249, 0.1500])
    elif _expr == 5:
        adj_rangem = array([0.45, 0.65, 0.80, 0.90, 1.0]) * adj_km * 1000.0
        wadj = array([0.2384, 0.2455, 0.2169, 0.1501, 0.1491])
    else:
        pass

    nadj_pix = int(adj_rangem / pixelsize + 0.5)
    return(nadj_pix, wadj)

#-----------------------------------------------------------------------
def load_wv_tables_summer():
# Output:
#    wv_av = 1000
#      uu1 = wv values (cm*1000)
#     suu1 = string array of wv values
#     nuu1 = n_elements(suu1)
#     nuu1, uu1, suu1 will later be updated according to actual number of atmospheric wv files (.atm)

    def _ret():  return (wv_av, uu1, suu1, nuu1, uu1_altit)

    wv_av = 1000  # default wv=1.0 cm, scaled with 1000

    uu1 = array([400, 1002, 2001, 2904, 4000, 5007])    #  scaled wv=0.4, 1.0, 2.0, 2.9, 4.0, 5.0 cm
    suu1 = array(['04', '10', '20', '29', '40', '50'])    #  (refers to sea level)
    nuu1 = uu1.size

    # 41 altitudes for wv interpolation: 0(0.1)4.0 km
    # These indices (0:40) correspond to 100 m DEM elevation classes
    # Currently the range 0-3.5 km is used in ATCOR.

    uu1_altit = zeros([41, nuu1], float32)

    # MOD43 values:
    # wv=400, 1000, 2000, 2900 is generated with model=2 and H2OSTR=0.137, 0.343, 0.685, 0.994
    # wv=4000 is calculated with model=2 and H2OSTR=1.36878
    # wv=5000 is calculated with model=2 and H2OSTR=1.71360

    # scaled wv column 400, all 46 altitude entries
    uu1_altit[0:10,0] = array([400, 381, 363, 346, 329, 313, 298, 283, 269, 256])
    uu1_altit[10:20,0] = array([243, 230, 219, 207, 196, 186, 176, 166, 157, 149])
    uu1_altit[20:30,0] = array([141, 133, 125, 118, 111, 105, 99, 94, 89, 84])
    uu1_altit[30:40,0] = array([79, 74, 70, 66, 63, 59, 44, 33, 25, 19])
    uu1_altit[40,0] = array([17])

    # scaled wv column 1000, all 46 altitude entries
    uu1_altit[0:10,1] = array([1002, 955, 909, 866, 824, 785, 746, 710, 675, 641])
    uu1_altit[10:20,1] = array([609, 578, 548, 519, 492, 466, 441, 417, 395, 373])
    uu1_altit[20:30,1] = array([353, 333, 314, 297, 280, 264, 249, 236, 222, 210])
    uu1_altit[30:40,1] = array([198, 187, 177, 167, 158, 149, 111, 83, 63, 48])
    uu1_altit[40,1] = array([43])

    # scaled wv column 2000, all 46 altitude entries
    uu1_altit[0:10,2] = array([2001, 1907, 1817, 1730, 1647, 1567, 1491, 1418, 1348, 1281])
    uu1_altit[10:20,2] = array([1216, 1154, 1095, 1038, 983, 930, 881, 834, 789, 746])
    uu1_altit[20:30,2] = array([705, 666, 628, 593, 559, 528, 499, 471, 445, 420])
    uu1_altit[30:40,2] = array([396, 374, 353, 334, 315, 298, 223, 167, 126, 96])
    uu1_altit[40,2] = array([87])

    # scaled wv column 2900, all 46 altitude entries
    uu1_altit[0:10,3] = array([2904, 2768, 2637, 2511, 2390, 2274, 2164, 2057, 1956, 1858])
    uu1_altit[10:20,3] = array([1765, 1675, 1589, 1506, 1426, 1350, 1278, 1210, 1145, 1083])
    uu1_altit[20:30,3] = array([1023, 966, 912, 861, 812, 767, 724, 683, 645, 609])
    uu1_altit[30:40,3] = array([575, 543, 513, 484, 457, 432, 323, 242, 183, 140])
    uu1_altit[40,3] = array([125])

    # scaled wv column 4000, all 46 altitude entries
    uu1_altit[0:10,4] = array([4000, 3811, 3631, 3458, 3291, 3132, 2980, 2833, 2693, 2559])
    uu1_altit[10:20,4] = array([2431, 2307, 2189, 2074, 1964, 1859, 1761, 1666, 1577, 1491])
    uu1_altit[20:30,4] = array([1409, 1331, 1256, 1185, 1118, 1056, 997, 941, 889, 839])
    uu1_altit[30:40,4] = array([792, 748, 706, 667, 630, 595, 445, 334, 253, 193])
    uu1_altit[40,4] = array([152])

    # scaled wv column 5000, all 46 altitude entries
    uu1_altit[0:10,5] = array([5007, 4772, 4546, 4329, 4121, 3921, 3730, 3547, 3372, 3204])
    uu1_altit[10:20,5] = array([3043, 2889, 2740, 2597, 2459, 2328, 2204, 2086, 1974, 1867])
    uu1_altit[20:30,5] = array([1765, 1666, 1573, 1484, 1400, 1322, 1248, 1179, 1113, 1051])
    uu1_altit[30:40,5] = array([992, 937, 884, 835, 789, 745, 558, 418, 317, 242])
    uu1_altit[40,5] = array([172])

    return _ret()

#-----------------------------------------------------------------------------------
def load_wv_tables_winter():
#
# Output:
#    wv_av =  400
#      uu1 = wv values (cm*1000)
#     suu1 = string array of wv values
#     nuu1 = n_elements(suu1)
#     nuu1, uu1, suu1 will later be updated according to actual number of atmospheric wv files (.atm or .bp7)

    def _ret():  return (wv_av, uu1, suu1, nuu1, uu1_altit)

    wv_av = 400  # default wv=0.4 cm, scaled with 1000

    uu1 = array([200, 400, 800, 1100])    #  scaled wv=0.2, 0.4, 0.8, 1.1 cm
    suu1 = array(['02', '04', '08', '11'])    #  (refers to sea level)
    nuu1 = uu1.size

    # 36 altitudes for wv interpolation: 0(0.1)3.5 km
    # These indices (0:35) correspond to 100 m DEM elevation classes

    uu1_altit = zeros([36, nuu1], float32)

    # MOD43 values:
    # wv= 200, 400, 800, 1100 is generated with model=3 and H2OSTR=0.23482, 0.46964,  0.93928, 1.29155
    # all are calulated with the same model to have the same trace gas columns for all wv levels

    # scaled wv column  200, all 36 altitude entries
    uu1_altit[0:10,0] = array([199, 191, 184, 176, 169, 162, 155, 148, 142, 136])
    uu1_altit[10:20,0] = array([130, 124, 119, 113, 108, 103, 98, 93, 89, 84])
    uu1_altit[20:30,0] = array([80, 76, 72, 68, 65, 61, 58, 54, 51, 48])
    uu1_altit[30:36,0] = array([45, 43, 40, 37, 35, 33])

    # scaled wv column  400, all 36 altitude entries
    uu1_altit[0:10,1] = array([399, 383, 368, 353, 338, 324, 311, 297, 285, 272])
    uu1_altit[10:20,1] = array([260, 249, 238, 227, 217, 207, 197, 187, 178, 169])
    uu1_altit[20:30,1] = array([160, 152, 144, 137, 130, 122, 116, 109, 103, 97])
    uu1_altit[30:36,1] = array([91, 86, 80, 75, 71, 66])

    # scaled wv column  800, all 36 altitude entries
    uu1_altit[0:10,2] = array([799, 767, 736, 706, 677, 649, 622, 595, 570, 545])
    uu1_altit[10:20,2] = array([521, 498, 476, 455, 434, 414, 394, 375, 356, 339])
    uu1_altit[20:30,2] = array([321, 305, 289, 274, 260, 245, 232, 218, 206, 194])
    uu1_altit[30:36,2] = array([182, 172, 161, 151, 142, 133])

    # scaled wv column 1100, all 36 altitude entries
    uu1_altit[0:10,3] = array([1100, 1055, 1012, 971, 931, 892, 855, 819, 784, 750])
    uu1_altit[10:20,3] = array([717, 685, 655, 625, 597, 569, 542, 515, 490, 466])
    uu1_altit[20:30,3] = array([442, 420, 398, 377, 357, 337, 319, 301, 283, 267])
    uu1_altit[30:36,3] = array([251, 236, 222, 208, 196, 184])

    return _ret()

#------------------------------------------------------------------------------------------
def read_wv_trans945_1375(gamma, solze, h1_cirrus):

# Purpose: calculate the two-way (sun-cirrus-sensor) wv transmittance
#             for the 945 nm band of Sentinel2
#
# Input:
#             gamma = transmittance of the 1.38 micron band  (sun-cirrus-sensor)
#             solze = solar zenith angle (deg)
# Output:  trans945  = two-way (sun-cirrus-sensor) wv transmittance 945 nm
#             h1_cirrus = approximate height of cirrus (km)
#                             based on the US standard temperature/humidity profile

#  trans for 1375 nm and R=t(945)/t1375  for US standard 1976 atmosphere
#  R is independent of selected atmosphere (for h> 6 km)

# updated for final S2 filter functions (Sept. 2011)

    arr = array([0.28593, 3.00129, 0.28436, 3.01551, 0.27951, 3.06023, 0.27091, 3.14325, \
                             0.25765, 3.28119, 0.23796, 3.51201, 0.20821, 3.93628, 0.15891, 4.95467, \
                             0.44296, 2.05593, 0.44135, 2.06239, 0.43633, 2.08275, 0.42738, 2.12018, \
                             0.41336, 2.18185, 0.39207, 2.28328, 0.35866, 2.46481, 0.29931, 2.87851, \
                             0.61264, 1.55084, 0.61129, 1.55381, 0.60704, 1.56317, 0.59942, 1.58030, \
                             0.58731, 1.60836, 0.56860, 1.65383, 0.53835, 1.73362, 0.48151, 1.90849, \
                             0.75786, 1.28681, 0.75690, 1.28824, 0.75389, 1.29269, 0.74846, 1.30085, \
                             0.73982, 1.31400, 0.72627, 1.33520, 0.70395, 1.37168, 0.66037, 1.44940, \
                             0.84854, 1.16355, 0.84790, 1.16433, 0.84587, 1.16678, 0.84221, 1.17125, \
                             0.83635, 1.17843, 0.82717, 1.18983, 0.81193, 1.20921, 0.78172, 1.24949, \
                             0.93536, 1.06474, 0.93505, 1.06505, 0.93410, 1.06603, 0.93240, 1.06778, \
                             0.92966, 1.07060, 0.92537, 1.07502, 0.91829, 1.08231, 0.90424, 1.09694, \
                             0.96321, 1.03605, 0.96304, 1.03622, 0.96248, 1.03675, 0.96148, 1.03773, \
                             0.95990, 1.03927, 0.95743, 1.04167, 0.95337, 1.04562, 0.94541, 1.05326])

    arr = reshape(arr, (7,8,2)) # two columns: t(1375), R=t(945)/t(1375)
    hh = array([6., 7., 8., 9., 10., 12., 14.])
    nh = hh.size

    arr1375 = (arr[:,:,0])
    ratio = arr[:,:,1]
    x = arange(nh, dtype = float32)
    y = ones(nh, dtype = float32) * 0.1 * solze
    xy = array([x,y])

    tr1375 = ravel(map_coordinates(arr1375, xy, order=1, mode='nearest'))
    ratio1 = ravel(map_coordinates(ratio, xy, order=1, mode='nearest'))

    h1_cirrus = interpol(hh, tr1375, gamma)
    r1 = interpol(ratio1, hh, h1_cirrus)

    tr945 = (minimum(r1 * gamma, 1.0))
    return tr945, h1_cirrus

#-----------------------------------------------------------------------
def image_cells(nrow, ncol, cell_length_x):

# Input:
#     nrow                = rows of image     (y direction)
#     ncol                = columns of image (x direction)
#     cell_length_x = cell length in x direction (pixels)
#     Comment:     cell_length_y = 2*cell_length_x  as Delta(SZA) much smaller than Delta(VZA)
#                     i.e. the cell length in y is twice the size of the length in x direction
#                     Advantage: increase in speed and more DDV reference pixels can be found per cell
#
# Output:
#     nx_cell =  number of cells in x direction
#     ny_cell =     "          "    y    "
#     xcell    = lonarr(2,nx_cell) left/right pixel coordinates of all x cells
#     ycell    = lonarr(2,ny_cell) left/right pixel coordinates of all y cells

    cell_length_y = 2 * cell_length_x

    nx_cell = int(maximum((ncol / cell_length_x)+0.5, 1))
    ny_cell = int(maximum((nrow / (cell_length_y))+0.5, 1))
    nx = 0
    ny = 0

    if (cell_length_x < ncol):
        nx = ncol % cell_length_x
    if (cell_length_y < nrow):
        ny = nrow % cell_length_y

    if (nx > 0):
        nx_cell = nx_cell + 1
    if (ny > 0):
        ny_cell = ny_cell + 1

    xcell = zeros((nx_cell,2), long)  # stores left/right pixel indices for each xcell (starting with 0)
    ycell = zeros((ny_cell,2), long)  #      "                  "          "                ycell

    # calculate x cell pixel coordinates
    p1 = -cell_length_x
    p2 = -1

    for i in arange(0, (nx_cell)):
        p1 = p1 + cell_length_x
        p2 = minimum((p2 + cell_length_x), (ncol - 1))
        xcell[i,0:2] = [p1, p2]

    # calculate y cell pixel coordinates
    p1 = -cell_length_y
    p2 = -1

    for i in arange(0, (ny_cell)):
        p1 = p1 + cell_length_y
        p2 = minimum((p2 + cell_length_y), (nrow - 1))
        ycell[i,0:2] = [p1, p2]

    return ny_cell, nx_cell, ycell, xcell


#-----------------------------------------------------------------------
def column_subset(vza_arr, ncol, itilt):

# Purpose:
# Check is performed for the close-to-nadir part of the scene with VZA < 10 deg
# In case of tilt the 10 degrees closest to nadir are selected.
# ncol is simply scaled assuming pixel 1 corresponds to VZA[0,0] = UL
#                                 and pixel ncol corresponds to VZA[1,0] = UR
# This is usually also sufficient to guarantee a small SZA difference.
# return the corresponding column subset nc1 to nc2 (starting from 0)
#
# Input; vza_arr = fltarr(2,2) : VZA for UL, UR and LL, LR
#          ncol     = number of image columns
#          itilt    = 0 no tilt, =1 tilt capability
# Output:  nc1, nc2

# calculate column pixels for close-to-nadir subset

    delta_vza = vza_arr[0,0] - vza_arr[0,1]                         # vza_arr( UL - UR)

    if (itilt == 0):
        # for S2 with VZA < 10.5 degrees the close-to-nadir region includes the full swath
        # (no tilt)
        nc1 = 0
        nc2 = ncol - 1
    else:
        # tilt capability
        fact = (minimum(10.0 / abs(delta_vza), 1.0))
        if (vza_arr[0,0] > vza_arr[0,1]):
            nc1 = array((1.0 - fact) * ncol, copy=0).astype(int32)  ;  nc2 = ncol - 1
        else:
            nc1 = 1  ;  nc2 = (minimum(array(ncol * fact, copy=0).astype(int32), (ncol - 1)))

    return array([nc1, nc2])


#-----------------------------------------------------------------------
def set_nadir_geometry(sza_arr, saa_arr, vza_arr, vaa_arr, itilt):
# Set nadir geometry or in case of tilt the nearest possible nadir region
# Input:
#     sza_arr  = fltarr(2,2) solar zenith angles (UL, UR, LL, LR)
#     saa_arr  =         "              azimuth  "
#     vza_arr  =         "        view zenith    "
#     vaa_arr  =         "        view azimuth  "
#     itilt     =0 no tilt, =1=tilt capability
# Output:
#      4-element vector [vza, vaa, sza, saa] with (close) nadir geometry

    if (itilt == 0):
        vza = 0.0
        vaa = 0.0  # irrelevant for vza=0
        saa = 0.0  # irrelevant for vza=0

        if (ravel(abs(vza_arr[0,0] - vza_arr[0,1]))[0] < 10.0):
            sza = float32(mean(sza_arr))
        else:
            sza = 0.5 * (bilinear_interpolation(sza_arr, 1.0 - vza_arr[0,1] / vza_arr[0,0], 0) + bilinear_interpolation(sza_arr, 1.0, 1.0 - vza_arr[1,1] / vza_arr[1,0]))
    else:
        delt = abs(vza_arr[0,0] - vza_arr[0,1])
        fact = (minimum(10.0 / delt, 1.0))                    # 10 deg closest to nadir
        if (vza_arr[0,0] < vza_arr[0,1]):
            ind = 0
        else:
            ind = 1
        vza = vza_arr[0,ind] + fact * delt
        vaa = vaa_arr[0,ind] + fact * abs(vaa_arr[0,0] - vaa_arr[0,1]) + fact * abs(vaa_arr[1,ind] - vaa_arr[0,ind])
        sza = sza_arr[0,ind] + fact * abs(sza_arr[0,0] - sza_arr[0,1]) + fact * abs(sza_arr[1,ind] - sza_arr[0,ind])
        saa = saa_arr[0,ind] + fact * abs(saa_arr[0,0] - saa_arr[0,1]) + fact * abs(saa_arr[1,ind] - saa_arr[0,ind])

    return array([vza, vaa, sza, saa],dtype=float32)

