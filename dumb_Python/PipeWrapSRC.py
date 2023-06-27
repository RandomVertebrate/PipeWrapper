from numpy import *
from numpy.linalg import *
import matplotlib.pyplot as plt
import pandas as pd

# Plotting helper functions for equal-scale axes in 3D
# Functions from @Mateen Ulhaq and @karlo
def set_axes_equal(ax: plt.Axes):

    limits = array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])
    origin = mean(limits, axis=1)
    radius = 0.5 * max(abs(limits[:, 1] - limits[:, 0]))
    _set_axes_radius(ax, origin, radius)

def _set_axes_radius(ax, origin, radius):
    x, y, z = origin
    ax.set_xlim3d([x - radius, x + radius])
    ax.set_ylim3d([y - radius, y + radius])
    ax.set_zlim3d([z - radius, z + radius])
# End plotting helper functions

def Triangle3dTo2d(A3d, B3d, C3d, A2d, B2d):
    BaseVec3d = B3d - A3d
    BaseDir3d = BaseVec3d/norm(BaseVec3d)
    BasePt3d = A3d + BaseDir3d*(dot((C3d - A3d), BaseDir3d))
    height = norm(C3d - BasePt3d)
    footdist = norm(BasePt3d - A3d)
    
    A2d = append(A2d, 0)
    B2d = append(B2d, 0)
    
    BaseVec2d = B2d - A2d
    BaseDir2d = BaseVec2d/norm(BaseVec2d)
    
    HeightDir2d = cross(BaseDir2d, array([0, 0, 1]))
    
    C2d = A2d + footdist*BaseDir2d + height*HeightDir2d
    
    return(C2d[0:2])

def axang2rotm(axang):
    axis = array(axang[0:3:1])
    axis = axis/norm(axis)
    angle = axang[3]

    S = array([[0, -axis[2], axis[1]],[axis[2], 0, -axis[0]],[-axis[1], axis[0], 0]])
    
    rotm = cos(angle)*eye(3) + (1-cos(angle))*matmul(transpose([axis]), [axis]) + sin(angle)*S

    return rotm

def GenPipeline(cornersx, cornersy, cornersz, bendradii, bendpoints, showplot = False):

    numbends = len(cornersx) - 2

    if bendpoints < 2:
        raise Exception("Bendpoints cannot be less than 2")
    
    if not len(cornersx) == len(cornersy) == len(cornersz):
        raise Exception("Inconsistent number of input coordinates")
    
    if len(bendradii) != numbends:
        raise Exception("Incorrect number of bend radii specified")
    
    corners = array([cornersx, cornersy, cornersz])
    
    pipeline = zeros((3, numbends*bendpoints + 2))
    pipeline[:,0] = [cornersx[0], cornersy[0], cornersz[0]]
    
    for i in range(numbends):
    
        oldDirection = (corners[:,i+1] - corners[:,i])/norm(corners[:,i+1] - corners[:,i])
    
        newDirection = (corners[:, i+2] - corners[:, i+1])/norm(corners[:, i+2] - corners[:, i+1])
    
        cornerangle = arccos(-(dot(newDirection, oldDirection)))
    
        bendAngle = pi - cornerangle
    
        AngInc = bendAngle/(bendpoints - 1)
    
        bendAxis = cross(oldDirection, newDirection)
        if norm(bendAxis) == 0:
            raise Exception("Bend angle cannot be zero")
        bendAxis = bendAxis/norm(bendAxis)
    
        startradvec = cross(oldDirection, bendAxis)*bendradii[i]
    
        bendCenter = corners[:, i+1] - oldDirection*bendradii[i]/tan(cornerangle/2) - startradvec
        
        for j in range(bendpoints):
        
            pipeline[:, i*bendpoints+j+1] = bendCenter + matmul(axang2rotm([bendAxis[0], bendAxis[1], bendAxis[2], AngInc*j]), transpose([startradvec]))[:,0]

    pipeline[:, numbends*bendpoints+1] = [cornersx[numbends+1], cornersy[numbends+1], cornersz[numbends+1]]

    if showplot:
        fig = plt.figure()
        ax = fig.add_subplot(projection = '3d')
    
        ax.plot(pipeline[0, :], pipeline[1, :], pipeline[2, :])

        plt.title("Pipe Centerline")
    
        plt.show(block = False)

    return pipeline

def WrapPipe(Profile, PipeRadius, TurnsPerMeter, Overlap, Resolution, PlotAngle, RefVector, OutputFile = None, fig = None):

    linein = Profile

    maxpts = len(transpose(linein))

    Spacing = 1/Resolution
    line = transpose([linein[:, 0]])

    partition = 0

    while norm(linein[:, -1] - line[:, -1]) > Spacing:
        for j in range(partition, maxpts):
            potentialIncrement = linein[:, j] - line[:, -1]
            if norm(potentialIncrement) > Spacing:
                incrementDirection = potentialIncrement/norm(potentialIncrement)
                actualIncrement = incrementDirection*Spacing
                line = append(line, transpose([line[:, -1] + actualIncrement]), axis = 1)
                partition = j
                break

    angrate = TurnsPerMeter*2*pi*Spacing
    
    ijump = int(round(2*pi/angrate))
    
    npts = len(transpose(line))
    
    if angrate*npts<2*pi:
        raise Exception("Wrap rate too low")
    
    directions = diff(line, axis = 1)
    directions = append(directions, transpose([directions[:, -1]]), axis = 1)
    for i in range(npts):
        directions[:, i] = directions[:, i]/norm(directions[:, i])
    
    nonParallel = zeros((3, npts))
    
    nonParallel[:, 0] = RefVector

    for i in range(1, npts-1):
        rotax = cross(directions[:, i-1], directions[:, i])
        if norm(rotax) == 0:
            nonParallel[:, i] = nonParallel[:, i-1]
        else:
            rotax = rotax/norm(rotax)
            rotang = real(arccos(dot(directions[:, i-1], directions[:, i])+0j))
            nonParallel[:, i] = matmul(axang2rotm(append(rotax, rotang)), transpose([nonParallel[:, i-1]]))[:, 0]
    
    nonParallel[:, -1] = nonParallel[:, -2]
    
    offsets = zeros((3, npts))
    
    for i in range(npts):
        refvec = cross(directions[:, i], nonParallel[:, i])
        rotang = i*angrate
        rotax = directions[:, i]
        offsets[:, i] = matmul(axang2rotm(append(rotax, rotang)), transpose([refvec]))[:, 0]
        offsets[:, i] = PipeRadius*offsets[:, i]/norm(offsets[:, i])
    
    seam = line + offsets

    if fig == None:
        fig = plt.figure()
    
    ax1 = fig.add_subplot(2, 2, 3, projection = "3d")
    ax1.clear()
    ax1.plot(line[0, :], line[1, :], line[2, :], 'g', linewidth = 2)
    ax1.plot(seam[0, :], seam[1, :], seam[2, :], 'r', linewidth = 0.7)
    
    for i in range(npts-ijump):
        ax1.plot([seam[0, i], seam[0, i+ijump]], [seam[1, i], seam[1, i+ijump]], [seam[2, i], seam[2, i+ijump]], 'k', linewidth = 0.2)
    
    plt.title("Wrapped Pipe")
    
    ax1.set_box_aspect([1, 1, 1])
    set_axes_equal(ax1)
    plt.show(block = False)

    line1_3d = seam[:, 0:npts-ijump]
    line2_3d = seam[:, ijump:npts]
    
    line1_2d = zeros((2, npts-ijump))
    line2_2d = zeros((2, npts-ijump))

    line1_2d[:, 0] = [0, 0];

    width0 = norm(line1_3d[:, 0] - line2_3d[:, 0])

    line2_2d[:, 0] = [width0*cos(PlotAngle), width0*sin(PlotAngle)]
    
    for i in range(npts-ijump-1):
        if i%2 == 0:
            A3d = line1_3d[:, i]
            B3d = line2_3d[:, i]
            C3d = line1_3d[:, i+1]
            A2d = line1_2d[:, i]
            B2d = line2_2d[:, i]
            C2d = Triangle3dTo2d(A3d, B3d, C3d, A2d, B2d)
            line1_2d[:, i+1] = C2d
            A3d = line1_3d[:, i+1]
            C3d = line2_3d[:, i+1]
            A2d = line1_2d[:, i+1]
            C2d = Triangle3dTo2d(A3d, B3d, C3d, A2d, B2d)
            line2_2d[:, i+1] = C2d
        else:
            A3d = line1_3d[:, i]
            B3d = line2_3d[:, i]
            C3d = line2_3d[:, i+1]
            A2d = line1_2d[:, i]
            B2d = line2_2d[:, i]
            C2d = Triangle3dTo2d(A3d, B3d, C3d, A2d, B2d)
            line2_2d[:, i+1] = C2d
            B3d = line2_3d[:, i+1]
            C3d = line1_3d[:, i+1]
            B2d = line2_2d[:, i+1]
            C2d = Triangle3dTo2d(A3d, B3d, C3d, A2d, B2d)
            line1_2d[:, i+1] = C2d
    
    dottedline = zeros((2, npts-ijump))
    
    for i in range(npts-ijump):
        dottedline[:, i] = line2_2d[:, i] + Overlap*(line2_2d[:, i] - line1_2d[:, i])/norm(line2_2d[:, i] - line1_2d[:, i])
    
    ax2 = fig.add_subplot(2, 2, 4)
    ax2.clear()
    ax2.plot(line1_2d[0, :], line1_2d[1,:], 'k', linewidth = 0.5)
    ax2.plot(line2_2d[0, :], line2_2d[1,:], 'k', linewidth = 0.5)
    ax2.plot(dottedline[0, :], dottedline[1, :], 'k--', linewidth = 0.5)
    
    numturns = int(ceil(npts/ijump))

    for i in range(numturns-1):
        idx = i*ijump
        pt1 = line1_2d[:, idx]
        pt2 = line2_2d[:, idx]
        ax2.plot([pt1[0], pt2[0]], [pt1[1], pt2[1]], 'k', linewidth = 0.5)
    
    pt1 = line1_2d[:, -1]
    pt2 = line2_2d[:, -1]
    ax2.plot([pt1[0], pt2[0]], [pt1[1], pt2[1]], 'k', linewidth = 0.5)
    
    ax2.set_aspect("equal")
    plt.title("Wrapper")
    plt.show(block = False)

    if not OutputFile == None:
        df = pd.DataFrame(transpose([line1_2d[0, :], line1_2d[1, :], line2_2d[0, :], line2_2d[1, :], dottedline[0, :], dottedline[1, :]]))
        df.to_excel(OutputFile, header = ["line 1 x", "line 1 y", "line 2 x", "line 2 y", "dotted line x", "dotted line y"], index = False)
        
    return
