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

# Triangle3dTo2d() Gives the third point in a 2D drawing of
# a triangle ABC defined in 3D space.
# (Transforms, rather than projects, the triangle).
# 3D points are vectors [x, y, z] and 2D points are [x, y].
def Triangle3dTo2d(A3d, B3d, C3d, A2d, B2d):
    BaseVec3d = B3d - A3d                                       # 3D vector AB
    BaseDir3d = BaseVec3d/norm(BaseVec3d)                       # Direction of 3D AB
    BasePt3d = A3d + BaseDir3d*(dot((C3d - A3d), BaseDir3d))    # Foot of perpendicular line dropped from C onto AB
    height = norm(C3d - BasePt3d)                               # Height of perpendicular line
    footdist = norm(BasePt3d - A3d)                             # Distance of foot from A
    # Defining the 2D points in 3D space
    A2d = append(A2d, 0)
    B2d = append(B2d, 0)
    
    BaseVec2d = B2d - A2d                                       # 2D vector AB
    BaseDir2d = BaseVec2d/norm(BaseVec2d)                       # Direction of 2D AB
    
    HeightDir2d = cross(BaseDir2d, array([0, 0, 1]))            # Vector perpendicular to AB
    
    C2d = A2d + footdist*BaseDir2d + height*HeightDir2d         # Finding 2D C
    
    return(C2d[0:2])                                            # Removing the (zero) z component of C

# Converts axis-angle [x, y, z, theta] to a rotation matrix
def axang2rotm(axang):
    axis = array(axang[0:3:1])                                  # The axis
    axis = axis/norm(axis)                                      # Normalizing for good measure
    
    angle = axang[3]                                            # The angle

    # Skew symmetric matrix corresponding to cross product with axis
    S = array([[0, -axis[2], axis[1]],[axis[2], 0, -axis[0]],[-axis[1], axis[0], 0]])
    # Rotation matrix formula
    rotm = cos(angle)*eye(3) + (1-cos(angle))*matmul(transpose([axis]), [axis]) + sin(angle)*S

    return rotm


# Given a sequence of points in 3D space, GenPipeline() connects them,
# rounds the corners, and returns the resulting curve (discretized, so a
# set of points)
# cornersx, cornersy, cornersz contain the coordinates of the input points
# bendradii contains the radii of the bends
# bendpoints is the resolution (number of subdivisions) of each bend
# showplot specifies whether to 3D-plot the generated pipeline (in a new window)
def GenPipeline(cornersx, cornersy, cornersz, bendradii, bendpoints, showplot = False):

    # Number of bends is 2 less than number of points (all points are bends except first and last)
    numbends = len(cornersx) - 2

    if bendpoints < 2:
        raise Exception("Bendpoints cannot be less than 2")
    
    if not len(cornersx) == len(cornersy) == len(cornersz):
        raise Exception("Inconsistent number of input coordinates")
    
    if len(bendradii) != numbends:
        raise Exception("Incorrect number of bend radii specified")

    corners = array([cornersx, cornersy, cornersz])     # Each column is an input point

    # pipeline[][] will store the output. Each column will be a point.
    pipeline = zeros((3, numbends*bendpoints + 2))
    # The first point will be unchanged.
    pipeline[:,0] = [cornersx[0], cornersy[0], cornersz[0]]
    
    for i in range(numbends):
        # Unit vector parallel to pipe before bend
        oldDirection = (corners[:,i+1] - corners[:,i])/norm(corners[:,i+1] - corners[:,i])
        # Unit vector parallel to pipe after bend
        newDirection = (corners[:, i+2] - corners[:, i+1])/norm(corners[:, i+2] - corners[:, i+1])
        # Corner angle (if there was no bend) is the inverse cosine of the dot
        # product of the unit vectors
        cornerangle = arccos(-(dot(newDirection, oldDirection)))
        # The angle through which an arc will be drawn. Equal to pi-cornerangle
        # for geometry reasons
        bendAngle = pi - cornerangle
        # The "arc angle" between each pair of subsequent points in the bend.
        AngInc = bendAngle/(bendpoints - 1)
        # The bend axis is perpendicular to the pipe before and after the bend.
        bendAxis = cross(oldDirection, newDirection)
        if norm(bendAxis) == 0:
            raise Exception("Bend angle cannot be zero")
        bendAxis = bendAxis/norm(bendAxis)
        # Vector pointing from center of the arc to the start of the bend (end
        # of the straight section of pipe)
        startradvec = cross(oldDirection, bendAxis)*bendradii[i]
        # This gives the position vector of the center of the arc.
        bendCenter = corners[:, i+1] - oldDirection*bendradii[i]/tan(cornerangle/2) - startradvec
        # calculating points along bend and storing to output array
        for j in range(bendpoints):
            # The jth point along the bend is obtained by rotating startradvec
            # through j*AngInc about bendAxis, then adding this rotated vector
            # to bendCenter
            pipeline[:, i*bendpoints+j+1] = bendCenter + matmul(axang2rotm([bendAxis[0], bendAxis[1], bendAxis[2], AngInc*j]), transpose([startradvec]))[:,0]

    # The last point is unchanged.
    pipeline[:, numbends*bendpoints+1] = [cornersx[numbends+1], cornersy[numbends+1], cornersz[numbends+1]]

    if showplot:
        # PLOT PIPE PROFILE
        fig = plt.figure()
        ax = fig.add_subplot(projection = '3d')
    
        ax.plot(pipeline[0, :], pipeline[1, :], pipeline[2, :])

        plt.title("Pipe Centerline")
    
        plt.show(block = False)

    return pipeline


# Function WrapPipe()
# Profile is a matrix in which each column is [x, y, z] coordinates of
# a point along the centerline of the pipe
# PipeRadius is the radius of the pipe in meters.
# TurnsPerMeter is the number of turns of wrapping per meter of pipe.
# Overlap is the overlap between successive wraps/turns. A dashed line is
# added to the output to mark it.
# Resolution is calculation spatial resolution in points per unit length
# Convergence of the output should be checked by considering progressively
# larger Resolution values.
# PlotAngle rotates the output.
# RefVector is used as a reference for the helix angle and effectively
# decides where the wrapping starts. It must be non-parallel to the pipe at
# the pipe's starting point.
def WrapPipe(Profile, PipeRadius, TurnsPerMeter, Overlap, Resolution, PlotAngle, RefVector, OutputFile = None, ax2d = None, ax3d = None):

    RefVector = RefVector/norm(RefVector)     # Normalizing reference vector. Not strictly necessary but doesn't hurt.

    linein = Profile                          # Input data. Each column is a point [x, y, z]
    maxpts = len(transpose(linein))           # Number of points read

    # GENERATING EVENLY SPACED INTERPOLATION OF PIPE CENTERLINE
    # We need the centerline to consist of (physically) evenly spaced points so
    # that we can use an array index to determine how far along the line we are.
    Spacing = 1/Resolution                    # The physical spacing corresponding to the chosen resolution

    if Spacing*TurnsPerMeter > 0.5:
        # Resolution too low: Nyquist, aliasing, something.
        raise Exception("Resolution too low")

    line = transpose([linein[:, 0]])          # Generated interpolation will be stored here. The starting point is unchanged.
    # The variable partition keeps track of how much of the input data has been processed.
    # (All up to the (partition-1)'th point have been processed and
    # the partition'th point is currently being processed)
    partition = 0
    # This loop goes through the input data, generates evenly spaced points, and adds them to the output data.
    while norm(linein[:, -1] - line[:, -1]) > Spacing:           # Keep going until the last generated point is within Spacing distance of the last input point
        for j in range(partition, maxpts):                       # Iterate through the unprocessed input data
            # If this point were added as-is to the generated data, that would
            # correspond to a length increment of
            potentialIncrement = linein[:, j] - line[:, -1]
            # Wait until the potential increment is longer than Spacing
            if norm(potentialIncrement) > Spacing:
                # Now add an increment in this direction but of magnitude Spacing
                incrementDirection = potentialIncrement/norm(potentialIncrement)
                actualIncrement = incrementDirection*Spacing
                # New point added to generated data
                line = append(line, transpose([line[:, -1] + actualIncrement]), axis = 1)
                # Update partition. Note: partition'th point will be checked
                # again now, just not points to its left. This means that
                # further increments towards the partition'th point are
                # possible if necessary.
                partition = j
                break

    angrate = TurnsPerMeter*2*pi*Spacing           # Converting wrap rate to radians per point
    
    ijump = int(round(2*pi/angrate))               # Index difference between subsequent wraps
    
    npts = len(transpose(line))                    # Number of points in the generated evenly spaced pipe centerline
    
    if angrate*npts<2*pi:
        raise Exception("Wrap rate too low")

    # Array of vectors along pipe centerline    
    directions = diff(line, axis = 1)
    # Duplicating the last vector so that directions has length npts
    directions = append(directions, transpose([directions[:, -1]]), axis = 1)
    # Normalizing directions
    for i in range(npts):
        directions[:, i] = directions[:, i]/norm(directions[:, i])

    # GENERATING NON-PARALLEL VECTORS THAT FOLLOW THE CURVES OF THE PIPE
    # The generation of a helix that follow the pipe will be done by taking
    # vectors perpendicular to the pipe centerline and rotating them about it
    # by appropriate angles.
    # To find these perpendicular vectors, the pipe directions have to be
    # crossed with something.
    # The following loop finds those somethings.
    # Each new non-parallel vector is found by rotating the previous one
    # through whatever angle the pipe direction changes by.
    nonParallel = zeros((3, npts))               # Initializing    
    nonParallel[:, 0] = RefVector                # The first non-parallel vector is RefVector

    for i in range(1, npts-1):
        # The axis of rotation is normal to the plane of the pipeline, i.e.,
        # perpendicular to two consecutive pipe directions.
        rotax = cross(directions[:, i-1], directions[:, i])
        # If the cross product is zero, there is no change in direction and the
        # non-parallel vector need not be changed.
        if norm(rotax) == 0:
            nonParallel[:, i] = nonParallel[:, i-1]
        else:
            rotax = rotax/norm(rotax)          # Normalizing, just in case
            rotang = real(arccos(dot(directions[:, i-1], directions[:, i])+0j))
            # Rotating
            nonParallel[:, i] = matmul(axang2rotm(append(rotax, rotang)), transpose([nonParallel[:, i-1]]))[:, 0]
    # Duplicating the last non-parallel vector for array length matching
    nonParallel[:, -1] = nonParallel[:, -2]

    # GENERATING A 'HELIX' THAT FOLLOWS THE PIPE CENTERLINE
    offsets = zeros((3, npts))           # This will store offset vectors between the pipe centerline and a helix that follows it
    # Iterating through centerline points and generating offsets to the helix
    for i in range(npts):                # for each point,
        # Creating a vector perpendicular to the pipe centerline
        # The vector that will be rotated will be perpendicular to the pipe centerline at that point.
        refvec = cross(directions[:, i], nonParallel[:, i])
        # The angle it will be rotated by is the rotation per index times the index
        rotang = i*angrate
        rotax = directions[:, i]
        # Now rotating the vector (applying the 'twist')
        offsets[:, i] = matmul(axang2rotm(append(rotax, rotang)), transpose([refvec]))[:, 0]
        offsets[:, i] = PipeRadius*offsets[:, i]/norm(offsets[:, i])    # Setting the pipe radius
    
    seam = line + offsets       # The helix

    # PLOTTING THE WRAPPED PIPE
    new_window_3d = False
    if ax3d == None:
        # If no axes are given, make new ones in a new figure
        ax3d = plt.figure().add_subplot(projection = "3d")
        new_window_3d = True

    ax3d.clear()
    ax3d.plot(line[0, :], line[1, :], line[2, :], 'g', linewidth = 2)
    ax3d.plot(seam[0, :], seam[1, :], seam[2, :], 'r', linewidth = 0.7)
    
    for i in range(npts-ijump):
        ax3d.plot([seam[0, i], seam[0, i+ijump]], [seam[1, i], seam[1, i+ijump]], [seam[2, i], seam[2, i+ijump]], 'k', linewidth = 0.2)
    
    ax3d.set_title("Wrapped Pipe")

    # Setting equal aspect ratio using helper functions by @Mateen Ulhaq and @karlo
    ax3d.set_box_aspect([1, 1, 1])
    set_axes_equal(ax3d)

    # If a new figure was created, show it
    if new_window_3d:
        plt.show(block = False)


    # UNWRAPPING
    # Making two copies of the helix offset by one wrap.
    # These correspond to the two edges of the wrapper.
    line1_3d = seam[:, 0:npts-ijump]
    line2_3d = seam[:, ijump:npts]
    # line1_2D and line2_2D will store the coordinates of
    # the edges of the wrapper, as unwrapped and drawn in 2D
    line1_2d = zeros((2, npts-ijump))
    line2_2d = zeros((2, npts-ijump))
    # The first point of the first edge of the wrapper can start anywhere.
    line1_2d[:, 0] = [0, 0];
    # The first point of the second edge of the wrapper is at the same distance
    # from the first point of the first edge of the wrapper in 2D as in 3D, and
    # it is at the arbitrary (function parameter) angle PlotAngle from the
    # horizontal in the 2D plot.
    width0 = norm(line1_3d[:, 0] - line2_3d[:, 0])
    line2_2d[:, 0] = [width0*cos(PlotAngle), width0*sin(PlotAngle)]
    # Now iterating through the "rungs of the ladder", i.e., lines joining
    # corresponding points on opposite 3D edges, and adding each new rung to
    # the 2D wrapper profile
    # Each consecutive pair of rungs forms a quadrilateral that is not
    # necessarily planar. To unwrap to 2D, this quadrilateral can be divided
    # into two triangles about either of its two diagonals. Alternating between
    # diagonals from one pair of rungs to the next presumably cancels some
    # cumulative errors.
    for i in range(npts-ijump-1):
        if i%2 == 0:
            # For even indices, divide along diagonal from
            # first-edge-next-rung to second-edge-this-rung
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
            # For odd indices, divide along diagonal from 
            # first-edge-this-rung to second-edge-next-rung 
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

    # DEFINING DOTTED LINE
    dottedline = zeros((2, npts-ijump))            # This will hold the dotted line
    for i in range(npts-ijump):
        # Each point of the dotted line is distance Overlap away from the
        # second line in the "ladder rung" direction.
        dottedline[:, i] = line2_2d[:, i] + Overlap*(line2_2d[:, i] - line1_2d[:, i])/norm(line2_2d[:, i] - line1_2d[:, i])

    # PLOTTING THE 2D WRAPPER PROFILE
    new_window_2d = False
    if ax2d == None:
        # If no axes, are given, make new ones in a new figure
        ax2d = plt.figure().add_subplot()
        new_window_2d = True
    
    ax2d.clear()
    ax2d.plot(line1_2d[0, :], line1_2d[1,:], 'k', linewidth = 0.5)
    ax2d.plot(line2_2d[0, :], line2_2d[1,:], 'k', linewidth = 0.5)
    ax2d.plot(dottedline[0, :], dottedline[1, :], 'k--', linewidth = 0.5)
    
    numturns = int(ceil(npts/ijump))

    for i in range(numturns-1):
        idx = i*ijump
        pt1 = line1_2d[:, idx]
        pt2 = line2_2d[:, idx]
        ax2d.plot([pt1[0], pt2[0]], [pt1[1], pt2[1]], 'k', linewidth = 0.5)
    
    pt1 = line1_2d[:, -1]
    pt2 = line2_2d[:, -1]
    ax2d.plot([pt1[0], pt2[0]], [pt1[1], pt2[1]], 'k', linewidth = 0.5)

    # Equal aspect ratio supported natively for 2D axes
    ax2d.set_aspect("equal")
    ax2d.set_title("Wrapper")

    # If a new figure was created, show it
    if new_window_2d:
        plt.show(block = False)

    # If output filename is given write output data to it
    if OutputFile != None and OutputFile != "":
        df = pd.DataFrame(transpose([line1_2d[0, :], line1_2d[1, :], line2_2d[0, :], line2_2d[1, :], dottedline[0, :], dottedline[1, :]]))
        df.to_excel(OutputFile, header = ["line 1 x", "line 1 y", "line 2 x", "line 2 y", "dotted line x", "dotted line y"], index = False)
        
    return
