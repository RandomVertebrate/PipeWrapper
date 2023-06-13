function WrapperProfile2D = WrapPipe(ProfileFile, PipeRadius, TurnsPerMeter, Overlap, Resolution, PlotAngle, normvec)

% ProfileFile is the neame of a file that contains 'x', 'y', and 'z'
% columns containing coordinates defining the centerline of the pipe.
% PipeRadius is the radius of the pipe in meters.
% TurnsPerMeter is the number of turns of wrapping per meter of pipe.
% Overlap is the overlap between successive wraps/turns. A dashed line is
% added to the output to mark it.
% Resolution is the spatial resolution (in meters) of the calculations.
% Convergence of the output should be checked by considering progressively
% smaller Resolution values.
% PlotAngle rotates the output.
% normvec is used as a reference for the helix angle.  It must be chosen 
% such that the pipe centerline is never parallel to it.


% READING PIPE CENTERLINE
% Expected contents of file:
% Three columns called 'x', 'y' and 'z' with coordinates of points along pipe centerline
Table = readtable(ProfileFile); % Read file contents
xin = Table.x';                 % Row vetor of x-coordinates
yin = Table.y';                 % Row vetor of y-coordinates
zin = Table.z';                 % Row vetor of z-coordinates
maxpts = length(xin);           % Number of points read

linein = [xin; yin; zin];       % Data read from file. Each column is a point [x; y; z]

% GENERATING EVENLY SPACED INTERPOLATION OF PIPE CENTERLINE
% We need the centerline to consist of (physically) evenly spaced points so
% that we can use an array index to determine how far along the line we are.
line = linein(:, 1);    % Generated interpolation will be stored here. The starting point is unchanged.
% partition keeps track of how much of the input data has been processed.
% (All up to the (partition-1)'th point have been processed and
% the partition'th point is currently being processed)
partition = 1;
% This loop goes through the input data, generates evenly spaced points, and adds them to the output data.
while norm(linein(:, end) - line(:, end)) > Resolution        % Keep going until the last generated point is within Resolution distance of the last input point
    for j = partition:maxpts                                  % Iterate through the unprocessed input data
        % If this point were added as-is to the generated data, that would
        % correspond to a length increment of
        potentialIncrement = linein(:, j) - line(:, end);
        % Wait until the potential increment is longer than Resolution
        if norm(potentialIncrement) > Resolution
            % Now add an increment in this direction but of magnitude Resolution
            incrementDirection = potentialIncrement/norm(potentialIncrement); 
            actualIncrement = incrementDirection*Resolution;
            line(:, end+1) = line(:, end) + actualIncrement;  % New point added to generated data
            % Update partition. Note: partition'th point will be checked
            % again now, just not points to its left. This means that
            % further increments towards the partition'th point are
            % possible if necessary.
            partition = j;
            break
        end
    end
end

angrate = TurnsPerMeter*2*pi*Resolution;        % Converting wrap rate to radians per point

ijump = round(2*pi/angrate);                    % Index difference between subsequent wraps

npts = length(line);                            % Number of points in the generated evenly spaced pipe centerline

directions = diff(line, 1, 2);                  % Array of vectors along pipe centerline
directions = [directions, directions(:, end)];  % Duplicating the last vector so that directions has length npts
% Normalizing directions
for i = 1:npts
    directions(:, i) = directions(:, i)/norm(directions(:, i));
end

% GENERATING A 'HELIX' THAT FOLLOWS THE PIPE CENTERLINE
offsets = zeros(3, npts);                       % This will store offset vectors between the pipe centerline and a helix that follows it
% Iterating through centerline points and generating offsets to the helix
for i = 1:npts                                                       % For each point,
    refvec = cross(directions(:, i), normvec);                       % The vector that will be rotated will be perpendicular to the pipe centerline at that point.
    rotang = (i-1)*angrate;                                          % The angle it will be rotated by is the rotation per index times the index
    rotax = directions(:, i);                                        % The axis of rotation is the pipe centerline direction
    Smat = [0, -rotax(3), rotax(2); rotax(3), 0, -rotax(1); -rotax(2), rotax(1), 0];   % Skew-symmetric matrix
    Rmat = cos(rotang)*eye(3) + (1-cos(rotang))*rotax*rotax' + sin(rotang)*Smat;       % Rotation matrix formula
    offsets(:, i) = Rmat*refvec;
    offsets(:, i) = PipeRadius*offsets(:, i)/norm(offsets(:, i));    % Setting the pipe radius
end

seam = line + offsets;       % The helix

% PLOTTING THE WRAPPED PIPE
views = [[0; 0; 1], [1; 1; 1], [0; 1; 0], [1; 0; 0]];
figure()
for curplot = 1:4
    subplot(2, 2, curplot);
    plot3(linein(1, :), linein(2, :), linein(3, :), 'b')
    hold on
    plot3(line(1, :), line(2, :), line(3, :), 'g', 'LineWidth', 5)
    plot3(seam(1, :), seam(2, :), seam(3, :), 'r', 'LineWidth', 5)
    for i = 1:npts-ijump
        plot3([seam(1, i), seam(1, i+ijump)], [seam(2, i), seam(2, i+ijump)], [seam(3, i), seam(3, i+ijump)], 'k')
    end
    view(views(:, curplot))
    axis equal
end
sgtitle("Wrapped Pipe")

% UNWRAPPING
% Making two copies of the helix offset by one wrap.
% These correspond to the two edges of the wrapper
line1_3d = seam(:, 1:npts-ijump);
line2_3d = seam(:, ijump:npts);
% line1_2D and line2_2D will store the coordinates of the edges of the wrapper, as unwrapped and drawn in 2D
line1_2d = zeros(2, npts-ijump);
line2_2d = zeros(2, npts-ijump);
% The first point of the first edge of the wrapper can start anywhere.
line1_2d(:, 1) = [0; 0];
% The first point of the second edge of the wrapper is at the same distance
% from the first point of the first edge of the wrapper in 2D as in 3D, and
% it is at the arbitrary (function parameter) angle PlotAngle from the
% horizontal in the 2D plot.
width0 = norm(line1_3d(:, 1)-line2_3d(:, 1));
line2_2d(:, 1) = [width0*cos(PlotAngle); width0*sin(PlotAngle)];
% Now iterating through the "rungs of the ladder", i.e., lines joining
% corresponding points on opposite 3D edges, and adding each new rung to
% the 2D wrapper profile
for i = 1:npts-ijump-1
    % The ends of the current rung form a triangle ABC with one end of the
    % next rung
    A3d = line1_3d(:, i);
    B3d = line2_3d(:, i);
    C3d = line1_3d(:, i+1);
    % We already have the current rung in 2D, so A and B are known in 2D
    A2d = line1_2d(:, i);
    B2d = line2_2d(:, i);
    % Finding C in 2D
    C2d = Triangle3dTo2d(A3d, B3d, C3d, A2d, B2d);
    line1_2d(:, i+1) = C2d;           % Now we know one end of the next rung in 2D
    % Repeat for the other end
    C3d = line2_3d(:, i+1);
    C2d = Triangle3dTo2d(A3d, B3d, C3d, A2d, B2d);
    line2_2d(:, i+1) = C2d;
end

% DEFINING DOTTED LINE
dottedline = zeros(2, npts-ijump);            % This will hold the dotted line
for i = 1:npts-ijump
    % Each point of the dotted line is distance Overlap away from the
    % second line in the "ladder rung" direction.
    dottedline(:, i) = line1_2d(:, i) + Overlap*(line1_2d(:, i) - line2_2d(:, i))/norm(line1_2d(:, i) - line2_2d(:, i));
end

% SETTING RETURN TABLE
WrapperProfile2D = array2table([line1_2d' line2_2d' dottedline'], 'VariableNames', ["Edge1_x", "Edge1_y", "Edge2_x", "Edge2_y", "Overlap_x", "Overlap_y"]);

% PLOTTING THE 2D WRAPPER PROFILE
figure();
clf
plot(line1_2d(1, :), line1_2d(2, :), 'k')
hold on
plot(line2_2d(1, :), line2_2d(2, :), 'k')
plot(dottedline(1, :), dottedline(2, :), '--k')
numturns = ceil(npts/ijump);
pt1 = line1_2d(:, 1);
pt2 = line2_2d(:, 1);
plot([pt1(1), pt2(1)], [pt1(2), pt2(2)], 'k')
for i = 0:numturns-3
    idx = floor((i+1)*ijump);
    pt1 = line1_2d(:, idx);
    pt2 = line2_2d(:, idx);
    plot([pt1(1), pt2(1)], [pt1(2), pt2(2)], 'k')
end
pt1 = line1_2d(:, end);
pt2 = line2_2d(:, end);
plot([pt1(1), pt2(1)], [pt1(2), pt2(2)], 'k')
title('Wrapper')
xlabel('Scale (meters)')
axis image
set(gcf, 'Units', 'centimeters');
end
