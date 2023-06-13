function pipeline = GenPipeline(outputfile, cornersx, cornersy, cornersz, bendradii, bendpoints)
% Given a sequence of points in 3D space, this function connects them,
% rounds the corners, and returns the resulting curve (discretized, so a
% set of points)

% outputfile is the file to write the resulting data to
% cornersx, cornersy, cornersz contain the coordinates of the input points
% bendradii contains the radii of the bends
% bendpoints is the resolution (number of subdivisions) of each bend

numbends = length(cornersx) - 2;                       % Number of bends is 2 less than number of points (all points are bends except first and last)

corners = [cornersx; cornersy; cornersz];              % Each column is an input point

pipeline = [cornersx(1); cornersy(1); cornersz(1)];    % This will store the output. Each column will be a point. Currently only contains the first point, which will be unchanged.

for i = 1:numbends                                     % For each bend
    
    thisbend = zeros(3, bendpoints);                   % Initialize array to store points for this bend
    
    % Unit vector parallel to pipe before bend
    oldDirection = (corners(:, i+1) - corners(:, i))/norm(corners(:, i+1) - corners(:, i));
    % Unit vector parallel to pipe after bend
    newDirection = (corners(:, i+2) - corners(:, i+1))/norm(corners(:, i+2) - corners(:, i+1));
    % Corner angle (if there was no bend) is the inverse cosine of the dot
    % product of the unit vectors
    cornerangle = acos(abs(dot(newDirection, oldDirection)));
    % The angle through which an arc will be drawn. Equal to pi-cornerangle
    % for geometry reasons
    bendAngle = pi-cornerangle;
    % The "arc angle" between each pair of subsequent points in the bend.
    AngInc = bendAngle/(bendpoints-1);
    % The bend axis is perpendicular to the pipe before and after the bend.
    bendAxis = cross(oldDirection, newDirection);
    bendAxis = bendAxis/norm(bendAxis);
    % Vector pointing from center of the arc to the start of the bend (end
    % of the straight section of pipe)
    startradvec = cross(oldDirection, bendAxis)*bendradii(i);
    % This gives the position vector of the center of the arc.
    bendCenter = corners(:, i+1) ...
        - oldDirection*bendradii(i)/tan((cornerangle)/2) ...
        - startradvec;
    % Rotate startradvec about bendAxis by multiples of AngInc and add to 
    % bendCenter to get the points along the bend
    for j = 1:bendpoints
        thisbend(:, j) = bendCenter + axang2rotm([bendAxis(1), bendAxis(2), bendAxis(3), AngInc*(j-1)])*startradvec;
    end
    % Add (concatenate) the calculated bend to the output profile
    pipeline = [pipeline, thisbend];
end

% Add the last point, which is unchanged.
pipeline = [pipeline, [cornersx(end); cornersy(end); cornersz(end)]];

% PLOT PIPE PROFILE
figure()
subplot(2, 1, 1)
plot3(cornersx, cornersy, cornersz);
title('Unrounded');
axis equal

subplot(2, 1, 2);
plot3(pipeline(1, :), pipeline(2, :), pipeline(3, :))
title('Rounded');
axis equal

% WRITE OUTPUT DATA PIPE LINE POINTS
writematrix(pipeline', outputfile, 'Range', 'A2')

end