function C2d = Triangle3dTo2d(A3d, B3d, C3d, A2d, B2d)

% Gives the third point in a 2D drawing of a triangle ABC defined in 3D space.
% (Transforms, rather than projects, the triangle).
% 3D points are column vectors [x; y; z] and 2D points are [x; y].

BaseVec3d = B3d - A3d;                                         % 3D vector AB
BaseDir3d = BaseVec3d/norm(BaseVec3d);                         % Direction of 3D AB
BasePt3d = A3d + BaseDir3d*dot((C3d - A3d), BaseDir3d);        % Foot of perpendicular line dropped from C onto AB
height = norm(C3d - BasePt3d);                                 % Height of perpendicular line
footdist = norm(BasePt3d - A3d);                               % Distance of foot from A

% Defining the 2D points in 3D space
A2d = [A2d; 0];
B2d = [B2d; 0];

BaseVec2d = B2d - A2d;                                         % 2D vector AB
BaseDir2d = BaseVec2d/norm(BaseVec2d);                         % Direction of 2D AB

HeightDir2d = cross(BaseDir2d, [0; 0; 1]);                     % Vector perpendicular to AB

C2d = A2d + footdist*BaseDir2d + height*HeightDir2d;           % Finding 2D C

C2d = C2d(1:2);                                                % Removing the (zero) z component of C

end
