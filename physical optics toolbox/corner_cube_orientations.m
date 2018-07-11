function [RotMatrix, cc] = corner_cube_orientations(corner_cube_no,tilt_x,tilt_y,varargin)
% [RotMatrix, ccube_edges] =
% corner_cube_orientations(corner_cube_no,tilt_x,tilt_y,options)
%
% corner_cube_no = 1 or 2, reference:
%%% Corner Cube Orientations
%%% from KimKubes_08-10-2002-#19D08A.xls, Fabien Nicaise, 3-6754
%
% tilt_x, tilt_y = tilt about the x-axis and tilt about the y-axis (rad)
%
% phi_beam = (optional) angle between the metrology beam in the flight
% system and the z-axis. The diffraction code always places the beam on the
% z-axis. The default value is the angle between the z-axis and the
% internal metrology beam according to the latest configuration.
%
% RotMatrix is a rotation matrix that can be used as an input to the corner
% cube object in the diffraction code. The rotation is that required to
% align the default corner cube in the diffraction code with the corner
% cube as oriented in the flight system. tilt_x and tilt_y are applied last
% to model the articulation of the corner cube over the field of view.
%
% ccube_edges is a matrix of three column vectors, where each vector is a
% direction vector indicating the direction of a roof line. The orientation
% is for SIM coordinates with the metrology beam rotated to the z-axis.
% That is, ccube_edges are NOT the input to Benson's code.
%
% options:
%   'simversion' = '2005-05-05'(default), '2004-11-15'
%   'phi_beam'   = override the default value specified by choice of
%       simversion

%normalize to unit vectors
normalize = inline('a./sqrt(a''*a)','a');

%%% constants:
constants;
global P;

%%% check for options
[phitmp, simversion] = CheckOptions(varargin{:});

% surface normal vectors for the flight system from the reference
% each surface normal vector (roof line vector if there is no dihedral) is
% an (x,y,z) column vector
switch simversion,
    %%%%%%% OLD NUMBERS, OBSOLETE 2004-07-30%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     case 3,
%         % corner cube #3 is primary science with NO DELAY LINE
%         cc(:,1) = normalize([
%             -8.42601
%             47.49128
%             46.79125
%             ]);
%         cc(:,2) = normalize([
%             -6.873631
%             -38.82749
%             38.17059
%             ]);
%         cc(:,3) = normalize([
%             54.12936
%             0
%             9.747433
%             ]);
%     case 4,
%         % corner cube #4 is primary science with DELAY LINE
%         cc(:,1) = normalize([
%             -8.42601
%             -47.49128
%             46.79125
%             ]);
%         cc(:,2) = normalize([
%             54.01126
%             0
%             9.726165
%             ]);
%         cc(:,3) = normalize([
%             -6.888661
%             38.91239
%             38.25405
%             ]);
%
%%%%%%%% NEW NUMBERS 2004-07-30%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apparently, these numbers are reversed for #3 and #4
%     case 3,
%         % new corner cube orientation, but don't know which one
%         cc(:,1) = normalize([-0.4355 -0.5741 0.6915]');
%         cc(:,2) = normalize([0.8848 -0.4144 0.2132]');
%         cc(:,3) = normalize([0.1644 0.7056 0.6893]');
%     case 4
%         cc(:,1) = normalize([-0.4355 0.5741 0.6915]');
%         cc(:,2) = normalize([0.1644 -0.7056 0.6893]');
%         cc(:,3) = normalize([0.8848 0.4144 0.2132]');
%
% %%%%%%% NEW NUMBERS 2004-08-16%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     case 3,
%         cc(:,1) = normalize([-0.4875  0.4701 0.7357]');
%         cc(:,2) = normalize([ 0.2666 -0.7223 0.6382]');
%         cc(:,3) = normalize([ 0.8314  0.5073 0.2268]');
%     case 4
%         cc(:,1) = normalize([-0.4875 -0.4701 0.7357]');
%         cc(:,2) = normalize([ 0.8314 -0.5073 0.2268]');
%         cc(:,3) = normalize([ 0.2666  0.7223 0.6382]');
%     PHI_BEAM = 10*P;
%%%%%%%%% NEW NUMBERS 2004-11-15 %%%%%%%%%%%%%%%%%%%%%%%%%
    case '2004-11-15',
        switch corner_cube_no,
            case {1, 3, 'ew', 'hh'},
                cc(:,1) = normalize([ 0.60594  0.76177 0.22920]');
                cc(:,2) = normalize([-0.46851  0.10889 0.87672]');
                cc(:,3) = normalize([ 0.64291 -0.63863 0.42288]');
            case {2, 4, 'ns', 'vv'},
                cc(:,1) = normalize([ 0.60594 -0.76177 0.22920]');
                cc(:,2) = normalize([ 0.64291  0.63863 0.42288]');
                cc(:,3) = normalize([-0.46851 -0.10889 0.87672]');
            otherwise,
                error(['unknown corner cube #: ' num2str(corner_cube_no)]);
        end
        PHI_BEAM = 13*P;
   
%%%%%%%%%%% NEW NUMBERS 2005-05-05 SIM Redesign Arch 1a, file# 27076
    case '2005-05-05',
        switch corner_cube_no,
            case {1, 3, 'ew', 'hh'},
                cc = [
                    0.5341497540369890	-0.5280513438161600	0.6601861999134580
                    0.7929796634220710	0.0422725064340014	-0.6077798027236560
                    0.2930312161890670	0.8481596627367600	0.4413138258031600
                    ];
            case {2, 4, 'ns', 'vv'},
                cc = [
                    0.5341497540369890	0.6601861999134580	-0.5280513438161600
                    -0.7929796634220710	0.6077798027236560	-0.0422725064340014
                    0.2930312161890670	0.4413138258031600	0.8481596627367600
                    ];
            otherwise,
                error(['unknown corner cube #: ' num2str(corner_cube_no)]);
        end
        PHI_BEAM = 9.5*P;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    otherwise,
        error(['unknown SIM version: ' simversion]);
        
end

if ~isempty(phitmp),
    PHI_BEAM = phitmp;
end

%% default corner cube orientation in diffraction program
%ccd = get_benson_ccube_edges;
ccd = [
    0.707106781186548      -0.707106781186548       0.000000000000000E+000
    0.408248290463863       0.408248290463863      -0.816496580927726
    0.577350269189626       0.577350269189626       0.577350269189626
    ];
%%%% step 1: rotate the flight corner cube so that the input beam is along
%%%% the z-axis, and articulate corner cube according to desired tilt_x and
%%%% tilt_y
cc = make_rotmatrix(tilt_x,tilt_y,0) * make_rotmatrix(0,-PHI_BEAM,0) * cc;

%%%% step 2: rotate default orientation to line up edge#1 x-y projection
% find the transformation required to tip/tilt ccd to agree with cc3 or cc4
% first rotate about z to make 1st edge agree in x-y plane projection
phid = atan2(ccd(2,1),ccd(1,1));
phi3 = atan2(cc(2,1), cc(1,1));
phi = phi3 - phid;
% ccdz = get_benson_ccube_edges(phi);
Rz = make_rotmatrix(0,0,phi);
ccdz = Rz * ccd; 

%%%% step 3: rotate about axis normal to z-axis and edge #1 to fully align
%%%% edge #1
% determine angle between the 1st edges, and apply the tilt around the axis
% normal to z and 1st edge
thetad = acos(ccdz(3,1));
theta3 = acos(cc(3,1));
theta = theta3 - thetad;
% axis of rotation:
vrot = normalize(cross(ccdz(:,1),[0 0 1]'));
trot = atan2(vrot(2),vrot(1));
Rt = make_rotmatrix(0,0,trot) * make_rotmatrix(-theta,0,0) * make_rotmatrix(0,0,-trot);
% ccdzz = get_benson_ccube_edges(phi,Rt)
% can also use ccdzz = get_benson_ccube_edges(0,Rt * Rz)
% edge #1 of ccdzz and desired cce are now identical
ccdzz = Rt * ccdz;

%%%% step 4: rotate about edge#1 to align the other edges.
% final rotation is about edge1 to align the other edges. use angle
% between edge #2, use cross product to determine the direction of the
% rotation.
% check that both systems are right handed:
if abs(ccdzz(:,2)'*cc(:,2) - ccdzz(:,3)'*cc(:,3)) < 2e-4,
    psi = acos(ccdzz(:,2)'*cc(:,2)) .*...
        sign(ccdzz(:,1)' * cross(ccdzz(:,2),cc(:,2)));
elseif abs(ccdzz(:,3)'*cc(:,2) - ccdzz(:,2)'*cc(:,2)) < 2e-4,
    warning('not a right-handed system?');
    psi = pi/2 + acos(ccdzz(:,3)'*cc(:,2));
else,
    disp(abs(ccdzz(:,2)'*cc(:,2) - ccdzz(:,3)'*cc(:,3)));
    error('final rotation check failed');
end
% axis of rotation is edge 1:
edge1theta = acos(ccdzz(3,1));
edge1phi   = atan2(ccdzz(2,1),ccdzz(1,1));
Rp = make_rotmatrix(0,0,edge1phi) * make_rotmatrix(0,edge1theta,0) * ...
    make_rotmatrix(0,0,psi) * ...
    make_rotmatrix(0,-edge1theta,0) * make_rotmatrix(0,0,-edge1phi);
% ccdzzz = get_benson_ccube_edges(phi,Rp*Rt)
% ccdzzz = get_benson_ccube_edges(0,Rp*Rt*Rz)
% ccdzzztest = Rp * ccdzz;

RotMatrix = Rp*Rt*Rz;

%%% check:
check = RotMatrix * ccd - cc;
if any( abs(check(:)) > 1e-4 ),
    warning('corner cube orientation: check failed');
end

return

function [phitmp, simversion] = CheckOptions(varargin)

%   'simversion' = '2005-05-05'(default), '2004-11-15'
%   'phi_beam'   = override the default value specified by choice of

% default values:
phitmp = [];
simversion = '2005-05-05';

for ii = 1:length(varargin)
    switch lower(varargin{ii}),
        case 'simversion',
            simversion = varargin{ii+1};
            ii = ii+1;
        case 'phi_beam',
            phitmp = varargin{ii+1};
            ii = ii+1;
        otherwise,
            % do nothing
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check surface normals for consistency
R = cc;

R(:,1)'*R(:,2)
R(:,2)'*R(:,3)
R(:,3)'*R(:,1)

for ii = 1:3,
    disp( cross(R(:,ii),R(:,mod(ii,3)+1)) - R(:, mod(ii+1,3)+1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% draw corner cube faces
face1 = [[0; 0; 0] cc(:,2) cc(:,2)+cc(:,3) cc(:,3)];
face2 = [[0; 0; 0] cc(:,1) cc(:,1)+cc(:,3) cc(:,3)];
face3 = [[0; 0; 0] cc(:,1) cc(:,1)+cc(:,2) cc(:,2)];
figure, fill3(face1(1,:),face1(2,:),face1(3,:),'b',...
    face2(1,:),face2(2,:),face2(3,:),'r',...
    face3(1,:),face3(2,:),face3(3,:),'g'), grid
xlabel('X'), ylabel('Y'), zlabel('Z')
