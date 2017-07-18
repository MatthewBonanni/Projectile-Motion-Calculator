function [r, range, maxPoint, ToF] = projMotion(y0, v0, theta0, varargin)
% projMotion Determine flight paths for projectiles
%   
%   Uses the Runge-Kutta method to solve ordinary differential equations of
%   motion for projectile flight trajectories.
%   
%   Required inputs:
%       y0          Initial height of the projectile
%       v0          Scalar initial velocity of the projectile
%       theta0      Launch angle, measured from horizontal
%   
%   [r, range, maxPoint, ToF] = projMotion(y0, v0, theta0) returns an array
%   of x and y coordinates r corresponding to the trajectory of the
%   flight, as well as the following data:
%       
%       range       Maximum x displacement of the projectile
%       maxPoint    [x, y] vector of the coordinates of the apex of the
%                   projectile's flight path
%       ToF         Total time of flight of the projectile
%   
%   [r, range, maxPoint, ToF] = projMotion(y0, v0, theta0, g) calculates
%   flight path using specified gravitational acceleration g. (Default: g =
%   9.807 m/s^2)
%   
%   [r, range, maxPoint, ToF] = projMotion(y0, v0, theta0, g, 'dragType', <>, 'mediumProps', <>, 'objectProps', <>)
%   
%   'dragType' Optional parameter. Takes string corresponding to air drag
%   model to be applied
%       
%       'none'      No drag (default)
%       'lin'       Linear (viscous) drag model
%       'quad'      Quadratic (Newtonian) drag model
%       'linquad'   Linear and Quadratic drag models
%   
%   'mediumProps' Optional parameter. Takes structure with physical
%   constants for the medium the projectile is traveling through. Defaults
%   to values for air at standard temperature and pressure. Expects the
%   following fields:
%   
%       beta        Linear drag coefficient (default 0.00016)
%       gamma       Quadratic drag coefficient (default 0.25)
%   
%   If dragType is not 'none', then the 'objectProps' parameter is
%   required:
%   
%   'objectProps' Takes structure with values for physical properties of
%   projectile. (Note: projectile is assumed to be spherical) Expects the
%   following fields:
%   
%       mass        Mass of the projectile
%       diameter    Diameter of the projectile's cross section normal to
%                   direction of motion
% 
% Matthew R. Bonanni
% 6-10-2016

%% Manage inputs
p = inputParser;

expectedDragTypes = {'none', 'lin', 'quad', 'linquad'};
defaultDragType = 'none';
gDefault = 9.807;

defaultMediumProps = struct('beta', 0.00016, 'gamma', 0.25);

addRequired(p, 'y0', @isnumeric);
addRequired(p, 'v0', @isnumeric);
addRequired(p, 'theta0', @isnumeric);
addOptional(p, 'g', gDefault, @isnumeric);
addParameter(p, 'dragType', defaultDragType, @(x) any(validatestring(x,expectedDragTypes)));
addParameter(p, 'mediumProps', defaultMediumProps, @(x) isstruct(x) && isfield(x, 'beta') && isfield(x, 'gamma'));
addParameter(p, 'objectProps', struct(), @isstruct); 

parse(p, y0, v0, theta0, varargin{:});
inputs = p.Results;

if ~strcmp(inputs.dragType, 'none') && (~isfield(inputs.objectProps, 'mass') || ~isfield(inputs.objectProps, 'diameter'))
    error('objectProps structure must exist and have fields mass, diameter');
end

%% Get functions and fit data

%Calculate initial conditions for differential equations
v0components = [inputs.v0*cosd(inputs.theta0), v0*sind(inputs.theta0)];

%Calculate necessary drag coefficients based on drag type
constants.g = inputs.g;
switch(inputs.dragType)
    case 'none'
        constants.b = 0;
        constants.c = 0;
        constants.mass = 1;
        
    case 'lin'
        constants.b = inputs.mediumProps.beta * inputs.objectProps.diameter;
        constants.c = 0;
        constants.mass = inputs.objectProps.mass;
        
    case 'quad'
        constants.b = 0;
        constants.c = inputs.mediumProps.gamma * inputs.objectProps.diameter^2;
        constants.mass = inputs.objectProps.mass;
        
    case 'linquad'
        constants.b = inputs.mediumProps.beta * inputs.objectProps.diameter;
        constants.c = inputs.mediumProps.gamma * inputs.objectProps.diameter^2;
        constants.mass = inputs.objectProps.mass;
end

%Approximate time span based on no air resistance
tEndApprox = 2 * (inputs.v0*sind(inputs.theta0) + sqrt((inputs.v0*sind(inputs.theta0))^2 + 2*inputs.g*inputs.y0)) / inputs.g;
if tEndApprox == 0
    error('Invalid input: no distance traveled');
end
tSpan = linspace(0,tEndApprox,1000);

%Solve velocity differential equations
[t, v] = ode45(@motionEq, tSpan, v0components, [], constants);

%Get Vx(t)
[xData, yData] = prepareCurveData(t, v(:,1));
fitType = 'splineinterp';
fittedVx_vs_t = fit( xData, yData, fitType, 'Normalize', 'on' );

%Get Vy(t)
[xData, yData] = prepareCurveData(t, v(:,2));
fitType = 'splineinterp';
fittedVy_vs_t = fit( xData, yData, fitType, 'Normalize', 'on' );

%Calculate acceleration and position based on velocity
a = [differentiate(fittedVx_vs_t, tSpan), differentiate(fittedVy_vs_t, tSpan)]; %#ok<NASGU>
r = [integrate(fittedVx_vs_t, tSpan, 0)', integrate(fittedVy_vs_t, tSpan, 0)'];
r(:,2) = r(:,2) + inputs.y0;

%Get trajectory function (negated for finding minimum)
[xData, yData] = prepareCurveData(r(:,1), r(:,2));
fitType = 'splineinterp';
fittedTrajectory = fit( xData, -yData, fitType, 'Normalize', 'on' );

%% Get important points

rPositive = r(r(:,2) >= 0, :);

%Find minimum, set to max height
[xMax, yMax] = fminbnd(fittedTrajectory, r(1,1), r(end,1));
yMax = -yMax;
maxPoint = [xMax yMax];

%Find zero of trajectory, set to range
range = fzero(fittedTrajectory, rPositive(end,1));

%Get height as function of time
[xData, yData] = prepareCurveData(t, r(:,2));
fitType = 'splineinterp';
fittedYvsTime = fit( xData, yData, fitType, 'Normalize', 'on' );

%Find zero of height, set to flight time
ToF = fzero(fittedYvsTime, t(length(rPositive)));

%% Plot results
if nargout == 0
    figure(1);
    plot(linspace(0,range,1000), -fittedTrajectory(linspace(0,range,1000)));
    xlabel('x');
    ylabel('y');
    axis equal;
    hold on
    plot([r(1,1) r(end,1)], [0 0], 'r--');
    plot(maxPoint(1), maxPoint(2), 'g*');
    plot(range, -fittedTrajectory(range), 'g*');
end