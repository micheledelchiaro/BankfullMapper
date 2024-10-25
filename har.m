function REM=har(DEM,S)

% HAR computing Height Above River (HAR) or Relative Elevation Model (REM)
%
% Syntax
%
%    REM=har(DEM,S)
%
% Description
%
%     This function allows to compute the Relative Elevation Model (REM) or
%     Height Above River (HAR) model of a river
%
%
% Input arguments
%
%     DEM       high resolution GRIDobj
%     S         river course, as a STREAMobj or as a structure array obtained by a
%               shapefile
%
% Output arguments
%
%     HAR       high resolution GRIDobj
%
% See also: GRIDobj/interp2GRIDobj
%
% Author: Michele Delchiaro (michele.delchiaro[at]uniroma1.it)
% Date: 08. November, 2023

% Parse inputs
p = inputParser;
p.FunctionName = 'har';
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
addRequired(p,'S',@(x) isa(x,'STREAMobj')|| isstruct(S));

%% computing REM
if isa(S,'STREAMobj')
    trunk_river=trunk(klargestconncomps(S,1)); % extract trunk river
    [x,y,z] = STREAMobj2XY(trunk_river,DEM); % extract x,y,z coordinates of the trunk river
    z(isnan(z))=[];
    x(isnan(x))=[];
    y(isnan(y))=[];
    A = interp2GRIDobj(DEM,x,y,z,'linear','nearest'); % interpolate a plane based on the coordinates
    REM=DEM-A; % remove the plane from the DEM

elseif isstruct(S)

    x=[S.X]; % extract x coordinate
    y=[S.Y]; % extract y coordinate
    x(isnan(x))=[];
    y(isnan(y))=[];
    z = double(interp(DEM,x,y)); % extract z coordinate from x and y coordinates
    A = interp2GRIDobj(DEM,x',y',z','linear','nearest'); % interpolate a plane based on the coordinates
    REM=DEM-A;% remove the plane from the DEM

end
