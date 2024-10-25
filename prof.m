function [d,z,z_dem,x,y,SW]=prof(DEM,HAR,S,varargin)

% PROF extracting profiles from SWATHobj of the river course
%
% Syntax
%
%    [d,z,z_dem,x,y,SW]=prof(DEM,HAR,S)
%
% Description
%
%     This function allows to extract profile from Relative Elevation
%     Model (REM) or Height Above River (HAR) model as well as from Digital
%     Elevation Model (DEM) 
%
%
% Input arguments
%
%     DEM       high resolution GRIDobj
%     HAR       high resolution GRIDobj
%     S         river course, as a STREAMobj or as a structure array obtained by a
%               shapefile or from an free hand drawn line
%
% Optional Input arguments
%     step      stepping in meters between transversal river profiles
%     width     width in meters of the profiles
%     smooth    smoothing of the trace from which the profiles are
%               extracted
%     plot      plotting the results
%
% Output arguments
%
%     d         an array of the distance of the profiles
%     z         z-values of the profiles (from HAR)
%     z_dem     z-values of the profiles (from DEM)
%     x         x-coordinates of the profiles
%     y         y-coordinates of the profiles
%     SW        SWATHobj of the river profile
%
% See also: SWATHobj/profiles
%
% Author: Michele Delchiaro (michele.delchiaro[at]uniroma1.it)
% Date: 08. November, 2023

% Parse inputs
p = inputParser;
p.FunctionName = 'prof';
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
addRequired(p,'HAR',@(x) isa(x,'GRIDobj'));
addRequired(p,'S',@(x) isa(x,'STREAMobj') || isstruct(S)||ischar(validatestring(x,{'interactive'})));
addParameter(p,'step',1e3,@(x) isnumeric(x))
addParameter(p,'width',1e3,@(x) isnumeric(x))
addParameter(p,'smooth',100,@(x) isnumeric(x))
addParameter(p,'plot',true);

parse(p,DEM,HAR,S,varargin{:});


S=p.Results.S;
plt=p.Results.plot;
step=p.Results.step;
%% extracting profiles

if isa(S,'STREAMobj')

    [X,Y] = STREAMobj2XY(S); % extract x and y coordinates of the river
    
    ix = ~isnan(X); % create a index array to remove nans from the coordinates

    SW = SWATHobj(DEM,X(ix),Y(ix),'width',p.Results.width,'smooth',p.Results.smooth); % compute the swath along the river

    [~,z_dem,~,~] = profiles(SW,'dist','x','step',step,'format','mat'); % extract profiles from DEM-based swath

    SW = mapswath(SW,HAR); % map swath with HAR

    [d,z,x,y] = profiles(SW,'dist','x','step',step,'format','mat'); % extract profiles from HAR-based swath

elseif isstruct(S)

    X=[S.X]; % extract x coordinate of the river

    Y=[S.Y]; % extract y coordinate of the river

    ix = ~isnan(X); % create a index array to remove nans from the coordinates

    SW = SWATHobj(DEM,X(ix),Y(ix),'width',p.Results.width,'smooth',p.Results.smooth); % compute the swath along the river

    [~,z_dem,~,~] = profiles(SW,'dist','x','step',step,'format','mat'); % extract profiles from DEM-based swath

    SW = mapswath(SW,HAR); % map swath with HAR% map swath with HAR

    [d,z,x,y] = profiles(SW,'dist','x','step',step,'format','mat'); % extract profiles from HAR-based swath
    

elseif ischar(S)

    SW= SWATHobj(DEM,'width',p.Results.width); % interactively digitize the river and compute the swath along the river

    [~,z_dem,~,~] = profiles(SW,'dist','x','step',step,'format','mat'); % extract profiles from DEM-based swath  

    SW = mapswath(SW,HAR); % map swath with HAR

    [d,z,x,y] = profiles(SW,'dist','x','step',step,'format','mat'); % extract profiles from HAR-based swath

end

%% plotting profiles
if plt

figure

imageschs(DEM,[],'colormap',[.9 .9 .9],'colorbar',false)
hold on

c=parula(size(z_dem,1));

X = reshape(x,[],size(z_dem,1));
Y = reshape(y,[],size(z_dem,1));

for i =1:size(z_dem,1)
    plot(X(:,i),Y(:,i),'Color',c(i,:))
    hold on
end

plot(SW)

% xlim([3.467e5 3.615e5])
% ylim([4.786e6 4.7925e6])

figure

for i =1:size(z_dem,1)
    plot(d,z_dem(i,:),'Color',c(i,:))
    hold on
end

ylabel('Elevation (m a.s.l.)')
xlabel('Distance from thalweg (m)')

end