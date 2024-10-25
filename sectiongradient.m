function [slope,zw,zc,zs]=sectiongradient(DEM,SW,step,varargin)

% SECTIONGRADIENT computes slope between the section of the centerline used for the bankfull extraction
%
% Syntax
%
%    slope=sectiongradient(DEM,SW,step)
%
% Description
%
%     This function allows to compute the river gradient between the
%     consecutive sections using a downstream elevation correction and a
%     smoothing of the profile applying a moving mean with a user defined
%     window size. Pay attention that the river course is supposed to go
%     downstream.
%
%
% Input arguments
%
%     DEM       high-resolution GRIDobj (used in the prof function)
%     SW        SWATHobj of the river profile (obtained from the prof function)
%     z         z-values of the profiles (from HAR)
%     step      distance in meters between the profiles (used in the prof function)
%
% Optional Input arguments
%
%     windowsize   window size for the moving mean of the corrected
%                  elevation profile
%
%
% Output arguments
%
%     slope      river gradient (with a negative sign)
%     zw         original elevation profile
%     zc         elevation profile after correction
%     zs         elevation profile after smoothing
%
% Author: Michele Delchiaro (michele.delchiaro[at]uniroma1.it)
% Date: 11. July, 2024
%
%
% Parse inputs

p = inputParser;
p.FunctionName = 'sectiongradient';
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
addRequired(p,'SW',@(x) isa(x,'SWATHobj'));
addRequired(p,'step',@(x) isnumeric(x))
addParameter(p,'windowsize',5,@(x) isnumeric(x));
addParameter(p,'plot',true)

parse(p,DEM,SW,step,varargin{:});

step=p.Results.step;
windowsize=p.Results.windowsize;
plt=p.Results.plot;

%%
SWw = SWATHobj(DEM,SW.xy0(:,1),SW.xy0(:,2),'width',0);
[~,zw,~,~] = profiles(SWw,'dist','x','step',step,'format','mat');

% Example elevation data (replace this with your actual data)
zw = zw(:,1)';
zc = zw;

% Ensure the elevation decreases downstream
for i = 2:length(zc)
    if zc(i) > zc(i-1)
        zc(i) = zc(i-1);
    end
end

% Apply a moving average filter to smooth the profile
zs = movmean(zc, windowsize);
slope=diff(zs(1,:)/step);
slope=[NaN,slope];

if plt
    figure
    plot(zw,'k-');

    hold on;

    plot(zc,'r-');
    plot(zs,'b-','LineWidth',1)

    hold off

    ylabel('Elevation (m a.s.l.)')

    yyaxis right
    plot(slope,'g-')

    ylabel('Slope (m/m)')
    xlabel ('Section number')

    legend('Original profile','Corrected profile','Smoothed profile','Slope')
end

end


