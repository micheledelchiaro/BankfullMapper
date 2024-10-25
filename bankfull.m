function bank=bankfull(d,z,varargin)

% BANKFULL extracting flow height above thalveg vs. hydraulic depth values
% from the profiles taken from the river course. 
% 
% Syntax
%
%    bank=bankfull(d,z)
%
% Description
%
%     This function allows to extract flow height above thalveg vs. 
%     hydraulic depth values from the profiles obtained by the function PROF.
%     For each profile, the hydraulic depth, defined as the ratio
%     between the flow area and the flow width, is computed for different
%     heights from the thalveg. In detail, the heigths range from 0 to
%     max_depth every 10 cm.
%
%
% Input arguments
%
%     d         an array of the distance of the profiles
%     z         z-values of the profiles (from HAR)
%
% Optional Input arguments
%
%     max_depth maximum height from the thalveg for which the routine is
%               computed
%     plot      plotting the results
%
% Output arguments
%
%     bank      the structure array with elevation above thalveg values
%               (elevation), hydraulic depth values (h), area
%               values (area), and width values (width) for
%               each profile
%
% 
% Author: Michele Delchiaro (michele.delchiaro[at]uniroma1.it)
% Date: 08. November, 2023
%
% Parse inputs
p = inputParser;
p.FunctionName = 'bankfull';
addRequired(p,'d',@(x) isnumeric(x))
addRequired(p,'z',@(x) isnumeric(x))
addParameter(p,'max_depth',6,@(x) isnumeric(x))
addParameter(p,'plot',true);

parse(p,d,z,varargin{:});

%% flow height above thalveg vs. hydraulic depth computation

plt=p.Results.plot;
max_depth=p.Results.max_depth; 
d_elev=0.1:0.1:max_depth;

sz=size(z);

% initialize the output variables

h=zeros(max_depth./0.1,sz(1,1));
a=zeros(max_depth./0.1,sz(1,1));
w=zeros(max_depth./0.1,sz(1,1));

verbose = waitbar(0, 'Please wait...'); % initialize the waitbar

tic

for p=1:sz(1,1)


    x1 = d;

    x2 = [min(d) max(d)];

    y1 = z(p,:);

    xq = linspace(min([x1,x2]), max([x1,x2]),1e6);     % New ,x1 Vector For Interpolation

    y1i = interp1(x1(~isnan(x1)), y1(~isnan(x1)), xq, 'linear','extrap');   % Interpolate To ‘xq’


    for dz=0.1:0.1:max_depth

        y2 = [dz dz];
        y2i = interp1(x2, y2, xq, 'linear','extrap');   % Interpolate To ‘xq’

        [xi,~] = polyxpoly(x1,y1,x2,y2);

        di=diff(xi);

        width=sum(di(1:2:end));
        area=trapz(xq,y2i)-trapz([min(xq) xq(y1i<dz) max(xq)], [dz y1i(y1i<dz) dz]);
        i=round(dz*10);
        h(i,p)=area./width;
        a(i,p)=area;
        w(i,p)=width;


    end

    waitbar(p/sz(1,1), verbose, sprintf('Processing section %d of %d', p, sz(1,1)));
end

close(verbose)

toc



bank.elevation=d_elev;
bank.h=h;
bank.area=a;
bank.width=w;

bank.area(isinf(bank.h))=nan;
bank.width(isinf(bank.h))=nan;
bank.h(isinf(bank.h))=nan;

%% plotting results

if plt

c=parula(sz(1,1));

figure

for p=1:sz(1,1)

    plot(d_elev,h(:,p),'Color',c(p,:))
    hold on
    scatter(d_elev,h(:,p),20,'MarkerFaceColor',c(p,:),'MarkerEdgeColor','none','Marker','o','SizeData',10,'LineWidth',.1)

end

set(gca, 'YScale', 'log')
ylabel('Hydraulic depth (area/width) (m)')
xlabel ('Elevation above thalweg (m)')

end


