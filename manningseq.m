function R=manningseq(SW,z,z_dem,step,d,zvar,varargin)

% MANNINGSEQ computes the Manning's Equation Q = 1/n * (A/P)^(2/3) * S^(1/2)
% * A
%
% Syntax
%
%    R=manningseq(SW,z,z_dem,step,d,H)
%
% Description
%
%     This function allows to compute the river discharge using the
%     Manning's equation
%
% Input arguments
%
%     SW        SWATHobj of the river profile
%     z         z-values of the profiles (from HAR)
%     z_dem     z-values of the profiles (from DEM)
%     step      distance in meters between the profiles
%     d         an array of the distance of the profiles
%     zvar      variable height above thalveg (from DETECT_PEAK function)
%     slope     river gradient
%
% Optional Input arguments
%     n         n coefficient
%     plot      plotting the results
%
% Output arguments
%
% The outputs are organized in a struct file where are organized as follows:
%
%     Q_zvar     river discharge for each section
%     Qi_zvar    most probable river discharge
%     area_zvar  flow area for each section
%     wet_p_zvar wet perimeter for each section
%     width_zvar width for each section
%     slope      river gradient for each section
%
% Author: Michele Delchiaro (michele.delchiaro[at]uniroma1.it)
% Date: 11. July, 2024
%
%
% Parse inputs

p = inputParser;
p.FunctionName = 'manningseq';
addRequired(p,'SW',@(x) isa(x,'SWATHobj'));
addRequired(p,'z',@(x) isnumeric(x))
addRequired(p,'z_dem',@(x) isnumeric(x))
addRequired(p,'step',@(x) isnumeric(x))
addRequired(p,'d',@(x) isnumeric(x))
addRequired(p,'zvar',@(x) isnumeric(x))
addParameter(p,'n',0.05,@(x) isnumeric(x))
addParameter(p,'slope',[],@(x) isnumeric(x));
addParameter(p,'plot',true)

parse(p,SW,z,z_dem,step,d,zvar,varargin{:});

n=p.Results.n;
slope=p.Results.slope;
plt=p.Results.plot;

%%

sz=size(z);

% initialize the variables

area_zvar=zeros(1,sz(1,1));
wet_p_zvar=zeros(1,sz(1,1));
width_zvar=zeros(1,sz(1,1));

verbose = waitbar(0, 'Please wait...'); % initialize the waitbar

tic

for p=1:sz(1,1)

    x1 = d;

    x2 = [min(d) max(d)];

    y1 = z(p,:);

    if isnan(zvar(1,p))
        % Skip when zvar(1,p) is NaN
        continue;
    else

        y2 = [zvar(1,p) zvar(1,p)];
        % make common x coordinates

        xq = linspace(min([x1,x2]), max([x1,x2]),1e6);     % New ,x1 and x2 Vector For Interpolation

        y1i = interp1(x1(~isnan(x1)), y1(~isnan(x1)), xq, 'linear','extrap');   % Interpolate To ‘xq’

        y2i = interp1(x2, y2, xq, 'linear','extrap');   % Interpolate To ‘xq’


        % Use polyxpoly to find intersection points for wet perimeter calculus
        [intersection_x, ~] = polyxpoly(x1, y1, x2, y2);

        % Combine x and y coordinates of both lines
        all_x_coordinates = [xq, x2];
        all_y_coordinates = [y1i, y2];

        % Initialize sum of lengths
        length_sum = 0;



        % Process segments below intersections
        for i = 1:length(intersection_x)-1
            % Find indices corresponding to the current segment
            indices = find(all_x_coordinates >= intersection_x(i) & all_x_coordinates <= intersection_x(i+1));

            % Calculate length of the current segment only if it's below the line2
            for j = 1:length(indices)-1
                if all_y_coordinates(indices(j)) <= y2(1)
                    segment_length = sqrt((all_x_coordinates(indices(j+1)) - all_x_coordinates(indices(j)))^2 + ...
                        (all_y_coordinates(indices(j+1)) - all_y_coordinates(indices(j)))^2);
                    length_sum = length_sum + segment_length;
                end
            end
        end


        wet_p_zvar(1,p)=length_sum;

        % area calculus
        di=diff(intersection_x);
        width_zvar(1,p)=sum(di(1:2:end));

        area_zvar(1,p)=trapz(xq,y2i)-trapz([min(xq) xq(y1i<zvar(1,p)) max(xq)], [zvar(1,p) y1i(y1i<zvar(1,p)) zvar(1,p)]);

    end



    waitbar(p/sz(1,1), verbose, sprintf('Processing section %d of %d', p, sz(1,1)));

end

close(verbose)

toc

if isempty(slope)

    slope=abs((SW.zd0(end,1)-SW.zd0(1,1))/SW.zd0(end,2));

end

Q_zvar=1/n.*((area_zvar/wet_p_zvar).^(2/3)).*(abs(slope).^(1/2)).*area_zvar;

Q_zvar(Q_zvar==0)=NaN;

[f, xi] = ksdensity(Q_zvar, 'Bandwidth', 10);

Qi_zvar=xi(f==max(f));

% make struct file of the results

R.Q_zvar=Q_zvar;
R.Qi_zvar=Qi_zvar;

R.area_zvar=area_zvar;
R.wet_p_zvar=wet_p_zvar;

R.width_zvar=width_zvar;

R.slope=slope;

if plt

    % plotting histogram Q
    figure

    histogram(Q_zvar,100)

    ylabel('Count')
    xlabel('River Discharge (m^{3} s^{-1})')

    yyaxis right

    plot(xi,f,'r-','LineWidth',1)
    ylabel ('Density probability')


    xline(Qi_zvar,':','LineWidth',1);
    

end

end

