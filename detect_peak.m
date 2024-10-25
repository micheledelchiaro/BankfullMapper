function lim=detect_peak(h,a,w,d_elev,d,z,varargin)

% DETECT_PEAK detecting peaks from the flow height above thalveg vs.
% hydraulic depth curves of each profile
%
% Syntax
%
%    lim=detect_peak(h,a,w,d_elev,d,z)
%
% Description
%
%     This function allows to detect the peaks from the flow heigh above thalveg vs. hydraulic depth
%     hydraulic depth curves of each profiles. It is possible to detect the
%     all the peaks (using 'all' option), most prominent peaks (using 'max' option),
%     and lowest peaks (using 'lowest' option). For the 'max' and 'lowest'
%     option, the most probable Z value is also computed.
%
%
% Input arguments
%
%     h         hydraulic depth values
%     a         flow area values
%     w         flow width values
%     d_elev    flow height values
%
% Optional Input arguments
%     d         an array of the distance of the profiles (from prof
%               function). It is necessary if 'peak' option is 'max' or 'lowest'
%     z         z-values of the profiles (from HAR). It is necessary if 'peak' option is 'max' or 'lowest'
%     'peak'    string
%        'all': all the peaks are extracted
%        'max': the most prominent peaks are extracted
%        'lowest': the lowest peaks are extracted
%      plot      plotting the results
%
% Output arguments
%
%     lim       the structure array with section numbers (sections field),
%               peak elevation above thalveg values (zlim), the most probable 
%               peak elevation above thalveg for each section based on the 
%               2d kernel distribution estimation of zlim values (zvar), peak hydraulic
%               depth values (hlim), peak hydraulic depth values using zvar 
%               (hvar only if 'peak' option is 'max' or 'lowest'), flow area values (alim),
%               flow area values using zvar (avar only if 'peak' option is 'max' or 'lowest'), 
%               flow width values (wlim), flow width values using zvar (wvar only if 'peak' option is 'max' or 'lowest'). 
%               For 'max' and 'lowest' peak detection, the most probable elevation above thalveg
%               value of the peaks is also reported (Z).
%
%
% Author: Michele Delchiaro (michele.delchiaro[at]uniroma1.it)
% Date: 11. July, 2024
%
% Parse inputs
p = inputParser;
p.FunctionName = 'detect_peak';

addRequired(p,'h',@(x) isnumeric(x))
addRequired(p,'a',@(x) isnumeric(x))
addRequired(p,'w',@(x) isnumeric(x))
addRequired(p,'d_elev',@(x) isnumeric(x))
addOptional(p,'d',@(x) isnumeric(x))
addOptional(p,'z',@(x) isnumeric(x))
addParameter(p,'peak','all',@(x) ischar(validatestring(x,{'all','max','lowest'})))
addParameter(p,'plot',true)

parse(p,h,a,w,d_elev,d,z,varargin{:});

%% peak detection
plt=p.Results.plot;

meth = validatestring(p.Results.peak,{'all','max','lowest'});

sz=size(h);

switch meth

    case 'all'

        % detect zlim and hlim values

        zl = cell(sz(1,2), 1);
        hl = cell(sz(1,2), 1);


        % Store the peak elevation and hydraulic depth values in zlim
        % and hlim vectors

        zlim = [];
        hlim = [];

        % Store the sections numbers in the sections vector
        sections = [];

        for p = 1:sz(1,2)

            [hl{p,1}, zl{p,1},~,~] = findpeaks(h(:,p),d_elev);

            if isnumeric(zl{p})
                % Concatenate the values to the zlim and hlim arrays

                zlim = [zlim, zl{p, 1}];

                hlim = [hlim; hl{p, 1}];

                % Store the sections numbers
                sections = [sections, repmat(p, 1, numel(zl{p, 1}))];

            end


        end

        % detect alim and wlim values

        alim=[];
        wlim=[];

        zlim(isinf(hlim))=[];
        sections(isinf(hlim))=[];
        hlim(isinf(hlim))=[];

        for p=1:sz(1,2)

            ht=h(:,p);
            hp=hlim(sections==p);
            wt=w(:,p);
            at=a(:,p);



            for m=1:numel(hp)

                A=at(ht==hp(m));
                W=wt(ht==hp(m));
                alim=[alim,A];
                wlim=[wlim,W];


            end

        end


        lim.sections=sections;
        lim.zlim=zlim;
        lim.hlim=hlim';
        lim.alim=alim;
        lim.wlim=wlim;



    case 'max'

        sections=1:1:sz(1,2);
        zlim=zeros(1,sz(1,2));
        hlim=zeros(1,sz(1,2));
        h_var=zeros(1,sz(1,2));
        alim=zeros(1,sz(1,2));
        area_zvar=zeros(1,sz(1,2));
        width_zvar=zeros(1,sz(1,2));
        wlim=zeros(1,sz(1,2));



        for p=1:sz(1,2)
            wt=w(:,p);
            at=a(:,p);

            % Find the indices of the local maxima
            [hl, zl,~,pr] = findpeaks(h(:,p),d_elev,'SortStr','descend');

            max_pr=max(pr(~isinf(pr)));

            if isempty(zl)
                zlim(1,p)=nan;
                hlim(1,p)=nan;
            else
                zlim(1,p)=zl(pr==max_pr);
                hlim(1,p)=hl(zl==zlim(1,p));

            end

            if isnan(hlim(1,p))||isinf(hlim(1,p))
                alim(1,p)=NaN;
                wlim(1,p)=NaN;

            else

                alim(1,p)=at(h(:,p)==hlim(1,p));
                wlim(1,p)=wt(h(:,p)==hlim(1,p));
            end


        end

        [f, xi] = ksdensity(zlim, 'Bandwidth', 0.2);
        Z=xi(f==max(f));

        % computing the zvar from the most probable kernel distribution
        n_points = 100;
        [xi, yi] = meshgrid(linspace(min(sections), max(sections), n_points), linspace(min(zlim), max(zlim), n_points));
        xy = [sections', zlim'];
        [f,~] = ksdensity(xy, [xi(:), yi(:)]);
        f = reshape(f, n_points, n_points);

        % Initialize array to store mode zvar for each section
        zvar = zeros(n_points, 1);

        for i = 1:n_points
            % Find the index of the maximum density value for the current section
            [~, max_idx] = max(f(:, i));

            % Get the corresponding zvar value
            zvar(i) = yi(max_idx, i);
        end

        % Interpolate mode_zvar values for the full length of sections
        zvar = interp1(linspace(min(sections), max(sections), n_points), zvar, sections, 'linear');

        %%
        for p=1:sz(1,2)

            x1 = d;

            x2 = [min(d) max(d)];

            y1 = z(p,:);

            y2 = [zvar(1,p) zvar(1,p)];

            % make common x coordinates

            xq = linspace(min([x1,x2]), max([x1,x2]),1e6);     % New ,x1 and x2 Vector For Interpolation

            y1i = interp1(x1(~isnan(x1)), y1(~isnan(x1)), xq, 'linear','extrap');   % Interpolate To ‘xq’

            y2i = interp1(x2, y2, xq, 'linear','extrap');   % Interpolate To ‘xq’


            % Use polyxpoly to find intersection points for wet perimeter calculus
            [intersection_x, ~] = polyxpoly(x1, y1, x2, y2);


            di=diff(intersection_x);
            width_zvar(1,p)=sum(di(1:2:end));

            area_zvar(1,p)=trapz(xq,y2i)-trapz([min(xq) xq(y1i<zvar(1,p)) max(xq)], [zvar(1,p) y1i(y1i<zvar(1,p)) zvar(1,p)]);
            h_var(1,p)=area_zvar(1,p)./width_zvar(1,p);

        end

        %%
        lim.sections=sections;
        lim.zlim=zlim;
        lim.zvar=zvar;
        lim.hlim=hlim;
        lim.hvar=h_var;
        lim.alim=alim;
        lim.avar=area_zvar;
        lim.wlim=wlim;
        lim.wvar=width_zvar;
        lim.Z=Z;





    case 'lowest'

        sections=1:1:sz(1,2);
        zlim=zeros(1,sz(1,2));
        hlim=zeros(1,sz(1,2));
        h_var=zeros(1,sz(1,2));
        alim=zeros(1,sz(1,2));
        area_zvar=zeros(1,sz(1,2));
        width_zvar=zeros(1,sz(1,2));
        wlim=zeros(1,sz(1,2));

        for p=1:sz(1,2)
            wt=w(:,p);
            at=a(:,p);
            % Find the indices of the local maxima
            [hl, zl,~,~] = findpeaks(h(:,p),d_elev,'SortStr','descend');

            if isempty(zl)
                zlim(1,p)=nan;
                hlim(1,p)=nan;
            else
                zlim(1,p)=min(zl);
                hlim(1,p)=hl(zl==min(zl));

            end

            if isnan(hlim(1,p))
                alim(1,p)=NaN;
                wlim(1,p)=NaN;

            else

                alim(1,p)=at(h(:,p)==hlim(1,p));
                wlim(1,p)=wt(h(:,p)==hlim(1,p));
            end


        end

        [f, xi] = ksdensity(zlim, 'Bandwidth', 0.2);
        Z=xi(f==max(f));

        % computing the zvar from the most probable kernel distribution
        n_points = 400;
        [xi, yi] = meshgrid(linspace(min(sections), max(sections), n_points), linspace(min(zlim), max(zlim), n_points));
        xy = [sections', zlim'];
        [f,~] = ksdensity(xy, [xi(:), yi(:)]);
        f = reshape(f, n_points, n_points);

        % Initialize array to store mode zvar for each section
        zvar = zeros(n_points, 1);

        for i = 1:n_points
            % Find the index of the maximum density value for the current section
            [~, max_idx] = max(f(:, i));

            % Get the corresponding zvar value
            zvar(i) = yi(max_idx, i);
        end

        % Interpolate mode_zvar values for the full length of sections
        zvar = interp1(linspace(min(sections), max(sections), n_points), zvar, sections, 'linear');


        %%
        for p=1:sz(1,2)

            x1 = d;

            x2 = [min(d) max(d)];

            y1 = z(p,:);

            y2 = [zvar(1,p) zvar(1,p)];

            % make common x coordinates

            xq = linspace(min([x1,x2]), max([x1,x2]),1e6);     % New ,x1 and x2 Vector For Interpolation

            y1i = interp1(x1(~isnan(x1)), y1(~isnan(x1)), xq, 'linear','extrap');   % Interpolate To ‘xq’

            y2i = interp1(x2, y2, xq, 'linear','extrap');   % Interpolate To ‘xq’


            % Use polyxpoly to find intersection points for wet perimeter calculus
            [intersection_x, ~] = polyxpoly(x1, y1, x2, y2);


            di=diff(intersection_x);
            width_zvar(1,p)=sum(di(1:2:end));

            area_zvar(1,p)=trapz(xq,y2i)-trapz([min(xq) xq(y1i<zvar(1,p)) max(xq)], [zvar(1,p) y1i(y1i<zvar(1,p)) zvar(1,p)]);
            h_var(1,p)=area_zvar(1,p)./width_zvar(1,p);

        end

        %%
        lim.sections=sections;
        lim.zlim=zlim;
        lim.zvar=zvar;
        lim.hlim=hlim;
        lim.hvar=h_var;
        lim.alim=alim;
        lim.avar=area_zvar;
        lim.wlim=wlim;
        lim.wvar=width_zvar;
        lim.Z=Z;

end

%% plotting

if plt
    switch meth

        case 'all'
            figure
            scatterhist(sections, zlim,'NBins',[sz(1, 2),round(max(zlim)/0.1)],'Location', 'SouthEast','Marker','.','Color','k','Kernel','overlay','Bandwidth',[2;.2]);
            xlabel('Section number');
            ylabel('Peak Elevation (m)');
            xlim([1 sz(1,2)])
            ylim([0 max(zlim)])

            figure

            subplot(4,1,1)

            n_points=100;

            [xi, yi] = meshgrid(linspace(min(sections), max(sections), n_points), linspace(min(zlim), max(zlim), n_points));
            xy = [sections', zlim'];
            [f,~] = ksdensity(xy, [xi(:), yi(:)]);
            f = reshape(f, n_points, n_points);

            contour(xi, yi, f,100, 'LineWidth', 1.5);
            hold on
            scatter(sections,zlim, 'k.')

            xlim([1 sz(1,2)])

            ylabel('Peak elevation (m)')
            xlabel ('Section number')



            subplot(4,1,2)

            scatter(sections,wlim,'k.')

            xlim([1 sz(1,2)])

            ylabel('Width (m)')
            xlabel ('Section number')

            subplot(4,1,3)

            scatter(sections,alim,'k.')

            xlim([1 sz(1,2)])

            ylabel('Area (m^{2})')
            xlabel ('Section number')


            subplot(4,1,4)

            scatter(sections,hlim,'k.')

            xlim([1 sz(1,2)])

            ylabel('Hydraulic depth (m)')
            xlabel ('Section number')



        case 'max'

            figure

            histogram(zlim,round(max(zlim)/0.1))

            ylabel('Count')
            xlabel('Elevation above thalveg (m)')

            yyaxis right

            [f, xi] = ksdensity(zlim, 'Bandwidth', 0.2);

            plot(xi,f,'r-','LineWidth',1)

            xline(Z,'k-','LineWidth',1);
            xline(median(zvar),'r-','LineWidth',1)
            xline(median(zvar)-std(zvar),'r--','LineWidth',1)
            xline(median(zvar)+std(zvar),'r--','LineWidth',1)

            ylabel ('Density probability')

            
            figure

            subplot(4,1,1)

            n_points=100;

            [xi, yi] = meshgrid(linspace(min(sections), max(sections), n_points), linspace(min(zlim), max(zlim), n_points));
            xy = [sections', zlim'];
            [f,~] = ksdensity(xy, [xi(:), yi(:)]);
            f = reshape(f, n_points, n_points);

            contour(xi, yi, f,100, 'LineWidth', 1.5);
            hold on
            scatter(sections,zlim, 'k.')
            plot(sections,zvar, 'red','LineWidth',1)

            xlim([1 sz(1,2)])

            ylabel('Peak elevation (m)')
            xlabel ('Section number')


            subplot(4,1,2)

            scatter(sections,wlim,'k.')
            hold on
            plot(sections,width_zvar, 'red','LineWidth',1)

            xlim([1 sz(1,2)])

            ylabel('Width (m)')
            xlabel ('Section number')

            subplot(4,1,3)

            scatter(sections,alim,'k.')
            hold on
            plot(sections,area_zvar, 'red','LineWidth',1)

            xlim([1 sz(1,2)])

            ylabel('Area (m^{2})')
            xlabel ('Section number')

            subplot(4,1,4)

            scatter(sections,hlim,'k.')
            hold on
            plot(sections,h_var, 'red','LineWidth',1)

            xlim([1 sz(1,2)])

            ylabel('Hydraulic depth (m)')
            xlabel ('Section number')


        case 'lowest'

            figure

            histogram(zlim,round(max(zlim)/0.1))

            ylabel('Count')
            xlabel('Elevation above thalveg (m)')

            yyaxis right

            [f, xi] = ksdensity(zlim, 'Bandwidth', 0.2);

            plot(xi,f,'r-','LineWidth',1)

            xline(Z,'k-','LineWidth',1);
            xline(median(zvar),'r-','LineWidth',1)
            xline(median(zvar)-std(zvar),'r--','LineWidth',1)
            xline(median(zvar)+std(zvar),'r--','LineWidth',1)

            ylabel ('Density probability')

            
            figure

            subplot(4,1,1)

            n_points=100;

            [xi, yi] = meshgrid(linspace(min(sections), max(sections), n_points), linspace(min(zlim), max(zlim), n_points));
            xy = [sections', zlim'];
            [f,~] = ksdensity(xy, [xi(:), yi(:)]);
            f = reshape(f, n_points, n_points);

            contour(xi, yi, f,100, 'LineWidth', 1.5);
            hold on
            scatter(sections,zlim, 'k.')
            plot(sections,zvar, 'red','LineWidth',1)

            xlim([1 sz(1,2)])

            ylabel('Peak elevation (m)')
            xlabel ('Section number')


            subplot(4,1,2)

            scatter(sections,wlim,'k.')
            hold on
            plot(sections,width_zvar, 'red','LineWidth',1)

            xlim([1 sz(1,2)])

            ylabel('Width (m)')
            xlabel ('Section number')

            subplot(4,1,3)

            scatter(sections,alim,'k.')
            hold on
            plot(sections,area_zvar, 'red','LineWidth',1)

            xlim([1 sz(1,2)])

            ylabel('Area (m^{2})')
            xlabel ('Section number')

            subplot(4,1,4)

            scatter(sections,hlim,'k.')
            hold on
            plot(sections,h_var, 'red','LineWidth',1)

            xlim([1 sz(1,2)])

            ylabel('Hydraulic depth (m)')
            xlabel ('Section number')



    end
end

