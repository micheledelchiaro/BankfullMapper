function [final_Zlim,final_Zvar,final_Qlim,final_Qvar,residual_Z,residual_Q,res_Z,res_Q]=visual(HAR,DEM,SW,step,zlim,zvar,qlim,qvar,varargin)

%     map visualization of the bankfull geometry and discharge 
%
% Syntax
%
%    [final_Zlim,final_Zvar,final_Qlim,final_Qvar,residual_Z,residual_Q,res_Z,res_Q]=visual(HAR,DEM,SW,step,zlim,zvar,qlim,qvar)
%
% Description
%
%     This function allows to map the bankfull geometry and discharge based
%     on the results obtained by DETECT_PEAK and MANNINGEQ functions. In
%     detail, both peak elevation values and most probable values for each
%     section are used for the bankfull geometry visualization. Regarding
%     the discharge visualization, both the discharge computed with the
%     peak elevation values and that with the most probable elevation
%     values are used. The bankfull elevation and discharge residuals are
%     also computed by subtracting the most probable values with the peak values.
%
%
%
% Input arguments
%
%     DEM       high resolution GRIDobj
%     HAR       high resolution GRIDobj
%     SW        SWATHobj of the river profile
%     step      stepping in meters of the profiles
%     zlim      peak elevation above thalveg values (see DETECT_PEAK results)
%     zvar      most probable elevation values above thalveg values (see DETECT_PEAK results)
%     qlim      discharge values (computed with peak elevation values zlim using MANNINGEQ function)
%     qvar      discharge values (computed with the most probable elevation values zvar using MANNINGEQ function)
%
% Output arguments
% 
%     final_Zlim        GRIDobj of the REM with the values below the most probable peak elevation 
%     final_Zvar        GRIDobj of the REM with the variable z values 
%     final_Qlim        GRIDobj of the REM with discharge values with constant z 
%     final_Qvar        GRIDobj of the REM with discharge values with variable z
%     residual_Z        GRIDobj of the residual between the residual error between
%                       the peak elevation of each section and the most probable 
%                       elevation 
%     residual_Q        GRIDobj of the residual between the residual error between
%                       the discharge of each section extracted with constant z and the most probable 
%                       discharge 
%     res_Z             array of the elevation residual values    
%     res_Q             array of the discharge values extracted with constant z     
%
%
% 
% Author: Michele Delchiaro (michele.delchiaro[at]uniroma1.it)
% Date: 11. July, 2024
%
% Parse input
p = inputParser;
p.FunctionName = 'visual';
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
addRequired(p,'HAR',@(x) isa(x,'GRIDobj'));
addRequired(p,'SW',@(x) isa(x,'SWATHobj'));
addParameter(p,'step',@(x) isnumeric(x))
addParameter(p,'zlim',@(x) isnumeric(x))
addParameter(p,'zvar',@(x) isnumeric(x))
addParameter(p,'qlim',@(x) isnumeric(x))
addParameter(p,'qvar',@(x) isnumeric(x))
addParameter(p,'plot',true);

parse(p,DEM,HAR,SW,varargin{:});

plt=p.Results.plot;

%% SECTIONING OF STREAM PROFILE
DX = SWATHobj2GRIDobj(SW,HAR,'distx');

dist_along=SW.distx;
dx = SW.dx;
dix = round(step/dx);
ix = 1:dix:length(dist_along);

ix=ix*DEM.cellsize;

DX = reclassify(DX,'definedintervals',ix);

% clip HAR in the SW area
T=SWATHobj2GRIDobj(SW,HAR);
T=T~=0;
REMc=clip(HAR,T);

%% reclassification and plotting of final_Zlim, final_Zvar, residual_Z and res_Z

%
final_Zlim=REMc;

for i =1:length(zlim)
    final_Zlim.Z(DX.Z==i+1)=final_Zlim.Z(DX.Z==i+1)<zlim(i);
end

M=final_Zlim==1;

for i =1:length(zlim)
    final_Zlim.Z(DX.Z==i+1)=zlim(i);
end

final_Zlim=clip(final_Zlim,M);

%
final_Zvar=REMc;

for i =1:length(zvar)
    final_Zvar.Z(DX.Z==i+1)=final_Zvar.Z(DX.Z==i+1)<zvar(i);
end

M=final_Zvar==1;

for i =1:length(zvar)
    final_Zvar.Z(DX.Z==i+1)=zvar(i);
end

final_Zvar=clip(final_Zvar,M);

% 
res_Z=zlim-zvar;

residual_Z=REMc;


for i =1:length(res_Z)
    residual_Z.Z(DX.Z==i+1)=res_Z(i);
    
end

M=~isnan(final_Zvar);
residual_Z=clip(residual_Z,M);

%% reclassification and plotting of final_Qlim, final_Qvar, residual_Q and res_Q

%
final_Qlim=final_Zlim;

for i =1:length(qlim)
    final_Qlim.Z(DX.Z==i+1)=qlim(i);
    
end

M=~isnan(final_Zlim);
final_Qlim=clip(final_Qlim,M);

%
final_Qvar=final_Zvar;

for i =1:length(qvar)
    final_Qvar.Z(DX.Z==i+1)=qvar(i);
    
end

M=~isnan(final_Zvar);
final_Qvar=clip(final_Qvar,M);

% 
res_Q=qlim-qvar;

residual_Q=REMc;


for i =1:length(res_Q)
    residual_Q.Z(DX.Z==i+1)=res_Q(i);
    
end

M=~isnan(final_Qvar);
residual_Q=clip(residual_Q,M);



%% plotting 
if plt

% plot final_Z
figure
imageschs(DEM,final_Zlim,'caxis',[0 max(final_Zlim)],'colorbarylabel','Elevation above thalveg (m)')

title('Final Zlim')

% plot final_Zvar
figure
imageschs(DEM,final_Zvar,'caxis',[0 max(final_Zvar)],'colorbarylabel','Elevation above thalveg (m)')

title('Final Zvar')

% plot final_Qlim
figure
imageschs(DEM,final_Qlim,'caxis',[0 max(final_Qlim)],'colorbarylabel','Discharge (m^{3} s^{-1})')

title('Final Q')

% plot final_Qvar
figure
imageschs(DEM,final_Qvar,'caxis',[0 max(final_Qvar)],'colorbarylabel','Discharge (m^{3} s^{-1})')

title('Final Q with variable z')

% plot residual_Z
figure
imageschs(DEM,residual_Z,'caxis',[-max(residual_Z) max(residual_Z)],'colormap','jet','colorbarylabel','Residual (m)')

title('Residual Z')   

% plot residual_Q
figure
imageschs(DEM,residual_Q,'caxis',[min(residual_Q) max(residual_Q)],'colormap','jet','colorbarylabel','Residual discharge (m^{3} s^{-1})')

title('Residual discharge') 


% plot res_Z and res_Q

figure
subplot(2,1,1)

plot(zlim,'k-')
hold on
plot(zvar,'r-')

ylabel('Peak elevation (m)')
xlabel('Section')

yyaxis right
plot(res_Z,'g')
yline(0)

ylabel('Residual (m)')

legend('Peak elevation','Most probable peak elevation','Residual')
  
subplot(2,1,2)
plot(qlim,'k-')
hold on
plot(qvar,'r-')

ylabel('Discharge (m^{3} s^{-1})')
xlabel('Section')

yyaxis right
plot(res_Q,'g')
yline(0)

ylabel('Residual discharge (m^{3} s^{-1})')

legend('Peak discharge','Most probable peak discharge','Residual')

end

end

