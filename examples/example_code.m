%% EXAMPLE1: POTENZA RIVER CASE STUDY (MARCHE, ITALY)

%% WORKING DIRECTORY
WD='path\to\downloaded\data\';

%% INPUT DATA
% load DEM
DEM=GRIDobj(WD+"DEM_Potenza.tif");

% load REM
REM=GRIDobj(WD+"REM_Potenza.tif");
REM=resample(REM,DEM,'nearest');

% load shapefile of the river course
S=shaperead(WD+"potenza.shp");

% Optional inputs
step=50;        % stepping in meters between transversal river profiles
width=200;     % width in meters of the transversal river profiles
smooth=0;     % smoothing of the planar trace from which the profiles are extracted
max_depth=6;    % maximum height from the thalveg for which the bankfull is computed
n=0.05;         % n riverbed roughness coefficient in Manning's equation (Manning, 1904)

%% PROFILING
[d,z,z_dem,x,y,SW]=prof(DEM,REM,S,'step',step,'width',width,'smooth',smooth,'plot',true);

%% HYDRAULIC DEPTH FUNCTION COMPUTATION
bank=bankfull(d,z,'max_depth',max_depth);

%% PEAKS EXTRACTION FROM HYDRAULIC DEPTH FUNCTION 
% ALL MODE: all the peaks are extracted
lim_all=detect_peak(bank.h,bank.area,bank.width,bank.elevation,'peak','all');
% LOWEST MODE: the lowest peaks are extracted
lim_low=detect_peak(bank.h,bank.area,bank.width,bank.elevation,d,z,'peak','lowest');
% MAX MODE: the most prominent peaks are extracted
lim_max=detect_peak(bank.h,bank.area,bank.width,bank.elevation,d,z,'peak','max');

%% SLOPE COMPUTATION
ws=10;          % cell windowsize for moving average computation on slope smoothing
[slope,~,~,~]=sectiongradient(DEM,SW,step,'windowsize',ws);

%% MANNING'S EQUATION COMPUTATION (for lowest and max mode extracted peaks)
% discharge computation using the measured peak elevation 
R_low=manningseq(SW,z,z_dem,step,d,lim_low.zlim,'n',n,'slope',slope);
R_max=manningseq(SW,z,z_dem,step,d,lim_max.zlim,'n',n,'slope',slope);
% discharge computation using the most probable elevation
R_low2=manningseq(SW,z,z_dem,step,d,lim_low.zvar,'n',n,'slope',slope);
R_max2=manningseq(SW,z,z_dem,step,d,lim_max.zvar,'n',n,'slope',slope);


%% BANKFULL GEOMETRY AND DISCHARGE MAPPING 
% for lowest mode extracted peaks
[f_Zlim_low,f_Zvar_low,f_Qlim_low,f_Qvar_low,residual_Z_low,residual_Q_low,res_Z_low,res_Q_low]=visual(REM,DEM,SW,step,lim_low.zlim,lim_low.zvar,R_low.Q_zvar,R_low2.Q_zvar,'plot',false);
% for max mode extracted peaks
[f_Zlim_max,f_Zvar_max,f_Qlim_max,f_Qvar_max,residual_Z_max,residual_Q_max,res_Z_max,res_Q_max]=visual(REM,DEM,SW,step,lim_max.zlim,lim_max.zvar,R_max.Q_zvar,R_max2.Q_zvar,'plot',false);


%% EXPORTING OUTPUTS 

GRIDobj2geotiff(f_Zlim_low,'final_Zlim_low.tif')
GRIDobj2geotiff(f_Zlim_max,'final_Zlim_max.tif')

GRIDobj2geotiff(f_Zvar_low,'final_Zvar_low.tif')
GRIDobj2geotiff(f_Zvar_max,'final_Zvar_max.tif')

GRIDobj2geotiff(f_Qlim_low,'final_Qlim_low.tif')
GRIDobj2geotiff(f_Qlim_max,'final_Qlim_max.tif')

GRIDobj2geotiff(f_Qvar_low,'final_Qvar_low.tif')
GRIDobj2geotiff(f_Qvar_max,'final_Qvar_max.tif')

GRIDobj2geotiff(residual_Z_max,'residual_Z_max.tif')
GRIDobj2geotiff(residual_Q_max,'residual_Q_max.tif')

GRIDobj2geotiff(residual_Z_low,'residual_Z_low.tif')
GRIDobj2geotiff(residual_Q_low,'residual_Q_low.tif')
