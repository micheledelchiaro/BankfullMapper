%% EXAMPLE2: MARECCHIA RIVER CASE STUDY (EMILIA ROMAGNA, ITALY)

%% WORKING DIRECTORY
WD='path\to\downloaded\data\';

%% INPUT DATA
% load DEMs 2009 and 2022
DEM=GRIDobj("DEM_Marecchia_2008.tif");
DEM2=GRIDobj("DEM_Marecchia_2022.tif");

% load REMs 2009 and 2022
REM=GRIDobj("REM_Marecchia_2008.tif");
REM2=GRIDobj("REM_Marecchia_2022.tif");

DEM2=resample(DEM2,DEM,'nearest');
REM=resample(REM,DEM,'nearest');
REM2=resample(REM2,DEM,'nearest');

% load shapefiles of river courses related to 2009 and 2022
S=shaperead("marecchia2008.shp");
S2=shaperead("marecchia2022.shp");

% Optional inputs
step=50;        % stepping in meters between transversal river profiles
width=2500;     % width in meters of the transversal river profiles
smooth=500;     % smoothing of the planar trace from which the profiles are extracted
max_depth=6;    % maximum height from the thalveg for which the bankfull is computed
n=0.05;         % n riverbed roughness coefficient in Manning's equation (Manning, 1904)

%% PROFILING
% 2009
[d,z,z_dem,x,y,SW]=prof(DEM,REM,S,'step',step,'width',width,'smooth',smooth,'plot',true);
% 2022
[d2,z2,z_dem2,x2,y2,SW2]=prof(DEM2,REM2,S,'step',step,'width',width,'smooth',smooth,'plot',true);

%% HYDRAULIC DEPTH FUNCTION COMPUTATION
% 2009
bank=bankfull(d,z,'max_depth',max_depth);
% 2022
bank2=bankfull(d2,z2,'max_depth',max_depth);

%% PEAKS EXTRACTION FROM HYDRAULIC DEPTH FUNCTION 
% 2009
% ALL MODE: all the peaks are extracted
lim_all=detect_peak(bank.h,bank.area,bank.width,bank.elevation,'peak','all');
% LOWEST MODE: the lowest peaks are extracted
lim_low=detect_peak(bank.h,bank.area,bank.width,bank.elevation,d,z,'peak','lowest');
% MAX MODE: the most prominent peaks are extracted
lim_max=detect_peak(bank.h,bank.area,bank.width,bank.elevation,d,z,'peak','max');

% 2022
% ALL MODE: all the peaks are extracted
lim_all2=detect_peak(bank2.h,bank2.area,bank2.width,bank2.elevation,'peak','all');
% LOWEST MODE: the lowest peaks are extracted
lim_low2=detect_peak(bank2.h,bank2.area,bank2.width,bank2.elevation,d2,z2,'peak','lowest');
% MAX MODE: the most prominent peaks are extracted
lim_max2=detect_peak(bank2.h,bank2.area,bank2.width,bank2.elevation,d2,z2,'peak','max');


%% SLOPE COMPUTATION
ws=10;          % cell windowsize for moving average computation on slope smoothing
% 2009
[slope1,zw1,~,~]=sectiongradient(DEM,SW,step,'windowsize',ws);
% 2022
[slope2,zw2,~,~]=sectiongradient(DEM2,SW2,step,'windowsize',ws);

%% MANNING'S EQUATION COMPUTATION (for lowest and max mode extracted peaks)
% 2009
% discharge computation using the measured peak elevation 
R_lowA=manningseq(SW,z,z_dem,step,d,lim_low.zlim,'n',n,'slope',slope1,'plot',false);
R_maxA=manningseq(SW,z,z_dem,step,d,lim_max.zlim,'n',n,'slope',slope1,'plot',false);
% discharge computation using the most probable elevation
R_lowB=manningseq(SW,z,z_dem,step,d,lim_low.zvar,'n',n,'slope',slope1,'plot',false);
R_maxB=manningseq(SW,z,z_dem,step,d,lim_max.zvar,'n',n,'slope',slope1,'plot',false);

% 2022
% discharge computation using the measured peak elevation 
R_low2A=manningseq(SW2,z2,z_dem2,step,d2,lim_low2.zlim,'n',n,'slope',slope2,'plot',false);
R_max2A=manningseq(SW2,z2,z_dem2,step,d2,lim_max2.zlim,'n',n,'slope',slope2,'plot',false);
% discharge computation using the most probable elevation
R_low2B=manningseq(SW2,z2,z_dem2,step,d2,lim_low2.zvar,'n',n,'slope',slope2,'plot',false);
R_max2B=manningseq(SW2,z2,z_dem2,step,d2,lim_max2.zvar,'n',n,'slope',slope2,'plot',false);


%% BANKFULL GEOMETRY AND DISCHARGE MAPPING 
% 2009
% for lowest mode extracted peaks
[f_Zlim_low,f_Zvar_low,f_Qlim_low,f_Qvar_low,residual_Z_low,residual_Q_low,res_Z_low,res_Q_low]=visual(REM,DEM,SW,step,lim_low.zlim,lim_low.zvar,R_lowA.Q_zvar,R_lowB.Q_zvar,'plot',false);
% for max mode extracted peaks
[f_Zlim_max,f_Zvar_max,f_Qlim_max,f_Qvar_max,residual_Z_max,residual_Q_max,res_Z_max,res_Q_max]=visual(REM,DEM,SW,step,lim_max.zlim,lim_max.zvar,R_maxA.Q_zvar,R_maxB.Q_zvar,'plot',false);

% 2022
% for lowest mode extracted peaks
[f_Zlim_low2,f_Zvar_low2,f_Qlim_low2,f_Qvar_low2,residual_Z_low2,residual_Q_low2,res_Z_low2,res_Q_low2]=visual(REM2,DEM2,SW2,step,lim_low2.zlim,lim_low2.zvar,R_low2A.Q_zvar,R_low2B.Q_zvar,'plot',false);
% for max mode extracted peaks
[f_Zlim_max2,f_Zvar_max2,f_Qlim_max2,f_Qvar_max2,residual_Z_max2,residual_Q_max2,res_Z_max2,res_Q_max2]=visual(REM2,DEM2,SW2,step,lim_max2.zlim,lim_max2.zvar,R_max2A.Q_zvar,R_max2B.Q_zvar,'plot',false);


%% EXPORTING OUTPUTS 
% 2009
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

% 2022
GRIDobj2geotiff(f_Zlim_low2,'final_Zlim_low2.tif')
GRIDobj2geotiff(f_Zlim_max2,'final_Zlim_max2.tif')

GRIDobj2geotiff(f_Zvar_low2,'final_Zvar_low2.tif')
GRIDobj2geotiff(f_Zvar_max2,'final_Zvar_max2.tif')

GRIDobj2geotiff(f_Qlim_low2,'final_Qlim_low2.tif')
GRIDobj2geotiff(f_Qlim_max2,'final_Qlim_max2.tif')

GRIDobj2geotiff(f_Qvar_low2,'final_Qvar_low2.tif')
GRIDobj2geotiff(f_Qvar_max2,'final_Qvar_max2.tif')

GRIDobj2geotiff(residual_Z_max2,'residual_Z_max2.tif')
GRIDobj2geotiff(residual_Q_max2,'residual_Q_max2.tif')

GRIDobj2geotiff(residual_Z_low2,'residual_Z_low2.tif')
GRIDobj2geotiff(residual_Q_low2,'residual_Q_low2.tif')


