%% EXAMPLE1: POTENZA RIVER CASE STUDY (MARCHE, ITALY)

%% WORKING DIRECTORY
WD='path\to\downloaded\data\';

%% LOAD DATA
DEM=GRIDobj(WD+"DEM_San_Severino_Dettaglio.tif");

REM=GRIDobj(WD+"REM_San_Severino_Dettaglio.tif");
REM=resample(REM,DEM,'nearest');
S=shaperead(WD+"San_Severino_Centerline2.shp");

step=50;

%% PROFILING
[d,z,z_dem,x,y,SW]=prof(DEM,REM,S,'step',step,'width',200,'smooth',0,'plot',true);

%% HYDRAULIC DEPTH FUNCTION
bank=bankfull(d,z,'max_depth',6);

%% PEAK EXTRACTION
lim_all=detect_peak(bank.h,bank.area,bank.width,bank.elevation,'peak','all');
lim_low=detect_peak(bank.h,bank.area,bank.width,bank.elevation,d,z,'peak','lowest');
lim_max=detect_peak(bank.h,bank.area,bank.width,bank.elevation,d,z,'peak','max');

%% SLOPE COMPUTATION
[slope,~,~,~]=sectiongradient(DEM,SW,step,'windowsize',10);

%% MANNING'S EQUATION 

R_low=manningseq(SW,z,z_dem,step,d,lim_low.zlim,'n',0.05,'slope',slope);
R_max=manningseq(SW,z,z_dem,step,d,lim_max.zlim,'n',0.05,'slope',slope);

R_low2=manningseq(SW,z,z_dem,step,d,lim_low.zvar,'n',0.05,'slope',slope);
R_max2=manningseq(SW,z,z_dem,step,d,lim_max.zvar,'n',0.05,'slope',slope);


%% VISUALIZATION

[f_Zlim_low,f_Zvar_low,f_Qlim_low,f_Qvar_low,residual_Z_low,residual_Q_low,res_Z_low,res_Q_low]=visual(REM,DEM,SW,step,lim_low.zlim,lim_low.zvar,R_low.Q_zvar,R_low2.Q_zvar,'plot',false);

[f_Zlim_max,f_Zvar_max,f_Qlim_max,f_Qvar_max,residual_Z_max,residual_Q_max,res_Z_max,res_Q_max]=visual(REM,DEM,SW,step,lim_max.zlim,lim_max.zvar,R_max.Q_zvar,R_max2.Q_zvar,'plot',false);


%% EXPORT2geotiff

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
