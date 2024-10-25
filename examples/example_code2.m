%% LOAD DATA

% DEM=GRIDobj('dem_final_2008.tif');
% DEM.Z(DEM.Z==0)=NaN;
% DEM2=GRIDobj('dem_final_2022.tif');
% DEM2.Z(DEM2.Z==0)=NaN;
% REM=GRIDobj('rem_final_2008.tif');
% REM.Z(REM.Z==0)=NaN;
% REM2=GRIDobj('rem_final_2022.tif');
% REM2.Z(REM2.Z==0)=NaN;
% 
% DEM2=resample(DEM2,DEM,'nearest');
% REM=resample(REM,DEM,'nearest');
% REM2=resample(REM2,DEM,'nearest');
% 
% 
% S=shaperead("marecchia2008.shp");
% 
% S2=shaperead("marecchia2022.shp");
% 
% step=5;

DEM=GRIDobj("D:\CARTELLA_IRAN\iran geology\DOTTORATO\VAL_MARECCHIA\DEM_2008_clipped_roi.tif");
DEM.Z(DEM.Z==0)=NaN;
DEM2=GRIDobj("D:\CARTELLA_IRAN\iran geology\DOTTORATO\VAL_MARECCHIA\DEM_2022_clipped_roi.tif");
DEM2.Z(DEM2.Z==0)=NaN;
REM=GRIDobj("D:\CARTELLA_IRAN\iran geology\DOTTORATO\VAL_MARECCHIA\REM\REM_Marecchia_2008.tif");
REM.Z(REM.Z==0)=NaN;
REM2=GRIDobj("D:\CARTELLA_IRAN\iran geology\DOTTORATO\VAL_MARECCHIA\REM\REM_Marecchia_2022.tif");
REM2.Z(REM2.Z==0)=NaN;

%%
DEM=resample(DEM,1);

DEM2=resample(DEM2,DEM,'nearest');
REM=resample(REM,DEM,'nearest');
REM2=resample(REM2,DEM,'nearest');


S=shaperead("D:\CARTELLA_IRAN\iran geology\DOTTORATO\VAL_MARECCHIA\ANALISI\test2\marecchia2008.shp");

S2=shaperead("D:\CARTELLA_IRAN\iran geology\DOTTORATO\VAL_MARECCHIA\ANALISI\test2\marecchia2022.shp");

step=50;

% associare i valori nan del dem1 ai valori reali del dem2
DEM.Z(isnan(DEM.Z)& ~isnan(DEM2.Z))=DEM2.Z(isnan(DEM.Z)& ~isnan(DEM2.Z));
REM.Z(isnan(REM.Z)& ~isnan(REM2.Z))=REM2.Z(isnan(REM.Z)& ~isnan(REM2.Z));

%% PROFILING
[d,z,z_dem,x,y,SW]=prof(DEM,REM,S,'step',step,'width',2500,'smooth',500,'plot',false);
[d2,z2,z_dem2,x2,y2,SW2]=prof(DEM2,REM2,S,'step',step,'width',2500,'smooth',500,'plot',false); %stessa traccia del time 1

%% HYDRAULIC DEPTH FUNCTION
bank=bankfull(d,z,'max_depth',6);
bank2=bankfull(d2,z2,'max_depth',6);

%% PEAK EXTRACTION
lim_all=detect_peak(bank.h,bank.area,bank.width,bank.elevation,'peak','all');
lim_low=detect_peak(bank.h,bank.area,bank.width,bank.elevation,d,z,'peak','lowest');
lim_max=detect_peak(bank.h,bank.area,bank.width,bank.elevation,d,z,'peak','max');

lim_all2=detect_peak(bank2.h,bank2.area,bank2.width,bank2.elevation,'peak','all');
lim_low2=detect_peak(bank2.h,bank2.area,bank2.width,bank2.elevation,d2,z2,'peak','lowest');
lim_max2=detect_peak(bank2.h,bank2.area,bank2.width,bank2.elevation,d2,z2,'peak','max');


%% SLOPE COMPUTATION
[slope1,zw1,~,~]=sectiongradient(DEM,SW,step,'windowsize',10);
[slope2,zw2,~,~]=sectiongradient(DEM2,SW2,step,'windowsize',10);

%% MANNING'S EQUATION 

R_lowA=manningseq(SW,z,z_dem,step,d,lim_low.zlim,'n',0.05,'slope',slope1,'plot',false);
R_maxA=manningseq(SW,z,z_dem,step,d,lim_max.zlim,'n',0.05,'slope',slope1,'plot',false);

R_lowB=manningseq(SW,z,z_dem,step,d,lim_low.zvar,'n',0.05,'slope',slope1,'plot',false);
R_maxB=manningseq(SW,z,z_dem,step,d,lim_max.zvar,'n',0.05,'slope',slope1,'plot',false);

R_low2A=manningseq(SW2,z2,z_dem2,step,d2,lim_low2.zlim,'n',0.05,'slope',slope2,'plot',false);
R_max2A=manningseq(SW2,z2,z_dem2,step,d2,lim_max2.zlim,'n',0.05,'slope',slope2,'plot',false);

R_low2B=manningseq(SW2,z2,z_dem2,step,d2,lim_low2.zvar,'n',0.05,'slope',slope2,'plot',false);
R_max2B=manningseq(SW2,z2,z_dem2,step,d2,lim_max2.zvar,'n',0.05,'slope',slope2,'plot',false);


%% VISUALIZATION

[f_Zlim_low,f_Zvar_low,f_Qlim_low,f_Qvar_low,residual_Z_low,residual_Q_low,res_Z_low,res_Q_low]=visual(REM,DEM,SW,step,lim_low.zlim,lim_low.zvar,R_lowA.Q_zvar,R_lowB.Q_zvar,'plot',false);

[f_Zlim_max,f_Zvar_max,f_Qlim_max,f_Qvar_max,residual_Z_max,residual_Q_max,res_Z_max,res_Q_max]=visual(REM,DEM,SW,step,lim_max.zlim,lim_max.zvar,R_maxA.Q_zvar,R_maxB.Q_zvar,'plot',false);

[f_Zlim_low2,f_Zvar_low2,f_Qlim_low2,f_Qvar_low2,residual_Z_low2,residual_Q_low2,res_Z_low2,res_Q_low2]=visual(REM2,DEM2,SW2,step,lim_low2.zlim,lim_low2.zvar,R_low2A.Q_zvar,R_low2B.Q_zvar,'plot',false);

[f_Zlim_max2,f_Zvar_max2,f_Qlim_max2,f_Qvar_max2,residual_Z_max2,residual_Q_max2,res_Z_max2,res_Q_max2]=visual(REM2,DEM2,SW2,step,lim_max2.zlim,lim_max2.zvar,R_max2A.Q_zvar,R_max2B.Q_zvar,'plot',false);


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

%
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


