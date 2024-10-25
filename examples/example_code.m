%% LOAD DATA
DEM=GRIDobj("D:\CARTELLA_IRAN\iran geology\DOTTORATO\POTENZA\DEM_San_Severino_Dettaglio.tif");
% DEM=resample(DEM,5);
REM=GRIDobj("D:\CARTELLA_IRAN\iran geology\DOTTORATO\POTENZA\REM_San_Severino_Dettaglio.tif");
REM=resample(REM,DEM,'nearest');
S=shaperead("D:\CARTELLA_IRAN\iran geology\DOTTORATO\POTENZA\San_Severino_Centerline2.shp");

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


%%  VALIDATION

TEST2008=shaperead('canale_attivo_2008.shp');

TEST2008=polygon2GRIDobj(DEM,TEST2008);

MODEL1=~isnan(f_Zvar_low);
MODEL2=~isnan(f_Zvar_max);
MODEL3=~isnan(f_Zlim_low);
MODEL4=~isnan(f_Zlim_low);

conf_matrix1=confusionmat(TEST2008.Z(:)',MODEL1.Z(:)');
conf_matrix2=confusionmat(TEST2008.Z(:)',MODEL2.Z(:)');
conf_matrix3=confusionmat(TEST2008.Z(:)',MODEL3.Z(:)');
conf_matrix4=confusionmat(TEST2008.Z(:)',MODEL4.Z(:)');

figure
subplot(2,2,1)
confusionchart(conf_matrix1,'Normalization','total-normalized')
subplot(2,2,2)
confusionchart(conf_matrix2,'Normalization','total-normalized')
subplot(2,2,3)
confusionchart(conf_matrix3,'Normalization','total-normalized')
subplot(2,2,4)
confusionchart(conf_matrix4,'Normalization','total-normalized')

%
TP1 = conf_matrix1(1,1);
FP1 = conf_matrix1(1,2);
FN1 = conf_matrix1(2,1);
TN1 = conf_matrix1(2,2);

TP2 = conf_matrix2(1,1);
FP2 = conf_matrix2(1,2);
FN2 = conf_matrix2(2,1);
TN2 = conf_matrix2(2,2);

TP3 = conf_matrix3(1,1);
FP3= conf_matrix3(1,2);
FN3 = conf_matrix3(2,1);
TN3 = conf_matrix3(2,2);

TP4 = conf_matrix4(1,1);
FP4 = conf_matrix4(1,2);
FN4 = conf_matrix4(2,1);
TN4 = conf_matrix4(2,2);

%
accuracy1 = (TP1 + TN1) / sum(conf_matrix1(:));
precision1 = TP1 / (TP1 + FP1);
recall1 = TP1 / (TP1 + FN1);
specificity1 = TN1 / (TN1 + FP1);
f1_score1 = 2 * (precision1 * recall1) / (precision1 + recall1);
mcc1 = (TP1 * TN1 - FP1 * FN1) / sqrt((TP1 + FP1) * (TP1 + FN1) * (TN1 + FP1) * (TN1 + FN1));

%
accuracy2 = (TP2 + TN2) / sum(conf_matrix2(:));
precision2 = TP2 / (TP2 + FP2);
recall2 = TP2 / (TP2 + FN2);
specificity2 = TN2 / (TN2 + FP2);
f1_score2 = 2 * (precision2 * recall2) / (precision2 + recall2);
mcc2 = (TP2 * TN2 - FP2 * FN2) / sqrt((TP2 + FP2) * (TP2 + FN2) * (TN2 + FP2) * (TN2 + FN2));

%
accuracy3 = (TP3 + TN3) / sum(conf_matrix3(:));
precision3 = TP3 / (TP3 + FP3);
recall3 = TP3 / (TP3 + FN3);
specificity3 = TN3 / (TN3 + FP3);
f1_score3 = 2 * (precision3 * recall3) / (precision3 + recall3);
mcc3 = (TP3 * TN3 - FP3 * FN3) / sqrt((TP3 + FP3) * (TP3 + FN3) * (TN3 + FP3) * (TN3 + FN3));

%
accuracy4 = (TP4 + TN4) / sum(conf_matrix4(:));
precision4 = TP4 / (TP4 + FP4);
recall4 = TP4 / (TP4 + FN4);
specificity4 = TN4 / (TN4 + FP4);
f1_score4 = 2 * (precision4 * recall4) / (precision4 + recall4);
mcc4 = (TP4 * TN4 - FP4 * FN4) / sqrt((TP4 + FP4) * (TP4 + FN4) * (TN4 + FP4) * (TN4 + FN4));
