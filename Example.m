%%%
I1.Location='C:\Users\Undergrads\Dropbox\Benninger Lab Shared Files\NikkiCalciumCode\6wkRAG_11mM_1.czi'
I1=OpenFilesV2(I1);
% figure
% dbdb_untreated_01_2=OpenFilesV2(dbdb_untreated_01);
% figure
% dbdb_untreated_01_3=OpenFilesV2(dbdb_untreated_01);

% dbdb_untreated_02.Location='C:\Users\Undergrads\Desktop\Nurin\JoshCaData\101216\BLK\dbdb_untreated_11mM_02.czi'
% figure
% dbdb_untreated_02_1=OpenFilesV2(dbdb_untreated_02);
% figure
% dbdb_untreated_02_2=OpenFilesV2(dbdb_untreated_02);
% % figure
% % dbdb_untreated_02_3=OpenFilesV2(dbdb_untreated_02);
% 
% dbdb_untreated_03.Location='C:\Users\Undergrads\Desktop\Nurin\JoshCaData\101216\BLK\dbdb_untreated_11mM_03.czi'
% figure
% dbdb_untreated_03_1=OpenFilesV2(dbdb_untreated_03);
% figure
% dbdb_untreated_03_2=OpenFilesV2(dbdb_untreated_03);
% figure
% dbdb_untreated_03_3=OpenFilesV2(dbdb_untreated_03);

%%


dbdb_untreatedData(1,1)=dbdb_untreated_01_1.RAA; dbdb_untreatedData(1,2)=dbdb_untreated_01_1.RMCA;  dbdb_untreatedData(1,3:5)=dbdb_untreated_01_1.percentCorr;
dbdb_untreatedData(2,1)=dbdb_untreated_01_2.RAA; dbdb_untreatedData(2,2)=dbdb_untreated_01_2.RMCA;  dbdb_untreatedData(2,3:5)=dbdb_untreated_01_2.percentCorr;
% dbdb_untreatedData(3,1)=dbdb_untreated_01_3.RAA; dbdb_untreatedData(3,2)=dbdb_untreated_01_3.RMCA;  dbdb_untreatedData(3,3:5)=dbdb_untreated_01_3.percentCorr;

dbdb_untreatedData(4,1)=dbdb_untreated_02_1.RAA; dbdb_untreatedData(4,2)=dbdb_untreated_02_1.RMCA;  dbdb_untreatedData(4,3:5)=dbdb_untreated_02_1.percentCorr
dbdb_untreatedData(5,1)=dbdb_untreated_02_2.RAA; dbdb_untreatedData(5,2)=dbdb_untreated_02_2.RMCA;  dbdb_untreatedData(5,3:5)=dbdb_untreated_02_2.percentCorr
% dbdb_untreatedData(6,1)=dbdb_untreated_02_3.RAA; dbdb_untreatedData(6,2)=dbdb_untreated_02_3.RMCA;  dbdb_untreatedData(6,3:5)=dbdb_untreated_02_3.percentCorr

dbdb_untreatedData(7,1)=dbdb_untreated_03_1.RAA; dbdb_untreatedData(7,2)=dbdb_untreated_03_1.RMCA;  dbdb_untreatedData(7,3:5)=dbdb_untreated_03_1.percentCorr
dbdb_untreatedData(8,1)=dbdb_untreated_03_2.RAA; dbdb_untreatedData(8,2)=dbdb_untreated_03_2.RMCA;  dbdb_untreatedData(8,3:5)=dbdb_untreated_03_2.percentCorr
dbdb_untreatedData(9,1)=dbdb_untreated_03_3.RAA; dbdb_untreatedData(9,2)=dbdb_untreated_03_3.RMCA;  dbdb_untreatedData(9,3:5)=dbdb_untreated_03_3.percentCorr
%%
filename='H:\JoshCaAnalysis2\101216\dbdb_untreated.mat';
save(filename)

%% Save

ppt=saveppt2('dbdb_untreated','init'); 
for i=1:21
    figure(i)
    saveppt2('ppt',ppt) 
end

%%





