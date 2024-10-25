close all
selectBackgroundCell = 1; %% Set = 1 for images where you want to subtract background and Set = 0 for images without backgroung

% RAG_2G_I1.Location='C:\Users\dwuletj\Documents\JennsFiles\CalciumCode\NikkiCalciumData\PracticeData\Nikon Data\RAG_2_2mM.nd2';
% RAG_20G_I1.Location = 'C:\Users\dwuletj\Documents\JennsFiles\CalciumCode\NikkiCalciumData\PracticeData\Nikon Data\RAG_1_20mM.nd2';
% NOD_20G_I1.Location = 'C:\Users\dwuletj\Documents\JennsFiles\CalciumCode\NikkiCalciumData\PracticeData\Nikon Data\NOD_3_20mM.nd2';
% NOD12wk_2G_I1.Location = 'C:\Users\dwuletj\Documents\JennsFiles\CalciumCode\NikkiCalciumData\PracticeData\LSM 780 Data\12wkNOD_2mM_2.lsm'; 
% NOD12wk_11G_I1.Location = 'C:\Users\dwuletj\Documents\JennsFiles\CalciumCode\NikkiCalciumData\PracticeData\LSM 780 Data\12wkNOD_11mM_2.lsm'; 
% RAG6wk_2mM_1.Location = 'C:\Users\dwuletj\Documents\JennsFiles\CalciumCode\NikkiCalciumData\PracticeData\MorePracticeData\6wkRAG_2mM_1.czi';
% RAG6wk_11mM_1_czi.Location = 'C:\Users\dwuletj\Documents\JennsFiles\CalciumCode\NikkiCalciumData\PracticeData\MorePracticeData\6wkRAG_11mM_1.czi';
% RAG6wk_11mM_1_lsm.Location = 'C:\Users\dwuletj\Documents\JennsFiles\CalciumCode\NikkiCalciumData\PracticeData\MorePracticeData\6wkRAG_11mM_1.lsm';
% NOD12wk_11mM_3.Location = 'C:\Users\dwuletj\Documents\JennsFiles\CalciumCode\NikkiCalciumData\PracticeData\MorePracticeData\12wkNOD_11mM_3.lsm';

fileName.Location = '/Users/jdwulet/Documents/Files from Lacie/Nikki/CalciumAnalysis/JaeAnnUpdate/Sample Data/LSM 780 Data/12wkNOD_2mM_2.lsm'; 
fileName.Location = '/Users/jdwulet/Documents/Files from Lacie/Nikki/CalciumAnalysis/JaeAnnUpdate/Sample Data/LSM 780 Data/12wkNOD_11mM_2.lsm';
% fileName.Location = '/Users/jdwulet/Documents/Files from Lacie/Nikki/CalciumAnalysis/JaeAnnUpdate/Sample Data/Nikon Data/NOD_3_20mM.nd2';
% fileName.Location = '/Users/jdwulet/Documents/Files from Lacie/Nikki/CalciumAnalysis/JaeAnnUpdate/Sample Data/Nikon Data/RAG_1_20mM.nd2';
% fileName.Location = '/Users/jdwulet/Documents/Files from Lacie/Nikki/CalciumAnalysis/JaeAnnUpdate/Sample Data/Nikon Data/RAG_2_2mM.nd2';
output = Run(fileName, selectBackgroundCell);
%save(output)
% RAG_2G_I1=Run(RAG_2G_I1);
% RAG_20G_I1=Run(RAG_20G_I1);
% NOD_20G_I1=Run(NOD_20G_I1);
% NOD12wk_11G_I1=Run(NOD12wk_11G_I1);
% NOD12wk_2G_I1=Run(NOD12wk_2G_I1);
% RAG6wk_2mM_1 = Run(RAG6wk_2mM_1);
% RAG6wk_11mM_1_czi = Run(RAG6wk_11mM_1_czi);
% RAG6wk_11mM_1_lsm = Run(RAG6wk_11mM_1_lsm);
% NOD12wk_11mM_3 = Run(NOD12wk_11mM_3);



%KO_11G_I1.Location='C:\Users\Undergrads\Desktop\Nurin\Videos\3MH_KO_11G_I1.czi';
%KO_11G_I1=Run(KO_11G_I1)


