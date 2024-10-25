%statistical test

T = readtable('/Users/jdwulet/Google Drive/LabComputer/Nurins Exper/Statistical_tests_for_manuscript.xlsx');
S = table2struct( T ); %Create a structure so that you can call the columns by fieldnames
F = fieldnames(S); 