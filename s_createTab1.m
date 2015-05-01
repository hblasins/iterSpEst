close all;
clear all;
clc;

if ~exist('./Results/sampleResults.mat','file');
    run('../s_analyzeTime.m');
else
    load('./Results/macbethResults.mat');
end


fName = fopen('./Figures/time.tex','w');
fprintf(fName,'\\begin{table}\n');
fprintf(fName,'\\centering\n');
fprintf(fName,'\\begin{tabular}{| l | c | c | c | c | c |}\n');
fprintf(fName,'\\hline\n');
fprintf(fName,'\\multirow{2}{*}{Alg.} & \\multicolumn{2}{c|}{Time, s} & \\multicolumn{2}{c|}{No. iter.} \\\\\n');
fprintf(fName,'& IP & ADMM \\\\\n');
fprintf(fName,'\\hline\n\\hline\n');


fprintf(fName,'Least-squares & %.2f & %.2f  & %i & %i \\\\\n',tlsCvx,tls,1,histBox.pcg.iter(1));
fprintf(fName,'Box & %.2f & %.2f  & %i & %i \\\\\n',tBoxCvx,tBoxAdmm,length(histBox.dualRes),sum(histBox.pcg.iter));
fprintf(fName,'Spatial & %.2f & %.2f  & %i & %i \\\\\n',tSpatialCvx,tSpatialAdmm,length(histSpatial.dualRes),sum(histSpatial.pcg.iter));
fprintf(fName,'Box Spatial & %.2f & %.2f  & %i & %i & %f & %f \\\\\n',tBoxSpatialCvx,tBoxSpatialAdmm,length(histBoxSpatial.dualRes),sum(histBoxSpatial.pcg.iter));


fprintf(fName,'\\hline\\end{tabular}\n\\end{table}\n');
fclose(fName);

