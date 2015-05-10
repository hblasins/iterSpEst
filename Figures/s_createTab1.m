%% Run this script to generate Tab. 1 from 
%
% H. Blasinski, J. Farrell and B. Wandell; 'An iterative algorithm for spectral
% estimation with spatial smoothing,' ICIP 2015, Quebec City
%
% Copyright, Henryk Blasinski, 2015


close all;
clear all;
clc;

fName = fullfile(iterSpEstRoot,'Results','macbethResults.mat');
if ~exist(fName,'file')
    sName = fullfile(iterSpEstRoot,'s_analyzeTime');
    run(sName);
else
    load(fName);
end


fName = fullfile(iterSpEstRoot,'Figures','time.tex');
fid = fopen(fName,'w');
fprintf(fid,'\\begin{table}\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\begin{tabular}{| l | c | c | c | c | c |}\n');
fprintf(fid,'\\hline\n');
fprintf(fid,'\\multirow{2}{*}{Algorithm} & \\multicolumn{2}{c|}{Time, s} & \\multicolumn{2}{c|}{No. iter.} \\\\\n');
fprintf(fid,'& \\texttt{cvx} & ADMM & ADMM & CG \\\\\n');
fprintf(fid,'\\hline\n\\hline\n');


fprintf(fid,'Least-squares & %.2f & %.2f  & %i & %i \\\\\n',tlsCvx,tls,1,histBox.pcg.iter(1));
fprintf(fid,'Box & %.2f & %.2f  & %i & %i \\\\\n',tBoxCvx,tBoxAdmm,length(histBox.dualRes),sum(histBox.pcg.iter));
fprintf(fid,'Spatial & %.2f & %.2f  & %i & %i \\\\\n',tSpatialCvx,tSpatialAdmm,length(histSpatial.dualRes),sum(histSpatial.pcg.iter));
fprintf(fid,'Box Spatial & %.2f & %.2f  & %i & %i & %f & %f \\\\\n',tBoxSpatialCvx,tBoxSpatialAdmm,length(histBoxSpatial.dualRes),sum(histBoxSpatial.pcg.iter));


fprintf(fid,'\\hline\\end{tabular}\n\\end{table}\n');
fclose(fid);

