%% Simulation environment for iGEM Chalmers 2018
%  Author: Francisco Zorrilla
%  2018-08-05
%
%  Community model of Human Gut + Colorectal Cancer + Engineering
%  S.boulardii.

cd('C:\Users\zorrilla\Desktop\iGEM simulation environment')
%% 1. Load Genome Scale Models

% Saccharomyces boulardii derived from homology using S.cerevisiae as template
load('modelSbo.mat'); 

% Colon tissue GEM derived from HMR2.0
modelColon =importModel('colon.xml');% Colon specific GEM obtained from http://www.metabolicatlas.org/downloads/tissue

