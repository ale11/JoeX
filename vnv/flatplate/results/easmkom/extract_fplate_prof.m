clear all; clc; close all;

addpath /Users/alejandrocampos/myjoe/trunk/src/tools;

geom = 'plate';
model = 'sst';
wallfile = 'wall-7.dat';

prof_and_wall_vals(geom,model,wallfile);

rmpath /Users/alejandrocampos/myjoe/trunk/src/tools;

