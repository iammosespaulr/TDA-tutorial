% This script demonstrates the use of plex on a simple "house" example

clc; clear; close all;
import edu.stanford.math.plex4.*;

max_dimension = 2;
max_filtration_value = 4;
num_divisions = 1000;

% create the set of points
point_cloud = examples.PointCloudExamples.getHouseExample();

% create a Vietoris-Rips stream 
stream = api.Plex4.createVietorisRipsStream(point_cloud, max_dimension + 1, max_filtration_value, num_divisions);

% get the default persistence algorithm
persistence = api.Plex4.getDefaultSimplicialAlgorithm(max_dimension);

% this initially creates a set of intervals which contains the filtration
% indices (which are integers).
filtration_value_intervals = persistence.computeIntervals(stream)

% create the barcode plots
plot_barcodes(filtration_value_intervals)