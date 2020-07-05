close all
clear
clc

load_javaplex;

%diary('trial1')

global delta alpha gamma beta omega

delta=0.3;
alpha=1;
beta=1;
omega=1.2;

% Use one of the lines given below, comment the other one out.
%GGAMM=[0.36]; % Single gamma value
GGAMM=[0.35:0.01:0.38]; % Array of Gamma values 

gm_size=size(GGAMM);

for i_g=1:gm_size(2)

j=33;

gamma=GGAMM(i_g);

[t x]=ode45(@duffing,0:2*pi/omega/1000:20000,[0 1]);
x_size=size(x);

low=10000;
high=20000;

%figure,
%subplot(1,2,1);

%Plots Phase space trajectory
%plot(x(low:high,2),x(low:high,1),'b') 
%axis tight
%axis equal
%title(strcat(num2str(GGAMM(i_g)),'phase space'));

%Creates an array that skips j-number of points in the solution set
a(:,1)= x(low:j:high,1);
a(:,2)= x(low:j:high,2);


total = high-low;  
count = floor((total)/j);   %Number of points
% disp(count);
% disp(total);

%PLots the points that have been used for sampling in Javaplex
%subplot(1,2,2)
%scatter(a(:,2),a(:,1),'.');
%axis tight
%axis equal
%title(strcat((num2str(count)), ' of ', (num2str(total))," ",'points'));

%section end


%Javaplex work Begins
import edu.stanford.math.plex4.*;

max_dimension = 2;
max_filtration_value = 1.5;
num_divisions = 2000;

% create the set of points

point_cloud = a;

% create a Vietoris-Rips stream 
stream = api.Plex4.createVietorisRipsStream(point_cloud, max_dimension, max_filtration_value, num_divisions);

% get persistence algorithm over Z/2Z
persistence = api.Plex4.getModularSimplicialAlgorithm(max_dimension, 2);

% compute the intervals
intervals = persistence.computeIntervals(stream);
disp(GGAMM(i_g));
disp(high);
disp(low);
disp(j);
disp(intervals);


val=gamma*10;

%Name for the title of the plot
name=strcat(num2str(GGAMM(i_g)*100),' Barcode ', num2str(count),' ','Points');

% create the barcode plots
options.filename = name;
options.max_filtration_value = max_filtration_value;
options.max_dimension = max_dimension - 1;

% options.fileformat = ;
plot_barcodes(intervals, options);

clear a


end

