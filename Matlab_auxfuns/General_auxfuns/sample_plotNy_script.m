close all;clear all;clc
t=linspace(0,2*pi,100);
y1=2*sin(t);
y2=2*sin(2*t);

y3=200*sin(t+.5);
y4=100*sin(t+.25);

y5=2000*sin(t+.15);
y6=3000*sin(t+.75);


acc=cos(t)+1;
vel=sin(t)+t;
pos=-cos(t)+t.^2./2+1;


%You can either not include an input arguement, or leave its answer blank
%to leave default setting. If you give an output arguement, it will give
%you all of the handles used to construct the plot. See 'help plotNy' for a
%full description of each property

%% Sets of three inputs (no figure given)

linecolors={'r' [0 .5 0] 'b'};

h3i=plotNy(t,acc,1,...
    t,vel,2,...
    t,pos,3,...
    'Linewidth',1,...
    'YAxisLabels',{'[m/s^2]' '[m/s]' '[m]'},...
    'XAxisLabel','Time [s]',...
    'TitleStr','Partical Movement with Time',...
    'LineColor',linecolors,...
    'FontSize',12,...
    'Fontname','TimesNewRoman',...
    'Grid','on',...
    'LegendString',{'Acceleration [m/s^2]' 'Velocity [m/s]' 'Position [m]'});


%Change the axis colors to match the requested line colors
for i=1:length(h3i.ax)
	set(h3i.ax(i),'ycolor',linecolors{i});	
end



%% Cell array inputs

%here no color or line/marker style is specified, so they matlab default
%color order is used

fh=figure('units','normalized','position',[0.25 0.25 .5 .5]);
hci=plotNy({t,t,t,t,t,t},{y1,y2,y3,y4,y5,y6},[1 1 2 2 4 4 ],...
    'XAxisLabel','Time [s]',...
    'TitleStr','Sin Waves - Triplet Input',...
    'Grid','on',...
    'Parent',fh,... 
    'Xlim',[1 2*pi],...
    'Ylim',[-5 5; -300 300; -4000 4000]);


%% Cell array y with shared x

%Here I have not specified line color and style for each, but have
%overridden the default order. plotNy will plot with unique combinations of
%the color order and line style order

fh2=figure;
plot(1:10);

hsx=plotNy(t,{y1,y2,y3,y4,y5,y6},[1 1 2 2 3 3],...
    'YAxisLabels',{'A1 [m]' 'A2 [cm]' 'A3 [mm]'},...
    'XAxisLabel','Time [s]',...
    'TitleStr','Sin Waves - Triplet Input',...
    'Grid','on',...
    'Parent',fh2,... % I have give same parent figure, so function clears old content
    'Xlim',[1 2*pi],...
    'Ylim',[-5 5; -300 300; -4000 4000],...
    'colorord',[1 0 0; 0 1 0; 0 0 1; .5 .5 .5],...
    'lineord',{'-' '-.'},...
    'LegendString',{'y1 [m]' 'y2 [m]' 'y3 [cm]' 'y4 [cm]' 'y5 [mm]' 'y6 [mm]'},...
	'LegendLoc','none');


