%% read in the data files

% example of reading in data
fileID = fopen('Ba133_1325s.usf');
spec = textscan(fileID,'%u %*[^\n]');
fclose(fileID);
thisSpectrum = spec{1};

fileIDAm = fopen('Am241_662s.usf');
specAm = textscan(fileIDAm,'%u %*[^\n]');
fclose(fileIDAm);
thisSpectrumAm = specAm{1};

fileIDnofilt = fopen('Spectrum-NoFilter.usf');
specnofilt = textscan(fileIDnofilt,'%u %*[^\n]');
fclose(fileIDnofilt);
thisSpectrumnofilt = specnofilt{1};

fileIDzn = fopen('zinc filter.usf');
speczn = textscan(fileIDzn,'%u %*[^\n]');
fclose(fileIDzn);
thisSpectrumzn = speczn{1};

fileIDSt = fopen('Steel_filter.usf');
specSt = textscan(fileIDSt,'%u %*[^\n]');
fclose(fileIDSt);
thisSpectrumSt = specSt{1};
figure 
plot(x_keV, thisSpectrumSt)

%% Tube settings (convert the measured values to kV and mA)

%% Calibrate ADU to Energy
% start by identifying peaks in your data and matching these to the known
% emissions

adu = [312, 542, 1023, 1389];
kev = [17.7, 30.9, 59.6, 80.9];

[p, s] = polyfit(adu, kev, 1);


x_adu = linspace(1, 4096, 4096)';
x_kev = p(1)*x_adu + p(2);
figure
plot(x_kev, thisSpectrumAm); %Am graph 

xlim([50 65]);


%% Energy resolution @ 60 keV
% You need to fit a particular type of function to a peak that is well
% separated from the others

%This will not run because curve-fitting toolbox (cftool) gui used to
%calculate a fitted model for the peak at 60keV, within the x range of
%X_peak.
%3 term Gaussian
%fittedmodel, goodness, output

y=feval(fittedmodel,X_peak);
plot(X_peak,y);
disp(max(y));
[m, I]=max(y);
disp(I)

%FWHM= Width of Curve at half Y-value of maximum
HM=0.5*max(y)


%Method: create linspace of values between X-peak min and max (i.e. along the entire curve of peak), and plot
%against fittedmodel (many different linspace values attempted before finding accuracy to
%within 4 decimal places)
%First check the range of X_values of the peak: 
diff= X_peak(401)-X_peak(1);
smalldiff=diff/100;
%Then create a linspace of 10000000 x-value points on the curve:
z=linspace(X_peak(1), X_peak(401),10000000);
%Evaluate the fitted gaussian curve at each of these x-values
y2=feval(fittedmodel,z);
%Determine the Half-maximum value:
HM2=0.5*max(y2)
%Find the points on the curve that are equal to the Half-Maximum:
idx= find(y2 ==HM2)
j2 = find (y2>= HM2);
disp(y2(j2(1)));
disp(y2(j2(end)));
X1=j2(1);
X2= j2(end);

%Now that we have found the index of the corresponding half maximum
%Y-value, we can find the energy resolution:
FWHM= ((X2-X1)/10000000)*diff
Resolution = (FWHM/I)*100



%% Read in your measured specta

%% Compton shift energy axis
% Based on the calibration you know the conversion from ADU to energy
% This energy axis needs to be shifted to account for the scatter using
% Compton kinematics
Ef = x_kev;
theta = 90;
a=(Ef/511)*(1-cosd(theta));
b= 1-a;
Ei = Ef./b;


%hold on
%plot(Ef, thisSpectrumSt,'b');
%plot(Ef, thisSpectrumnofilt);
%plot(Ei, thisSpectrumSt,'r');
%plot(Ei, thisSpectrumnofilt,'g');; %all three graphs
%xlabel('Energy (keV)');
%ylabel('No. of counts')
%legend('Uncorrected X-ray spectrum with ANSI 300 filter','Uncorrected X-ray spectrum without filter', 'Corrected X-ray spectrum with ANSI 300 filter','Corrected X-ray spectrum (no filter)')
%title('X-ray spectra collected')
%saveas(gcf,'Allspectracomparison.png')


hold on
plot(Ef, thisSpectrumSt,'b');
plot(Ei, thisSpectrumSt,'r');
plot(Ef, thisSpectrumnofilt);
plot(Ei, thisSpectrumnofilt,'g');
xlim([7 70])
xlabel('Energy (keV)');
ylabel('No. of counts')
legend('Uncorrected X-ray spectrum w/ filter', 'Corrected X-ray spectrum w/filter','Uncorrected X-ray spectrum with no filter','Corrected X-ray spectrum with no filter', '65.68keV')
xline(65.68,'.')
xline(59.321,'.')
title('Spectra before & after Compton kinematics corrections')
%saveas(gcf,'Zoomedspectracomparison3.png');


%% Simulate spectrum 
% you know the kV from above
%Simulating spectra
simSpk = spektrSpectrum(65.68,[0,0],'TASMIP',0);
simSpk_int = spektrBeers(simSpk,[4 0.12]);




%adding Iron filter 74.675%
simSpk_Fe = spektrBeers(simSpk_int,[26 0.373375]);
%adding Chromium filter 16%
simSpk_Cr = spektrBeers(simSpk_Fe,[24 0.08]);
%adding Nickel filter 6%
simSpk_Ni=spektrBeers(simSpk_Cr,[28 0.03]);
%Adding Mn filter 2%
simSpk_Mn=spektrBeers(simSpk_Ni,[25 0.01]);
%Adding Si filter 1% 
simSpk_Si=spektrBeers(simSpk_Mn,[14  0.005]);
%Adding Nitrogen 0.1%
simSpk_N=spektrBeers(simSpk_Si,[7 0.0005]);
%Adding Carbon 0.15%
simSpk_C=spektrBeers(simSpk_N,[6  0.00075]);
%Adding Phosphorous 0.045%
simSpk_P=spektrBeers(simSpk_C,[15  0.000225]);
%Adding Sulfur filter 0.03%
simSpk_ANSI=spektrBeers(simSpk_P,[16 0.00015]);
figure
plot(simSpk_int);
hold on
plot(simSpk_ANSI);
plot(Ei, thisSpectrumnofilt);
plot(Ei, thisSpectrumSt);
xlabel('Energy (keV)')
ylabel('No. of counts')
xlim([0 100]);
legend('Simulated spectrum with inherent filtration', 'Simulated spectrum with Steel Filter', 'Acquired spectrum with inherent filtration','Acquired spectra with Steel Filter');
title('SPEKTR simulated spectra');
saveas(gcf,'SPEKTRsimszoom.png');

%% Energy resolution
% adjust the simulated spectrum to account for the energy resolution of the
% detector - you proabaly need a loop and some kind of function
%Idea: Simulation currently putsenergy values into 1keV bins, to adjust to
%the detector resolution these bins would need to be 4.19keV
%Modifying the energy resolution:
%Need to change 1keV resolution to 4.19keV resolution:
newxvalues=linspace(0,150,(150/4.19));
no=length(newxvalues);

newrows=zeros(35,4);
newrows(1,:)=simSpk_ANSI(1:end/35,1);
for i=1:34
    x=simSpk_ANSI(i*(end/35)+1:(i+1)*(end/35));
    newrows(i+1,:)=x;
    
end
% Summing up the total counts in the new energy bins;
newrowssum=sum(newrows,2);
ave=newrowssum./4;
plot(newxvalues,ave,'m');

%For the unfilteres spectrum:
newrows2=zeros(35,4);
newrows2(1,:)=simSpk_int(1:end/35,1);
for i=1:34
    x=simSpk_int(i*(end/35)+1:(i+1)*(end/35));
    newrows2(i+1,:)=x;
    
end
% Summing up the total counts in the new energy bins;
newrowssum2=sum(newrows2,2);
ave2=newrowssum2./4;
plot(newxvalues,ave2,'c');
legend('Simulated spectrum with inherent filtration', 'Simulated spectrum with Steel Filter', 'Acquired spectrum with inherent filtration','Acquired spectrum with Steel Filter', 'Sim. spectrum w/ filter, corrected for resolution','Sim. spectrum no filter, corrected for resolution' );
saveas(gcf,'SPEKTRsimszoomall2.png');


%%

55.6

