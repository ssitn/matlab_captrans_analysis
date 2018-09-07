%%%%%%%%% Choose your CapTrans files here %%%%%%%

[flnm,path] = uigetfile('*abf','Multiselect','on');
cd(path)
% script won't work with just one file, since variable flnm has to be
% variable type cell. If you chose just one file then flnm will be variable
% type char and script will throw an error
%% Actuall fitting happens here

for kk = 1:size(flnm,2)
    Vpulse = -5;                                                            % voltage pulse size in mV 
    filename = [path flnm{kk}];                                             % creates full path to the file
[a,b,h] = abfload2(filename);                                               % loads data (requires function pvpmod.m, make sure it is added to Matlab PATH), a - matrix containing recorded data NxMxL, where N - number of samples in each sweep, M usually equals 2 (data for I and data for V), L - number of sweeps. Example if you want to see I from the second sweep than it will be in a(:,1,2). b - frame duration in microseconds, h - info about the file 
I = a(:,1);                                                                 % creates separate variable for I from the current file
V = a(:,2);                                                                 % creates separate variable for V from the current file
mul = 20/b;                                                                 % multiplicator which is used to correct for difference in the sampling rate. For example if you want to get averaged data over 1 ms from the recording at 50kHz you'll need to average 50 sampling points (sampling interval is 20us and mul will be equal 1), if your data is recorded at 100kHz your sampling rate is 10us and mul is 2 in this case, you will average 50x2=100 sampling points, i.e. exactly 1 ms as in previous case.
t = (cumsum(ones(size(a,1),1))-1)*b/1000;                                   % creates an array of time poins (in ms) corresponding to each data sampling point. It always starts from 0ms
I = I -mean(I(1:mul*20));                                                   % subtracts baseline current so your current pulse always starts from 0pA
[peak(1), locs(1)] = findpeaks(-I, 'MINPEAKHEIGHT', 0.9*max(-I), 'NPEAKS', 1); % finds the first peak of on the CapTrans recording
[peak(2), locs(2)] = findpeaks(I, 'MINPEAKHEIGHT', 0.9*max(I), 'NPEAKS', 1);   % finds the second peak of the recording
if size(a,2) >1 
    [peak(3), locs(3)] = findpeaks(-diff(V), 'MINPEAKHEIGHT', 0.9*max(-diff(V)), 'NPEAKS', 1); % finds the begining of the voltage pulse
    locs(3) = locs(3)+1;
else
    [peak(3), locs(3)] = findpeaks(-diff(I), 'MINPEAKHEIGHT', 0.9*max(-diff(I)), 'NPEAKS', 1); % finds the begining of the pulse if voltage was not recorded 
end

I = I - mean(I(1:locs(1)-mul*20));                                          % more precise baseline subtraction
Isteady = mean(I(locs(2)-mul*500:locs(2)-mul*10));                          % measures strady-state current to calculate Rm+Rs

fitType = 'A*exp(-(x-x1)/tau)+C';                                           % set a function type which will be fitted to the data
rb=mul*20;                                                                  % set the right border of fitting (counting from position of the first peak)
lb = 1;                                                                     % set the left border of fitting (counting from position of the first peak)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here script will do series of fittings at different right border and
% record ggodness of fit data. Right border which leads to highest goodness
% of fit will be used to aqcuire the final result

 ii = 1;
while rb<mul*800                                                            % loop will stop when the right border is 16 ms away from the peak (typical decay time of the peak is around 1 ms)
[fitcoef1,gof1] = fit(t(locs(1)+lb:locs(1)+rb),I(locs(1)+lb:locs(1)+rb),fitType, 'StartPoint', [-peak(1) Isteady 0.5 t(locs(1))]); % fit the function defined previously, record parameters of fit and goodness of fit
rb = rb+mul*10;                                                             % increase right border by 10 sampling intervals
gofdata1(ii) = gof1.rsquare;                                                % create the array of goodness of fits and right border values
rbdata(ii) = rb;                                                            %
ii = ii+1;
end
%gofdata1(1)=0;
gofdata_sm = smooth(gofdata1,5);                                            % smoothens goodness of fit data to prevent picking max goodness of fit when it is actually an outlier
gofdata_sm(1)=0;                                                            % brings the first value to 0 as outliers often corresponds to few first rb values
[maxgof, ind] = max(gofdata_sm);                                            % determines max of goodness of fit
rb_max = rbdata(ind);                                                       % finds corresponding rb value

%%%%%%%%%%%%%%%%%%%%%%%%%%

[fitcoef1,gof1] = fit(t(locs(1)+lb:locs(1)+rb_max),I(locs(1)+lb:locs(1)+rb_max),fitType, 'StartPoint', [-peak(1) Isteady 0.5 t(locs(1))]); % actuall final fitting happens here
coeff=coeffvalues(fitcoef1);                                                % gets fitting coefficients: coeff(1) - A, coeff(2) - C, coeff(3) - tau, coeff(4) - x1
%%%%% this part records differet fitting parameters to variable fitparams
fitparams{kk}.filename      = flnm{kk};
fitparams{kk}.fit_function  = fitType;
fitparams{kk}.A_fit         = coeff(1);
fitparams{kk}.tau           = coeff(3);
fitparams{kk}.C_fit         = coeff(2);
fitparams{kk}.A_extrapolate = coeff(1)*exp(-(t(locs(3)-2)-coeff(4))/coeff(3))+coeff(2);
fitparams{kk}.t1            = t(locs(1)+lb);                                % left border in ms
fitparams{kk}.t2            = t(locs(1)+rb_max);                            % right border n ms
fitparams{kk}.gofdata       = gofdata1;
fitparams{kk}.rbdata        = rbdata;
fitparams{kk}.Isteady       = Isteady;
fitparams{kk}.Rs_peak       = -peak(1);
% here is calculation of Rs, Rm and Cm, which are also saved in fitparams
% variable
Rs = Vpulse/(0.001*fitparams{kk}.A_extrapolate);
Rm = Vpulse/(0.001*fitparams{kk}.Isteady);
Cm = 1000*fitparams{kk}.tau*(Rs+Rm)/(Rs*Rm);
fitparams{kk}.Rs_Rm_Cm(:,1) = Rs;
fitparams{kk}.Rs_Rm_Cm(:,2) = abs(Rm);
fitparams{kk}.Rs_Rm_Cm(:,3) = Cm;
fitparams{kk}.x1            = coeff(4);
%%%%%%%%%%%%%%%%%%

% Here is visualisation of the fitting, which is saved as .fig file so you
% can check it manually if needed
figure(123);
plot(t,I)                                                                   % actual recorded data
hold on
plot(t,coeff(1)*exp(-(t-coeff(4))/coeff(3))+coeff(2),'r')                   % fitted exponent

line([t(locs(3)-2) t(locs(3)-2)],[50 -1.25*peak(1)], 'Color', [0 0 0])      % start of the pulse to which exponent was extrapolated for Rs calculation
title(fitparams{kk}.filename)

axis([0 t(mul*800) -1.25*peak(1) 100])
savefig([fitparams{kk}.filename '.fig'])                                    % saving the figure with the same name as the recording file
close 123

end
%% here you can save fitted parameters for further reference
save('fitparams.mat', 'fitparams')

%% here you can create excel file containing filename of the recording and Rs Rm Cm data
start_ind = 0;

for jj= 1:size(fitparams,2)
   xlswrite([pwd '\Rs_Rm_Cm.xls'], {fitparams{jj}.filename}, 1, ['A' num2str(jj+start_ind)]);
   xlswrite([pwd '\Rs_Rm_Cm.xls'], fitparams{jj}.Rs_Rm_Cm, 1, ['B' num2str(jj+start_ind)]);
end
%% Here you can display goodness of fit vs right borrder for the required recording file
jj=2;
figure
plot(fitparams{jj}.rbdata, fitparams{jj}.gofdata)


