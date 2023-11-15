% A Matlab script for processing TRAST monitoring data (File 2 of 2)
% -------------------------------------------
% A. Karampatzakis and T. Wohland
% -------------------------------------------
% Centre for Bioimaging Sciences, National University of Singapore, Singapore 117557
%
% Downloaded from: http://www.dbs.nus.edu.sg/lab/BFL/
% September 2016
%
% This script loads the .mat file created by "prepare_stack.m" and
% calculates the triplet relaxation and intersystem crossing times in a
% pixel-wise manner. The user has the option to bin bixels and investigate
% specific ROIs. The value of interesystem crossing can either be held or
% free-fitted. Experimental parameters and starting values for the fitting
% are defined by the user.
%
% Disclaimer: The software and data are provided for personal or academic
% use only and may not be used in any commercial venture or distributions.
% All files have been virus scanned, however, for your own protection; you
% should scan these files again. You assume the entire risk related to your
% use of this software and data. By using the software and data on this
% site your expressly assume all risks of data loss or damage alleged to
% have been caused by the software and data. The Biophysical Fluorescence
% Laboratory at NUS is providing this data "as is," and disclaims any and
% all warranties, whether express or implied, including (without
% limitation) any implied warranties of merchantability or fitness for a
% particular purpose. In no event will the Biophysical Fluorescence
% Laboratory at NUS and/or NUS be liable to you or to any third party for
% any direct, indirect, incidental, consequential, special or exemplary
% damages or lost profit resulting from any use or misuse of this software
% and data.

close all
clear

sample = 'myMeas';      %Set the name of your sample. Must be same as your folder.

% Do you want to evaluate tt and tisc in a ROI? If set to 1, at the end of fitting
% there will be a popup that lets you create a ROI by making a mask.
% This functionality uses 'roipoly' - type: help roipoly for more info
WantROI = 1;

% Do you want to bin the image? This is optional. You can set to =1 for no reduction.
% This functionality uses function 'imresize'. A value of 0.03 will reduce pixels
% to 3% of original etc. - type: help imresize for more info
binscale = 1;

% The problem parameters
G = 2E-6;                   %Local concentration of fluorophores, in M.

% The below are used to calculate power density for SPIM.
% You can override I_density_laser for other types of setups.
I_laser=23E-3;              %Laser power in W
h = 6.62606957E-34;         %Planck constant
c = 299792458;              %Speed of light in m/s
lambda = 488E-9;            %Laser wavelength in m
sigma=1.7E-20;              %Abs. cross section in m2
beamwidth = 1.8E-6;         %Beam width (1/e2 value) in m
beamheight = 200E-6;        %Beam height in m
A = beamwidth*beamheight;   %Beam crosssection in m2

I_density_laser = I_laser/A; % Laser power density in W/m2

% The starting (or fixed) values for fitting
tt = 4E-6; %Triplet relaxation time (s)
tisc = 3E-9; %Intersystem crossing time (s)
t10 = 1.55E-9; %Fluorescence lifetime (s)
tex = h*c/(lambda*I_density_laser*sigma); %Excitation time (s)
eta = 1; %Conversion factor between fluorescence and digital readout

% Do you want to hold t_isc? // 1 means "hold it" to the starting value set above
hold_tisc = 0;

% Load the previously saved .mat
cd ([sample]);
eval(['load ' sample ]);

% Working parameters
Tp = period;
w = w(:);

% Set the range of values for fitting
eta_low=0.5;
eta_high=3;

tisc_low=1E-9;
tisc_high=5E-8;

tt_low=1.75E-6;
tt_high=1E-4;

%% --------- No editing is required in this section

eta_start=eta;
tt_start=tt;
tisc_start=tisc;

% Re-sampling the fit curve to make it smooth
ww = linspace(w(1),w(end),200);

% Re-shape the array
stacksize=length(w);

for i=1:stacksize
    trast_stack_3d(:,:,i)=imresize(traststack(:,:,i),binscale,'Method','nearest') ;
end

trast_stack_2d=reshape(trast_stack_3d,[size(trast_stack_3d,1)*size(trast_stack_3d,2),stacksize]);

% Fit for each pixel in the reshaped stack
for i=1:size(trast_stack_2d,1)
    
    pixel_int = im2double(trast_stack_2d(i,:)');
    pixel_int = pixel_int/pixel_int(1); %Normalise the plots
    
    if (hold_tisc==0)
        
        lowervalues=[eta_low  tisc_low  tt_low];
        uppervalues=[eta_high  tisc_high tt_high];
        startingpoints=[eta_start tisc_start tt_start];
        
        fo_ = fitoptions('MaxIter',1000,'MaxFunEvals',2000, 'method','NonlinearLeastSquares','Lower',lowervalues,'Upper',uppervalues);
        
        st_ = startingpoints;
        set(fo_,'Startpoint',st_);
        ft_ = fittype('(eta * G * (((1/t10)/((1/t10)+(1/tisc))) / (tex+t10+(t10/tisc)*tt)))*(1 - (tt * (((t10/tisc)*tt)/(tex+t10+(t10/tisc)*tt)) * (1-exp(-(Tp-tp)/tt)) / (1-exp(-(Tp-tp)/tt + (-((tex+t10+(t10/tisc)*tt)/(tt*(tex+t10))))*tp))) * ((exp((-((tex+t10+(t10/tisc)*tt)/(tt*(tex+t10))))*tp)-1)/tp))',...
            'dependent',{'Iout'},'independent',{'tp'},...
            'coefficients',{'eta','tisc','tt'},'problem',{'Tp', 'G','tex','t10'});
        
        [cf_,gof, output] = fit(w,pixel_int,ft_,fo_,'problem',{Tp,G,tex,t10});
        
    else
        
        lowervalues=[eta_low    tt_low];
        uppervalues=[eta_high    tt_high];
        startingpoints=[eta_start    tt_start];
        
        fo_ = fitoptions('MaxIter',1000,'MaxFunEvals',2000, 'method','NonlinearLeastSquares','Lower',lowervalues,'Upper',uppervalues);
        
        st_ = startingpoints;
        set(fo_,'Startpoint',st_);
        ft_ = fittype('(eta * G * (((1/t10)/((1/t10)+(1/tisc))) / (tex+t10+(t10/tisc)*tt)))*(1 - (tt * (((t10/tisc)*tt)/(tex+t10+(t10/tisc)*tt)) * (1-exp(-(Tp-tp)/tt)) / (1-exp(-(Tp-tp)/tt + (-((tex+t10+(t10/tisc)*tt)/(tt*(tex+t10))))*tp))) * ((exp((-((tex+t10+(t10/tisc)*tt)/(tt*(tex+t10))))*tp)-1)/tp))',...
            'dependent',{'Iout'},'independent',{'tp'},...
            'coefficients',{'eta','tt'},'problem',{'Tp', 'G','tex', 'tisc','t10'});
        
        [cf_,gof, output] = fit(w,pixel_int,ft_,fo_,'problem',{Tp,G,tex, tisc, t10});
        
    end
    
    tt_array(i)=cf_.tt;
    tisc_array(i)=cf_.tisc;
    rsq_array(i)=gof.rsquare;
    
    %Show the fitting results
    
    display(['Currently fitting pixel: ' num2str(i) ' of ' num2str(size(trast_stack_2d,1)) ' (' num2str(100*(i/size(trast_stack_2d,1))) '%)'])
    
    figure(1)
    fontsize=13;
    
    subplot(4,1,[1, 3]); %Show the raw data and fit of each pixel
    plot(ww,feval(cf_,ww),'-b',w,pixel_int','.r','MarkerSize',5,'LineWidth',2)
    legend('Fit function','Measurements')
    hold on
    h_ = plot(w, pixel_int, 'o r');
    xlim([min(w) max(w)])
    xlabel('Pulse duration (s)','fontsize',fontsize)
    ylabel('Normalised Intensity','fontsize',fontsize)
    grid on
    set(gca,'FontSize',fontsize)
    
    subplot(4,1,4); %Calculate and show residuals
    res_ = pixel_int - cf_(w);
    h_ = line(w,res_,'Color',[1 0 0],...
        'LineStyle','-', 'LineWidth',1,...
        'Marker','.', 'MarkerSize',6);
    xlim([0 1.0001*Tp])
    
end

%% Visualisation of data - plots

%Below what values you want to exclude from plotting?
%To bypass, set equal to 0
rsq_threshold = 0.9;

figure(2)
s1=subplot(2,2,1);
tt_array(rsq_array<rsq_threshold)=NaN;
tt_array_reshaped=reshape(tt_array,(size(trast_stack_3d(:,:,1))));
b=imagesc(tt_array_reshaped);
set(b,'AlphaData',~isnan(tt_array_reshaped))
h=colorbar;
title('Triplet time (s)')
axis image
colormap(s1,jet)

if hold_tisc == 0
    s1=subplot(2,2,2);
    tisc_array(rsq_array<rsq_threshold)=NaN;
    tisc_array_reshaped=reshape(tisc_array,(size(trast_stack_3d(:,:,1))));
    b=imagesc(tisc_array_reshaped);
    set(b,'AlphaData',~isnan(tt_array_reshaped))
    h=colorbar;
    title('Intersystem crossing time (s)')
    colormap(s1,jet)
    axis image
    %set(h, 'ylim', [2E-9 7.5E-9]);
end

s3=subplot(2,2,3);
rsq_array(rsq_array<rsq_threshold)=NaN;
rsq_array_reshaped=reshape(rsq_array,(size(trast_stack_3d(:,:,1))));
b=imagesc(rsq_array_reshaped);
set(b,'AlphaData',~isnan(tt_array_reshaped))
h=colorbar;
title('R^2')
%set(h, 'ylim', [0.92 1])
colormap(s3,winter)
axis image

subplot(2,2,4);
binranges=linspace(min(min(tt_array)), max(max(tt_array)), 20);
[bincounts,ind] = histc(tt_array,binranges);
bar(binranges,bincounts,'histc')
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0 .5 .5],'EdgeColor','w')
xlabel('Triplet Time (s)');

%% Save the results

save ([sample '_results.mat']);

%% This section lets you define a ROI and gives you the mean and std

if WantROI == 1
    
    figure()
    
    BW = roipoly(mat2gray(tt_array_reshaped));
    BW = double(BW);
    BW(BW==0) = NaN;
    
    masked_tt = tt_array_reshaped.*BW;
    
    display(['Triplet relaxation time in ROI: ' num2str(nanmean(nanmean(masked_tt))) ' +/- ' num2str(nanmean(nanstd(masked_tt))) ]);
    
    if hold_tisc == 0
        masked_tisc = tisc_array_reshaped.*BW;
        display(['Intersystem crossing time in ROI: ' num2str(nanmean(nanmean(masked_tisc))) ' +/- ' num2str(nanmean(nanstd(masked_tisc))) ]) ;
    end
   
end