% A Matlab script for processing TRAST monitoring data (File 1 of 2)
% ------------------------------------------- 
% A. Karampatzakis and T. Wohland
% -------------------------------------------
% Centre for Bioimaging Sciences, National University of Singapore, Singapore 117557 
%
% Downloaded from: http://www.dbs.nus.edu.sg/lab/BFL/ 
% September 2016
% 
% This script reads the contents of a folder where the images (.tif) taken
% with increasing duty cycles are saved. The filename of each image needs
% to resemble the duty cycle under which was taken. For example, if 4% duty
% cycle was used, the filename should be 4.tiff. The script visualises the
% raw intensity curves of each pixel and saves the variables as [sample].mat
%
% The saved .mat file can be processed by process_data.m
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

clear
close all

sample = 'myMeas';      %Set the name of your sample. Must be same as your folder.

period=50E-6;           %Set here the pulse train period used in your TRAST experiment, in seconds.

c_dir=[sample];         %The script enters the folder where files are saved
cd(c_dir)

files=dir('*.tif');     %Change accordingly for .tiff etc
filenames={};

%The script reads all .tif files into a list and builds a stack (traststack)

for i=1:length(files)   
    filenames{i}=files(i).name;
end

filelist=[];
sortedfilenames={};

for i=1:length(files);  
    current_file=filenames{i};
    filelist(i)=str2num(filenames{i}(1:end-3)); 
    %NOTE: For reading .tif files, above should be end-3. For .tiff, end-4
end

sortedfilelist=sort(filelist);

for i=1:length(files);
    sortedfilenames{i}=strcat(num2str(sortedfilelist(i)),'.tif');
end

no_of_tp=length(sortedfilenames);

for i=1:no_of_tp                        
    
    info = imfinfo(sortedfilenames{i});
    current_file=sortedfilenames{i};
    
    A = imread(current_file, 1, 'Info', info);
  
    x=size(A,1);
    y=size(A,2);
    
    traststack(:,:,i)=A;                %Build the stack
    
end
  
w=[];                                   %Build the x-axis
duty_cycle=sortedfilelist/100;          %Converts file names to [0-1] range
frequency=1/period;                     %In Hz

for i=1:length(sortedfilenames)         %Convert axis in w (s)
    cur_w=(1./frequency)*duty_cycle(i);
    f(i)=frequency;
    T(i)=1./frequency;
    w(i)=cur_w;
end

% Visualising the data

figure(1)
fontsize=15;
plot(w(1:end),reshape(traststack,size(traststack,1)*size(traststack,2),size(traststack,3))','x--','markersize',15)
xlabel('Pulse duration (s)','fontsize',fontsize)
ylabel('Fluorescence intensity (AU)','fontsize',fontsize)
grid on
set(gca,'FontSize',fontsize)
title([sample], 'Interpreter','none');
xlim([0, 1.01*T(1)]);

% Save the data that will need to be imported in the next script

save([sample],'w','traststack','period')

