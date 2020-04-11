clear all
close all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
target_range_init = 110;
target_velocity = -20; 


%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
range_max = 200;        
range_resolution = 1;
velocity_max = 100;
c = 299792458;      % speed of light

B = c/2/range_resolution;
Tchirp = 5.5*2*range_max/c;
slope = B/Tchirp;
disp(['Slope is ', num2str(slope)]);

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq

                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));
target_range = zeros(1,length(t));  % range 

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity.
    if i > 1
        target_range(i) = target_range(i-1) + target_velocity*(t(i) - t(i-1));
    else
        target_range(i) = target_range_init;
    end
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    delay = 2*target_range(i)/c;
    Tx(i) = cos(2*pi*(fc*t(i) + 0.5*slope*t(i)^2 ) );
    Rx(i) = cos(2*pi*(fc*(t(i) - delay) + 0.5*slope*(t(i)-delay)^2 ) );
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);
    
end

%% RANGE MEASUREMENT


 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix = reshape(Mix, Nr, Nd);

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
chirp_fft_1 = fft(Mix(:,1));

 % *%TODO* :
% Take the absolute value of FFT output
chirp_spec_double_1 = abs(chirp_fft_1/Nr);

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
chirp_spec_single_1 = chirp_spec_double_1(1:Nr/2+1);
chirp_spec_single_1(2:end-1) = 2*chirp_spec_single_1(2:end-1);

%plotting the range
figure ('Name','Range from First FFT')
% subplot(2,1,1)

 % *%TODO* :
 % plot FFT output 
Fs = 1/Tchirp*Nr;       % sampling frequency
axis_f = Fs*(0:(Nr/2))/Nr;  % x-axis in frequency (Hz)
axis_range = axis_f*c*Tchirp/2/B; % convert x-axis in frequency (Hz) to range estimate (m)

% plot(axis_f, chirp_spec_single_1) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')

plot(axis_range, chirp_spec_single_1);
xlabel('range (m)')
ylabel('amplitude')
axis ([0 200 0 1]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr = 10;
Td = 10;

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 4;
Gd = 4;

% *%TODO* :
% offset the threshold by SNR value in dB
offset_dB = 10;

% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);


% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.

n_train_cells = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1) - (2*Gr+1)*(2*Gd+1);


   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR
   
dim_range = size(RDM,1);
dim_dop = size(RDM,2);
detection_map = zeros(dim_range,dim_dop);    % initialize the detection map

for r = 1:dim_range
    for d = 1:dim_dop
        
        if (r <= Tr+Gr || r >= dim_range-Tr-Gr || ... % CUT cannot be located at the edges of matrix thus set to zero
            d <= Td+Gd || d >= dim_dop-Td-Gd )
            detection_map(r,d) = 0;
        else
            % compute the sum of FFT singals at training cells 
            sig_fft2_train_sum = sum(db2pow(RDM(r-Tr-Gr:r+Tr+Gr, d-Td-Gd:d+Td+Gd)),'all') -...
                sum(db2pow(RDM(r-Gr:r+Gr, d-Gd:d+Gd)),'all');
            
            % compute the average signal over the training cell
            sig_fft2_train_avg = sig_fft2_train_sum/n_train_cells;
            
            % compute the threshold in dB
            threshold = pow2db(sig_fft2_train_avg) + offset_dB;
            
            % detect the target
            if (RDM(r,d) >= threshold)
                detection_map(r,d) = 1;
            else
                detection_map(r,d) = 0;
            end
        end
    end   
end



% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 








% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,detection_map);
colorbar;


 
 