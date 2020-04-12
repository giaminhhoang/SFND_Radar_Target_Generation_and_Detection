# SFND_Radar_Target_Generation_and_Detection

The main idea of the CFAR algorithm is to dynamically ajust the threshold in order to detect the signals of interest in the presence of noise. 

For each cell under test (CUT), the signals in the neigboring cells are averged to find the threshold. In order to prevent the impact of the target on the noise estimation, guard cells are introduced. The neighboring cells are used to estimate the noise are called traning cells.

The number of traning cells should be large enough to estimate the noise level.

The number of guard cells are chosen based on the resolution of Range Dopler Map (RDM) whose 1 range cell has resolution of about 1 m whereas 1 dopler cell has resolution of 1.5748 m/s. 

Futhermore, an offset is also introduced to mitigate false alarm. The offset is derived from the SNR. From the Range Dopler Map (RDM), the SNR is more or less 15 dB. Therefore, 10 dB offset is chosen in order not to miss detection.

```
%Select the number of Training Cells in both the dimensions. 
Tr = 10;
Td = 10;

%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation

Gr = 4;
Gd = 4;

% offset the threshold by SNR value in dB
offset_dB = 10;
```

The 2D CFAR process is implemented in the following code block.
It takes RDM as the input and processes as follows:

Lines 31 and 32 show the nested loop to check all the cells the 2D RDM.

Lines 35 - 37 set the cells at the edge to zeros as they do have enough tranning and guard cells to apply CFAR.

Lines 41 and 42 supports to calcualte the average signal value over the training cells in line 45. I need to calculate the sum of all values in the trainng cells. Here, since the sum of all training cells, guard cells and CUT and the sum of only guard cells and CUT are easily computed, I then perform the subtraction to get the sum of all trainning cells only. 

Line 49 implements a threshold for a CUT with a predefined offset in order to mitigate false alarm. The offset is derived from the SNR. Note that as it is in dB we need to convert the averaged signal in trainning cells to dB.

Finally, lines 53 - 57 implements simple detection by comaparing the value in CUT with the calculated threshold.

```
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
```
