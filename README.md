# SFND_Radar_Target_Generation_and_Detection

```
%Select the number of Training Cells in both the dimensions.
Tr = 10;
Td = 10;

%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 4;
Gd = 4;

% offset the threshold by SNR value in dB
% From the Range Dopler Map (RDM), the SNR is more or less 15 dB. Therefore, 10 dB offset is chosen in order not to miss detection.
offset_dB = 10;
```

The 2D CFAR process is implemented in the following code block.
It takes RDM as the input. 

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
