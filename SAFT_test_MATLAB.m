% SAFT
% based on Thomsom 1984
  
clear all
format long

file = load("DadosEnsaio.mat");

idx_0 = 501; %gate start
idx_1 = 900; %gate end
b_scan = file.ptAco40dB_1.AscanValues(idx_0:idx_1,:); %bscan data
cl = 5836.575875486382; %material velocity
t = file.ptAco40dB_1.timeScale(idx_0:idx_1)*1e-6; %time grid
z = (cl*t)/2; %distance travelled
x = file.ptAco40dB_1.CscanData.X*1e-3; %transducer x axis position

f = zeros(size(b_scan));
for i=1:size(b_scan,1)
    for j=1:size(b_scan,2)
        coord_x = x(j); 
        coord_z = z(i);
        for pulse_echo=1:(size(b_scan,2))
            coord_trans_x = x(pulse_echo);
            coord_trans_z = 0;

            delta_x = coord_trans_x - coord_x;
            delta_z = coord_trans_z - coord_z;
            delta_x_pow = delta_x*delta_x;
            delta_z_pow = delta_z*delta_z;

            dist = sqrt(delta_x_pow + delta_z_pow);
        
            [subt, idx] = min(abs(z-dist));
            f(i,j) = f(i,j) + b_scan(idx,pulse_echo);
        end
    end  
end
f = f/size(b_scan,2);  

subplot(1,2,1);
image(b_scan,'CDataMapping','scaled')
title ('B-Scan')
subplot (1,2,2);
image(f,'CDataMapping','scaled')
title ('SAFT')
