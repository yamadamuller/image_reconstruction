% TFM
% Based on Holmes 2005

format long

file = load("dados_tfm.mat");

ascan = file.ascan;
Cl = file.cl;
f_sampl = double(file.f_sampling*1e6);
t_init = double(file.period*1e-6);
X = transpose(file.x*1e-3);

t = zeros(size(ascan,1),1);
for linha=1:size(t,1)
    t(linha) = (t_init + (linha/f_sampl));
end

Z = (Cl*t)/2;

b_scan = hilbert(ascan(:,:,:));

f = zeros(869,51);
for i=37:905
    for j=14:64
        coord_x = X(j); 
        coord_z = Z(i);
        for emissor=1:(size(b_scan,2))
            for receptor=1:(size(b_scan,3))
                coord_emissor_x = X(emissor);
                coord_emissor_z = 0;
                coord_receptor_x = X(receptor);
                coord_receptor_z = 0;

                delta_e_x = coord_emissor_x - coord_x;
                delta_e_z = coord_emissor_z - coord_z;
                delta_r_x = coord_receptor_x - coord_x;
                delta_r_z = coord_receptor_z - coord_z;

                delta_e_x_pow = delta_e_x^2;
                delta_e_z_pow = delta_e_z^2;
                delta_r_x_pow = delta_r_x^2;
                delta_r_z_pow = delta_r_z^2;

                dist1 = sqrt(delta_e_x_pow + delta_e_z_pow);
                dist2 = sqrt(delta_r_x_pow + delta_r_z_pow);
                dist = dist1 + dist2;
                t_e_r = dist/Cl;
                
                [subt, idx] = min(abs(t-t_e_r));
                f(i-36,j-13) = f(i-36,j-13) + b_scan(idx,emissor,receptor);
            end
        end
    end  
end
f = abs(f/size(b_scan,2));  

image(f,'CDataMapping','scaled')
title ('TFM')
