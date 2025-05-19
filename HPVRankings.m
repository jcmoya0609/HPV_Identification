function [RSS_HPV_list, RSS_values, MW_HPV_list, MW_values] = HPVRankings(RS2C,m,b, sigma)
    % RS2C:     - 3x3 matrix Sample to crystal grain averaged orientation

    % m:        - n x 3 array with n being the number of possible HPV. m is the 
                % plane normal vector associated with each HPV.  eg. for a
                % CuAlNi cubic to orthorhombic there are 96 HPVs

    % b:        - n x 3 array with n being the number of possible HPV. b is the 
                % shape strain vector associated with each HPV.

    % sigma:    - 3 x 3 x 3 3D matrix that contains the macroscopic stress tensor applied to
                % the grain in the lab frame 

F_sample = cell(1,96);
F_crystal = cell(1,96);    
schmid = cell(1,96);
schmid2 = cell(1,96);
transform_strain = cell(1,96);
work_transform=zeros(1,96);
tau_r=zeros(1,96);
sigma_prime = RS2C*sigma*RS2C';
HPVNum=1:96;
for ii = 1 : 96

    btmp = b(ii,:)';
    mtmp = m(ii,:)';
    btmpunit = (btmp) / (norm(b(ii,:)));

    schmid{ii} = btmp * transpose(mtmp);
    F_crystal{ii} = eye(3) + schmid{ii};
    F_sample{ii} = RS2C' * F_crystal{ii} * RS2C;

    transform_strain{ii} = 1/2 * (F_sample{ii}' * F_sample{ii} - eye(3));
    work_transform(ii) = trace(sigma * transform_strain{ii}');
    % work_transform(ii) = transform_strain{ii}(3,3);
    
    schmid2{ii} = mtmp * transpose(btmpunit);
    tau_r(ii) = trace(schmid2{ii} * sigma_prime');

end

[MW_values, ind] = sort(work_transform, 'descend');
MW_HPV_list = HPVNum(ind);

[RSS_values, ind2] = sort(tau_r, 'descend');
RSS_HPV_list = HPVNum(ind2);


end