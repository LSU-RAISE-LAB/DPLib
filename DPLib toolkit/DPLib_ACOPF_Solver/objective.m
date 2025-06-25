function f = objective(x, auxdata)
    % Base generation cost
    mpc = auxdata{1};
    nbuses = size(mpc.bus, 1);
    ngens = size(mpc.gen, 1);
    baseMVA = mpc.baseMVA;
    Pg = 2*nbuses + (1:ngens);
    if ngens~=0
    actgen = mpc.gen(:, 8);

    if mpc.gencost(1, 4) == 3
        f = sum(actgen .* (baseMVA^2 * mpc.gencost(:, 5) .* x(Pg).^2 + baseMVA * mpc.gencost(:, 6) .* x(Pg) + mpc.gencost(:, 7)));
    elseif mpc.gencost(1, 4) == 2
        f = sum(actgen .* (baseMVA * mpc.gencost(:, 5) .* x(Pg) + mpc.gencost(:, 6)));
    end
    else
        f = 0;
    end
    % Additional terms from tie-line consistency enforcement
    global rho PuScale landa lambdav regnum loadedData
    regionTag = ['R' num2str(regnum)];

    % Extract interregional tie-lines table
    tieTableKey = ['interregional_tielines' regionTag];
    dataMatrix = loadedData.(tieTableKey);

    % Sort tie-line pairs by alphabetical region order
for i = 1:size(dataMatrix, 1)
    r1 = dataMatrix{i, 1};  % e.g., 'mpc_regionR3'
    r2 = dataMatrix{i, 2};  % e.g., 'mpc_regionR11'

    % Extract full numeric region numbers
    r1_num = sscanf(r1, 'mpc_regionR%d');
    r2_num = sscanf(r2, 'mpc_regionR%d');

    % Reorder based on numeric comparison
    if r1_num > r2_num
        dataMatrix{i, 1} = r2;
        dataMatrix{i, 2} = r1;
    end
end

    Tablek = cell2table(dataMatrix);

    % Load region-specific values from loadedData
    con      = loadedData.(['con' regionTag]);
    Yon      = loadedData.(['Yon' regionTag]);
    region_buses = loadedData.(['region_' regionTag]);
    nitr      = loadedData.(['nit' regionTag]);
    nitall   = loadedData.('nit');

    for i = 1:nitr
        str1 = char(Tablek{i, 1});
        str2 = char(Tablek{i, 2});

region1_num = sscanf(str1, 'mpc_regionR%d');
region2_num = sscanf(str2, 'mpc_regionR%d');
localRegion_num = regnum;

if localRegion_num == region1_num
    remoteRegion_num = region2_num;
else
    remoteRegion_num = region1_num;
end

remoteRegionStr = ['R' num2str(remoteRegion_num)];


remotePhase = loadedData.(['phase_angles_' remoteRegionStr]);
remoteMag   = loadedData.(['voltage_mag_' remoteRegionStr]);
remote_nit  = loadedData.(['nit' remoteRegionStr]);


        y = Yon(i);  % index mapping

        if localRegion_num == region1_num
            grad_theta      = x(region_buses{i} + nbuses)              - remotePhase(remote_nit + con(i));
            grad_theta_dual = x(2*nbuses + 2*ngens + 2*i)           - remotePhase(con(i));
            grad_mag        = x(region_buses{i})                   - remoteMag(remote_nit + con(i));
            grad_mag_dual   = x(2*nbuses + 2*ngens + 2*i - 1)       - remoteMag(con(i));
        else
            grad_theta      = remotePhase(con(i))                  - x(2*nbuses + 2*ngens + 2*i);
            grad_theta_dual = remotePhase(remote_nit + con(i))     - x(region_buses{i} + nbuses);
            grad_mag        = remoteMag(con(i))                    - x(2*nbuses + 2*ngens + 2*i - 1);
            grad_mag_dual   = remoteMag(remote_nit + con(i))       - x(region_buses{i});
        end
if ngens~=0
        % Augmented Lagrangian accumulation
        f = f + PuScale * ( ...
              landa(abs(y))            * grad_theta      + (rho / 2) * grad_theta^2 + ...
              landa(abs(y) + nitall)      * grad_theta_dual + (rho / 2) * grad_theta_dual^2 + ...
              lambdav(abs(y))          * grad_mag        + (rho / 2) * grad_mag^2 + ...
              lambdav(abs(y) + nitall)    * grad_mag_dual   + (rho / 2) * grad_mag_dual^2 );
else
        f = f + ( ...
              landa(abs(y))            * grad_theta      + (rho / 2) * grad_theta^2 + ...
              landa(abs(y) + nitall)      * grad_theta_dual + (rho / 2) * grad_theta_dual^2 + ...
              lambdav(abs(y))          * grad_mag        + (rho / 2) * grad_mag^2 + ...
              lambdav(abs(y) + nitall)    * grad_mag_dual   + (rho / 2) * grad_mag_dual^2 );
end
    end

