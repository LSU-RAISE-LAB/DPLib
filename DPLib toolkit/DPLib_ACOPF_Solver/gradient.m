function g = gradient(x, auxdata)
    mpc = auxdata{1};
    nbuses = size(mpc.bus, 1);
    ngens = size(mpc.gen, 1);
    baseMVA = mpc.baseMVA;
    Pg = 2 * nbuses + (1:ngens);
    if ngens~=0
    actgen = mpc.gen(:, 8);
    end
    % Initialize gradient vector
    global regnum loadedData rho PuScale landa lambdav 
    regionTag = ['R' num2str(regnum)];
    nitr = loadedData.(['nit', regionTag]);nitall   = loadedData.('nit');
    g = zeros(2 * nbuses + 2 * ngens + 2 * nitr, 1);
if ngens~=0
    % Gradient of generation cost
    if mpc.gencost(1, 4) == 3
        g(Pg) = actgen .* (2 * baseMVA^2 * mpc.gencost(:, 5) .* x(Pg) + baseMVA * mpc.gencost(:, 6));
    elseif mpc.gencost(1, 4) == 2
        g(Pg) = actgen .* (baseMVA * mpc.gencost(:, 5));
    end
end
    % Interregional tie-line processing
    tieTableKey = ['interregional_tielines', regionTag];
    dataMatrix = loadedData.(tieTableKey);

    % Sort region pairs by numeric order
    for i = 1:size(dataMatrix, 1)
        r1 = dataMatrix{i, 1};
        r2 = dataMatrix{i, 2};
        r1_num = sscanf(r1, 'mpc_regionR%d');
        r2_num = sscanf(r2, 'mpc_regionR%d');
        if r1_num > r2_num
            dataMatrix{i, 1} = r2;
            dataMatrix{i, 2} = r1;
        end
    end
    Tablek = cell2table(dataMatrix);

    % Region-specific values
    con = loadedData.(['con', regionTag]);
    Yon = loadedData.(['Yon', regionTag]);
    region_buses = loadedData.(['region_', regionTag]);

    for i = 1:nitr
        str1 = char(Tablek{i, 1});
        str2 = char(Tablek{i, 2});
        region1_num = sscanf(str1, 'mpc_regionR%d');
        region2_num = sscanf(str2, 'mpc_regionR%d');
        localRegion_num = regnum;

        if localRegion_num == region1_num
            remoteRegion_num = region2_num;
            sign_y = 1;
        else
            remoteRegion_num = region1_num;
            sign_y = -1;
        end

        remoteTag = ['R' num2str(remoteRegion_num)];
        remotePhase = loadedData.(['phase_angles_', remoteTag]);
        remoteMag = loadedData.(['voltage_mag_', remoteTag]);
        nit_remote = loadedData.(['nit', remoteTag]);

        y = sign_y * Yon(i);

        if sign_y == 1
            grad_theta     = x(region_buses{i} + nbuses)            - remotePhase(nit_remote + con(i));
            grad_theta_dual = x(2 * nbuses + 2 * ngens + 2 * i) - remotePhase(con(i));
            grad_mag       = x(region_buses{i})                 - remoteMag(nit_remote + con(i));
            grad_mag_dual  = x(2 * nbuses + 2 * ngens + 2 * i - 1) - remoteMag(con(i));
        else
            grad_theta     = remotePhase(con(i))                - x(2 * nbuses + 2 * ngens + 2 * i);
            grad_theta_dual = remotePhase(nit_remote + con(i))   - x(region_buses{i} + nbuses);
            grad_mag       = remoteMag(con(i))                  - x(2 * nbuses + 2 * ngens + 2 * i - 1);
            grad_mag_dual  = remoteMag(nit_remote + con(i))      - x(region_buses{i});
        end

        % Apply updates to gradient vector
        idx1 = region_buses{i};
        idx2 = region_buses{i} + nbuses;
        idx3 = 2 * nbuses + 2 * ngens + 2 * i - 1;
        idx4 = 2 * nbuses + 2 * ngens + 2 * i;
if ngens~=0
        if y >= 0
            g(idx2) = g(idx2) + sign_y * PuScale * (landa(abs(y)) + rho * grad_theta);
            g(idx4) = g(idx4) + sign_y * PuScale * (landa(abs(y) + nitall) + rho * grad_theta_dual);
            g(idx1) = g(idx1) + sign_y * PuScale * (lambdav(abs(y)) + rho * grad_mag);
            g(idx3) = g(idx3) + sign_y * PuScale * (lambdav(abs(y) + nitall) + rho * grad_mag_dual);
        else
            g(idx2) = g(idx2) + sign_y * PuScale * (landa(abs(y) + nitall) + rho * grad_theta_dual);
            g(idx4) = g(idx4) + sign_y * PuScale * (landa(abs(y)) + rho * grad_theta);
            g(idx1) = g(idx1) + sign_y * PuScale * (lambdav(abs(y) + nitall) + rho * grad_mag_dual);
            g(idx3) = g(idx3) + sign_y * PuScale * (lambdav(abs(y)) + rho * grad_mag);
        end
else
        if y >= 0
            g(idx2) = g(idx2) + sign_y * (landa(abs(y)) + rho * grad_theta);
            g(idx4) = g(idx4) + sign_y * (landa(abs(y) + nitall) + rho * grad_theta_dual);
            g(idx1) = g(idx1) + sign_y * (lambdav(abs(y)) + rho * grad_mag);
            g(idx3) = g(idx3) + sign_y * (lambdav(abs(y) + nitall) + rho * grad_mag_dual);
        else
            g(idx2) = g(idx2) + sign_y * (landa(abs(y) + nitall) + rho * grad_theta_dual);
            g(idx4) = g(idx4) + sign_y * (landa(abs(y)) + rho * grad_theta);
            g(idx1) = g(idx1) + sign_y * (lambdav(abs(y) + nitall) + rho * grad_mag_dual);
            g(idx3) = g(idx3) + sign_y * (lambdav(abs(y)) + rho * grad_mag);
        end
end
    end

