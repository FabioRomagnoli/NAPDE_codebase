%% REDIMENSIONALIZATION

% concentration_plot = 1;
% potential_plot = 1;
% current_plot = 1;

% extracting vectors 
vAll = X(1:lr,:);
nAll = X(lr + 1:2*lr,:);
pAll = X(2*lr + 1:end,:);

% Descaling procedure
vAll = vAll*Vbar;
nAll = nAll*nbar;
pAll = pAll*nbar;

% Final values
vEnd = vAll(:,end);
nEnd = nAll(:,end);
pEnd = pAll(:,end);

%% CONCENTRATIONS ---------------------------------------------------------
if concentration_plot 
    kPlot = K;
    figure;
    title('Electron and hole concentrations')
    for k=kPlot:K
        clf; % Clear figure before plotting new data
        hold on;

        vk = vAll(:,k);
        nk = nAll(:,k);
        pk = pAll(:,k);

        % semilogy(r,p,'o-',r,n,'x-')
        % semilogy(r,pk, "DisplayName", "p");
        semilogy(r,nk, "DisplayName", "n");
        hold off;

        grid on; 
        legend();
        set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')

        drawnow;
    end
end
%% ELECTRIC POTENTIAL -----------------------------------------------------
if potential_plot 
    kPlot = K;
    figure;
    title('Electric potential')
    for k=kPlot:K
        clf; % Clear figure before plotting new data

        vk = vAll(:,k);
        plot(r,vk, "DisplayName", "v");

        grid on; 
        legend();
        drawnow;
    end

end

%% CURRENT PLOT -----------------------------------------------------------
if current_plot 
    kPlot = K;
    x_medi = (r(1:end-1)+r(2:end))/2;
    figure;

    title('Corrente elettrica')
    for k=kPlot:K
        clf; % Clear figure before plotting new data
        vk = vAll(:,k);
        nk = nAll(:,k);
        pk = pAll(:,k);

        Jn = Comp_current(r,mun,q,vk,Vth,-1,nk);
        Jp = Comp_current(r,mup,q,vk,Vth, 1,pk);
        JJ = Jn+Jp;

        hold on;
        plot(x_medi,Jn, "DisplayName","Jn");
        plot(x_medi,Jp, "DisplayName","Jp");
        plot(x_medi,JJ, "DisplayName","JJ");
        hold off;

        grid on; 
        legend();
        set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')

        drawnow;
    end

end