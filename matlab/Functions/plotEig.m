function [] = plotEig(eigA,zita_min,zita_lim,string_values)
    figure 
    % Group the complex eigenvalues with the same real part and color them similarly
    real_parts = unique(real(eigA(imag(eigA)~=0)),'stable'); %the real parts of the poles that have imaginary components ( numel = number of complex poles) 
    comp_real  = unique(real(eigA(imag(eigA)==0)),'stable'); %the real parts of the poles that don't have imaginary components ( numel = number of real poles)
    colors = jet(numel(real_parts)+numel(comp_real));
    hold on;
    for i = 1:numel(real_parts)
        groups = eigA(real(eigA) == real_parts(i) & imag(eigA)~=0);
        scatter(real(groups), imag(groups)/(2*pi), [], colors(i,:), 'filled');
        legend_entries{i} = sprintf('[\\lambda_{%d} \\lambda_{%d}] [%s %s %s]', find(real(eigA) == real_parts(i) & imag(eigA)~=0, 1), find(real(eigA) == real_parts(i) & imag(eigA)~=0, 1, 'last'),string_values(find(real(eigA) == real_parts(i) & imag(eigA)~=0, 1),:));
    end

    for i = 1:numel(comp_real)
        gr = eigA(real(eigA) == comp_real(i) & imag(eigA)==0);
        scatter(real(gr), imag(gr)/(2*pi), [], colors(i+numel(real_parts),:), 'filled');
        legend_entries{i+numel(real_parts)} = sprintf('\\lambda_{%d} [%s %s %s]', find(real(eigA) == comp_real(i),1), string_values(find(real(eigA) == comp_real(i),1),:));
    end
    hold off;
    % legend('boxoff')
%     xlim([-6500 1]);
%     ylim([-400 400]);
%     xL = xlim;
%     yL = ylim;
    legend(legend_entries,'FontSize',6,'AutoUpdate','off')
    hold on
    slope = sqrt((1-(zita_min)^2)/zita_min^2);
    xz = linspace(zita_lim, 0, 10000);
    yz = slope * xz/(2*pi);
    plot(xz, yz, '--r'); % 'b' for blue color
    hold on
    plot(xz, -yz, '--r'); % 'b' for blue color
    hold on
    slope = sqrt((1-(0.707)^2)/(0.707)^2);
    xz = linspace(zita_lim, 0, 10000);
    yz = slope * xz /(2*pi);
    plot(xz, yz, '--g'); % 'b' for blue color
    hold on
    plot(xz, -yz, '--g'); % 'b' for blue color

    xline(0, 'Color', 'k');  %x-axis
    yline(0, 'Color', 'k');  %y-axis
    ylabel('\Im [Hz]');
    xlabel('\Re');
    title('Eigenvalues of A');
end