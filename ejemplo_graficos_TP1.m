cotas = load('cotas_esquinas.dat'); % cota/x/y
escuelas = load('primarias_estatales.dat'); % escuela/partido/x/y

scatter(cotas(:,2), cotas(:,3), 4, cotas(:,1), 'filled')
% scatter(vector x, vector y, tamaño de marcador, color, otros argumentos)
hold on
scatter(escuelas(:,3), escuelas(:,4), 10, 'r')
daspect([1 1])
colormap(gca, 'winter') % otros colormaps en https://octave.sourceforge.io/octave/function/colormap.html

% Grilla
grid

% colorbar
colorbar
ylabel(colorbar, 'cota [m]')

% sacar ticklabels
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

% Rótulos de ejes
xlabel('x', 'fontsize', 16)
ylabel('y','fontsize', 16)

print -djpg ejemplo.jpg

hold off;
