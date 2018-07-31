function ContourPlot

    x = linspace(-5, 5);
    y = linspace(-5, 5);

    [X,Y] = meshgrid(x,y);
    f = (X.^2 + Y - 11).^2 + (X + Y.^2 - 7).^2;   
    
    a = 0.01;
    
    F = log(a + f);

    figure
    contour(X, Y, F)
    hold on
    plot(-2.8, 3.1, '*', 'Color', 'black') 
    hold on
    plot(-3.8, -3.3, '*', 'Color', 'black')
    hold on
    plot(3.6, -1.9, '*', 'Color', 'black')
    hold on
    plot(3.0, 2.0, '*', 'Color', 'black')   
    
    title('f(x,y)')
    xlabel('x');
    ylabel('y');
    
end