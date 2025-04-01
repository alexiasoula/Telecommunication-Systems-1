function DrawSeminology(x, y, titleName, xLabel, yLabel)
    semilogy(x,y);
    grid on;
    title(titleName);
    xlabel(xLabel);
    ylabel(yLabel);
end

