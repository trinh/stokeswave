function [X, Y, h] = ginputmouse

    but = 1;
    tf = ishold;
    hold on
    X = [];
    Y = [];
    h = [];
    while but == 1
        [X(end+1), Y(end+1), but] = ginput(1);
        h(end+1) = plot(X(end), Y(end), 'b.', 'MarkerSize', 25);
    end
    if tf == 0 % Was not originally held
        hold off
    end

    keyboard
end

