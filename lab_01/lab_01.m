function lab_1()
    function [amnt] = count_less(t, X)
        tolerance = 1e-4;
        amnt = 0;
        for x_index = 1:length(X)
            if (X(x_index) - t < tolerance)
                amnt = amnt + 1;
            end
        end
    end
    
    X = [11.89,9.60,9.29,10.06,9.50,8.93,9.58,6.81,8.69,9.62,...
        9.01,10.59,10.50,11.53,9.94,8.84,8.91,6.90,9.76,7.09,...
        11.29,11.25,10.84,10.76,7.42,8.49,10.10,8.79,11.87,...
        8.77,9.43,12.41,9.75,8.53,9.72,9.45,7.20,9.23,8.93,...
        9.15,10.19,9.57,11.09,9.97,8.81,10.73,9.57,8.53,9.21,...
        10.08,9.10,11.03,10.10,9.47,9.72,9.60,8.21,7.78,10.21,...
        8.99,9.14,8.60,9.14,10.95,9.33,9.98,9.09,10.35,8.61,...
        9.35,10.04,7.85,9.64,9.99,9.65,10.89,9.08,8.60,7.56,...
        9.27,10.33,10.09,8.51,9.86,9.24,9.63,8.67,8.85,11.57,...
        9.85,9.27,9.69,10.90,8.84,11.10,8.19,9.26,9.93,10.15,...
        8.42,9.36,9.93,9.11,9.07,7.21,8.22,9.08,8.88,8.71,...
        9.93,12.04,10.41,10.80,7.17,9.00,9.46,10.42,10.43,8.38,9.01];
   
    X = sort(X);
    
    m_max = max(X);
    m_min = min(X);

    fprintf("\n1)M_max = %f.\n  M_min = %f.\n\n", m_max, m_min);
    
    R = m_max - m_min;
    fprintf("\n2)Размах выборки: R = %f.\n\n", R);

    
    n = length(X);
    q = sum(X) / n;
    S2 = sum((X - q).^2) / (n - 1);
    fprintf("\n3)Оценка математического ожидания: q = %f.\n", q);
    fprintf("   Оценка дисперсии: S^2 = %f.\n", S2);
    
    m = floor(log2(n)) + 2;
    groups = [];
    curr = m_min;
    d = R / m;

    for i = 1:(m + 1)
        groups(i) = curr;
        curr = curr + d;
    end

    eps = 1e-6;
    amounts = [];

    tot = 0;

    for i = 1:(m - 1)
        curr = 0;

        for j = 1:n
            if ((X(j) - groups(i)) > eps || abs(groups(i) - X(j)) < eps) && (X(j) - groups(i + 1) < eps)
                curr = curr + 1;
            end
        end

        amounts(i) = curr;
        tot = tot + curr;
    end

    amounts(m) = n - tot;

    fprintf("\n4)Группировка значений выборки в %d интервалов:\n", m);
    for i = 1:(m)
        fprintf("Интервал %d [%f : %f) - %d значений из выборки.\n", i, groups(i), groups(i + 1), amounts(i));
    end


    fprintf("\n5)Построение гистограммы и графика функции плотности распределения нормальной СВ.\n");
    figure;
    hold on;
    grid on;

    cent = zeros(1, m);
    heights = zeros(1, m);

    for i = 1:m
        heights(i) = amounts(i) / (n * d);
    end

    for i = 1:m
        cent(i) = groups(i + 1) - (d / 2);
    end

    set(gca, "xtick", groups);
    set(gca, "ytick", unique(sort(heights)));
    set(gca, "xlim", [min(groups) - 1, max(groups) + 1]);

    bar(cent, heights, 1);
    points = (m_min - 5):(S2 / 250):(m_max + 5);
    X_pdf = normpdf(points, q, sqrt(S2));
    plot(points, X_pdf, "r", 'LineWidth', 2);

    xlabel('x')
    ylabel('f(x)')
    hold off;

    fprintf("\n6)Построение графика эмпирической функции распределения и функции распределения нормальной случайной величины.\n");
    figure;
    hold on;
    grid on;
    
    x = unique(X);
    x(end + 1) = x(end) + 1;
    heights = zeros([1 length(x)]);
    for i = 1:length(x)
        heights(i) = count_less(x(i), X) / n;
    end

    heightss = unique(heights);
    
    set(gca, "ylim", [0, max(heightss) + 0.3]);
    stairs(x, heights);

    points = (m_min - 2):(S2 / 250):(m_max + 2);
    x_cdf = normcdf(points, q, sqrt(S2));
    plot(points, x_cdf, "r");

    xlabel('x')
    ylabel('F(x)')
    hold off;
end

