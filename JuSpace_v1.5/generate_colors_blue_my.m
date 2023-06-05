function colors = generate_colors_blue_my(N)
% function colors = generate_colors_nice_my(N)
% from pa_LCH2RGB

red = (0:1:N-1)'./N;

blue = 0.5+((N-1:-1:0)'./N)./2;

colors = [zeros(length(red),1)+0.3 blue./1.2 blue];
