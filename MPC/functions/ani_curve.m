classdef ani_curve < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        color
        linewidth
        plt
    end
    
    methods
        function obj = ani_curve(color, linewidth)
            %UNTITLED3 Construct an instance of this class
            obj.color = color;
            obj.linewidth = linewidth;
        end
        
        function plot_curve(obj, x, y)
            obj.plt = plot(x,y, 'color', obj.color, 'linewidth', obj.linewidth);
        end
        
        function update_curve(obj, x, y)
            set(obj.plt, 'xdata', x);
            set(obj.plt, 'ydata', y);
        end
            
    end
end