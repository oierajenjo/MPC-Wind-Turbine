classdef ani_cart < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        W
        L
        color
        plt
        linewidth
    end
    
    methods
        function obj = ani_cart(L, W, color, linewidth)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            %obj.x = P(1);
            %obj.y = P(2);
            %obj.theta = angle;
            obj.L = L;
            obj.W = W;
            obj.linewidth = linewidth;
            %obj.v = v;
            obj.color = color;
            %obj.fig = fig;
            %obj.rot_mat = @(theta) [cos(theta) -sin(theta); sin(theta)  cos(theta)];
        end
        
        function plot_cart(obj, P, v, theta)
                
            cart = obj.cart_placement(P,v,theta);
            
            obj.plt = plot(cart(1,:), cart(2,:),'color', obj.color, 'linewidth', obj.linewidth);
        end
        
        function update_cart(obj, P, v, theta)
            cart = obj.cart_placement(P,v,theta);
            set(obj.plt, 'xdata', cart(1,:));
            set(obj.plt, 'ydata', cart(2,:));
        end
        
        function [R] = rotate(obj, theta)
            R = [cos(theta) -sin(theta); sin(theta)  cos(theta)];
        end
        
        function [pos] = position(obj,P)
            ones_v = ones(1, 8);
            pos = [ones_v*P(1); ones_v*P(2)];
        end
        
        function [cart] = cart(obj,v)
            cart = [obj.L/2 0; % tip of velocity
                     obj.L/2 obj.W/2; % north east
                    -obj.L/2 obj.W/2; % north west
                    -obj.L/2 -obj.W/2; % south west 
                    obj.L/2 -obj.W/2; % south east
                    obj.L/2 0; % 
                    0 0;
                    v 0]'; 
        end
        
        function [cart_place] = cart_placement(obj,P, v, theta)
            cart = obj.cart(v);
            R =  obj.rotate(theta);
            offset = obj.position(P);
            
            cart_place = offset + R*cart;
            
        end
        
        
    end
end

