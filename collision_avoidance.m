function [force_val]=collision_avoidance(r,r_min, L_f)
%% collision avoidance function
    p = r;
    d = r_min;
    
    if (norm(p) < d)
        if ( dot(r(1:3), L_f) < 0)
            L_f = - L_f;
        end
    end
    
    force_val = L_f;
end