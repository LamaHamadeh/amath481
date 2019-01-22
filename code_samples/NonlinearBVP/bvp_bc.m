function bc = bvp_bc(yl,yr,beta)

bc = [yl(1); yl(2)-1; yr(1)];

%bc = [yl(1); yl(2)-0.1; yr(1)];