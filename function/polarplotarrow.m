function polarplotarrow(orientation,resultant_length,arrowhead_angle,arrowhead_scale,num_arrowlines)

if nargin == 2
    arrowhead_angle = 15;
    arrowhead_scale = 0.15;
    num_arrowlines  = 100;
end


resultant_direction = deg2rad(orientation);

%%%% arrow head %%%%
arrowhead_length = resultant_length*arrowhead_scale; % arrow head length relative to resultant_length
arrowhead_radian  = deg2rad(arrowhead_angle);

%%%% arrow tip coordinates %%%%
t1 = repmat(resultant_direction,1,num_arrowlines);
r1 = repmat(resultant_length,1,num_arrowlines);

%%%% arrow base coordinates %%%%
b = arrowhead_length.*tan(linspace(0,arrowhead_radian,num_arrowlines/2));
theta = atan(b./(resultant_length-arrowhead_length));
pre_t2 = [theta, -theta];
r2 = (resultant_length-arrowhead_length)./cos(pre_t2);
t2 = t1(1)+pre_t2;

%%%% plot %%%%
polarplot([t1(1) t1(1)],[0 r1(1)-0.9*arrowhead_length],'r','linewidth',2); % plot arrow head
hold on
polarplot([t1; t2],[r1; r2],'r');

return