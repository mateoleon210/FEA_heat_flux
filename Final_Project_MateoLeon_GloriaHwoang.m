%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEA final Project
% December 14th, 2021
%
%
% Mateo Leon
% Gloria Hwoang
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% n_mat :       node matrix: contains the coordinates for each of the nodes
% c_mat :       conectivity matrix: contains information on how the triangles for
%               the mesh are formed 
% counter:      Counts the number of essential nodes
% x     :       Location of the boundary
% h     :       element number: this is used to go through each element in c_mat
% i, j, and k:  used to reference each of the three nodes on each element.
%               This is done by using h and 1, 2, and 3 to go to the appropiate node 

% CALL INITILIZING FUNCTION

[n_mat, c_mat, counter, x] = input('please input one of the 4 initializing functions: "coarse_circle", "fine_circle", "coarse_rounded", or "fine_rounded" ');
c_or_r = input('for circle enter 1, for rounded enter 0 :');


%IINITIAL VARIABLES
conductivity = 1;    % conductivity
q= -1;               % prescribed flux


figure(1)  % INIITAL MESH PLOT
for h = 1:length(c_mat)
i=c_mat(h,1);
j=c_mat(h,2);
k=c_mat(h,3);
plot(        [n_mat(i,1) n_mat(j,1)]  , [n_mat(i,2) n_mat(j,2) ], 'b' )
hold on
plot(        [n_mat(j,1) n_mat(k,1)]  , [n_mat(j,2) n_mat(k,2) ], 'b' )
hold on
plot(        [n_mat(i,1) n_mat(k,1)]  , [n_mat(i,2) n_mat(k,2) ], 'b' )
hold on
end 
axis off
axis equal

% display the boundary conditions 
if c_or_r == 1
text(15.4,  0,   '$\bar{T}$=0',    'Interpreter','latex')
text(-7.6,  0,   '$q_x=-1$',       'Interpreter','latex')
text(5,     5.5, '$\bar{q} = 0$',  'Interpreter','latex')
text(5,    -5.5, '$\bar{q} = 0$',  'Interpreter','latex')
else
text(10.5,   0,   '$\bar{T}$=0',    'Interpreter','latex')
text(-12.6,  0,   '$q_x=-1$',       'Interpreter','latex')
text(0,      5.5, '$\bar{q} = 0$',  'Interpreter','latex')
text(0,     -5.5, '$\bar{q} = 0$',  'Interpreter','latex')
end



% INITIALIZING MATRICES for speed
K_matrix_global = zeros(length(n_mat));                  % k global matrix
f_gamma_local = zeros(3,1);                              % f gamma local
f_gamma_global = zeros(length(n_mat),1);                 % f gamma global
temp_global = zeros(length(n_mat),1);                    % matrix of temperatures
b_tensor_local = zeros(2,3,length(c_mat));

for h = 1:length(c_mat)

i=c_mat(h,1);
j=c_mat(h,2);
k=c_mat(h,3);

% COMPUTE AREAS

two_area= (n_mat(j,1)*n_mat(k,2)-n_mat(k,1)*n_mat(j,2))-(n_mat(i,1)*n_mat(k,2)-n_mat(k,1)*n_mat(i,2))+(n_mat(i,1)*n_mat(j,2)-n_mat(j,1)*n_mat(i,2));

% B MATRIX

b_mat =  (1/two_area)* [(n_mat(j,2)-n_mat(k,2)),(n_mat(k,2)-n_mat(i,2)), (n_mat(i,2)-n_mat(j,2));
                        (n_mat(k,1)-n_mat(j,1)),(n_mat(i,1)-n_mat(k,1)), (n_mat(j,1)-n_mat(i,1))];

% storing b matrices to find fluxes later
b_tensor_local(:,:,h)= b_mat;


% local stiffness matrix
k_mat_local= (two_area/2) *conductivity* transpose(b_mat) * b_mat;

    for a = 1:3
        for b = 1:3
    K_matrix_global(c_mat(h,a) , c_mat(h,b)) = K_matrix_global(c_mat(h,a) , c_mat(h,b)) + k_mat_local(a,b);
        end
    end    

% ELEMENT BOUNDARY FLUX MATRIX
    
    % shape functions 
    N1= @(y) (1/two_area) * (n_mat(j,1)*n_mat(k,2) - n_mat(k,1)*n_mat(j,2) + (n_mat(j,2)-n_mat(k,2))*x + (n_mat(k,1)-n_mat(j,1))*y);
    N2= @(y) (1/two_area) * (n_mat(k,1)*n_mat(i,2) - n_mat(i,1)*n_mat(k,2) + (n_mat(k,2)-n_mat(i,2))*x + (n_mat(i,1)-n_mat(k,1))*y);
    N3= @(y) (1/two_area) * (n_mat(i,1)*n_mat(j,2) - n_mat(j,1)*n_mat(i,2) + (n_mat(i,2)-n_mat(j,2))*x + (n_mat(j,1)-n_mat(i,1))*y);

    if       n_mat(i,1) ==x  && n_mat(j,1) ==x
        y_array=[n_mat(i,2), n_mat(j,2)];
        ymin= min(y_array);
        ymax= max(y_array);
        f_gamma_global(i) = f_gamma_global(i)  -q* integral(N1, ymin, ymax);
        f_gamma_global(j) = f_gamma_global(j)  -q* integral(N2, ymin, ymax);  

    elseif   n_mat(j,1) ==x  && n_mat(k,1) ==x
        y_array=[n_mat(j,2), n_mat(k,2)];
        ymin= min(y_array);
        ymax= max(y_array); 
        f_gamma_global(j) = f_gamma_global(j) -q* integral(N2, ymin, ymax);
        f_gamma_global(k) = f_gamma_global(k) -q* integral(N3, ymin, ymax);  

    elseif   n_mat(i,1) ==x  && n_mat(k,1) ==x
        y_array=[n_mat(i,2), n_mat(k,2)];
        ymin= min(y_array);
        ymax= max(y_array);
        f_gamma_global(i) = f_gamma_global(i) -q* integral(N1, ymin, ymax);
        f_gamma_global(k) = f_gamma_global(k) -q* integral(N3, ymin, ymax);  
    end
end


% MATRICES
    % PARTITION
        % STIFFNESS
K_E    = K_matrix_global(1:counter,1:counter);
K_EF   = K_matrix_global(1:counter,counter+1 :length(n_mat));
K_EF_T = K_matrix_global(counter+1:length(n_mat) , 1:counter);
K_F    = K_matrix_global(counter+1:length(n_mat) , counter+1 :length(n_mat));

        % FLUX
f_E    = f_gamma_global(1:counter);
f_F    = f_gamma_global(counter+1 :length(n_mat));

        % TEMPERATURE
T_E    =  zeros(counter,1);

   
    % SOLVE THE TEMPERATURE MATRIX
T_F  = K_F\(f_F - K_EF_T * T_E);
T =  cat(1,T_E,T_F);


% FLUX MATRIX GRAPH
element_midpoint_mat=zeros(length(c_mat),2);
q_transpose_global=zeros(length(c_mat),2);

for h= 1:length(c_mat)
local_temp_mat = [ T(c_mat(h,1)); T(c_mat(h,2)); T(c_mat(h,3))];  
q= - conductivity * b_tensor_local(:,:,h) * local_temp_mat;
q_transpose= transpose(q);
q_transpose_global(h,:) = q_transpose;

midpoint_mat_x= (n_mat(c_mat(h,1),1) + n_mat(c_mat(h,2),1)+ n_mat(c_mat(h,3),1))/3;
midpoint_mat_y= (n_mat(c_mat(h,1),2) + n_mat(c_mat(h,2),2)+ n_mat(c_mat(h,3),2))/3;
midpoint = [midpoint_mat_x,midpoint_mat_y];
element_midpoint_mat(h,:) = midpoint;
end

% ARROW FLUX PLOT
figure(2) 
for h = 1:length(c_mat)
i=c_mat(h,1);
j=c_mat(h,2);
k=c_mat(h,3);
plot(        [n_mat(i,1) n_mat(j,1)]  , [n_mat(i,2) n_mat(j,2) ], 'b' )
hold on
plot(        [n_mat(j,1) n_mat(k,1)]  , [n_mat(j,2) n_mat(k,2) ], 'b' )
hold on
plot(        [n_mat(i,1) n_mat(k,1)]  , [n_mat(i,2) n_mat(k,2) ], 'b' )
hold on
end 
axis off
axis equal
quiver(element_midpoint_mat(:,1), element_midpoint_mat(:,2), q_transpose_global(:,1), q_transpose_global(:,2),1,'r')

% display the boundary conditions 
if c_or_r == 1
text(15.4,  0,   '$\bar{T}$=0',    'Interpreter','latex')
text(-7.6,  0,   '$q_x=-1$',       'Interpreter','latex')
text(5,     5.5, '$\bar{q} = 0$',  'Interpreter','latex')
text(5,    -5.5, '$\bar{q} = 0$',  'Interpreter','latex')
else
text(10.5,   0,   '$\bar{T}$=0',    'Interpreter','latex')
text(-12.6,  0,   '$q_x=-1$',       'Interpreter','latex')
text(0,      5.5, '$\bar{q} = 0$',  'Interpreter','latex')
text(0,     -5.5, '$\bar{q} = 0$',  'Interpreter','latex')
end


% TEMPERATURE COLOR PLOT
figure(3) 
for h = 1:length(c_mat)
i=c_mat(h,1);
j=c_mat(h,2);
k=c_mat(h,3);
x_vec = [n_mat(i,1); n_mat(j,1); n_mat(k,1)];
y_vec = [n_mat(i,2); n_mat(j,2); n_mat(k,2)];
T_vec = [T(i,1); T(j,1); T(k,1)];

patch(x_vec, y_vec, T_vec)
hold on
end 

axis equal
colorbar('Location','eastoutside')
axis off


% display the boundary conditions 
if c_or_r == 1

text(15.1,  0,   '$\bar{T}$=',    'Interpreter','latex', FontSize=9)
text(15.3,  -0.8, '0',     'Interpreter','latex', FontSize=10)
text(-7.6,  0,   '$q_x=-1$',       'Interpreter','latex')
text(5,     5.5, '$\bar{q} = 0$',  'Interpreter','latex')
text(5,    -5.5, '$\bar{q} = 0$',  'Interpreter','latex')
else
text(10.1,   0,   '$\bar{T}$=',    'Interpreter','latex', FontSize=9)
text(10.3,  -0.8, '0',     'Interpreter','latex', FontSize=10)
text(-12.6,  0,   '$q_x=-1$',       'Interpreter','latex')
text(0,      5.5, '$\bar{q} = 0$',  'Interpreter','latex')
text(0,     -5.5, '$\bar{q} = 0$',  'Interpreter','latex')
end
