%% ISOPARAMETRIC
[n_mat, c_mat,counter, x ]=coarse_circle; 

% gauss quadrature
s = 0;           % boundary condition
q = -1;          % prescribed flux
conductivity = 1;    % conductivity

 GN = [1 0 -1;
       0 1 -1];


K_matrix_global_iso= zeros(length(n_mat),length(n_mat));
f_matrix = zeros(length(n_mat),1);


for h= 1: length(c_mat)
    i=c_mat(h,1);
    j=c_mat(h,2);
    k=c_mat(h,3);
    
    J = [n_mat(i,1)-n_mat(k,1),   n_mat(i,2)-n_mat(k,2);
         n_mat(j,1)-n_mat(k,1),   n_mat(j,2)-n_mat(k,2)];

    b_mat= J\GN;
    b_tensor_local(:,:,h)= b_mat;
    J_tensor_local(h,1)= det(J);
    b_mat_t=transpose(b_mat);
    
    k_mat_local =  0.5 * conductivity * det(J) * b_mat_t * b_mat;

     for a = 1:3
        for b = 1:3
            K_matrix_global_iso(c_mat(h,a) , c_mat(h,b)) = K_matrix_global_iso(c_mat(h,a) , c_mat(h,b)) + k_mat_local(a,b);
        end
     end    


    %b_tensor_local(:,:,h)= b_mat;


    % FLUX

    t= 1/2;

    if       n_mat(i,1) ==x  && n_mat(j,1) ==x
  
    f_matrix(i) = f_matrix(i)  -q* 0.5 * (s)        * 2;
    f_matrix(j) = f_matrix(j)  -q* 0.5 * (t)        * 2;
    f_matrix(k) = f_matrix(k)  -q* 0.5 * (1-s-t)    * 2; 

    elseif   n_mat(j,1) ==x  && n_mat(k,1) ==x
  
    f_matrix(i) = f_matrix(i)  -q* 0.5 * (s)        * 2;
    f_matrix(j) = f_matrix(j)  -q* 0.5 * (t)        * 2;
    f_matrix(k) = f_matrix(k)  -q* 0.5 * (1-s- t)   * 2; 

    elseif   n_mat(i,1) ==x  && n_mat(k,1) ==x

    f_matrix(i) = f_matrix(i)  -q* 0.5 *  (s)       * 2;
    f_matrix(j) = f_matrix(j)  -q* 0.5 *   (t)      * 2;
    f_matrix(k) = f_matrix(k)  -q* 0.5 * (1-s-t)    * 2; 
    end
end


disp('flux matrix is:')
sum(f_matrix, 'all')


% MATRICES
    % PARTITION
        % STIFFNESS
K_E    = K_matrix_global_iso(1:counter,1:counter);
K_EF   = K_matrix_global_iso(1:counter,counter+1 :length(n_mat));
K_EF_T = K_matrix_global_iso(counter+1:length(n_mat) , 1:counter);
K_F    = K_matrix_global_iso(counter+1:length(n_mat) , counter+1 :length(n_mat));

        % FLUX
f_E    =  f_matrix(1:counter);
f_F    =  f_matrix(counter+1 :length(n_mat));

        % TEMPERATURE
T_E    =  zeros(counter,1);

   
    % SOLVE THE TEMPERATURE MATRIX

T_F  = K_F\(f_F - K_EF_T * T_E);

T_iso =  cat(1,T_E,T_F);


% PLOT THE TEMPERATURE DISTRIBUTION 

for h = 1:length(c_mat)
i=c_mat(h,1);
j=c_mat(h,2);
k=c_mat(h,3);
x_vec = [n_mat(i,1); n_mat(j,1); n_mat(k,1)];
y_vec = [n_mat(i,2); n_mat(j,2); n_mat(k,2)];
T_vec_iso = [T_iso(i,1); T_iso(j,1); T_iso(k,1)];

patch(x_vec, y_vec, T_vec_iso)
hold on
end 

axis equal
colorbar
axis off

text(15.1,  0,   '$\bar{T}$=',    'Interpreter','latex', FontSize=9)
text(15.3,  -0.8,   '0',          'Interpreter','latex', FontSize=10)
text(-7.6,  0,   '$q_x=-1$',       'Interpreter','latex')
text(5,     5.5, '$\bar{q} = 0$',  'Interpreter','latex')
text(5,    -5.5, '$\bar{q} = 0$',  'Interpreter','latex')

% PLOT THE TEMPERATURE DISTRIBUTION


% FLUX MATRIX GRAPH
element_midpoint_mat=zeros(length(c_mat),2);
q_transpose_global=zeros(length(c_mat),2);

for h= 1:length(c_mat)

local_temp_mat = [ T_iso(c_mat(h,1)); T_iso(c_mat(h,2)); T_iso(c_mat(h,3))];  

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
quiver(element_midpoint_mat(:,1), element_midpoint_mat(:,2), q_transpose_global(:,1), q_transpose_global(:,2), 1, 'r')

% display the boundary conditions 
text(15.4,  0,   '$\bar{T}$=0',    'Interpreter','latex')
text(-7.6,  0,   '$q_x=-1$',       'Interpreter','latex')
text(5,     5.5, '$\bar{q} = 0$',  'Interpreter','latex')
text(5,    -5.5, '$\bar{q} = 0$',  'Interpreter','latex')
