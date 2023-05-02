A = [4,3,0;
     3,4,-1;
     0,-1,4];

b = [24,30,-24];

#Gauss Seidel
function x = GS(A,b)
  x = zeros(length(b),1);
  normVal = Inf;
  tol = 0.01;
  itr = 0;
  while normVal > tol
    x_old = x;
    
    for i=1:n
      suma = 0;
      for j=1:i-1
  n = size(x,1);
        suma = suma+A(i,j)*x(j);
      endfor 
      for j=i+1:n
        suma = suma+A(i,j)*x_old(j);
      endfor
      
      x(i) =(b(i)-suma)*(1/A(i,i));
    endfor
    itr = itr + 1;
    normVal = norm(x_old-x);
  endwhile
  
  disp(itr)
endfunction

x_GS = GS(A,b);

function x = Jacobi(A,b)
  x = zeros(length(b),1);
  n = size(x,1);
  normVal = Inf;
  tol = 0.01;
  itr = 0;
  while normVal > tol
    x_old = x;
    
    for i=1:n
      suma = 0;
      for j=1:n
        if j!=i
         suma = suma + A(i,j)*x_old(j);
        endif
      endfor 
      x(i) =(b(i)-suma)*(1/A(i,i));
    endfor
    
    itr = itr + 1;
    normVal = norm(x_old-x);
  endwhile
  
  disp(itr)
  
endfunction

vec_x_Jacobi = Jacobi(A,b)