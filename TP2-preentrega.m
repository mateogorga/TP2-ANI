#TP2
#a) Armado del SEL
function resultado = armado_matriz(n)
  resultado = [];
  for fila = 1:(n-1)
    if (fila == 1)#fila 1
      resultado(fila,1) = -4;
      resultado(fila,2) = 5;
      resultado(fila,3) = -4;
      resultado(fila,4) = 1;
      for j = 5:(n+1)#completo con ceros
        resultado(fila,j) = 0;
      endfor
    endif
    
    if (fila > 1)&&(fila < (n-1))#filas 1 < i < n-1
      resultado(fila,(fila-1)) = 1;
      resultado(fila,(fila)) = -4;
      resultado(fila,(fila+1)) = 6;
      resultado(fila,(fila+2)) = -4;
      resultado(fila,(fila+3)) = 1;
      
      if (fila == 2)
        for j = 6:(n+1)
          resultado(fila,j) = 0;
        endfor
      else
        for j = 1:(fila-2)
          resultado(fila,j) = 0;
        endfor
        for j = (fila+4):(n+1)
          resultado(fila,j) = 0;
        endfor
      endif
    endif
    
    if (fila == (n-1))#fila n-1
      resultado(fila,n+1) = -4;
      resultado(fila,n) = 5;
      resultado(fila,n-1) = -4;
      resultado(fila,n-2) = 1;

      for j = 1:(n-3)
        resultado(fila,j) = 0;
      endfor
    endif
    
  endfor
endfunction

A10 = armado_matriz(10);
A50 = armado_matriz(50);
A100 = armado_matriz(100); 

A10 = A10(:,2:10);
A50 = A50(:,2:50);
A100 = A100(:,2:100);
  
function vector = armado_b(n)
  j = 1
  for j = 1:(n-1)
    Fi = (2 + 4*((j/n) - ((j/n)^2))) * ((1/n)^4);
    vector(j,1) = Fi
    j = j + 1
  endfor
endfunction
  
b10 = armado_b(10);
b50 = armado_b(50);
b100 = armado_b(100);

Ab10 = [A10 b10];
Ab50 = [A50 b50];
Ab100 = [A100 b100];


#b) 1) Eliminacion Gaussiana

function x = EG(Ab,n)
    #Triangulacion de la matriz

  for i = 1:n
    for j = (i+1):(n-1)
      Ab(j,:) = Ab(j,:) - (Ab(j,i)/Ab(i,i))*Ab(i,:);
    endfor
  endfor

   #Vector solucion

  x = zeros(n-1,1);

  for i = (n-1):-1:1
    suma = 0;
    if i<=(n-2)
      for j = i+1:n-1
        suma = suma + Ab(i,j)*x(j)
      endfor
    endif
    x(i) = (Ab(i,n) - suma)/Ab(i,i);
  endfor

endfunction

x10_EG = EG(Ab10,10);
x50_EG = EG(Ab50,50);
x100_EG = EG(Ab100,100);



# 2) Jacobi
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

x10_Jacobi = Jacobi(A10,b10);
x50_Jacobi = Jacobi(A50,b50);
x100_Jacobi = Jacobi(A100,b100);

# 3) Gauss Seidel
function x = GS(A,b)
  x = zeros(length(b),1);
  n = size(x,1);
  normVal = Inf;
  tol = 0.01;
  itr = 0;
  while normVal > tol
    x_old = x;
    
    for i=1:n
      suma = 0;
      for j=1:i-1
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

x10_GS = GS(A10,b10);
x50_GS = GS(A50,b50);
x100_GS = GS(A100,b100);

#Defino matrices D, L y U tal que A = D + L + U

#Matriz de iteracion para el metodo de Jacobi

function radio_espectral_Jacobi = r(A,b)
  radio_espectral_Jacobi = 0;
  D = diag(diag(A));
  L = tril(A)-D;
  U = triu(A)-D;
  T_j = -inv(D)*(L+U); #matriz de iteracion T de Jacobi

  #Chequeo condicion de convergencia
  radio_espectral_Jacobi = max(abs(eig(T_j)));
  if radio_espectral_Jacobi<1
    disp('El metodo converge.')
  else disp('El metodo no converge')
  endif

endfunction

r10_Jacobi = r(A10,b10);
r50_Jacobi = r(A50,b50);
r100_Jacobi = r(A100,b100);

#Matriz de iteracion para el metodo de Gauss-Seidel

function radio_espectral_GS = s(A,b)
  radio_espectral_GS = 0;
  D = diag(diag(A));
  L = tril(A)-D;
  U = triu(A)-D;
  T_gs = -inv(D + L)*U; #matriz de iteracion T de Gauss-Seidel

  #Chequeo condicion de convergencia
  radio_espectral_GS = max(abs(eig(T_gs)));
  if radio_espectral_GS<1
    disp('El metodo converge.')
  else disp('El metodo no converge')
  endif

endfunction

r10_GS = s(A10,b10);
r50_GS = s(A50,b50);
r100_GS = s(A100,b100);
