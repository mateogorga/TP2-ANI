#TP2
#a) Armado del SEL
function resultado = armado_matriz(n)
  resultado = [];
  for fila = 1:(n+1)

    if (fila == 1)              #fila 1
      resultado(fila,1) = 1;
        for j = 2:(n+1)         #completo con ceros
         resultado(fila,j) = 0;
        endfor
    endif

    if (fila == 2)              #fila 2
      resultado(fila,1) = -4;
      resultado(fila,2) = 5;
      resultado(fila,3) = -4;
      resultado(fila,4) = 1;
      for j = 5:(n+1)           #completo con ceros
        resultado(fila,j) = 0;
      endfor
    endif

    if (fila > 2)&&(fila < n)    #filas 2 < i < n
      resultado(fila,(fila-2)) = 1;
      resultado(fila,(fila-1)) = -4;
      resultado(fila,(fila)) = 6;
      resultado(fila,(fila+1)) = -4;
      resultado(fila,(fila+2)) = 1;

      if (fila == 3)
        for j = 7:(n+1)
          resultado(fila,j) = 0;
        endfor
        for j = (fila+4):(n+1)
          resultado(fila,j) = 0;
        endfor
      endif
    endif

    if (fila == n)                #fila n
      resultado(fila,n+1) = -4;
      resultado(fila,n)= 5;
      resultado(fila,n-1) = -4;
      resultado(fila,n-2) = 1;
      for j = 1:(n-3)
        resultado(fila,j) = 0;
      endfor
    endif

    if (fila == n+1)              #fila n + 1
      resultado(fila,n+1) = 1;

      for j = 1:n
        resultado(fila,j) = 0;
      endfor
    endif

  endfor
endfunction

A10 = armado_matriz(10);
A50 = armado_matriz(50);
A100 = armado_matriz(100);

function vector = armado_b(n)
  i = 1;
  for j = 1:(n+1)

    if (j == 1)
      vector(j,1) = 0;
    endif

    if (j > 1) && (j < n+1)
      Fi = (2 + 4*((i/n) - ((i/n)^2))) * ((1/n)^4);
      vector(j,1) = Fi;
      i = i + 1;
    endif

    if (j == n+1)
      vector(j,1) = 0;
    endif
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

#Condicion de convergencia para ambos metodos

#Matriz de iteracion para el metodo de Jacobi

function radio_espectral_Jacobi = r_Jacobi(A,b)
  radio_espectral_Jacobi = 0;
#(matrices D, L y U tal que A = D + L + U)
  D = diag(diag(A));
  L = tril(A)-D;
  U = triu(A)-D;
  T_Jacobi = -inv(D)*(L+U); #matriz de iteracion T de Jacobi

    radio_espectral_Jacobi = max(abs(eig(T_Jacobi)));
  if (radio_espectral_Jacobi<1)
    disp('El metodo converge.')
  else disp('El metodo no converge')
  endif

endfunction

r10_Jacobi = r_Jacobi(A10,b10);
r50_Jacobi = r_Jacobi(A50,b50);
r100_Jacobi = r_Jacobi(A100,b100);

#Matriz de iteracion para el metodo de Gauss-Seidel

function radio_espectral_GS = r_GS(A,b)
  radio_espectral_GS = 0;
#(matrices D, L y U tal que A = D + L + U)
  D = diag(diag(A));
  L = tril(A)-D;
  U = triu(A)-D;
  T_GS = -inv(D + L)*U; #matriz de iteracion T de Gauss-Seidel

  radio_espectral_GS = max(abs(eig(T_GS)));
  if (radio_espectral_GS<1)
    disp('El metodo converge.')
  else disp('El metodo no converge')
  endif

endfunction

r10_GS = r_GS(A10,b10);
r50_GS = r_GS(A50,b50);
r100_GS = r_GS(A100,b100);
