n = 10;
resultado = [];


for fila = 1:(n+1)
  if (fila == 1)#fila 1
    resultado(fila,1) = 1;
    for j = 2:(n+2)
      resultado(fila,j) = 0;
    endfor
  elseif (fila == 2)
    resultado(fila,1) = -4;
    resultado(fila,2) = 5;
    resultado(fila,3) = -4;
    resultado(fila,4) = 1;
    for j = 5:(n+1)#completo con ceros
      resultado(fila,j) = 0;
    endfor
  endif
  
  if (fila > 2)&&(fila < (n))#filas 1 < i < n-1
    resultado(fila,(fila-2)) = 1;
    resultado(fila,(fila-1)) = -4;
    resultado(fila,(fila)) = 6;
    resultado(fila,(fila+1)) = -4;
    resultado(fila,(fila+2)) = 1;
    
    if (fila == 3)
      for j = 6:(n+1)
        resultado(fila,j) = 0;
      endfor
    else
      for j = 2:(fila-2)
        resultado(fila,j-1) = 0;
      endfor
      for j = (fila+4):(n+1)
        resultado(fila,j) = 0;
      endfor
    endif
  endif
  
  if (fila == (n))#fila n-1
    resultado(fila,n) = -4;
    resultado(fila,n-1) = 5;
    resultado(fila,n-2) = -4;
    resultado(fila,n-3) = 1;

    for j = 1:(n-4)
      resultado(fila,j) = 0;
    endfor
  endif
  
  if (fila == n+1)
    for j = 1:n
      resultado(fila,j) = 0;
    endfor
    resultado(fila,n+1) = 1;
    resultado(fila,n+2) = 0;
  endif
    
  
endfor

for j = 1:(n-1)
  Fi = (2 + 4*((j/n) - ((j/n)^2))) * ((1/n)^4);
  resultado(j+1,(n+2)) = Fi;
endfor


disp(resultado);

for fila = 1:n+1
  for columna = 1:n+1
    A(fila,columna) = resultado(fila,columna);
  endfor
endfor  
for fila = 1:n+1
  b(fila,1) = resultado(fila,n+2);
endfor
  
disp(A);
disp(b);


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

disp(itr);
  
disp(x);