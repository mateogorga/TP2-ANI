n = 10;
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
    resultado(fila,n) = -4;
    resultado(fila,n-1) = 5;
    resultado(fila,n-2) = -4;
    resultado(fila,n-3) = 1;

    for j = 1:(n-4)
      resultado(fila,j) = 0;
    endfor
  endif
  
endfor

for j = 1:(n-1)
  Fi = (2 + 4*((j/n) - ((j/n)^2))) * ((1/n)^4);
  resultado(j,(n+1)) = Fi;
endfor


disp(resultado);

#Eliminacion Gaussiana

for i = 1:n
  for j = (i+1):(n-1)
    resultado(j,:) = resultado(j,:) - (resultado(j,i)/resultado(i,i))*resultado(i,:);
  endfor
endfor
