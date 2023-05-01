n = 10;
resultado = [];
tol = 0.01;
k = 0;
er = 1;

for fila = 1:(n-1)
  if (fila == 1)#fila 1
    resultado(fila,1) = -4;
    resultado(fila,2) = 5;
    resultado(fila,3) = 6;
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
  resultado(j,(n+2)) = Fi;
endfor


#disp(resultado);

#Metodo de Jacobi
#1° iteracion:
rta_parciales(1,1) = resultado(1,(n+2))/resultado(1,2);

for fila = 2:(n-2)
  rta_parciales(fila,1) = resultado(fila,(n+2))/resultado(fila,fila);
endfor

rta_parciales(n-1,1) = resultado(n-1,n+2)/resultado(n-1,n);
k = k + 1;



while er!=0
  er = 0;
  for fila = 1:n-1
    if fila == 1
      rta_nuevas(1,1) = rta_parciales(1,1)*5 - rta_parciales(2,1)*4 + rta_parciales(3,1);
    elseif fila>=2&&fila<=n-3
      rta_nuevas(fila,1) = rta_parciales(fila-1,1)*(-4) + rta_parciales(fila,1)*6 + rta_parciales(fila+1,1)*(-4) + rta_parciales(fila+2,1);
    elseif fila == n-2
      rta_nuevas(fila,1) = rta_parciales(fila-1,1)*(-4) + rta_parciales(fila,1)*6 + rta_parciales(fila+1,1)*(-4);
    elseif fila == n-1
      rta_nuevas(fila,1) = rta_parciales(fila-3,1) - rta_parciales(fila-2,1)*4 + rta_parciales(fila-1,1)*5 - rta_parciales(fila,1)*4;
    endif
  endfor
  
  for r = 1:n-1#actualizo el vector
    er = abs(rta_nuevas(r,1) - rta_parciales(r,1));
    if er > tol
      er = er + 1;
    endif
    
    rta_parciales(r,1) = rta_nuevas(r,1);
  endfor
  
  k = k+1;

endwhile  
  
disp(rta_parciales);
disp(k);