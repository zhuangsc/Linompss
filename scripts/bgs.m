#! /usr/bin/octave -qf

###
## Blocked Gauss-Seidel
###

if (length(argv) != 3)
	printf("./bgs.m [Matrix A] [Vector B] [Block size]\n");
	quit
endif

MAX_ITER = 20;
MT = spconvert(load(argv(){1}));
B = load(argv(){2});
BS = str2num(argv(){3});
DIM = rows(MT);
##Calculate the true result
X = linsolve(MT, B);

printf("Dim = %d\n", DIM);

X0 = zeros(length(X), 1);

for iter = 1 : MAX_ITER
	Xtmp = zeros(length(X), 1);
	for i = 1 : BS : DIM
		for j = 1 : BS : DIM
			if ( i != j )
				X0(i:i+BS-1,:) += MT(i:i+BS-1,j:j+BS-1) * X0(j:j+BS-1,:);
			endif
		endfor
		X0(i:i+BS-1,:) = B(i:i+BS-1,:) - X0(i:i+BS-1,:);
		LL = chol(MT(i:i+BS-1,i:i+BS-1),"lower");
		LU = chol(MT(i:i+BS-1,i:i+BS-1),"upper");
		X0(i:i+BS-1,:) = linsolve(LL, X0(i:i+BS-1,:));
		X0(i:i+BS-1,:) = linsolve(LU, X0(i:i+BS-1,:));
	endfor
endfor

#printf("norm x: %e\t norm x0: %e\n", norm(X), norm(X0));
#printf("(n-n0)/n0: %e\n", abs(norm(X0)-norm(X))/norm(X));
printf("B-A*x: %e\n", norm(B-MT*X0));
